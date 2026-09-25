! #################################################################################################################################
! Begin MIT license text.
! _______________________________________________________________________________________________________

! Copyright 2022 Dr William R Case, Jr (mystransolver@gmail.com)

! Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
! associated documentation files (the "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to
! the following conditions:

! The above copyright notice and this permission notice shall be included in all copies or substantial
! portions of the Software and documentation.

! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
! OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
! THE SOFTWARE.
! _______________________________________________________________________________________________________

! End MIT license text.
SUBROUTINE MITC4_BMBS ( R, S, BM, BB, BS )

  USE PENTIUM_II_KIND, ONLY   : DOUBLE
  USE MODEL_STUF, ONLY    : ELGP, EPROP
  USE CONSTANTS_1, ONLY       : ZERO, HALF, ONE

  USE MITC4_B_Interface
  USE MITC4_CARTESIAN_LOCAL_BASIS_Interface
  USE MITC_COVARIANT_BASIS_Interface
  USE MITC_TRANSFORM_B_Interface

  IMPLICIT NONE

  REAL(DOUBLE), INTENT(IN)    :: R, S
  REAL(DOUBLE), INTENT(OUT)   :: BM(3, 6*ELGP)
  REAL(DOUBLE), INTENT(OUT)   :: BB(3, 6*ELGP)
  REAL(DOUBLE), INTENT(OUT)   :: BS(2, 6*ELGP)

  REAL(DOUBLE)        :: BMEM(6, 6*ELGP)
  REAL(DOUBLE)        :: BBOT(6, 6*ELGP)
  REAL(DOUBLE)        :: BTOP(6, 6*ELGP)
  REAL(DOUBLE)        :: BSHR(6, 6*ELGP)
  REAL(DOUBLE)        :: G(3,3)            ! Covariant basis at the mid-surface, in element coordinates
  REAL(DOUBLE)        :: M1, M2            ! In-plane slope of the director relative to the element facet

! **********************************************************************************************************************************
! Subr MITC4_B returns the strain-displacement matrix in the CARTESIAN LOCAL basis, whose x axis lies along the covariant g_r
! direction. Everything that uses the operators returned here - the elasticity matrices built from EM, EB and ET, the stress,
! strain and engineering force output, and the extrapolation of Gauss point values to the corners - works in the ELEMENT
! coordinate system. So the rows have to be rotated out of the cartesian local basis before they are handed back.
!
! For a rectangular element the two systems differ by 180 degrees about z, because the element x axis starts along side 1-2
! while the cartesian local x axis lies along g_r, and side 1-2 runs in the negative r direction in Bathe's node ordering.
!
! MITC4_B doubles rows 4 to 6 on the way out to make them engineering shear strains. MITC_TRANSFORM_B rotates tensor
! components, so those rows are halved before the rotation and doubled again afterwards. This matches what subr MITC4 on the
! dev branch does at each of these points.

! **********************************************************************************************************************************
! Pure midsurface membrane operator.

  CALL MITC4_B( R, S, ZERO, .TRUE., .FALSE., .FALSE., BMEM )
  CALL TO_ELEMENT_BASIS( ZERO, BMEM )

  BM(1,:) = BMEM(1,:)  ! xx
  BM(2,:) = BMEM(2,:)  ! yy
  BM(3,:) = BMEM(4,:)  ! xy

! **********************************************************************************************************************************
! Pure curvature operator.
! Use bending-only calls so membrane terms cannot leak into curvature.
! The difference removes even-in-T bending terms.

  CALL MITC4_B( R, S, -ONE, .FALSE., .TRUE., .FALSE., BBOT )
  CALL TO_ELEMENT_BASIS( -ONE, BBOT )

  CALL MITC4_B( R, S, +ONE, .FALSE., .TRUE., .FALSE., BTOP )
  CALL TO_ELEMENT_BASIS( +ONE, BTOP )

  BB(1,:) = (BTOP(1,:) - BBOT(1,:)) / EPROP(1)  ! kxx-ish
  BB(2,:) = (BTOP(2,:) - BBOT(2,:)) / EPROP(1)  ! kyy-ish
  BB(3,:) = (BTOP(4,:) - BBOT(4,:)) / EPROP(1)  ! kxy-ish

! **********************************************************************************************************************************
! Pure transverse shear operator.
! MITC shear is already evaluated at tying points and is effectively midsurface shear here.
! Keep row order as zx, yz to match SHELL_T convention used by MITC4.f90.

  CALL MITC4_B( R, S, ZERO, .FALSE., .FALSE., .TRUE., BSHR )
  CALL TO_ELEMENT_BASIS( ZERO, BSHR )

  BS(1,:) = BSHR(6,:)  ! zx
  BS(2,:) = BSHR(5,:)  ! yz

! **********************************************************************************************************************************
! Remove the spurious membrane to transverse shear coupling that a director which is not normal to the element facet produces.
!
! The reference geometry of the degenerated shell is X(r,s,t) = Xbar(r,s) + (t*h/2) * V, so when the director V is not normal to
! the facet the through-thickness fibre is slanted and a point at height z sits in-plane offset by z*m, where m is the in-plane
! slope of V in element coordinates. Straining the mid-surface then drags the top of the fibre relative to the bottom by
! (du/dx) * m * z, which the covariant strain registers as transverse shear even though the fibre has not rotated relative to the
! facet at all. In element coordinates that spurious shear is exactly the in-plane strain contracted with the slope,
!
!    gamma_xz = eps_xx * m1 + eps_xy * m2 ,   gamma_yz = eps_xy * m1 + eps_yy * m2
!
! written in the sign convention of BS as it is returned above, and it is removed below.
!
! For a shell whose geometry is modelled exactly the director is the surface normal, m is zero and the term does not exist. It is
! produced purely by the mismatch between the flat facet and the nodal normals, so it is faceting error, not physics, and it grows
! linearly with the tilt. Subtracting it leaves the element free of membrane-shear coupling. Only the symmetric (strain) part of
! the in-plane displacement gradient is removed. The antisymmetric part is the in-plane rigid rotation, whose contribution is
! cancelled by the corresponding fibre rotation, so removing it as well would destroy rigid body invariance.

  CALL MITC_COVARIANT_BASIS( R, S, ZERO, G )

  IF (DABS(G(3,3)) > 1.0D-12 * DSQRT(DOT_PRODUCT(G(:,3), G(:,3)))) THEN

     M1 = G(1,3) / G(3,3)
     M2 = G(2,3) / G(3,3)

     BS(1,:) = BS(1,:) - ( BM(1,:) * M1 + HALF * BM(3,:) * M2 )
     BS(2,:) = BS(2,:) - ( HALF * BM(3,:) * M1 + BM(2,:) * M2 )

  ENDIF

  RETURN

! **********************************************************************************************************************************

CONTAINS

! **********************************************************************************************************************************

  SUBROUTINE TO_ELEMENT_BASIS ( T, B )

  ! Rotate one strain-displacement matrix from the cartesian local basis to the element coordinate system.

  REAL(DOUBLE), INTENT(IN)    :: T                 ! Isoparametric thickness coordinate the matrix was evaluated at
  REAL(DOUBLE), INTENT(INOUT) :: B(6, 6*ELGP)
  REAL(DOUBLE)                :: TRANSFORM(3,3)    ! Cartesian local basis vectors, in element coordinates

  TRANSFORM = MITC4_CARTESIAN_LOCAL_BASIS( R, S, T )

  B(4:6,:) = B(4:6,:) / 2                          ! Remove the engineering shear factor of 2 to rotate as a tensor
  CALL MITC_TRANSFORM_B( TRANSFORM, B )
  B(4:6,:) = B(4:6,:) * 2                          ! Reinstate it

  RETURN

  END SUBROUTINE TO_ELEMENT_BASIS

END SUBROUTINE MITC4_BMBS
