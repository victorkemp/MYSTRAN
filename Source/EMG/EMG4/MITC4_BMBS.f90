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
  USE CONSTANTS_1, ONLY       : ZERO, ONE

  USE MITC4_B_Interface
  USE MITC4_CARTESIAN_LOCAL_BASIS_Interface
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

! **********************************************************************************************************************************
! Subr MITC4_B returns the strain-displacement matrix in the CARTESIAN LOCAL basis, whose x axis lies along the covariant g_r
! direction. Everything that uses the operators returned here - the elasticity matrices built from EM, EB and ET, the stress,
! strain and engineering force output, and the extrapolation of Gauss point values to the corners - works in the ELEMENT
! coordinate system. So the rows have to be rotated out of the cartesian local basis before they are handed back.
!
! For a rectangular element the two systems differ by 180 degrees about z, because the element x axis starts along side 1-2
! while the cartesian local x axis lies along g_r, and side 1-2 runs in the negative r direction in Bathe's node ordering. A
! 180 degree rotation about z leaves the in-plane rows xx, yy and xy untouched and negates both transverse shear rows, which is
! why omitting this step produced correctly signed membrane and bending output alongside transverse shear strain, stress and
! engineering force of the wrong sign. For a skewed or warped element the angle is not 180 degrees and the in-plane rows are
! wrong as well, so this is a rotation and not a sign correction.
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
