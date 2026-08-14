! ##################################################################################################################################
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

      SUBROUTINE RBE3_PROC ( RTYPE, REC_NO, IERR )

! Processes a single RBE3 "rigid" element, per call, to get terms for the RMG constraint matrix. When the Bulk data was read, the
! RBE3 input data was written to file LINK1F. In this subr, file LINK1F is read and RBE3 terms for array RMG are calculated and
! written to file LINK1J.  Later, in subr SPARSE_RMG, LINK1J will be read to create the sparse array RMG (of all rigid element and
! MPC coefficients) which will be used in LINK2 to reduce the G-set mass, stiffness and load matrices to the N-set.

! The derivation of the equations for the RBE3 are shown in Appendix E to the MYSTRAN User's Reference Manual

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, L1F, LINK1F, L1F_MSG, L1J
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MRBE3, NCORD, NGRID, NTERM_RMG
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE MODEL_STUF, ONLY            :  CORD, GRID_ID, GRID, RCORD, RGRID
      USE PARAMS, ONLY                :  EPSIL
      USE DOF_TABLES, ONLY            :  TDOF, TDOF_ROW_START

      USE RBE3_PROC_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'RBE3_PROC'
      CHARACTER( 8*BYTE), INTENT(IN)  :: RTYPE             ! The type of rigid element being processed (RBE2)
      CHARACTER( 1*BYTE)              :: CDOF_D(6)         ! An output from subr RDOF (= 1 if a displ comp is in COMPS_D)
      CHARACTER( 1*BYTE)              :: CDOF_I(6)         ! An output from subr RDOF (= 1 if a displ comp 1-6 is in COMPS_I)

      INTEGER(LONG), INTENT(INOUT)    :: IERR              ! Count of errors in RIGID_ELEM_PROC
      INTEGER(LONG), INTENT(INOUT)    :: REC_NO            ! Record number when reading file L1F
      INTEGER(LONG)                   :: AGRID_D           ! Dep   grid ID (actual) read from a record of file LINK1F
      INTEGER(LONG)                   :: AGRID_I(MRBE3)    ! Indep grid ID (actual) read from a record of file LINK1F
      INTEGER(LONG)                   :: COMPS_D           ! Dipsl components associated with dep   grid, AGRID_D
      INTEGER(LONG)                   :: COMPS_I(MRBE3)    ! Dipsl components associated with indep grid, AGRID_I
      INTEGER(LONG)                   :: ECORD_D           ! Global coord ID (actual) for grid AGRID_D
      INTEGER(LONG)                   :: ECORD_I           ! Global coord ID (actual) for grid AGRID_I
      INTEGER(LONG)                   :: GRID_ID_ROW_NUM_D ! Row number in array GRID_ID where AGRID_D is found
      INTEGER(LONG)                   :: GRID_ID_ROW_NUM_I ! Row number in array GRID_ID where AGRID_I is found
      INTEGER(LONG)                   :: G_SET_COL_NUM     ! Col no., in TDOF array, of the G-set DOF list
      INTEGER(LONG)                   :: I,J,K,L           ! DO loop indices
      INTEGER(LONG)                   :: ICORD_D           ! Internal coord ID corresponding to ECORD_D
      INTEGER(LONG)                   :: ICORD_I           ! Internal coord ID corresponding to ECORD_I
      INTEGER(LONG)                   :: IGRID             ! Internal grid ID
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening/reading a file
      INTEGER(LONG)                   :: IRBE3             ! Number of triplets of grid/comp/weight on L1F for RBE3 elems
      INTEGER(LONG)                   :: IROW              ! A row number in matrix RDI_GLOBAL
      INTEGER(LONG)                   :: JERR              ! Local error count
      INTEGER(LONG)                   :: ITERM_RMG         ! Countof number of records written to L1J (should be NTERM_RMG at end)
      INTEGER(LONG)                   :: M_SET_COL_NUM     ! Col no., in TDOF array, of the M-set DOF list
      INTEGER(LONG)                   :: NUM_COMPS         ! Number of displ components for a grid
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN
      INTEGER(LONG)                   :: REID              ! RBE2 elem ID read from file LINK1F
      INTEGER(LONG)                   :: RMG_COL_NUM_D(6)  ! Col no's. in RMG for 6 components of dep DOF at ref pt (if they exist)
      INTEGER(LONG)                   :: RMG_ROW_NUM       ! Row no. of a term in array RMG
      INTEGER(LONG)                   :: ROW_NUM           ! A row number in array TDOF
      INTEGER(LONG)                   :: ROW_NUM_START_D   ! DOF number where TDOF data begins for the ref grid


      REAL(DOUBLE)                    :: EPS1              ! Small number
      REAL(DOUBLE)                    :: DX_BAR            ! Wgt'd avg diff in x dist from indep pt i to ref pt A (in ref pt global)
      REAL(DOUBLE)                    :: DY_BAR            ! Wgt'd avg diff in y dist from indep pt i to ref pt A (in ref pt global)
      REAL(DOUBLE)                    :: DZ_BAR            ! Wgt'd avg diff in z dist from indep pt i to ref pt A (in ref pt global)
      REAL(DOUBLE)                    :: SX_DY_BAR, SX_DZ_BAR   ! X-weight applied to Y, Z offsets
      REAL(DOUBLE)                    :: SY_DX_BAR, SY_DZ_BAR   ! Y-weight applied to X, Z offsets
      REAL(DOUBLE)                    :: SZ_DX_BAR, SZ_DY_BAR   ! Z-weight applied to X, Y offsets
      REAL(DOUBLE)                    :: DX0(3)            ! Differences in coords of one indep pt and ref pt in basic coord system
      REAL(DOUBLE)                    :: DXI(MRBE3)        ! Differences in X coords of indep pt and ref pt in ref pt global system
      REAL(DOUBLE)                    :: DYI(MRBE3)        ! Differences in Y coords of indep pt and ref pt in ref pt global system
      REAL(DOUBLE)                    :: DZI(MRBE3)        ! Differences in Z coords of indep pt and ref pt in ref pt global system
      REAL(DOUBLE)                    :: PHID, THETAD      ! Angles output from subr GEN_T0L, called herein but not needed here
      REAL(DOUBLE)                    :: DUM3(3)           ! Intermediate result in a calc
      REAL(DOUBLE)                    :: EBAR_YZ           ! Sum of weights times radii squared divided by WT for rotation about x
      REAL(DOUBLE)                    :: EBAR_ZX           ! Sum of weights times radii squared divided by WT for rotation about y
      REAL(DOUBLE)                    :: EBAR_XY           ! Sum of weights times radii squared divided by WT for rotation about z
      REAL(DOUBLE)                    :: T0D(3,3)          ! Transform a vector to basic coords from one in global coords at AGRID_D
      REAL(DOUBLE)                    :: TDI(3,3)          ! TOD'*T0I
      REAL(DOUBLE)                    :: T0I(3,3)          ! Transform a vector to basic coords from one in global coords at AGRID_I
      REAL(DOUBLE)                    :: X0_D(3)           ! Basic coords of AGRID_D reference point
      REAL(DOUBLE)                    :: X0_I(3)           ! Basic coords of AGRID_D reference point
      REAL(DOUBLE)                    :: WTi(MRBE3)        ! Weight value for an indep grid
      REAL(DOUBLE)                    :: WT                ! Sum of weights on this RBE3
      REAL(DOUBLE)                    :: WT6(6)            ! WT6(i) = Sum of weights in comp i of an indep grid NB *** new 10/03/21

      REAL(DOUBLE)                    :: SXY,SZX,SYZ       ! new Rdd terms according to victor
      REAL(DOUBLE)                    :: WTi6(MRBE3,6)     ! per-DoF grid weights

! FIX (partial-REFC coupling bug): the full 6x6 coupled rigid-body least-squares system and the
! full coefficient table (one column per active independent-grid component) are now ALWAYS built,
! regardless of which components REFC selects as dependent. Whenever REFC excludes some
! components, they are eliminated via a Schur complement (not simply dropped), so the retained
! (REFC-selected) rows correctly account for their coupling to the excluded ones -- matching MSC
! Nastran's actual behavior (confirmed: MSC always solves the full 6-DOF system internally and
! REFC only controls what's output, not what's computed). When REFC=123456 (nothing excluded),
! this reduces identically to the original per-row computation -- verified by regression test.
      REAL(DOUBLE)                    :: A6(6,6)              ! full coupled system matrix
      REAL(DOUBLE)                    :: B6(6,6*MRBE3)        ! full coefficient table (row, indep grid*comp)
      INTEGER(LONG)                   :: B6_COL_RMG(6*MRBE3)  ! actual RMG column number for each B6 column (0=inactive)
      INTEGER(LONG)                   :: NB6COLS              ! number of columns actually used in B6/B6_COL_RMG
      LOGICAL                         :: IS_R(6)              ! TRUE if component is REFC-selected (retained)
      INTEGER(LONG)                   :: R_IDX(6), D_IDX(6)   ! lists of retained / discarded component indices
      INTEGER(LONG)                   :: NR, ND               ! counts of retained / discarded components
      REAL(DOUBLE)                    :: ADD(5,5), ADD_SAVE(5,5) ! discarded-discarded block (max Nd=5) and a working copy
      REAL(DOUBLE)                    :: RHS_MAT(5,6+6*MRBE3) ! combined RHS for the Gauss solve: [A_dr | B_d]
      REAL(DOUBLE)                    :: A_EFF(6,6)           ! Schur-reduced system matrix (only r,r entries meaningful)
      REAL(DOUBLE)                    :: B_EFF(6,6*MRBE3)     ! Schur-reduced coefficient table (only r rows meaningful)
      INTEGER(LONG)                   :: II, JJ, KK, PIVROW   ! loop/pivot indices for the Gauss elimination
      REAL(DOUBLE)                    :: PIVVAL, FACTOR, TMPSWAP


! **********************************************************************************************************************************
! File LINK1F contains data from the logical RBE3 cards in the input B.D. deck. For each logical RBE3 card, LINK1F has:
!     1st record          :'RBE3' the element type. This record was read in subr RIGID_ELEM_PROC before calling this subr

!     2nd record          : REID   : elem ID
!                           AGRID_D: reference (or dependent) grid
!                           COMPS_D: dependent displ comps
!                           IRBE3  : number of independent sets of grid/components/weight in the element
!                           WT     : total of all of the WTi weights on the RBE3 entry (calc'd when RBE3 bdf entry read in BD_RBE3

!     3rd record          : GRID(1), COMP(1), WTi(1): 1st independent grid, the independent components and the weight for this grid

!     4th record          : GRID(2), COMP(2), WTi(2): 2nd independent grid, the independent components and the weight for this grid

!     5th record, and on, : GRID(3), COMP(3), WTi(3): 3rd independent grid, the independent components and the weight for this grid

! The above record structure is repeated for each RBE3 logical card in the data deck (in the order in which they were read from the
! B.D. deck).

! Make units for writing errors the error file and output file

      OUNT(1) = ERR
      OUNT(2) = F06

      EPS1 = EPSIL(1)

      JERR = 0

! Init weight totals in each of the 6 components           ! NB *** new 10/03/21

      DO I=1,6                                             ! NB *** new 10/03/21
         WT6(I) = ZERO                                     ! NB *** new 10/03/21
      ENDDO                                                ! NB *** new 10/03/21

      ! zero-init WTi6
      WTi6 = ZERO

! Start reading at the 2nd record of L1F for this RBE3 (first record, RYPE, was read above in calling subr, RIGID_ELEM_PROC):
                                                           ! Read 2nd record from L1F for this RBE3
      READ(L1F,IOSTAT=IOCHK) REID, AGRID_D, COMPS_D, IRBE3, WT

      CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID_D, GRID_ID_ROW_NUM_D )

      REC_NO = REC_NO + 1
      IF (IOCHK == 0) THEN
         CALL GET_GRID_NUM_COMPS ( GRID_ID_ROW_NUM_D, NUM_COMPS, SUBR_NAME )
         IF (NUM_COMPS /= 6) THEN
            IERR  = IERR + 1
            JERR = JERR + 1
            WRITE(ERR,1951) 'RBE3', REID, NUM_COMPS
            WRITE(F06,1951) 'RBE3', REID, NUM_COMPS
         ENDIF
      ELSE
         CALL READERR ( IOCHK, LINK1F, L1F_MSG, REC_NO, OUNT )
         IERR = IERR + 1
         JERR = JERR + 1
      ENDIF

      DO I=1,IRBE3                                         ! Read remaining records from L1F for this RBE3
         READ(L1F,IOSTAT=IOCHK) AGRID_I(I), COMPS_I(I), WTi(I)
         REC_NO = REC_NO + 1
         IF (IOCHK /= 0) THEN
            CALL READERR ( IOCHK, LINK1F, L1F_MSG, REC_NO, OUNT )
            IERR = IERR + 1
            JERR = JERR + 1
         ENDIF
         CALL RDOF ( COMPS_I(I), CDOF_I )
         DO J=1,6
            IF (CDOF_I(J) == '1') THEN
               WTi6(I,J) = WTi(I)
               WT6(J) = WT6(J) + WTi(I)
            END IF
         END DO
      ENDDO

! Return if error

      IF (JERR /= 0) THEN
         FATAL_ERR = FATAL_ERR + 1
         RETURN
      ENDIF

! Get T0D (transforms global vector at AGRID_D to basic)

      ECORD_D = GRID(GRID_ID_ROW_NUM_D,3)
      IF (ECORD_D /= 0) THEN
         DO I=1,NCORD
            IF (ECORD_D == CORD(I,2)) THEN
               ICORD_D = I
               EXIT
            ENDIF
         ENDDO
         CALL GEN_T0L ( GRID_ID_ROW_NUM_D, ICORD_D, THETAD, PHID, T0D )
      ELSE
         DO I=1,3
            DO J=1,3
               T0D(I,J) = ZERO
            ENDDO
            T0D(I,I) = ONE
         ENDDO
      ENDIF

! Get coords of the reference grid (AGRID_D) in basic coord system

      DO I=1,3
         X0_D(I) = RGRID(GRID_ID_ROW_NUM_D,I)
      ENDDO

! Calc DXI, DYI, DZI, DX_BAR, DY_BAR, DZ_BAR

      DX_BAR = ZERO
      DY_BAR = ZERO
      DZ_BAR = ZERO
      SX_DY_BAR = ZERO;  SX_DZ_BAR = ZERO
      SY_DX_BAR = ZERO;  SY_DZ_BAR = ZERO
      SZ_DX_BAR = ZERO;  SZ_DY_BAR = ZERO

      DO J=1,IRBE3

         CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID_I(J), GRID_ID_ROW_NUM_I )
         DO K=1,3
            X0_I(K) = RGRID(GRID_ID_ROW_NUM_I,K)
            DX0(K)  = X0_I(K) - X0_D(K)
         ENDDO
                                                           ! Transform rel coords from basic to the coord sys at the ref pt
         CALL MATMULT_FFF_T ( T0D, DX0, 3, 3, 1, DUM3 )
                                                           ! Calc radius and in-plane angle from ref pt to indep points
         DXI(J) = DUM3(1)
         DYI(J) = DUM3(2)
         DZI(J) = DUM3(3)

         DX_BAR = DX_BAR + WTi6(J,1)*DXI(J)
         DY_BAR = DY_BAR + WTi6(J,2)*DYI(J)
         DZ_BAR = DZ_BAR + WTi6(J,3)*DZI(J)
         SX_DY_BAR = SX_DY_BAR + WTi6(J,1)*DYI(J)
         SX_DZ_BAR = SX_DZ_BAR + WTi6(J,1)*DZI(J)
         SY_DX_BAR = SY_DX_BAR + WTi6(J,2)*DXI(J)
         SY_DZ_BAR = SY_DZ_BAR + WTi6(J,2)*DZI(J)
         SZ_DX_BAR = SZ_DX_BAR + WTi6(J,3)*DXI(J)
         SZ_DY_BAR = SZ_DY_BAR + WTi6(J,3)*DYI(J)


      ENDDO


! Calc the EBAR's

      EBAR_YZ = ZERO
      EBAR_ZX = ZERO
      EBAR_XY = ZERO

      DO J=1,IRBE3
         EBAR_YZ = EBAR_YZ + WTi6(J,3)*DYI(J)*DYI(J) + WTi6(J,2)*DZI(J)*DZI(J)
         EBAR_ZX = EBAR_ZX + WTi6(J,1)*DZI(J)*DZI(J) + WTi6(J,3)*DXI(J)*DXI(J)
         EBAR_XY = EBAR_XY + WTi6(J,2)*DXI(J)*DXI(J) + WTi6(J,1)*DYI(J)*DYI(J)
         EBAR_YZ = EBAR_YZ + WTi6(J,4)
         EBAR_ZX = EBAR_ZX + WTi6(J,5)
         EBAR_XY = EBAR_XY + WTi6(J,6)
      ENDDO

! Calc the S-terms
      SXY = ZERO
      SZX = ZERO
      SYZ = ZERO

      DO J=1,IRBE3
         SXY = SXY + WTi6(J,3) * DXI(J) * DYI(J)
         SZX = SZX + WTi6(J,2) * DZI(J) * DXI(J)
         SYZ = SYZ + WTi6(J,1) * DYI(J) * DZI(J)
      END DO

      CALL TDOF_COL_NUM ( 'G ', G_SET_COL_NUM )
      CALL TDOF_COL_NUM ( 'M ', M_SET_COL_NUM )
      CALL RDOF ( COMPS_D, CDOF_D )
! Calc RMG_COL_NUM_D's no's for up to 6 DOF's for ref pt

!xx   CALL CALC_TDOF_ROW_NUM ( AGRID_D, ROW_NUM_START_D, 'N' )
      CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID_D, IGRID )
      ROW_NUM_START_D = TDOF_ROW_START(IGRID)
      DO I=1,6
         RMG_COL_NUM_D(I) = 0
         IF (CDOF_D(I) == '1') THEN
            IROW = I
            RMG_COL_NUM_D(I) = TDOF(ROW_NUM_START_D+I-1, G_SET_COL_NUM)
         ENDIF
      ENDDO

! Write terms to L1J for the constraint equations. The outer loop is for each of the RBE3 equations and the inner loop cycles over
! the IRBE3 grids in the "independent" set (the i points). There are up to 6 constraint eqns per RBE3 (1 each for T1, T2, T3,
! R1, R2 R3 comps in the indep set for the RBE3)

! FIX (partial-REFC coupling bug): build the FULL 6x6 coupled system A6 first, regardless of REFC.
! This is exactly the same math as before (WT6 diagonal for translation rows, EBAR for rotation
! rows, the S-terms and DX_BAR-family cross terms) -- just assembled into an explicit matrix
! instead of being written straight to RMG one REFC-selected row at a time.
      DO II=1,6
         DO JJ=1,6
            A6(II,JJ) = ZERO
         ENDDO
      ENDDO
      A6(1,1) = WT6(1);  A6(1,5) =  SX_DZ_BAR;  A6(1,6) = -SX_DY_BAR
      A6(2,2) = WT6(2);  A6(2,4) = -SY_DZ_BAR;  A6(2,6) =  SY_DX_BAR
      A6(3,3) = WT6(3);  A6(3,4) =  SZ_DY_BAR;  A6(3,5) = -SZ_DX_BAR
      A6(4,4) = EBAR_YZ; A6(4,5) = -SXY;        A6(4,6) = -SZX
      A6(5,5) = EBAR_ZX; A6(5,4) = -SXY;        A6(5,6) = -SYZ
      A6(6,6) = EBAR_XY; A6(6,4) = -SZX;        A6(6,5) = -SYZ
      A6(4,2) = -SY_DZ_BAR;  A6(4,3) =  SZ_DY_BAR
      A6(5,1) =  SX_DZ_BAR;  A6(5,3) = -SZ_DX_BAR
      A6(6,1) = -SX_DY_BAR;  A6(6,2) =  SY_DX_BAR
                                                           ! Mirror to guarantee exact symmetry
      DO II=1,6
         DO JJ=II+1,6
            A6(JJ,II) = A6(II,JJ)
         ENDDO
      ENDDO
! Build the full coefficient table B6: one column per (independent grid, active component),
! for ALL 6 rows -- again regardless of REFC. TDI is computed once per independent grid here
! (previously computed inside the row loop, redundantly, once per REFC-selected row).

      NB6COLS = 0
      DO JJ=1,6*MRBE3
         B6_COL_RMG(JJ) = 0
         DO II=1,6
            B6(II,JJ) = ZERO
         ENDDO
      ENDDO
      DO J=1,IRBE3
         CALL RDOF ( COMPS_I(J), CDOF_I )
         CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID_I(J), GRID_ID_ROW_NUM_I )
         CALL GET_GRID_NUM_COMPS ( GRID_ID_ROW_NUM_I, NUM_COMPS, SUBR_NAME )
         IF (NUM_COMPS /= 6) THEN
            IERR  = IERR + 1
            JERR  = JERR + 1
            WRITE(ERR,1951) 'RBE3', REID, NUM_COMPS
            WRITE(F06,1951) 'RBE3', REID, NUM_COMPS
            FATAL_ERR = FATAL_ERR + 1
            RETURN
         ENDIF
         ECORD_I= GRID(GRID_ID_ROW_NUM_I,3)
         IF (ECORD_I /= 0) THEN
            DO K=1,NCORD
               IF (ECORD_I == CORD(K,2)) THEN
                  ICORD_I = K
                  EXIT
               ENDIF
            ENDDO
            CALL GEN_T0L ( GRID_ID_ROW_NUM_I, ICORD_I, THETAD, PHID, T0I )
         ELSE
            DO K=1,3
               DO L=1,3
                  T0I(K,L) = ZERO
               ENDDO
               T0I(K,K) = ONE
            ENDDO
         ENDIF
         CALL MATMULT_FFF_T ( T0D, T0I, 3, 3, 3, TDI )
! Resolve the RMG column number for each of this grid's 6 components once, up front.
         CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID_I(J), IGRID )
         ROW_NUM_START_D = TDOF_ROW_START(IGRID)             ! reused as "ROW_NUM_START_I" here
         IF (TDOF(ROW_NUM_START_D,G_SET_COL_NUM) <= 0) THEN
            WRITE(ERR,'(A,I8,A)') ' *ERROR: RBE3_PROC found no valid G-set column for grid ', AGRID_I(J), ' (independent grid)'
            WRITE(F06,'(A,I8,A)') ' *ERROR: RBE3_PROC found no valid G-set column for grid ', AGRID_I(J), ' (independent grid)'
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ENDIF

         DO K=1,3
            NB6COLS = NB6COLS + 1
            IF (CDOF_I(K) == '1') THEN
               B6_COL_RMG(NB6COLS) = (TDOF(ROW_NUM_START_D,G_SET_COL_NUM)-1) + K
                                                              ! translation column K: rows 1-3 (WRITE_L1J_123-style)
               DO II=1,3
                  B6(II,NB6COLS) = -WTi6(J,K)*TDI(II,K)
               ENDDO
                                                              ! and rows 4-6 (WRITE_L1J_456-style translation terms)
               B6(4,NB6COLS) = B6(4,NB6COLS) + WTi6(J,K)*(DZI(J)*TDI(2,K) - DYI(J)*TDI(3,K))
               B6(5,NB6COLS) = B6(5,NB6COLS) + WTi6(J,K)*(-DZI(J)*TDI(1,K) + DXI(J)*TDI(3,K))
               B6(6,NB6COLS) = B6(6,NB6COLS) + WTi6(J,K)*(DYI(J)*TDI(1,K) - DXI(J)*TDI(2,K))
            ENDIF
         ENDDO

         DO K=4,6                                            ! rotation column K: gap-3 direct rotation coupling
            NB6COLS = NB6COLS + 1
            IF (CDOF_I(K) == '1') THEN
               B6_COL_RMG(NB6COLS) = (TDOF(ROW_NUM_START_D,G_SET_COL_NUM)-1) + K
               KK = K - 3
               B6(4,NB6COLS) = B6(4,NB6COLS) - WTi6(J,K)*TDI(1,KK)
               B6(5,NB6COLS) = B6(5,NB6COLS) - WTi6(J,K)*TDI(2,KK)
               B6(6,NB6COLS) = B6(6,NB6COLS) - WTi6(J,K)*TDI(3,KK)
            ENDIF
         ENDDO

      ENDDO

! Split the 6 components into "retained" (R, = REFC-selected/dependent) and "discarded" (D, not
! part of REFC) sets, then eliminate the D set from A6/B6 via a Schur complement, so the retained
! rows correctly account for their coupling to the discarded ones instead of simply ignoring it.

      NR = 0;  ND = 0
      DO II=1,6
         IS_R(II) = (CDOF_D(II) == '1')
         IF (IS_R(II)) THEN
            NR = NR + 1
            R_IDX(NR) = II
         ELSE
            ND = ND + 1
            D_IDX(ND) = II
         ENDIF
      ENDDO

      IF (ND == 0) THEN                                     ! Nothing to eliminate -- exact fast path, byte-identical
         DO II=1,6                                          ! to the original (pre-fix) computation.
            DO JJ=1,6
               A_EFF(II,JJ) = A6(II,JJ)
            ENDDO
         ENDDO
         DO II=1,6
            DO JJ=1,NB6COLS
               B_EFF(II,JJ) = B6(II,JJ)
            ENDDO
         ENDDO
      ELSE
                                                              ! Build A_dd and the combined RHS [A_dr | B_d]
         DO II=1,ND
            DO JJ=1,ND
               ADD(II,JJ) = A6(D_IDX(II),D_IDX(JJ))
            ENDDO
            DO JJ=1,NR
               RHS_MAT(II,JJ) = A6(D_IDX(II),R_IDX(JJ))
            ENDDO
            DO JJ=1,NB6COLS
               RHS_MAT(II,NR+JJ) = B6(D_IDX(II),JJ)
            ENDDO
         ENDDO
                                                              ! Gauss elimination with partial pivoting: solve
                                                              ! ADD * X = RHS_MAT for X (overwrite RHS_MAT with X)
         DO KK=1,ND
            PIVROW = KK
            PIVVAL = DABS(ADD(KK,KK))
            DO II=KK+1,ND
               IF (DABS(ADD(II,KK)) > PIVVAL) THEN
                  PIVROW = II
                  PIVVAL = DABS(ADD(II,KK))
               ENDIF
            ENDDO
            IF (PIVROW /= KK) THEN
               DO JJ=1,ND
                  TMPSWAP = ADD(KK,JJ); ADD(KK,JJ) = ADD(PIVROW,JJ); ADD(PIVROW,JJ) = TMPSWAP
               ENDDO
               DO JJ=1,NR+NB6COLS
                  TMPSWAP = RHS_MAT(KK,JJ); RHS_MAT(KK,JJ) = RHS_MAT(PIVROW,JJ); RHS_MAT(PIVROW,JJ) = TMPSWAP
               ENDDO
            ENDIF
            IF (DABS(ADD(KK,KK)) > EPS1) THEN
               DO II=KK+1,ND
                  FACTOR = ADD(II,KK)/ADD(KK,KK)
                  DO JJ=KK,ND
                     ADD(II,JJ) = ADD(II,JJ) - FACTOR*ADD(KK,JJ)
                  ENDDO
                  DO JJ=1,NR+NB6COLS
                     RHS_MAT(II,JJ) = RHS_MAT(II,JJ) - FACTOR*RHS_MAT(KK,JJ)
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
                                                              ! Back-substitution
         DO KK=ND,1,-1
            IF (DABS(ADD(KK,KK)) > EPS1) THEN
               DO JJ=1,NR+NB6COLS
                  DO II=KK+1,ND
                     RHS_MAT(KK,JJ) = RHS_MAT(KK,JJ) - ADD(KK,II)*RHS_MAT(II,JJ)
                  ENDDO
                  RHS_MAT(KK,JJ) = RHS_MAT(KK,JJ)/ADD(KK,KK)
               ENDDO
            ELSE                                             ! degenerate discarded block: no correction from this DOF
               DO JJ=1,NR+NB6COLS
                  RHS_MAT(KK,JJ) = ZERO
               ENDDO
            ENDIF
         ENDDO
                                                              ! A_eff = A_rr - A_rd*X ;  B_eff = B_r - A_rd*Y
         DO II=1,NR
            DO JJ=1,NR
               A_EFF(II,JJ) = A6(R_IDX(II),R_IDX(JJ))
               DO KK=1,ND
                  A_EFF(II,JJ) = A_EFF(II,JJ) - A6(R_IDX(II),D_IDX(KK))*RHS_MAT(KK,JJ)
               ENDDO
            ENDDO
            DO JJ=1,NB6COLS
               B_EFF(II,JJ) = B6(R_IDX(II),JJ)
               DO KK=1,ND
                  B_EFF(II,JJ) = B_EFF(II,JJ) - A6(R_IDX(II),D_IDX(KK))*RHS_MAT(KK,NR+JJ)
               ENDDO
            ENDDO
         ENDDO
                                                              ! Re-expand A_EFF/B_EFF back to full 1..6 row indexing
                                                              ! (rows for indices in R_IDX only; others unused/ignored)
         DO II=NR,1,-1
            DO JJ=1,NR
               A_EFF(R_IDX(II),R_IDX(JJ)) = A_EFF(II,JJ)
            ENDDO
            DO JJ=1,NB6COLS
               B_EFF(R_IDX(II),JJ) = B_EFF(II,JJ)
            ENDDO
         ENDDO
      ENDIF

! Write terms to L1J for the constraint equations, now using the (possibly Schur-reduced) A_EFF
! and B_EFF instead of the raw pivots/S-terms/DX_BAR-family sums and per-row WRITE_L1J_123/456
! calls. There are up to 6 constraint eqns per RBE3 (1 each for T1, T2, T3, R1, R2, R3 comps).

      ITERM_RMG = 0
do_i1:DO I=1,6
cdof_dep:IF (CDOF_D(I) == '1') THEN                        ! The I-th component is in DDOF so write this row to RMG
            IROW = I
            CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID_D, IGRID )
            ROW_NUM_START_D = TDOF_ROW_START(IGRID)
            ROW_NUM = ROW_NUM_START_D + I - 1
            RMG_ROW_NUM = TDOF(ROW_NUM, M_SET_COL_NUM)

            IF ((RMG_ROW_NUM > 0) .AND. (RMG_COL_NUM_D(I) > 0)) THEN

               IF ((I == 1) .OR. (I == 2) .OR. (I == 3)) THEN
                  IF (DABS(WT) <= EPS1) THEN
                     WRITE(L1J) RMG_ROW_NUM, RMG_COL_NUM_D(I), A_EFF(I,I)
                     ITERM_RMG = ITERM_RMG + 1
                     CYCLE do_i1
                  ENDIF
               ENDIF
                                                           ! Pivot term (own column), with the same near-zero fallback
                                                           ! to a unit pivot the original code used for rows 4-6
               IF ((I == 4) .OR. (I == 5) .OR. (I == 6)) THEN
                  IF (DABS(A_EFF(I,I)) > EPS1) THEN
                     WRITE(L1J) RMG_ROW_NUM, RMG_COL_NUM_D(I), A_EFF(I,I)
                  ELSE
                     WRITE(L1J) RMG_ROW_NUM, RMG_COL_NUM_D(I), ONE
                  ENDIF
               ELSE
                  WRITE(L1J) RMG_ROW_NUM, RMG_COL_NUM_D(I), A_EFF(I,I)
               ENDIF
               ITERM_RMG = ITERM_RMG + 1
                                                           ! Cross terms to the OTHER retained (REFC-selected) comps
               DO JJ=1,6
                  IF ((JJ /= I) .AND. (CDOF_D(JJ) == '1')) THEN
                     IF (A_EFF(I,JJ) /= ZERO) THEN
                        WRITE(L1J) RMG_ROW_NUM, RMG_COL_NUM_D(JJ), A_EFF(I,JJ)
                        ITERM_RMG = ITERM_RMG + 1
                     ENDIF
                  ENDIF
               ENDDO
                                                           ! Independent-grid terms
               DO JJ=1,NB6COLS
                  IF ((B6_COL_RMG(JJ) > 0) .AND. (B_EFF(I,JJ) /= ZERO)) THEN
                     WRITE(L1J) RMG_ROW_NUM, B6_COL_RMG(JJ), B_EFF(I,JJ)
                     ITERM_RMG = ITERM_RMG + 1
                  ENDIF
               ENDDO

            ELSE
               IF (RMG_ROW_NUM  == 0) THEN
                  WRITE(ERR,1509) SUBR_NAME,RTYPE,REID,AGRID_D,IROW
                  WRITE(F06,1509) SUBR_NAME,RTYPE,REID,AGRID_D,IROW
                  FATAL_ERR = FATAL_ERR + 1
                  IERR = IERR + 1
                  JERR = JERR + 1
               ENDIF
               IF (RMG_COL_NUM_D(I) == 0) THEN
                  WRITE(ERR,1510) SUBR_NAME,RTYPE,REID,AGRID_D,IROW
                  WRITE(F06,1510) SUBR_NAME,RTYPE,REID,AGRID_D,IROW
                  FATAL_ERR = FATAL_ERR + 1
                  IERR = IERR + 1
                  JERR = JERR + 1
               ENDIF
            ENDIF

         ENDIF cdof_dep

      ENDDO do_i1

      NTERM_RMG = NTERM_RMG + ITERM_RMG

! Return if JERR > 0

      IF (JERR > 0) THEN
         RETURN
      ENDIF



      RETURN

! **********************************************************************************************************************************

 1509 FORMAT(' *ERROR  1509: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,15X,A8,' RIGID ELEMENT NUMBER ',I8,', DEPENDENT GRID NUMBER ',I8,', COMPONENT ',I2                          &
                    ,/,14X,' IS NOT A M-SET DOF IN TABLE TDOFI')

 1510 FORMAT(' *ERROR  1510: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,15X,A8,' RIGID ELEMENT NUMBER ',I8,', INDEPENDENT GRID NUMBER ',I8,', COMPONENT ',I2                        &
                    ,/,14X,' IS NOT A G-SET DOF IN TABLE TDOFI')

 1951 FORMAT(' *ERROR  1951: ',A,I8,' USES GRID ',I8,' WHICH IS A SCALAR POINT. SCALAR POINTS NOT ALLOWED FOR THIS ELEM TYPE')


! **********************************************************************************************************************************


      END SUBROUTINE RBE3_PROC
