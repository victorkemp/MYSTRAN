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

! Form the RMG constraint equations for one RBE3 element. Each active scalar independent DOF contributes one residual to a weighted
! rigid-body least-squares fit. If H is the kinematic vector for that DOF, its contribution is
!
!                  A6 = A6 + WEIGHT*H*TRANSPOSE(H)       and       B6(:,COL) = -WEIGHT*H .
!
! The six reference-grid components satisfy A6*Q + B6*U = 0. Components omitted by REFC are eliminated with a rank-revealing LAPACK
! solve before the retained equations are written to LINK1J. See Appendix E of the MYSTRAN User's Reference Manual.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, L1F, LINK1F, L1F_MSG, L1J
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MRBE3, NCORD, NGRID, NTERM_RMG
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE MODEL_STUF, ONLY            :  CORD, GRID_ID, GRID, RGRID
      USE PARAMS, ONLY                :  EPSIL
      USE DOF_TABLES, ONLY            :  TDOF, TDOF_ROW_START

      USE RBE3_PROC_USE_IFs

      IMPLICIT NONE

      INTEGER(LONG), PARAMETER        :: NUM_RIGID_DOF = 6
      INTEGER(LONG), PARAMETER        :: NUM_VECTOR_DOF = 3
      INTEGER(LONG), PARAMETER        :: ALL_DOF(NUM_RIGID_DOF) = [1_LONG, 2_LONG, 3_LONG, 4_LONG, 5_LONG, 6_LONG]

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'RBE3_PROC'
      CHARACTER( 8*BYTE), INTENT(IN)  :: RTYPE

      INTEGER(LONG), INTENT(INOUT)    :: IERR
      INTEGER(LONG), INTENT(INOUT)    :: REC_NO

      CHARACTER( 1*BYTE)              :: INDEP_DOF_ACTIVE(NUM_RIGID_DOF)
      CHARACTER( 1*BYTE)              :: REF_DOF_ACTIVE(NUM_RIGID_DOF)

      INTEGER(LONG)                   :: AGRID_D
      INTEGER(LONG)                   :: AGRID_I(MRBE3)
      INTEGER(LONG)                   :: B6_COL_RMG(NUM_RIGID_DOF*MRBE3)
      INTEGER(LONG)                   :: COMPS_D
      INTEGER(LONG)                   :: COMPS_I(MRBE3)
      INTEGER(LONG)                   :: DISCARDED_IDX(NUM_RIGID_DOF)
      INTEGER(LONG)                   :: G_SET_COL_NUM
      INTEGER(LONG)                   :: GRID_ID_ROW_NUM_D
      INTEGER(LONG)                   :: IRBE3
      INTEGER(LONG)                   :: ITERM_RMG
      INTEGER(LONG)                   :: JERR
      INTEGER(LONG)                   :: M_SET_COL_NUM
      INTEGER(LONG)                   :: NB6COLS
      INTEGER(LONG)                   :: NUM_DISCARDED
      INTEGER(LONG)                   :: NUM_RETAINED
      INTEGER(LONG)                   :: REID
      INTEGER(LONG)                   :: RETAINED_IDX(NUM_RIGID_DOF)

      REAL(DOUBLE)                    :: A6(NUM_RIGID_DOF,NUM_RIGID_DOF)
      REAL(DOUBLE)                    :: A_REDUCED(NUM_RIGID_DOF,NUM_RIGID_DOF)
      REAL(DOUBLE)                    :: B6(NUM_RIGID_DOF,NUM_RIGID_DOF*MRBE3)
      REAL(DOUBLE)                    :: B_REDUCED(NUM_RIGID_DOF,NUM_RIGID_DOF*MRBE3)
      REAL(DOUBLE)                    :: EPS1
      REAL(DOUBLE)                    :: REFERENCE_POSITION(NUM_VECTOR_DOF)
      REAL(DOUBLE)                    :: T0D(NUM_VECTOR_DOF,NUM_VECTOR_DOF)
      REAL(DOUBLE)                    :: WEIGHT(MRBE3)
      REAL(DOUBLE)                    :: WEIGHT_SUM_ON_FILE

      INTERFACE
         SUBROUTINE DGELSY ( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, WORK, LWORK, INFO )
            IMPORT LONG, DOUBLE
            INTEGER(LONG), INTENT(IN)    :: M, N, NRHS, LDA, LDB, LWORK
            INTEGER(LONG), INTENT(INOUT) :: JPVT(*)
            INTEGER(LONG), INTENT(OUT)   :: RANK, INFO
            REAL(DOUBLE), INTENT(INOUT)  :: A(LDA,*), B(LDB,*), WORK(*)
            REAL(DOUBLE), INTENT(IN)     :: RCOND
         END SUBROUTINE DGELSY
      END INTERFACE

! **********************************************************************************************************************************

      EPS1 = EPSIL(1)
      JERR = 0

      CALL READ_RBE3_INPUT
      IF (JERR /= 0) THEN
         FATAL_ERR = FATAL_ERR + 1
         RETURN
      ENDIF

      CALL TDOF_COL_NUM ( 'G ', G_SET_COL_NUM )
      CALL TDOF_COL_NUM ( 'M ', M_SET_COL_NUM )

      CALL GET_GRID_TRANSFORM ( GRID_ID_ROW_NUM_D, T0D )
      IF (JERR /= 0) THEN
         FATAL_ERR = FATAL_ERR + 1
         RETURN
      ENDIF
      REFERENCE_POSITION = RGRID(GRID_ID_ROW_NUM_D,1:NUM_VECTOR_DOF)

      CALL ASSEMBLE_LEAST_SQUARES_SYSTEM
      IF (JERR /= 0) THEN
         FATAL_ERR = FATAL_ERR + 1
         RETURN
      ENDIF

      CALL SELECT_REFERENCE_COMPONENTS
      CALL REDUCE_TO_REFERENCE_COMPONENTS
      IF (JERR /= 0) THEN
         FATAL_ERR = FATAL_ERR + 1
         RETURN
      ENDIF

      CALL WRITE_RMG_TERMS ( ITERM_RMG )
      NTERM_RMG = NTERM_RMG + ITERM_RMG

      IF (JERR /= 0) FATAL_ERR = FATAL_ERR + 1

      RETURN

! **********************************************************************************************************************************

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE READ_RBE3_INPUT

      INTEGER(LONG)                   :: I
      INTEGER(LONG)                   :: IOCHK
      INTEGER(LONG)                   :: NUM_COMPS
      INTEGER(LONG)                   :: OUNT(2)
      INTEGER(LONG)                   :: REID_LOCAL

      OUNT = [ERR, F06]
      AGRID_I = 0
      COMPS_I = 0
      WEIGHT = ZERO

! The element type record was read by RIGID_ELEM_PROC. Read the header and only use its values after IOSTAT has been checked.

      READ(L1F,IOSTAT=IOCHK) REID_LOCAL, AGRID_D, COMPS_D, IRBE3, WEIGHT_SUM_ON_FILE
      REC_NO = REC_NO + 1
      IF (IOCHK /= 0) THEN
         CALL READERR ( IOCHK, LINK1F, L1F_MSG, REC_NO, OUNT )
         IERR = IERR + 1
         JERR = JERR + 1
         RETURN
      ENDIF
      REID = REID_LOCAL

      IF ((IRBE3 < 1) .OR. (IRBE3 > MRBE3)) THEN
         WRITE(ERR,1952) REID, IRBE3, MRBE3
         WRITE(F06,1952) REID, IRBE3, MRBE3
         IERR = IERR + 1
         JERR = JERR + 1
         RETURN
      ENDIF

      CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID_D, GRID_ID_ROW_NUM_D )
      CALL GET_GRID_NUM_COMPS ( GRID_ID_ROW_NUM_D, NUM_COMPS, SUBR_NAME )
      IF (NUM_COMPS /= NUM_RIGID_DOF) THEN
         WRITE(ERR,1951) 'RBE3', REID, AGRID_D
         WRITE(F06,1951) 'RBE3', REID, AGRID_D
         IERR = IERR + 1
         JERR = JERR + 1
      ENDIF

      DO I = 1, IRBE3
         READ(L1F,IOSTAT=IOCHK) AGRID_I(I), COMPS_I(I), WEIGHT(I)
         REC_NO = REC_NO + 1
         IF (IOCHK /= 0) THEN
            CALL READERR ( IOCHK, LINK1F, L1F_MSG, REC_NO, OUNT )
            IERR = IERR + 1
            JERR = JERR + 1
            RETURN
         ENDIF
      ENDDO

 1951 FORMAT(' *ERROR  1951: ',A,I8,' USES GRID ',I8,' WHICH IS A SCALAR POINT. SCALAR POINTS NOT ALLOWED FOR THIS ELEM TYPE')
 1952 FORMAT(' *ERROR  1952: RBE3 ',I8,' HAS ',I8,' INDEPENDENT GRID RECORDS; THE SUPPORTED RANGE IS 1 THROUGH ',I8)

      END SUBROUTINE READ_RBE3_INPUT

! ##################################################################################################################################

      SUBROUTINE GET_GRID_TRANSFORM ( GRID_ROW, TRANSFORM )

      INTEGER(LONG), INTENT(IN)      :: GRID_ROW
      REAL(DOUBLE), INTENT(OUT)      :: TRANSFORM(NUM_VECTOR_DOF,NUM_VECTOR_DOF)

      INTEGER(LONG)                  :: COORD_ID
      INTEGER(LONG)                  :: COORD_ROW
      REAL(DOUBLE)                   :: PHI
      REAL(DOUBLE)                   :: THETA

      COORD_ID = GRID(GRID_ROW,3)
      IF (COORD_ID == 0) THEN
         TRANSFORM = IDENTITY_3()
         RETURN
      ENDIF

      COORD_ROW = FINDLOC(CORD(1:NCORD,2), COORD_ID, DIM=1)
      IF (COORD_ROW == 0) THEN
         WRITE(ERR,1953) REID, GRID_ID(GRID_ROW), COORD_ID
         WRITE(F06,1953) REID, GRID_ID(GRID_ROW), COORD_ID
         IERR = IERR + 1
         JERR = JERR + 1
         TRANSFORM = ZERO
         RETURN
      ENDIF

      CALL GEN_T0L ( GRID_ROW, COORD_ROW, THETA, PHI, TRANSFORM )

 1953 FORMAT(' *ERROR  1953: RBE3 ',I8,' USES GRID ',I8,' WITH UNRESOLVED DISPLACEMENT COORDINATE SYSTEM ',I8)

      END SUBROUTINE GET_GRID_TRANSFORM

! ##################################################################################################################################

      SUBROUTINE ASSEMBLE_LEAST_SQUARES_SYSTEM

      INTEGER(LONG)                  :: COMPONENT
      INTEGER(LONG)                  :: GRID_INDEX
      INTEGER(LONG)                  :: GRID_ROW
      INTEGER(LONG)                  :: NUM_COMPS
      INTEGER(LONG)                  :: RMG_COLUMN
      INTEGER(LONG)                  :: ROW_START
      REAL(DOUBLE)                   :: H(NUM_RIGID_DOF)
      REAL(DOUBLE)                   :: RELATIVE_POSITION(NUM_VECTOR_DOF)
      REAL(DOUBLE)                   :: T0I(NUM_VECTOR_DOF,NUM_VECTOR_DOF)
      REAL(DOUBLE)                   :: TDI(NUM_VECTOR_DOF,NUM_VECTOR_DOF)

      A6 = ZERO
      B6 = ZERO
      B6_COL_RMG = 0
      NB6COLS = 0

      DO GRID_INDEX = 1, IRBE3
         CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID_I(GRID_INDEX), GRID_ROW )
         CALL GET_GRID_NUM_COMPS ( GRID_ROW, NUM_COMPS, SUBR_NAME )
         IF (NUM_COMPS /= NUM_RIGID_DOF) THEN
            WRITE(ERR,1951) 'RBE3', REID, AGRID_I(GRID_INDEX)
            WRITE(F06,1951) 'RBE3', REID, AGRID_I(GRID_INDEX)
            IERR = IERR + 1
            JERR = JERR + 1
            CYCLE
         ENDIF

         CALL GET_GRID_TRANSFORM ( GRID_ROW, T0I )
         IF (JERR /= 0) CYCLE

         RELATIVE_POSITION = MATMUL(TRANSPOSE(T0D), RGRID(GRID_ROW,1:NUM_VECTOR_DOF) - REFERENCE_POSITION)
         TDI = MATMUL(TRANSPOSE(T0D), T0I)
         ROW_START = TDOF_ROW_START(GRID_ROW)

         CALL RDOF ( COMPS_I(GRID_INDEX), INDEP_DOF_ACTIVE )
         DO COMPONENT = 1, NUM_RIGID_DOF
            IF (INDEP_DOF_ACTIVE(COMPONENT) /= '1') CYCLE

            RMG_COLUMN = TDOF(ROW_START + COMPONENT - 1,G_SET_COL_NUM)
            IF (RMG_COLUMN <= 0) THEN
               WRITE(ERR,1513) SUBR_NAME, AGRID_I(GRID_INDEX), COMPONENT, RMG_COLUMN
               WRITE(F06,1513) SUBR_NAME, AGRID_I(GRID_INDEX), COMPONENT, RMG_COLUMN
               IERR = IERR + 1
               JERR = JERR + 1
               CYCLE
            ENDIF

            H = KINEMATIC_VECTOR(COMPONENT, TDI, RELATIVE_POSITION)

            NB6COLS = NB6COLS + 1
            B6_COL_RMG(NB6COLS) = RMG_COLUMN
            B6(:,NB6COLS) = -WEIGHT(GRID_INDEX)*H
            A6 = A6 + WEIGHT(GRID_INDEX)*OUTER_PRODUCT_6(H)
         ENDDO
      ENDDO

 1513 FORMAT(' *ERROR  1513: PROGRAMMING ERROR IN SUBROUTINE ',A,/,14X,'G-SET COLUMN FOR GRID ',I8,', COMPONENT ',I2,              &
                    ' MUST BE > 0 BUT IS ',I8)
 1951 FORMAT(' *ERROR  1951: ',A,I8,' USES GRID ',I8,' WHICH IS A SCALAR POINT. SCALAR POINTS NOT ALLOWED FOR THIS ELEM TYPE')

      END SUBROUTINE ASSEMBLE_LEAST_SQUARES_SYSTEM

! ##################################################################################################################################

      SUBROUTINE SELECT_REFERENCE_COMPONENTS

      LOGICAL                        :: IS_RETAINED(NUM_RIGID_DOF)

      CALL RDOF ( COMPS_D, REF_DOF_ACTIVE )
      IS_RETAINED = (REF_DOF_ACTIVE == '1')
      NUM_RETAINED = COUNT(IS_RETAINED)
      NUM_DISCARDED = NUM_RIGID_DOF - NUM_RETAINED

      RETAINED_IDX = 0
      DISCARDED_IDX = 0
      RETAINED_IDX(1:NUM_RETAINED) = PACK(ALL_DOF, IS_RETAINED)
      IF (NUM_DISCARDED > 0) DISCARDED_IDX(1:NUM_DISCARDED) = PACK(ALL_DOF, .NOT. IS_RETAINED)

      END SUBROUTINE SELECT_REFERENCE_COMPONENTS

! ##################################################################################################################################

      SUBROUTINE REDUCE_TO_REFERENCE_COMPONENTS

      INTEGER(LONG)                  :: INFO
      INTEGER(LONG)                  :: NUM_RIGHT_HAND_SIDES
      INTEGER(LONG)                  :: RANK
      REAL(DOUBLE)                   :: A_DD(5,5)
      REAL(DOUBLE)                   :: A_DD_ORIGINAL(5,5)
      REAL(DOUBLE)                   :: A_RD(NUM_RIGID_DOF,5)
      REAL(DOUBLE)                   :: RANK_CHECK_A(NUM_RIGID_DOF,NUM_RIGID_DOF)
      REAL(DOUBLE)                   :: RANK_CHECK_B(NUM_RIGID_DOF,1)
      REAL(DOUBLE)                   :: RESIDUAL_SCALE
      REAL(DOUBLE)                   :: RHS(5,NUM_RIGID_DOF+NUM_RIGID_DOF*MRBE3)
      REAL(DOUBLE)                   :: RHS_ORIGINAL(5,NUM_RIGID_DOF+NUM_RIGID_DOF*MRBE3)

      A_REDUCED = ZERO
      B_REDUCED = ZERO

      IF (NUM_RETAINED == 0) THEN
         WRITE(ERR,1954) REID, COMPS_D
         WRITE(F06,1954) REID, COMPS_D
         IERR = IERR + 1
         JERR = JERR + 1
         RETURN
      ENDIF

      IF (NUM_DISCARDED == 0) THEN
         A_REDUCED(1:NUM_RETAINED,1:NUM_RETAINED) = A6(RETAINED_IDX(1:NUM_RETAINED),RETAINED_IDX(1:NUM_RETAINED))
         B_REDUCED(1:NUM_RETAINED,1:NB6COLS) = B6(RETAINED_IDX(1:NUM_RETAINED),1:NB6COLS)
      ELSE
         NUM_RIGHT_HAND_SIDES = NUM_RETAINED + NB6COLS
         A_DD = ZERO
         RHS = ZERO
         A_RD = ZERO

         A_DD(1:NUM_DISCARDED,1:NUM_DISCARDED) = &
            A6(DISCARDED_IDX(1:NUM_DISCARDED),DISCARDED_IDX(1:NUM_DISCARDED))
         A_RD(1:NUM_RETAINED,1:NUM_DISCARDED) = &
            A6(RETAINED_IDX(1:NUM_RETAINED),DISCARDED_IDX(1:NUM_DISCARDED))
         RHS(1:NUM_DISCARDED,1:NUM_RETAINED) = &
            A6(DISCARDED_IDX(1:NUM_DISCARDED),RETAINED_IDX(1:NUM_RETAINED))
         RHS(1:NUM_DISCARDED,NUM_RETAINED+1:NUM_RIGHT_HAND_SIDES) = &
            B6(DISCARDED_IDX(1:NUM_DISCARDED),1:NB6COLS)

         A_DD_ORIGINAL = A_DD
         RHS_ORIGINAL = RHS
         CALL SOLVE_MINIMUM_NORM ( A_DD, RHS, NUM_DISCARDED, NUM_RIGHT_HAND_SIDES, RANK, INFO )
         IF (INFO /= 0) THEN
            WRITE(ERR,1955) REID, INFO
            WRITE(F06,1955) REID, INFO
            IERR = IERR + 1
            JERR = JERR + 1
            RETURN
         ENDIF

! A rank-deficient discarded block is permitted only when all elimination right-hand sides are in its range. DGELSY then supplies
! the unique minimum-norm solution needed for the generalized Schur complement.

         RESIDUAL_SCALE = MAX(ONE, MAXVAL(ABS(RHS_ORIGINAL(1:NUM_DISCARDED,1:NUM_RIGHT_HAND_SIDES))))
         IF (MAXIMUM_RESIDUAL(A_DD_ORIGINAL(1:NUM_DISCARDED,1:NUM_DISCARDED), &
                              RHS(1:NUM_DISCARDED,1:NUM_RIGHT_HAND_SIDES), &
                              RHS_ORIGINAL(1:NUM_DISCARDED,1:NUM_RIGHT_HAND_SIDES)) > SQRT(EPS1)*RESIDUAL_SCALE) THEN
            WRITE(ERR,1956) REID, RANK, NUM_DISCARDED
            WRITE(F06,1956) REID, RANK, NUM_DISCARDED
            IERR = IERR + 1
            JERR = JERR + 1
            RETURN
         ENDIF

         A_REDUCED(1:NUM_RETAINED,1:NUM_RETAINED) = SCHUR_REDUCTION( &
            A6(RETAINED_IDX(1:NUM_RETAINED),RETAINED_IDX(1:NUM_RETAINED)), &
            A_RD(1:NUM_RETAINED,1:NUM_DISCARDED), RHS(1:NUM_DISCARDED,1:NUM_RETAINED))
         B_REDUCED(1:NUM_RETAINED,1:NB6COLS) = SCHUR_REDUCTION( &
            B6(RETAINED_IDX(1:NUM_RETAINED),1:NB6COLS), A_RD(1:NUM_RETAINED,1:NUM_DISCARDED), &
            RHS(1:NUM_DISCARDED,NUM_RETAINED+1:NUM_RIGHT_HAND_SIDES))
      ENDIF

! A singular retained block cannot define every requested dependent DOF. Diagnose it here instead of inventing unit pivot terms.

      RANK_CHECK_A = ZERO
      RANK_CHECK_B = ZERO
      RANK_CHECK_A(1:NUM_RETAINED,1:NUM_RETAINED) = A_REDUCED(1:NUM_RETAINED,1:NUM_RETAINED)
      CALL SOLVE_MINIMUM_NORM ( RANK_CHECK_A, RANK_CHECK_B, NUM_RETAINED, 1_LONG, RANK, INFO )
      IF ((INFO /= 0) .OR. (RANK < NUM_RETAINED)) THEN
         WRITE(ERR,1957) REID, NUM_RETAINED, RANK
         WRITE(F06,1957) REID, NUM_RETAINED, RANK
         IERR = IERR + 1
         JERR = JERR + 1
      ENDIF

 1954 FORMAT(' *ERROR  1954: RBE3 ',I8,' HAS NO VALID REFC COMPONENTS IN VALUE ',I8)
 1955 FORMAT(' *ERROR  1955: LAPACK DGELSY FAILED WHILE REDUCING RBE3 ',I8,'; INFO = ',I8)
 1956 FORMAT(' *ERROR  1956: RBE3 ',I8,' HAS AN INCONSISTENT RANK-DEFICIENT DISCARDED-COMPONENT SYSTEM (RANK ',I2,' OF ',I2,')')
 1957 FORMAT(' *ERROR  1957: RBE3 ',I8,' CANNOT DETERMINE ALL ',I2,' REQUESTED REFC COMPONENTS; REDUCED RANK IS ',I2)

      END SUBROUTINE REDUCE_TO_REFERENCE_COMPONENTS

! ##################################################################################################################################

      SUBROUTINE SOLVE_MINIMUM_NORM ( MATRIX, RHS, N, NRHS, RANK, INFO )

      REAL(DOUBLE), INTENT(INOUT)    :: MATRIX(:,:)
      REAL(DOUBLE), INTENT(INOUT)    :: RHS(:,:)
      INTEGER(LONG), INTENT(IN)      :: N
      INTEGER(LONG), INTENT(IN)      :: NRHS
      INTEGER(LONG), INTENT(OUT)     :: RANK
      INTEGER(LONG), INTENT(OUT)     :: INFO

      INTEGER(LONG)                  :: JPVT(NUM_RIGID_DOF)
      INTEGER(LONG)                  :: LWORK
      REAL(DOUBLE), ALLOCATABLE      :: WORK(:)
      REAL(DOUBLE)                   :: WORK_QUERY(1)

      JPVT = 0
      CALL DGELSY ( N, N, NRHS, MATRIX, SIZE(MATRIX,1,KIND=LONG), RHS, SIZE(RHS,1,KIND=LONG), JPVT, EPS1, RANK, &
                    WORK_QUERY, -1_LONG, INFO )
      IF (INFO /= 0) RETURN

      LWORK = MAX(1_LONG, INT(WORK_QUERY(1),KIND=LONG))
      ALLOCATE(WORK(LWORK))
      JPVT = 0
      CALL DGELSY ( N, N, NRHS, MATRIX, SIZE(MATRIX,1,KIND=LONG), RHS, SIZE(RHS,1,KIND=LONG), JPVT, EPS1, RANK, WORK, LWORK, INFO )
      DEALLOCATE(WORK)

      END SUBROUTINE SOLVE_MINIMUM_NORM

! ##################################################################################################################################

      SUBROUTINE WRITE_RMG_TERMS ( TERM_COUNT )

      INTEGER(LONG), INTENT(OUT)     :: TERM_COUNT

      INTEGER(LONG)                  :: DEP_COMPONENT
      INTEGER(LONG)                  :: I
      INTEGER(LONG)                  :: J
      INTEGER(LONG)                  :: REF_COL(NUM_RIGID_DOF)
      INTEGER(LONG)                  :: REF_ROW_START
      INTEGER(LONG)                  :: RMG_ROW
      REAL(DOUBLE)                   :: COEFFICIENT_SCALE
      REAL(DOUBLE)                   :: WRITE_TOLERANCE

      TERM_COUNT = 0
      REF_COL = 0
      REF_ROW_START = TDOF_ROW_START(GRID_ID_ROW_NUM_D)

      DO I = 1, NUM_RETAINED
         DEP_COMPONENT = RETAINED_IDX(I)
         REF_COL(I) = TDOF(REF_ROW_START + DEP_COMPONENT - 1,G_SET_COL_NUM)
      ENDDO

      DO I = 1, NUM_RETAINED
         DEP_COMPONENT = RETAINED_IDX(I)
         RMG_ROW = TDOF(REF_ROW_START + DEP_COMPONENT - 1,M_SET_COL_NUM)
         IF ((RMG_ROW <= 0) .OR. (REF_COL(I) <= 0)) THEN
            IF (RMG_ROW <= 0) THEN
               WRITE(ERR,1509) SUBR_NAME, RTYPE, REID, AGRID_D, DEP_COMPONENT
               WRITE(F06,1509) SUBR_NAME, RTYPE, REID, AGRID_D, DEP_COMPONENT
            ENDIF
            IF (REF_COL(I) <= 0) THEN
               WRITE(ERR,1510) SUBR_NAME, RTYPE, REID, AGRID_D, DEP_COMPONENT
               WRITE(F06,1510) SUBR_NAME, RTYPE, REID, AGRID_D, DEP_COMPONENT
            ENDIF
            IERR = IERR + 1
            JERR = JERR + 1
            CYCLE
         ENDIF

         COEFFICIENT_SCALE = MAX(MAXVAL(ABS(A_REDUCED(I,1:NUM_RETAINED))), MAXVAL(ABS(B_REDUCED(I,1:NB6COLS))))
         WRITE_TOLERANCE = 10.0_DOUBLE*EPS1*COEFFICIENT_SCALE

         DO J = 1, NUM_RETAINED
            IF (ABS(A_REDUCED(I,J)) <= WRITE_TOLERANCE) CYCLE
            WRITE(L1J) RMG_ROW, REF_COL(J), A_REDUCED(I,J)
            TERM_COUNT = TERM_COUNT + 1
         ENDDO

         DO J = 1, NB6COLS
            IF (ABS(B_REDUCED(I,J)) <= WRITE_TOLERANCE) CYCLE
            WRITE(L1J) RMG_ROW, B6_COL_RMG(J), B_REDUCED(I,J)
            TERM_COUNT = TERM_COUNT + 1
         ENDDO
      ENDDO

 1509 FORMAT(' *ERROR  1509: PROGRAMMING ERROR IN SUBROUTINE ',A,/,15X,A8,' RIGID ELEMENT NUMBER ',I8,                           &
                    ', DEPENDENT GRID NUMBER ',I8,', COMPONENT ',I2,/,14X,' IS NOT AN M-SET DOF IN TABLE TDOFI')
 1510 FORMAT(' *ERROR  1510: PROGRAMMING ERROR IN SUBROUTINE ',A,/,15X,A8,' RIGID ELEMENT NUMBER ',I8,                           &
                    ', DEPENDENT GRID NUMBER ',I8,', COMPONENT ',I2,/,14X,' IS NOT A G-SET DOF IN TABLE TDOFI')

      END SUBROUTINE WRITE_RMG_TERMS

! ##################################################################################################################################

      PURE FUNCTION KINEMATIC_VECTOR ( COMPONENT, TDI, RELATIVE_POSITION ) RESULT ( H )

! Return the reference-motion coefficients for one scalar independent DOF. For translation, H = [D; R cross D]. For rotation,
! H = [0; D]. D is the independent component direction expressed in the reference grid's displacement coordinate system.

      INTEGER(LONG), INTENT(IN)      :: COMPONENT
      REAL(DOUBLE), INTENT(IN)       :: TDI(NUM_VECTOR_DOF,NUM_VECTOR_DOF)
      REAL(DOUBLE), INTENT(IN)       :: RELATIVE_POSITION(NUM_VECTOR_DOF)
      REAL(DOUBLE)                   :: DIRECTION(NUM_VECTOR_DOF)
      REAL(DOUBLE)                   :: H(NUM_RIGID_DOF)

      H = ZERO
      IF (COMPONENT <= NUM_VECTOR_DOF) THEN
         DIRECTION = TDI(:,COMPONENT)
         H(1:NUM_VECTOR_DOF) = DIRECTION
         H(4:NUM_RIGID_DOF) = CROSS_PRODUCT_3(RELATIVE_POSITION, DIRECTION)
      ELSE
         DIRECTION = TDI(:,COMPONENT - NUM_VECTOR_DOF)
         H(4:NUM_RIGID_DOF) = DIRECTION
      ENDIF

      END FUNCTION KINEMATIC_VECTOR

! ##################################################################################################################################

      PURE FUNCTION SCHUR_REDUCTION ( RETAINED, COUPLING, ELIMINATED_SOLUTION ) RESULT ( REDUCED )

      REAL(DOUBLE), INTENT(IN)       :: RETAINED(:,:)
      REAL(DOUBLE), INTENT(IN)       :: COUPLING(:,:)
      REAL(DOUBLE), INTENT(IN)       :: ELIMINATED_SOLUTION(:,:)
      REAL(DOUBLE)                   :: REDUCED(SIZE(RETAINED,1),SIZE(RETAINED,2))

      REDUCED = RETAINED - MATMUL(COUPLING, ELIMINATED_SOLUTION)

      END FUNCTION SCHUR_REDUCTION

! ##################################################################################################################################

      PURE FUNCTION MAXIMUM_RESIDUAL ( MATRIX, SOLUTION, RHS ) RESULT ( MAX_RESIDUAL )

      REAL(DOUBLE), INTENT(IN)       :: MATRIX(:,:)
      REAL(DOUBLE), INTENT(IN)       :: SOLUTION(:,:)
      REAL(DOUBLE), INTENT(IN)       :: RHS(:,:)
      REAL(DOUBLE)                   :: MAX_RESIDUAL
      REAL(DOUBLE)                   :: RESIDUAL(SIZE(RHS,1),SIZE(RHS,2))

      RESIDUAL = MATMUL(MATRIX, SOLUTION) - RHS
      MAX_RESIDUAL = MAXVAL(ABS(RESIDUAL))

      END FUNCTION MAXIMUM_RESIDUAL

! ##################################################################################################################################

      PURE FUNCTION CROSS_PRODUCT_3 ( LEFT, RIGHT ) RESULT ( PRODUCT )

      REAL(DOUBLE), INTENT(IN)       :: LEFT(NUM_VECTOR_DOF)
      REAL(DOUBLE), INTENT(IN)       :: RIGHT(NUM_VECTOR_DOF)
      REAL(DOUBLE)                   :: PRODUCT(NUM_VECTOR_DOF)

      PRODUCT = [ LEFT(2)*RIGHT(3) - LEFT(3)*RIGHT(2), &
                  LEFT(3)*RIGHT(1) - LEFT(1)*RIGHT(3), &
                  LEFT(1)*RIGHT(2) - LEFT(2)*RIGHT(1) ]

      END FUNCTION CROSS_PRODUCT_3

! ##################################################################################################################################

      PURE FUNCTION IDENTITY_3 () RESULT ( IDENTITY )

      REAL(DOUBLE)                   :: IDENTITY(NUM_VECTOR_DOF,NUM_VECTOR_DOF)

      IDENTITY = ZERO
      IDENTITY(1,1) = ONE
      IDENTITY(2,2) = ONE
      IDENTITY(3,3) = ONE

      END FUNCTION IDENTITY_3

! ##################################################################################################################################

      PURE FUNCTION OUTER_PRODUCT_6 ( VECTOR ) RESULT ( PRODUCT )

      REAL(DOUBLE), INTENT(IN)       :: VECTOR(NUM_RIGID_DOF)
      REAL(DOUBLE)                   :: PRODUCT(NUM_RIGID_DOF,NUM_RIGID_DOF)

      PRODUCT = SPREAD(VECTOR,DIM=2,NCOPIES=NUM_RIGID_DOF)*SPREAD(VECTOR,DIM=1,NCOPIES=NUM_RIGID_DOF)

      END FUNCTION OUTER_PRODUCT_6

! **********************************************************************************************************************************

      END SUBROUTINE RBE3_PROC
