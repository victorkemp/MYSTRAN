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
 
      SUBROUTINE FILE_CLOSE ( UNIT, FILNAM, CLOSE_STAT, WRITE_F04 )
 
! Closes files and writes message if the close fails
 
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  FILE_NAM_MAXLEN, WRT_ERR, WRT_LOG, F04, SC1
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  STIME, TSEC
      USE SUBR_BEGEND_LEVELS, ONLY    :  FILE_OPEN_BEGEND

      USE FILE_CLOSE_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'FILE_CLOSE'

!xx   CHARACTER(FILE_NAM_MAXLEN*BYTE), INTENT(IN) :: FILNAM            ! File name
      CHARACTER(LEN=*)   , INTENT(IN) :: FILNAM            ! File name

      CHARACTER(LEN=*)   , INTENT(IN) :: CLOSE_STAT        ! Status for close
      CHARACTER(LEN=*)   , INTENT(IN) :: WRITE_F04         ! If 'Y' write to F04, otherwise do not



      INTEGER(LONG), INTENT(IN)       :: UNIT              ! File unit number
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when closing a file
      INTEGER(LONG), PARAMETER        :: SUBR_BEGEND = FILE_OPEN_BEGEND



! **********************************************************************************************************************************
      CLOSE ( UNIT,STATUS=CLOSE_STAT,IOSTAT=IOCHK )

      IF (IOCHK /= 0) THEN
         WRITE(SC1,903) IOCHK
         CALL WRITE_FILNAM ( FILNAM, SC1, 15 )
         WRITE(SC1,9232)
         STOP
      ENDIF


      RETURN

! **********************************************************************************************************************************
  903 FORMAT(' *ERROR   903: ERROR ENCOUNTERED WITH IOSTAT = ',I8,' CLOSING FILE:')

 9232 FORMAT('               IT MAY BE OPEN IN ANOTHER PROGRAM. IF SO, CLOSE IT & START AGAIN.')

! **********************************************************************************************************************************

      END SUBROUTINE FILE_CLOSE
