! Begin MIT license text.                                                                                    
! _______________________________________________________________________________________________________
                                                                                                         
! Copyright 2019 Dr William R Case, Jr (dbcase29@gmail.com)                                              
                                                                                                         
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

      MODULE OUTPUT4_PROC_USE_IFs

! USE Interface statements for all subroutines called by SUBROUTINE OUTPUT4_PROC

      USE OURTIM_Interface
      USE GET_OU4_MAT_STATS_Interface
      USE OU4_PARTVEC_PROC_Interface
      USE WRITE_PARTNd_MAT_HDRS_Interface
      USE PARTITION_SS_NTERM_Interface
      USE ALLOCATE_SCR_CRS_MAT_Interface
      USE PARTITION_SS_Interface
      USE WRITE_OU4_SPARSE_MAT_Interface
      USE ALLOCATE_FULL_MAT_Interface
      USE PARTITION_FF_Interface
      USE WRITE_OU4_FULL_MAT_Interface
      USE OUTA_HERE_Interface
      USE DEALLOCATE_SCR_MAT_Interface
      USE DEALLOCATE_FULL_MAT_Interface

      END MODULE OUTPUT4_PROC_USE_IFs
