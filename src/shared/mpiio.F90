#ifdef OPT_IO_COCOMPI
! --- information -----------------------------------------------------
!
!  HISTORY
!     '15.04.07  M.Kurogi: MPI-IO
!
! ---------------------------------------------------------------------

      SUBROUTINE INFO_SEQ(FH, OFFSET, NSIZE0)
      use zocdim
      IMPLICIT NONE
#include "mpif.h"
      INTEGER NSIZE0, NSIZE
      INTEGER   FH
      INTEGER   (KIND = MPI_OFFSET_KIND):: OFFSET

      NSIZE=NSIZE0
      CALL REVERSE_INT4(NSIZE) !swap endian

      CALL MPI_FILE_SET_VIEW(                     &
     &     FH, OFFSET,                            &
     &     MPI_INTEGER4, MPI_INTEGER4,"native",   &
     &     MPI_INFO_NULL, IERR)
      IF (MYRANK .EQ. IROOT) THEN
         CALL MPI_FILE_WRITE(                     &
     &     FH, NSIZE, 1, MPI_INTEGER4,            &
     &     MPI_STATUS_IGNORE, IERR)
      END IF
      OFFSET = OFFSET + 4
      END



      SUBROUTINE MPI_WRITE_HEADER(CHEAD, FH, OFFSET)
      use zocdim
      IMPLICIT NONE
#include "mpif.h"

      CHARACTER CHEAD(64)*16
      INTEGER   FH
      INTEGER   (KIND = MPI_OFFSET_KIND):: OFFSET
      INTEGER NSIZE
      NSIZE=64*16
#ifdef OPT_IO_SEQUENTIAL
      CALL INFO_SEQ(FH, OFFSET, NSIZE)
#endif
!========= header
       CALL MPI_FILE_SET_VIEW(                   &
     &     FH, OFFSET,                           &
     &     MPI_CHARACTER,MPI_CHARACTER,"native", &
     &     MPI_INFO_NULL,IERR)

       IF (MYRANK .EQ. IROOT) THEN
       CALL MPI_FILE_WRITE(                      &
     &     FH,CHEAD,  NSIZE, MPI_CHARACTER,      &
     &     MPI_STATUS_IGNORE, IERR)
       END IF
      OFFSET = OFFSET + NSIZE
#ifdef OPT_IO_SEQUENTIAL
      CALL INFO_SEQ(FH, OFFSET, NSIZE)
#endif

      RETURN
      END SUBROUTINE MPI_WRITE_HEADER



!===============================================================
!     REVERSE_REAL4
!===============================================================

      SUBROUTINE REVERSE_REAL4(REAL4)
      IMPLICIT NONE

      REAL(4) :: REAL4, VALUE
      INTEGER(1) :: REVERSE(4), TMPVAL
      EQUIVALENCE(VALUE, REVERSE)

#ifdef OPT_IO_BYTESWAP
      REVERSE=0
      VALUE = REAL4

      TMPVAL = REVERSE(1)
      REVERSE(1) = REVERSE(4)
      REVERSE(4) = TMPVAL
      TMPVAL = REVERSE(2)
      REVERSE(2) = REVERSE(3)
      REVERSE(3) = TMPVAL

      REAL4 = VALUE
#endif
      RETURN
      END SUBROUTINE REVERSE_REAL4

!===============================================================
!     REVERSE_REAL8
!===============================================================

      SUBROUTINE REVERSE_REAL8(REAL8)
      IMPLICIT NONE

      real(8) :: REAL8, VALUE
      INTEGER(1) :: REVERSE(8), TMPVAL
      EQUIVALENCE(VALUE, REVERSE)

#ifdef OPT_IO_BYTESWAP
      REVERSE=0
      VALUE = REAL8

      TMPVAL = REVERSE(1)
      REVERSE(1) = REVERSE(8)
      REVERSE(8) = TMPVAL
      TMPVAL = REVERSE(2)
      REVERSE(2) = REVERSE(7)
      REVERSE(7) = TMPVAL
      TMPVAL = REVERSE(3)
      REVERSE(3) = REVERSE(6)
      REVERSE(6) = TMPVAL
      TMPVAL = REVERSE(4)
      REVERSE(4) = REVERSE(5)
      REVERSE(5) = TMPVAL

      REAL8 = VALUE
#endif
      RETURN
      END SUBROUTINE REVERSE_REAL8

!===============================================================
!     REVERSE_INT4
!===============================================================

      SUBROUTINE REVERSE_INT4(INT4)
      IMPLICIT NONE

      INTEGER   INT4, VALUE
      INTEGER(1) :: REVERSE(4), TMPVAL
      EQUIVALENCE(VALUE, REVERSE)
#ifdef OPT_IO_BYTESWAP
      REVERSE=0
      VALUE = INT4

      TMPVAL = REVERSE(1)
      REVERSE(1) = REVERSE(4)
      REVERSE(4) = TMPVAL
      TMPVAL = REVERSE(2)
      REVERSE(2) = REVERSE(3)
      REVERSE(3) = TMPVAL

      INT4 = VALUE
#endif
      RETURN
      END SUBROUTINE REVERSE_INT4



      SUBROUTINE MPI_READ_CHEAD(CHEAD, FH, DISP, ICREAD)
      use zocdim
      IMPLICIT NONE
#include "mpif.h"
      CHARACTER CHEAD(64)*16
      INTEGER FH
      INTEGER (KIND = MPI_OFFSET_KIND):: DISP
      INTEGER :: MPISTAT(MPI_STATUS_SIZE)
      INTEGER :: ICREAD
      INTEGER IFPAR, JFPAR

#ifdef OPT_IO_SEQUENTIAL
      DISP=DISP+4 
#endif
      CALL MPI_FILE_SET_VIEW(                      &
     &     FH, DISP,                               &
     &     MPI_CHARACTER, MPI_CHARACTER,"native",  &
     &     MPI_INFO_NULL,IERR)

      IF (MYRANK .EQ. IROOT) THEN
          CALL MPI_FILE_READ(                      &
     &    FH, CHEAD, 1024,                         &
     &    MPI_CHARACTER, MPISTAT, IERR)
          CALL MPI_GET_COUNT(MPISTAT,MPI_CHARACTER, ICREAD,IERR)     
      END IF
         CALL MPI_BCAST(CHEAD, 1024, MPI_CHARACTER,  &
     &                  IROOT, MPI_COMM_WORLD, IERR)
         CALL MPI_BCAST(ICREAD, 1, MPI_INTEGER4,     &
     &                  IROOT, MPI_COMM_WORLD, IERR)
      DISP=DISP+ 1024
#ifdef OPT_IO_SEQUENTIAL
      DISP=DISP + 4
#endif
      IF(ICREAD .NE. 1024) THEN
      DISP=DISP-1024
#ifdef OPT_IO_SEQUENTIAL
      DISP=DISP-8
#endif
      END IF

      RETURN
      END SUBROUTINE MPI_READ_CHEAD



      SUBROUTINE MPI_READ_SFC(BUF, FH, DISP)
      use zocdim
      IMPLICIT NONE

#include "mpif.h"

      real(8) ::  BUF(NX, NY)
      INTEGER FH
      INTEGER (KIND = MPI_OFFSET_KIND):: DISP
      INTEGER   I, J, K
      INTEGER :: IFILE
      INTEGER ISTART(2), IGSIZE(2), ISIZE(2)
      INTEGER IFPAR, JFPAR

#ifdef OPT_IO_SEQUENTIAL
      DISP=DISP+4 
#endif
                  ISTART=(/IRANK*NX, JRANK*NY/)
                  IGSIZE=(/NXG, NYG/)
                  ISIZE =(/NX , NY /)

                  CALL MPI_TYPE_CREATE_SUBARRAY(  &
     &                 2,      & !array dimension
     &                 IGSIZE, & !global size 
     &                 ISIZE,  & !subarray size
     &                 ISTART, & !subarray start index (start from 0)
     &                 MPI_ORDER_FORTRAN, &
     &                 MPI_REAL8,         &
     &                 IFILE,             &
     &                 IERR)
                  CALL MPI_TYPE_COMMIT(IFILE, IERR)

                  CALL MPI_FILE_SET_VIEW(        &
     &                 FH, DISP,                 &
     &                 MPI_REAL8,IFILE,"native", &
     &                 MPI_INFO_NULL,IERR)

                  CALL MPI_FILE_READ_ALL(        &
     &                 FH, BUF, NX*NY,           &
     &                 MPI_REAL8, MPI_STATUS_IGNORE, IERR)


      DO J=1,NY
      DO I=1,NX
            CALL REVERSE_REAL8(BUF(I,J))
      END DO
      END DO
      DISP=DISP+ NXG*NYG*8
#ifdef OPT_IO_SEQUENTIAL
      DISP=DISP+4 
#endif
      RETURN
      END SUBROUTINE MPI_READ_SFC



      SUBROUTINE MPI_READ_BDY(BUF, FH, DISP)
      use zocdim
      IMPLICIT NONE
#include "mpif.h"

      real(8) ::  BUF(NX, NY, NZ)
      INTEGER FH
      INTEGER (KIND = MPI_OFFSET_KIND):: DISP
      INTEGER   I, J, K
      INTEGER :: IFILE
      INTEGER ISTART(3), IGSIZE(3), ISIZE(3)
      INTEGER IFPAR, JFPAR
      INTEGER(8) :: INT1, INT2, INT3, INT4

#ifdef OPT_IO_SEQUENTIAL
      DISP=DISP+4 
#endif
                  ISTART=(/IRANK*NX, JRANK*NY, 0/)
                  IGSIZE=(/NXG, NYG, NZ/)
                  ISIZE =(/NX , NY , NZ /)

                  CALL MPI_TYPE_CREATE_SUBARRAY( &
     &                 3,       & !array dimension
     &                 IGSIZE,  & !global size 
     &                 ISIZE,   & !subarray size
     &                 ISTART,  & !subarray start index (start from 0)
     &                 MPI_ORDER_FORTRAN, &
     &                 MPI_REAL8, &
     &                 IFILE, &
     &                 IERR)
                  CALL MPI_TYPE_COMMIT(IFILE, IERR)

                  CALL MPI_FILE_SET_VIEW(        &
     &                 FH, DISP,                 &
     &                 MPI_REAL8,IFILE,"native", &
     &                 MPI_INFO_NULL,IERR)

                  CALL MPI_FILE_READ_ALL(        &
     &                 FH, BUF, NX*NY*NZ,        &
     &                 MPI_REAL8, MPI_STATUS_IGNORE, IERR)

      DO K=1,NZ
      DO J=1,NY
      DO I=1,NX
            CALL REVERSE_REAL8(BUF(I,J,K))
      END DO
      END DO
      END DO


#ifdef OPT_IO_SEQUENTIAL
      DISP=DISP+ NXG*NYG*NZ*8
      DISP=DISP+4 
#else
      INT1=NXG
      INT2=NYG
      INT3=NZ
      INT4=8
      DISP=DISP+ INT1*INT2*INT3*INT4
#endif

      RETURN
      END SUBROUTINE MPI_READ_BDY




      SUBROUTINE MPI_READ_ROOT(BUF,NBUF, FH, DISP)
      use zocdim
      IMPLICIT NONE
#include "mpif.h"

      real(8) ::  BUF(NBUF)
      INTEGER NBUF, FH, I
      INTEGER (KIND = MPI_OFFSET_KIND):: DISP
#ifdef OPT_IO_SEQUENTIAL
      DISP=DISP+4 
#endif
      CALL MPI_FILE_SET_VIEW( FH, DISP,    &
     &       MPI_REAL8,MPI_REAL8,"native", MPI_INFO_NULL,IERR)

      IF (MYRANK .EQ. IROOT) THEN
         call MPI_FILE_READ(FH, BUF, NBUF, &
     &                      MPI_REAL8, MPI_STATUS_IGNORE, IERR)
         DO I=1,NBUF
         call REVERSE_REAL8(BUF(I))
         END DO
      END IF
#ifdef OPT_IO_SEQUENTIAL
      DISP=DISP+ NBUF*8 + 4
#else
      DISP=DISP+ NBUF*8
#endif
      RETURN
      END SUBROUTINE MPI_READ_ROOT


      SUBROUTINE MPI_READ_2D_INTX(BUF, FH, DISP)
      use zocdim
      IMPLICIT NONE

#include "mpif.h"

      INTEGER  BUF(NXDIM, NYDIM)
      INTEGER  TMP(NX, NY)
      INTEGER FH
      INTEGER (KIND = MPI_OFFSET_KIND):: DISP
      INTEGER   I, J, K
      INTEGER :: IFILE
      INTEGER ISTART(2), IGSIZE(2), ISIZE(2)
#ifdef OPT_IO_SEQUENTIAL
      DISP=DISP+4 
#endif
                  ISTART=(/IRANK*NX +IGSTR-1, JRANK*NY+JGSTR-1/)
                  IGSIZE=(/NXGDIM, NYGDIM/)
                  ISIZE =(/NX    ,    NY /)

                  CALL MPI_TYPE_CREATE_SUBARRAY( &
     &                 2,      &   !array dimension
     &                 IGSIZE, & !global size 
     &                 ISIZE,  & !subarray size
     &                 ISTART, &  !subarray start index (start from 0)
     &                 MPI_ORDER_FORTRAN, &
     &                 MPI_INTEGER4, &
     &                 IFILE, &
     &                 IERR)
                  CALL MPI_TYPE_COMMIT(IFILE, IERR)

                  CALL MPI_FILE_SET_VIEW( &
     &                 FH, DISP,          &
     &                 MPI_INTEGER4,IFILE,"native", &
     &                 MPI_INFO_NULL,IERR)

                  CALL MPI_FILE_READ_ALL( &
     &                 FH, TMP, NX*NY,    &
     &                 MPI_INTEGER4, MPI_STATUS_IGNORE, IERR)



      DO J=1,NY
      DO I=1,NX
            CALL REVERSE_INT4(TMP(I,J))
            BUF(I+ISTR-1,J+JSTR-1)=TMP(I,J)
      END DO
      END DO
#ifdef OPT_IO_SEQUENTIAL
      DISP=DISP+ NXGDIM*NYGDIM*4 + 4
#else
      DISP=DISP+ NXGDIM*NYGDIM*4
#endif
      RETURN
      END SUBROUTINE MPI_READ_2D_INTX



      SUBROUTINE MPI_READ_2D_DIMX(BUF, FH, DISP)
      use zocdim
      IMPLICIT NONE
#include "mpif.h"

      real(8) ::  BUF(NXDIM, NYDIM)
      real(8) ::  TMP(NX, NY)
      INTEGER FH
      INTEGER (KIND = MPI_OFFSET_KIND):: DISP
      INTEGER   I, J, K
      INTEGER :: IFILE
      INTEGER ISTART(2), IGSIZE(2), ISIZE(2)

#ifdef OPT_IO_SEQUENTIAL
      DISP=DISP+4 
#endif
                  ISTART=(/IRANK*NX +IGSTR-1, JRANK*NY+JGSTR-1/)
                  IGSIZE=(/NXGDIM, NYGDIM/)
                  ISIZE =(/NX    ,    NY /)

                  CALL MPI_TYPE_CREATE_SUBARRAY(  &
     &                 2,      & !array dimension
     &                 IGSIZE, &  !global size 
     &                 ISIZE,  & !subarray size
     &                 ISTART, & !subarray start index (start from 0)
     &                 MPI_ORDER_FORTRAN, &
     &                 MPI_REAL8,         &
     &                 IFILE,             &
     &                 IERR)
                  CALL MPI_TYPE_COMMIT(IFILE, IERR)

                  CALL MPI_FILE_SET_VIEW( &
     &                 FH, DISP,          &
     &                 MPI_REAL8,IFILE,"native", &
     &                 MPI_INFO_NULL,IERR)

                  CALL MPI_FILE_READ_ALL( &
     &                 FH, TMP, NX*NY,    &
     &                 MPI_REAL8, MPI_STATUS_IGNORE, IERR)



      DO J=1,NY
      DO I=1,NX
            CALL REVERSE_REAL8(TMP(I,J))
            BUF(I+ISTR-1,J+JSTR-1)=TMP(I,J)
      END DO
      END DO

#ifdef OPT_IO_SEQUENTIAL
      DISP=DISP+ NXGDIM*NYGDIM*8 + 4
#else
      DISP=DISP+ NXGDIM*NYGDIM*8
#endif
      RETURN
      END SUBROUTINE MPI_READ_2D_DIMX




      SUBROUTINE MPI_READ_3D_DIMX(BUF, FH, DISP)
      use zocdim
      IMPLICIT NONE

#include "mpif.h"

      real(8) ::  BUF(NXDIM, NYDIM, NZDIM)
      real(8) ::  TMP(NX, NY, NZDIM)
      INTEGER FH
      INTEGER (KIND = MPI_OFFSET_KIND):: DISP
      INTEGER   I, J, K
      INTEGER :: IFILE
      INTEGER ISTART(3), IGSIZE(3), ISIZE(3)
      INTEGER(8) :: INT1, INT2
#ifdef OPT_IO_SEQUENTIAL
      DISP=DISP+4 
#endif
                  ISTART=(/IRANK*NX +IGSTR-1, JRANK*NY+JGSTR-1, 0/)
                  IGSIZE=(/NXGDIM, NYGDIM, NZDIM/)
                  ISIZE =(/NX    ,    NY , NZDIM/)

                  CALL MPI_TYPE_CREATE_SUBARRAY( &
     &                 3,       & !array dimension
     &                 IGSIZE,  & !global size 
     &                 ISIZE,   & !subarray size
     &                 ISTART,  & !subarray start index (start from 0)
     &                 MPI_ORDER_FORTRAN, &
     &                 MPI_REAL8, &
     &                 IFILE,     &
     &                 IERR)
                  CALL MPI_TYPE_COMMIT(IFILE, IERR)

                  CALL MPI_FILE_SET_VIEW( &
     &                 FH, DISP,          &
     &                 MPI_REAL8,IFILE,"native", &
     &                 MPI_INFO_NULL,IERR)

                  CALL MPI_FILE_READ_ALL(    &
     &                 FH, TMP, NX*NY*NZDIM, &
     &                 MPI_REAL8, MPI_STATUS_IGNORE, IERR)



      DO K=1,NZDIM
      DO J=1,NY
      DO I=1,NX
            CALL REVERSE_REAL8(TMP(I,J,K))
            BUF(I+ISTR-1,J+JSTR-1,K)=TMP(I,J,K)
      END DO
      END DO
      END DO
#ifdef OPT_IO_SEQUENTIAL
      DISP=DISP+ NXYZGD*8 + 4
#else
      INT1=NXYZGD
      INT2=8
      DISP=DISP+ INT1*INT2
#endif
      RETURN
      END SUBROUTINE MPI_READ_3D_DIMX



      SUBROUTINE MPI_READ_2D(BUF, FH, DISP)
      use zocdim
      IMPLICIT NONE
#include "mpif.h"

      real(8) ::  BUF(NXDIM, NYDIM)
      real(8) ::  TMP(NX, NY)
      INTEGER FH
      INTEGER (KIND = MPI_OFFSET_KIND):: DISP
      INTEGER   I, J, K
      INTEGER :: IFILE
      INTEGER ISTART(2), IGSIZE(2), ISIZE(2)
      INTEGER  IFPAR,  JFPAR

#ifdef OPT_IO_SEQUENTIAL
      DISP=DISP+4 
#endif
                  ISTART=(/IRANK*NX, JRANK*NY/)
                  IGSIZE=(/NXG, NYG/)
                  ISIZE =(/NX , NY /)

                  CALL MPI_TYPE_CREATE_SUBARRAY( &
     &                 2,      &  !array dimension
     &                 IGSIZE, & !global size 
     &                 ISIZE,  & !subarray size
     &                 ISTART, & !subarray start index (start from 0)
     &                 MPI_ORDER_FORTRAN, &
     &                 MPI_REAL8,         &
     &                 IFILE,             &
     &                 IERR)
  CALL MPI_TYPE_COMMIT(IFILE, IERR)
  CALL MPI_FILE_SET_VIEW(FH, DISP, MPI_REAL8,IFILE,"native", MPI_INFO_NULL,IERR)
  CALL MPI_FILE_READ_ALL(FH, TMP, NX*NY, MPI_REAL8, MPI_STATUS_IGNORE, IERR)

      DO J=1,NY
      DO I=1,NX
            CALL REVERSE_REAL8(TMP(I,J))
      END DO
      END DO
      BUF(ISTR:IEND,JSTR:JEND)=TMP(:,:)

      DISP=DISP+ NXG*NYG*8
#ifdef OPT_IO_SEQUENTIAL
      DISP=DISP+4 
#endif
      RETURN
      END SUBROUTINE MPI_READ_2D


      SUBROUTINE MPI_READ_ID(BUF, FH, DISP)
      use zocdim
      IMPLICIT NONE
#include "mpif.h"

      real(8) ::  TMP(NX, NY, NIC)
      real(8) ::  BUF(NXDIM, NYDIM, 0:NIC)
      INTEGER FH
      INTEGER (KIND = MPI_OFFSET_KIND):: DISP
      INTEGER   I, J, K
      INTEGER :: IFILE
      INTEGER ISTART(3), IGSIZE(3), ISIZE(3)

#ifdef OPT_IO_SEQUENTIAL
      DISP=DISP+4 
#endif
                  ISTART=(/IRANK*NX, JRANK*NY, 0/)
                  IGSIZE=(/NXG, NYG, NIC/)
                  ISIZE =(/NX , NY , NIC/)

  CALL MPI_TYPE_CREATE_SUBARRAY(3, IGSIZE, ISIZE, ISTART, &
     &   MPI_ORDER_FORTRAN, MPI_REAL8, IFILE, IERR)
  CALL MPI_TYPE_COMMIT(IFILE, IERR)
  CALL MPI_FILE_SET_VIEW(FH, DISP, MPI_REAL8,IFILE,"native",MPI_INFO_NULL,IERR)
  CALL MPI_FILE_READ_ALL(FH, TMP, NX*NY*NIC, MPI_REAL8, MPI_STATUS_IGNORE, IERR)



      DO K=1,NIC
      DO J=1,NY
      DO I=1,NX
            CALL REVERSE_REAL8(TMP(I,J,K))
            BUF(I+ISTR-1,J+JSTR-1,K)=TMP(I,J,K)
      END DO
      END DO
      END DO
      DISP=DISP+ NXG*NYG*NIC*8
#ifdef OPT_IO_SEQUENTIAL
      DISP=DISP+4 
#endif
      RETURN
      END SUBROUTINE MPI_READ_ID



      SUBROUTINE MPI_READ_3D(BUF, FH, DISP)
      use zocdim
      IMPLICIT NONE
#include "mpif.h"

      real(8) ::  BUF(NXDIM, NYDIM, NZDIM)
      real(8) ::  TMP(NX, NY, NZ)
      INTEGER FH
      INTEGER (KIND = MPI_OFFSET_KIND):: DISP
      INTEGER   I, J, K
      INTEGER :: IFILE
      INTEGER ISTART(3), IGSIZE(3), ISIZE(3)
      INTEGER(8) :: INT1, INT2, INT3, INT4

#ifdef OPT_IO_SEQUENTIAL
      DISP=DISP+4 
#endif
                  ISTART=(/IRANK*NX, JRANK*NY, 0/)
                  IGSIZE=(/NXG, NYG, NZ/)
                  ISIZE =(/NX , NY , NZ/)

  CALL MPI_TYPE_CREATE_SUBARRAY(3,IGSIZE, ISIZE, ISTART, &
     &        MPI_ORDER_FORTRAN, MPI_REAL8, IFILE, IERR)
  CALL MPI_TYPE_COMMIT(IFILE, IERR)
  CALL MPI_FILE_SET_VIEW(FH, DISP, MPI_REAL8,IFILE,"native", MPI_INFO_NULL,IERR)
  CALL MPI_FILE_READ_ALL(FH, TMP, NX*NY*NZ, MPI_REAL8, MPI_STATUS_IGNORE, IERR)

      DO K=1,NZ
      DO J=1,NY
      DO I=1,NX
            CALL REVERSE_REAL8(TMP(I,J,K))
            BUF(I+ISTR-1, J+JSTR-1, K+KSTR-1)=TMP(I,J,K)
      END DO
      END DO
      END DO

#ifdef OPT_IO_SEQUENTIAL
      DISP=DISP+ NXG*NYG*NZ*8
      DISP=DISP+4 
#else
      INT1=NXG
      INT2=NYG
      INT3=NZ
      INT4=8
      DISP=DISP+ INT1*INT2*INT3*INT4
#endif
      RETURN
      END SUBROUTINE MPI_READ_3D





      SUBROUTINE MPI_WRITE_2D(BUF, FH, DISP)
      use zocdim
      IMPLICIT NONE
#include "mpif.h"

      real(8) ::  BUF(NXDIM, NYDIM)
      real(8) ::  TMP(NX, NY)
      INTEGER FH
      INTEGER (KIND = MPI_OFFSET_KIND):: DISP
      INTEGER   I, J, K, KDIM
      INTEGER :: IFILE
      INTEGER ISTART(2), IGSIZE(2), ISIZE(2)
      INTEGER NSIZE

      NSIZE=8*NXG*NYG
#ifdef OPT_IO_SEQUENTIAL
      CALL INFO_SEQ(FH, DISP, NSIZE)
#endif

      DO J=1,NY
      DO I=1,NX
            TMP(I,J)=BUF(I+ISTR-1,J+JSTR-1)
            CALL REVERSE_REAL8(TMP(I,J))
      END DO
      END DO

                  ISTART=(/IRANK*NX, JRANK*NY/)
                  IGSIZE=(/NXG, NYG/)
                  ISIZE =(/NX , NY /)

  CALL MPI_TYPE_CREATE_SUBARRAY(2, IGSIZE, ISIZE, ISTART, &
     &       MPI_ORDER_FORTRAN, MPI_REAL8, IFILE, IERR)
  CALL MPI_TYPE_COMMIT(IFILE, IERR)
  CALL MPI_FILE_SET_VIEW(FH, DISP, MPI_REAL8, IFILE, "native", MPI_INFO_NULL, IERR)
  CALL MPI_FILE_WRITE_ALL(FH, TMP, NX*NY, MPI_REAL8, MPI_STATUS_IGNORE, IERR)

      DISP=DISP+NSIZE
#ifdef OPT_IO_SEQUENTIAL
      CALL INFO_SEQ(FH, DISP, NSIZE)
#endif
      RETURN
      END SUBROUTINE MPI_WRITE_2D


      SUBROUTINE MPI_WRITE_ID(BUF, FH, DISP)
      use zocdim
      IMPLICIT NONE
#include "mpif.h"

      real(8) ::  BUF(NXDIM, NYDIM, 0:NIC)
      real(8) ::  TMP(NX, NY, NIC)
      INTEGER FH
      INTEGER (KIND = MPI_OFFSET_KIND):: DISP
      INTEGER   I, J, K, KDIM
      INTEGER :: IFILE
      INTEGER ISTART(3), IGSIZE(3), ISIZE(3)
      INTEGER NSIZE
      INTEGER(8) :: INT1, INT2, INT3, INT4,  NSIZE2

      INT1=8
      INT2=NXG
      INT3=NYG
      INT4=NIC
      NSIZE2=INT1*INT2*INT3*INT4      

      NSIZE= 8*NXG*NYG*NIC
#ifdef OPT_IO_SEQUENTIAL
      CALL INFO_SEQ(FH, DISP, NSIZE)
#endif
      DO K=1,NIC
      DO J=1,NY
      DO I=1,NX
            TMP(I,J,K)=BUF(I+ISTR-1,J+JSTR-1,K)
            CALL REVERSE_REAL8(TMP(I,J,K))
      END DO
      END DO
      END DO

                  ISTART=(/IRANK*NX, JRANK*NY, 0/)
                  IGSIZE=(/NXG, NYG, NIC/)
                  ISIZE =(/NX , NY , NIC/)

  CALL MPI_TYPE_CREATE_SUBARRAY(3, IGSIZE, ISIZE, &
     &    ISTART, MPI_ORDER_FORTRAN, MPI_REAL8, IFILE, IERR)
  CALL MPI_TYPE_COMMIT(IFILE, IERR)
  CALL MPI_FILE_SET_VIEW(FH, DISP,MPI_REAL8,IFILE,"native",MPI_INFO_NULL,IERR)
  CALL MPI_FILE_WRITE_ALL(FH, TMP, NX*NY*NIC, MPI_REAL8, MPI_STATUS_IGNORE, IERR)
      DISP=DISP+NSIZE2
#ifdef OPT_IO_SEQUENTIAL
      CALL INFO_SEQ(FH, DISP, NSIZE)
#endif
      RETURN
      END SUBROUTINE MPI_WRITE_ID


      SUBROUTINE MPI_WRITE_3D(BUF, FH, DISP)
      use zocdim
      IMPLICIT NONE
#include "mpif.h"

      real(8) ::  BUF(NXDIM, NYDIM, NZDIM)
      real(8) ::  TMP(NX, NY, NZ)
      INTEGER FH
      INTEGER (KIND = MPI_OFFSET_KIND):: DISP
      INTEGER   I, J, K, KDIM
      INTEGER :: IFILE
      INTEGER ISTART(3), IGSIZE(3), ISIZE(3)
      INTEGER NSIZE
      INTEGER(8) :: INT1, INT2, INT3, INT4,  NSIZE2

      INT1=8
      INT2=NXG
      INT3=NYG
      INT4=NZ
      NSIZE2=INT1*INT2*INT3*INT4

      NSIZE=8*NXG*NYG*NZ
#ifdef OPT_IO_SEQUENTIAL
      CALL INFO_SEQ(FH, DISP, NSIZE)
#endif

      DO K=1,NZ
      DO J=1,NY
      DO I=1,NX
            TMP(I,J,K)=BUF(I+ISTR-1,J+JSTR-1,K+KSTR-1)
            CALL REVERSE_REAL8(TMP(I,J,K))
      END DO
      END DO
      END DO

                  ISTART=(/IRANK*NX, JRANK*NY, 0/)
                  IGSIZE=(/NXG, NYG, NZ/)
                  ISIZE =(/NX , NY , NZ/)

  CALL MPI_TYPE_CREATE_SUBARRAY( 3, IGSIZE, ISIZE, ISTART, &
     & MPI_ORDER_FORTRAN, MPI_REAL8, IFILE, IERR)
  CALL MPI_TYPE_COMMIT(IFILE, IERR)
  CALL MPI_FILE_SET_VIEW(FH, DISP, MPI_REAL8,IFILE,"native",MPI_INFO_NULL,IERR)
  CALL MPI_FILE_WRITE_ALL(FH, TMP, NX*NY*NZ, MPI_REAL8, MPI_STATUS_IGNORE, IERR)
      DISP=DISP+NSIZE2
#ifdef OPT_IO_SEQUENTIAL
      CALL INFO_SEQ(FH, DISP, NSIZE)
#endif
      RETURN
      END SUBROUTINE MPI_WRITE_3D

#else
      SUBROUTINE MPI_IO
      RETURN
      END
#endif
