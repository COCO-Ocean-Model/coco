#ifdef OPT_PARALLEL
      SUBROUTINE SHFT1I(
     M                      IQ,
     D                    IDIM,   JDIM)
      use zocdim
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  One-on-one communication routines for a 2-dimensional integer
c variable.
c
c  HISTORY
c     '99.10.06  H.Hasumi
c     '07.04.23  H.Hasumi
c     '10.04.14  M.Kurogi: (COCO4.4 tripolar code by Dr. Suzuki)
c
c ---------------------------------------------------------------------

#include "mpif.h"

      INTEGER   IDIM,   JDIM
      INTEGER     IQ(IDIM, JDIM)

      INTEGER ISDBX1(ICOMM, NYDIM), ISDBX2(ICOMM, NYDIM)
      INTEGER ISDBY1(NXDIM, JCOMM), ISDBY2(NXDIM, JCOMM)
      INTEGER IRVBX1(ICOMM, NYDIM), IRVBX2(ICOMM, NYDIM)
      INTEGER IRVBY1(NXDIM, JCOMM), IRVBY2(NXDIM, JCOMM)
      INTEGER      I,      J
      INTEGER NBFDIM

      IF (IDOWN .NE. MPI_PROC_NULL) THEN
         DO J = 1, NYDIM
            DO I = 1, ICOMM
               ISDBX1(I, J) = IQ(I+ISTR-1, J)
            END DO
         END DO
      END IF
      IF (IUP .NE. MPI_PROC_NULL) THEN
         DO J = 1, NYDIM
            DO I = 1, ICOMM
               ISDBX2(I, J) = IQ(I+IEND-ICOMM, J)
            END DO
         END DO
      END IF
      NBFDIM = NYDIM * ICOMM
      CALL SHFTXI(
     O            IRVBX1, IRVBX2,
     I            ISDBX1, ISDBX2,
     D            NBFDIM)
      IF (IDOWN .NE. MPI_PROC_NULL) THEN
         DO J = 1, NYDIM
            DO I = 1, ICOMM
               IQ(I, J) = IRVBX1(I, J)
            END DO
         END DO
      END IF
      IF (IUP .NE. MPI_PROC_NULL) THEN
         DO J = 1, NYDIM
            DO I = 1, ICOMM
               IQ(IEND+I, J) = IRVBX2(I, J)
            END DO
         END DO
      END IF

      IF (JDOWN .NE. MPI_PROC_NULL) THEN
         DO J = 1, JCOMM
            DO I = 1, NXDIM
               ISDBY1(I, J) = IQ(I, J+JSTR-1)
            END DO
         END DO
      END IF
      IF (JUP .NE. MPI_PROC_NULL) THEN
         DO J = 1, JCOMM
            DO I = 1, NXDIM
               ISDBY2(I, J) = IQ(I, J+JEND-JCOMM)
            END DO
         END DO
      END IF
      NBFDIM = NXDIM * JCOMM
      CALL SHFTYI(
     O            IRVBY1, IRVBY2,
     I            ISDBY1, ISDBY2,
     D            NBFDIM)
      IF (JDOWN .NE. MPI_PROC_NULL) THEN
         DO J = 1, JCOMM
            DO I = 1, NXDIM
               IQ(I, J) = IRVBY1(I, J)
            END DO
         END DO
      END IF
      IF (JUP .NE. MPI_PROC_NULL) THEN
         DO J = 1, JCOMM
            DO I = 1, NXDIM
               IQ(I, JEND+J) = IRVBY2(I, J)
            END DO
         END DO
      END IF

      RETURN
      END

c *********************************************************************

      SUBROUTINE SHFTXI(
     O            IRVBX1, IRVBX2,
     I            ISDBX1, ISDBX2,
     D            NBFDIM)
      use zocdim
      IMPLICIT NONE

#include "mpif.h"

      INTEGER ISDBX1(ICOMM, NYDIM), ISDBX2(ICOMM, NYDIM)
      INTEGER IRVBX1(ICOMM, NYDIM), IRVBX2(ICOMM, NYDIM)
      INTEGER NBFDIM

      INTEGER ISRQX1, ISRQX2
      INTEGER IRRQX1, IRRQX2
      INTEGER ISTMPI(MPI_STATUS_SIZE)

      CALL MPI_ISEND(
     &               ISDBX1, NBFDIM, MPI_INTEGER,
     &                IDOWN,      1, MPI_COMM_WORLD,
     &               ISRQX1,   IERR)
      CALL MPI_ISEND(
     &               ISDBX2, NBFDIM, MPI_INTEGER,
     &                  IUP,      2, MPI_COMM_WORLD,
     &               ISRQX2,   IERR)
      CALL MPI_IRECV(
     &               IRVBX2, NBFDIM, MPI_INTEGER,
     &                  IUP,      1, MPI_COMM_WORLD,
     &               IRRQX1,   IERR)
      CALL MPI_IRECV(
     &               IRVBX1, NBFDIM, MPI_INTEGER,
     &                IDOWN,      2, MPI_COMM_WORLD,
     &               IRRQX2,   IERR)
      CALL MPI_WAIT(ISRQX1, ISTMPI,   IERR)
      CALL MPI_WAIT(ISRQX2, ISTMPI,   IERR)
      CALL MPI_WAIT(IRRQX1, ISTMPI,   IERR)
      CALL MPI_WAIT(IRRQX2, ISTMPI,   IERR)

      RETURN
      END

c *********************************************************************

      SUBROUTINE SHFTYI(
     O            IRVBY1, IRVBY2,
     I            ISDBY1, ISDBY2,
     D            NBFDIM)
      use zocdim
      IMPLICIT NONE

#include "mpif.h"

      INTEGER ISDBY1(NXDIM, JCOMM), ISDBY2(NXDIM, JCOMM)
      INTEGER IRVBY1(NXDIM, JCOMM), IRVBY2(NXDIM, JCOMM)
      INTEGER NBFDIM

      INTEGER ISRQY1, ISRQY2
      INTEGER IRRQY1, IRRQY2
      INTEGER ISTMPI(MPI_STATUS_SIZE)

      CALL MPI_ISEND(
     &               ISDBY1, NBFDIM, MPI_INTEGER,
     &                JDOWN,      3, MPI_COMM_WORLD,
     &               ISRQY1,   IERR)
      CALL MPI_ISEND(
     &               ISDBY2, NBFDIM, MPI_INTEGER,
     &                  JUP,      4, MPI_COMM_WORLD,
     &               ISRQY2,   IERR)
      CALL MPI_IRECV(
     &               IRVBY2, NBFDIM, MPI_INTEGER,
     &                  JUP,      3, MPI_COMM_WORLD,
     &               IRRQY1,   IERR)
      CALL MPI_IRECV(
     &               IRVBY1, NBFDIM, MPI_INTEGER,
     &                JDOWN,      4, MPI_COMM_WORLD,
     &               IRRQY2,   IERR)
      CALL MPI_WAIT(ISRQY1, ISTMPI,   IERR)
      CALL MPI_WAIT(ISRQY2, ISTMPI,   IERR)
      CALL MPI_WAIT(IRRQY1, ISTMPI,   IERR)
      CALL MPI_WAIT(IRRQY2, ISTMPI,   IERR)

      RETURN
      END

#ifdef OPT_TRIPOLE
      SUBROUTINE SHFTINT(
     M                      IQ,
     D                    IDIM,   JDIM)
      use zocdim
      IMPLICIT NONE

#include "mpif.h"

      INTEGER   IDIM,   JDIM
      INTEGER     IQ(IDIM, JDIM)

      INTEGER ISDBX1(ICOMM, NY), ISDBX2(ICOMM, NY)
      INTEGER ISDBY1(NXDIM, JCOMM), ISDBY2(NXDIM, JCOMM)
      INTEGER IRVBX1(ICOMM, NY), IRVBX2(ICOMM, NY)
      INTEGER IRVBY1(NXDIM, JCOMM), IRVBY2(NXDIM, JCOMM)
      INTEGER      I,      J
      INTEGER NBFDIM

      IF (JUPE .NE. MPI_PROC_NULL) THEN
         DO J = 1, JCOMM
            DO I = 1, NXDIM
               ISDBY1(I, J) = IQ(NXDIM-I+1,JEND-J+1)
            END DO
         END DO
      END IF
      IF (JUPW .NE. MPI_PROC_NULL) THEN
         DO J = 1, JCOMM
            DO I = 1, NXDIM
               ISDBY2(I, J) = IQ(NXDIM-I+1,JEND-J+1)
            END DO
         END DO
      END IF
      NBFDIM = NXDIM * JCOMM
      CALL SHFTYIN(
     O            IRVBY1, IRVBY2,
     I            ISDBY1, ISDBY2,
     D            NBFDIM)
      IF (JUPE .NE. MPI_PROC_NULL) THEN
         DO J = 1, JCOMM
            DO I = 1, NXDIM
               IQ(I, JEND+J) = IRVBY1(I, J)
            END DO
         END DO
      END IF
      IF (JUPW .NE. MPI_PROC_NULL) THEN
         DO J = 1, JCOMM
            DO I = 1, NXDIM
               IQ(I, JEND+J) = IRVBY2(I, J)
            END DO
         END DO
      END IF

      RETURN
      END
CC
CC
      SUBROUTINE SHFTINV(
     M                      IQ,
     D                    IDIM,   JDIM)
      use zocdim
      IMPLICIT NONE

#include "mpif.h"

      INTEGER   IDIM,   JDIM
      INTEGER     IQ(IDIM, JDIM)

      INTEGER ISDBX1(ICOMM, NY), ISDBX2(ICOMM, NY)
      INTEGER ISDBY1(NXDIM, JCOMM), ISDBY2(NXDIM, JCOMM)
      INTEGER IRVBX1(ICOMM, NY), IRVBX2(ICOMM, NY)
      INTEGER IRVBY1(NXDIM, JCOMM), IRVBY2(NXDIM, JCOMM)
      INTEGER      I,      J
      INTEGER NBFDIM

      IF (JUPE .NE. MPI_PROC_NULL) THEN
         DO J = 1, JCOMM
            DO I = 1, NXDIM-1
               ISDBY1(I, J) = IQ(NXDIM-I,JEND-J)
            END DO
               ISDBY1(NXDIM, J) = IQ(1,JEND-J)
         END DO
      END IF
      IF (JUPW .NE. MPI_PROC_NULL) THEN
         DO J = 1, JCOMM
            DO I = 1, NXDIM-1
               ISDBY2(I, J) = IQ(NXDIM-I,JEND-J)
            END DO
               ISDBY2(NXDIM, J) = IQ(1,JEND-J)
         END DO
      END IF
      NBFDIM = NXDIM * JCOMM
      CALL SHFTYIN(
     O            IRVBY1, IRVBY2,
     I            ISDBY1, ISDBY2,
     D            NBFDIM)
      IF (JUPE .NE. MPI_PROC_NULL) THEN
         DO J = 1, JCOMM
            DO I = 1, NXDIM
               IQ(I, JEND+J) = IRVBY1(I, J)
            END DO
         END DO
      END IF
      IF (JUPW .NE. MPI_PROC_NULL) THEN
         DO J = 1, JCOMM
            DO I = 1, NXDIM
               IQ(I, JEND+J) = IRVBY2(I, J)
            END DO
         END DO
      END IF

      RETURN
      END
c *********************************************************************

      SUBROUTINE SHFTYIN(
     O            IRVBY1, IRVBY2,
     I            ISDBY1, ISDBY2,
     D            NBFDIM)
      use zocdim
      IMPLICIT NONE

#include "mpif.h"

      INTEGER ISDBY1(NXDIM, JCOMM), ISDBY2(NXDIM, JCOMM)
      INTEGER IRVBY1(NXDIM, JCOMM), IRVBY2(NXDIM, JCOMM)
      INTEGER NBFDIM

      INTEGER ISRQY1, ISRQY2
      INTEGER IRRQY1, IRRQY2
      INTEGER ISTMPI(MPI_STATUS_SIZE)

      CALL MPI_ISEND(
     &               ISDBY1, NBFDIM, MPI_INTEGER,
     &                 JUPE,      3, MPI_COMM_WORLD,
     &               ISRQY1,   IERR)
      CALL MPI_ISEND(
     &               ISDBY2, NBFDIM, MPI_INTEGER,
     &                 JUPW,      4, MPI_COMM_WORLD,
     &               ISRQY2,   IERR)
      CALL MPI_IRECV(
     &               IRVBY2, NBFDIM, MPI_INTEGER,
     &                 JUPW,      3, MPI_COMM_WORLD,
     &               IRRQY1,   IERR)
      CALL MPI_IRECV(
     &               IRVBY1, NBFDIM, MPI_INTEGER,
     &                 JUPE,      4, MPI_COMM_WORLD,
     &               IRRQY2,   IERR)
      CALL MPI_WAIT(ISRQY1, ISTMPI,   IERR)
      CALL MPI_WAIT(ISRQY2, ISTMPI,   IERR)
      CALL MPI_WAIT(IRRQY1, ISTMPI,   IERR)
      CALL MPI_WAIT(IRRQY2, ISTMPI,   IERR)

      RETURN
      END
#endif

#else
c *********************************************************************
      SUBROUTINE SHFIST
      RETURN
      END
#endif
