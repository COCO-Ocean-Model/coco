#ifdef OPT_PARALLEL
      SUBROUTINE GS3DST
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Group communication (GATHER and SCATTER) routines for 3-dimensional
c variables.
c
c  HISTORY
c     '99.10.05  H.Hasumi
c     '07.04.23  H.Hasumi
c
c ---------------------------------------------------------------------

#include "mpif.h"
#include "zocdim.F"

      REAL*8       Q(NXDIM , NYDIM , NZDIM)
      REAL*8      QG(NXGDIM, NYGDIM, NZDIM)
      REAL*8    QCHK(NX         , NY         , KSTR:KEND)
      REAL*8   QCHKG(IGSTR:IGEND, JGSTR:JGEND, KSTR:KEND)

      INTEGER, ALLOCATABLE :: IRSTR(:), JRSTR(:)
      INTEGER, ALLOCATABLE :: IREND(:), JREND(:)
      INTEGER, ALLOCATABLE :: NSNDC(:), NRCVC(:), NDISP(:)
      INTEGER  NGSND,  NSRCV
      SAVE IRSTR, JRSTR, IREND, JREND, NSNDC, NRCVC, NDISP,
     &     NGSND, NSRCV

      REAL*8,  ALLOCATABLE :: SNDBUF(:, :, :), RCVBUF(:, :, :)
      REAL*8,  ALLOCATABLE :: DUMMY(:)
      SAVE SNDBUF, RCVBUF, DUMMY

      INTEGER      I,      J,   IJKG,      K,      N

      ALLOCATE(IRSTR(0:NPROCS-1))
      ALLOCATE(JRSTR(0:NPROCS-1))
      ALLOCATE(IREND(0:NPROCS-1))
      ALLOCATE(JREND(0:NPROCS-1))
      ALLOCATE(NSNDC(0:NPROCS-1))
      ALLOCATE(NRCVC(0:NPROCS-1))
      ALLOCATE(NDISP(0:NPROCS-1))
      ALLOCATE(DUMMY(NXYZG+NPROCS))
      IF (MYRANK .GE. IJNODE) THEN
         ALLOCATE(SNDBUF(1, 1, 1))
         ALLOCATE(RCVBUF(1, 1, 1))
      ELSE
         ALLOCATE(SNDBUF(NX, NY, KSTR:KEND))
         ALLOCATE(RCVBUF(NX, NY, KSTR:KEND))
      END IF

      SNDBUF = 0.D0
      RCVBUF = 0.D0
      DUMMY  = 0.D0

      DO J = 0, JNODES-1
         DO I = 0, INODES-1
            N = J * INODES + I
            IRSTR(N) = IGSTR + I * NX
            IREND(N) = IRSTR(N) + NX - 1
            JRSTR(N) = JGSTR + J * NY
            JREND(N) = JRSTR(N) + NY - 1
            NSNDC(N) = NXYZ
            NRCVC(N) = NXYZ
         END DO
      END DO
      DO N = IJNODE, NPROCS-1
         IRSTR(N) = 2
         IREND(N) = 1
         JRSTR(N) = 2
         JREND(N) = 1
         NSNDC(N) = 1
         NRCVC(N) = 1
      END DO
      NDISP(0) = 0
      DO N = 1, NPROCS-1
         NDISP(N) = NDISP(N-1) + NSNDC(N-1)
      END DO

      IF (MYRANK .GE. IJNODE) THEN
         NGSND = 1
         NSRCV = 1
      ELSE
         NGSND = NXYZ
         NSRCV = NXYZ
      END IF

      RETURN

c =====================================================================

      ENTRY GATHER_3D(
     O                    QG,
     I                     Q)

      IF (MYRANK .LT. IJNODE) THEN
         DO K = KSTR, KEND
            DO J = 1, NY
               DO I = 1, NX
                  SNDBUF(I, J, K) = Q(I+ISTR-1, J+JSTR-1, K)
               END DO
            END DO
         END DO
      END IF

      CALL MPI_GATHERV(
     &                 SNDBUF,  NGSND,         MPI_REAL8,
     &                  DUMMY,  NRCVC,  NDISP, MPI_REAL8,
     &                  IROOT, MPI_COMM_WORLD, IERR)

      IF (MYRANK .EQ. IROOT) THEN
         DO N = 0, NPROCS-1
            IJKG = NDISP(N)
            DO K = KSTR, KEND
               DO J = JRSTR(N), JREND(N)
                  DO I = IRSTR(N), IREND(N)
                     IJKG = IJKG + 1
                     QG(I, J, K) = DUMMY(IJKG)
                  END DO
               END DO
            END DO
         END DO
      END IF

      RETURN

c =====================================================================

      ENTRY GATHER_CHK(
     O                  QCHKG,
     I                   QCHK)

      IF (MYRANK .LT. IJNODE) THEN
         DO K = KSTR, KEND
            DO J = 1, NY
               DO I = 1, NX
                  SNDBUF(I, J, K) = QCHK(I, J, K)
               END DO
            END DO
         END DO
      END IF

      CALL MPI_GATHERV(
     &                 SNDBUF,  NGSND,         MPI_REAL8,
     &                  DUMMY,  NRCVC,  NDISP, MPI_REAL8,
     &                  IROOT, MPI_COMM_WORLD, IERR)

      IF (MYRANK .EQ. IROOT) THEN
         DO N = 0, NPROCS-1
            IJKG = NDISP(N)
            DO K = KSTR, KEND
               DO J = JRSTR(N), JREND(N)
                  DO I = IRSTR(N), IREND(N)
                     IJKG = IJKG + 1
                     QCHKG(I, J, K) = DUMMY(IJKG)
                  END DO
               END DO
            END DO
         END DO
      END IF

      RETURN

c =====================================================================

      ENTRY SCATTER_3D(
     O                      Q,
     I                     QG)

      IF (MYRANK .EQ. IROOT) THEN
         DO N = 0, NPROCS-1
            IJKG = NDISP(N)
            DO K = KSTR, KEND
               DO J = JRSTR(N), JREND(N)
                  DO I = IRSTR(N), IREND(N)
                     IJKG = IJKG + 1
                     DUMMY(IJKG) = QG(I, J, K)
                  END DO
               END DO
            END DO
         END DO
      END IF

      CALL MPI_SCATTERV(
     &                   DUMMY,  NSNDC,  NDISP, MPI_REAL8,
     &                  RCVBUF,  NSRCV,         MPI_REAL8,
     &                   IROOT, MPI_COMM_WORLD, IERR)

      IF (MYRANK .LT. IJNODE) THEN
         DO K = KSTR, KEND
            DO J = 1, NY
               DO I = 1, NX
                  Q(I+ISTR-1, J+JSTR-1, K) = RCVBUF(I, J, K)
               END DO
            END DO
         END DO
      END IF

      RETURN

c =====================================================================

      ENTRY SCATTER_BDY(
     O                    QCHK,
     I                   QCHKG)

      IF (MYRANK .EQ. IROOT) THEN
         DO N = 0, NPROCS-1
            IJKG = NDISP(N)
            DO K = KSTR, KEND
               DO J = JRSTR(N), JREND(N)
                  DO I = IRSTR(N), IREND(N)
                     IJKG = IJKG + 1
                     DUMMY(IJKG) = QCHKG(I, J, K)
                  END DO
               END DO
            END DO
         END DO
      END IF

      CALL MPI_SCATTERV(
     &                   DUMMY,  NSNDC,  NDISP, MPI_REAL8,
     &                  RCVBUF,  NSRCV,         MPI_REAL8,
     &                   IROOT, MPI_COMM_WORLD, IERR)

      IF (MYRANK .LT. IJNODE) THEN
         DO K = KSTR, KEND
            DO J = 1, NY
               DO I = 1, NX
                  QCHK(I, J, K) = RCVBUF(I, J, K)
               END DO
            END DO
         END DO
      END IF

      RETURN
      END
#else
c *********************************************************************
      SUBROUTINE GS3DST
      RETURN
      END
#endif
