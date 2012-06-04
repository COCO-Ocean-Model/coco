#ifdef OPT_PARALLEL
      SUBROUTINE GSIDST
      use zocdim
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Group communication (GATHER and SCATTER) routines for sea ice
c variables (multi-category thickness representation).
c
c  HISTORY
c     '07.09.25  H.Hasumi: from multi-category icedcoco3
c
c ---------------------------------------------------------------------

#include "mpif.h"

      REAL*8       Q(NXDIM , NYDIM , 0:NIC)
      REAL*8      QG(NXGDIM, NYGDIM, 0:NIC)

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
      ALLOCATE(DUMMY(NXYG*NIC+NPROCS))
      IF (MYRANK .GE. IJNODE) THEN
         ALLOCATE(SNDBUF(1, 1, 1))
         ALLOCATE(RCVBUF(1, 1, 1))
      ELSE
         ALLOCATE(SNDBUF(NX, NY, NIC))
         ALLOCATE(RCVBUF(NX, NY, NIC))
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
            NSNDC(N) = NXY*NIC
            NRCVC(N) = NXY*NIC
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
         NGSND = NXY*NIC
         NSRCV = NXY*NIC
      END IF

      RETURN

c =====================================================================

      ENTRY GATHER_ID(
     O                    QG,
     I                     Q)

      IF (MYRANK .LT. IJNODE) THEN
         DO K = 1, NIC
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
            DO K = 1, NIC
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

      ENTRY SCATTER_ID(
     O                      Q,
     I                     QG)

      IF (MYRANK .EQ. IROOT) THEN
         DO N = 0, NPROCS-1
            IJKG = NDISP(N)
            DO K = 1, NIC
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
         DO K = 1, NIC
            DO J = 1, NY
               DO I = 1, NX
                  Q(I+ISTR-1, J+JSTR-1, K) = RCVBUF(I, J, K)
               END DO
            END DO
         END DO
      END IF

      RETURN
      END
#else
c *********************************************************************
      SUBROUTINE GSIDST
      RETURN
      END
#endif
