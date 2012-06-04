#ifdef OPT_PARALLEL
      SUBROUTINE GS2DST
      use zocdim
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Group communication (GATHER and SCATTER) routines for 2-dimensional
c variables.
c
c  HISTORY
c     '99.10.05  H.Hasumi
c     '07.04.23  H.Hasumi
c
c ---------------------------------------------------------------------

#include "mpif.h"

      REAL*8       Q(NXDIM, NYDIM),     QG(NXGDIM, NYGDIM)
      REAL*8      QS(NX   , NY   ),    QSG(IGSTR:IGEND, JGSTR:JGEND)
      INTEGER     IQ(NXDIM, NYDIM),    IQG(NXGDIM, NYGDIM)

      INTEGER, ALLOCATABLE :: IRSTR(:), JRSTR(:)
      INTEGER, ALLOCATABLE :: IREND(:), JREND(:)
      INTEGER, ALLOCATABLE :: NSNDC(:), NRCVC(:), NDISP(:)
      INTEGER  NGSND,  NSRCV
      SAVE IRSTR, JRSTR, IREND, JREND, NSNDC, NRCVC, NDISP,
     &     NGSND, NSRCV

      REAL*8,  ALLOCATABLE :: SNDBUF(:, :), RCVBUF(:, :)
      INTEGER, ALLOCATABLE :: IRCBUF(:, :)
      REAL*8,  ALLOCATABLE :: DUMMY(:)
      INTEGER, ALLOCATABLE :: IDUMY(:)
      SAVE SNDBUF, RCVBUF, IRCBUF, DUMMY, IDUMY

      INTEGER      I,      J,    IJG,      N

      ALLOCATE(IRSTR(0:NPROCS-1))
      ALLOCATE(JRSTR(0:NPROCS-1))
      ALLOCATE(IREND(0:NPROCS-1))
      ALLOCATE(JREND(0:NPROCS-1))
      ALLOCATE(NSNDC(0:NPROCS-1))
      ALLOCATE(NRCVC(0:NPROCS-1))
      ALLOCATE(NDISP(0:NPROCS-1))
      ALLOCATE(DUMMY(NXYG+NPROCS))
      ALLOCATE(IDUMY(NXYG+NPROCS))
      IF (MYRANK .GE. IJNODE) THEN
         ALLOCATE(SNDBUF(1, 1))
         ALLOCATE(RCVBUF(1, 1))
         ALLOCATE(IRCBUF(1, 1))
      ELSE
         ALLOCATE(SNDBUF(NX, NY))
         ALLOCATE(RCVBUF(NX, NY))
         ALLOCATE(IRCBUF(NX, NY))
      END IF

      SNDBUF = 0.D0
      RCVBUF = 0.D0
      DUMMY  = 0.D0
      IRCBUF = 0
      IDUMY  = 0

      DO J = 0, JNODES-1
         DO I = 0, INODES-1
            N = J * INODES + I
            IRSTR(N) = IGSTR + I * NX
            IREND(N) = IRSTR(N) + NX - 1
            JRSTR(N) = JGSTR + J * NY
            JREND(N) = JRSTR(N) + NY - 1
            NSNDC(N) = NXY
            NRCVC(N) = NXY
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
         NGSND = NXY
         NSRCV = NXY
      END IF

      RETURN

c =====================================================================

      ENTRY GATHER_2D(
     O                    QG,
     I                     Q)

      IF (MYRANK .LT. IJNODE) THEN
         DO J = 1, NY
            DO I = 1, NX
               SNDBUF(I, J) = Q(I+ISTR-1, J+JSTR-1)
            END DO
         END DO
      END IF

      CALL MPI_GATHERV(
     &                 SNDBUF,  NGSND,         MPI_REAL8,
     &                  DUMMY,  NRCVC,  NDISP, MPI_REAL8,
     &                  IROOT, MPI_COMM_WORLD, IERR)

      IF (MYRANK .EQ. IROOT) THEN
         DO N = 0, NPROCS-1
            IJG = NDISP(N)
            DO J = JRSTR(N), JREND(N)
               DO I = IRSTR(N), IREND(N)
                  IJG = IJG + 1
                  QG(I, J) = DUMMY(IJG)
               END DO
            END DO
         END DO
      END IF

      RETURN

c =====================================================================

      ENTRY SCATTER_2D(
     O                      Q,
     I                     QG)

      IF (MYRANK .EQ. IROOT) THEN
         DO N = 0, NPROCS-1
            IJG = NDISP(N)
            DO J = JRSTR(N), JREND(N)
               DO I = IRSTR(N), IREND(N)
                  IJG = IJG + 1
                  DUMMY(IJG) = QG(I, J)
               END DO
            END DO
         END DO
      END IF

      CALL MPI_SCATTERV(
     &                   DUMMY,  NSNDC,  NDISP, MPI_REAL8,
     &                  RCVBUF,  NSRCV,         MPI_REAL8,
     &                   IROOT, MPI_COMM_WORLD, IERR)

      IF (MYRANK .LT. IJNODE) THEN
         DO J = 1, NY
            DO I = 1, NX
               Q(I+ISTR-1, J+JSTR-1) = RCVBUF(I, J)
            END DO
         END DO
      END IF

      RETURN

c =====================================================================

      ENTRY SCATTER_2D_INT(
     O                         IQ,
     I                        IQG)

      IF (MYRANK .EQ. IROOT) THEN
         DO N = 0, NPROCS-1
            IJG = NDISP(N)
            DO J = JRSTR(N), JREND(N)
               DO I = IRSTR(N), IREND(N)
                  IJG = IJG + 1
                  IDUMY(IJG) = IQG(I, J)
               END DO
            END DO
         END DO
      END IF

      CALL MPI_SCATTERV(
     &                   IDUMY,  NSNDC,  NDISP, MPI_INTEGER,
     &                  IRCBUF,  NSRCV,         MPI_INTEGER,
     &                   IROOT, MPI_COMM_WORLD, IERR)

      IF (MYRANK .LT. IJNODE) THEN
         DO J = 1, NY
            DO I = 1, NX
               IQ(I+ISTR-1, J+JSTR-1) = IRCBUF(I, J)
            END DO
         END DO
      END IF

      RETURN

c =====================================================================

      ENTRY SCATTER_SFC(
     O                         QS,
     I                        QSG)

      IF (MYRANK .EQ. IROOT) THEN
         DO N = 0, NPROCS-1
            IJG = NDISP(N)
            DO J = JRSTR(N), JREND(N)
               DO I = IRSTR(N), IREND(N)
                  IJG = IJG + 1
                  DUMMY(IJG) = QSG(I, J)
               END DO
            END DO
         END DO
      END IF

      CALL MPI_SCATTERV(
     &                   DUMMY,  NSNDC,  NDISP, MPI_REAL8,
     &                  RCVBUF,  NSRCV,         MPI_REAL8,
     &                   IROOT, MPI_COMM_WORLD, IERR)

      IF (MYRANK .LT. IJNODE) THEN
         DO J = 1, NY
            DO I = 1, NX
               QS(I, J) = RCVBUF(I, J)
            END DO
         END DO
      END IF

      RETURN
      END
#else
c *********************************************************************
      SUBROUTINE GS2DST
      RETURN
      END
#endif
