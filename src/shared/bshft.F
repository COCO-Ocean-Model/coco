#ifdef OPT_PARALLEL
      SUBROUTINE SHFTST
      use zocdim
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  One-on-one communication routines for real*8 variables.
c
c  HISTORY
c     '99.10.06  H.Hasumi
c     '02.06.02  H.Hasumi
c     '07.04.23  H.Hasumi
c     '10.04.14  M.Kurogi: (COCO4.4 tripolar code by Dr. Suzuki)
c
c ---------------------------------------------------------------------

#include "mpif.h"

      INTEGER   IDIM,   JDIM,   KDIM
      REAL*8      Q1(IDIM, JDIM, KDIM)
      REAL*8      Q2(IDIM, JDIM, KDIM)
      REAL*8      Q3(IDIM, JDIM, KDIM)

      REAL*8  SDBFX1(ICOMM, NY, NZTDIM+NZDIM)
      REAL*8  SDBFX2(ICOMM, NY, NZTDIM+NZDIM)
      REAL*8  SDBFY1(NXDIM, JCOMM, NZTDIM+NZDIM)
      REAL*8  SDBFY2(NXDIM, JCOMM, NZTDIM+NZDIM)
      REAL*8  RVBFX1(ICOMM, NY, NZTDIM+NZDIM)
      REAL*8  RVBFX2(ICOMM, NY, NZTDIM+NZDIM)
      REAL*8  RVBFY1(NXDIM, JCOMM, NZTDIM+NZDIM)
      REAL*8  RVBFY2(NXDIM, JCOMM, NZTDIM+NZDIM)
      INTEGER      I,      J,      K,      N
      INTEGER NBFDIM
#ifdef OPT_TRIPOLE
      REAL*8  SDBFX1N(ICOMM, JCOMM+1, NZTDIM+NZDIM)
      REAL*8  SDBFX2N(ICOMM, JCOMM+1, NZTDIM+NZDIM)
      REAL*8  RVBFX1N(ICOMM, JCOMM+1, NZTDIM+NZDIM)
      REAL*8  RVBFX2N(ICOMM, JCOMM+1, NZTDIM+NZDIM)

      REAL*8  SDBFN1(NXDIM, 0:JCOMM, NZTDIM+NZDIM)
      REAL*8  SDBFN2(NXDIM, JCOMM,   NZTDIM+NZDIM)
      REAL*8  RVBFN1(NXDIM, JCOMM,   NZTDIM+NZDIM)
      REAL*8  RVBFN2(NXDIM, 0:JCOMM, NZTDIM+NZDIM)

      REAL*8  FACT        ! set 1 for scalar, set -1 for vector
      INTEGER NBFDM0
      INTEGER IOFF, JOFF  !(IOFF,JOFF) = (0,0) for T-cell, (-1, -1) for V-cell 
      INTEGER ISTV
#endif
c =====================================================================

      ENTRY SHIFT1(
     M                 Q1,
#ifndef OPT_TRIPOLE
     D               IDIM,   JDIM,   KDIM )
#else
     D               IDIM,   JDIM,   KDIM,
     I               FACT,   IOFF,   JOFF )
#endif

      IF (IDOWN .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, NY
               DO I = 1, ICOMM
                  SDBFX1(I, J, K) = Q1(I+ISTR-1, J+JSTR-1, K)
               END DO
            END DO
         END DO
      END IF
      IF (IUP .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, NY
               DO I = 1, ICOMM
                  SDBFX2(I, J, K) = Q1(I+IEND-ICOMM, J+JSTR-1, K)
               END DO
            END DO
         END DO
      END IF
      NBFDIM = KDIM * NY * ICOMM
      CALL SHIFTX(
     O            RVBFX1, RVBFX2,
     I            SDBFX1, SDBFX2,
     D            NBFDIM)
      IF (IDOWN .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, NY
               DO I = 1, ICOMM
                  Q1(I, J+JSTR-1, K) = RVBFX1(I, J, K)
               END DO
            END DO
         END DO
      END IF
      IF (IUP .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, NY
               DO I = 1, ICOMM
                  Q1(IEND+I, J+JSTR-1, K) = RVBFX2(I, J, K)
               END DO
            END DO
         END DO
      END IF

      IF (JDOWN .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM
               DO I = 1, NXDIM
                  SDBFY1(I, J, K) = Q1(I, J+JSTR-1, K)
               END DO
            END DO
         END DO
      END IF
      IF (JUP .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM
               DO I = 1, NXDIM
                  SDBFY2(I, J, K) = Q1(I, J+JEND-JCOMM, K)
               END DO
            END DO
         END DO
      END IF
      NBFDIM = KDIM * NXDIM * JCOMM
      CALL SHIFTY(
     O            RVBFY1, RVBFY2,
     I            SDBFY1, SDBFY2,
     D            NBFDIM)
      IF (JDOWN .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM
               DO I = 1, NXDIM
                  Q1(I, J, K) = RVBFY1(I, J, K)
               END DO
            END DO
         END DO
      END IF
      IF (JUP .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM
               DO I = 1, NXDIM
                  Q1(I, JEND+J, K) = RVBFY2(I, J, K)
               END DO
            END DO
         END DO
      END IF
#ifdef OPT_TRIPOLE
      IF (JOFF. EQ. -1) THEN
C ------
      IF (JUPE .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 0, JCOMM
               DO I = 1+IABS(IOFF), NXDIM-IABS(IOFF)
         SDBFN1(I, J, K) = FACT*Q1(NXDIM-I+1+IOFF, JEND-J, K)
               END DO
            END DO
         END DO
      END IF
      IF (JUPW .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM
               DO I = 1+IABS(IOFF), NXDIM-IABS(IOFF)
         SDBFN2(I, J, K) = FACT*Q1(NXDIM-I+1+IOFF, JEND-J, K)
               END DO
            END DO
         END DO
      END IF
      NBFDIM = KDIM * NXDIM * JCOMM
      NBFDM0 = KDIM * NXDIM * (JCOMM+1)
      CALL SHIFTNV(
     O            RVBFN1, RVBFN2,
     I            SDBFN1, SDBFN2,
     D            NBFDIM, NBFDM0)
      IF (JUPE .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM
               DO I = 1, NXDIM
                  Q1(I, JEND+J, K) = RVBFN1(I, J, K)
               END DO
            END DO
         END DO
      END IF

      ISTV=1
      IF(INODES .EQ. 1) ISTV=NXDIM/2+1
      IF (JUPW .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 0, JCOMM
C               DO I = 1, NXDIM
               DO I = ISTV, NXDIM
                  Q1(I, JEND+J, K) = RVBFN2(I, J, K)
               END DO
            END DO
         END DO
      END IF
C-----------
      ELSE
      IF (JUPE .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM
               DO I = 1+IABS(IOFF), NXDIM-IABS(IOFF)
         SDBFY1(I, J, K) = FACT*Q1(NXDIM-I+1+IOFF,JEND-J+1+JOFF,K)
               END DO
            END DO
         END DO
      END IF
      IF (JUPW .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM
               DO I = 1+IABS(IOFF), NXDIM-IABS(IOFF)
         SDBFY2(I, J, K) = FACT*Q1(NXDIM-I+1+IOFF,JEND-J+1+JOFF,K)
               END DO
            END DO
         END DO
      END IF
      NBFDIM = KDIM * NXDIM * JCOMM
      CALL SHIFTYN(
     O            RVBFY1, RVBFY2,
     I            SDBFY1, SDBFY2,
     D            NBFDIM)
      IF (JUPE .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM
               DO I = 1, NXDIM
                  Q1(I, JEND+J, K) = RVBFY1(I, J, K)
               END DO
            END DO
         END DO
      END IF
      IF (JUPW .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM
               DO I = 1, NXDIM
                  Q1(I, JEND+J, K) = RVBFY2(I, J, K)
               END DO
            END DO
         END DO
      END IF
      ENDIF
CC
CC
      IF(IOFF.NE.0 .and. JRANK.EQ.JNODES-1) THEN
      IF (IDOWN .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM+1
               DO I = 1, ICOMM
            SDBFX1N(I, J, K) = Q1(I+ISTR-1, J+JEND-1, K)
               END DO
            END DO
         END DO
      END IF
      IF (IUP .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM+1
               DO I = 1, ICOMM
            SDBFX2N(I, J, K) = Q1(I+IEND-ICOMM, J+JEND-1, K)
               END DO
            END DO
         END DO
      END IF
      NBFDIM = KDIM * (JCOMM+1) * ICOMM
      CALL SHIFTX(
     O            RVBFX1N, RVBFX2N,
     I            SDBFX1N, SDBFX2N,
     D            NBFDIM)
      IF (IDOWN .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM+1
               DO I = 1, ICOMM
            Q1(I, J+JEND-1, K) = RVBFX1N(I, J, K)
               END DO
            END DO
         END DO
      END IF
      IF (IUP .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM+1
               DO I = 1, ICOMM
            Q1(IEND+I, J+JEND-1, K) = RVBFX2N(I, J, K)
               END DO
            END DO
         END DO
      END IF
      ENDIF
#endif
      RETURN

c =====================================================================

      ENTRY SHIFT2(
     M                 Q1,     Q2,
#ifndef OPT_TRIPOLE
     D               IDIM,   JDIM,   KDIM )
#else
     D               IDIM,   JDIM,   KDIM,
     I               FACT,   IOFF,   JOFF )
#endif

      IF (IDOWN .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, NY
               DO I = 1, ICOMM
                  SDBFX1(I, J, K) = Q1(I+ISTR-1, J+JSTR-1, K)
                  SDBFX1(I, J, K+KDIM) = Q2(I+ISTR-1, J+JSTR-1, K)
               END DO
            END DO
         END DO
      END IF
      IF (IUP .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, NY
               DO I = 1, ICOMM
                  SDBFX2(I, J, K) = Q1(I+IEND-ICOMM, J+JSTR-1, K)
                  SDBFX2(I, J, K+KDIM) = Q2(I+IEND-ICOMM, J+JSTR-1, K)
               END DO
            END DO
         END DO
      END IF
      NBFDIM = 2 * KDIM * NY * ICOMM
      CALL SHIFTX(
     O            RVBFX1, RVBFX2,
     I            SDBFX1, SDBFX2,
     D            NBFDIM)
      IF (IDOWN .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, NY
               DO I = 1, ICOMM
                  Q1(I, J+JSTR-1, K) = RVBFX1(I, J, K)
                  Q2(I, J+JSTR-1, K) = RVBFX1(I, J, K+KDIM)
               END DO
            END DO
         END DO
      END IF
      IF (IUP .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, NY
               DO I = 1, ICOMM
                  Q1(IEND+I, J+JSTR-1, K) = RVBFX2(I, J, K)
                  Q2(IEND+I, J+JSTR-1, K) = RVBFX2(I, J, K+KDIM)
               END DO
            END DO
         END DO
      END IF

      IF (JDOWN .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM
               DO I = 1, NXDIM
                  SDBFY1(I, J, K) = Q1(I, J+JSTR-1, K)
                  SDBFY1(I, J, K+KDIM) = Q2(I, J+JSTR-1, K)
               END DO
            END DO
         END DO
      END IF
      IF (JUP .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM
               DO I = 1, NXDIM
                  SDBFY2(I, J, K) = Q1(I, J+JEND-JCOMM, K)
                  SDBFY2(I, J, K+KDIM) = Q2(I, J+JEND-JCOMM, K)
               END DO
            END DO
         END DO
      END IF
      NBFDIM = 2 * KDIM * NXDIM * JCOMM
      CALL SHIFTY(
     O            RVBFY1, RVBFY2,
     I            SDBFY1, SDBFY2,
     D            NBFDIM)
      IF (JDOWN .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM
               DO I = 1, NXDIM
                  Q1(I, J, K) = RVBFY1(I, J, K)
                  Q2(I, J, K) = RVBFY1(I, J, K+KDIM)
               END DO
            END DO
         END DO
      END IF
      IF (JUP .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM
               DO I = 1, NXDIM
                  Q1(I, JEND+J, K) = RVBFY2(I, J, K)
                  Q2(I, JEND+J, K) = RVBFY2(I, J, K+KDIM)
               END DO
            END DO
         END DO
      END IF
#ifdef OPT_TRIPOLE
      IF (JOFF.EQ.-1) THEN
C---
      IF (JUPE .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 0, JCOMM
               DO I = 1+IABS(IOFF), NXDIM-IABS(IOFF)
         SDBFN1(I, J, K) = FACT*Q1(NXDIM-I+1+IOFF, JEND-J, K)
         SDBFN1(I, J, K+KDIM) = FACT*Q2(NXDIM-I+1+IOFF, JEND-J, K)
               END DO
            END DO
         END DO
      END IF
      IF (JUPW .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM
               DO I = 1+IABS(IOFF), NXDIM-IABS(IOFF)
         SDBFN2(I, J, K) = FACT*Q1(NXDIM-I+1+IOFF, JEND-J, K)
         SDBFN2(I, J, K+KDIM) = FACT*Q2(NXDIM-I+1+IOFF, JEND-J, K)
               END DO
            END DO
         END DO
      END IF
      NBFDIM = 2 * KDIM * NXDIM * JCOMM
      NBFDM0 = 2 * KDIM * NXDIM * (JCOMM+1)
      CALL SHIFTNV(
     O            RVBFN1, RVBFN2,
     I            SDBFN1, SDBFN2,
     D            NBFDIM, NBFDM0)
      IF (JUPE .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM
               DO I = 1, NXDIM
                  Q1(I, JEND+J, K) = RVBFN1(I, J, K)
                  Q2(I, JEND+J, K) = RVBFN1(I, J, K+KDIM)
               END DO
            END DO
         END DO
      END IF

      ISTV=1
      IF(INODES .EQ. 1) ISTV=NXDIM/2+1
      IF (JUPW .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 0, JCOMM
C               DO I = 1, NXDIM
               DO I = ISTV, NXDIM
                  Q1(I, JEND+J, K) = RVBFN2(I, J, K)
                  Q2(I, JEND+J, K) = RVBFN2(I, J, K+KDIM)
               END DO
            END DO
         END DO
      END IF
C---
      ELSE
      IF (JUPE .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM
               DO I = 1+IABS(IOFF), NXDIM-IABS(IOFF)
         SDBFY1(I, J, K) = FACT*Q1(NXDIM-I+1+IOFF,JEND-J+1+JOFF,K)
         SDBFY1(I, J, K+KDIM) = FACT*Q2(NXDIM-I+1+IOFF,JEND-J+1+JOFF,K)
               END DO
            END DO
         END DO
      END IF
      IF (JUPW .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM
               DO I = 1+IABS(IOFF), NXDIM-IABS(IOFF)
         SDBFY2(I, J, K) = FACT*Q1(NXDIM-I+1+IOFF,JEND-J+1+JOFF,K)
         SDBFY2(I, J, K+KDIM) = FACT*Q2(NXDIM-I+1+IOFF,JEND-J+1+JOFF,K)
               END DO
            END DO
         END DO
      END IF
      NBFDIM = 2 * KDIM * NXDIM * JCOMM
      CALL SHIFTYN(
     O            RVBFY1, RVBFY2,
     I            SDBFY1, SDBFY2,
     D            NBFDIM)
      IF (JUPE .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM
               DO I = 1, NXDIM
                  Q1(I, JEND+J, K) = RVBFY1(I, J, K)
                  Q2(I, JEND+J, K) = RVBFY1(I, J, K+KDIM)
               END DO
            END DO
         END DO
      END IF
      IF (JUPW .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM
               DO I = 1, NXDIM
                  Q1(I, JEND+J, K) = RVBFY2(I, J, K)
                  Q2(I, JEND+J, K) = RVBFY2(I, J, K+KDIM)
               END DO
            END DO
         END DO
      END IF
      ENDIF
CC
CC
      IF(IOFF.NE.0 .and. JRANK.EQ.JNODES-1) THEN
      IF (IDOWN .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM+1
               DO I = 1, ICOMM
            SDBFX1N(I, J, K) = Q1(I+ISTR-1, J+JEND-1, K)
            SDBFX1N(I, J, K+KDIM) = Q2(I+ISTR-1, J+JEND-1, K)
               END DO
            END DO
         END DO
      END IF
      IF (IUP .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM+1
               DO I = 1, ICOMM
            SDBFX2N(I, J, K) = Q1(I+IEND-ICOMM, J+JEND-1, K)
            SDBFX2N(I, J, K+KDIM) = Q2(I+IEND-ICOMM, J+JEND-1, K)
               END DO
            END DO
         END DO
      END IF
      NBFDIM = 2 * KDIM * (JCOMM+1) * ICOMM
      CALL SHIFTX(
     O            RVBFX1N, RVBFX2N,
     I            SDBFX1N, SDBFX2N,
     D            NBFDIM)
      IF (IDOWN .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM+1
               DO I = 1, ICOMM
            Q1(I, J+JEND-1, K) = RVBFX1N(I, J, K)
            Q2(I, J+JEND-1, K) = RVBFX1N(I, J, K+KDIM)
               END DO
            END DO
         END DO
      END IF
      IF (IUP .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM+1
               DO I = 1, ICOMM
            Q1(IEND+I, J+JEND-1, K) = RVBFX2N(I, J, K)
            Q2(IEND+I, J+JEND-1, K) = RVBFX2N(I, J, K+KDIM)
               END DO
            END DO
         END DO
      END IF
      ENDIF
#endif
      RETURN

c =====================================================================

      ENTRY SHIFT3(
     M                 Q1,     Q2,     Q3,
#ifndef OPT_TRIPOLE
     D               IDIM,   JDIM,   KDIM )
#else
     D               IDIM,   JDIM,   KDIM,
     I               FACT,   IOFF,   JOFF )
#endif

      IF (IDOWN .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, NY
               DO I = 1, ICOMM
                  SDBFX1(I, J, K) = Q1(I+ISTR-1, J+JSTR-1, K)
                  SDBFX1(I, J, K+KDIM) = Q2(I+ISTR-1, J+JSTR-1, K)
                  SDBFX1(I, J, K+KDIM*2) = Q3(I+ISTR-1, J+JSTR-1, K)
               END DO
            END DO
         END DO
      END IF
      IF (IUP .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, NY
               DO I = 1, ICOMM
                  SDBFX2(I, J, K) = Q1(I+IEND-ICOMM, J+JSTR-1, K)
                  SDBFX2(I, J, K+KDIM) = Q2(I+IEND-ICOMM, J+JSTR-1, K)
                  SDBFX2(I, J, K+KDIM*2) = Q3(I+IEND-ICOMM, J+JSTR-1, K)
               END DO
            END DO
         END DO
      END IF
      NBFDIM = 3 * KDIM * NY * ICOMM
      CALL SHIFTX(
     O            RVBFX1, RVBFX2,
     I            SDBFX1, SDBFX2,
     D            NBFDIM)
      IF (IDOWN .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, NY
               DO I = 1, ICOMM
                  Q1(I, J+JSTR-1, K) = RVBFX1(I, J, K)
                  Q2(I, J+JSTR-1, K) = RVBFX1(I, J, K+KDIM)
                  Q3(I, J+JSTR-1, K) = RVBFX1(I, J, K+KDIM*2)
               END DO
            END DO
         END DO
      END IF
      IF (IUP .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, NY
               DO I = 1, ICOMM
                  Q1(IEND+I, J+JSTR-1, K) = RVBFX2(I, J, K)
                  Q2(IEND+I, J+JSTR-1, K) = RVBFX2(I, J, K+KDIM)
                  Q3(IEND+I, J+JSTR-1, K) = RVBFX2(I, J, K+KDIM*2)
               END DO
            END DO
         END DO
      END IF

      IF (JDOWN .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM
               DO I = 1, NXDIM
                  SDBFY1(I, J, K) = Q1(I, J+JSTR-1, K)
                  SDBFY1(I, J, K+KDIM) = Q2(I, J+JSTR-1, K)
                  SDBFY1(I, J, K+KDIM*2) = Q3(I, J+JSTR-1, K)
               END DO
            END DO
         END DO
      END IF
      IF (JUP .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM
               DO I = 1, NXDIM
                  SDBFY2(I, J, K) = Q1(I, J+JEND-JCOMM, K)
                  SDBFY2(I, J, K+KDIM) = Q2(I, J+JEND-JCOMM, K)
                  SDBFY2(I, J, K+KDIM*2) = Q3(I, J+JEND-JCOMM, K)
               END DO
            END DO
         END DO
      END IF
      NBFDIM = 3 * KDIM * NXDIM * JCOMM
      CALL SHIFTY(
     O            RVBFY1, RVBFY2,
     I            SDBFY1, SDBFY2,
     D            NBFDIM)
      IF (JDOWN .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM
               DO I = 1, NXDIM
                  Q1(I, J, K) = RVBFY1(I, J, K)
                  Q2(I, J, K) = RVBFY1(I, J, K+KDIM)
                  Q3(I, J, K) = RVBFY1(I, J, K+KDIM*2)
               END DO
            END DO
         END DO
      END IF
      IF (JUP .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM
               DO I = 1, NXDIM
                  Q1(I, JEND+J, K) = RVBFY2(I, J, K)
                  Q2(I, JEND+J, K) = RVBFY2(I, J, K+KDIM)
                  Q3(I, JEND+J, K) = RVBFY2(I, J, K+KDIM*2)
               END DO
            END DO
         END DO
      END IF
#ifdef OPT_TRIPOLE
      IF (JOFF .EQ. -1) THEN
C--
      IF (JUPE .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 0, JCOMM
               DO I = 1+IABS(IOFF), NXDIM-IABS(IOFF)
         SDBFN1(I, J, K) = FACT*Q1(NXDIM-I+1+IOFF, JEND-J, K)
         SDBFN1(I, J, K+KDIM) = FACT*Q2(NXDIM-I+1+IOFF, JEND-J, K)
         SDBFN1(I, J, K+2*KDIM) = FACT*Q3(NXDIM-I+1+IOFF, JEND-J, K)
                END DO
            END DO
         END DO
      END IF
      IF (JUPW .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM
               DO I = 1+IABS(IOFF), NXDIM-IABS(IOFF)
         SDBFN2(I, J, K) = FACT*Q1(NXDIM-I+1+IOFF, JEND-J, K)
         SDBFN2(I, J, K+KDIM) = FACT*Q2(NXDIM-I+1+IOFF, JEND-J, K)
         SDBFN2(I, J, K+2*KDIM) = FACT*Q3(NXDIM-I+1+IOFF, JEND-J, K)
               END DO
            END DO
         END DO
      END IF
      NBFDIM = 3 * KDIM * NXDIM * JCOMM
      NBFDM0 = 3 * KDIM * NXDIM * (JCOMM+1)
      CALL SHIFTNV(
     O            RVBFN1, RVBFN2,
     I            SDBFN1, SDBFN2,
     D            NBFDIM, NBFDM0)
      IF (JUPE .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM
               DO I = 1, NXDIM
                  Q1(I, JEND+J, K) = RVBFN1(I, J, K)
                  Q2(I, JEND+J, K) = RVBFN1(I, J, K+KDIM)
                  Q3(I, JEND+J, K) = RVBFN1(I, J, K+2*KDIM)
               END DO
            END DO
         END DO
      END IF

      ISTV=1
      IF(INODES .EQ. 1) ISTV=NXDIM/2+1
      IF (JUPW .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 0, JCOMM
C               DO I = 1, NXDIM
               DO I = ISTV, NXDIM
                  Q1(I, JEND+J, K) = RVBFN2(I, J, K)
                  Q2(I, JEND+J, K) = RVBFN2(I, J, K+KDIM)
                  Q3(I, JEND+J, K) = RVBFN2(I, J, K+2*KDIM)
               END DO
            END DO
         END DO
      END IF
C--
      ELSE
      IF (JUPE .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM
               DO I = 1+IABS(IOFF), NXDIM-IABS(IOFF)
         SDBFY1(I,J,K) = FACT*Q1(NXDIM-I+1+IOFF,JEND-J+1+JOFF,K)
         SDBFY1(I,J,K+KDIM) = FACT*Q2(NXDIM-I+1+IOFF,JEND-J+1+JOFF,K)
         SDBFY1(I,J,K+2*KDIM) = FACT*Q3(NXDIM-I+1+IOFF,JEND-J+1+JOFF,K)
               END DO
            END DO
         END DO
      END IF
      IF (JUPW .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM
               DO I = 1+IABS(IOFF), NXDIM-IABS(IOFF)
         SDBFY2(I,J,K) = FACT*Q1(NXDIM-I+1+IOFF,JEND-J+1+JOFF,K)
         SDBFY2(I,J,K+KDIM) = FACT*Q2(NXDIM-I+1+IOFF,JEND-j+1+JOFF,K)
         SDBFY2(I,J,K+2*KDIM) = FACT*Q3(NXDIM-I+1+IOFF,JEND-j+1+JOFF,K)
               END DO
            END DO
         END DO
      END IF
      NBFDIM = 3 * KDIM * NXDIM * JCOMM
      CALL SHIFTYN(
     O            RVBFY1, RVBFY2,
     I            SDBFY1, SDBFY2,
     D            NBFDIM)
      IF (JUPE .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM
               DO I = 1, NXDIM
                  Q1(I, JEND+J, K) = RVBFY1(I, J, K)
                  Q2(I, JEND+J, K) = RVBFY1(I, J, K+KDIM)
                  Q3(I, JEND+J, K) = RVBFY1(I, J, K+2*KDIM)
               END DO
            END DO
         END DO
      END IF
      IF (JUPW .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM
               DO I = 1, NXDIM
                  Q1(I, JEND+J, K) = RVBFY2(I, J, K)
                  Q2(I, JEND+J, K) = RVBFY2(I, J, K+KDIM)
                  Q3(I, JEND+J, K) = RVBFY2(I, J, K+2*KDIM)
               END DO
            END DO
         END DO
      END IF
      ENDIF
CC
CC
      IF( IOFF.NE.0 .and. JRANK.EQ.JNODES-1 ) THEN
      IF (IDOWN .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM+1
               DO I = 1, ICOMM
            SDBFX1N(I, J, K) = Q1(I+ISTR-1, J+JEND-1, K)
            SDBFX1N(I, J, K+KDIM) = Q2(I+ISTR-1, J+JEND-1, K)
            SDBFX1N(I, J, K+KDIM*2) = Q3(I+ISTR-1, J+JEND-1, K)
               END DO
            END DO
         END DO
      END IF
      IF (IUP .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM+1
               DO I = 1, ICOMM
            SDBFX2N(I, J, K) = Q1(I+IEND-ICOMM, J+JEND-1, K)
            SDBFX2N(I, J, K+KDIM) = Q2(I+IEND-ICOMM, J+JEND-1, K)
            SDBFX2N(I, J, K+KDIM*2) = Q3(I+IEND-ICOMM, J+JEND-1, K)
               END DO
            END DO
         END DO
      END IF
      NBFDIM = 3 * KDIM * (JCOMM+1) * ICOMM
      CALL SHIFTX(
     O            RVBFX1N, RVBFX2N,
     I            SDBFX1N, SDBFX2N,
     D            NBFDIM)
      IF (IDOWN .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM+1
               DO I = 1, ICOMM
            Q1(I, J+JEND-1, K) = RVBFX1N(I, J, K)
            Q2(I, J+JEND-1, K) = RVBFX1N(I, J, K+KDIM)
            Q3(I, J+JEND-1, K) = RVBFX1N(I, J, K+KDIM*2)
               END DO
            END DO
         END DO
      END IF
      IF (IUP .NE. MPI_PROC_NULL) THEN
         DO K = 1, KDIM
            DO J = 1, JCOMM+1
               DO I = 1, ICOMM
            Q1(IEND+I, J+JEND-1, K) = RVBFX2N(I, J, K)
            Q2(IEND+I, J+JEND-1, K) = RVBFX2N(I, J, K+KDIM)
            Q3(IEND+I, J+JEND-1, K) = RVBFX2N(I, J, K+KDIM*2)
               END DO
            END DO
         END DO
      END IF
      ENDIF
#endif
      RETURN
      END

c *********************************************************************

      SUBROUTINE SHIFTX(
     O            RVBFX1, RVBFX2,
     I            SDBFX1, SDBFX2,
     D            NBFDIM)
      use zocdim
      IMPLICIT NONE

#include "mpif.h"

      REAL*8  SDBFX1(ICOMM, NY, NZTDIM+NZDIM)
      REAL*8  SDBFX2(ICOMM, NY, NZTDIM+NZDIM)
      REAL*8  RVBFX1(ICOMM, NY, NZTDIM+NZDIM)
      REAL*8  RVBFX2(ICOMM, NY, NZTDIM+NZDIM)
      INTEGER NBFDIM

      INTEGER ISRQX1, ISRQX2
      INTEGER IRRQX1, IRRQX2
      INTEGER ISTMPI(MPI_STATUS_SIZE)

      CALL MPI_ISEND(
     &               SDBFX1, NBFDIM, MPI_REAL8,
     &                IDOWN,      1, MPI_COMM_WORLD,
     &               ISRQX1,   IERR)
      CALL MPI_ISEND(
     &               SDBFX2, NBFDIM, MPI_REAL8,
     &                  IUP,      2, MPI_COMM_WORLD,
     &               ISRQX2,   IERR)
      CALL MPI_IRECV(
     &               RVBFX2, NBFDIM, MPI_REAL8,
     &                  IUP,      1, MPI_COMM_WORLD,
     &               IRRQX1,   IERR)
      CALL MPI_IRECV(
     &               RVBFX1, NBFDIM, MPI_REAL8,
     &                IDOWN,      2, MPI_COMM_WORLD,
     &               IRRQX2,   IERR)
      CALL MPI_WAIT(ISRQX1, ISTMPI,   IERR)
      CALL MPI_WAIT(ISRQX2, ISTMPI,   IERR)
      CALL MPI_WAIT(IRRQX1, ISTMPI,   IERR)
      CALL MPI_WAIT(IRRQX2, ISTMPI,   IERR)

      RETURN
      END

c *********************************************************************

      SUBROUTINE SHIFTY(
     O            RVBFY1, RVBFY2,
     I            SDBFY1, SDBFY2,
     D            NBFDIM)
      use zocdim
      IMPLICIT NONE

#include "mpif.h"

      REAL*8  SDBFY1(NXDIM, JCOMM, NZTDIM+NZDIM)
      REAL*8  SDBFY2(NXDIM, JCOMM, NZTDIM+NZDIM)
      REAL*8  RVBFY1(NXDIM, JCOMM, NZTDIM+NZDIM)
      REAL*8  RVBFY2(NXDIM, JCOMM, NZTDIM+NZDIM)
      INTEGER NBFDIM

      INTEGER ISRQY1, ISRQY2
      INTEGER IRRQY1, IRRQY2
      INTEGER ISTMPI(MPI_STATUS_SIZE)

      CALL MPI_ISEND(
     &               SDBFY1, NBFDIM, MPI_REAL8,
     &                JDOWN,      3, MPI_COMM_WORLD,
     &               ISRQY1,   IERR)
      CALL MPI_ISEND(
     &               SDBFY2, NBFDIM, MPI_REAL8,
     &                  JUP,      4, MPI_COMM_WORLD,
     &               ISRQY2,   IERR)
      CALL MPI_IRECV(
     &               RVBFY2, NBFDIM, MPI_REAL8,
     &                  JUP,      3, MPI_COMM_WORLD,
     &               IRRQY1,   IERR)
      CALL MPI_IRECV(
     &               RVBFY1, NBFDIM, MPI_REAL8,
     &                JDOWN,      4, MPI_COMM_WORLD,
     &               IRRQY2,   IERR)
      CALL MPI_WAIT(ISRQY1, ISTMPI,   IERR)
      CALL MPI_WAIT(ISRQY2, ISTMPI,   IERR)
      CALL MPI_WAIT(IRRQY1, ISTMPI,   IERR)
      CALL MPI_WAIT(IRRQY2, ISTMPI,   IERR)

      RETURN
      END
#ifdef OPT_TRIPOLE
      SUBROUTINE SHIFTYN(
     O            RVBFY1, RVBFY2,
     I            SDBFY1, SDBFY2,
     D            NBFDIM)
      use zocdim
      IMPLICIT NONE

#include "mpif.h"

      REAL*8  SDBFY1(NXDIM, JCOMM, NZTDIM+NZDIM)
      REAL*8  SDBFY2(NXDIM, JCOMM, NZTDIM+NZDIM)
      REAL*8  RVBFY1(NXDIM, JCOMM, NZTDIM+NZDIM)
      REAL*8  RVBFY2(NXDIM, JCOMM, NZTDIM+NZDIM)
      INTEGER NBFDIM

      INTEGER ISRQY1, ISRQY2
      INTEGER IRRQY1, IRRQY2
      INTEGER ISTMPI(MPI_STATUS_SIZE)

      CALL MPI_ISEND(
     &               SDBFY1, NBFDIM, MPI_REAL8,
     &                JUPE,       3, MPI_COMM_WORLD,
     &               ISRQY1,   IERR)
      CALL MPI_ISEND(
     &               SDBFY2, NBFDIM, MPI_REAL8,
     &                 JUPW,      4, MPI_COMM_WORLD,
     &               ISRQY2,   IERR)
      CALL MPI_IRECV(
     &               RVBFY2, NBFDIM, MPI_REAL8,
     &                 JUPW,      3, MPI_COMM_WORLD,
     &               IRRQY1,   IERR)
      CALL MPI_IRECV(
     &               RVBFY1, NBFDIM, MPI_REAL8,
     &                 JUPE,      4, MPI_COMM_WORLD,
     &               IRRQY2,   IERR)
      CALL MPI_WAIT(ISRQY1, ISTMPI,   IERR)
      CALL MPI_WAIT(ISRQY2, ISTMPI,   IERR)
      CALL MPI_WAIT(IRRQY1, ISTMPI,   IERR)
      CALL MPI_WAIT(IRRQY2, ISTMPI,   IERR)

      RETURN
      END

      SUBROUTINE SHIFTNV(
     O            RVBFN1, RVBFN2,
     I            SDBFN1, SDBFN2,
     D            NBFDIM, NBFDM0)
      use zocdim
      IMPLICIT NONE

#include "mpif.h"

      REAL*8  SDBFN1(NXDIM, JCOMM+1, NZTDIM+NZDIM)
      REAL*8  SDBFN2(NXDIM, JCOMM,   NZTDIM+NZDIM)
      REAL*8  RVBFN1(NXDIM, JCOMM,   NZTDIM+NZDIM)
      REAL*8  RVBFN2(NXDIM, JCOMM+1, NZTDIM+NZDIM)
      INTEGER NBFDIM, NBFDM0

      INTEGER ISRQY1, ISRQY2
      INTEGER IRRQY1, IRRQY2
      INTEGER ISTMPI(MPI_STATUS_SIZE)

      CALL MPI_ISEND(
     &               SDBFN1, NBFDM0, MPI_REAL8,
     &                JUPE,       3, MPI_COMM_WORLD,
     &               ISRQY1,   IERR)
      CALL MPI_ISEND(
     &               SDBFN2, NBFDIM, MPI_REAL8,
     &                 JUPW,      4, MPI_COMM_WORLD,
     &               ISRQY2,   IERR)
      CALL MPI_IRECV(
     &               RVBFN2, NBFDM0, MPI_REAL8,
     &                 JUPW,      3, MPI_COMM_WORLD,
     &               IRRQY1,   IERR)
      CALL MPI_IRECV(
     &               RVBFN1, NBFDIM, MPI_REAL8,
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
      SUBROUTINE SHFTST
      RETURN
      END
#endif
