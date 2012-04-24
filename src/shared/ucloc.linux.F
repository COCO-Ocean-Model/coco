      SUBROUTINE CLCOUT

c --- information -----------------------------------------------------
c
c  Monitoring CPU time on Sun OS (Sun FORTRAN, using the intrinsic
c function "ETIME")
c
c  HISTORY
c     '97.03.18  H.Hasumi: From AGCM5.4 developed by A.Numaguti
c     '07.04.23  H.Hasumi
c
c ---------------------------------------------------------------------

      INTEGER   NCLMAX
      PARAMETER(NCLMAX = 100)

      CHARACTER HTTL*(*)

      INTEGER IC
      REAL*8  CPUT, VPUT
      INTEGER IFPAR, JFPAR

      INTEGER NCLOCK
      REAL*8  CPUTIM(NCLMAX), VPUTIM(NCLMAX)
      REAL*8  CPUOLD(NCLMAX), VPUOLD(NCLMAX)
      LOGICAL OFIRST
      CHARACTER HTITLE(NCLMAX)*16
      SAVE

      DATA HTITLE / NCLMAX*' ' /
      DATA OFIRST / .TRUE. /
      DATA NCLOCK / 0 /

      CALL REWNML(IFPAR, JFPAR)

      DO 100 IC = 1, NCLOCK
         IF (HTITLE(IC) .NE. ' ') THEN
            WRITE(JFPAR, 51) HTITLE(IC), CPUTIM(IC), VPUTIM(IC)
         END IF
  100 CONTINUE
   51 FORMAT(' ', A16, 2F15.6)

      CALL YCLOCK(CPUT, VPUT)
      WRITE(JFPAR, 51) ' TOTAL TIME = ', CPUT, VPUT

      RETURN

c =====================================================================

      ENTRY CLCSTR(HTTL)

      IF (OFIRST) THEN
         OFIRST = .FALSE.
         CALL YCLOCL
      END IF

      CALL YCLOCK(CPUT, VPUT)
      DO 200 IC = 1, NCLOCK
         IF (HTITLE(IC) .EQ. HTTL) THEN
            CPUOLD(IC) = CPUT
            VPUOLD(IC) = VPUT
            RETURN
         END IF
  200 CONTINUE
      IF (NCLOCK .LE. NCLMAX) THEN
         NCLOCK = NCLOCK + 1
         HTITLE(NCLOCK) = HTTL
         CPUOLD(NCLOCK) = CPUT
         VPUOLD(NCLOCK) = VPUT
      END IF

      RETURN

c =====================================================================

      ENTRY CLCEND(HTTL)

      CALL YCLOCK(CPUT, VPUT)
      DO 300 IC = 1, NCLOCK
         IF (HTITLE(IC) .EQ. HTTL) THEN
            CPUTIM(IC) = CPUTIM(IC) + CPUT - CPUOLD(IC)
            VPUTIM(IC) = VPUTIM(IC) + VPUT - VPUOLD(IC)
            RETURN
         END IF
  300 CONTINUE

      RETURN
      END

c *********************************************************************

      SUBROUTINE YCLOCK(CPUTIM, VPUTIM)

      REAL*8 CPUTIM, VPUTIM

      REAL*4 TARRAY(2)
      REAL*4 ETIME

      REAL*8 TICKS, TICK0, TUSR0
      SAVE   TICKS, TICK0, TUSR0

      DATA TICK0, TUSR0 / 0., 0. /

      TICKS = ETIME(TARRAY)
      CPUTIM = TICKS - TICK0
      VPUTIM = TARRAY(1) - TUSR0

      RETURN

c =====================================================================

      ENTRY YCLOCL

      TICK0 = ETIME(TARRAY)
      TUSR0 = TARRAY(1)

      RETURN
      END
