      SUBROUTINE FWATER(
     M                AX,    HIX,    HSX,
     M              PREC,   SNOW,
     O              EVAP,   SUBI, ADJLAT,
     I               WEV,    WSB,   SOFF)
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  HISTORY
c     '02.10.10  H.Hasumi: from MIROC3.1
c     '07.04.23  H.Hasumi
c     '07.09.25  H.Hasumi: multi-category sea ice
c     '08.06.11  H.Hasumi: initial/final processing
c     '08.07.10  H.Hasumi: initial/final processing
c
c ---------------------------------------------------------------------

#include "zocdim.F"
#include "zocgrd.F"
#include "zocmsk.F"
#include "zocphy.F"

      REAL*8      AX(NXYDIM, 0:NIC),    HIX(NXYDIM, 0:NIC)
      REAL*8     HSX(NXYDIM, 0:NIC)
      REAL*8     WEV(NXYDIM),    WSB(NXYDIM, NIC)
      REAL*8    EVAP(NXYDIM), ADJLAT(NXYDIM)
      REAL*8    SUBI(NXYDIM, NIC)
      REAL*8    PREC(NXYDIM),   SNOW(NXYDIM),   SOFF(NXYDIM)

      REAL*8     RRI,    RRS
      LOGICAL OFIRST
      SAVE RRI, RRS, OFIRST

      REAL*8    WDIF
      REAL*8     DHS,    DHI
      INTEGER     IJ,      K
      INTEGER  IFPAR,  JFPAR

      REAL*8      AZ(NXYDIM, 0:NIC),    HIZ(NXYDIM, 0:NIC)
      REAL*8     HSZ(NXYDIM, 0:NIC)
      COMMON /WORK/ AZ, HIZ, HSZ

      DATA OFIRST / .TRUE. /

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF

      IF (OFIRST) THEN
         OFIRST = .FALSE.
         RRI = rhoo / rhoi
         RRS = rhoo / rhos
      END IF

      DO K = 0, NIC
         DO IJ = 1, NXYDIM
            AZ (IJ, K) = AX(IJ, K)
            HIZ(IJ, K) = HIX(IJ, K)
            HSZ(IJ, K) = HSX(IJ, K)
         END DO
      END DO

      DO K = 1, NIC
         DO IJ = IJTSTR, IJTEND
            IF (AZ(IJ, K) .GT. 0.D0) THEN
               IF (HSZ(IJ, K) .GT. 0.D0) THEN
                  DHS = TS * RRS * WSB(IJ, K) / AX(IJ, K)
                  HSX(IJ, K) = MAX(HSZ(IJ, K) - DHS, 0.D0)
                  WSB(IJ, K) = WSB(IJ, K)
     &                       + AX(IJ, K) * (HSX(IJ, K) - HSZ(IJ, K))
     &                         / RRS / TS
               END IF
               DHI = TS * RRI * WSB(IJ, K) / AX(IJ, K)
               HIX(IJ, K) = MAX(HIZ(IJ, K) - DHI, 0.D0)
               SUBI(IJ, K) = AX(IJ, K) * (HIZ(IJ, K) - HIX(IJ, K))
     &                       / RRI / TS
               WSB(IJ, K) = WSB(IJ, K) - SUBI(IJ, K)
            END IF
         END DO
      END DO

      DO K = 1, NIC
         DO IJ = IJTSTR, IJTEND
            EVAP(IJ) = EVAP(IJ) + WSB(IJ, K)
         END DO
      END DO
      DO IJ = IJTSTR, IJTEND
         ADJLAT(IJ) = EVAP(IJ) - WEV(IJ)
      END DO

      DO IJ = IJTSTR, IJTEND
c         PREC(IJ) = PREC(IJ) - SNOW(IJ)
         SNOW(IJ) = SNOW(IJ) + SOFF(IJ)
      END DO

      RETURN
      END
