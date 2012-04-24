      SUBROUTINE PRIDGE(
     O               PICE,
     M                 AX,    HIX,    EIX,    HSX,    TIX,
     I                 AZ,    HIZ,    EIZ,    HSZ,
     I                 UI,     VI)
      IMPLICIT NONE

c ---- information ----------------------------------------------------
c
c  Mechanical redistribution of sea ice among thickness categories.
c
c  HISTORY
c     '03.08.06  H.Hasumi
c     '07.09.25  H.Hasumi: for COCO4
c     '07.10.03  H.Hasumi: 1-layer sea ice thermodynamics
c     '08.02.??  Y.Komuro
c     '08.06.11  H.Hasumi: initial/final processing
c     '08.07.10  H.Hasumi: initial/final processing
c     '08.11.22  Y.Komuro: coping with the situation HRIDGE = HIX
c
c ---------------------------------------------------------------------

#include "zocdim.F"
#include "zocgrd.F"
#include "zocmsk.F"
#include "zocphy.F"

      REAL*8    PICE(NXYDIM)
      REAL*8      AX(NXYDIM, 0:NIC),     AZ(NXYDIM, 0:NIC)
      REAL*8     HIX(NXYDIM, 0:NIC),    HIZ(NXYDIM, 0:NIC)
      REAL*8     EIX(NXYDIM, 0:NIC),    EIZ(NXYDIM, 0:NIC)
      REAL*8     TIX(NXYDIM, 0:NIC)
      REAL*8     HSX(NXYDIM, 0:NIC),    HSZ(NXYDIM, 0:NIC)
      REAL*8      UI(NXYDIM),     VI(NXYDIM)

      REAL*8   AXHIX(NXYDIM, 0:NIC),  AXHSX(NXYDIM, 0:NIC)
      REAL*8   AXEIX(NXYDIM, 0:NIC)
      REAL*8    DIVV(NXYDIM),  EDIS(NXYDIM)
      REAL*8      WA(NXYDIM, 0:NIC),     WN(NXYDIM, 0:NIC)
      REAL*8      WW(NXYDIM)
      REAL*8       G(NXYDIM, -1:NIC)
      REAL*8       Y(NXYDIM, -1:NIC)
      REAL*8      DA(NXYDIM, NIC)
      REAL*8    DAHI(NXYDIM, NIC),   DAHS(NXYDIM, NIC)
      REAL*8    DAEI(NXYDIM, NIC)
      COMMON /WORK/ AXHIX, AXHSX,
     &              DIVV, EDIS, WA, WN, WW,
     &              G, Y, DA, DAHI, DAHS, DAEI

      REAL*8      C1,     C2,     C3
      REAL*8     GAM(NXYDIM, NIC, NIC)
      REAL*8   PIFCT
      REAL*8     TMI
      LOGICAL OFIRST
      SAVE C1, C2, C3, GAM, PIFCT, TMI, OFIRST

      REAL*8     EXX,    EYY,    EXY,    DEL
      REAL*8   HRMAX,  HRMIN
      REAL*8    HIKL,   HSKL,   HEKL
      INTEGER     IJ,      K,      L
      INTEGER   IJLW,   IJLS,  IJLSW
      INTEGER  IFPAR,  JFPAR

      REAL*8     ECC,   DMIN,  FLOSS
      REAL*8      CS, GRIDGE, HRIDGE
      REAL*8      SI
      SAVE ECC, DMIN, FLOSS, CS, GRIDGE, HRIDGE, SI
      NAMELIST /NMIDYN/ ECC, DMIN, FLOSS
      NAMELIST /NMIRDG/ CS, GRIDGE, HRIDGE
      NAMELIST /NMISLT/ SI

#include "zocite.F"

      DATA OFIRST / .TRUE. /
      DATA ECC, DMIN, FLOSS / 2.D0, 2.0D-7, 17.D0 /
      DATA CS, GRIDGE, HRIDGE / 0.5D0, 0.15D0, 1.D4 /
      DATA SI / 5.D0 /

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF

      IF (OFIRST) THEN
         OFIRST = .FALSE.
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMIDYN, END=99)
 99      WRITE(JFPAR, NMIDYN)
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMIRDG, END=98)
 98      WRITE(JFPAR, NMIRDG)
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMISLT, END=97)
 97      WRITE(JFPAR, NMISLT)

         TMI = dtds * SI

         C1 = 1.D0 + 1.D0 / ECC / ECC
         C2 = 4.D0 / ECC / ECC
         C3 = 2.D0 * (1.D0 - 1.D0 / ECC / ECC)

         PIFCT = 0.5D0 * rhoi / rhoo * gravit * (rhoo - rhoi) * FLOSS

         DO L = 1, NIC
            DO K = 1, NIC
               DO IJ = 1, NXYDIM
                  GAM(IJ, K, L) = 0.D0
               END DO
            END DO
         END DO
      END IF

      DO K = 0, NIC
         DO IJ = 1, NXYDIM
            AXHIX(IJ, K) = AX(IJ, K) * HIX(IJ, K)
            AXHSX(IJ, K) = AX(IJ, K) * HSX(IJ, K)
            AXEIX(IJ, K) = AX(IJ, K) * EIX(IJ, K)
            G(IJ, K) = 0.D0
            WN(IJ, K) = 0.D0
         END DO
      END DO
      DO IJ = 1, NXYDIM
         PICE(IJ) = 0.D0
         G(IJ, -1) = 0.D0
      END DO

      DO IJ = IJTSTR, IJTEND
         IJLW = IJ + LW
         IJLS = IJ + LS
         IJLSW = IJ + LSW
         EXX = (UI(IJ) + UI(IJLS) - UI(IJLW) - UI(IJLSW)) *
     &         RX * RXT(IJ) * 0.5D0
     &       + (VI(IJ) + VI(IJLS) + VI(IJLW) + VI(IJLSW)) *
     &         HXYT(IJ) * 0.25D0
         EYY = (VI(IJ) + VI(IJLW) - VI(IJLS) - VI(IJLSW)) *
     &         RY(IJ) * RYT(IJ) * 0.5D0
     &       + (UI(IJ) + UI(IJLW) + UI(IJLS) + UI(IJLSW)) *
     &         HYXT(IJ) * 0.25D0
         EXY = (UI(IJ) + UI(IJLW) - UI(IJLS) - UI(IJLSW)) *
     &         RY(IJ) * RYT(IJ) * 0.25D0
     &       - (VI(IJ) + VI(IJLW) + VI(IJLS) + VI(IJLSW)) *
     &         HYXT(IJ) * 0.125D0
     &       + (VI(IJ) + VI(IJLS) - VI(IJLW) - VI(IJLSW)) *
     &         RX * RXT(IJ) * 0.25D0
     &       - (UI(IJ) + UI(IJLS) + UI(IJLW) + UI(IJLSW)) *
     &         HXYT(IJ) * 0.125D0
         DIVV(IJ) = EXX + EYY
         DEL = SQRT(  C1 * (EXX * EXX + EYY * EYY)
     &              + C2 * EXY * EXY + C3 * EXX * EYY)
         EDIS(IJ) = 0.5D0 * CS * (DEL - ABS(DIVV(IJ)))
     &            - MIN(DIVV(IJ), 0.D0)
      END DO

      DO IJ = IJTSTR, IJTEND
         G(IJ, 0) = max(AX(IJ, 0), 0.0D0)
      END DO
      DO K = 1, NIC
         DO IJ = IJTSTR, IJTEND
            G(IJ, K) = G(IJ, K-1) + AX(IJ, K)
         END DO
      END DO

      DO K = -1, NIC
         DO IJ = IJTSTR, IJTEND
            Y(IJ, K) = MAX(0.D0, 1.D0 - G(IJ, K) / GRIDGE)**2
         END DO
      END DO

      DO K = 0, NIC
         DO IJ = IJTSTR, IJTEND
            WA(IJ, K) = Y(IJ, K-1) - Y(IJ, K)
         END DO
      END DO

      DO L = 1, NIC
         DO K = 1, L
            DO IJ = IJTSTR, IJTEND
               HRMAX = 2.D0 * SQRT(HRIDGE * HIX(IJ, K))
               HRMIN = 2.D0 * HIX(IJ, K)
               IF (     (HIC(L+1) .LT. HRMIN)
     &             .OR. (HIC(L) .GT. HRMAX)
     &             .OR. (HRMIN .GE. HRMAX) ) THEN
                  GAM(IJ, K, L) = 0.D0
               ELSE
                  GAM(IJ, K, L) = (  MIN(HRMAX, HIC(L+1))
     &                             - MAX(HRMIN, HIC(L))) * 0.5D0
     &                            / (HRIDGE - HIX(IJ, K))
               END IF
            END DO
         END DO
      END DO

      DO L = 1, NIC
         DO K = 1, L
            DO IJ = IJTSTR, IJTEND
               WN(IJ, L) = WN(IJ, L) + WA(IJ, K) * GAM(IJ, K, L)
           END DO
         END DO
      END DO

      DO IJ = 1, NXYDIM
         WW(IJ) = 0.D0
      END DO
      DO K = 0, NIC
         DO IJ = IJTSTR, IJTEND
            WW(IJ) = WW(IJ) + WA(IJ, K) - WN(IJ, K)
         END DO
      END DO
      DO K = 0, NIC
         DO IJ = IJTSTR, IJTEND
            WA(IJ, K) = WA(IJ, K) / WW(IJ)
            WN(IJ, K) = WN(IJ, K) / WW(IJ)
         END DO
      END DO

      DO K = 1, NIC
         DO IJ = IJTSTR, IJTEND
            DA(IJ, K) = EDIS(IJ) * (WN(IJ, K) - WA(IJ, K))
         END DO
      END DO
      DO L = 1, NIC
         DO IJ = IJTSTR, IJTEND
            DAHI(IJ, L) = - HIX(IJ, L) * WA(IJ, L) * EDIS(IJ)
            DAHS(IJ, L) = - HSX(IJ, L) * WA(IJ, L) * EDIS(IJ)
            DAEI(IJ, L) = - EIX(IJ, L) * WA(IJ, L) * EDIS(IJ)
            PICE(IJ) = PICE(IJ)
     &               - HIX(IJ, L) * HIX(IJ, L) * WA(IJ, L) * PIFCT
            DO K = 1, L
               HRMAX = 2.D0 * SQRT(HRIDGE * HIX(IJ, K))
               HRMIN = 2.D0 * HIX(IJ, K)
               HIKL = (  MIN(HRMAX, HIC(L+1))
     &                 + MAX(HRMIN, HIC(L))) * 0.5D0
               HSKL = HSX(IJ, K) * (HIX(IJ, K) + HRMAX * 0.5D0)
     &                / HIX(IJ, K)
               IF (TIX(IJ, K) .LT. TMI) THEN
                  HEKL = HIKL * EI(TIX(IJ, K), SI)
               ELSE
                  HEKL = 0.D0
               END IF
               DAHI(IJ, L) = DAHI(IJ, L)
     &                     + HIKL * WA(IJ, K) * GAM(IJ, K, L) *
     &                       EDIS(IJ)
               DAHS(IJ, L) = DAHS(IJ, L)
     &                     + HSKL * WA(IJ, K) * GAM(IJ, K, L) *
     &                       EDIS(IJ)
               DAEI(IJ, L) = DAEI(IJ, L)
     &                     + HEKL * WA(IJ ,K) * GAM(IJ, K, L) *
     &                       EDIS(IJ)
               PICE(IJ) = PICE(IJ)
     &                  + HIKL * HIKL * WA(IJ, K) * GAM(IJ, K, L) *
     &                    PIFCT
            END DO
         END DO
      END DO

      DO K = 1, NIC
         DO IJ = IJTSTR, IJTEND
            AX(IJ, K) = AX(IJ, K) + TS * DA(IJ, K)
            AX(IJ, K) = MIN(1.D0, MAX(0.D0, AX(IJ, K)))
         END DO
      END DO
      DO IJ = IJTSTR, IJTEND
         AX(IJ, 0) = 1.D0
      END DO
      DO K = 1, NIC
         DO IJ = IJTSTR, IJTEND
            AX(IJ, 0) = AX(IJ, 0) - AX(IJ, K)
         END DO
      END DO
      DO IJ = IJTSTR, IJTEND
         IF (AX(IJ, 0) .LT. 0.D0) THEN
            DO K = 1, NIC
               AX(IJ, K) = AX(IJ, K) / (1.D0 - AX(IJ, 0))
            END DO
            AX(IJ, 0) = 0.D0
         END IF
      END DO

      DO IJ = IJTSTR, IJTEND
         HIX(IJ, 0) = 0.D0
         HSX(IJ, 0) = 0.D0
         EIX(IJ, 0) = 0.D0
         IF (AX(IJ, 0) .EQ. 1.D0) THEN
            DO K = 1, NIC
               AX(IJ, K) = 0.D0
               HIX(IJ, K) = HIC(K)
               HSX(IJ, K) = 0.D0
               EIX(IJ, K) = 0.D0
               HIX(IJ, 0) = HIX(IJ, 0) + AXHIX(IJ, K)
               HSX(IJ, 0) = HSX(IJ, 0) + AXHSX(IJ, K)
               EIX(IJ, 0) = EIX(IJ, 0) + AXEIX(IJ, K)
               DAHI(IJ, K) = 0.D0
               DAHS(IJ, K) = 0.D0
               DAEI(IJ, K) = 0.D0
            END DO
         END IF
      END DO

      DO K = 1, NIC
         DO IJ = IJTSTR, IJTEND
            AXHIX(IJ, K) = AXHIX(IJ, K) + TS * DAHI(IJ, K)
            AXHSX(IJ, K) = AXHSX(IJ, K) + TS * DAHS(IJ, K)
            AXEIX(IJ, K) = AXEIX(IJ, K) + TS * DAEI(IJ, K)
         END DO
      END DO

*** Ridging process conserves the mass of ice and snow, and the
*** amount of the thickest category ice is never reduced by ridging.
*** Thus, the situation that AXHIX(IJ, NIC) becomes less than zero
*** needs not be considered in the adjustment below. 
      DO K = 1, NIC-1
         DO IJ = IJTSTR, IJTEND
            IF (AX(IJ, K) .EQ. 0.D0) THEN
               AXHIX(IJ, K+1) = AXHIX(IJ, K+1) + AXHIX(IJ, K)
               AXHSX(IJ, K+1) = AXHSX(IJ, K+1) + AXHSX(IJ, K)
               AXEIX(IJ, K+1) = AXEIX(IJ, K+1) + AXEIX(IJ, K)
               AXHIX(IJ, K) = 0.D0
               AXHSX(IJ, K) = 0.D0
               AXEIX(IJ, K) = 0.D0
               HIX(IJ, K) = HIC(K)
               HSX(IJ, K) = 0.D0
               EIX(IJ, K) = 0.D0
            ELSE IF (AXHIX(IJ, K) .LT. 0.D0) THEN
               AXHIX(IJ, K+1) = AXHIX(IJ, K+1) + AXHIX(IJ, K)
               AXHSX(IJ, K+1) = AXHSX(IJ, K+1) + AXHSX(IJ, K)
               AXEIX(IJ, K+1) = AXEIX(IJ, K+1) + AXEIX(IJ, K)
               AX(IJ, K+1) = AX(IJ, K+1) + AX(IJ, K)
               AX(IJ, K) = 0.D0
               AXHIX(IJ, K) = 0.D0
               AXHSX(IJ, K) = 0.D0
               AXEIX(IJ, K) = 0.D0
               HIX(IJ, K) = HIC(K)
               HSX(IJ, K) = 0.D0
               EIX(IJ, K) = 0.D0
            ELSE
               IF (AXHSX(IJ, K) .LT. 0.D0) THEN
                  AXHSX(IJ, K+1) = AXHSX(IJ, K+1) + AXHSX(IJ, K)
                  AXHSX(IJ, K) = 0.D0
                  HSX(IJ, K) = 0.D0
               END IF
               IF (AXEIX(IJ, K) .LT. 0.D0) THEN
                  AXEIX(IJ, K+1) = AXEIX(IJ, K+1) + AXEIX(IJ, K)
                  AXEIX(IJ, K) = 0.D0
                  EIX(IJ, K) = 0.D0
               END IF
            END IF
         END DO
      END DO

      DO K = 1, NIC
         DO IJ = IJTSTR, IJTEND
            IF (AX(IJ, K) .GT. 0.D0) THEN
               HIX(IJ, K) = AXHIX(IJ, K) / AX(IJ, K)
               HSX(IJ, K) = AXHSX(IJ, K) / AX(IJ, K)
               EIX(IJ, K) = AXEIX(IJ, K) / AX(IJ, K)
               TIX(IJ, K) = TI(EIX(IJ, K)/HIX(IJ, K), SI)
            ELSE
               TIX(IJ, K) = TMI
            END IF
         END DO
      END DO

      RETURN
      END
