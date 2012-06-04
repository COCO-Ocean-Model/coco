      SUBROUTINE PTHERM(
     M                AX,    HIX,    HSX,
     M               EIX,    TIX,
     M              PREC,   SNOW,     FT,     FS,
     O             FTITD, IGRFRA, IGRCON, IGRSNI,
     O            INRLAT, IMRSNO, IMRISF, IMRIBS,
     I                TX,
     I               WIO,    WAO,    WAS,    WIL,
     I              EVAP,   SUBI,   ROFF, ADJLAT,
     I               QIO)
      use zocdim
      use zocgrd
      use zocmsk
      use zocphy
c      use zocite
      IMPLICIT NONE

c ---- information ----------------------------------------------------
c
c  HISTORY
c     '02.10.10  H.Hasumi: from MIROC3.1
c     '03.07.30  H.Hasumi: multi-category sea ice thickness
c     '07.09.25  H.Hasumi: for COCO4
c     '07.10.03  H.Hasumi: 1-layer sea ice thermodynamics
c     '08.04.07  H.Hasumi: condition is modified for new ice formation
c                          on open water
c     '08.06.11  H.Hasumi: initial/final processing
c     '08.07.10  H.Hasumi: initial/final processing
c     '09.05.25  Y.Komuro: CMIP5 output code included
c                          (basal/lateral melting processes separated)
c
c ---------------------------------------------------------------------

      REAL*8      AX(NXYDIM, 0:NIC),    HIX(NXYDIM, 0:NIC)
      REAL*8     HSX(NXYDIM, 0:NIC)
      REAL*8     EIX(NXYDIM, 0:NIC),    TIX(NXYDIM, 0:NIC)
      REAL*8     WIO(NXYDIM, NIC),    WAS(NXYDIM, NIC)
      REAL*8     WIL(NXYDIM, NIC)
      REAL*8     WAO(NXYDIM)
      REAL*8    PREC(NXYDIM),   SNOW(NXYDIM)
      REAL*8    EVAP(NXYDIM),   SUBI(NXYDIM, NIC)
      REAL*8    ROFF(NXYDIM), ADJLAT(NXYDIM)
      REAL*8      FT(NXYDIM, NTDIM),     FS(NXYDIM)
      REAL*8      TX(NXYDIM, NZDIM, NTDIM)
      REAL*8     QIO(NXYDIM, NIC)
      REAL*8   FTITD(NXYDIM)
      REAL*8  IGRFRA(NXYDIM), IGRCON(NXYDIM), IGRSNI(NXYDIM)
      REAL*8  INRLAT(NXYDIM)
      REAL*8  IMRSNO(NXYDIM), IMRISF(NXYDIM), IMRIBS(NXYDIM)

      REAL*8      AZ(NXYDIM, 0:NIC)
      REAL*8   AXHIX(NXYDIM, 0:NIC)
      REAL*8   AXHSX(NXYDIM, 0:NIC), AXHSXN(NXYDIM, 0:NIC)
      REAL*8   AXEIX(NXYDIM, 0:NIC), AXEIXN(NXYDIM, 0:NIC)
      REAL*8     WAI(NXYDIM, NIC)
      REAL*8      WI(NXYDIM),     WS(NXYDIM)
      REAL*8     WEN(NXYDIM),    WSN(NXYDIM)
      REAL*8  LAXHIX(NXYDIM, 0:NIC), LAXEIX(NXYDIM, 0:NIC)
      REAL*8  DAXEIT(NXYDIM, 0:NIC), DAXEIB(NXYDIM, 0:NIC)
c      COMMON /WORK/ AZ, AXHIX, AXHSX, AXHSXN, WAI,
c     &              WI, WS, WEN, WSN, AXEIX, AXEIXN,
c     &              LAXHIX, LAXEIX, DAXEIT, DAXEIB

      REAL*8     RRI,    RRS, RORIRS
      REAL*8   RSFUS
      REAL*8     TMI
      REAL*8  EPSAEI
      LOGICAL OFIRST
      SAVE RRI, RRS, RORIRS, RSFUS, TMI, OFIRST, EPSAEI

      REAL*8    WRES
      REAL*8    HSXO,    DHS
      REAL*8  AX1MAX
      REAL*8  DAXHIX, DAXEIX, DAXHIT, DAXHIB,   EIEQ
      INTEGER     IJ,      K
      INTEGER  IFPAR,  JFPAR

      REAL*8     EPS, EPSL
      SAVE EPS, EPSL

      REAL*8    AMIN,   AMAX,     SI
      INTEGER    MIC
      NAMELIST /NMAMIN/ AMIN, AMAX, MIC
      NAMELIST /NMISLT/ SI
      SAVE AMIN, SI

#include "zocite.F"

      DATA OFIRST / .TRUE. /
      DATA SI / 5.D0 /
      DATA AMIN, AMAX / 1.0D-6, 1.D0 /
      DATA MIC / NIC /
      DATA EPS, EPSL / 1.D-3, 1.0D-6 /

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF

      IF (OFIRST) THEN
         OFIRST = .FALSE.
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMAMIN, END=98)
 98      WRITE(JFPAR, NMAMIN)
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMISLT, END=97)
 97      WRITE(JFPAR, NMISLT)

         TMI = dtds * SI
         RRI    = rhoo / rhoi
         RRS    = rhoo / rhos
         RORIRS = (rhoo - rhoi) / rhos
         RSFUS  = rhos * hfus
         EPSAEI = AMIN * ABS(EI(TMI+EPS, SI))
      END IF

      DO IJ = 1, NXYDIM
         IGRFRA(IJ) = 0.0D0
         IGRCON(IJ) = 0.0D0
         IGRSNI(IJ) = 0.0D0
         INRLAT(IJ) = 0.0D0
         IMRSNO(IJ) = 0.0D0
         IMRISF(IJ) = 0.0D0
         IMRIBS(IJ) = 0.0D0
      END DO

      DO K = 0, NIC
         DO IJ = 1, NXYDIM
            AZ(IJ, K) = AX(IJ, K)
            DAXEIT(IJ, K) = 0.0D0
            DAXEIB(IJ, K) = 0.0D0
         END DO
      END DO

      DO K = 0, NIC
         DO IJ = 1, NXYDIM
            AXHIX(IJ, K) = AX(IJ, K) * HIX(IJ, K)
            AXHSX(IJ, K) = AX(IJ, K) * HSX(IJ, K)
            AXEIX(IJ, K) = AX(IJ, K) * EIX(IJ, K)
         END DO
      END DO

*** snowfall ***
      DO K = 1, NIC
         DO IJ = IJTSTR, IJTEND
            IF (AX(IJ, K) .GT. 0.D0) THEN
               HSX(IJ, K) = AXHSX(IJ, K) / AX(IJ, K)
     &                    + TS * RRS * SNOW(IJ)
               AXHSX(IJ, K) = AX(IJ, K) * HSX(IJ, K)
            ELSE
               HSX(IJ, K) = 0.D0
               AXHSX(IJ, K) = 0.D0
            END IF
         END DO
      END DO
      DO IJ = IJTSTR, IJTEND
         SNOW(IJ) = AX(IJ, 0) * SNOW(IJ)
         PREC(IJ) = PREC(IJ) + SNOW(IJ)
      END DO

*** snow melting ***
      DO K = 1, NIC
         DO IJ = IJTSTR, IJTEND
            WRES = AXHSX(IJ, K) * RSFUS / TS + WAS(IJ, K)
            IF (AX(IJ, K) .GT. 0.D0) THEN
               IF (WRES .LT. 0.D0) THEN
                  HSX(IJ, K) = 0.D0
                  WAI(IJ, K) = WRES
               ELSE
                  HSX(IJ, K) = WRES * TS / RSFUS / AX(IJ, K)
                  WAI(IJ, K) = 0.D0
               END IF
            ELSE
               WAI(IJ, K) = WAS(IJ, K)
            END IF
            AXHSXN(IJ, K) = AX(IJ, K) * HSX(IJ, K) * AMSKT(IJ, KSTR)
            IMRSNO(IJ) = IMRSNO(IJ) -
     &                   ( AXHSXN(IJ, K) - AXHSX(IJ, K) )
         END DO
      END DO

*** ice top melting ***
      DO K = 1, NIC
         DO IJ = IJTSTR, IJTEND
            WRES = rhoi * AXEIX(IJ, K) / TS + WAI(IJ, K)
            IF (AX(IJ, K) .GT. 0.D0) THEN
               IF (WRES .LT. 0.D0) THEN
                  AXEIXN(IJ, K) = 0.D0
                  WIO(IJ, K) = WIO(IJ, K) + WRES
               ELSE
                  AXEIXN(IJ, K) = WRES / rhoi * TS
               END IF
            ELSE
               WIO(IJ, K) = WIO(IJ, K) + WAI(IJ, K)
               AXEIXN(IJ, K) = AXEIX(IJ, K)
            END IF
            DAXEIT(IJ, K) = AXEIXN(IJ, K) - AXEIX(IJ, K)
         END DO
      END DO

*** new ice formation on open water ***
      DO IJ = IJTSTR, IJTEND
         AXEIXN(IJ, 0) = WAO(IJ) * TS / rhoi * AMSKT(IJ, KSTR)
         AXHSXN(IJ, 0) = 0.D0
         TIX(IJ, 0) = MIN(TX(IJ, KSTR, 1)*AMSKT(IJ, KSTR), TMI)
         IF (      (AX(IJ, 0) .GT. 0.D0)
     &       .AND. (AXEIXN(IJ, 0) .GT. 0.D0)
c     &       .AND. (TIX(IJ, 0) .LT. TMI)) THEN
     &       .AND. (TIX(IJ, 0) .LT. TMI-EPS)) THEN
            EIX(IJ, 0) = AXEIXN(IJ, 0) / AX(IJ, 0)
            HIX(IJ, 0) = EIX(IJ, 0) / EI(TIX(IJ, 0), SI)
            IGRFRA(IJ) = AX(IJ, 0) * HIX(IJ, 0) * rhoi
         ELSE
            EIX(IJ, 0) = 0.D0
            HIX(IJ, 0) = 0.D0
         END IF
      END DO

*** basal ice formation/melting
      DO K = 1, NIC
         DO IJ = IJTSTR, IJTEND
            DAXEIB(IJ, K) = AXEIXN(IJ, K)
            AXEIXN(IJ, K) = AXEIXN(IJ, K)
     &                    + WIO(IJ, K) * TS / rhoi *
     &                      AMSKT(IJ, KSTR)
         END DO
      END DO
      DO K = 1, NIC
         DO IJ = IJTSTR, IJTEND
            IF (AXEIXN(IJ, K) .LE. 0.D0) THEN
               AX(IJ, K) = 0.D0
               HIX(IJ, K) = HIC(K)
               HSX(IJ, K) = 0.D0
               TIX(IJ, K) = TMI
               EIX(IJ, K) = 0.D0
            ELSE
               EIX(IJ, K) = AXEIXN(IJ, K) / AX(IJ, K)
               HIX(IJ, K) = EIX(IJ, K) / EI(TIX(IJ, K), SI)
            END IF
            DAXEIB(IJ, K) = AX(IJ, K) * EIX(IJ, K) - DAXEIB(IJ, K)
         END DO
      END DO

      DO K = 1, NIC
         DO IJ = IJTSTR, IJTEND
            IF (AZ(IJ, K) .GT. 0.0D0) THEN
               DAXHIX = AX(IJ, K) * HIX(IJ, K) - AXHIX(IJ, K)
               DAXEIX = DAXEIT(IJ, K) + DAXEIB(IJ, K)
               IF ( ABS(DAXEIX) .LT. EPSAEI ) THEN
                  DAXHIT = DAXEIT(IJ, K) / EI(TIX(IJ, K), SI)
                  DAXHIB = DAXEIB(IJ, K) / EI(TIX(IJ, K), SI)
               ELSE
                  DAXHIT = DAXHIX * DAXEIT(IJ, K) / DAXEIX
                  DAXHIB = DAXHIX * DAXEIB(IJ, K) / DAXEIX
               ENDIF
               IMRISF(IJ) = IMRISF(IJ) - rhoi * DAXHIT
               IGRCON(IJ) = IGRCON(IJ) + rhoi * MAX(0.0D0, DAXHIB)
               IMRIBS(IJ) = IMRIBS(IJ) - rhoi * MIN(0.0D0, DAXHIB)
            END IF
         END DO
      END DO

*** snow-ice formation ***
      DO K = 1, NIC
         DO IJ = IJTSTR, IJTEND
            HSXO = HSX(IJ, K)
            HSX(IJ, K) = MIN(HSX(IJ, K), RORIRS * HIX(IJ, K))
            DHS = HSXO - HSX(IJ, K)
            HIX(IJ, K) = HIX(IJ, K) + DHS * rhos / rhoi
            EIX(IJ, K) = EIX(IJ, K) + RSFUS / rhoi * DHS
            IGRSNI(IJ) = IGRSNI(IJ) + AX(IJ, K) * DHS * rhos
            IF (AX(IJ, K) .GT. 0.D0) THEN
               TIX(IJ, K) = TI(EIX(IJ, K)/HIX(IJ, K), SI)
            ELSE
               TIX(IJ, K) = TMI
            END IF
         END DO
      END DO

*** lateral ice formation/melting
      DO K = 1, NIC
         DO IJ = IJTSTR, IJTEND
            LAXHIX(IJ, K) = AX(IJ, K) * HIX(IJ, K)
            LAXEIX(IJ, K) = AX(IJ, K) * EIX(IJ, K)
     &                    + WIL(IJ, K) * TS / rhoi *
     &                      AMSKT(IJ, KSTR)
            AXEIXN(IJ, K) = AXEIXN(IJ, K)
     &                    + WIL(IJ, K) * TS / rhoi *
     &                      AMSKT(IJ, KSTR)
         END DO
      END DO
      DO K = 1, NIC
         DO IJ = IJTSTR, IJTEND
            IF (AX(IJ, K) .GT. 0.D0) THEN
               AX(IJ, K) = AX(IJ, K)
     &                   + TS * WIL(IJ, K) / EIX(IJ, K) / rhoi *
     &                     AMSKT(IJ, KSTR)
            END IF
         END DO
      END DO
      DO K = 1, NIC
         DO IJ = IJTSTR, IJTEND
            IF (LAXEIX(IJ, K) .LE. 0.D0) THEN
               AX(IJ, K) = 0.D0
               HIX(IJ, K) = HIC(K)
               HSX(IJ, K) = 0.D0
               TIX(IJ, K) = TMI
               EIX(IJ, K) = 0.D0
            ELSE IF (AX(IJ, K) .LE. 0.D0) THEN
               AX(IJ, K) = 0.D0
               HIX(IJ, K) = HIC(K)
               HSX(IJ, K) = 0.D0
               TIX(IJ, K) = TMI
               EIX(IJ, K) = 0.D0
            ELSE
               EIX(IJ, K) = LAXEIX(IJ, K) / AX(IJ, K)
               TIX(IJ, K) = TI(EIX(IJ, K)/HIX(IJ, K), SI)
            END IF
            INRLAT(IJ) = INRLAT(IJ) + rhoi * 
     &                   ( AX(IJ, K)*HIX(IJ, K) - LAXHIX(IJ, K) )
         END DO
      END DO

*** heat and freshwater budget ***
      DO IJ = 1, NXYDIM
         FTITD(IJ) = 0.D0
      END DO
      DO IJ = IJTSTR, IJTEND
         WI(IJ) = (AX(IJ, 0) * HIX(IJ, 0) - AXHIX(IJ, 0)) / RRI / TS
         WS(IJ) = (AX(IJ, 0) * HSX(IJ, 0) - AXHSX(IJ, 0)) / RRS / TS
         WSN(IJ) = (AX(IJ, 0) * HSX(IJ, 0) - AXHSXN(IJ, 0)) / RRS / TS
         WEN(IJ) = (AX(IJ, 0) * EIX(IJ, 0) - AXEIXN(IJ, 0)) / RRI / TS
         FS(IJ) = 0.D0
      END DO
      DO K = 1, NIC
         DO IJ = IJTSTR, IJTEND
            WI(IJ) = WI(IJ)
     &             + (AX(IJ, K) * HIX(IJ, K) - AXHIX(IJ, K))
     &               / RRI / TS
            WS(IJ) = WS(IJ)
     &             + (AX(IJ, K) * HSX(IJ, K) - AXHSX(IJ, K))
     &               / RRS / TS
            WSN(IJ) = WSN(IJ)
     &              + (AX(IJ, K) * HSX(IJ, K) - AXHSXN(IJ, K))
     &                / RRS / TS
            WEN(IJ) = WEN(IJ)
     &              + (AX(IJ, K) * EIX(IJ, K) - AXEIXN(IJ, K))
     &                / RRI / TS
            FS(IJ) = FS(IJ) - SUBI(IJ, K) * SI * AMSKT(IJ, KSTR)
            FTITD(IJ) = FTITD(IJ)
     &                - QIO(IJ, K) * AZ(IJ, K) * AMSKT(IJ, KSTR)
         END DO
      END DO
      DO IJ = IJTSTR, IJTEND
         FT(IJ, 2) = (  EVAP(IJ) - PREC(IJ) - ROFF(IJ)
     &                + WS(IJ) + WI(IJ)) * AMSKT(IJ, KSTR)
         FS(IJ) = FS(IJ) + WI(IJ) * SI * AMSKT(IJ, KSTR)
         FT(IJ, 1) = - FT(IJ, 1)
     &               + hfus / cpo * (WSN(IJ) - SNOW(IJ))
     &               + WEN(IJ) / cpo
     &               + ADJLAT(IJ) * hfus / cpo
         FT(IJ, 1) = FT(IJ, 1) * AMSKT(IJ, KSTR)
         FTITD(IJ) = FTITD(IJ) + AMSKT(IJ, KSTR) *
     &               ( rhoo * hfus * WI(IJ) 
     &               + rhoo * hfus * WSN(IJ) )
      END DO

*** merging newly formed ice into the category 1 ***
      DO IJ = IJTSTR, IJTEND
         AXHIX(IJ, 1) = AX(IJ, 1) * HIX(IJ, 1)
     &                + AX(IJ, 0) * HIX(IJ, 0)
         AXEIX(IJ, 1) = AX(IJ, 1) * EIX(IJ, 1)
     &                + AX(IJ, 0) * EIX(IJ, 0)
         AXHSX(IJ, 1) = AX(IJ, 1) * HSX(IJ, 1)
         HIX(IJ, 1) = MAX(HIX(IJ, 1), HIC(1))
         HIX(IJ, 0) = 0.D0
         EIX(IJ, 0) = 0.D0
         TIX(IJ, 0) = TMI
         HSX(IJ, 0) = 0.D0
      END DO
      DO IJ = IJTSTR, IJTEND
         AX(IJ, 1) = AXHIX(IJ, 1) / HIX(IJ, 1)
         AX1MAX = AZ(IJ, 0) + AZ(IJ, 1)
         IF (AX(IJ, 1) .GT. AX1MAX) THEN
            AX(IJ, 1) = AX1MAX
            HIX(IJ, 1) = AXHIX(IJ, 1) / AX1MAX
            EIX(IJ, 1) = AXEIX(IJ, 1) / AX1MAX
            TIX(IJ, 1) = TI(EIX(IJ, 1) / HIX(IJ, 1), SI)
            HSX(IJ, 1) = AXHSX(IJ, 1) / AX1MAX
         ELSE IF (AX(IJ, 1) .GT. 0.D0) THEN
            HSX(IJ, 1) = AXHSX(IJ, 1) / AX(IJ, 1)
            EIX(IJ, 1) = AXEIX(IJ, 1) / AX(IJ, 1)
            TIX(IJ, 1) = TI(EIX(IJ, 1) / HIX(IJ, 1), SI)
         END IF
      END DO

c     for CMIP5 output: unit conversion
      DO IJ = IJTSTR, IJTEND
         IGRFRA(IJ) = IGRFRA(IJ) / TS
         IGRCON(IJ) = IGRCON(IJ) / TS
         IGRSNI(IJ) = IGRSNI(IJ) / TS
         INRLAT(IJ) = INRLAT(IJ) / TS
         IMRSNO(IJ) = IMRSNO(IJ) / TS
         IMRISF(IJ) = IMRISF(IJ) / TS
         IMRIBS(IJ) = IMRIBS(IJ) / TS
      END DO

      RETURN
      END
