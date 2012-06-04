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
c     '08.12.12  Y.Komuro: linear remapping of Lipscomb(2001)
c     '09.05.25  Y.Komuro: CMIP5 output code included
c
c ---------------------------------------------------------------------

#ifdef OPT_PARALLEL
#include "mpif.h"
#endif

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

      REAL*8     HIZ(NXYDIM, 0:NIC)
      REAL*8   AFLRM(NXYDIM, 0:NIC), AFLRMC(NXYDIM)
      REAL*8    FDTN(NXYDIM, 0:NIC),  FDTCN(NXYDIM, 0:NIC)
      REAL*8    HICN(NXYDIM, 0:NIC)
      REAL*8      G0(NXYDIM, 0:NIC),     G1(NXYDIM, 0:NIC)
      REAL*8     HIL(NXYDIM, 0:NIC),    HIR(NXYDIM, 0:NIC)
      REAL*8      DA(NXYDIM, 0:NIC),   DAHI(NXYDIM, 0:NIC)
      REAL*8    DAHS(NXYDIM, 0:NIC),   DAEI(NXYDIM, 0:NIC)
      REAL*8  LAXHIX(NXYDIM, 0:NIC), LAXHSX(NXYDIM, 0:NIC) 
      REAL*8  LAXEIX(NXYDIM, 0:NIC)
      REAL*8  DAXEIT(NXYDIM, 0:NIC), DAXEIB(NXYDIM, 0:NIC)
c      COMMON /WORK/ AZ, AXHIX, AXHSX, AXHSXN, WAI,
c     &              WI, WS, WEN, WSN, AXEIX, AXEIXN,
c     &              HIZ, AFLRM, AFLRMC, FDTN, FDTCN, HICN,
c     &              G0, G1, HIL, HIR, DA, DAHI,
c     &              DAHS, DAEI, LAXHIX, LAXHSX, LAXEIX,
c     &              DAXEIT, DAXEIB


      REAL*8     RRI,    RRS, RORIRS
      REAL*8   RSFUS
      REAL*8     TMI
      REAL*8  EPSAEI
      LOGICAL OFIRST
      REAL*8    HIC0(0:NIC+1)
      SAVE RRI, RRS, RORIRS, RSFUS, TMI, OFIRST, HIC0, EPSAEI

      REAL*8    WRES
      REAL*8    HSXO,    DHS
      REAL*8  AX1MAX
      REAL*8    ETAN,   ETAR, ETANRR,   GIL,    GIR
      REAL*8      X0,     X1,   GINT
      REAL*8    FAHI,   FAHS,   FAEI,   PVOL
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
         HIC0(0) = 0.0D0
         DO K = 1, NIC+1
            HIC0(K) = HIC(K)
         END DO
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
            HIZ(IJ, K) = HIX(IJ, K)
            AFLRM(IJ, K) = 0.0D0
            DAXEIT(IJ,K) = 0.0D0
            DAXEIB(IJ,K) = 0.0D0
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

*** basal and lateral ice formation/melting processes are divided
*** in order to apply linear-remapping method
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

CCCCCC linear remapping of Lipscomb(2001)
*** setting flags
      DO IJ = IJTSTR, IJTEND
         AFLRMC(IJ) = AMSKT(IJ, KSTR)
         AFLRM(IJ, 0) = 1.0D0
      END DO
 
      DO K = 1, NIC
         DO IJ = IJTSTR, IJTEND
            IF ((AX(IJ, K) .GT. 0.0D0).AND.
     &          (AZ(IJ, K) .GT. 0.0D0) ) THEN
               AFLRM(IJ, K) = 1.0D0
            END IF
         END DO
      END DO

*** growth rate of each categories and category boundaries
      DO IJ = IJTSTR, IJTEND
         FDTN(IJ, 0) = HIX(IJ, 0)
      END DO
      DO K = 1, NIC
         DO IJ = IJTSTR, IJTEND
            FDTN(IJ, K) = ( HIX(IJ, K) - HIZ(IJ, K) ) * AFLRM(IJ, K)
         END DO
      END DO

      DO IJ = IJTSTR, IJTEND
         IF (FDTN(IJ, 0) .GT. 0.0D0) THEN
            FDTCN(IJ, 1) = FDTN(IJ, 0)
         ELSEIF (AFLRM(IJ, 1) .EQ. 1.0D0) THEN
            FDTCN(IJ, 1) = FDTN(IJ, 1)
         ELSE
            FDTCN(IJ, 1) = 0.0D0
         ENDIF
      END DO
      DO K = 2, NIC
         DO IJ = IJTSTR, IJTEND
            IF ((AFLRM(IJ, K-1) .EQ. 1.0D0)
     &          .AND.(AFLRM(IJ, K) .EQ. 1.0D0)) THEN
               FDTCN(IJ, K) = FDTN(IJ, K-1) +
     &              ( FDTN(IJ, K) - FDTN(IJ, K-1) ) *
     &              ( HIC0(K) - HIZ(IJ, K-1) )
     &              / ( HIZ(IJ, K) - HIZ(IJ, K-1) )
            ELSEIF (AFLRM(IJ, K-1) .EQ. 1.0D0) THEN
               FDTCN(IJ, K) = FDTN(IJ, K-1)
            ELSEIF (AFLRM(IJ, K) .EQ. 1.0D0) THEN
               FDTCN(IJ, K) = FDTN(IJ, K)
            ELSE
               FDTCN(IJ, K) = 0.0D0
            ENDIF
         END DO
      END DO

c      CALL CHEKIN(FDTN,  'FDTN',
c     &            NX,     NY,    NIC, NXYIDM, 'ICE')
c      CALL CHEKIN(FDTCN,   'FDTCN',
c     &            NX,     NY,    NIC, NXYIDM, 'ICE')
c      CALL CHEKIN(HIZ,  'HIZ',
c     &            NX,     NY,    NIC, NXYIDM, 'ICE')
c      CALL CHEKIN(HIX,  'HIX1',
c     &            NX,     NY,    NIC, NXYIDM, 'ICE')

*** temporally shift category boundaries
      DO K = 1, NIC
         DO IJ = IJTSTR, IJTEND
            HICN(IJ, K) = HIC(K) + FDTCN(IJ,K)
         END DO
      END DO

*** validation check: will not execute remapping when...
***  - HICN(IJ,K) does not lie between HIX(IJ,K-1) and HIX(IJ,K)
***  - HICN(IJ,K) does not lie between HIC0(K-1) and HIC0(K+1)
      DO K = 1, NIC
         DO IJ = IJTSTR, IJTEND
            HICN(IJ, K) = HIC(K) + FDTCN(IJ,K)
            IF ( (AFLRM(IJ, K-1) .EQ. 1.0D0) .AND.
     &           (HICN(IJ, K) .LE. HIX(IJ, K-1)) ) THEN
               AFLRMC(IJ) = 0.0D0
            ENDIF
            IF ( (AFLRM(IJ, K) .EQ. 1.0D0) .AND.
     &           (HICN(IJ, K) .GE. HIX(IJ, K)) ) THEN
               AFLRMC(IJ) = 0.0D0
            ENDIF
            IF ( (HICN(IJ, K) .LE. HIC0(K-1)) .OR.
     &           (HICN(IJ, K) .GE. HIC0(K+1)) ) THEN
               AFLRMC(IJ) = 0.0D0
            ENDIF
         END DO
      END DO

      DO K = 0, NIC
         DO IJ = IJTSTR, IJTEND
            AFLRM(IJ, K) = AFLRM(IJ, K) * AFLRMC(IJ)
         END DO
      END DO

c      CALL CHEKIN(HICN,  'HICN',
c     &            NX,     NY,    NIC, NXYIDM, 'ICE')

*** determine linear distribution function within each category
      DO K = 1, NIC-1
         DO IJ = IJTSTR, IJTEND
            HIL(IJ, K) = MAX( HICN(IJ, K),
     &                        3.0D0*HIX(IJ, K) - 2.0D0*HICN(IJ, K+1) )
            HIR(IJ, K) = MIN( HICN(IJ, K+1),
     &                        3.0D0*HIX(IJ, K) - 2.0D0*HICN(IJ, K) )
            ETAN = HIX(IJ, K) - HIL(IJ, K)
c           function MAX would not applied when AX(IJ, K) > 0
            ETAR = MAX(HIR(IJ, K) - HIL(IJ, K), EPSL)
            ETANRR = ETAN / ETAR
            G1(IJ, K) = 12.0D0 * AX(IJ, K) / (ETAR*ETAR)
     &                  * ( ETANRR - 0.5D0)
     &                  * AFLRM(IJ, K)
            G0(IJ, K) = 6.0D0 * AX(IJ, K) / ETAR
     &                  * ( 2.0D0 / 3.0D0 - ETANRR)
     &                  * AFLRM(IJ, K)
         END DO
      END DO
      DO IJ = IJTSTR, IJTEND
         HIL(IJ, NIC) = HICN(IJ, NIC)
         HIR(IJ, NIC) = 3.0D0*HIX(IJ, NIC) - 2.0D0*HICN(IJ, NIC)
c         HIR(IJ, NIC) = MAX(3.0D0*HIX(IJ, NIC) - 2.0D0*HICN(IJ, NIC),
c     &                      HIC(NIC))
         ETAN = HIX(IJ, NIC) - HIL(IJ, NIC)
c        function MAX would not applied when AX(IJ, NIC) > 0
         ETAR = MAX(HIR(IJ, NIC) - HIL(IJ, NIC), EPSL)
         ETANRR = ETAN / ETAR
         G1(IJ, NIC) = 12.0D0 * AX(IJ, NIC) / (ETAR*ETAR)
     &                 * ( ETANRR - 0.5D0 )
     &                 * AFLRM(IJ, NIC)
         G0(IJ, NIC) = 6.0D0 * AX(IJ, NIC) / ETAR
     &                 * ( 2.0D0 / 3.0D0 - ETANRR)
     &                 * AFLRM(IJ, NIC)
      END DO

c      CALL CHEKIN(HIL,  'HIL',
c     &            NX,     NY,    NIC, NXYIDM, 'ICE')
c      CALL CHEKIN(HIR,  'HIR',
c     &            NX,     NY,    NIC, NXYIDM, 'ICE')
c      CALL CHEKIN(G0,  'G0',
c     &            NX,     NY,    NIC, NXYIDM, 'ICE')
c      CALL CHEKIN(G1,  'G1',
c     &            NX,     NY,    NIC, NXYIDM, 'ICE')


*** compute transfer fluxes
*** area and volume fluxes are calculated from the linear distribution
*** snow and enthalpy fluxes are propotional to the volume flux
      DO K = 1, NIC
         DO IJ = IJTSTR, IJTEND
            LAXHIX(IJ, K) = AX(IJ, K) * HIX(IJ, K)
            LAXHSX(IJ, K) = AX(IJ, K) * HSX(IJ, K)
            LAXEIX(IJ, K) = AX(IJ, K) * EIX(IJ, K)
            DA(IJ, K) = 0.0D0
            DAHI(IJ, K) = 0.0D0
            DAHS(IJ, K) = 0.0D0
            DAEI(IJ, K) = 0.0D0
         END DO
      END DO

      DO IJ = IJTSTR, IJTEND
         DA(IJ, 0) = 0.0D0
c        for the flux across the lowest boundary,
c        consider only melting situation at the boundary
         IF ( ( HICN(IJ, 1) .LT. HIC(1) ) .AND.
     &        ( AFLRM(IJ, 1) .EQ. 1.0D0 ) ) THEN
c         IF ( HICN(IJ, 1) .LT. HIC(1) ) THEN
            IF (HIR(IJ, 1) .LE. HIC(1)) THEN
c              all sea ice transfered out of the lowest category
               GINT = AX(IJ, 1)
            ELSE
c               GIL = HIL(IJ, 1) - HIL(IJ, 1)
               GIL = 0.0D0
               GIR = MAX(HIC(1) - HIL(IJ, 1), GIL)
               X0 = GIR - GIL
               X1 = 0.5D0 * (GIR*GIR - GIL*GIL)
               GINT = G0(IJ, 1)*X0 + G1(IJ, 1)*X1
            END IF
c           area loss does not cause sea-ice to thicken
            GINT = MIN( GINT,
     &                  AX(IJ, 1)*(1.0D0 - HIX(IJ, 1)/HIZ(IJ, 1)) )
c           only area flux can across the lowest boundary
            DA(IJ, 1) = DA(IJ, 1) - GINT
            DA(IJ, 0) = DA(IJ, 0) + GINT
         END IF
      END DO

      DO K = 2, NIC
         DO IJ = IJTSTR, IJTEND
            IF ( ( HICN(IJ, K) .GE. HIC(K) ) .AND.
     &           ( AFLRM(IJ, K-1) .EQ. 1.0D0 ) ) THEN
c            IF ( HICN(IJ, K) .GE. HIC(K) ) THEN
               IF (HIL(IJ, K-1) .GE. HIC(K)) THEN
c                 all sea ice transfered to the upper category
                  GINT = AX(IJ, K-1)
                  FAHI = LAXHIX(IJ, K-1)
                  PVOL = 1.0D0
                  FAHS = LAXHSX(IJ, K-1)
                  FAEI = LAXEIX(IJ, K-1)
               ELSE
                  GIL = HIC(K) - HIL(IJ, K-1)
                  GIR = MAX(HIR(IJ, K-1) - HIL(IJ, K-1), GIL)
                  X0 = GIR - GIL
                  X1 = 0.5D0 * (GIR*GIR - GIL*GIL)
                  GINT = G0(IJ, K-1)*X0 + G1(IJ, K-1)*X1
                  X0 = 0.5D0 * (GIR*GIR - GIL*GIL)
                  X1 = (GIR*GIR*GIR - GIL*GIL*GIL) / 3.0D0
                  FAHI = HIL(IJ, K-1) * GINT +
     &                     G0(IJ, K-1)*X0 + G1(IJ, K-1)*X1
                  PVOL = FAHI / LAXHIX(IJ, K-1)
                  FAHS = LAXHSX(IJ, K-1) * PVOL
                  IF (TIX(IJ, K-1) .LT. TMI) THEN
                     FAEI = LAXEIX(IJ, K-1) * PVOL
c                     FAEI = FAHI * EI(TIX(IJ, K-1), SI)
                  ELSE
                     FAEI = 0.D0
                  END IF
c               FAEI = AX(IJ, K-1) * EIX(IJ, K-1) * PVOL
               END IF
               DA(IJ, K-1) = DA(IJ, K-1) - GINT
               DAHI(IJ, K-1) = DAHI(IJ, K-1) - FAHI
               DAHS(IJ, K-1) = DAHS(IJ, K-1) - FAHS
               DAEI(IJ, K-1) = DAEI(IJ, K-1) - FAEI
               DA(IJ, K) = DA(IJ, K) + GINT
               DAHI(IJ, K) = DAHI(IJ, K) + FAHI
               DAHS(IJ, K) = DAHS(IJ, K) + FAHS
               DAEI(IJ, K) = DAEI(IJ, K) + FAEI
            ELSEIF ( ( HICN(IJ, K) .LT. HIC(K) ) .AND.
     &               ( AFLRM(IJ, K) .EQ. 1.0D0 ) ) THEN
               IF (HIR(IJ, K) .LE. HIC(K)) THEN
c                 all sea ice transfered to the lower category
                  GINT = AX(IJ, K)
                  FAHI = LAXHIX(IJ, K)
                  PVOL = 1.0D0
                  FAHS = LAXHSX(IJ, K)
                  FAEI = LAXEIX(IJ, K)
               ELSE
c                  GIL = HIL(IJ, K) - HIL(IJ, K)
                  GIL = 0.0D0
                  GIR = MAX(HIC(K) - HIL(IJ, K), GIL)
                  X0 = GIR - GIL
                  X1 = 0.5D0 * (GIR*GIR - GIL*GIL)
                  GINT = G0(IJ, K)*X0 + G1(IJ, K)*X1
                  GINT = MIN(GINT, AX(IJ, K))
                  X0 = 0.5D0 * (GIR*GIR - GIL*GIL)
                  X1 = (GIR*GIR*GIR - GIL*GIL*GIL) / 3.0D0
                  FAHI = HIL(IJ, K) * GINT +
     &                     G0(IJ, K)*X0 + G1(IJ, K)*X1
                  FAHI = MIN(FAHI, LAXHIX(IJ, K))
                  PVOL = FAHI / LAXHIX(IJ, K)
                  FAHS = LAXHSX(IJ, K) * PVOL
                  IF (TIX(IJ, K) .LT. TMI) THEN
                     FAEI = LAXEIX(IJ, K) * PVOL
c                     FAEI = FAHI * EI(TIX(IJ, K), SI)
                  ELSE
                     FAEI = 0.D0
                  END IF
c                  FAEI = AX(IJ, K) * EIX(IJ, K) * PVOL
               ENDIF
               DA(IJ, K-1) = DA(IJ, K-1) + GINT
               DAHI(IJ, K-1) = DAHI(IJ, K-1) + FAHI
               DAHS(IJ, K-1) = DAHS(IJ, K-1) + FAHS
               DAEI(IJ, K-1) = DAEI(IJ, K-1) + FAEI
               DA(IJ, K) = DA(IJ, K) - GINT
               DAHI(IJ, K) = DAHI(IJ, K) - FAHI
               DAHS(IJ, K) = DAHS(IJ, K) - FAHS
               DAEI(IJ, K) = DAEI(IJ, K) - FAEI
            END IF
         END DO
      END DO

c      CALL CHEKIN(DA,  'DA',
c     &            NX,     NY,    NIC, NXYIDM, 'ICE')
c      CALL CHEKIN(DAHI,  'DAHI',
c     &            NX,     NY,    NIC, NXYIDM, 'ICE')
c      CALL CHEKIN(DAHS,  'DAHS',
c     &            NX,     NY,    NIC, NXYIDM, 'ICE')
c      CALL CHEKIN(DAEI,  'DAEI',
c     &            NX,     NY,    NIC, NXYIDM, 'ICE')

*** update prediction variables.
      DO K = 1, NIC
         DO IJ = IJTSTR, IJTEND
            AX(IJ, K) = AX(IJ, K) + DA(IJ, K) * AMSKT(IJ, KSTR)
            LAXHIX(IJ, K) = LAXHIX(IJ, K) +
     &                        DAHI(IJ, K) * AMSKT(IJ, KSTR)
            LAXHSX(IJ, K) = LAXHSX(IJ, K) +
     &                        DAHS(IJ, K) * AMSKT(IJ, KSTR)
            LAXEIX(IJ, K) = LAXEIX(IJ, K) +
     &                        DAEI(IJ, K) * AMSKT(IJ, KSTR)
            AXHSXN(IJ, K) = AXHSXN(IJ, K) +
     &                        DAHS(IJ, K) * AMSKT(IJ, KSTR)
            AXEIXN(IJ, K) = AXEIXN(IJ, K) +
     &                        DAEI(IJ, K) * AMSKT(IJ, KSTR)
            IF (LAXEIX(IJ, K) .LE. 0.D0) THEN
               AX(IJ, K) = 0.D0
               HIX(IJ, K) = HIC(K)
               HSX(IJ, K) = 0.D0
               TIX(IJ, K) = TMI
               EIX(IJ, K) = 0.D0
            ELSEIF (AX(IJ, K) .LE. 0.D0) THEN
               AX(IJ, K) = 0.D0
               HIX(IJ, K) = HIC(K)
               HSX(IJ, K) = 0.D0
               TIX(IJ, K) = TMI
               EIX(IJ, K) = 0.D0
            ELSE
               HIX(IJ, K) = LAXHIX(IJ, K) / AX(IJ, K)
               HSX(IJ, K) = LAXHSX(IJ, K) / AX(IJ, K)
               EIX(IJ, K) = LAXEIX(IJ, K) / AX(IJ, K)
c               TIX(IJ, K) = TI(EIX(IJ, K)/HIX(IJ, K), SI)
            ENDIF
         END DO
      END DO

*** check if variables are in valid range
c      DO K = 1, NIC
c         DO IJ = IJTSTR, IJTEND
c            IF (AX(IJ, K) .LT. 0.D0) THEN
c               CALL REWNML(IFPAR, JFPAR)
c               WRITE(JFPAR, *) '### NEGATIVE AREA (linrmp) ###',
c     &         IJ, K, AX(IJ, K)
c               stop
c            END IF
cc            IF (     (LAXHIX(IJ, K) .LT. 0.D0)
cc     &          .OR. (AXHSXN(IJ, K) .LT. 0.D0)) THEN
cc               CALL REWNML(IFPAR, JFPAR)
cc               WRITE(JFPAR, *) '### NEGATIVE ICE/SNOW (linrmp) ###',
cc     &         IJ, K, AX(IJ, K), LAXHIX(IJ, K), AXHSXN(IJ, K)
cc               stop
cc            END IF
cc            IF (AXEIXN(IJ, K) .LT. 0.D0) THEN
cc               CALL REWNML(IFPAR, JFPAR)
cc               WRITE(JFPAR, *) '### NEGATIVE ENTHALPY (linrmp) ###',
cc     &         IJ, K, AX(IJ, K), AXEIXN(IJ, K)
cc               stop
cc            END IF
c            IF (HIX(IJ, K) .LE. 0.D0) THEN
c               CALL REWNML(IFPAR, JFPAR)
c               WRITE(JFPAR, *) '### INVALID ICE THICKNESS (linrmp) ###',
c     &         IJ, K, AX(IJ, K), HIX(IJ, K)
c               stop
c            END IF
c         END DO
c      END DO

CCCCCC end of linear remapping of Lipscomb(2001)

*** lateral ice formation/melting
      DO K = 1, NIC
         DO IJ = IJTSTR, IJTEND
            LAXHIX(IJ, K) = AX(IJ, K) * HIX(IJ, K)
            LAXEIX(IJ, K) = LAXEIX(IJ, K)
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
c     &               + rhoo * WEN(IJ) )
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
c         AX1MAX = AZ(IJ, 0) + AZ(IJ, 1)
         AX1MAX = AZ(IJ, 0) + AZ(IJ, 1) + DA(IJ, 0) + DA(IJ, 1) 
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
