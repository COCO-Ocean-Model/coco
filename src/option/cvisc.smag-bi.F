      SUBROUTINE VSCVEL(
     M                 GX,     GY,     XX,     YY,
     I                 UY,     UX,     VY,     VX,
     I                 HX,     HY,
     I                AMV,   TAUX,   TAUY)
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Biharmonic smagorinsky viscosity
c
c  HISTORY
c     '03.05.13  H.Hasumi: from COCO3.4
c     '07.04.18  H.Hasumi: Formulation changed
c     '07.11.16  H.Hasumi: bug fix
c     '08.06.11  H.Hasumi: initial/final processing
c     '08.07.10  H.Hasumi: initial/final processing
c     '09.01.06  T.Suzuki: bug fix (FUZ/FVZ added to GX/GY in VSCVLB)
c
c ---------------------------------------------------------------------

#include "zocdim.F"
#include "zocgrd.F"
#include "zocmsk.F"
#include "zocphy.F"

      REAL*8      UY(NXYDIM, NZDIM),     UX(NXYDIM, NZDIM)
      REAL*8      VY(NXYDIM, NZDIM),     VX(NXYDIM, NZDIM)
      REAL*8      HY(NXYDIM)       ,     HX(NXYDIM)
      REAL*8     AMV(NXYDIM, NZDIM)
      REAL*8    TAUX(NXYDIM),   TAUY(NXYDIM)
      REAL*8      GX(NXYDIM, NZDIM),     GY(NXYDIM, NZDIM)
      REAL*8      XX(NXYDIM, NZDIM),     YY(NXYDIM, NZDIM)

      REAL*8     ASM(NXYDIM)
      LOGICAL OFIRST
      SAVE ASM, OFIRST

      REAL*8     FUX(NXYDIM, NZDIM),    FVX(NXYDIM, NZDIM)
      REAL*8     FUY(NXYDIM, NZDIM),    FVY(NXYDIM, NZDIM)
      REAL*8     FUZ(NXYDIM, NZDIM),    FVZ(NXYDIM, NZDIM)
      REAL*8      RZ(NXYDIM, NZDIM),    RZM(NXYDIM, NZDIM)
      REAL*8     SXX(NXYDIM),    SYY(NXYDIM)
      REAL*8     SXY(NXYDIM)
      REAL*8     SZX(NXYDIM, NZDIM),    SZY(NXYDIM, NZDIM)
      COMMON /WORK/ FUX, FUY, FUZ, FVX, FVY, FVZ, RZ, RZM
      EQUIVALENCE(FUZ(1, 1), SZX(1, 1))
      EQUIVALENCE(FVZ(1, 1), SZY(1, 1))

      REAL*8   HVBOT(NXYDIM)
      REAL*8     EXX,    EYY,    EXY
      REAL*8      EZ
      REAL*8      DD(NXYDIM, NZDIM)
      REAL*8      HU(NXYDIM, NZDIM),     HV(NXYDIM, NZDIM)
      REAL*8      PI
      INTEGER     IJ,      K
      INTEGER   IJLN,   IJLS,   IJLE,   IJLW
      INTEGER  IJLNE,  IJLSW
      INTEGER  IFPAR,  JFPAR

      REAL*8    CSMB
      NAMELIST /NMSMAB/ CSMB

      DATA OFIRST / .TRUE. /
      DATA CSMB / 0.D0 /

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF

      IF (OFIRST) THEN
         OFIRST = .FALSE.
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMSMAB, END=99)
 99      CONTINUE
         WRITE(JFPAR, NMSMAB)

         PI = 4.D0 * ATAN(1.D0)
         DO IJ = 1, NXYDIM
            ASM(IJ) = CSMB * MIN(DX*HXT(IJ), DY(IJ)*HYT(IJ))**2
     &                / PI / SQRT(8.D0)
         END DO
      END IF

      DO 110 K = 1, NZDIM
         DO 100 IJ = 1, NXYDIM
            FUX(IJ, K) = 0.D0
            FVX(IJ, K) = 0.D0
            FUY(IJ, K) = 0.D0
            FVY(IJ, K) = 0.D0
            FUZ(IJ, K) = 0.D0
            FVZ(IJ, K) = 0.D0
  100    CONTINUE
  110 CONTINUE
      DO 120 IJ = IJVSTR, IJVEND
         FUZ(IJ, KSTR) = TAUX(IJ)
         FVZ(IJ, KSTR) = TAUY(IJ)
  120 CONTINUE
      DO IJ = IJVSTR, IJVEND
         HVBOT(IJ) = (  (HX(IJ)    + HX(IJ+LE) ) * DY(IJ)
     &                + (HX(IJ+LN) + HX(IJ+LNE)) * DY(IJ+LN)) *
     &               RYM(IJ) * 0.25D0
     &             + ZBOT
         HVBOT(IJ) = 1.D0 / HVBOT(IJ)
      END DO
      DO K = KSTR, KSTR+KZ-1
         DO IJ = IJVSTR, IJVEND
            RZ (IJ, K) = 1.D0 * RS (K) * HVBOT(IJ)
            RZM(IJ, K) = 1.D0 * RSM(K) * HVBOT(IJ)
         END DO
      END DO
      DO K = KSTR+KZ, KEND
         DO IJ = IJVSTR, IJVEND
            RZ (IJ, K) = 1.D0 / DZV(IJ, K)
            RZM(IJ, K) = 1.D0 / DZM(IJ, K)
         END DO
      END DO

      DO 210 K = KSTR+1, KEND
         DO 200 IJ = IJVSTR, IJVEND
            FUZ(IJ, K) = AMV(IJ, K) * RZM(IJ, K) *
     &                   (UX(IJ, K-1) - UX(IJ, K))
            FVZ(IJ, K) = AMV(IJ, K) * RZM(IJ, K) *
     &                   (VX(IJ, K-1) - VX(IJ, K))
  200    CONTINUE
  210 CONTINUE
      DO K = KSTR, KEND
         DO IJ = IJVSTR, IJVEND
            GX(IJ, K) = GX(IJ, K)
     &                + (FUZ(IJ, K) - FUZ(IJ, K+1)) * RZ(IJ, K)
            GY(IJ, K) = GY(IJ, K)
     &                + (FVZ(IJ, K) - FVZ(IJ, K+1)) * RZ(IJ, K)
         END DO
      END DO

      DO K = KSTR, KEND
         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+NXDIM+2
            IJLS = IJ + LS
            IJLW = IJ + LW
            IJLSW = IJ + LSW
            EXX = (  UX(IJ  , K) + UX(IJLS , K)
     &             - UX(IJLW, K) - UX(IJLSW, K)) * 0.5D0 *
     &            RX * RXT(IJ)
     &          + (  VX(IJ  , K) + VX(IJLS , K)
     &             + VX(IJLW, K) + VX(IJLSW, K)) * 0.25D0 *
     &            HXYT(IJ)
            EYY = (  VX(IJ  , K) + VX(IJLW , K)
     &             - VX(IJLS, K) - VX(IJLSW, K)) * 0.5D0 *
     &            RY(IJ) * RYT(IJ)
     &          + (  UX(IJ  , K) + UX(IJLW , K)
     &             + UX(IJLS, K) + UX(IJLSW, K)) * 0.25D0 *
     &            HYXT(IJ)
            EXY = (  UX(IJ  , K) + UX(IJLW , K)
     &             - UX(IJLS, K) - UX(IJLSW, K)) * 0.5D0 *
     &            RY(IJ) * RYT(IJ)
     &          + (  VX(IJ  , K) + VX(IJLS , K)
     &             - VX(IJLW, K) - VX(IJLSW, K)) * 0.5D0 *
     &            RX * RXT(IJ)
     &          - (  UX(IJ  , K) + UX(IJLW , K)
     &             + UX(IJLS, K) + UX(IJLSW, K)) * 0.25D0 *
     &            HXYT(IJ)
     &          - (  VX(IJ  , K) + VX(IJLW , K)
     &             + VX(IJLS, K) + VX(IJLSW, K)) * 0.25D0 *
     &            HYXT(IJ)

            DD(IJ, K) = SQRT(SQRT(
     &                  (EXX - EYY) * (EXX - EYY) + EXY * EXY))

            SXX(IJ) = ASM(IJ) * DD(IJ, K) * (EXX - EYY)
            SYY(IJ) = ASM(IJ) * DD(IJ, K) * (EYY - EXX)
            SXY(IJ) = ASM(IJ) * DD(IJ, K) * EXY
         END DO

         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+2
            IJLN = IJ + LN
            FUX(IJ, K) = (  SXX(IJ) * HYT(IJ) * HYT(IJ)
     &                    + SXX(IJLN) * HYT(IJLN) * HYT(IJLN)) *
     &                   0.5D0
            FVX(IJ, K) = (  SXY(IJ) * HYT(IJ) * HYT(IJ)
     &                    + SXY(IJLN) * HYT(IJLN) * HYT(IJLN)) *
     &                   0.5D0
         END DO
         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+NXDIM+1
            IJLE = IJ + LE
            FUY(IJ, K) = (  SXY(IJ) * HXT(IJ) * HXT(IJ)
     &                    + SXY(IJLE) * HXT(IJLE) * HXT(IJLE)) *
     &                   0.5D0
            FVY(IJ, K) = (  SYY(IJ) * HXT(IJ) * HXT(IJ)
     &                    + SYY(IJLE) * HXT(IJLE) * HXT(IJLE)) *
     &                   0.5D0
         END DO
      END DO

      DO K = KSTR, KEND
         DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
            HU(IJ, K) = - (  (FUX(IJ+LE, K) - FUX(IJ, K)) *
     &                       RX * RYU(IJ)
     &                     + (FUY(IJ+LN, K) - FUY(IJ, K)) *
     &                       RYM(IJ) * RXU(IJ)) *
     &                    RXU(IJ) * RYU(IJ) * AMSKV(IJ, K)
            HV(IJ, K) = - (  (FVX(IJ+LE, K) - FVX(IJ, K)) *
     &                       RX * RYU(IJ)
     &                     + (FVY(IJ+LN, K) - FVY(IJ, K)) *
     &                       RYM(IJ) * RXU(IJ)) *
     &                    RXU(IJ) * RYU(IJ) * AMSKV(IJ, K)
         END DO
      END DO

      DO K = KSTR, KEND
         DO IJ = IJSTR, IJEND+NXDIM+1
            IJLS = IJ + LS
            IJLW = IJ + LW
            IJLSW = IJ + LSW
            EXX = (  HU(IJ  , K) + HU(IJLS , K)
     &             - HU(IJLW, K) - HU(IJLSW, K)) * 0.5D0 *
     &            RX * RXT(IJ)
     &          + (  HV(IJ  , K) + HV(IJLS , K)
     &             + HV(IJLW, K) + HV(IJLSW, K)) * 0.25D0 *
     &            HXYT(IJ)
            EYY = (  HV(IJ  , K) + HV(IJLW , K)
     &             - HV(IJLS, K) - HV(IJLSW, K)) * 0.5D0 *
     &            RY(IJ) * RYT(IJ)
     &          + (  HU(IJ  , K) + HU(IJLW , K)
     &             + HU(IJLS, K) + HU(IJLSW, K)) * 0.25D0 *
     &            HYXT(IJ)
            EXY = (  HU(IJ  , K) + HU(IJLW , K)
     &             - HU(IJLS, K) - HU(IJLSW, K)) * 0.5D0 *
     &            RY(IJ) * RYT(IJ)
     &          + (  HV(IJ  , K) + HV(IJLS , K)
     &             - HV(IJLW, K) - HV(IJLSW, K)) * 0.5D0 *
     &            RX * RXT(IJ)
     &          - (  HU(IJ  , K) + HU(IJLW , K)
     &             + HU(IJLS, K) + HU(IJLSW, K)) * 0.25D0 *
     &            HXYT(IJ)
     &          - (  HV(IJ  , K) + HV(IJLW , K)
     &             + HV(IJLS, K) + HV(IJLSW, K)) * 0.25D0 *
     &            HYXT(IJ)

            SXX(IJ) = ASM(IJ) * DD(IJ, K) * (EXX - EYY)
            SYY(IJ) = ASM(IJ) * DD(IJ, K) * (EYY - EXX)
            SXY(IJ) = ASM(IJ) * DD(IJ, K) * EXY
         END DO
         DO IJ = IJVSTR, IJVEND+1
            IJLN = IJ + LN
            FUX(IJ, K) = (  SXX(IJ) * HYT(IJ) * HYT(IJ)
     &                    + SXX(IJLN) * HYT(IJLN) * HYT(IJLN)) *
     &                   0.5D0 * AMFVX(IJ, K)
            FVX(IJ, K) = (  SXY(IJ) * HYT(IJ) * HYT(IJ)
     &                    + SXY(IJLN) * HYT(IJLN) * HYT(IJLN)) *
     &                   0.5D0 * AMFVX(IJ, K)
         END DO
         DO IJ = IJVSTR, IJVEND+NXDIM
            IJLE = IJ + LE
            FUY(IJ, K) = (  SXY(IJ) * HXT(IJ) * HXT(IJ)
     &                    + SXY(IJLE) * HXT(IJLE) * HXT(IJLE)) *
     &                   0.5D0 * AMFVY(IJ, K)
            FVY(IJ, K) = (  SYY(IJ) * HXT(IJ) * HXT(IJ)
     &                    + SYY(IJLE) * HXT(IJLE) * HXT(IJLE)) *
     &                   0.5D0 * AMFVY(IJ, K)
         END DO
      END DO

      DO K = KSTR, KEND
         DO IJ = IJVSTR, IJVEND
            EZ = (  (AMV(IJ, K) + AMV(IJ, K+1)) * 0.5D0 / REA
     &            + (AMV(IJ, K) - AMV(IJ, K+1)) * RZ(IJ, K))
            SZX(IJ, K) = - UX(IJ, K) * EZ
            SZY(IJ, K) = - VX(IJ, K) * EZ
         END DO
      END DO

      DO 510 K = KSTR, KSTR+KZ-1
         DO 500 IJ = IJVSTR, IJVEND
            GX(IJ, K) = (  GX(IJ, K)
     &                   + (  (FUX(IJ+LE, K) - FUX(IJ, K)) *
     &                        RX * RYU(IJ)
     &                      + (FUY(IJ+LN, K) - FUY(IJ, K)) *
     &                        RYM(IJ) * RXU(IJ)) * RXU(IJ) * RYU(IJ)
     &                   + SZX(IJ, K) / REA) * AMSKV(IJ, K)
            GY(IJ, K) = (  GY(IJ, K)
     &                   + (  (FVX(IJ+LE, K) - FVX(IJ, K)) *
     &                        RX * RYU(IJ)
     &                      + (FVY(IJ+LN, K) - FVY(IJ, K)) *
     &                        RYM(IJ) * RXU(IJ)) * RXU(IJ) * RYU(IJ)
     &                   + SZY(IJ, K) / REA) * AMSKV(IJ, K)
  500    CONTINUE
  510 CONTINUE
      DO 530 K = KSTR+KZ, KEND
         DO 520 IJ = IJVSTR, IJVEND
            GX(IJ, K) = (  GX(IJ, K)
     &                   + (  (FUX(IJ+LE, K) - FUX(IJ, K)) *
     &                        RX * RYU(IJ)
     &                      + (FUY(IJ+LN, K) - FUY(IJ, K)) *
     &                        RYM(IJ) * RXU(IJ)
     &                     ) * RXU(IJ) * RYU(IJ) * RZ(IJ, K)
     &                   + SZX(IJ, K) / REA) * AMSKV(IJ, K)
            GY(IJ, K) = (  GY(IJ, K)
     &                   + (  (FVX(IJ+LE, K) - FVX(IJ, K)) *
     &                        RX * RYU(IJ)
     &                      + (FVY(IJ+LN, K) - FVY(IJ, K)) *
     &                        RYM(IJ) * RXU(IJ)
     &                     ) * RXU(IJ) * RYU(IJ) * RZ(IJ, K)
     &                   + SZY(IJ, K) / REA) * AMSKV(IJ, K)
  520    CONTINUE
  530 CONTINUE

      DO 610 K = KSTR, KEND
         DO 600 IJ = IJVSTR, IJVEND
            XX(IJ, K) = GX(IJ, K)
            YY(IJ, K) = GY(IJ, K)
 600     CONTINUE
 610  CONTINUE

      RETURN
      END
#ifdef OPT_BBL
c *********************************************************************

      SUBROUTINE VSCVLB(
     M                 GX,     GY,     XX,     YY,
     I                 UY,     UX,     VY,     VX,
     I                AMV)
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Visocity term of the BBL momentum eqs.
c
c  HISTORY
c     '01.12.03  H.Hasumi: from cflxvb.F of COCO3.3.05
c
c ---------------------------------------------------------------------

#include "zocdim.F"
#include "zocgrd.F"
#include "zocmsk.F"
#include "zocphy.F"

      REAL*8      UY(NXYDIM, NZDIM),     UX(NXYDIM, NZDIM)
      REAL*8      VY(NXYDIM, NZDIM),     VX(NXYDIM, NZDIM)
      REAL*8     AMV(NXYDIM, NZDIM)
      REAL*8      GX(NXYDIM, NZDIM),     GY(NXYDIM, NZDIM)
      REAL*8      XX(NXYDIM, NZDIM),     YY(NXYDIM, NZDIM)

      REAL*8      RZ(NXYDIM),    RZM(NXYDIM)
      REAL*8     ASM(NXYDIM)
      LOGICAL OFIRST
      SAVE RZ, RZM, ASM

      REAL*8     FUX(NXYDIM),    FVX(NXYDIM)
      REAL*8     FUY(NXYDIM),    FVY(NXYDIM)
      REAL*8     FUZ(NXYDIM),    FVZ(NXYDIM)
      REAL*8     SXX(NXYDIM),    SYY(NXYDIM)
      REAL*8     SXY(NXYDIM)
      REAL*8     SZX(NXYDIM),    SZY(NXYDIM)
      REAL*8      HU(NXYDIM),     HV(NXYDIM)
      REAL*8      DD(NXYDIM)
      COMMON /WORK/ FUX, FUY, FUZ, FVX, FVY, FVZ,
     &              SXX, SYY, SXY, SZX, SZY,
     &               HU,  HV,  DD

      REAL*8     EXX,    EYY,    EXY
      REAL*8      PI
      INTEGER     IJ,    KUP,      K
      INTEGER   IJLN,   IJLS,   IJLE,   IJLW
      INTEGER  IJLNE,  IJLSW
      INTEGER  IFPAR,  JFPAR

      REAL*8     CSMBBB
      NAMELIST /NMSMBB/ CSMBBB

      DATA OFIRST / .TRUE. /
      DATA CSMBBB / 0.D0 /

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF

      IF (OFIRST) THEN
         OFIRST = .FALSE.
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMSMBB, END=99)
 99      CONTINUE
         WRITE(JFPAR, NMSMBB)

         PI = 4.D0 * ATAN(1.D0)
         DO IJ = 1, NXYDIM
            ASM(IJ) = CSMBBB * MIN(DX*HXT(IJ), DY(IJ)*HYT(IJ))**2
     &                / PI / SQRT(8.D0)
            RZ (IJ) = 1.D0 / DZV(IJ, KEND)
            RZM(IJ) = 1.D0 / DZM(IJ, KEND)
         END DO
      END IF

      K = KEND

      DO IJ = IJVSTR, IJVEND
         KUP = MAX(NBOTV(IJ)-1, 1)
         FUZ(IJ) = AMV(IJ, KEND) * RZM(IJ) *
     &             (UX(IJ, KUP) - UX(IJ, KEND))
         FVZ(IJ) = AMV(IJ, KEND) * RZM(IJ) *
     &             (VX(IJ, KUP) - VX(IJ, KEND))
      END DO

      DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+NXDIM+2
         IJLS = IJ + LS
         IJLW = IJ + LW
         IJLSW = IJ + LSW
         EXX = (  UX(IJ  , K) + UX(IJLS , K)
     &          - UX(IJLW, K) - UX(IJLSW, K)) * 0.5D0 *
     &         RX * RXT(IJ)
     &       + (  VX(IJ  , K) + VX(IJLS , K)
     &          + VX(IJLW, K) + VX(IJLSW, K)) * 0.25D0 *
     &         HXYT(IJ)
         EYY = (  VX(IJ  , K) + VX(IJLW , K)
     &          - VX(IJLS, K) - VX(IJLSW, K)) * 0.5D0 *
     &         RY(IJ) * RYT(IJ)
     &       + (  UX(IJ  , K) + UX(IJLW , K)
     &          + UX(IJLS, K) + UX(IJLSW, K)) * 0.25D0 *
     &         HYXT(IJ)
         EXY = (  UX(IJ  , K) + UX(IJLW , K)
     &          - UX(IJLS, K) - UX(IJLSW, K)) * 0.5D0 *
     &         RY(IJ) * RYT(IJ)
     &       + (  VX(IJ  , K) + VX(IJLS , K)
     &          - VX(IJLW, K) - VX(IJLSW, K)) * 0.5D0 *
     &         RX * RXT(IJ)
     &       - (  UX(IJ  , K) + UX(IJLW , K)
     &          + UX(IJLS, K) + UX(IJLSW, K)) * 0.25D0 *
     &         HXYT(IJ)
     &       - (  VX(IJ  , K) + VX(IJLW , K)
     &          + VX(IJLS, K) + VX(IJLSW, K)) * 0.25D0 *
     &         HYXT(IJ)

         DD(IJ) = SQRT(SQRT((EXX - EYY) * (EXX - EYY) + EXY * EXY))

         SXX(IJ) = ASM(IJ) * DD(IJ) * (EXX - EYY)
         SYY(IJ) = ASM(IJ) * DD(IJ) * (EYY - EXX)
         SXY(IJ) = ASM(IJ) * DD(IJ) * EXY
      END DO

      DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+2
         IJLN = IJ + LN
         FUX(IJ) = (  SXX(IJ) * HYT(IJ) * HYT(IJ)
     &              + SXX(IJLN) * HYT(IJLN) * HYT(IJLN)) * 0.5D0
         FVX(IJ) = (  SXY(IJ) * HYT(IJ) * HYT(IJ)
     &              + SXY(IJLN) * HYT(IJLN) * HYT(IJLN)) * 0.5D0
      END DO
      DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+NXDIM+1
         IJLE = IJ + LE
         FUY(IJ) = (  SXY(IJ) * HXT(IJ) * HXT(IJ)
     &              + SXY(IJLE) * HXT(IJLE) * HXT(IJLE)) * 0.5D0
         FVY(IJ) = (  SYY(IJ) * HXT(IJ) * HXT(IJ)
     &              + SYY(IJLE) * HXT(IJLE) * HXT(IJLE)) * 0.5D0
      END DO

      DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
         HU(IJ) = - (  (FUX(IJ+LE) - FUX(IJ)) *
     &                 RX * RYU(IJ)
     &               + (FUY(IJ+LN) - FUY(IJ)) *
     &                 RYM(IJ) * RXU(IJ)) *
     &              RXU(IJ) * RYU(IJ) * AMSKV(IJ, K)
         HV(IJ) = - (  (FVX(IJ+LE) - FVX(IJ)) *
     &                 RX * RYU(IJ)
     &               + (FVY(IJ+LN) - FVY(IJ)) *
     &                 RYM(IJ) * RXU(IJ)) *
     &              RXU(IJ) * RYU(IJ) * AMSKV(IJ, K)
      END DO

      DO IJ = IJSTR, IJEND+NXDIM+1
         IJLS = IJ + LS
         IJLW = IJ + LW
         IJLSW = IJ + LSW
         EXX = (  HU(IJ  ) + HU(IJLS )
     &          - HU(IJLW) - HU(IJLSW)) * 0.5D0 *
     &         RX * RXT(IJ)
     &       + (  HV(IJ  ) + HV(IJLS )
     &          + HV(IJLW) + HV(IJLSW)) * 0.25D0 *
     &         HXYT(IJ)
         EYY = (  HV(IJ  ) + HV(IJLW )
     &          - HV(IJLS) - HV(IJLSW)) * 0.5D0 *
     &         RY(IJ) * RYT(IJ)
     &       + (  HU(IJ  ) + HU(IJLW )
     &          + HU(IJLS) + HU(IJLSW)) * 0.25D0 *
     &         HYXT(IJ)
         EXY = (  HU(IJ  ) + HU(IJLW )
     &          - HU(IJLS) - HU(IJLSW)) * 0.5D0 *
     &         RY(IJ) * RYT(IJ)
     &       + (  HV(IJ  ) + HV(IJLS )
     &          - HV(IJLW) - HV(IJLSW)) * 0.5D0 *
     &         RX * RXT(IJ)
     &       - (  HU(IJ  ) + HU(IJLW )
     &          + HU(IJLS) + HU(IJLSW)) * 0.25D0 *
     &         HXYT(IJ)
     &       - (  HV(IJ  ) + HV(IJLW )
     &          + HV(IJLS) + HV(IJLSW)) * 0.25D0 *
     &         HYXT(IJ)

         SXX(IJ) = ASM(IJ) * DD(IJ) * (EXX - EYY)
         SYY(IJ) = ASM(IJ) * DD(IJ) * (EYY - EXX)
         SXY(IJ) = ASM(IJ) * DD(IJ) * EXY
      END DO
      DO IJ = IJVSTR, IJVEND+1
         IJLN = IJ + LN
         FUX(IJ) = (  SXX(IJ) * HYT(IJ) * HYT(IJ)
     &              + SXX(IJLN) * HYT(IJLN) * HYT(IJLN)) *
     &             0.5D0 * AMFVX(IJ, K)
         FVX(IJ) = (  SXY(IJ) * HYT(IJ) * HYT(IJ)
     &              + SXY(IJLN) * HYT(IJLN) * HYT(IJLN)) *
     &             0.5D0 * AMFVX(IJ, K)
      END DO
      DO IJ = IJVSTR, IJVEND+NXDIM
         IJLE = IJ + LE
         FUY(IJ) = (  SXY(IJ) * HXT(IJ) * HXT(IJ)
     &              + SXY(IJLE) * HXT(IJLE) * HXT(IJLE)) *
     &             0.5D0 * AMFVY(IJ, K)
         FVY(IJ) = (  SYY(IJ) * HXT(IJ) * HXT(IJ)
     &              + SYY(IJLE) * HXT(IJLE) * HXT(IJLE)) *
     &             0.5D0 * AMFVY(IJ, K)
      END DO

      DO IJ = IJVSTR, IJVEND
         SZX(IJ) = - UX(IJ, K) * AMV(IJ, K) / REA
         SZY(IJ) = - VX(IJ, K) * AMV(IJ, K) / REA
      END DO

c     '09.01.20 ???
c      DO IJ = IJVSTR, IJVEND
c         GX(IJ, K) = (  GX(IJ, K)
c     &                + (  (FUX(IJ+LE) - FUX(IJ)) *
c     &                     RX * RYU(IJ)
c     &                   + (FUY(IJ+LN) - FUY(IJ)) *
c     &                     RYM(IJ) * RXU(IJ)) * RXU(IJ) * RYU(IJ)
c     &                + SZX(IJ) / REA) * AMSKV(IJ, K)
c         GY(IJ, K) = (  GY(IJ, K)
c     &                + (  (FVX(IJ+LE) - FVX(IJ)) *
c     &                     RX * RYU(IJ)
c     &                   + (FVY(IJ+LN) - FVY(IJ)) *
c     &                     RYM(IJ) * RXU(IJ)) * RXU(IJ) * RYU(IJ)
c     &                + SZY(IJ) / REA) * AMSKV(IJ, K)
      END DO
      DO IJ = IJVSTR, IJVEND
         GX(IJ, K) = (  GX(IJ, K)
     &                + FUZ(IJ) * RZ(IJ)
     &                + (  (FUX(IJ+LE) - FUX(IJ)) *
     &                     RX * RYU(IJ)
     &                   + (FUY(IJ+LN) - FUY(IJ)) *
     &                     RYM(IJ) * RXU(IJ)
     &                  ) * RXU(IJ) * RYU(IJ) * RZ(IJ)
     &                + SZX(IJ) / REA) * AMSKV(IJ, K)
         GY(IJ, K) = (  GY(IJ, K)
     &                + FUZ(IJ) * RZ(IJ)
     &                + (  (FVX(IJ+LE) - FVX(IJ)) *
     &                     RX * RYU(IJ)
     &                   + (FVY(IJ+LN) - FVY(IJ)) *
     &                     RYM(IJ) * RXU(IJ)
     &                  ) * RXU(IJ) * RYU(IJ) * RZ(IJ)
     &                + SZY(IJ) / REA) * AMSKV(IJ, K)
      END DO

      DO IJ = IJVSTR, IJVEND
         XX(IJ, KEND) = GX(IJ, KEND)
         YY(IJ, KEND) = GY(IJ, KEND)
      END DO

      RETURN
      END
#endif
