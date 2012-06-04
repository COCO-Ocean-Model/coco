      SUBROUTINE TNDUVD(
     O                GXX,    GYY,
     M                 GX,     GY,
     I                 XX,     YY,
     I                 UY,     VY,     HY,      R,
     I                 UX,     VX,
     I                FUX,    FUY,   FUNE,   FUSE,
     I                FVX,    FVY,   FVNE,   FVSE)
      use zocdim
      use zocgrd
      use zocmsk
      use zocphy
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Estimate the pressure gradient and Coriolis terms of the equation
c of motion.
c
c  HISTORY
c     '03.05.13  H.Hasumi: from COCO3.4
c     '07.04.23  H.Hasumi
c     '07.05.01  H.Hasumi: for McDougall et al. (2003) eq. of state
c     '08.07.10  H.Hasumi: initial/final processing
c     '09.01.15  T.Suzuki: bug fix vertical integration thanks to Hasumi
c
c ---------------------------------------------------------------------

      REAL*8       R(NXYDIM, NZDIM),     HY(NXYDIM)
      REAL*8      UY(NXYDIM, NZDIM),     VY(NXYDIM, NZDIM)
      REAL*8      UX(NXYDIM, NZDIM),     VX(NXYDIM, NZDIM)
      REAL*8      GX(NXYDIM, NZDIM),    GXX(NXYDIM)
      REAL*8      GY(NXYDIM, NZDIM),    GYY(NXYDIM)
      REAL*8      XX(NXYDIM, NZDIM),     YY(NXYDIM, NZDIM)
      REAL*8     FUX(NXYDIM, NZDIM),    FVX(NXYDIM, NZDIM)
      REAL*8     FUY(NXYDIM, NZDIM),    FVY(NXYDIM, NZDIM)
      REAL*8    FUNE(NXYDIM, NZDIM),   FVNE(NXYDIM, NZDIM)
      REAL*8    FUSE(NXYDIM, NZDIM),   FVSE(NXYDIM, NZDIM)

      REAL*8     COF
      SAVE COF

      REAL*8      PX(NXYDIM, NZDIM),    PXM(NXYDIM, NZDIM)
      REAL*8      PY(NXYDIM, NZDIM),    PYM(NXYDIM, NZDIM)
      REAL*8      UU(NXYDIM)       ,     UV(NXYDIM)
      REAL*8      VV(NXYDIM)
      REAL*8   HVBOT(NXYDIM)
      REAL*8   DZSIG(NXYDIM, NZDIM)
c      COMMON /WORK/ PX, PXM,  PY, PYM, UU, UV, VV,
c     &              HVBOT, DZSIG
      COMMON /WORK2/ PX, PXM,  PY, PYM

      INTEGER     IJ,      K
      INTEGER  IFPAR,  JFPAR

      LOGICAL OFIRST
      SAVE OFIRST
      DATA OFIRST / .TRUE. /

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF

      IF (OFIRST) THEN
         OFIRST = .FALSE.
         COF = -0.25D0 * gravit / rhoo * 1.D-3
      END IF

      DO IJ = IJVSTR, IJVEND+NXDIM+1
         HVBOT(IJ) = (  (HY(IJ) + HY(IJ+LE)) * DY(IJ)
     &                + (HY(IJ+LN) + HY(IJ+LNE)) * DY(IJ+LN)) *
     &               RYM(IJ) * 0.25D0
     &             + ZBOT
      END DO
      DO K = KSTR, KSTR+KZ-1
         DO IJ = IJVSTR, IJVEND+NXDIM+1
            DZSIG(IJ, K) = DS(K) * HVBOT(IJ)
         END DO
      END DO

      DO 100 IJ = IJVSTR, IJVEND
         PXM(IJ, KSTR) = (  R(IJ+LNE, KSTR) + R(IJ+LE, KSTR)
     &                    - R(IJ+LN , KSTR) - R(IJ   , KSTR)) *
     &                   COF * RX * RXU(IJ) * DZSIG(IJ, KSTR) *
     &                   AMSKV(IJ, KSTR)
         PX (IJ, KSTR) = PXM(IJ, KSTR)
         GX (IJ, KSTR) = GX(IJ, KSTR)
     &                 + (COR(IJ) * VX(IJ, KSTR) + PX(IJ, KSTR)) *
     &                   AMSKV(IJ, KSTR)
         XX (IJ, KSTR) = XX(IJ, KSTR)
     &                 + PX(IJ, KSTR) * AMSKV(IJ, KSTR)
         GXX(IJ)       = XX(IJ, KSTR) * DZSIG(IJ, KSTR) *
     &                   AMSKV(IJ, KSTR)
         PYM(IJ, KSTR) = (  R(IJ+LNE, KSTR) + R(IJ+LN, KSTR)
     &                    - R(IJ+LE , KSTR) - R(IJ   , KSTR)) *
     &                   COF * RYM(IJ) * RYU(IJ) * DZSIG(IJ, KSTR) *
     &                   AMSKV(IJ, KSTR)
         PY (IJ, KSTR) = PYM(IJ, KSTR)
         GY (IJ, KSTR) = GY(IJ, KSTR)
     &                 - (COR(IJ) * UX(IJ, KSTR) - PY(IJ, KSTR)) *
     &                   AMSKV(IJ, KSTR)
         YY (IJ, KSTR) = YY(IJ, KSTR)
     &                 + PY(IJ, KSTR) * AMSKV(IJ, KSTR)
         GYY(IJ)       = YY(IJ, KSTR) * DZSIG(IJ, KSTR) *
     &                   AMSKV(IJ, KSTR)
  100 CONTINUE

      DO 210 K = KSTR+1, KSTR+KZ-1
         DO 200 IJ = IJVSTR, IJVEND
            PXM(IJ, K) = (  R(IJ+LNE, K) + R(IJ+LE, K)
     &                    - R(IJ+LN , K) - R(IJ   , K)) *
     &                   COF * RX * RXU(IJ) * DZSIG(IJ, K) *
     &                   AMSKV(IJ, K)
            PX (IJ, K) = (PX(IJ, K-1) + PXM(IJ, K-1) + PXM(IJ, K)) *
     &                   AMSKV(IJ, K)
            GX (IJ, K) = GX(IJ, K)
     &                 + (COR(IJ) * VX(IJ, K) + PX(IJ, K)) *
     &                   AMSKV(IJ, K)
            XX (IJ, K) = XX(IJ, K)
     &                 + PX(IJ, K) * AMSKV(IJ, K)
            GXX(IJ)    = GXX(IJ)
     &                 + XX(IJ, K) * DZSIG(IJ, K) * AMSKV(IJ, K)
            PYM(IJ, K) = (  R(IJ+LNE, K) + R(IJ+LN, K)
     &                    - R(IJ+LE , K) - R(IJ   , K)) *
     &                   COF * RYM(IJ) * RYU(IJ) * DZSIG(IJ, K) *
     &                   AMSKV(IJ, K)
            PY (IJ, K) = (PY(IJ, K-1) + PYM(IJ, K-1) + PYM(IJ, K)) *
     &                   AMSKV(IJ, K)
            GY (IJ, K) = GY(IJ, K)
     &                 - (COR(IJ) * UX(IJ, K) - PY(IJ, K)) *
     &                   AMSKV(IJ, K)
            YY (IJ, K) = YY(IJ, K)
     &                 + PY(IJ, K) * AMSKV(IJ, K)
            GYY(IJ)    = GYY(IJ) + YY(IJ, K) * DZSIG(IJ, K) *
     &                   AMSKV(IJ, K)
  200    CONTINUE
  210 CONTINUE

      DO 240 K = KSTR+KZ, KEND
         DO 230 IJ = IJVSTR, IJVEND
            PXM(IJ, K) = (  R(IJ+LNE, K) - R(IJ+LN, K)
     &                    + R(IJ+LE , K) - R(IJ   , K)) *
     &                   COF * RX * RXU(IJ) * DZV(IJ, K) *
     &                   AMSKV(IJ, K)
            PX (IJ, K) = (PX(IJ, K-1) + PXM(IJ, K-1) + PXM(IJ, K)) *
     &                   AMSKV(IJ, K)
            GX (IJ, K) = GX(IJ, K)
     &                 + (COR(IJ) * VX(IJ, K) + PX(IJ, K)) *
     &                   AMSKV(IJ, K)
            XX (IJ, K) = XX(IJ, K)
     &                 + PX(IJ, K) * AMSKV(IJ, K)
            GXX(IJ)    = GXX(IJ) + XX(IJ, K) * DZV(IJ, K) * AMSKV(IJ, K)
            PYM(IJ, K) = (  R(IJ+LNE, K) + R(IJ+LN, K)
     &                    - R(IJ+LE , K) - R(IJ   , K)) *
     &                   COF * RYM(IJ) * RYU(IJ) * DZV(IJ, K) *
     &                   AMSKV(IJ, K)
            PY (IJ, K) = (PY(IJ, K-1) + PYM(IJ, K-1) + PYM(IJ, K)) *
     &                   AMSKV(IJ, K)
            GY (IJ, K) = GY(IJ, K)
     &                 - (COR(IJ) * UX(IJ, K) - PY(IJ, K)) *
     &                   AMSKV(IJ, K)
            YY (IJ, K) = YY(IJ, K)
     &                 + PY(IJ, K) * AMSKV(IJ, K)
            GYY(IJ)    = GYY(IJ) + YY(IJ, K) * DZV(IJ, K) * AMSKV(IJ, K)
  230    CONTINUE
  240 CONTINUE
C -----------------------------------------------------------------------
      DO K = KSTR, KSTR+KZ-1
         DO IJ = IJVSTR, IJVEND
         GXX(IJ) = GXX(IJ)
     &           + (  (FUX(IJ+LE, K) - FUX(IJ, K)) * RX
     &              + (FUY(IJ+LN, K) - FUY(IJ, K)) * RYM(IJ)) *
     &             RXU(IJ) * RYU(IJ) * DZSIG(IJ,K)
     &           + (  FUNE(IJ+LNE, K) - FUNE(IJ   , K)
     &              + FUSE(IJ+LE , K) - FUSE(IJ+LN, K)) *
     &             RX * RYM(IJ) * RXU(IJ) * RYU(IJ) * DZSIG(IJ,K)
         GYY(IJ) = GYY(IJ)
     &           + (  (FVX(IJ+LE, K) - FVX(IJ, K)) * RX
     &              + (FVY(IJ+LN, K) - FVY(IJ, K)) * RYM(IJ)) *
     &             RXU(IJ) * RYU(IJ) * DZSIG(IJ,K)
     &           + (  FVNE(IJ+LNE, K) - FVNE(IJ   , K)
     &              + FVSE(IJ+LE , K) - FVSE(IJ+LN, K)) *
     &             RX * RYM(IJ) * RXU(IJ) * RYU(IJ) * DZSIG(IJ,K)
         END DO
      END DO
      DO K = KSTR+KZ, KEND
         DO IJ = IJVSTR, IJVEND
         GXX(IJ) = GXX(IJ)
     &           + (  (FUX(IJ+LE, K) - FUX(IJ, K)) * RX
     &              + (FUY(IJ+LN, K) - FUY(IJ, K)) * RYM(IJ)) *
     &             RXU(IJ) * RYU(IJ)
     &           + (  FUNE(IJ+LNE, K) - FUNE(IJ   , K)
     &              + FUSE(IJ+LE , K) - FUSE(IJ+LN, K)) *
     &             RX * RYM(IJ) * RXU(IJ) * RYU(IJ)
         GYY(IJ) = GYY(IJ)
     &           + (  (FVX(IJ+LE, K) - FVX(IJ, K)) * RX
     &              + (FVY(IJ+LN, K) - FVY(IJ, K)) * RYM(IJ)) *
     &             RXU(IJ) * RYU(IJ)
     &           + (  FVNE(IJ+LNE, K) - FVNE(IJ   , K)
     &              + FVSE(IJ+LE , K) - FVSE(IJ+LN, K)) *
     &             RX * RYM(IJ) * RXU(IJ) * RYU(IJ)
         END DO
      END DO
C -----------------------------------------------------------------------
C      DO 300 IJ = IJVSTR, IJVEND+NXDIM+1
C         FUX (IJ, KSTR) = FUX (IJ, KSTR) * DZSIG(IJ, KSTR)
C         FUY (IJ, KSTR) = FUY (IJ, KSTR) * DZSIG(IJ, KSTR)
C         FUNE(IJ, KSTR) = FUNE(IJ, KSTR) * DZSIG(IJ, KSTR)
C         FUSE(IJ, KSTR) = FUSE(IJ, KSTR) * DZSIG(IJ, KSTR)
C         FVX (IJ, KSTR) = FVX (IJ, KSTR) * DZSIG(IJ, KSTR)
C         FVY (IJ, KSTR) = FVY (IJ, KSTR) * DZSIG(IJ, KSTR)
C         FVNE(IJ, KSTR) = FVNE(IJ, KSTR) * DZSIG(IJ, KSTR)
C         FVSE(IJ, KSTR) = FVSE(IJ, KSTR) * DZSIG(IJ, KSTR)
C 300  CONTINUE
C      DO 330 K = KSTR+1, KSTR+KZ-1
C         DO 320 IJ = IJVSTR, IJVEND+NXDIM+1
C            FUX (IJ, KSTR) = FUX (IJ, KSTR)
C     &                     + FUX (IJ, K) * DZSIG(IJ, K)
C            FUY (IJ, KSTR) = FUY (IJ, KSTR)
C     &                     + FUY (IJ, K) * DZSIG(IJ, K)
C            FUNE(IJ, KSTR) = FUNE(IJ, KSTR)
C     &                     + FUNE(IJ, K) * DZSIG(IJ, K)
C            FUSE(IJ, KSTR) = FUSE(IJ, KSTR)
C     &                     + FUSE(IJ, K) * DZSIG(IJ, K)
C            FVX (IJ, KSTR) = FVX (IJ, KSTR)
C     &                     + FVX (IJ, K) * DZSIG(IJ, K)
C            FVY (IJ, KSTR) = FVY (IJ, KSTR)
C     &                     + FVY (IJ, K) * DZSIG(IJ, K)
C            FVNE(IJ, KSTR) = FVNE(IJ, KSTR)
C     &                     + FVNE(IJ, K) * DZSIG(IJ, K)
C            FVSE(IJ, KSTR) = FVSE(IJ, KSTR)
C     &                     + FVSE(IJ, K) * DZSIG(IJ, K)
C 320     CONTINUE
C 330  CONTINUE
C      DO 410 K = KSTR+KZ, KEND
C         DO 400 IJ = IJVSTR, IJVEND+NXDIM+1
C            FUX (IJ, KSTR) = FUX (IJ, KSTR) + FUX (IJ, K)
C            FUY (IJ, KSTR) = FUY (IJ, KSTR) + FUY (IJ, K)
C            FUNE(IJ, KSTR) = FUNE(IJ, KSTR) + FUNE(IJ, K)
C            FUSE(IJ, KSTR) = FUSE(IJ, KSTR) + FUSE(IJ, K)
C            FVX (IJ, KSTR) = FVX (IJ, KSTR) + FVX (IJ, K)
C            FVY (IJ, KSTR) = FVY (IJ, KSTR) + FVY (IJ, K)
C            FVNE(IJ, KSTR) = FVNE(IJ, KSTR) + FVNE(IJ, K)
C            FVSE(IJ, KSTR) = FVSE(IJ, KSTR) + FVSE(IJ, K)
C  400    CONTINUE
C  410 CONTINUE

      DO IJ = 1, NXYDIM
         UU(IJ) = 0.D0
         VV(IJ) = 0.D0
         UV(IJ) = 0.D0
      END DO

      DO 500 IJ = IJVSTR, IJVEND
         UU(IJ) = UY(IJ, KSTR) * UY(IJ, KSTR) *
     &            DZSIG(IJ, KSTR) * AMSKV(IJ, KSTR)
         VV(IJ) = VY(IJ, KSTR) * VY(IJ, KSTR) *
     &            DZSIG(IJ, KSTR) * AMSKV(IJ, KSTR)
         UV(IJ) = UY(IJ, KSTR) * VY(IJ, KSTR) *
     &            DZSIG(IJ, KSTR) * AMSKV(IJ, KSTR)
  500 CONTINUE
      DO 530 K = KSTR+1, KSTR+KZ-1
         DO 520 IJ = IJVSTR, IJVEND
            UU(IJ) = UU(IJ)
     &             + UY(IJ, K) * UY(IJ, K) *
     &               DZSIG(IJ, K) * AMSKV(IJ, K)
            VV(IJ) = VV(IJ)
     &             + VY(IJ, K) * VY(IJ, K) *
     &               DZSIG(IJ, K) * AMSKV(IJ, K)
            UV(IJ) = UV(IJ)
     &             + UY(IJ, K) * VY(IJ, K) *
     &               DZSIG(IJ, K) * AMSKV(IJ, K)
  520    CONTINUE
  530 CONTINUE
      DO 610 K = KSTR+KZ, KEND
         DO 600 IJ = IJVSTR, IJVEND
            UU(IJ) = UU(IJ)
     &             + UY(IJ, K) * UY(IJ, K) * DZV(IJ, K) * AMSKV(IJ, K)
            VV(IJ) = VV(IJ)
     &             + VY(IJ, K) * VY(IJ, K) * DZV(IJ, K) * AMSKV(IJ, K)
            UV(IJ) = UV(IJ)
     &             + UY(IJ, K) * VY(IJ, K) * DZV(IJ, K) * AMSKV(IJ, K)
  600    CONTINUE
  610 CONTINUE

      DO 700 IJ = IJVSTR, IJVEND
         GXX(IJ) = GXX(IJ) + VV(IJ) * HYXU(IJ) - UV(IJ) * HXYU(IJ)
         GYY(IJ) = GYY(IJ) + UU(IJ) * HXYU(IJ) - UV(IJ) * HYXU(IJ)
  700 CONTINUE

C      DO 700 IJ = IJVSTR, IJVEND
C         GXX(IJ) = GXX(IJ)
C     &           + (  (FUX(IJ+LE, KSTR) - FUX(IJ, KSTR)) * RX
C     &              + (FUY(IJ+LN, KSTR) - FUY(IJ, KSTR)) * RYM(IJ)) *
C     &             RXU(IJ) * RYU(IJ)
C     &           + (  FUNE(IJ+LNE, KSTR) - FUNE(IJ   , KSTR)
C     &              + FUSE(IJ+LE , KSTR) - FUSE(IJ+LN, KSTR)) *
C     &             RX * RYM(IJ) * RXU(IJ) * RYU(IJ)
C     &           + VV(IJ) * HYXU(IJ) - UV(IJ) * HXYU(IJ)
C         GYY(IJ) = GYY(IJ)
C     &           + (  (FVX(IJ+LE, KSTR) - FVX(IJ, KSTR)) * RX
C     &              + (FVY(IJ+LN, KSTR) - FVY(IJ, KSTR)) * RYM(IJ)) *
C     &             RXU(IJ) * RYU(IJ)
C     &           + (  FVNE(IJ+LNE, KSTR) - FVNE(IJ   , KSTR)
C     &              + FVSE(IJ+LE , KSTR) - FVSE(IJ+LN, KSTR)) *
C     &             RX * RYM(IJ) * RXU(IJ) * RYU(IJ)
C     &           + UU(IJ) * HXYU(IJ) - UV(IJ) * HYXU(IJ)
C  700 CONTINUE

      RETURN
      END
#ifdef OPT_BBL
c *********************************************************************

      SUBROUTINE TNDUVB(
     M                GXX,    GYY,
     M                 GX,     GY,
     I                 XX,     YY,
     I                 UY,     VY,     TY,
     I                 UX,     VX,
     I                FUX,    FUY,   FUNE,   FUSE,
     I                FVX,    FVY,   FVNE,   FVSE)
      use zocdim
      use zocgrd
      use zocmsk
      use zocphy
      IMPLICIT NONE

      REAL*8     GXX(NXYDIM)       ,    GYY(NXYDIM)
      REAL*8      GX(NXYDIM, NZDIM),     GY(NXYDIM, NZDIM)
      REAL*8      XX(NXYDIM, NZDIM),     YY(NXYDIM, NZDIM)
      REAL*8      UX(NXYDIM, NZDIM),     VX(NXYDIM, NZDIM)
      REAL*8      UY(NXYDIM, NZDIM),     VY(NXYDIM, NZDIM)
      REAL*8      TY(NXYDIM, NZDIM, NTDIM)
      REAL*8     FUX(NXYDIM, NZDIM),    FVX(NXYDIM, NZDIM)
      REAL*8     FUY(NXYDIM, NZDIM),    FVY(NXYDIM, NZDIM)
      REAL*8    FUNE(NXYDIM, NZDIM),   FVNE(NXYDIM, NZDIM)
      REAL*8    FUSE(NXYDIM, NZDIM),   FVSE(NXYDIM, NZDIM)

      REAL*8      PX(NXYDIM, NZDIM),    PXM(NXYDIM, NZDIM)
      REAL*8      PY(NXYDIM, NZDIM),    PYM(NXYDIM, NZDIM)
      COMMON /WORK2/ PX, PXM,  PY, PYM

      REAL*8  C0B(NXYDIM), C1B(NXYDIM), C2B(NXYDIM), C3B(NXYDIM)
      REAL*8  C4B(NXYDIM), C5B(NXYDIM), C6B(NXYDIM)
      REAL*8  D0B(NXYDIM), D1B(NXYDIM), D2B(NXYDIM), D3B(NXYDIM)
      REAL*8  D4B(NXYDIM), D5B(NXYDIM), D6B(NXYDIM), D7B(NXYDIM)
      REAL*8  D8B(NXYDIM), D9B(NXYDIM)
      SAVE C0B, C1B, C2B, C3B, C4B, C5B, C6B,
     &     D0B, D1B, D2B, D3B, D4B, D5B, D6B, D7B, D8B, D9B

      REAL*8  DEPT0(NXYDIM)
      REAL*8    COF
      LOGICAL OFIRST
      SAVE DEPT0, COF, OFIRST

      REAL*8    RBBL,     TB,     SB
      REAL*8  C0, C1, C2, C3, C4, C5, C6
      REAL*8  D0, D1, D2, D3, D4, D5, D6, D7, D8, D9
      REAL*8      TL,     SL,   RTMP
      REAL*8       P,     PE,     PN,    PNE
      INTEGER     IJ,      K,     KK
      INTEGER   IJLE,   IJLN,   IJNE

      DATA OFIRST / .TRUE. /

      RBBL(TB, SB, C0, C1, C2, C3, C4, C5, C6,
     &     D0, D1, D2, D3, D4, D5, D6, D7, D8, D9) = 
     &    (C0 + (C1 + (C2 + C3 * TB) * TB) * TB
     &        + (C4 + C5 * TB + C6 * SB) * SB)
     &  / (D0 + (D1 + (D2 + (D3 + D4 * TB) * TB) * TB) * TB
     &        + (D5 + (D6 + D7 * TB * TB) * TB
     &              + (D8 + D9 * TB * TB) * SQRT(SB)) * SB)
     &  - 1.D3

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF

      IF (OFIRST) THEN
         OFIRST = .FALSE.
         CALL SECFBB(
     O               C0B,    C1B,    C2B,    C3B,
     O               C4B,    C5B,    C6B,
     O               D0B,    D1B,    D2B,    D3B,    D4B,
     O               D5B,    D6B,    D7B,    D8B,    D9B)

         COF = -0.5D0 * gravit / rhoo * 1.D-3
         DO IJ = 1, NXYDIM
            DEPT0(IJ) = 0.D0
            DO K = KSTR, NBOT(IJ)-1
               DEPT0(IJ) = DEPT0(IJ) + DZ0(K)
            END DO
         END DO
      END IF

      DO IJ = IJVSTR, IJVEND

         K = MAX(NBOTV(IJ)-1, 1)
         IJLE = IJ + LE
         IJLN = IJ + LN
         IJNE = IJ + LNE

         TL = TY(IJ, K, 1)
         SL = TY(IJ, K, 2)
         RTMP = RBBL(TL, SL,
     &               C0B(IJ), C1B(IJ), C2B(IJ), C3B(IJ),
     &               C4B(IJ), C5B(IJ), C6B(IJ),
     &               D0B(IJ), D1B(IJ), D2B(IJ), D3B(IJ), D4B(IJ),
     &               D5B(IJ), D6B(IJ), D7B(IJ), D8B(IJ), D9B(IJ))
         P = RTMP * DZ0(K) * 0.5D0
         DO KK = K+1, NBOT(IJ)-1
            TL = TY(IJ, KK, 1)
            SL = TY(IJ, KK, 2)
            RTMP = RBBL(TL, SL,
     &               C0B(IJ), C1B(IJ), C2B(IJ), C3B(IJ),
     &               C4B(IJ), C5B(IJ), C6B(IJ),
     &               D0B(IJ), D1B(IJ), D2B(IJ), D3B(IJ), D4B(IJ),
     &               D5B(IJ), D6B(IJ), D7B(IJ), D8B(IJ), D9B(IJ))
            P = P + RTMP * DZ0(KK)
         END DO
         TL = TY(IJ, KEND, 1)
         SL = TY(IJ, KEND, 2)
         RTMP = RBBL(TL, SL,
     &               C0B(IJ), C1B(IJ), C2B(IJ), C3B(IJ),
     &               C4B(IJ), C5B(IJ), C6B(IJ),
     &               D0B(IJ), D1B(IJ), D2B(IJ), D3B(IJ), D4B(IJ),
     &               D5B(IJ), D6B(IJ), D7B(IJ), D8B(IJ), D9B(IJ))
         P = P + RTMP * DZ0(KEND) * 0.5D0

         TL = TY(IJLE, K, 1)
         SL = TY(IJLE, K, 2)
         RTMP = RBBL(TL, SL,
     &               C0B(IJ), C1B(IJ), C2B(IJ), C3B(IJ),
     &               C4B(IJ), C5B(IJ), C6B(IJ),
     &               D0B(IJ), D1B(IJ), D2B(IJ), D3B(IJ), D4B(IJ),
     &               D5B(IJ), D6B(IJ), D7B(IJ), D8B(IJ), D9B(IJ))
         PE = RTMP * DZ0(K) * 0.5D0
         DO KK = K+1, NBOT(IJLE)-1
            TL = TY(IJLE, KK, 1)
            SL = TY(IJLE, KK, 2)
            RTMP = RBBL(TL, SL,
     &               C0B(IJ), C1B(IJ), C2B(IJ), C3B(IJ),
     &               C4B(IJ), C5B(IJ), C6B(IJ),
     &               D0B(IJ), D1B(IJ), D2B(IJ), D3B(IJ), D4B(IJ),
     &               D5B(IJ), D6B(IJ), D7B(IJ), D8B(IJ), D9B(IJ))
            PE = PE + RTMP * DZ0(KK)
         END DO
         TL = TY(IJLE, KEND, 1)
         SL = TY(IJLE, KEND, 2)
         RTMP = RBBL(TL, SL,
     &               C0B(IJ), C1B(IJ), C2B(IJ), C3B(IJ),
     &               C4B(IJ), C5B(IJ), C6B(IJ),
     &               D0B(IJ), D1B(IJ), D2B(IJ), D3B(IJ), D4B(IJ),
     &               D5B(IJ), D6B(IJ), D7B(IJ), D8B(IJ), D9B(IJ))
         PE = PE + RTMP * DZ0(KEND) * 0.5D0

         TL = TY(IJLN, K, 1)
         SL = TY(IJLN, K, 2)
         RTMP = RBBL(TL, SL,
     &               C0B(IJ), C1B(IJ), C2B(IJ), C3B(IJ),
     &               C4B(IJ), C5B(IJ), C6B(IJ),
     &               D0B(IJ), D1B(IJ), D2B(IJ), D3B(IJ), D4B(IJ),
     &               D5B(IJ), D6B(IJ), D7B(IJ), D8B(IJ), D9B(IJ))
         PN = RTMP * DZ0(K) * 0.5D0
         DO KK = K+1, NBOT(IJLN)-1
            TL = TY(IJLN, KK, 1)
            SL = TY(IJLN, KK, 2)
            RTMP = RBBL(TL, SL,
     &               C0B(IJ), C1B(IJ), C2B(IJ), C3B(IJ),
     &               C4B(IJ), C5B(IJ), C6B(IJ),
     &               D0B(IJ), D1B(IJ), D2B(IJ), D3B(IJ), D4B(IJ),
     &               D5B(IJ), D6B(IJ), D7B(IJ), D8B(IJ), D9B(IJ))
            PN = PN + RTMP * DZ0(KK)
         END DO
         TL = TY(IJLN, KEND, 1)
         SL = TY(IJLN, KEND, 2)
         RTMP = RBBL(TL, SL,
     &               C0B(IJ), C1B(IJ), C2B(IJ), C3B(IJ),
     &               C4B(IJ), C5B(IJ), C6B(IJ),
     &               D0B(IJ), D1B(IJ), D2B(IJ), D3B(IJ), D4B(IJ),
     &               D5B(IJ), D6B(IJ), D7B(IJ), D8B(IJ), D9B(IJ))
         PN = PN + RTMP * DZ0(KEND) * 0.5D0

         TL = TY(IJNE, K, 1)
         SL = TY(IJNE, K, 2)
         RTMP = RBBL(TL, SL,
     &               C0B(IJ), C1B(IJ), C2B(IJ), C3B(IJ),
     &               C4B(IJ), C5B(IJ), C6B(IJ),
     &               D0B(IJ), D1B(IJ), D2B(IJ), D3B(IJ), D4B(IJ),
     &               D5B(IJ), D6B(IJ), D7B(IJ), D8B(IJ), D9B(IJ))
         PNE = RTMP * DZ0(K) * 0.5D0
         DO KK = K+1, NBOT(IJNE)-1
            TL = TY(IJNE, KK, 1)
            SL = TY(IJNE, KK, 2)
            RTMP = RBBL(TL, SL,
     &               C0B(IJ), C1B(IJ), C2B(IJ), C3B(IJ),
     &               C4B(IJ), C5B(IJ), C6B(IJ),
     &               D0B(IJ), D1B(IJ), D2B(IJ), D3B(IJ), D4B(IJ),
     &               D5B(IJ), D6B(IJ), D7B(IJ), D8B(IJ), D9B(IJ))
            PNE = PNE + RTMP * DZ0(KK)
         END DO
         TL = TY(IJNE, KEND, 1)
         SL = TY(IJNE, KEND, 2)
         RTMP = RBBL(TL, SL,
     &               C0B(IJ), C1B(IJ), C2B(IJ), C3B(IJ),
     &               C4B(IJ), C5B(IJ), C6B(IJ),
     &               D0B(IJ), D1B(IJ), D2B(IJ), D3B(IJ), D4B(IJ),
     &               D5B(IJ), D6B(IJ), D7B(IJ), D8B(IJ), D9B(IJ))
         PNE = PNE + RTMP * DZ0(KEND) * 0.5D0

         TL = 0.25D0 * (  TY(IJ  , KEND, 1) + TY(IJLE, KEND, 1)
     &                  + TY(IJLN, KEND, 1) + TY(IJNE, KEND, 1))
         SL = 0.25D0 * (  TY(IJ  , KEND, 2) + TY(IJLE, KEND, 2)
     &                  + TY(IJLN, KEND, 2) + TY(IJNE, KEND, 2))
         RTMP = RBBL(TL, SL,
     &               C0B(IJ), C1B(IJ), C2B(IJ), C3B(IJ),
     &               C4B(IJ), C5B(IJ), C6B(IJ),
     &               D0B(IJ), D1B(IJ), D2B(IJ), D3B(IJ), D4B(IJ),
     &               D5B(IJ), D6B(IJ), D7B(IJ), D8B(IJ), D9B(IJ))

         PX(IJ, KEND) = PX(IJ, K)
     &                + (  PNE + PE - PN - P
     &                   - RTMP * (  DEPT0(IJNE) + DEPT0(IJLE)
     &                             - DEPT0(IJLN) - DEPT0(IJ  ))) *
     &                  COF * RX * RXU(IJ) * AMSKVB(IJ)
         GX(IJ, KEND) = GX(IJ, KEND)
     &                + (COR(IJ) * VX(IJ, KEND) + PX(IJ, KEND)) *
     &                  AMSKVB(IJ)
         XX(IJ, KEND) = XX(IJ, KEND) + PX(IJ, KEND) * AMSKVB(IJ)
         GXX(IJ)      = GXX(IJ)
     &                + XX(IJ, KEND) * DZV(IJ, KEND) * AMSKVB(IJ)

         PY(IJ, KEND) = PY(IJ, K)
     &                + (  PNE + PN - PE - P
     &                   - RTMP * (  DEPT0(IJNE) + DEPT0(IJLN)
     &                             - DEPT0(IJLE) - DEPT0(IJ  ))) *
     &                  COF * RYM(IJ) * RYU(IJ) * AMSKVB(IJ)
         GY(IJ, KEND) = GY(IJ, KEND)
     &                - (COR(IJ) * UX(IJ, KEND) - PY(IJ, KEND)) *
     &                  AMSKVB(IJ)
         YY(IJ, KEND) = YY(IJ, KEND) + PY(IJ, KEND) * AMSKVB(IJ)
         GYY(IJ)      = GYY(IJ)
     &                + YY(IJ, KEND) * DZV(IJ, KEND) * AMSKVB(IJ)

      END DO

      DO IJ = IJVSTR, IJVEND
         GXX(IJ) = GXX(IJ)
     &           + (  VY(IJ, KEND) * VY(IJ, KEND) * HYXU(IJ)
     &              - UY(IJ, KEND) * VY(IJ, KEND) * HXYU(IJ)) *
     &             DZV(IJ, KEND) * AMSKVB(IJ)
         GYY(IJ) = GYY(IJ)
     &           + (  UY(IJ, KEND) * UY(IJ, KEND) * HXYU(IJ)
     &              - UY(IJ, KEND) * VY(IJ, KEND) * HYXU(IJ)) *
     &             DZV(IJ, KEND) * AMSKVB(IJ)
      END DO

      RETURN
      END
#endif
