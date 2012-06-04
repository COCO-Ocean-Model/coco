      SUBROUTINE ADVVEL(
     M                 GX,     GY,     XX,     YY,
     I                 UY,     UX,     VY,     VX,
     I                 HX,     HY,
     I               UADV,   VADV,   WADV)
      use zocdim
      use zocgrd
      use zocmsk
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Advection and metric terms of the equation of motion.
c
c  HISTORY
c     '03.04.22  H.Hasumi: from COCO3.4
c     '07.04.23  H.Hasumi
c     '08.06.11  H.Hasumi: initial/final processing
c     '08.07.10  H.Hasumi: initial/final processing
c
c ---------------------------------------------------------------------

      REAL*8      UY(NXYDIM, NZDIM),     UX(NXYDIM, NZDIM)
      REAL*8      VY(NXYDIM, NZDIM),     VX(NXYDIM, NZDIM)
      REAL*8      HY(NXYDIM)       ,     HX(NXYDIM)
      REAL*8    UADV(NXYDIM, NZDIM),   VADV(NXYDIM, NZDIM)
      REAL*8    WADV(NXYDIM, NZDIM)
      REAL*8      GX(NXYDIM, NZDIM),     GY(NXYDIM, NZDIM)
      REAL*8      XX(NXYDIM, NZDIM),     YY(NXYDIM, NZDIM)

      REAL*8     FUX(NXYDIM, NZDIM),    FVX(NXYDIM, NZDIM)
      REAL*8     FUY(NXYDIM, NZDIM),    FVY(NXYDIM, NZDIM)
      REAL*8     FUZ(NXYDIM, NZDIM),    FVZ(NXYDIM, NZDIM)
      REAL*8      RZ(NXYDIM, NZDIM),    RZM(NXYDIM, NZDIM)
c      COMMON /WORK/ FUX, FUY, FUZ, FVX, FVY, FVZ, RZ, RZM
*** RZ is taken over from the previous subroutine VSCVEL ***

      REAL*8     DIV(NXYDIM, KSTR:KSTR+KZ-1)
      REAL*8   HVBOT(NXYDIM), HVBOTX(NXYDIM)
      INTEGER     IJ,      K

      IF (OINIT .OR. OFINAL) THEN
         RETURN
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

      DO IJ = IJVSTR, IJVEND
         HVBOT(IJ) = (  (HY(IJ)    + HY(IJ+LE) ) * DY(IJ)
     &                + (HY(IJ+LN) + HY(IJ+LNE)) * DY(IJ+LN)) *
     &               RYM(IJ) * 0.25D0
     &             + ZBOT
         HVBOT(IJ) = 1.D0 / HVBOT(IJ)
         HVBOTX(IJ) = (  (HX(IJ)    + HX(IJ+LE) ) * DY(IJ)
     &                + (HX(IJ+LN) + HX(IJ+LNE)) * DY(IJ+LN)) *
     &               RYM(IJ) * 0.25D0
     &             + ZBOT
         HVBOTX(IJ) = 1.D0 / HVBOTX(IJ)
      END DO

      DO K = KSTR, KSTR+KZ-1
         DO IJ = IJVSTR, IJVEND
            RZ (IJ, K) = 1.D0 * RS (K) * HVBOTX(IJ)
         END DO
      END DO
      DO K = KSTR+KZ, KEND
         DO IJ = IJVSTR, IJVEND
            RZ (IJ, K) = 1.D0 / DZV(IJ, K)
         END DO
      END DO

      DO K = KSTR, KSTR+KZ-2
         DO IJ = IJVSTR, IJVEND
            DIV(IJ, K) = (  (  (  UADV(IJ+LE, K) *
     &                            (HYU(IJ+LE) + HYU(IJ))
     &                          - UADV(IJ, K) *
     &                            (HYU(IJ) + HYU(IJ+LW))) *
     &                         0.5D0 * RX
     &                       + (  VADV(IJ+LN, K) *
     &                            (HXU(IJ+LN) + HXU(IJ))
     &                          - VADV(IJ, K) *
     &                            (HXU(IJ) + HXU(IJ+LS))) *
     &                         0.5D0 * RYM(IJ)) * RXU(IJ) * RYU(IJ)
     &                    + (WADV(IJ, K) - WADV(IJ, K+1)) * RS(K)) *
     &                   AMSKV(IJ, K)
         END DO
      END DO
      K = KSTR+KZ-1
         DO IJ = IJVSTR, IJVEND
            DIV(IJ, K) = (  (  (  UADV(IJ+LE, K) *
     &                            (HYU(IJ+LE) + HYU(IJ))
     &                          - UADV(IJ, K) *
     &                            (HYU(IJ) + HYU(IJ+LW))) *
     &                         0.5D0 * RX
     &                       + (  VADV(IJ+LN, K) *
     &                            (HXU(IJ+LN) + HXU(IJ))
     &                          - VADV(IJ, K) *
     &                            (HXU(IJ) + HXU(IJ+LS))) *
     &                         0.5D0 * RYM(IJ)) * RXU(IJ) * RYU(IJ)
     &                    + (  WADV(IJ, K)
     &                       - WADV(IJ, K+1) * HVBOT(IJ)) * RS(K)) *
     &                   AMSKV(IJ, K)
         END DO

      DO IJ = IJVSTR, IJVEND
         FUZ(IJ, KSTR) = 0.D0
         FVZ(IJ, KSTR) = 0.D0
      END DO
      DO 1210 K = KSTR+1, KEND
         DO 1200 IJ = IJVSTR, IJVEND
            FUZ(IJ, K) =
     &                 - WADV(IJ, K) * 0.5D0 *
     &                   (UY(IJ, K-1) + UY(IJ, K))
            FVZ(IJ, K) =
     &                 - WADV(IJ, K) * 0.5D0 *
     &                   (VY(IJ, K-1) + VY(IJ, K))
 1200    CONTINUE
 1210 CONTINUE

      DO 1310 K = KSTR, KEND
         DO 1300 IJ = IJVSTR, IJVEND+NXDIM
            FUY(IJ, K) =
     &                 - VADV(IJ, K) * 0.5D0 *
     &                   (  UY(IJ, K) * HXU(IJ)
     &                    + UY(IJ+LS, K) * HXU(IJ+LS))
            FVY(IJ, K) =
     &                 - VADV(IJ, K) * 0.5D0 *
     &                   (  VY(IJ, K) * HXU(IJ)
     &                    + VY(IJ+LS, K) * HXU(IJ+LS))
 1300    CONTINUE
 1310 CONTINUE

      DO 1410 K = KSTR, KEND
         DO 1400 IJ = IJVSTR, IJVEND+1
            FUX(IJ, K) =
     &                 - UADV(IJ, K) * 0.5D0 *
     &                   (  UY(IJ, K) * HYU(IJ)
     &                    + UY(IJ+LW, K) * HYU(IJ+LW))
            FVX(IJ, K) =
     &                 - UADV(IJ, K) * 0.5D0 *
     &                   (  VY(IJ, K) * HYU(IJ)
     &                    + VY(IJ+LW, K) * HYU(IJ+LW))
 1400    CONTINUE
 1410 CONTINUE

      DO 1510 K = KSTR, KSTR+KZ-2
         DO 1500 IJ = IJVSTR, IJVEND
            GX(IJ, K) = (  GX(IJ, K)
     &                   + (  (FUX(IJ+LE, K) - FUX(IJ, K)) * RX
     &                      + (FUY(IJ+LN, K) - FUY(IJ, K)) * RYM(IJ)) *
     &                     RXU(IJ) * RYU(IJ)
     &                   + (FUZ(IJ, K) - FUZ(IJ, K+1)) * RS(K)
     &                   + UY(IJ, K) * DIV(IJ, K)
     &                   + VY(IJ, K) * VY(IJ, K) * HYXU(IJ)
     &                   - UY(IJ, K) * VY(IJ, K) * HXYU(IJ)) *
     &                  AMSKV(IJ, K)
            GY(IJ, K) = (  GY(IJ, K)
     &                   + (  (FVX(IJ+LE, K) - FVX(IJ, K)) * RX
     &                      + (FVY(IJ+LN, K) - FVY(IJ, K)) * RYM(IJ)) *
     &                     RXU(IJ) * RYU(IJ)
     &                   + (FVZ(IJ, K) - FVZ(IJ, K+1)) * RS(K)
     &                   + VY(IJ, K) * DIV(IJ, K)
     &                   + UY(IJ, K) * UY(IJ, K) * HXYU(IJ)
     &                   - UY(IJ, K) * VY(IJ, K) * HYXU(IJ)) *
     &                  AMSKV(IJ, K)
 1500    CONTINUE
 1510 CONTINUE
      K = KSTR+KZ-1
      DO 1520 IJ = IJVSTR, IJVEND
            GX(IJ, K) = (  GX(IJ, K)
     &                   + (  (FUX(IJ+LE, K) - FUX(IJ, K)) * RX
     &                      + (FUY(IJ+LN, K) - FUY(IJ, K)) * RYM(IJ)) *
     &                     RXU(IJ) * RYU(IJ)
     &                   + (  FUZ(IJ, K)
     &                      - FUZ(IJ, K+1) * HVBOT(IJ)) * RS(K)
     &                   + UY(IJ, K) * DIV(IJ, K)
     &                   + VY(IJ, K) * VY(IJ, K) * HYXU(IJ)
     &                   - UY(IJ, K) * VY(IJ, K) * HXYU(IJ)) *
     &                  AMSKV(IJ, K)
            GY(IJ, K) = (  GY(IJ, K)
     &                   + (  (FVX(IJ+LE, K) - FVX(IJ, K)) * RX
     &                      + (FVY(IJ+LN, K) - FVY(IJ, K)) * RYM(IJ)) *
     &                     RXU(IJ) * RYU(IJ)
     &                   + (  FVZ(IJ, K)
     &                      - FVZ(IJ, K+1) * HVBOT(IJ)) * RS(K)
     &                   + VY(IJ, K) * DIV(IJ, K)
     &                   + UY(IJ, K) * UY(IJ, K) * HXYU(IJ)
     &                   - UY(IJ, K) * VY(IJ, K) * HYXU(IJ)) *
     &                  AMSKV(IJ, K)
 1520 CONTINUE
      DO 1540 K = KSTR+KZ, KEND
         DO 1530 IJ = IJVSTR, IJVEND
            GX(IJ, K) = (  GX(IJ, K)
     &                   + (  (  (  FUX(IJ+LE, K)
     &                            - FUX(IJ   , K)) * RX
     &                         + (  FUY(IJ+LN, K)
     &                            - FUY(IJ   , K)) * RYM(IJ)) *
     &                        RXU(IJ) * RYU(IJ)
     &                      + FUZ(IJ, K) - FUZ(IJ, K+1)) * RZ(IJ, K)
     &                   + VY(IJ, K) * VY(IJ, K) * HYXU(IJ)
     &                   - UY(IJ, K) * VY(IJ, K) * HXYU(IJ)) *
     &                  AMSKV(IJ, K)
            GY(IJ, K) = (  GY(IJ, K)
     &                   + (  (  (  FVX(IJ+LE, K)
     &                            - FVX(IJ   , K)) * RX
     &                         + (  FVY(IJ+LN, K)
     &                            - FVY(IJ   , K)) * RYM(IJ)) *
     &                        RXU(IJ) * RYU(IJ)
     &                      + FVZ(IJ, K) - FVZ(IJ, K+1)) * RZ(IJ, K)
     &                   + UY(IJ, K) * UY(IJ, K) * HXYU(IJ)
     &                   - UY(IJ, K) * VY(IJ, K) * HYXU(IJ)) *
     &                  AMSKV(IJ, K)
 1530    CONTINUE
 1540 CONTINUE

      RETURN
      END
#ifdef OPT_BBL
c *********************************************************************

      SUBROUTINE ADVVLB(
     M                 GX,     GY,     XX,     YY,
     I                 UY,     UX,     VY,     VX,
     I               UADV,   VADV,   WADV)
      use zocdim
      use zocgrd
      use zocmsk
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Advection term of the BBL momentum eqs.
c
c ---------------------------------------------------------------------

      REAL*8      UY(NXYDIM, NZDIM),     UX(NXYDIM, NZDIM)
      REAL*8      VY(NXYDIM, NZDIM),     VX(NXYDIM, NZDIM)
      REAL*8    UADV(NXYDIM, NZDIM),   VADV(NXYDIM, NZDIM)
      REAL*8    WADV(NXYDIM, NZDIM)
      REAL*8      GX(NXYDIM, NZDIM),     GY(NXYDIM, NZDIM)
      REAL*8      XX(NXYDIM, NZDIM),     YY(NXYDIM, NZDIM)

      REAL*8      RZ(NXYDIM)
      SAVE RZ

      REAL*8     FUX(NXYDIM),    FVX(NXYDIM)
      REAL*8     FUY(NXYDIM),    FVY(NXYDIM)
      REAL*8     FUZ(NXYDIM),    FVZ(NXYDIM)
c      COMMON /WORK/ FUX, FUY, FUZ, FVX, FVY, FVZ

      INTEGER     IJ,    KUP

      LOGICAL OFIRST
      SAVE OFIRST

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF

      IF (OFIRST) THEN
         OFIRST = .FALSE.
         DO IJ = 1, NXYDIM
            RZ   (IJ) = 1.D0 / DZV(IJ, KEND)
         END DO
      END IF

      DO IJ = 1, NXYDIM
         FUX(IJ) = 0.D0
         FVX(IJ) = 0.D0
         FUY(IJ) = 0.D0
         FVY(IJ) = 0.D0
      END DO

      DO IJ = IJVSTR, IJVEND
         KUP = MAX(NBOTV(IJ)-1, 1)
         FUZ(IJ) = - WADV(IJ, KEND) * 0.5D0 *
     &               (UY(IJ, KUP) + UY(IJ, KEND))
         FVZ(IJ) = - WADV(IJ, KEND) * 0.5D0 *
     &               (VY(IJ, KUP) + VY(IJ, KEND))
      END DO

      DO IJ = IJVSTR, IJVEND+NXDIM
         FUY(IJ) = - VADV(IJ, KEND) * 0.5D0 *
     &               (  UY(IJ, KEND) * HXU(IJ)
     &                + UY(IJ+LS, KEND) * HXU(IJ+LS))
         FVY(IJ) = - VADV(IJ, KEND) * 0.5D0 *
     &               (  VY(IJ, KEND) * HXU(IJ)
     &                + VY(IJ+LS, KEND) * HXU(IJ+LS))
      END DO

      DO IJ = IJVSTR, IJVEND+1
         FUX(IJ) = - UADV(IJ, KEND) * 0.5D0 *
     &               (  UY(IJ, KEND) * HYU(IJ)
     &                + UY(IJ+LW, KEND) * HYU(IJ+LW))
         FVX(IJ) = - UADV(IJ, KEND) * 0.5D0 *
     &               (  VY(IJ, KEND) * HYU(IJ)
     &                + VY(IJ+LW, KEND) * HYU(IJ+LW))
      END DO

      DO IJ = IJVSTR, IJVEND
         GX(IJ, KEND) = (  GX(IJ, KEND)
     &                   + (  (  (FUX(IJ+LE) - FUX(IJ)) * RX
     &                         + (FUY(IJ+LN) - FUY(IJ)) * RYM(IJ)) *
     &                        RXU(IJ) * RYU(IJ)
     &                      + FUZ(IJ)) * RZ(IJ)
     &                   + VY(IJ, KEND) * VY(IJ, KEND) * HYXU(IJ)
     &                   - UY(IJ, KEND) * VY(IJ, KEND) * HXYU(IJ)) *
     &                  AMSKVB(IJ)
         GY(IJ, KEND) = (  GY(IJ, KEND)
     &                   + (  (  (FVX(IJ+LE) - FVX(IJ)) * RX
     &                         + (FVY(IJ+LN) - FVY(IJ)) * RYM(IJ)) *
     &                        RXU(IJ) * RYU(IJ)
     &                      + FVZ(IJ)) * RZ(IJ)
     &                   + UY(IJ, KEND) * UY(IJ, KEND) * HXYU(IJ)
     &                   - UY(IJ, KEND) * VY(IJ, KEND) * HYXU(IJ)) *
     &                  AMSKVB(IJ)
      END DO

      RETURN
      END
#endif
