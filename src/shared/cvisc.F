      SUBROUTINE VSCVEL(
     M                 GX,     GY,     XX,     YY,
     I                 UY,     UX,     VY,     VX,
     I                 HX,     HY,
     I                AMV,   TAUX,   TAUY)
      use zocdim
      use zocgrd
      use zocmsk
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Viscosity term of the equation of motion.
c
c  HISTORY
c     '03.04.21  H.Hasumi: from COCO3.4
c     '07.04.16  H.Hasumi: Formulation changed
c     '07.04.23  H.Hasumi
c     '07.06.11  H.Hasumi: initial/final processing
c     '07.07.10  H.Hasumi: initial/final processing
c
c ---------------------------------------------------------------------

      REAL*8      UY(NXYDIM, NZDIM),     UX(NXYDIM, NZDIM)
      REAL*8      VY(NXYDIM, NZDIM),     VX(NXYDIM, NZDIM)
      REAL*8      HY(NXYDIM)       ,     HX(NXYDIM)
      REAL*8     AMV(NXYDIM, NZDIM)
      REAL*8    TAUX(NXYDIM),   TAUY(NXYDIM)
      REAL*8      GX(NXYDIM, NZDIM),     GY(NXYDIM, NZDIM)
      REAL*8      XX(NXYDIM, NZDIM),     YY(NXYDIM, NZDIM)

      REAL*8     FUX(NXYDIM, NZDIM),    FVX(NXYDIM, NZDIM)
      REAL*8     FUY(NXYDIM, NZDIM),    FVY(NXYDIM, NZDIM)
      REAL*8     FUZ(NXYDIM, NZDIM),    FVZ(NXYDIM, NZDIM)
      REAL*8      RZ(NXYDIM, NZDIM),    RZM(NXYDIM, NZDIM)
      REAL*8     SMX(NXYDIM, NZDIM),    SMY(NXYDIM, NZDIM)
      REAL*8     SXX(NXYDIM),    SYY(NXYDIM)
      REAL*8     SXY(NXYDIM),    SYX(NXYDIM)
      REAL*8     SZX(NXYDIM),    SZY(NXYDIM)
c      COMMON /WORK/ FUX, FUY, FUZ, FVX, FVY, FVZ, RZ, RZM
c      EQUIVALENCE(FUZ(1, 1), SMY(1, 1))
c      EQUIVALENCE(FVZ(1, 1), SMX(1, 1))

      REAL*8   HVBOT(NXYDIM)
      INTEGER     IJ,      K
      INTEGER   IJLN,   IJLS,   IJLE,   IJLW
      INTEGER   IJNW,   IJSE,   IJSW
      INTEGER  IFPAR,  JFPAR

      LOGICAL OFIRST
      SAVE OFIRST

      REAL*8     AMH
      SAVE AMH
      NAMELIST /NMVISH/ AMH

      DATA OFIRST / .TRUE. /
      DATA AMH / 0.D0 /

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF

      IF (OFIRST) THEN
         OFIRST = .FALSE.
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMVISH, END=99)
 99      CONTINUE
         WRITE(JFPAR, NMVISH)
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
         DO IJ = IJSTR, IJEND+NXDIM+1
            IJLS = IJ + LS
            IJLW = IJ + LW
            IJLN = IJ + LN
            IJLE = IJ + LE
            IJNW = IJ + LNW
            IJSE = IJ + LSE
            IJSW = IJ + LSW
            SXX(IJ) = (UX(IJ, K) - UX(IJLW, K)) * AMH *
     &                RX / (HXU(IJ) + HXU(IJLW)) * 2.D0
     &              + (VX(IJ, K) + VX(IJLW, K)) * AMH *
     &                (HXYU(IJ) + HXYU(IJLW)) * 0.25D0
     &              - (  VX(IJLN, K) + VX(IJNW, K)
     &                 - VX(IJLS, K) - VX(IJSW, K)) * AMH
     &                / (HYU(IJ) + HYU(IJLW)) / (DY(IJ) + DY(IJLN))
     &              - (UX(IJ, K) + UX(IJLW, K)) * AMH *
     &                (HYXU(IJ) + HYXU(IJLW)) * 0.25D0
            SXY(IJ) = (UX(IJ, K) - UX(IJLS, K)) * AMH *
     &                RY(IJ) / (HYU(IJ) + HYU(IJLS)) * 2.D0
     &              - (UX(IJ, K) + UX(IJLS, K)) * AMH *
     &                (HXYU(IJ) + HXYU(IJLS)) * 0.25D0
     &              + (  VX(IJLE, K) + VX(IJSE, K)
     &                 - VX(IJLW, K) - VX(IJSW, K)) * AMH
     &                / (HXU(IJ) + HXU(IJLS)) * RX * 0.5D0
     &              - (VX(IJ, K) + VX(IJLS, K)) * AMH *
     &                (HYXU(IJ) + HYXU(IJLS)) * 0.25D0
            SYY(IJ) = (VX(IJ, K) - VX(IJLS, K)) * AMH *
     &                RY(IJ) / (HYU(IJ) + HYU(IJLS)) * 2.D0
     &              + (UX(IJ, K) + UX(IJLS, K)) * AMH *
     &                (HYXU(IJ) + HYXU(IJLS)) * 0.25D0
     &              - (  UX(IJLE, K) + UX(IJSE, K)
     &                 - UX(IJLW, K) - UX(IJSW, K)) * AMH
     &                / (HXU(IJ) + HXU(IJLS)) * RX * 0.5D0
     &              - (VX(IJ, K) + VX(IJLS, K)) * AMH *
     &                (HXYU(IJ) + HXYU(IJLS)) * 0.25D0
            SYX(IJ) = (VX(IJ, K) - VX(IJLW, K)) * AMH *
     &                RX / (HXU(IJ) + HXU(IJLW)) * 2.D0
     &              - (VX(IJ, K) + VX(IJLW, K)) * AMH *
     &                (HYXU(IJ) + HYXU(IJLW)) * 0.25D0
     &              + (  UX(IJLN, K) + UX(IJNW, K)
     &                 - UX(IJLS, K) - UX(IJSW, K)) * AMH
     &                / (HYU(IJ) + HYU(IJLW)) / (DY(IJ) + DY(IJLN))
     &              - (UX(IJ, K) + UX(IJLW, K)) * AMH *
     &                (HXYU(IJ) + HXYU(IJLW)) * 0.25D0
            SZX(IJ) = - UX(IJ, K) *
     &                  (  (AMV(IJ, K) + AMV(IJ, K+1)) * 0.5D0 / REA
     &                   + (AMV(IJ, K) - AMV(IJ, K+1)) * RZ(IJ, K))
            SZY(IJ) = - VX(IJ, K) *
     &                  (  (AMV(IJ, K) + AMV(IJ, K+1)) * 0.5D0 / REA
     &                   + (AMV(IJ, K) - AMV(IJ, K+1)) * RZ(IJ, K))
         END DO

         DO IJ = IJVSTR, IJVEND+1
            IJLW = IJ + LW
            FUX(IJ, K) = SXX(IJ) * (HYU(IJ) + HYU(IJLW)) *
     &                             (HYU(IJ) + HYU(IJLW)) * 0.25D0 *
     &                   AMFVX(IJ, K)
            FVX(IJ, K) = SYX(IJ) * (HYU(IJ) + HYU(IJLW)) *
     &                             (HYU(IJ) + HYU(IJLW)) * 0.25D0 *
     &                   AMFVX(IJ, K)
         END DO
         DO IJ = IJVSTR, IJVEND+NXDIM
            IJLS = IJ + LS
            FUY(IJ, K) = SXY(IJ) * (HXU(IJ) + HXU(IJLS)) *
     &                             (HXU(IJ) + HXU(IJLS)) * 0.25D0 *
     &                   AMFVY(IJ, K)
            FVY(IJ, K) = SYY(IJ) * (HXU(IJ) + HXU(IJLS)) *
     &                             (HXU(IJ) + HXU(IJLS)) * 0.25D0 *
     &                   AMFVY(IJ, K)
         END DO

         DO IJ = IJVSTR, IJVEND
            SMX(IJ, K) = SZX(IJ) / REA
            SMY(IJ, K) = SZY(IJ) / REA
         END DO

      END DO

      DO K = KSTR, KSTR+KZ-1
         DO IJ = IJVSTR, IJVEND
            GX(IJ, K) = (  GX(IJ, K)
     &                   + (  (FUX(IJ+LE, K) - FUX(IJ, K)) *
     &                        RX * RYU(IJ)
     &                      + (FUY(IJ+LN, K) - FUY(IJ, K)) *
     &                        RYM(IJ) * RXU(IJ)) * RXU(IJ) * RYU(IJ)
     &                   + SMX(IJ, K)) * AMSKV(IJ, K)
            GY(IJ, K) = (  GY(IJ, K)
     &                   + (  (FVX(IJ+LE, K) - FVX(IJ, K)) *
     &                        RX * RYU(IJ)
     &                      + (FVY(IJ+LN, K) - FVY(IJ, K)) *
     &                        RYM(IJ) * RXU(IJ)) * RXU(IJ) * RYU(IJ)
     &                   + SMY(IJ, K)) * AMSKV(IJ, K)
         END DO
      END DO
      DO K = KSTR+KZ, KEND
         DO IJ = IJVSTR, IJVEND
            GX(IJ, K) = (  GX(IJ, K)
     &                   + (  (FUX(IJ+LE, K) - FUX(IJ, K)) *
     &                        RX * RYU(IJ)
     &                      + (FUY(IJ+LN, K) - FUY(IJ, K)) *
     &                        RYM(IJ) * RXU(IJ)
     &                     ) * RXU(IJ) * RYU(IJ) * RZ(IJ, K)
     &                   + SMX(IJ, K)) * AMSKV(IJ, K)
            GY(IJ, K) = (  GY(IJ, K)
     &                   + (  (FVX(IJ+LE, K) - FVX(IJ, K)) *
     &                        RX * RYU(IJ)
     &                      + (FVY(IJ+LN, K) - FVY(IJ, K)) *
     &                        RYM(IJ) * RXU(IJ)
     &                     ) * RXU(IJ) * RYU(IJ) * RZ(IJ, K)
     &                   + SMY(IJ, K)) * AMSKV(IJ, K)
         END DO
      END DO

      DO K = KSTR, KEND
         DO IJ = IJVSTR, IJVEND
            XX(IJ, K) = GX(IJ, K)
            YY(IJ, K) = GY(IJ, K)
         END DO
      END DO

      RETURN
      END
#ifdef OPT_BBL
c *********************************************************************

      SUBROUTINE VSCVLB(
     M                 GX,     GY,     XX,     YY,
     I                 UY,     UX,     VY,     VX,
     I                AMV)
      use zocdim
      use zocgrd
      use zocmsk
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Visocity term of the BBL momentum eqs.
c
c  HISTORY
c     '09.01.06  T.Suzuki: bug fix (FUZ/FVZ added to GX/GY in VSCVLB)
c
c ---------------------------------------------------------------------

      REAL*8      UY(NXYDIM, NZDIM),     UX(NXYDIM, NZDIM)
      REAL*8      VY(NXYDIM, NZDIM),     VX(NXYDIM, NZDIM)
      REAL*8     AMV(NXYDIM, NZDIM)
      REAL*8      GX(NXYDIM, NZDIM),     GY(NXYDIM, NZDIM)
      REAL*8      XX(NXYDIM, NZDIM),     YY(NXYDIM, NZDIM)

      REAL*8      RZ(NXYDIM),    RZM(NXYDIM)
      SAVE RZ, RZM

      REAL*8     FUX(NXYDIM),    FVX(NXYDIM)
      REAL*8     FUY(NXYDIM),    FVY(NXYDIM)
      REAL*8     FUZ(NXYDIM),    FVZ(NXYDIM)
      REAL*8     SMX(NXYDIM),    SMY(NXYDIM)
      REAL*8     SXX(NXYDIM),    SYY(NXYDIM)
      REAL*8     SXY(NXYDIM),    SYX(NXYDIM)
      REAL*8     SZX(NXYDIM),    SZY(NXYDIM)
c      COMMON /WORK/ FUX, FUY, FUZ, FVX, FVY, FVZ, SMX, SMY,
c     &              SXX, SYY, SXY, SYX, SZX, SZY

      INTEGER     IJ,    KUP,      K
      INTEGER   IJLN,   IJLS,   IJLE,   IJLW
      INTEGER   IJNW,   IJSE,   IJSW
      INTEGER  IFPAR,  JFPAR

      LOGICAL OFIRST
      SAVE OFIRST

      REAL*8  AMHBBL
      SAVE AMHBBL
      NAMELIST /NMBBVH/ AMHBBL

      DATA OFIRST / .TRUE. /
      DATA AMHBBL / 0.D0 /

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF

      IF (OFIRST) THEN
         OFIRST = .FALSE.
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMBBVH, END=99)
 99      CONTINUE
         WRITE(JFPAR, NMBBVH)

         DO IJ = 1, NXYDIM
            RZ   (IJ) = 1.D0 / DZV(IJ, KEND)
            RZM  (IJ) = 1.D0 / DZM(IJ, KEND)
         END DO
      END IF

      K = KEND

      DO IJ = 1, NXYDIM
         FUX(IJ) = 0.D0
         FVX(IJ) = 0.D0
         FUY(IJ) = 0.D0
         FVY(IJ) = 0.D0
      END DO

      DO IJ = IJVSTR, IJVEND
         KUP = MAX(NBOTV(IJ)-1, 1)
         FUZ(IJ) = AMV(IJ, KEND) * RZM(IJ) *
     &             (UX(IJ, KUP) - UX(IJ, KEND))
         FVZ(IJ) = AMV(IJ, KEND) * RZM(IJ) *
     &             (VX(IJ, KUP) - VX(IJ, KEND))
      END DO

      DO IJ = IJSTR, IJEND+NXDIM+1
         IJLS = IJ + LS
         IJLW = IJ + LW
         IJLN = IJ + LN
         IJLE = IJ + LE
         IJNW = IJ + LNW
         IJSE = IJ + LSE
         IJSW = IJ + LSW
         SXX(IJ) = (UX(IJ, K) - UX(IJLW, K)) * AMHBBL *
     &             RX / (HXU(IJ) + HXU(IJLW)) * 2.D0
     &           + (VX(IJ, K) + VX(IJLW, K)) * AMHBBL *
     &             (HXYU(IJ) + HXYU(IJLW)) * 0.25D0
     &           - (  VX(IJLN, K) + VX(IJNW, K)
     &              - VX(IJLS, K) - VX(IJSW, K)) * AMHBBL
     &             / (HYU(IJ) + HYU(IJLW)) / (DY(IJ) + DY(IJLN))
     &           - (UX(IJ, K) + UX(IJLW, K)) * AMHBBL *
     &             (HYXU(IJ) + HYXU(IJLW)) * 0.25D0
         SXY(IJ) = (UX(IJ, K) - UX(IJLS, K)) * AMHBBL *
     &             RY(IJ) / (HYU(IJ) + HYU(IJLS)) * 2.D0
     &           - (UX(IJ, K) + UX(IJLS, K)) * AMHBBL *
     &             (HXYU(IJ) + HXYU(IJLS)) * 0.25D0
     &           + (  VX(IJLE, K) + VX(IJSE, K)
     &              - VX(IJLW, K) - VX(IJSW, K)) * AMHBBL
     &             / (HXU(IJ) + HXU(IJLS)) * RX * 0.5D0
     &           - (VX(IJ, K) + VX(IJLS, K)) * AMHBBL *
     &             (HYXU(IJ) + HYXU(IJLS)) * 0.25D0
         SYY(IJ) = (VX(IJ, K) - VX(IJLS, K)) * AMHBBL *
     &             RY(IJ) / (HYU(IJ) + HYU(IJLS)) * 2.D0
     &           + (UX(IJ, K) + UX(IJLS, K)) * AMHBBL *
     &             (HYXU(IJ) + HYXU(IJLS)) * 0.25D0
     &           - (  UX(IJLE, K) + UX(IJSE, K)
     &              - UX(IJLW, K) - UX(IJSW, K)) * AMHBBL
     &             / (HXU(IJ) + HXU(IJLS)) * RX * 0.5D0
     &           - (VX(IJ, K) + VX(IJLS, K)) * AMHBBL *
     &             (HXYU(IJ) + HXYU(IJLS)) * 0.25D0
         SYX(IJ) = (VX(IJ, K) - VX(IJLW, K)) * AMHBBL *
     &             RX / (HXU(IJ) + HXU(IJLW)) * 2.D0
     &           - (VX(IJ, K) + VX(IJLW, K)) * AMHBBL *
     &             (HYXU(IJ) + HYXU(IJLW)) * 0.25D0
     &           + (  UX(IJLN, K) + UX(IJNW, K)
     &              - UX(IJLS, K) - UX(IJSW, K)) * AMHBBL
     &             / (HYU(IJ) + HYU(IJLW)) / (DY(IJ) + DY(IJLN))
     &           - (UX(IJ, K) + UX(IJLW, K)) * AMHBBL *
     &             (HXYU(IJ) + HXYU(IJLW)) * 0.25D0
         SZX(IJ) = - UX(IJ, K) * AMV(IJ, K) / REA
         SZY(IJ) = - VX(IJ, K) * AMV(IJ, K) / REA
      END DO

      DO IJ = IJVSTR, IJVEND+1
         IJLW = IJ + LW
         FUX(IJ) = SXX(IJ) * (HYU(IJ) + HYU(IJLW)) *
     &                       (HYU(IJ) + HYU(IJLW)) * 0.25D0 *
     &             AMFVX(IJ, K)
         FVX(IJ) = SYX(IJ) * (HYU(IJ) + HYU(IJLW)) *
     &                       (HYU(IJ) + HYU(IJLW)) * 0.25D0 *
     &             AMFVX(IJ, K)
      END DO
      DO IJ = IJVSTR, IJVEND+NXDIM
         IJLS = IJ + LS
         FUY(IJ) = SXY(IJ) * (HXU(IJ) + HXU(IJLS)) *
     &                       (HXU(IJ) + HXU(IJLS)) * 0.25D0 *
     &             AMFVY(IJ, K)
         FVY(IJ) = SYY(IJ) * (HXU(IJ) + HXU(IJLS)) *
     &                       (HXU(IJ) + HXU(IJLS)) * 0.25D0 *
     &             AMFVY(IJ, K)
      END DO

      DO IJ = IJVSTR, IJVEND
         SMX(IJ) = SZX(IJ) / REA
         SMY(IJ) = SZY(IJ) / REA
      END DO

      DO IJ = IJVSTR, IJVEND
         GX(IJ, K) = (  GX(IJ, K)
     &                + FUZ(IJ) * RZ(IJ)
     &                + (  (FUX(IJ+LE) - FUX(IJ)) *
     &                     RX * RYU(IJ)
     &                   + (FUY(IJ+LN) - FUY(IJ)) *
     &                     RYM(IJ) * RXU(IJ)
     &                  ) * RXU(IJ) * RYU(IJ) * RZ(IJ)
     &                + SMX(IJ)) * AMSKVB(IJ)
         GY(IJ, K) = (  GY(IJ, K)
     &                + FVZ(IJ) * RZ(IJ)
     &                + (  (FVX(IJ+LE) - FVX(IJ)) *
     &                     RX * RYU(IJ)
     &                   + (FVY(IJ+LN) - FVY(IJ)) *
     &                     RYM(IJ) * RXU(IJ)
     &                  ) * RXU(IJ) * RYU(IJ) * RZ(IJ)
     &                + SMY(IJ)) * AMSKVB(IJ)
      END DO

      DO IJ = IJVSTR, IJVEND
         XX(IJ, KEND) = GX(IJ, KEND)
         YY(IJ, KEND) = GY(IJ, KEND)
      END DO

      RETURN
      END
#endif
