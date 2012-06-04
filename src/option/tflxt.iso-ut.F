      SUBROUTINE FLXTRC(
     O                ADT,  DIFFZ,
     I                 TX,     HX,
     I                 TY,     HZ,
     I                 UY,     VY,
     I                  W,    AHV)
      use zocdim
      use zocgrd
      use zocmsk
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  HISTORY
c     '03.05.12  H.Hasumi: from COCO3.4
c     '07.04.23  H.Hasumi
c     '07.05.01  H.Hasumi: McDougall et al. (2003) eq. of state
c     '07.09.25  H.Hasumi: arguments of CHEKIN
c     '08.06.11  H.Hasumi: initial/final setup
c     '08.07.10  H.Hasumi: initial/final setup
c
c ---------------------------------------------------------------------

      REAL*8     ADT(NXYDIM, NZDIM, NTDIM)    
      REAL*8   DIFFZ(NXYDIM, NZDIM)
      REAL*8      TX(NXYDIM, NZDIM, NTDIM)
      REAL*8      TY(NXYDIM, NZDIM, NTDIM)
      REAL*8      HX(NXYDIM),     HZ(NXYDIM)
      REAL*8      UY(NXYDIM, NZDIM),     VY(NXYDIM, NZDIM)
      REAL*8       W(NXYDIM, NZDIM),    AHV(NXYDIM, NZDIM)

      REAL*8  RYMLLS(NXYDIM), RYMRYM(NXYDIM)
      REAL*8     DX2,    RX2
      REAL*8     EPS
      LOGICAL OFIRST
      SAVE RYMLLS, RYMRYM, DX2, RX2, EPS, OFIRST

      REAL*8     FTX(NXYDIM, NZDIM, NTDIM)
      REAL*8     FTY(NXYDIM, NZDIM, NTDIM)
      REAL*8     FTZ(NXYDIM, NZDIM, NTDIM)
      REAL*8     WZC(NXYDIM, NZDIM),    RZM(NXYDIM, NZDIM)
c      COMMON /WORK/ FTX, FTY, FTZ, WZC, RZM
c      COMMON /WORK/ WZC, RZM
      COMMON /WORKFT/ FTX, FTY, FTZ

      REAL*8   HZBOT(NXYDIM)

      REAL*8   XDZDX(NXYDIM, NZDIM),  YDZDY(NXYDIM, NZDIM)
      REAL*8   ZDZDX(NXYDIM, NZDIM),  ZDZDY(NXYDIM, NZDIM)
      REAL*8   XDTDZ(NXYDIM, NZDIM, NTDIM),  YDTDZ(NXYDIM, NZDIM, NTDIM)
      REAL*8   ZDTDX(NXYDIM, NZDIM, NTDIM),  ZDTDY(NXYDIM, NZDIM, NTDIM)
      REAL*8     TNT(NXYDIM),    TET(NXYDIM)
      REAL*8     TST(NXYDIM),    TS2(NXYDIM)
      REAL*8     TWT(NXYDIM),    TWW(NXYDIM)
      REAL*8     TNW(NXYDIM),    TSW(NXYDIM)
      REAL*8     TSE(NXYDIM)
      REAL*8       U,      V
      REAL*8     CXL,    CYL,    CXD,    CYD,    WUT
      REAL*8      A0,     A1,     A2
      REAL*8     C00,    C01,    C02,    C10,    C11,    C20
      REAL*8     D00,    D01,    D02,    D10,    D11,    D20
      INTEGER     IJ,      K,   IJLS,   IJLW,      N
      INTEGER    KUU,     KU,     KD
      INTEGER  IFPAR,  JFPAR

      REAL*8    TADV
      REAL*8    VPOS,   VNEG
      REAL*8   RZMUP,   DZCN,  RZMCN
      REAL*8   RYMUP,   DYCN,  RYMCN
      REAL*8     TUP,    TCN,    TDN
      REAL*8    TCNE,   TCNW,   TCNN,   TCNS
      REAL*8   TDIFF,  TCURV
      REAL*8    TREF,  TREF1,  TREF2
      REAL*8    TMIN,   TMAX
      REAL*8    SLP1,   SLP2

      REAL*8   TADVX(NXYDIM, NTDIM),  TADVY(NXYDIM, NTDIM)
      COMMON /UTADV/ TADVX, TADVY

c      REAL*8  PSIGMX(NXYDIM, NZDIM), PSIGMY(NXYDIM, NZDIM)

      REAL*8     AHH,    AHI,    AHG
      SAVE       AHH,    AHI,    AHG
      NAMELIST /NMDIFH/ AHH
      NAMELIST /NMDIFI/ AHI
      NAMELIST /NMDIFG/ AHG

      DATA OFIRST / .TRUE. /
      DATA AHH / 0.D0 /
      DATA AHI / 0.D0 /
      DATA AHG / 0.D0 /

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF

      IF (OFIRST) THEN
         OFIRST = .FALSE.
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMDIFH, END=97)
 97      CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMDIFI, END=96)
 96      CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMDIFG, END=95)
 95      CONTINUE
         WRITE(JFPAR, NMDIFH)
         WRITE(JFPAR, NMDIFI)
         WRITE(JFPAR, NMDIFG)

         DO 40 IJ = IJSTR-NXDIM, IJEND+NXDIM
            RYMLLS(IJ) = 1.D0 / DYM(IJ) / DYM(IJ-NXDIM)
            RYMRYM(IJ) = RYMLLS(IJ) / (DYM(IJ) + DYM(IJ-NXDIM))
 40      CONTINUE
         DX2 = DX * DX
         RX2 = 1.D0 / DX2
         EPS = 1.D-10
      END IF

      CALL DNSGRD(
     O             XDZDX,  YDZDY,  ZDZDX,  ZDZDY,
     O             XDTDZ,  YDTDZ,  ZDTDX,  ZDTDY,
     I                TY,     TX,     
     W               FTX,    FTY,    FTZ)

      DO 120 N = 1, NTDIM
         DO 110 K = 1, NZDIM
            DO 100 IJ = 1, NXYDIM
               ADT(IJ, K, N) = 0.D0
               FTX(IJ, K, N) = 0.D0
               FTY(IJ, K, N) = 0.D0
               FTZ(IJ, K, N) = 0.D0
  100       CONTINUE
  110    CONTINUE
  120 CONTINUE
      DO K = 1, NZDIM
         DO IJ = 1, NXYDIM
            DIFFZ(IJ, K) = 0.D0
         END DO
      END DO

      DO IJ = 1, NXYDIM
         HZBOT(IJ) = HZ(IJ) + ZBOT
      END DO
      DO K = KSTR, KSTR+KZ-1
         DO IJ = 1, NXYDIM
            WZC(IJ, K) = W(IJ, K) * HZBOT(IJ)
            RZM(IJ, K) = 1.D0 / DSM(K) / HZBOT(IJ)
         END DO
      END DO
      DO K = KSTR+KZ, KEND
         DO IJ = 1, NXYDIM
            WZC(IJ, K) = W(IJ, K)
            RZM(IJ, K) = 1.D0 / DZM(IJ, K)
         END DO
      END DO

      CALL CHEKIN(WZC, 'WZC', NX, NY, NZ, NXYZDM, 'OCN')

      DO 230 N = 1, NTDIM
         DO 220 K = KSTR+1, KEND
            KUU = K - 2
            KU  = K - 1
            KD  = K + 1
            DO 210 IJ = IJTSTR, IJTEND
               DIFFZ(IJ, K) =
     &         (  AHV(IJ, K)
     &          + AHI * (  ZDZDX(IJ, K) * ZDZDX(IJ, K)
     &                   + ZDZDY(IJ, K) * ZDZDY(IJ, K))) *
     &         RZM(IJ, K) * AMFTZ(IJ, K)
               WUT = WZC(IJ, K) * TS
               
               VPOS = 0.5D0 + SIGN(0.5D0, WUT)
               VNEG = 0.5D0 - SIGN(0.5D0, WUT)
               RZMUP = VPOS / DZM(IJ, KD) + VNEG / DZM(IJ, KU)
               DZCN  = VPOS * DZ(IJ, K) + VNEG * DZ(IJ, KU)
               RZMCN = VPOS / (DZM(IJ, K) + DZM(IJ, KD))
     &               + VNEG / (DZM(IJ, K) + DZM(IJ, KU))

               TUP = VPOS * TX(IJ, KD, N) + VNEG * TX(IJ, KUU, N)
               TCN = VPOS * TX(IJ, K , N) + VNEG * TX(IJ, KU,  N)
               TDN = VPOS * TX(IJ, KU, N) + VNEG * TX(IJ, K ,  N)
               TDIFF = TDN - TUP
               TCURV = (TDN - TCN) / DZM(IJ, K) - (TCN - TUP) * RZMUP

               A2 = TCURV * RZMCN
               A1 = (TX(IJ, KU, N) - TX(IJ, K, N)) / DZM(IJ, K)
     &            + 0.5D0 * (DZ(IJ, K) - DZ(IJ, KU)) * A2
               A0 = (  DZ(IJ, K ) * TX(IJ, KU, N)
     &               + DZ(IJ, KU) * TX(IJ, K , N))
     &              / DZM(IJ, K) * 0.5D0
     &            - DZ(IJ, KU) * DZ(IJ, K) * 0.25D0 * A2
               TADV = (  WUT * WUT / 3.D0
     &                 - DZ(IJ, K) * DZ(IJ, K) / 12.D0) * A2
     &              - WUT * 0.5D0 * A1 + A0
               TCURV = ABS(TCURV) * DZCN

               IF (     (TCURV .GT. ABS(TDIFF))
     &              .OR. (ABS(WUT) .LT. EPS)) THEN
                  TADV = TCN
               ELSE
                  TREF = TUP + (TCN - TUP) / ABS(WUT) * DZM(IJ, K)
                  VPOS = 0.5D0 + SIGN(0.5D0, TDIFF)
                  VNEG = 0.5D0 - SIGN(0.5D0, TDIFF)
                  TMIN = VPOS * TCN + VNEG * MAX(TDN, TREF)
                  TMAX = VPOS * MIN(TDN, TREF) + VNEG * TCN
                  TADV = MIN(MAX(TADV, TMIN), TMAX)
               END IF

               FTZ(IJ, K, N) = 
     &         (  DIFFZ(IJ, K) * (TX(IJ, KU, N) - TX(IJ, K, N))
     &          - (AHI + AHG) *
     &            (  ZDZDX(IJ, K) * ZDTDX(IJ, K, N)
     &             + ZDZDY(IJ, K) * ZDTDY(IJ, K, N))
     &          - WZC(IJ, K) * TADV) * AMFTZ(IJ, K)
  210       CONTINUE
  220    CONTINUE
  230 CONTINUE

      DO 340 N = 1, NTDIM
         DO 330 K = KSTR, KEND
            DO 310 IJ = IJTSTR, IJTEND+NXDIM
               TS2(IJ) = TX(IJ+LSS, K, N) * AMSKT(IJ+LSS, K)
     &                 + TX(IJ+LS , K, N) * (1.D0 - AMSKT(IJ+LSS, K))
               TSE(IJ) = TX(IJ+LSE, K, N) * AMSKT(IJ+LSE, K)
     &                 + TX(IJ+LS , K, N) * (1.D0 - AMSKT(IJ+LSE, K))
               TSW(IJ) = TX(IJ+LSW, K, N) * AMSKT(IJ+LSW, K)
     &                 + TX(IJ+LS , K, N) * (1.D0 - AMSKT(IJ+LSW, K))
               TNT(IJ) = TX(IJ+LN , K, N) * AMSKT(IJ+LN, K)
     &                 + TX(IJ    , K, N) * (1.D0 - AMSKT(IJ+LN, K))
               TET(IJ) = TX(IJ+LE , K, N) * AMSKT(IJ+LE, K)
     &                 + TX(IJ    , K, N) * (1.D0 - AMSKT(IJ+LE, K))
               TWT(IJ) = TX(IJ+LW , K, N) * AMSKT(IJ+LW, K)
     &                 + TX(IJ    , K, N) * (1.D0 - AMSKT(IJ+LW, K))
  310       CONTINUE
            
            DO 320 IJ = IJTSTR, IJTEND+NXDIM

               IJLS = IJ + LS
               U = 0.5D0 *
     &             (  UY(IJ+LSW, K) * DZV(IJ+LSW, K) * HYU(IJ+LSW)
     &              + UY(IJLS  , K) * DZV(IJLS  , K) * HYU(IJLS  ))
               V = 0.5D0 *
     &             (  VY(IJ+LSW, K) * DZV(IJ+LSW, K) * HXU(IJ+LSW)
     &              + VY(IJLS  , K) * DZV(IJLS  , K) * HXU(IJLS  ))
               CXD = U * TS * 4.D0
     &               / (HXT(IJ) + HXT(IJLS)) / (HYT(IJ) + HYT(IJLS))
     &               / MAX(EPS, AMFTY(IJ, K))
               CYD = V * TS * 4.D0
     &               / (HXT(IJ) + HXT(IJLS)) / (HYT(IJ) + HYT(IJLS))
     &               / MAX(EPS, AMFTY(IJ, K))
               
               VPOS = 0.5D0 + SIGN(0.5D0, CYD)
               VNEG = 0.5D0 - SIGN(0.5D0, CYD)
               RYMUP = VPOS / DYM(IJ+LSS) + VNEG / DYM(IJ)
               RYMCN = VPOS / (DYM(IJLS) + DYM(IJ+LSS))
     &               + VNEG / (DYM(IJLS) + DYM(IJ))
               DYCN  = VPOS * DY(IJLS) + VNEG * DY(IJ)
               
               TUP  = VPOS * TS2(IJ)        + VNEG * TNT(IJ)
               TCN  = VPOS * TX(IJLS, K, N) + VNEG * TX(IJ, K, N)
               TDN  = VPOS * TX(IJ, K, N)   + VNEG * TX(IJLS, K, N)
               TCNE = VPOS * TSE(IJ)        + VNEG * TET(IJ)
               TCNW = VPOS * TSW(IJ)        + VNEG * TWT(IJ)
               TDIFF = TDN - TUP
               TCURV = (TDN - TCN) * RYM(IJLS) - (TCN - TUP) * RYMUP

               D02 = TCURV * RYMCN
               D00 = (  TX(IJ, K, N) * DY(IJLS)
     &                + TX(IJLS, K, N) * DY(IJ)) * RYM(IJLS) * 0.5D0
     &             - D02 * DY(IJ) * DY(IJLS) * 0.25D0
               D01 = (TX(IJ, K, N) - TX(IJLS, K, N)) * RYM(IJLS)
     &             - D02 * (DY(IJ) - DY(IJLS)) * 0.5D0
               D20 = (TCNE - 2.D0 * TCN + TCNW) * RX2 * 0.5D0
               D11 = (0.5D0 + SIGN(0.5D0, U)) *
     &               (  TSW(IJ) - TX(IJLS, K, N)
     &                - TWT(IJ) + TX(IJ, K, N)) * RX * RYM(IJLS)
     &              + (0.5D0 - SIGN(0.5D0, U)) *
     &                (  TET(IJ) - TX(IJ, K, N)
     &                 - TSE(IJ) + TX(IJLS, K, N)) * RX * RYM(IJLS)
               D10 = (0.5D0 + SIGN(0.5D0, U * V)) *
     &               (  (TX(IJ, K, N) - TWT(IJ)) * RX
     &                + D20 * DX
     &                - D11 * DY(IJ) * 0.5D0)
     &             + (0.5D0 - SIGN(0.5D0, U * V)) *
     &               (  (TET(IJ) - TX(IJ, K, N)) * RX
     &                - D20 * DX
     &                - D11 * DY(IJ) * 0.5D0)
               TADVY(IJ, N) = D00
     &              - D01 * 0.5D0 * CYD / 12.D0
     &              + D02 * (CYD * CYD / 3.D0 - DY(IJ) * DY(IJ) / 12.D0)
     &              - D10 * 0.5D0 * CXD
     &              + D20 * 0.5D0 * CXD * CXD
     &              + D11 * 0.25D0 * CXD * CYD
               TCURV = ABS(TCURV) * DYCN
               
               CXD = CXD * RX
               CYD = CYD * RYM(IJLS)
               
               D10 = (  (0.5D0 + SIGN(0.5D0, CXD)) * (TCN - TCNW)
     &                + (0.5D0 - SIGN(0.5D0, CXD)) * (TCNE - TCN))
               TADVY(IJ, N) = TADVY(IJ, N) + 0.5D0 * CXD * D10
               IF (     (TCURV .GT. ABS(TDIFF))
     &             .OR. (      (ABS(CXD) .LT. EPS)
     &                   .AND. (ABS(CYD) .LT. EPS))) THEN
                  TADVY(IJ, N) = TCN
               ELSE
                  SLP1 = (1.D0 + ABS(CXD * CYD)) 
     &                   / (ABS(CXD) + ABS(CYD))
                  SLP2 = ABS(CXD)
                  TREF1 = TUP + (TCN - TUP) * SLP1
                  TREF2 = TDN + (TCN - TDN) * SLP2
                  VPOS = 0.5D0 + SIGN(0.5D0, TDIFF)
                  VNEG = 0.5D0 - SIGN(0.5D0, TDIFF)
                  TMIN = VPOS * TCN + VNEG * MAX(TDN, TREF2)
                  TMAX = VPOS * MIN(TDN, TREF1) + VNEG * TCN
                  TADVY(IJ, N) = MIN(MAX(TADVY(IJ, N), TMIN), TMAX)
                  TADVY(IJ, N) = TADVY(IJ, N) - 0.5D0 * CXD * D10
               END IF

               FTY(IJ, K, N) =
     &           (  (AHH + AHI) * RYM(IJLS) / (HYT(IJ) + HYT(IJLS)) *
     &              (TX(IJ, K, N) - TX(IJLS, K, N)) * 2.D0
     &            - (AHI - AHG) * YDZDY(IJ, K) * YDTDZ(IJ, K, N)) *
     &           (HXU(IJLS) + HXU(IJ+LSW)) * 0.5D0 * AMFTY(IJ, K)
     &         - V * TADVY(IJ, N)
  320       CONTINUE
  330    CONTINUE

  340 CONTINUE

      DO 440 N = 1, NTDIM
         DO 430 K = KSTR, KEND
            DO 410 IJ = IJTSTR, IJTEND+1
               TWW(IJ) = TX(IJ+LWW, K, N) * AMSKT(IJ+LWW, K)
     &                 + TX(IJ+LW , K, N) * (1.D0 - AMSKT(IJ+LWW, K))
               TNW(IJ) = TX(IJ+LNW, K, N) * AMSKT(IJ+LNW, K)
     &                 + TX(IJ+LW , K, N) * (1.D0 - AMSKT(IJ+LNW, K))
               TSW(IJ) = TX(IJ+LSW, K, N) * AMSKT(IJ+LSW, K)
     &                 + TX(IJ+LW , K, N) * (1.D0 - AMSKT(IJ+LSW, K))
               TET(IJ) = TX(IJ+LE , K, N) * AMSKT(IJ+LE, K)
     &                 + TX(IJ    , K, N) * (1.D0 - AMSKT(IJ+LE, K))
               TNT(IJ) = TX(IJ+LN , K, N) * AMSKT(IJ+LN, K)
     &                 + TX(IJ    , K, N) * (1.D0 - AMSKT(IJ+LN, K))
               TST(IJ) = TX(IJ+LS , K, N) * AMSKT(IJ+LS, K)
     &                 + TX(IJ    , K, N) * (1.D0 - AMSKT(IJ+LS, K))
  410       CONTINUE

            DO 420 IJ = IJTSTR, IJTEND+1

               IJLW = IJ + LW
               IJLS = IJ + LS
               U   = 0.5D0 *
     &               (  UY(IJLW  , K) * DZV(IJLW  , K) * HYU(IJLW)
     &                + UY(IJ+LSW, K) * DZV(IJ+LSW, K) * HYU(IJ+LSW))
               V   = 0.5D0 *
     &               (  VY(IJLW  , K) * DZV(IJLW  , K) * HXU(IJLW)
     &                + VY(IJ+LSW, K) * DZV(IJ+LSW, K) * HXU(IJ+LSW))
               CXL = U * TS * 4.D0
     &               / (HXT(IJ) + HXT(IJLW)) / (HYT(IJ) + HYT(IJLW))
     &               / MAX(EPS, AMFTX(IJ, K))
               CYL = V * TS * 4.D0
     &               / (HXT(IJ) + HXT(IJLW)) / (HYT(IJ) + HYT(IJLW))
     &               / MAX(EPS, AMFTX(IJ, K))

               VPOS = 0.5D0 + SIGN(0.5D0, CXL)
               VNEG = 0.5D0 - SIGN(0.5D0, CXL)

               TUP  = VPOS * TWW(IJ)        + VNEG * TET(IJ)
               TCN  = VPOS * TX(IJLW, K, N) + VNEG * TX(IJ, K, N)
               TDN  = VPOS * TX(IJ, K, N)   + VNEG * TX(IJLW, K, N)
               TCNN = VPOS * TNW(IJ)        + VNEG * TNT(IJ)
               TCNS = VPOS * TSW(IJ)        + VNEG * TST(IJ)
               TDIFF = TDN - TUP
               TCURV = TDN - 2.D0 * TCN + TUP
               
               C10 = (TX(IJ, K, N) - TX(IJLW, K, N)) * RX
               C20 = 0.5D0 * TCURV
               C00 = (TX(IJ, K, N) + TX(IJLW, K, N)) * 0.5D0
     &             - C20 * 0.25D0
               C20 = C20 * RX2
               C11 = (0.5D0 + SIGN(0.5D0, V)) *
     &               (  TSW(IJ) - TX(IJLW, K, N)
     &                - TST(IJ) + TX(IJ, K, N)) * RX * RYM(IJLS)
     &             + (0.5D0 - SIGN(0.5D0, V)) *
     &               (  TNT(IJ) - TX(IJ, K, N)
     &                - TNW(IJ) + TX(IJLW, K, N)) * RX * RYM(IJ)
               C02 = (  (TCNN - TCN) * DYM(IJLS)
     &                + (TCNS - TCN) * DYM(IJ)) * RYMRYM(IJ)
               C01 = (0.25D0 + SIGN(0.25D0, U * V)) *
     &               (  (  (TNW(IJ) - TX(IJLW, K, N)) * DYM(IJLS)
     &                   - (TST(IJ) - TX(IJ, K, N)) * DYM(IJ)) *
     &                  RYMLLS(IJ)
     &                - C02 * (DYM(IJ) - DYM(IJLS)))
     &              + (0.25D0 - SIGN(0.25D0, U * V)) *
     &                (  (  (TNT(IJ) - TX(IJ, K, N)) * DYM(IJLS)
     &                    - (TSW(IJ) - TX(IJLW, K, N)) * DYM(IJ)) *
     &                   RYMLLS(IJ)
     &                 - C02 * (DYM(IJ) - DYM(IJLS)))
               TADVX(IJ, N) = C00
     &              - C10 * 0.5D0 * CXL
     &              + C20 * (CXL * CXL / 3.D0 - DX2 / 12.D0)
     &              - C01 * 0.5D0 * CYL
     &              + C02 * 0.5D0 * CYL * CYL
     &              + C11 * 0.25D0 * CXL * CYL
               
               CXL = CXL * RX
               CYL = CYL * RY(IJ)
               
               C01 = (0.5D0 + SIGN(0.5D0, CYL)) * (TCN - TCNS)
     &             + (0.5D0 - SIGN(0.5D0, CYL)) * (TCNN - TCN)
               TADVX(IJ, N) = TADVX(IJ, N) + 0.5D0 * CYL * C01
               IF (     (ABS(TCURV) .GT. ABS(TDIFF))
     &             .OR. (      (ABS(CXL) .LT. EPS)
     &                   .AND. (ABS(CYL) .LT. EPS))) THEN
                  TADVX(IJ, N) = TCN
               ELSE
                  SLP1 = (1.D0 + ABS(CXL * CYL))
     &                   / (ABS(CXL) + ABS(CYL))
                  SLP2 = ABS(CYL)
                  TREF1 = TUP + (TCN - TUP) * SLP1
                  TREF2 = TDN + (TCN - TDN) * SLP2
                  VPOS = 0.5D0 + SIGN(0.5D0, TDIFF)
                  VNEG = 0.5D0 - SIGN(0.5D0, TDIFF)
                  TMIN = VPOS * TCN + VNEG * MAX(TDN, TREF2)
                  TMAX = VPOS * MIN(TDN, TREF1) + VNEG * TCN
                  TADVX(IJ, N) = MIN(MAX(TADVX(IJ, N), TMIN), TMAX)
                  TADVX(IJ, N) = TADVX(IJ, N) - 0.5D0 * CYL * C01
               END IF
               
               FTX(IJ, K, N) =
     &           (  (AHH + AHI) * RX / (HXT(IJ) + HXT(IJLW)) *
     &              (TX(IJ, K, N) - TX(IJLW, K, N)) * 2.D0
     &            - (AHI - AHG) * XDZDX(IJ, K) * XDTDZ(IJ, K, N)) *
     &           (HYU(IJLW) + HYU(IJ+LSW)) * 0.5D0 * AMFTX(IJ, K)
     &         - U * TADVX(IJ, N)
  420       CONTINUE
  430    CONTINUE
         
  440 CONTINUE

      DO 630 N = 1, NTDIM
         DO 620 K = KSTR, KEND
            DO 610 IJ = IJTSTR, IJTEND
               ADT(IJ, K, N) =
     &         (  (  (FTX(IJ+LE, K, N) - FTX(IJ, K, N)) * RX
     &             + (FTY(IJ+LN, K, N) - FTY(IJ, K, N)) * RY(IJ)) *
     &            RXT(IJ) * RYT(IJ)
     &          + FTZ(IJ, K, N) - FTZ(IJ, K+1, N)) / DZ(IJ, K)
  610       CONTINUE
  620    CONTINUE
  630 CONTINUE

c      DO K = 1, NZDIM
c         DO IJ = 1, NXYDIM
c            PSIGMX(IJ, K) = 0.0D0
c            PSIGMY(IJ, K) = 0.0D0
c         END DO
c      END DO
c
c      DO K = KSTR, KEND
c         DO IJ = IJTSTR, IJTEND+NXDIM
c            PSIGMX(IJ, K) = AHG * YDZDY(IJ, K)
c            PSIGMY(IJ, K) = - AHG * XDZDX(IJ, K)
c         END DO
c      END DO
c
c      CALL CHEKIN(PSIGMX, 'PSIGMX', NX, NY, NZ, NXYZDM, 'OCN')
c      CALL CHEKIN(PSIGMY, 'PSIGMY', NX, NY, NZ, NXYZDM, 'OCN')

      RETURN
      END

c *********************************************************************

      SUBROUTINE DNSGRD(
     O              XDZDX,  YDZDY,  ZDZDX,  ZDZDY,
     O              XDTDZ,  YDTDZ,  ZDTDX,  ZDTDY,
     I                 TY,     TX,     
     W               DTDX,   DTDY,   DTDZ)
      use zocdim
      use zocgrd
      use zocmsk
      IMPLICIT NONE

      REAL*8   XDZDX(NXYDIM, NZDIM),  YDZDY(NXYDIM, NZDIM)
      REAL*8   ZDZDX(NXYDIM, NZDIM),  ZDZDY(NXYDIM, NZDIM)
      REAL*8   XDTDZ(NXYDIM, NZDIM, NTDIM),  YDTDZ(NXYDIM, NZDIM, NTDIM)
      REAL*8   ZDTDX(NXYDIM, NZDIM, NTDIM),  ZDTDY(NXYDIM, NZDIM, NTDIM)
      REAL*8      TY(NXYDIM, NZDIM, NTDIM)
      REAL*8      TX(NXYDIM, NZDIM, NTDIM)
      REAL*8    DTDX(NXYDIM, NZDIM, NTDIM),   DTDY(NXYDIM, NZDIM, NTDIM)
      REAL*8    DTDZ(NXYDIM, NZDIM, NTDIM)

      REAL*8 C0(NZDIM), C1(NZDIM), C2(NZDIM)
      REAL*8 C3(NZDIM), C4(NZDIM), C5(NZDIM), C6(NZDIM)
      REAL*8 D0(NZDIM), D1(NZDIM), D2(NZDIM), D3(NZDIM), D4(NZDIM)
      REAL*8 D5(NZDIM), D6(NZDIM), D7(NZDIM), D8(NZDIM), D9(NZDIM)
      REAL*8 EPS
      LOGICAL OFIRST
      SAVE C0, C1, C2, C3, C4, C5, C6,
     &     D0, D1, D2, D3, D4, D5, D6, D7, D8, D9
      SAVE EPS, OFIRST

      REAL*8      P1,     P2
      REAL*8      TL,     SL
      REAL*8      RL,    RLW,    RLS,    RLU
      REAL*8    DZDX,   DZDY
      INTEGER     IJ,      K,      N
      INTEGER  IFPAR,  JFPAR

      REAL*8  SLPMAX
      SAVE SLPMAX
      NAMELIST /NMSLPM/ SLPMAX

      DATA OFIRST / .TRUE. /
      DATA SLPMAX / 1.D-2 /

      IF (OFIRST) THEN
         OFIRST = .FALSE.
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMSLPM, END=99)
 99      CONTINUE
         WRITE(JFPAR, NMSLPM)

         CALL SECOEF(
     O            C0(KSTR), C1(KSTR), C2(KSTR), C3(KSTR),
     O            C4(KSTR), C5(KSTR), C6(KSTR),
     O            D0(KSTR), D1(KSTR), D2(KSTR), D3(KSTR),  D4(KSTR),
     O            D5(KSTR), D6(KSTR), D7(KSTR), D8(KSTR),  D9(KSTR))
         EPS = 1.D-20
      END IF

      DO 115 N = 1, NTDIM
         DO 110 K = 1, NZDIM
            DO 100 IJ = 1, NXYDIM
               DTDX(IJ, K, N) = 0.D0
               DTDY(IJ, K, N) = 0.D0
               DTDZ(IJ, K, N) = 0.D0
  100       CONTINUE
  110    CONTINUE
  115 CONTINUE

      DO 200 K = KSTR, KEND
         DO 210 IJ = IJTSTR, IJTEND+NXDIM
            TL = TY(IJ, K, 1)
            SL = TY(IJ, K, 2)
            P1 = C0(K)
     &         + (C1(K) + (C2(K) + C3(K) * TL) * TL) * TL
     &         + (C4(K) + C5(K) * TL + C6(K) * SL) * SL
            P2 = D0(K)
     &         + (D1(K) + (D2(K) + (D3(K) + D4(K) * TL) * TL) * TL) * TL
     &         + (D5(K) + (D6(K) + D7(K) * TL * TL) * TL
     &                  + (D8(K) + D9(K) * TL * TL) * SQRT(SL)) * SL
            RL = P1 / P2

            TL = TY(IJ+LW, K, 1)
            SL = TY(IJ+LW, K, 2)
            P1 = C0(K)
     &         + (C1(K) + (C2(K) + C3(K) * TL) * TL) * TL
     &         + (C4(K) + C5(K) * TL + C6(K) * SL) * SL
            P2 = D0(K)
     &         + (D1(K) + (D2(K) + (D3(K) + D4(K) * TL) * TL) * TL) * TL
     &         + (D5(K) + (D6(K) + D7(K) * TL * TL) * TL
     &                  + (D8(K) + D9(K) * TL * TL) * SQRT(SL)) * SL
            RLW = P1 / P2

            TL = TY(IJ+LS, K, 1)
            SL = TY(IJ+LS, K, 2)
            P1 = C0(K)
     &         + (C1(K) + (C2(K) + C3(K) * TL) * TL) * TL
     &         + (C4(K) + C5(K) * TL + C6(K) * SL) * SL
            P2 = D0(K)
     &         + (D1(K) + (D2(K) + (D3(K) + D4(K) * TL) * TL) * TL) * TL
     &         + (D5(K) + (D6(K) + D7(K) * TL * TL) * TL
     &                  + (D8(K) + D9(K) * TL * TL) * SQRT(SL)) * SL
            RLS = P1 / P2

            DTDX(IJ, K, 1) = (RL - RLW) * RX * 2.D0
     &                       / (HXT(IJ) + HXT(IJ+LW)) *
     &                       AMSKT(IJ, K) * AMSKT(IJ+LW, K)
            DTDY(IJ, K, 1) = (RL - RLS) * RYM(IJ+LS) * 2.D0
     &                       / (HYT(IJ) + HYT(IJ+LS)) *
     &                       AMSKT(IJ, K) * AMSKT(IJ+LS, K)
  210    CONTINUE
  200 CONTINUE

      DO 230 K = KSTR+1, KEND
         DO 240 IJ = IJTSTR-NXDIM, IJTEND+NXDIM
            TL = TY(IJ, K, 1)
            SL = TY(IJ, K, 2)
            P1 = C0(K)
     &         + (C1(K) + (C2(K) + C3(K) * TL) * TL) * TL
     &         + (C4(K) + C5(K) * TL + C6(K) * SL) * SL
            P2 = D0(K)
     &         + (D1(K) + (D2(K) + (D3(K) + D4(K) * TL) * TL) * TL) * TL
     &         + (D5(K) + (D6(K) + D7(K) * TL * TL) * TL
     &                  + (D8(K) + D9(K) * TL * TL) * SQRT(SL)) * SL
            RL = P1 / P2

            TL = TY(IJ, K-1, 1)
            SL = TY(IJ, K-1, 2)
            P1 = C0(K)
     &         + (C1(K) + (C2(K) + C3(K) * TL) * TL) * TL
     &         + (C4(K) + C5(K) * TL + C6(K) * SL) * SL
            P2 = D0(K)
     &         + (D1(K) + (D2(K) + (D3(K) + D4(K) * TL) * TL) * TL) * TL
     &         + (D5(K) + (D6(K) + D7(K) * TL * TL) * TL
     &                  + (D8(K) + D9(K) * TL * TL) * SQRT(SL)) * SL
            RLU = P1 / P2

            DTDZ(IJ, K, 1) = MIN((RLU - RL) / DZM(IJ, K), 0.D0) *
     &                       AMFTZ(IJ, K)
  240    CONTINUE
  230 CONTINUE

      DO 300 K = KSTR, KEND
         DO 310 IJ = IJTSTR, IJTEND+NXDIM
            DZDX = DTDX(IJ, K, 1) * 4.D0
     &             / (  DTDZ(IJ, K  , 1) + DTDZ(IJ+LW, K ,  1)
     &                + DTDZ(IJ, K+1, 1) + DTDZ(IJ+LW, K+1, 1) - EPS)
            XDZDX(IJ, K) = MIN(SLPMAX, MAX(-SLPMAX, DZDX))
            DZDY = DTDY(IJ, K, 1) * 4.D0 * DYM(IJ+LS)
     &          / (  (DTDZ(IJ, K, 1) + DTDZ(IJ, K+1, 1)) * DY(IJ)
     &          + (DTDZ(IJ+LS, K, 1) + DTDZ(IJ+LS, K+1, 1)) * DY(IJ+LS)
     &                - EPS)
            YDZDY(IJ, K) = MIN(SLPMAX, MAX(-SLPMAX, DZDY))
  310    CONTINUE
  300 CONTINUE

      DO 330 K = KSTR+1, KEND
         DO 340 IJ = IJTSTR, IJTEND
            DZDX = ((DTDX(IJ, K-1, 1)+DTDX(IJ+LE, K-1, 1))*DZ(IJ, K-1)
     &            + (DTDX(IJ, K  , 1)+DTDX(IJ+LE, K  , 1))*DZ(IJ, K)) *
     &             0.25D0 / DZM(IJ, K) / (DTDZ(IJ, K, 1) - EPS)
            ZDZDX(IJ, K) = MIN(SLPMAX, MAX(-SLPMAX, DZDX)) *
     &                     AMFTZ(IJ, K)
            DZDY = 
     &        (  (DTDY(IJ, K-1, 1)+DTDY(IJ+LN, K-1, 1))*DZ(IJ, K-1)
     &         + (DTDY(IJ, K  , 1)+DTDY(IJ+LN, K  , 1))*DZ(IJ, K)) *
     &          0.25D0 / DZM(IJ, K) / (DTDZ(IJ, K, 1) - EPS)
            ZDZDY(IJ, K) = MIN(SLPMAX, MAX(-SLPMAX, DZDY)) *
     &                     AMFTZ(IJ, K)
  340    CONTINUE
  330 CONTINUE

      DO 165 N = 1, NTDIM
         DO 160 K = 1, NZDIM
            DO 150 IJ = 1, NXYDIM
               DTDX(IJ, K, N) = 0.D0
               DTDY(IJ, K, N) = 0.D0
               DTDZ(IJ, K, N) = 0.D0
  150      CONTINUE
  160   CONTINUE
  165 CONTINUE

      DO 420 N = 1, NTDIM
         DO 410 K = KSTR, KEND
            DO 400 IJ = IJTSTR, IJTEND+NXDIM
            DTDX(IJ, K, N) = (TX(IJ, K, N) - TX(IJ+LW, K, N)) * RX *
     &                       2.D0 / (HXT(IJ) + HXT(IJ+LW)) *
     &                       AMSKT(IJ, K) * AMSKT(IJ+LW, K)
            DTDY(IJ, K, N) = (TX(IJ, K, N) - TX(IJ+LS, K, N)) *
     &                       RYM(IJ+LS) *
     &                       2.D0 / (HYT(IJ) + HYT(IJ+LS)) *
     &                       AMSKT(IJ, K) * AMSKT(IJ+LS, K)
  400       CONTINUE
  410    CONTINUE
  420 CONTINUE

      DO 450 N = 1, NTDIM
         DO 440 K = KSTR+1, KEND
            DO 430 IJ = IJTSTR-NXDIM, IJTEND+NXDIM
               DTDZ(IJ, K, N) = (TX(IJ, K-1, N) - TX(IJ, K, N)) 
     &                          / DZM(IJ, K) * AMFTZ(IJ, K)
  430       CONTINUE
  440    CONTINUE
  450 CONTINUE

      DO 520 N = 1, NTDIM
         DO 510 K = KSTR, KEND
            DO 500 IJ = IJTSTR, IJTEND+NXDIM
               XDTDZ(IJ, K, N) =
     &         (  DTDZ(IJ, K, N)   + DTDZ(IJ+LW, K, N)
     &          + DTDZ(IJ, K+1, N) + DTDZ(IJ+LW, K+1, N)) *  0.25D0
               YDTDZ(IJ, K, N) =
     &         (  (DTDZ(IJ, K, N)    + DTDZ(IJ, K+1, N)) *
     &            DY(IJ)
     &          + (DTDZ(IJ+LS, K, N) + DTDZ(IJ+LS, K+1, N)) *
     &            DY(IJ+LS)) *
     &         RYM(IJ+LS) * 0.25D0
  500       CONTINUE
  510    CONTINUE
  520 CONTINUE

      DO 550 N = 1, NTDIM
         DO 540 K = KSTR+1, KEND
            DO 530 IJ = IJTSTR, IJTEND
               ZDTDX(IJ, K, N) = 
     &              (  (DTDX(IJ, K-1, N) + DTDX(IJ+LE, K-1, N)) *
     &                 DZ(IJ, K-1)
     &               + (DTDX(IJ, K, N) + DTDX(IJ+LE, K, N)) *
     &                 DZ(IJ, K))
     &              / DZM(IJ, K) * 0.25D0 * AMFTZ(IJ, K)
               ZDTDY(IJ, K, N) =
     &              (  (DTDY(IJ, K-1, N) + DTDY(IJ+LN, K-1, N)) *
     &                 DZ(IJ, K-1)
     &               + (DTDY(IJ, K, N) + DTDY(IJ+LN, K, N)) *
     &                 DZ(IJ, K))
     &              / DZM(IJ, K) * 0.25D0 * AMFTZ(IJ, K)
  530       CONTINUE
  540    CONTINUE
  550 CONTINUE

      RETURN
      END
#ifdef OPT_BBL
c *********************************************************************

      SUBROUTINE FLXTRB(
     O                ADT,  DIFFZ,
     I                 TX,
     I                 TY,
     I                 UY,     VY,      W,
     I                AHV)
      use zocdim
      use zocgrd
      use zocmsk
      IMPLICIT NONE
      
      REAL*8     ADT(NXYDIM, NZDIM, NTDIM)
      REAL*8   DIFFZ(NXYDIM, NZDIM)
      REAL*8      TX(NXYDIM, NZDIM, NTDIM),     TY(NXYDIM, NZDIM, NTDIM)
      REAL*8      UY(NXYDIM, NZDIM),     VY(NXYDIM, NZDIM)
      REAL*8       W(NXYDIM, NZDIM),    AHV(NXYDIM, NZDIM)

      REAL*8     FTX(NXYDIM, NZDIM, NTDIM)
      REAL*8     FTY(NXYDIM, NZDIM, NTDIM)
      REAL*8     FTZ(NXYDIM, NZDIM, NTDIM)
      COMMON /WORKFT/ FTX, FTY, FTZ

      REAL*8   TADVX(NXYDIM, NTDIM),  TADVY(NXYDIM, NTDIM)
      COMMON /UTADV/ TADVX, TADVY

      REAL*8       U,      V
      INTEGER     IJ,      K,   IJLS,   IJLW,      N
      INTEGER  IFPAR,  JFPAR

      LOGICAL OFIRST
      SAVE OFIRST

      REAL*8  AHHBBL
      SAVE AHHBBL
      NAMELIST /NMBBDH/ AHHBBL

      DATA OFIRST / .TRUE. /
      DATA AHHBBL / 0.D0 /

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF

      IF (OFIRST) THEN
         OFIRST = .FALSE.
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMBBDH, END=98)
 98      CONTINUE
         WRITE(JFPAR, NMBBDH)
      END IF

      DO N = 1, NTDIM
         DO IJ = IJTSTR, IJTEND
            K = NBOT(IJ)
            FTZ(IJ, KEND, N) = FTZ(IJ, K, N)
         END DO
      END DO

      DO N = 1, NTDIM
         DO IJ = IJTSTR, IJTEND+NXDIM
            IJLS = IJ + LS
            V = 0.5D0 *
     &          (  VY(IJ+LSW, KEND) * DZV(IJ+LSW, KEND) * HXU(IJ+LSW)
     &           + VY(IJLS  , KEND) * DZV(IJLS  , KEND) * HXU(IJLS  ))
            FTY(IJ, KEND, N) =
     &        AHHBBL * (TX(IJ, KEND, N) - TX(IJ+LS, KEND, N)) *
     &        RYM(IJ) * 
     &        (HXU(IJ+LS) + HXU(IJ+LSW)) / (HYT(IJ) + HYT(IJ+LS)) *
     &        AMFTY(IJ, KEND)
     &      - V * TADVY(IJ, N)
         END DO

         DO IJ = IJTSTR, IJTEND+1
            IJLW = IJ + LW
            U   = 0.5D0 *
     &            (  UY(IJLW  , KEND) * DZV(IJLW  , KEND) * HYU(IJLW)
     &             + UY(IJ+LSW, KEND) * DZV(IJ+LSW, KEND) * HYU(IJ+LSW))
            FTX(IJ, KEND, N) = 
     &        AHHBBL * (TX(IJ, KEND, N) - TX(IJ+LW, KEND, N)) * RX *
     &        (HYU(IJ+LW) + HYU(IJ+LSW)) / (HXT(IJ) + HXT(IJ+LW)) *
     &        AMFTX(IJ, KEND)
     &      - U * TADVX(IJ, N)
         END DO

         DO IJ = IJTSTR, IJTEND
            ADT(IJ, KEND, N) =
     &      (  (  (  FTX(IJ+LE, KEND, N) - FTX(IJ, KEND, N)) * RX
     &          + (  FTY(IJ+LN, KEND, N) - FTY(IJ, KEND, N)) * RY(IJ)) *
     &         RXT(IJ) * RYT(IJ)
     &       + FTZ(IJ, KEND, N)) / DZ(IJ, KEND)
         END DO
      END DO

      RETURN
      END
#endif
c *********************************************************************
      SUBROUTINE CHKFTX
      use zocdim
      implicit none

      REAL*8     FTX(NXYDIM, NZDIM, NTDIM)
      REAL*8     FTY(NXYDIM, NZDIM, NTDIM)
      REAL*8     FTZ(NXYDIM, NZDIM, NTDIM)
      COMMON /WORKFT/ FTX, FTY, FTZ

      IF (OINIT .OR. OFINAL) THEN
         RETURN
      END IF

      CALL CHEKIN(   FTX,  'FTX',
     &                NX,     NY,     NZ, NXYZDM, 'OCN')
      CALL CHEKIN(   FTY,  'FTY',
     &                NX,     NY,     NZ, NXYZDM, 'OCN')
      CALL CHEKIN(   FTZ,  'FTZ',
     &                NX,     NY,     NZ, NXYZDM, 'OCN')
      CALL CHEKIN(   FTX(1, 1, 2),  'FSX',
     &                NX,     NY,     NZ, NXYZDM, 'OCN')
      CALL CHEKIN(   FTY(1, 1, 2),  'FSY',
     &                NX,     NY,     NZ, NXYZDM, 'OCN')
      CALL CHEKIN(   FTZ(1, 1, 2),  'FSZ',
     &                NX,     NY,     NZ, NXYZDM, 'OCN')

      RETURN
      END 
