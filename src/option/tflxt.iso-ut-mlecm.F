      SUBROUTINE FLXTRC(
     O                ADT,  DIFFZ,
     I                 TX,     HX,
     I                 TY,     HZ,
     I                 UY,     VY,
     I                  W,    AHV)
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
c     '08.08.??  Y.Komuro: Fox-Kemper et al. (2008) parameterization
c     '09.07.30  Y.Komuro: bug fix (sign of mle fluxes)
c     '10.11.29  Y.Komuro: Coarse-resolution modification applied
c                          following Fox-Kemper et al. (2011, OM).
c
c ---------------------------------------------------------------------

#include "zocdim.F"
#include "zocgrd.F"
#include "zocmsk.F"

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
      COMMON /WORK/ FTX, FTY, FTZ, WZC, RZM

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

      REAL*8  PSIGMX(NXYDIM, NZDIM), PSIGMY(NXYDIM, NZDIM)

c     for mixed layer eddy parameterization
      REAL*8   XPSIY(NXYDIM, NZDIM),  YPSIX(NXYDIM, NZDIM)
      REAL*8   ZPSIX(NXYDIM, NZDIM),  ZPSIY(NXYDIM, NZDIM)
      REAL*8    IGSY(NXYDIM, NZDIM),   IGSX(NXYDIM, NZDIM)

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
     O             XPSIY,  YPSIX,
     O             ZPSIX,  ZPSIY,
     I                TY,     TX,     HZ,
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
     &          - ( ( AHI + AHG ) *
     &             ZDZDX(IJ, K) - ZPSIY(IJ, K) ) * ZDTDX(IJ, K, N)
     &          - ( ( AHI + AHG ) *
     &             ZDZDY(IJ, K) + ZPSIX(IJ, K) ) * ZDTDY(IJ, K, N)
     &          - WZC(IJ, K) * TADV ) * AMFTZ(IJ, K)
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
     &           (  ( AHH + AHI ) * RYM(IJLS) 
     &            / ( HYT(IJ) + HYT(IJLS) ) *
     &              ( TX(IJ, K, N) - TX(IJLS, K, N) ) * 2.D0
     &            - ( ( AHI - AHG ) 
     &              * YDZDY(IJ, K) - YPSIX(IJ, K) ) * YDTDZ(IJ, K, N) )
     &            * ( HXU(IJLS) + HXU(IJ+LSW) ) * 0.5D0 * AMFTY(IJ, K)
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
     &           (  ( AHH + AHI ) * RX 
     &            / ( HXT(IJ) + HXT(IJLW) ) 
     &            * ( TX(IJ, K, N) - TX(IJLW, K, N) ) * 2.D0
     &            - ( ( AHI - AHG )  
     &            * XDZDX(IJ, K) + XPSIY(IJ, K) ) * XDTDZ(IJ, K, N) ) 
     &            * ( HYU(IJLW) + HYU(IJ+LSW) ) * 0.5D0 * AMFTX(IJ, K)
     &         - U * TADVX(IJ, N)

  420       CONTINUE
  430    CONTINUE
         
  440 CONTINUE

      DO K = KSTR, KEND
         DO IJ=1, NXYDIM
            IGSX(IJ, K) = ( AHI- AHG )  
     &            * XDZDX(IJ, K)
            IGSY(IJ, K) = ( AHI - AHG ) 
     &              * YDZDY(IJ, K)
         END DO
      END DO
      CALL CHEKIN(IGSX, 'IGSX', NX, NY, NZ, NXYZDM, 'OCN')
      CALL CHEKIN(IGSY, 'IGSY', NX, NY, NZ, NXYZDM, 'OCN')

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

      DO K = 1, NZDIM
         DO IJ = 1, NXYDIM
            PSIGMX(IJ, K) = 0.0D0
            PSIGMY(IJ, K) = 0.0D0
         END DO
      END DO

      DO K = KSTR, KEND
         DO IJ = IJTSTR, IJTEND+NXDIM
            PSIGMX(IJ, K) = AHG * YDZDY(IJ, K)
            PSIGMY(IJ, K) = - AHG * XDZDX(IJ, K)
         END DO
      END DO

      CALL CHEKIN(PSIGMX, 'PSIGMX', NX, NY, NZ, NXYZDM, 'OCN')
      CALL CHEKIN(PSIGMY, 'PSIGMY', NX, NY, NZ, NXYZDM, 'OCN')

      RETURN
      END

c *********************************************************************

      SUBROUTINE DNSGRD(
     O              XDZDX,  YDZDY,  ZDZDX,  ZDZDY,
     O              XDTDZ,  YDTDZ,  ZDTDX,  ZDTDY,
     O              XPSIY,  YPSIX,
     O              ZPSIX,  ZPSIY,
     I                 TY,     TX,     HZ,
     W               DTDX,   DTDY,  DTFDZ)
      IMPLICIT NONE

#include "zocdim.F"
#include "zocgrd.F"
#include "zocmsk.F"
#include "zocphy.F"

      REAL*8   XDZDX(NXYDIM, NZDIM),  YDZDY(NXYDIM, NZDIM)
      REAL*8   ZDZDX(NXYDIM, NZDIM),  ZDZDY(NXYDIM, NZDIM)
      REAL*8   XDTDZ(NXYDIM, NZDIM, NTDIM),  YDTDZ(NXYDIM, NZDIM, NTDIM)
      REAL*8   ZDTDX(NXYDIM, NZDIM, NTDIM),  ZDTDY(NXYDIM, NZDIM, NTDIM)
      REAL*8      TY(NXYDIM, NZDIM, NTDIM)
      REAL*8      TX(NXYDIM, NZDIM, NTDIM)
      REAL*8    DTDX(NXYDIM, NZDIM, NTDIM),   DTDY(NXYDIM, NZDIM, NTDIM)
      REAL*8    DTFDZ(NXYDIM, NZDIM, NTDIM)

      REAL*8   XPSIY(NXYDIM, NZDIM),  YPSIX(NXYDIM, NZDIM)
      REAL*8   ZPSIX(NXYDIM, NZDIM),  ZPSIY(NXYDIM, NZDIM)
      REAL*8      HZ(NXYDIM)

      REAL*8 C0(NZDIM), C1(NZDIM), C2(NZDIM)
      REAL*8 C3(NZDIM), C4(NZDIM), C5(NZDIM), C6(NZDIM)
      REAL*8 D0(NZDIM), D1(NZDIM), D2(NZDIM), D3(NZDIM), D4(NZDIM)
      REAL*8 D5(NZDIM), D6(NZDIM), D7(NZDIM), D8(NZDIM), D9(NZDIM)
      REAL*8 EPS
      LOGICAL OFIRST
      SAVE C0, C1, C2, C3, C4, C5, C6,
     &     D0, D1, D2, D3, D4, D5, D6, D7, D8, D9
      SAVE EPS, OFIRST

      REAL*8 CXPSY(NXYDIM), CYPSX(NXYDIM) 
      REAL*8 CZPSX(NXYDIM), CZPSY(NXYDIM) 
      INTEGER KZMIN
      SAVE CXPSY, CYPSX, CZPSX, CZPSY, KZMIN

      REAL*8       R(NXYDIM, NZDIM)
      REAL*8    HMLD(NXYDIM), HMLD1(NXYDIM)
      REAL*8  RMAVEZ(NXYDIM), RMAV1(NXYDIM)
      REAL*8   DZSIG(NXYDIM, NZDIM), DZMSIG(NXYDIM, NZDIM)
      REAL*8      ZT(NXYDIM, NZDIM),    ZTM(NXYDIM, NZDIM)
      REAL*8  RSIGTH(NXYDIM, NZDIM)
      REAL*8     NBV(NXYDIM),            LF(NXYDIM)
      REAL*8  XPSIY1(NXYDIM, NZDIM), YPSIX1(NXYDIM, NZDIM)
      REAL*8  ZPSIX1(NXYDIM, NZDIM), ZPSIY1(NXYDIM, NZDIM)
      REAL*8   ZMLD0(NXYDIM),         ZHMLD(NXYDIM, NZDIM)
      REAL*8   HMLDX(NXYDIM),         HMLDY(NXYDIM)
      INTEGER   KMLD(NXYDIM)

      REAL*8  RMAVDX(NXYDIM), RMAVDY(NXYDIM)
c      REAL*8    MUZX(NXYDIM, NZDIM),   MUZY(NXYDIM, NZDIM)

      REAL*8      P1,     P2
      REAL*8      TL,     SL
      REAL*8      RL,    RLW,    RLS,    RLU
      REAL*8    DZDX,   DZDY
      INTEGER     IJ,      K,      N
      INTEGER  IFPAR,  JFPAR

      REAL*8    MUZH,  IN2DZ, N2MIN, N2L, HMLDT
      REAL*8  RSIGDF, RSIGBT,  DHMLD
      REAL*8      PI,  OMEGA, CORMIN

      REAL*8  SLPMAX
      SAVE SLPMAX
      NAMELIST /NMSLPM/ SLPMAX

      DATA OFIRST / .TRUE. /
      DATA SLPMAX / 1.D-2 /

      REAL*8 CM, CE, FMINLT, LFMIN, TAUMLE, VSCL, DRSIG
      INTEGER MZ, NFLTDM, NFLTPS, NFLTRM, MZMIN
      LOGICAL OFLTDM, OFLTPS, OFLTRM, OCOAMP
      SAVE CM, CE, FMINLT, MZ, OFLTDM, NFLTDM, LFMIN, TAUMLE, VSCL,
     &     OFLTPS, NFLTPS, MZMIN, DRSIG, OCOAMP, OFLTRM, NFLTRM
      NAMELIST /NMMLEP/ CM, CE, FMINLT, MZ, OFLTDM, NFLTDM,
     &                  LFMIN, TAUMLE, VSCL, OFLTPS, NFLTPS, MZMIN,
     &                  DRSIG, OCOAMP, OFLTRM, NFLTRM
      DATA CM, CE, FMINLT / 8.0D0, 0.06D0, 10.0D0 /
      DATA LFMIN, TAUMLE, VSCL, DRSIG / 1.0D5, 10.0D0, 50.0D0, 0.1D0 /
      DATA MZ, MZMIN / NZ, 1 /
      DATA OFLTDM, OFLTPS, OFLTRM, OCOAMP
     &               / .FALSE., .FALSE., .FALSE., .TRUE. /
      DATA NFLTDM, NFLTPS, NFLTRM / 0, 0, 0 /

      IF (OFIRST) THEN
         OFIRST = .FALSE.
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMSLPM, END=99)
 99      CONTINUE
         WRITE(JFPAR, NMSLPM)
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMMLEP, END=98)
 98      WRITE(JFPAR, NMMLEP)

         CALL SECOEF(
     O            C0(KSTR), C1(KSTR), C2(KSTR), C3(KSTR),
     O            C4(KSTR), C5(KSTR), C6(KSTR),
     O            D0(KSTR), D1(KSTR), D2(KSTR), D3(KSTR),  D4(KSTR),
     O            D5(KSTR), D6(KSTR), D7(KSTR), D8(KSTR),  D9(KSTR))
         EPS = 1.D-20

c     === Mixed layer eddy parameterization
c         from Fox-Kemper and Ferrari(2008) ===
c     *** Coarse-resolution modification is applied.
c         The last row in CXPSY/CYPSX corresponds to ds,
c         and is cancelled out by the row just before it. *** 
         PI = ATAN( 1.D0 )*4.D0
         OMEGA = 2.D0 * PI / 86400.D0
         CORMIN = 2.D0 * OMEGA * SIN( PI*ABS(FMINLT)/180.D0 )
         KZMIN = MZMIN + KSTR - 1

         DO IJ = NXDIM+2, NXYDIM
            IF (OCOAMP) THEN
c               DSCM = max(DX * 0.5D0 * (HXT(IJ) + HXT(IJ+LW)),
c     &                    DYM(IJ+LS) * 0.5D0 * (HYT(IJ) + HYT(IJ+LS)))
               CXPSY(IJ) = CE
     &               / sqrt(
     &                 max( abs(0.5D0*(COR(IJ+LSW) + COR(IJ+LW))),
     &                      CORMIN )**2.0D0
     &                 + 1.0D0 / (8.64D4 * TAUMLE)**2.0D0 )
c     &               * DSCM
     &               * DX * 0.5D0 * (HXT(IJ) + HXT(IJ+LW))
               CYPSX(IJ) = CE
     &               / sqrt(
     &                 max( abs(0.5D0*(COR(IJ+LSW) + COR(IJ+LS))),
     &                      CORMIN )**2.0D0
     &                 + 1.0D0 / (8.64D4 * TAUMLE)**2.0D0 )
c     &               * DSCM
     &               * DYM(IJ+LS) * 0.5D0 * (HYT(IJ) + HYT(IJ+LS))        
               CZPSY(IJ) = CE
     &               / sqrt(
     &                 max( abs(0.25D0*
     &                      ( COR(IJ    ) + COR(IJ+LW )
     &                      + COR(IJ+LS ) + COR(IJ+LSW))), CORMIN )
     &                                                      **2.0D0
     &                 + 1.0D0 / (8.64D4 * TAUMLE)**2.0D0 )
     &                 * DX * HXT(IJ)
               CZPSX(IJ) = CE
     &               / sqrt(
     &                 max( abs(0.25D0*
     &                      ( COR(IJ    ) + COR(IJ+LW )
     &                      + COR(IJ+LS ) + COR(IJ+LSW))), CORMIN )
     &                                                      **2.0D0
     &                 + 1.0D0 / (8.64D4 * TAUMLE)**2.0D0 )
     &                 * DY(IJ) * HYT(IJ)        
            ELSE
               CXPSY(IJ) = CE
     &               / sqrt(
     &                 max( abs(0.5D0*(COR(IJ+LSW) + COR(IJ+LW))),
     &                      CORMIN )**2.0D0 )
               CYPSX(IJ) = CE
     &               / sqrt(
     &                 max( abs(0.5D0*(COR(IJ+LSW) + COR(IJ+LS))),
     &                      CORMIN )**2.0D0 )
               CZPSY(IJ) = CE
     &               / sqrt(
     &                 max( abs(0.25D0*
     &                      ( COR(IJ    ) + COR(IJ+LW )
     &                      + COR(IJ+LS ) + COR(IJ+LSW))), CORMIN )
     &                                                      **2.0D0 )
               CZPSX(IJ) = CE
     &               / sqrt(
     &                 max( abs(0.25D0*
     &                      ( COR(IJ    ) + COR(IJ+LW )
     &                      + COR(IJ+LS ) + COR(IJ+LSW))), CORMIN )
     &                                                      **2.0D0 )
            END IF
         END DO
      END IF

      DO 115 N = 1, NTDIM
         DO 110 K = 1, NZDIM
            DO 100 IJ = 1, NXYDIM
               DTDX(IJ, K, N) = 0.D0
               DTDY(IJ, K, N) = 0.D0
               DTFDZ(IJ, K, N) = 0.D0
  100       CONTINUE
  110    CONTINUE
  115 CONTINUE

      DO K = KSTR, KEND
         DO IJ = 1, NXYDIM
            TL = TY(IJ, K, 1)
            SL = TY(IJ, K, 2)
            P1 = C0(K)
     &         + (C1(K) + (C2(K) + C3(K) * TL) * TL) * TL
     &         + (C4(K) + C5(K) * TL + C6(K) * SL) * SL
            P2 = D0(K)
     &         + (D1(K) + (D2(K) + (D3(K) + D4(K) * TL) * TL) * TL) * TL
     &         + (D5(K) + (D6(K) + D7(K) * TL * TL) * TL
     &                  + (D8(K) + D9(K) * TL * TL) * SQRT(SL)) * SL
            R(IJ, K) = P1 / P2
            P1 = C0(KSTR)
     &         + (C1(KSTR) + (C2(KSTR) + C3(KSTR) * TL) * TL) * TL
     &         + (C4(KSTR) + C5(KSTR) * TL + C6(KSTR) * SL) * SL
            P2 = D0(KSTR)
     &         + (D1(KSTR) + (D2(KSTR) +
     &              (D3(KSTR) + D4(KSTR) * TL) * TL) * TL) * TL
     &         + (D5(KSTR) + (D6(KSTR) + D7(KSTR) * TL * TL) * TL
     &                + (D8(KSTR) + D9(KSTR) * TL * TL) * SQRT(SL)) * SL
            RSIGTH(IJ, K) = P1 / P2
         END DO
      END DO

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

            DTFDZ(IJ, K, 1) = MIN((RLU - RL) / DZM(IJ, K), 0.D0) *
     &                        AMFTZ(IJ, K)
  240    CONTINUE
  230 CONTINUE

c     === Mixed layer eddy parameterization
c         from Fox-Kemper and Ferrari(2008) ===

      DO K = KSTR, KSTR+KZ-1
         DO IJ = 1, NXYDIM
c            DZSIG (IJ, K) = (HZ(IJ) + ZBOT) * DS(K)
c            DZMSIG(IJ, K) = (HZ(IJ) + ZBOT) * DSM(K)
            DZSIG (IJ, K) = ZBOT * DS(K)
            DZMSIG(IJ, K) = ZBOT * DSM(K)
         END DO
      END DO
      DO K = KSTR+KZ, KEND
         DO IJ = 1, NXYDIM
            DZSIG (IJ, K) = DZ(IJ, K)
            DZMSIG(IJ, K) = DZM(IJ, K)
         END DO
      END DO
      DO K = KEND+1, NZDIM
         DO IJ = 1, NXYDIM
            DZSIG (IJ, K) = DZSIG(IJ, KEND)
            DZMSIG(IJ, K) = DZMSIG(IJ, KEND)
         END DO
      END DO

      DO IJ = 1, NXYDIM
         ZT(IJ, KSTR) = 0.D0
        ZTM(IJ, KSTR) = 0.5D0 * DZMSIG(IJ, KSTR)
      END DO
      DO K = KSTR, KEND
         DO IJ = 1, NXYDIM
            ZT(IJ, K+1) = ZT(IJ, K) + DZSIG(IJ, K)
           ZTM(IJ, K+1) = ZTM(IJ, K) + DZMSIG(IJ, K+1)
         END DO
      END DO

      DO K = KSTR, KEND
         DO IJ = 1, NXYDIM
            ZHMLD(IJ, K) = 0.0D0
         END DO
      END DO

      DO IJ = 1, NXYDIM
         RMAVDX(IJ) = 0.0D0
         RMAVDY(IJ) = 0.0D0
         HMLDX(IJ) = 0.0D0
         HMLDY(IJ) = 0.0D0
      END DO

c     calculating mixed layer depth HMLD and MLD-averaged buoyancy RMAVEZ
c      DO IJ = IJTSTR-NXDIM, IJTEND+NXDIM
cc      DO IJ = 1, NXYDIM
c         IN2DZ = (-1.0D-3) * gravit / rhoo
c     &           * DTFDZ(IJ, KSTR+1, 1) * ZTM(IJ, KSTR+1)
c         N2MIN = (-1.0D-3) * gravit / rhoo
c     &           * DTFDZ(IJ, KSTR+1, 1)
c         HMLD(IJ) = ZTM(IJ, KSTR+1)
c         RMAVEZ(IJ) = R(IJ, KSTR) * DZSIG(IJ, KSTR)
c     &                + 0.5D0 * R(IJ, KSTR+1) * DZSIG(IJ, KSTR+1)
c         DO K = KSTR+2, MIN(NBOT(IJ), MZ)
c            N2L = (-1.0D-3) * gravit / rhoo * DTFDZ(IJ, K, 1)
c            IN2DZ = IN2DZ + N2L * DZMSIG(IJ, K)
c            N2MIN = MIN(N2MIN, N2L)
c            IF ( ( (N2L - N2MIN) * ZTM(IJ, K) )
c     &           .GE. (CM * IN2DZ) ) THEN
cc            IF ( ( N2L * ZTM(IJ, K) )
cc     &           .GE. (CM * IN2DZ) ) THEN
c               EXIT
c            END IF
c            HMLD(IJ) = ZTM(IJ, K)
c            RMAVEZ(IJ) = RMAVEZ(IJ) + 0.5D0 * 
c     &           (  R(IJ, K-1) * DZSIG(IJ, K-1)
c     &                + R(IJ, K) * DZSIG(IJ, K) )
c         END DO
c         RMAVEZ(IJ) = RMAVEZ(IJ) / HMLD(IJ)
c      END DO

c      DO IJ = IJTSTR-NXDIM, IJTEND+NXDIM
      DO IJ = 1, NXYDIM
         HMLD(IJ) = ZTM(IJ, KSTR)
         ZHMLD(IJ, KSTR) = ZTM(IJ, KSTR)
         KMLD(IJ) = KSTR
         RMAVEZ(IJ) = RSIGTH(IJ, KSTR) * ZTM(IJ, KSTR)
         RSIGBT = RSIGTH(IJ, KSTR) + DRSIG
         RSIGDF = 0.D0
         DO K = KSTR+1, MIN(NBOT(IJ), MZ)
            IF ( RSIGTH(IJ, K) .GE. RSIGBT ) THEN
               DHMLD = DZMSIG(IJ, K) * (RSIGBT - RSIGTH(IJ, K-1))
     &                 / (RSIGTH(IJ, K) - RSIGTH(IJ, K-1))
               HMLD(IJ) = HMLD(IJ) + DHMLD
               IF (DHMLD .LE. 0.5D0*DZSIG(IJ, K-1)) THEN
                  ZHMLD(IJ, K-1) = ZHMLD(IJ, K-1) + DHMLD
               ELSE
                  ZHMLD(IJ, K-1) = DZSIG(IJ, K-1)
                  ZHMLD(IJ, K) = DHMLD - 0.5D0 * DZSIG(IJ, K-1)
               END IF
               RMAVEZ(IJ) = RMAVEZ(IJ) +
     &              0.5D0 * (RSIGTH(IJ, K-1)+RSIGBT) * DHMLD
               RSIGDF = DRSIG
               EXIT
            END IF
            HMLD(IJ) = HMLD(IJ) + DZMSIG(IJ, K)
            ZHMLD(IJ, K-1) = DZSIG(IJ, K-1)
            ZHMLD(IJ, K) = 0.5D0 * DZSIG(IJ, K)
            KMLD(IJ) = K
            RMAVEZ(IJ) = RMAVEZ(IJ) + 0.5D0 *
     &                 ( RSIGTH(IJ, K-1) * DZSIG(IJ, K-1)
     &                 + RSIGTH(IJ, K  ) * DZSIG(IJ, K  ) )
            RSIGDF = RSIGTH(IJ, K) - RSIGTH(IJ, KSTR)
         END DO
         RMAVEZ(IJ) = RMAVEZ(IJ) / HMLD(IJ)
c         RMAVEZ(IJ) = RSIGTH(IJ, KSTR)
         NBV(IJ) = sqrt(max(1.0D-3*gravit/rhoo*RSIGDF/HMLD(IJ),0.0D0))
      END DO

      DO IJ = IJTSTR-NXDIM, IJTEND+NXDIM
         IF (OCOAMP) THEN
            LF(IJ) = max(4.0D0*NBV(IJ)*HMLD(IJ)/
     &                max( abs( COR(IJ   ) + COR(IJ+LW )
     &                        + COR(IJ+LS) + COR(IJ+LSW) ), EPS), LFMIN)
         ELSE
            LF(IJ) = 1.0D0
         END IF
      END DO

      IF (OFLTDM) THEN
#ifdef OPT_PARALLEL
         CALL SHIFT1(
     M                 HMLD,
     D                NXDIM,  NYDIM,      1)
#else
         CALL STBCSH(
     M                 HMLD)
#endif
         DO N = 1, NFLTDM
c            DO IJ = IJTSTR-NXDIM, IJTEND+NXDIM
            DO IJ = IJSTR, IJEND
               IF (AMSKT(IJ, KSTR) .EQ. 1.0D0) THEN
                  HMLD1(IJ) = ( 4.0D0 * HMLD(IJ) * AMSKT(IJ, KSTR)
     &                       + HMLD(IJ+LS) * AMSKT(IJ+LS, KSTR)
     &                       + HMLD(IJ+LN) * AMSKT(IJ+LN, KSTR)
     &                       + HMLD(IJ+LW) * AMSKT(IJ+LW, KSTR)
     &                       + HMLD(IJ+LE) * AMSKT(IJ+LE, KSTR) )
     &                    / ( 4.0D0 * AMSKT(IJ, KSTR)
     &                       + AMSKT(IJ+LS, KSTR) + AMSKT(IJ+LN, KSTR)
     &                       + AMSKT(IJ+LW, KSTR) + AMSKT(IJ+LE, KSTR) )
               ELSE
                  HMLD1(IJ) = HMLD(IJ)
               END IF
            END DO
c            DO IJ = IJTSTR-NXDIM, IJTEND+NXDIM
            DO IJ = IJSTR, IJEND
               HMLD(IJ) = HMLD1(IJ)
            END DO
#ifdef OPT_PARALLEL
            CALL SHIFT1(
     M                    HMLD,
     D                   NXDIM,  NYDIM,      1)
#else
            CALL STBCSH(
     M                    HMLD)
#endif
         END DO
      END IF

      IF (OFLTRM) THEN
#ifdef OPT_PARALLEL
         CALL SHIFT1(
     M               RMAVEZ,
     D                NXDIM,  NYDIM,      1)
#else
         CALL STBCSH(
     M               RMAVEZ)
#endif
         DO N = 1, NFLTRM
c            DO IJ = IJTSTR-NXDIM, IJTEND+NXDIM
            DO IJ = IJSTR, IJEND
               IF (AMSKT(IJ, KSTR) .EQ. 1.0D0) THEN
                  RMAV1(IJ) = ( 4.0D0 * RMAVEZ(IJ) * AMSKT(IJ, KSTR)
     &                       + RMAVEZ(IJ+LS) * AMSKT(IJ+LS, KSTR)
     &                       + RMAVEZ(IJ+LN) * AMSKT(IJ+LN, KSTR)
     &                       + RMAVEZ(IJ+LW) * AMSKT(IJ+LW, KSTR)
     &                       + RMAVEZ(IJ+LE) * AMSKT(IJ+LE, KSTR) )
     &                    / ( 4.0D0 * AMSKT(IJ, KSTR)
     &                       + AMSKT(IJ+LS, KSTR) + AMSKT(IJ+LN, KSTR)
     &                       + AMSKT(IJ+LW, KSTR) + AMSKT(IJ+LE, KSTR) )
               ELSE
                  RMAV1(IJ) = RMAVEZ(IJ)
               END IF
            END DO
c            DO IJ = IJTSTR-NXDIM, IJTEND+NXDIM
            DO IJ = IJSTR, IJEND
               RMAVEZ(IJ) = RMAV1(IJ)
            END DO
#ifdef OPT_PARALLEL
            CALL SHIFT1(
     M                  RMAVEZ,
     D                   NXDIM,  NYDIM,      1)
#else
            CALL STBCSH(
     M                  RMAVEZ)
#endif
         END DO
      END IF

      DO IJ = 1, NXYDIM
         IF (KMLD(IJ) .LT. KZMIN) THEN
            HMLD(IJ) = 0.0D0
         END IF
         ZMLD0(IJ) = KMLD(IJ) - KSTR + 1
      END DO
        
      CALL CHEKIN(HMLD, 'HMLD', NX, NY, 1, NXYDIM, 'SFC')
      CALL CHEKIN(ZMLD0, 'KMLD', NX, NY, 1, NXYDIM, 'SFC')
      CALL CHEKIN(RMAVEZ, 'RMAVEZ', NX, NY, 1, NXYDIM, 'SFC')
      CALL CHEKIN(R, 'R', NX, NY, NZ, NXYZDM, 'OCN')
      CALL CHEKIN(LF, 'LF', NX, NY, 1, NXYDIM, 'SFC')

c      DO IJ = IJTSTR, IJTEND+NXDIM
c         RMAVDX(IJ) = (-1.0D-3) * gravit / rhoo 
c     &               * RX * 2.D0 / (HXT(IJ) + HXT(IJ+LW))
c     &               * (RMAVEZ(IJ) - RMAVEZ(IJ+LW))
c     &               * AMSKT(IJ, KSTR) * AMSKT(IJ+LW, KSTR)
c         RMAVDY(IJ) = (-1.0D-3) * gravit / rhoo 
c     &               * RYM(IJ+LS) * 2.D0 / (HYT(IJ) + HYT(IJ+LS))
c     &               * (RMAVEZ(IJ) - RMAVEZ(IJ+LS))
c     &               * AMSKT(IJ, KSTR) * AMSKT(IJ+LS, KSTR)
c      END DO

      DO K = KSTR, MZ
         DO IJ = IJTSTR, IJTEND+NXDIM
c            DHMLD = MIN(ZHMLD(IJ, K), ZHMLD(IJ+LW, K))
            DHMLD = MAX(ZHMLD(IJ, K), ZHMLD(IJ+LW, K))
            RMAVDX(IJ) = RMAVDX(IJ) + 
     &               DHMLD * (-1.0D-3) * gravit / rhoo 
     &               * RX * 2.D0 / (HXT(IJ) + HXT(IJ+LW))
     &               * (RSIGTH(IJ, K) - RSIGTH(IJ+LW, K))
            HMLDX(IJ) = HMLDX(IJ) + DHMLD
c            DHMLD = MIN(ZHMLD(IJ, K), ZHMLD(IJ+LS, K))
            DHMLD = MAX(ZHMLD(IJ, K), ZHMLD(IJ+LS, K))
            RMAVDY(IJ) = RMAVDY(IJ) +
     &               DHMLD * (-1.0D-3) * gravit / rhoo 
     &               * RYM(IJ+LS) * 2.D0 / (HYT(IJ) + HYT(IJ+LS))
     &               * (RSIGTH(IJ, K) - RSIGTH(IJ+LS, K))
            HMLDY(IJ) = HMLDY(IJ) + DHMLD
         END DO
      END DO

      DO IJ = IJTSTR, IJTEND+NXDIM
         RMAVDX(IJ) = RMAVDX(IJ) / MAX(HMLDX(IJ), EPS)
     &              * AMSKT(IJ, KSTR) * AMSKT(IJ+LW, KSTR)
         RMAVDY(IJ) = RMAVDY(IJ) / MAX(HMLDY(IJ), EPS)
     &              * AMSKT(IJ, KSTR) * AMSKT(IJ+LS, KSTR)
      END DO

c     calculating XPSIY and YPSIX
      DO K = KSTR, KEND
         DO IJ = IJTSTR, IJTEND+NXDIM
c         DO IJ = NXDIM+2, NXYDIM
            HMLDT = MAX(HMLD(IJ), HMLD(IJ+LW))
c            HMLDT = MIN(HMLD(IJ), HMLD(IJ+LW))
c            HMLDT = 0.5D0 * (HMLD(IJ) + HMLD(IJ+LW))
            MUZH = ( 1.0D0 - 2.0D0 *
     &               MIN(0.5D0*(ZTM(IJ, K)+ZTM(IJ+LW, K))
     &                   / MAX(HMLDT, EPS), 1.0D0)
     &               ) ** 2.0D0
c            MUZX(IJ, K) = ( 1.0D0 - MUZH )
c     &                  * ( 1.0D0 + 5.0D0 / 21.0D0 * MUZH )
            XPSIY(IJ, K) = CXPSY(IJ)
     &           * 2.0D0 / (LF(IJ) + LF(IJ+LW))
     &           * HMLDT * HMLDT
     &           * (-1.0D0) * RMAVDX(IJ)
     &           * ( 1.0D0 - MUZH ) * ( 1.0D0 + 5.0D0 / 21.0D0 * MUZH )     
     &           * AMSKT(IJ, K) * AMSKT(IJ+LW, K)
            XPSIY(IJ, K) = sign(1.0D0, XPSIY(IJ, K)) *
     &                     min(abs(XPSIY(IJ, K)), VSCL*DZSIG(IJ, K))
         END DO
      END DO
      DO K = KSTR, KEND
         DO IJ = IJTSTR, IJTEND+NXDIM
c         DO IJ = NXDIM+2, NXYDIM
            HMLDT = MAX(HMLD(IJ), HMLD(IJ+LS))
c            HMLDT = MIN(HMLD(IJ), HMLD(IJ+LS))
c            HMLDT = 0.5D0 * (HMLD(IJ) + HMLD(IJ+LS))
            MUZH = ( 1.0D0 - 2.0D0 *
     &               MIN(0.5D0*(ZTM(IJ, K)+ZTM(IJ+LS, K))
     &                   / MAX(HMLDT, EPS), 1.0D0)
     &               ) ** 2.0D0
c            MUZY(IJ, K) = ( 1.0D0 - MUZH )
c     &                  * ( 1.0D0 + 5.0D0 / 21.0D0 * MUZH )
            YPSIX(IJ, K) = CYPSX(IJ)
     &           * 2.0D0 / (LF(IJ) + LF(IJ+LS))
     &           * HMLDT * HMLDT
     &           * RMAVDY(IJ)
     &           * ( 1.0D0 - MUZH ) * ( 1.0D0 + 5.0D0 / 21.0D0 * MUZH )
     &           * AMSKT(IJ, K) * AMSKT(IJ+LS, K)
            YPSIX(IJ, K) = sign(1.0D0, YPSIX(IJ, K)) *
     &                     min(abs(YPSIX(IJ, K)), VSCL*DZSIG(IJ, K))
         END DO
      END DO

c      DO K = KSTR+1, KEND
c         DO IJ = IJTSTR, IJTEND
c            HMLDT = HMLD(IJ)
c            MUZH = ( 1.0D0 - 2.0D0 *
c     &               MIN(ZT(IJ, K) / MAX(HMLDT, EPS), 1.0D0)
c     &               ) ** 2.0D0
cc            MUZX(IJ, K) = ( 1.0D0 - MUZH )
cc     &                  * ( 1.0D0 + 5.0D0 / 21.0D0 * MUZH )
c            ZPSIX(IJ, K) = CZPSX(IJ)
c     &           / LF(IJ)
c     &           * HMLDT * HMLDT
c     &           * 0.5D0 * (RMAVDY(IJ) + RMAVDY(IJ+LN))
c     &           * ( 1.0D0 - MUZH ) * ( 1.0D0 + 5.0D0 / 21.0D0 * MUZH )
c     &           * AMFTZ(IJ, K)
c            ZPSIY(IJ, K) = CZPSY(IJ)
c     &           / LF(IJ)
c     &           * HMLDT * HMLDT
c     &           * (-0.5D0) * (RMAVDX(IJ) + RMAVDX(IJ+LE))
c     &           * ( 1.0D0 - MUZH ) * ( 1.0D0 + 5.0D0 / 21.0D0 * MUZH )     
c     &           * AMFTZ(IJ, K)
c            ZPSIX(IJ, K) = sign(1.0D0, ZPSIX(IJ, K)) *
c     &                     min(abs(ZPSIX(IJ, K)), 
c     &                         VSCL*0.5D0*(DZSIG(IJ, K-1)+DZSIG(IJ, K)))
c            ZPSIY(IJ, K) = sign(1.0D0, ZPSIY(IJ, K)) *
c     &                     min(abs(ZPSIY(IJ, K)), 
c     &                         VSCL*0.5D0*(DZSIG(IJ, K-1)+DZSIG(IJ, K)))
c         END DO
c      END DO


      IF (OFLTPS) THEN
#ifdef OPT_PARALLEL
         CALL SHIFT2(
     M                XPSIY,  YPSIX,
     D                NXDIM,  NYDIM,  NZDIM)
c         CALL SHIFT2(
c     M                ZPSIX,  ZPSIY,
c     D                NXDIM,  NYDIM,  NZDIM)
#else
         CALL STBCAV(
     M                XPSIY,  YPSIX)
c         CALL STBCAV(
c     M                ZPSIX,  ZPSIY)
#endif
         DO N = 1, NFLTPS
c            DO IJ = IJTSTR-NXDIM, IJTEND+NXDIM
            DO K = KSTR, KEND
               DO IJ = IJSTR, IJEND
                  IF (AMFTX(IJ, K) .EQ. 1.0D0) THEN
                     XPSIY1(IJ, K) =
     &                    ( 4.0D0 * XPSIY(IJ, K) * AMFTX(IJ, K)
     &                       + XPSIY(IJ+LS, K) * AMFTX(IJ+LS, K)
     &                       + XPSIY(IJ+LN, K) * AMFTX(IJ+LN, K)
     &                       + XPSIY(IJ+LW, K) * AMFTX(IJ+LW, K)
     &                       + XPSIY(IJ+LE, K) * AMFTX(IJ+LE, K) )
     &                    / ( 4.0D0 * AMFTX(IJ, K)
     &                       + AMFTX(IJ+LS, K) + AMFTX(IJ+LN, K)
     &                       + AMFTX(IJ+LW, K) + AMFTX(IJ+LE, K) )
                  ELSE
                     XPSIY1(IJ, K) = XPSIY(IJ, K)
                  END IF
                  IF (AMFTY(IJ, K) .EQ. 1.0D0) THEN
                     YPSIX1(IJ, K) =
     &                    ( 4.0D0 * YPSIX(IJ, K) * AMFTY(IJ, K)
     &                       + YPSIX(IJ+LS, K) * AMFTY(IJ+LS, K)
     &                       + YPSIX(IJ+LN, K) * AMFTY(IJ+LN, K)
     &                       + YPSIX(IJ+LW, K) * AMFTY(IJ+LW, K)
     &                       + YPSIX(IJ+LE, K) * AMFTY(IJ+LE, K) )
     &                    / ( 4.0D0 * AMFTY(IJ, K)
     &                       + AMFTY(IJ+LS, K) + AMFTY(IJ+LN, K)
     &                       + AMFTY(IJ+LW, K) + AMFTY(IJ+LE, K) )
                  ELSE
                     YPSIX1(IJ, K) = YPSIX(IJ, K)
                  END IF
c                  IF (AMFTZ(IJ, K) .EQ. 1.0D0) THEN
c                     ZPSIX1(IJ, K) =
c     &                    ( 4.0D0 * ZPSIX(IJ, K) * AMFTZ(IJ, K)
c     &                       + ZPSIX(IJ+LS, K) * AMFTZ(IJ+LS, K)
c     &                       + ZPSIX(IJ+LN, K) * AMFTZ(IJ+LN, K)
c     &                       + ZPSIX(IJ+LW, K) * AMFTZ(IJ+LW, K)
c     &                       + ZPSIX(IJ+LE, K) * AMFTZ(IJ+LE, K) )
c     &                    / ( 4.0D0 * AMFTZ(IJ, K)
c     &                       + AMFTZ(IJ+LS, K) + AMFTZ(IJ+LN, K)
c     &                       + AMFTZ(IJ+LW, K) + AMFTZ(IJ+LE, K) )
c                     ZPSIY1(IJ, K) =
c     &                    ( 4.0D0 * ZPSIY(IJ, K) * AMFTZ(IJ, K)
c     &                       + ZPSIY(IJ+LS, K) * AMFTZ(IJ+LS, K)
c     &                       + ZPSIY(IJ+LN, K) * AMFTZ(IJ+LN, K)
c     &                       + ZPSIY(IJ+LW, K) * AMFTZ(IJ+LW, K)
c     &                       + ZPSIY(IJ+LE, K) * AMFTZ(IJ+LE, K) )
c     &                    / ( 4.0D0 * AMFTZ(IJ, K)
c     &                       + AMFTZ(IJ+LS, K) + AMFTZ(IJ+LN, K)
c     &                       + AMFTZ(IJ+LW, K) + AMFTZ(IJ+LE, K) )
c                  ELSE
c                     ZPSIX1(IJ, K) = ZPSIX(IJ, K)
c                     ZPSIY1(IJ, K) = ZPSIY(IJ, K)
c                  END IF
               END DO
            END DO
c            DO IJ = IJTSTR-NXDIM, IJTEND+NXDIM
            DO K = KSTR, KEND
               DO IJ = IJSTR, IJEND
                  XPSIY(IJ, K) = XPSIY1(IJ, K)
                  YPSIX(IJ, K) = YPSIX1(IJ, K)
c                  ZPSIX(IJ, K) = ZPSIX1(IJ, K)
c                  ZPSIY(IJ, K) = ZPSIY1(IJ, K)
               END DO 
            END DO
#ifdef OPT_PARALLEL
            CALL SHIFT2(
     M                   XPSIY,  YPSIX,
     D                   NXDIM,  NYDIM,  NZDIM)
c            CALL SHIFT2(
c     M                   ZPSIX,  ZPSIY,
c     D                   NXDIM,  NYDIM,  NZDIM)
#else
            CALL STBCAV(
     M                   XPSIY,  YPSIX)
c            CALL STBCAV(
c     M                   ZPSIX,  ZPSIY)
#endif
         END DO
      END IF
        
      DO K = KSTR+1, KEND
         DO IJ = IJTSTR, IJTEND
            ZPSIX(IJ, K) = 0.25D0 * 
     &         ( ( YPSIX(IJ, K-1) + YPSIX(IJ+LN, K-1) ) * DZ(IJ, K-1)
     &         + ( YPSIX(IJ, K  ) + YPSIX(IJ+LN, K  ) ) * DZ(IJ, K  ) )
     &           / DZM(IJ, K) * AMFTZ(IJ, K)
            ZPSIY(IJ, K) = 0.25D0 * 
     &         ( ( XPSIY(IJ, K-1) + XPSIY(IJ+LE, K-1) ) * DZ(IJ, K-1)
     &         + ( XPSIY(IJ, K  ) + XPSIY(IJ+LE, K  ) ) * DZ(IJ, K  ) )
     &           / DZM(IJ, K) * AMFTZ(IJ, K)
c            ZPSIX(IJ, K) = sign(1.0D0, ZPSIX(IJ, K)) *
c     &                     min(abs(ZPSIX(IJ, K)), 
c     &                         VSCL*0.5D0*(DZSIG(IJ, K-1)+DZSIG(IJ, K)))
c            ZPSIY(IJ, K) = sign(1.0D0, ZPSIY(IJ, K)) *
c     &                     min(abs(ZPSIY(IJ, K)), 
c     &                         VSCL*0.5D0*(DZSIG(IJ, K-1)+DZSIG(IJ, K)))
         END DO
      END DO

      CALL CHEKIN(CXPSY, 'CXPSY', NX, NY, 1, NXYDIM, 'SFC')
      CALL CHEKIN(CYPSX, 'CYPSX', NX, NY, 1, NXYDIM, 'SFC')
      CALL CHEKIN(XPSIY, 'XPSIY', NX, NY, NZ, NXYZDM, 'OCN')
      CALL CHEKIN(YPSIX, 'YPSIX', NX, NY, NZ, NXYZDM, 'OCN')
      CALL CHEKIN(ZPSIX, 'ZPSIX', NX, NY, NZ, NXYZDM, 'OCN')
      CALL CHEKIN(ZPSIY, 'ZPSIY', NX, NY, NZ, NXYZDM, 'OCN')
      CALL CHEKIN(RMAVDX, 'RMAVDX', NX, NY, 1, NXYDIM, 'SFC')
      CALL CHEKIN(RMAVDY, 'RMAVDY', NX, NY, 1, NXYDIM, 'SFC')
c      CALL CHEKIN(MUZX, 'MUZX', NX, NY, NZ, NXYZDM, 'OCN')
c      CALL CHEKIN(MUZY, 'MUZY', NX, NY, NZ, NXYZDM, 'OCN')

      DO 300 K = KSTR, KEND
         DO 310 IJ = IJTSTR, IJTEND+NXDIM
            DZDX = DTDX(IJ, K, 1) * 4.D0
     &             / (  DTFDZ(IJ, K  , 1) + DTFDZ(IJ+LW, K ,  1)
     &                + DTFDZ(IJ, K+1, 1) + DTFDZ(IJ+LW, K+1, 1) - EPS)
            XDZDX(IJ, K) = MIN(SLPMAX, MAX(-SLPMAX, DZDX))
            DZDY = DTDY(IJ, K, 1) * 4.D0 * DYM(IJ+LS)
     &        / (  (DTFDZ(IJ, K, 1) + DTFDZ(IJ, K+1, 1)) * DY(IJ)
     &        + (DTFDZ(IJ+LS, K, 1) + DTFDZ(IJ+LS, K+1, 1)) * DY(IJ+LS)
     &                - EPS)
            YDZDY(IJ, K) = MIN(SLPMAX, MAX(-SLPMAX, DZDY))
  310    CONTINUE
  300 CONTINUE

      DO 330 K = KSTR+1, KEND
         DO 340 IJ = IJTSTR, IJTEND
            DZDX = ((DTDX(IJ, K-1, 1)+DTDX(IJ+LE, K-1, 1))*DZ(IJ, K-1)
     &            + (DTDX(IJ, K  , 1)+DTDX(IJ+LE, K  , 1))*DZ(IJ, K)) *
     &             0.25D0 / DZM(IJ, K) / (DTFDZ(IJ, K, 1) - EPS)
            ZDZDX(IJ, K) = MIN(SLPMAX, MAX(-SLPMAX, DZDX)) *
     &                     AMFTZ(IJ, K)
            DZDY = 
     &        (  (DTDY(IJ, K-1, 1)+DTDY(IJ+LN, K-1, 1))*DZ(IJ, K-1)
     &         + (DTDY(IJ, K  , 1)+DTDY(IJ+LN, K  , 1))*DZ(IJ, K)) *
     &          0.25D0 / DZM(IJ, K) / (DTFDZ(IJ, K, 1) - EPS)
            ZDZDY(IJ, K) = MIN(SLPMAX, MAX(-SLPMAX, DZDY)) *
     &                     AMFTZ(IJ, K)
  340    CONTINUE
  330 CONTINUE

      DO 165 N = 1, NTDIM
         DO 160 K = 1, NZDIM
            DO 150 IJ = 1, NXYDIM
               DTDX(IJ, K, N) = 0.D0
               DTDY(IJ, K, N) = 0.D0
               DTFDZ(IJ, K, N) = 0.D0
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
               DTFDZ(IJ, K, N) = (TX(IJ, K-1, N) - TX(IJ, K, N)) 
     &                          / DZM(IJ, K) * AMFTZ(IJ, K)
  430       CONTINUE
  440    CONTINUE
  450 CONTINUE

      DO 520 N = 1, NTDIM
         DO 510 K = KSTR, KEND
            DO 500 IJ = IJTSTR, IJTEND+NXDIM
               XDTDZ(IJ, K, N) =
     &         (  DTFDZ(IJ, K, N)   + DTFDZ(IJ+LW, K, N)
     &          + DTFDZ(IJ, K+1, N) + DTFDZ(IJ+LW, K+1, N)) *  0.25D0
               YDTDZ(IJ, K, N) =
     &         (  (DTFDZ(IJ, K, N)    + DTFDZ(IJ, K+1, N)) *
     &            DY(IJ)
     &          + (DTFDZ(IJ+LS, K, N) + DTFDZ(IJ+LS, K+1, N)) *
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
      IMPLICIT NONE

#include "zocdim.F"
#include "zocgrd.F"
#include "zocmsk.F"

      REAL*8     ADT(NXYDIM, NZDIM, NTDIM)
      REAL*8   DIFFZ(NXYDIM, NZDIM)
      REAL*8      TX(NXYDIM, NZDIM, NTDIM),     TY(NXYDIM, NZDIM, NTDIM)
      REAL*8      UY(NXYDIM, NZDIM),     VY(NXYDIM, NZDIM)
      REAL*8       W(NXYDIM, NZDIM),    AHV(NXYDIM, NZDIM)

      REAL*8     FTX(NXYDIM, NZDIM, NTDIM)
      REAL*8     FTY(NXYDIM, NZDIM, NTDIM)
      REAL*8     FTZ(NXYDIM, NZDIM, NTDIM)
      COMMON /WORK/ FTX, FTY, FTZ

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
