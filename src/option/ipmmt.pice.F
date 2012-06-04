      SUBROUTINE PMOMNT(
     M                UIX,    VIX,
     O               TAUX,   TAUY,
     I                 AX,     AY,     AZ,
     I                HIX,    HIY,    HIZ,
     I                HSX,    HSY,    HSZ,
     I                UIY,    VIY,   PICE,
     I                 UX,     VX,     HY,   PTOP,
     I             TAUAIX, TAUAIY, TAUAOX, TAUAOY)
      use zocdim
      use zocgrd
      use zocmsk
      use zocphy
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Elastic-visco-plastic rheology of Hunke and Dukowicz (1997, JPO).
c
c  HISTORY
c     '03.06.05  H.Hasumi: from IcedCOCO3
c     '07.04.23  H.Hasumi
c     '07.05.01  H.Hasumi: sea surface slope term neglected
c     '07.09.25  H.Hasumi: multi-category sea ice
c     '08.06.11  H.Hasumi: initial/final processing
c     '08.07.10  H.Hasumi: initial/final processing
c     '09.01.09  Y.Komuro: bug fix (dimension of SGMXX etc. at FINADD)
c                          masking UX and VX, which are referred to
c                          in calculating ice-ocean stress
c     '09.02.20  Y.Komuro: parallel forward time march
c     '10.04.14  M.Kurogi: (COCO4.4 tripolar code by Dr. Suzuki)
c
c ---------------------------------------------------------------------

      REAL*8      AX(NXYDIM, 0:NIC),     AY(NXYDIM, 0:NIC)
      REAL*8      AZ(NXYDIM, 0:NIC)
      REAL*8     HIX(NXYDIM, 0:NIC),    HIY(NXYDIM, 0:NIC)
      REAL*8     HIZ(NXYDIM, 0:NIC)
      REAL*8     HSX(NXYDIM, 0:NIC),    HSY(NXYDIM, 0:NIC)
      REAL*8     HSZ(NXYDIM, 0:NIC)
      REAL*8     UIX(NXYDIM),    UIY(NXYDIM)
      REAL*8     VIX(NXYDIM),    VIY(NXYDIM)
      REAL*8      UX(NXYDIM, NZDIM),     VX(NXYDIM, NZDIM)
      REAL*8      HY(NXYDIM)
      REAL*8  TAUAIX(NXYDIM), TAUAIY(NXYDIM)
      REAL*8  TAUAOX(NXYDIM), TAUAOY(NXYDIM)
      REAL*8    TAUX(NXYDIM),   TAUY(NXYDIM)
      REAL*8    PTOP(NXYDIM)
      REAL*8    PICE(NXYDIM)

      REAL*8   SGMXX(NXYDIM),  SGMYY(NXYDIM),  SGMXY(NXYDIM)
      REAL*8    RDTS
      SAVE SGMXX, SGMYY, SGMXY, RDTS

      REAL*8    ZETA(NXYDIM),    ETA(NXYDIM)
      REAL*8     EMZ(NXYDIM),    EPZ(NXYDIM)
      REAL*8    ECOF(NXYDIM)
      REAL*8     EXX(NXYDIM),    EYY(NXYDIM),    EXY(NXYDIM)
      REAL*8    MICE(NXYDIM),   AICE(NXYDIM)
      REAL*8  AVRMSX(NXYDIM),   AVRA(NXYDIM)
      REAL*8  ACCELU(NXYDIM), ACCELV(NXYDIM)
      REAL*8    CTAU(NXYDIM),     HH(NXYDIM)
c      COMMON /WORK/ ZETA, ETA, EMZ, EPZ, ECOF,
c     &              EXX, EYY, EXY,
c     &              MICE, AICE, AVRMSX, AVRA,
c     &              ACCELU, ACCELV, CTAU, HH

      REAL*8    CAIC,   CAIS,   CIOC,   CIOS
      REAL*8    ELST(NXYDIM)
      REAL*8    CEMZ,   CEPZ,  RECC2
      SAVE CAIC, CAIS, CIOC, CIOS, ELST, CEMZ, CEPZ, RECC2

      INTEGER     IJ,      K, ISPLIT, NNSPLT
      INTEGER   IJLE,   IJLN,  IJLNE
      REAL*8     RAD
      REAL*8    ALPH,   BETA,   GAMM,   DELT
      REAL*8   RAABB
      INTEGER  IFPAR,  JFPAR

      LOGICAL OFIRST,   OEOF
      SAVE OFIRST, OEOF

      REAL*8  THETAA, THETAO, CWDRAG, ELAST0
      INTEGER  KGLEV
      INTEGER NSPLIT
      REAL*8     ECC,   DMIN,  FLOSS
      SAVE CWDRAG, KGLEV
      NAMELIST /NMIDYN/ ECC, DMIN, FLOSS
      NAMELIST /NMIPRM/ THETAA, THETAO, CWDRAG, ELAST0, KGLEV
      NAMELIST /NMITSP/ NSPLIT

      DATA ECC, DMIN, FLOSS / 2.D0, 2.0D-7, 17.D0 /
      DATA THETAA, THETAO / 0.D0, 25.D0 /
      DATA KGLEV / 1 /
      DATA CWDRAG / 5.D-3 /
      DATA ELAST0 / 0.25D0 /
      DATA NSPLIT / 60 /
      DATA OFIRST / .TRUE. /

      IF (OINIT) THEN
#ifdef OPT_TRIPOLE
         CALL RSTADD(SGMXX, OEOF, NXDIM, NYDIM, 1, 'SGMXX', 'SFC',
     &                                              1.D0,  0,  0 )
         CALL RSTADD(SGMYY, OEOF, NXDIM, NYDIM, 1, 'SGMYY', 'SFC',
     &                                              1.D0,  0,  0 )
         CALL RSTADD(SGMXY, OEOF, NXDIM, NYDIM, 1, 'SGMXY', 'SFC',
     &                                              1.D0,  0,  0 )
#else
         CALL RSTADD(SGMXX, OEOF, NXDIM, NYDIM, 1, 'SGMXX', 'SFC')
         CALL RSTADD(SGMYY, OEOF, NXDIM, NYDIM, 1, 'SGMYY', 'SFC')
         CALL RSTADD(SGMXY, OEOF, NXDIM, NYDIM, 1, 'SGMXY', 'SFC')
#endif
         RETURN
      END IF

      IF (OFINAL) THEN
         CALL FINADD(SGMXX, NXDIM, NYDIM, 1, 'SGMXX', 'SFC')
         CALL FINADD(SGMYY, NXDIM, NYDIM, 1, 'SGMYY', 'SFC')
         CALL FINADD(SGMXY, NXDIM, NYDIM, 1, 'SGMXY', 'SFC')
         RETURN
      END IF

      IF (OFIRST) THEN
ccc         OFIRST = .FALSE.
         CALL REWNML(IFPAR, JFPAR)
         READ (IFPAR, NMIPRM, END=99)
 99      WRITE(JFPAR, NMIPRM)
         CALL REWNML(IFPAR, JFPAR)
         READ (IFPAR, NMITSP, END=98)
 98      WRITE(JFPAR, NMITSP)
         CALL REWNML(IFPAR, JFPAR)
         READ (IFPAR, NMIDYN, END=97)
 97      WRITE(JFPAR, NMIDYN)

         CEMZ = (1.D0 / ECC / ECC - 1.D0) * 0.5D0
         CEPZ = (1.D0 / ECC / ECC + 1.D0) * 0.5D0
         RECC2 = 0.5D0 / ECC / ECC
c         RDTS = DBLE(NSPLIT) / DT
         RAD = ATAN(1.D0) * 4.D0 / 180.D0
         KGLEV = KGLEV + KSTR - 1
         CAIC = COS(RAD * THETAA)
         CAIS = SIN(RAD * THETAA)
         CIOC = COS(RAD * THETAO)
         CIOS = SIN(RAD * THETAO)
         CWDRAG = rhoo * CWDRAG
         DO 40 IJ = 1, NXYDIM
            ELST(IJ) = 2.D0 * rhoi * ELAST0 *
     &                 MIN(DX*DX*HXT(IJ)*HXT(IJ),
     &                     DY(IJ)*DY(IJ)*HYT(IJ)*HYT(IJ)) *
     &                 AMSKT(IJ, KSTR)
 40      CONTINUE
      END IF
      
      IF (ITS .EQ. 2) THEN
         NNSPLT = NSPLIT
      ELSE
         NNSPLT = NSPLIT / 2
      ENDIF
      RDTS = DBLE(NNSPLT) / TS

      DO 100 IJ = 1, NXYDIM
         ZETA  (IJ) = 0.D0
         ETA   (IJ) = 0.D0
         EMZ   (IJ) = 0.D0
         EPZ   (IJ) = 0.D0
         EXX   (IJ) = 0.D0
         EYY   (IJ) = 0.D0
         EXY   (IJ) = 0.D0
         AVRMSX(IJ) = 0.D0
         AVRA  (IJ) = 0.D0
         ACCELU(IJ) = 0.D0
         ACCELV(IJ) = 0.D0
         CTAU  (IJ) = 0.D0
         ECOF  (IJ) = 0.D0
         AICE  (IJ) = 1.D0 - AX(IJ, 0)
         MICE  (IJ) = 0.D0
         TAUAIX(IJ) = TAUAIX(IJ) * AICE(IJ)
         TAUAIY(IJ) = TAUAIY(IJ) * AICE(IJ)
         TAUAOX(IJ) = TAUAOX(IJ) * AX(IJ, 0)
         TAUAOY(IJ) = TAUAOY(IJ) * AX(IJ, 0)
c         HH    (IJ) = HY(IJ)
c     &              + AY(IJ) * (rhoi * HIY(IJ) + rhos * HSY(IJ)) / rhoo
         HH    (IJ) = 0.D0
  100 CONTINUE
      DO K = 1, NIC
         DO IJ = 1, NXYDIM
            MICE(IJ) = MICE(IJ) + AX(IJ, K) * HIX(IJ, K)
         END DO
      END DO

c      DO 500 IJ = 1, NXYDIM
c         UIX(IJ) = UIY(IJ)
c         VIX(IJ) = VIY(IJ)
c  500 CONTINUE
      CALL STRAIN(
     O               EXX,    EYY,    EXY,
     I               UIX,    VIX)
      CALL RHEOLO(
     O              ZETA,    ETA,    EMZ,    EPZ,   PICE,
     I              AICE,   MICE,
     I               EXX,    EYY,    EXY)

      IF (OFIRST) THEN
         OFIRST = .FALSE.
         IF (OEOF) THEN
            DO 50 IJ = 1, NXYDIM
               SGMXX(IJ) = 2.D0 * ETA(IJ) * EXX(IJ)
     &                - EMZ(IJ) * (EXX(IJ) + EYY(IJ))
     &                - 0.5D0 * PICE(IJ)
               SGMYY(IJ) = 2.D0 * ETA(IJ) * EYY(IJ)
     &                - EMZ(IJ) * (EXX(IJ) + EYY(IJ))
     &                - 0.5D0 * PICE(IJ)
               SGMXY(IJ) = 2.D0 * ETA(IJ) * EXY(IJ)
 50         CONTINUE
         END IF
      END IF

      DO 400 IJ = IJVSTR, IJVEND
         AVRMSX(IJ) = (  MICE(IJ) + MICE(IJ+LE)
     &                 + MICE(IJ+LN) + MICE(IJ+LNE)) * 0.25D0 * rhoi *
     &                AMSKV(IJ, KSTR)
         AVRA(IJ)   = (  AICE(IJ) + AICE(IJ+LE)
     &                 + AICE(IJ+LN) + AICE(IJ+LNE)) * 0.25D0
         ACCELU(IJ) = (  CAIC * TAUAIX(IJ)
     &                 - SIGN(CAIS, COR(IJ)) * TAUAIY(IJ)
     &                 - AVRMSX(IJ) * RX * RXU(IJ) * 0.5D0 *
     &                   (  HH(IJ+LNE) + HH(IJ+LE)
     &                    - HH(IJ+LN) - HH(IJ)) * gravit
     &                ) * AMSKV(IJ, KSTR)
         ACCELV(IJ) = (  CAIC * TAUAIY(IJ)
     &                 + SIGN(CAIS, COR(IJ)) * TAUAIX(IJ)
     &                 - AVRMSX(IJ) * RYM(IJ) * RYU(IJ) * 0.5D0 *
     &                   (  HH(IJ+LNE) + HH(IJ+LN)
     &                    - HH(IJ+LE) - HH(IJ)) * gravit
     &                ) * AMSKV(IJ, KSTR)
  400 CONTINUE

c      DO 1000 ISPLIT = 1, NSPLIT
      DO 1000 ISPLIT = 1, NNSPLT

         DO 600 IJ = IJTSTR, IJTEND
            IF (AICE(IJ) .EQ. 0.D0) THEN
               SGMXX(IJ) = 0.D0
               SGMYY(IJ) = 0.D0
               SGMXY(IJ) = 0.D0
            ELSE
               ECOF(IJ) = ELST(IJ) * MICE(IJ) * RDTS * RDTS
               ALPH = 2.D0 * ETA(IJ) * RDTS + CEPZ * ECOF(IJ)
               BETA = CEMZ * ECOF(IJ)
               GAMM = 2.D0 * ETA(IJ) *
     &                (SGMXX(IJ) * RDTS + ECOF(IJ) * EXX(IJ))
     &              - ECOF(IJ) * PICE(IJ) * RECC2
               DELT = 2.D0 * ETA(IJ) *
     &                (SGMYY(IJ) * RDTS + ECOF(IJ) * EYY(IJ))
     &              - ECOF(IJ) * PICE(IJ) * RECC2
               RAABB = 1.D0 / (ALPH * ALPH - BETA * BETA)
               SGMXX(IJ) = (ALPH * GAMM - BETA * DELT) * RAABB
               SGMYY(IJ) = (ALPH * DELT - BETA * GAMM) * RAABB
               SGMXY(IJ) = 2.D0 * ETA(IJ) *
     &                     (SGMXY(IJ) * RDTS + ECOF(IJ) * EXY(IJ))
     &                     / (2.D0 * ETA(IJ) * RDTS + ECOF(IJ))
            END IF
  600    CONTINUE

#ifdef OPT_PARALLEL
#ifdef OPT_TRIPOLE
         CALL SHIFT3( SGMXX,  SGMYY,  SGMXY,   
     D                NXDIM,  NYDIM,      1,
     I                 1.D0,      0,      0 )
#else
         CALL SHIFT3(
     M                SGMXX,  SGMYY,  SGMXY,
     D                NXDIM,  NYDIM,      1)
#endif
#endif

         DO 700 IJ = IJVSTR, IJVEND
            IF (AVRMSX(IJ) .EQ. 0.D0) THEN
               UIX(IJ) = 0.D0
               VIX(IJ) = 0.D0
            ELSE
               IJLE = IJ + LE
               IJLN = IJ + LN
               IJLNE = IJ + LNE
               CTAU(IJ) = AVRA(IJ) * CWDRAG *
     &                    SQRT(  (UIX(IJ)
     &                           - UX(IJ, KGLEV)*AMSKV(IJ, KGLEV) )**2
     &                         + (VIX(IJ)
     &                           - VX(IJ, KGLEV)*AMSKV(IJ, KGLEV) )**2)
               ALPH = AVRMSX(IJ) * RDTS
     &              + CTAU(IJ) * CIOC
               BETA = AVRMSX(IJ) * COR(IJ)
     &              + CTAU(IJ) * SIGN(CIOS, COR(IJ))
               GAMM = (  (  SGMXX(IJLNE) * HYT(IJLNE)
     &                    + SGMXX(IJLE ) * HYT(IJLE )
     &                    - SGMXX(IJLN ) * HYT(IJLN )
     &                    - SGMXX(IJ   ) * HYT(IJ   )) * RX
     &                 + (  SGMXY(IJLNE) * HXT(IJLNE)
     &                    + SGMXY(IJLN ) * HXT(IJLN )
     &                    - SGMXY(IJLE ) * HXT(IJLE )
     &                    - SGMXY(IJ   ) * HXT(IJ   )) * RYM(IJ)) *
     &                0.5D0 * RXU(IJ) * RYU(IJ)
     &              + (  SGMXY(IJLNE) + SGMXY(IJLE)
     &                 + SGMXY(IJLN) + SGMXY(IJ)) * 0.25D0 * HXYU(IJ)
     &              - (  SGMYY(IJLNE) + SGMYY(IJLE)
     &                 + SGMYY(IJLN) + SGMYY(IJ)) * 0.25D0 * HYXU(IJ)
     &              + ACCELU(IJ)
     &              + CTAU(IJ) *
     &                (  UX(IJ, KGLEV)*AMSKV(IJ, KGLEV) * CIOC
     &                 - VX(IJ, KGLEV)*AMSKV(IJ, KGLEV)
     &                     * SIGN(CIOS, COR(IJ)))
     &              + AVRMSX(IJ) * RDTS * UIX(IJ)
               DELT = (  (  SGMXY(IJLNE) * HYT(IJLNE)
     &                    + SGMXY(IJLE ) * HYT(IJLE )
     &                    - SGMXY(IJLN ) * HYT(IJLN )
     &                    - SGMXY(IJ   ) * HYT(IJ   )) * RX
     &                 + (  SGMYY(IJLNE) * HXT(IJLNE)
     &                    + SGMYY(IJLN ) * HXT(IJLN )
     &                    - SGMYY(IJLE ) * HXT(IJLE )
     &                    - SGMYY(IJ   ) * HXT(IJ   )) * RYM(IJ)) *
     &                0.5D0 * RXU(IJ) * RYU(IJ)
     &              + (  SGMXY(IJLNE) + SGMXY(IJLE)
     &                 + SGMXY(IJLN) + SGMXY(IJ)) * 0.25D0 * HYXU(IJ)
     &              - (  SGMXX(IJLNE) + SGMXX(IJLE)
     &                 + SGMXX(IJLN) + SGMXX(IJ)) * 0.25D0 * HXYU(IJ)
     &              + ACCELV(IJ)
     &              + CTAU(IJ) *
     &                (  VX(IJ, KGLEV)*AMSKV(IJ, KGLEV) * CIOC
     &                 + UX(IJ, KGLEV)*AMSKV(IJ, KGLEV)
     &                     * SIGN(CIOS, COR(IJ)))
     &              + AVRMSX(IJ) * RDTS * VIX(IJ)
               RAABB = 1.D0 / (ALPH * ALPH + BETA * BETA)
               UIX(IJ) = (ALPH * GAMM + BETA * DELT) * RAABB
               VIX(IJ) = (ALPH * DELT - BETA * GAMM) * RAABB
            END IF
  700    CONTINUE

#ifdef OPT_PARALLEL
#ifdef OPT_TRIPOLE
         CALL SHIFT2(   UIX,    VIX,
     D                NXDIM,  NYDIM,      1,
     I                -1.D0,     -1,     -1 )
#else
         CALL SHIFT2(
     M                  UIX,    VIX,
     D                NXDIM,  NYDIM,      1)
#endif
#else
         CALL STBCVI(
     M               UIX,    VIX)
#endif

         CALL STRAIN(
     O                  EXX,    EYY,    EXY,
     I                  UIX,    VIX)
         CALL RHEOLO(
     O                 ZETA,    ETA,    EMZ,    EPZ,   PICE,
     I                 AICE,   MICE,
     I                  EXX,    EYY,    EXY)

 1000 CONTINUE

      DO 800 IJ = IJVSTR, IJVEND
         TAUX(IJ) = (  CTAU(IJ) *
     &                 (  (UIX(IJ)
     &                      - UX(IJ, KGLEV)*AMSKV(IJ, KGLEV)) * CIOC
     &                  - (VIX(IJ)
     &                      - VX(IJ, KGLEV)*AMSKV(IJ, KGLEV)) *
     &                    SIGN(CIOS, COR(IJ)))
     &               + TAUAOX(IJ) * CAIC
     &               - TAUAOY(IJ) * SIGN(CAIS, COR(IJ))) *
     &              AMSKV(IJ, KSTR)
         TAUY(IJ) = (  CTAU(IJ) *
     &                 (  (VIX(IJ)
     &                      - VX(IJ, KGLEV)*AMSKV(IJ, KGLEV)) * CIOC
     &                  + (UIX(IJ)
     &                      - UX(IJ, KGLEV)*AMSKV(IJ, KGLEV)) *
     &                    SIGN(CIOS, COR(IJ)))
     &               + TAUAOY(IJ) * CAIC
     &               + TAUAOX(IJ) * SIGN(CAIS, COR(IJ))) *
     &              AMSKV(IJ, KSTR)
  800 CONTINUE

      RETURN
      END

c *********************************************************************

      SUBROUTINE STRAIN(
     O                EXX,    EYY,    EXY,
     I                 UI,     VI)
      use zocdim
      use zocgrd
      use zocmsk
      IMPLICIT NONE

      REAL*8     EXX(NXYDIM),    EYY(NXYDIM),    EXY(NXYDIM)
      REAL*8      UI(NXYDIM),     VI(NXYDIM)

      INTEGER     IJ,   IJLW,   IJLS,  IJLSW

      DO 100 IJ = IJTSTR, IJTEND+NXDIM+1
         IJLW = IJ + LW
         IJLS = IJ + LS
         IJLSW = IJ + LSW
         EXX(IJ) = (UI(IJ) + UI(IJLS) - UI(IJLW) - UI(IJLSW)) *
     &             RX * RXT(IJ) * 0.5D0
     &           + (VI(IJ) + VI(IJLS) + VI(IJLW) + VI(IJLSW)) *
     &             HXYT(IJ) * 0.25D0
         EYY(IJ) = (VI(IJ) + VI(IJLW) - VI(IJLS) - VI(IJLSW)) *
     &             RY(IJ) * RYT(IJ) * 0.5D0
     &           + (UI(IJ) + UI(IJLW) + UI(IJLS) + UI(IJLSW)) *
     &             HYXT(IJ) * 0.25D0
         EXY(IJ) = (UI(IJ) + UI(IJLW) - UI(IJLS) - UI(IJLSW)) *
     &             RY(IJ) * RYT(IJ) * 0.25D0
     &           - (VI(IJ) + VI(IJLW) + VI(IJLS) + VI(IJLSW)) *
     &             HYXT(IJ) * 0.125D0
     &           + (VI(IJ) + VI(IJLS) - VI(IJLW) - VI(IJLSW)) *
     &             RX * RXT(IJ) * 0.25D0
     &           - (UI(IJ) + UI(IJLS) + UI(IJLW) + UI(IJLSW)) *
     &             HXYT(IJ) * 0.125D0
  100 CONTINUE

      RETURN
      END

c *********************************************************************

      SUBROUTINE RHEOLO(
     O               ZETA,    ETA,    EMZ,    EPZ,   PICE,
     I               AICE,   MICE,
     I                EXX,    EYY,    EXY)
      use zocdim
      use zocmsk
      IMPLICIT NONE

      REAL*8    ZETA(NXYDIM),    ETA(NXYDIM)
      REAL*8     EMZ(NXYDIM),    EPZ(NXYDIM)
      REAL*8     EXX(NXYDIM),    EYY(NXYDIM),    EXY(NXYDIM)
      REAL*8    AICE(NXYDIM),   MICE(NXYDIM)
      REAL*8    PICE(NXYDIM)

      REAL*8   DELTA(NXYDIM)

      REAL*8      C1,     C2,     C3,     C4
      LOGICAL OFIRST
      SAVE C1, C2, C3, C4, OFIRST

      REAL*8     DEL
      INTEGER     IJ
      INTEGER  IFPAR,  JFPAR

      REAL*8     ECC,   DMIN,  FLOSS
      REAL*8      P0,     CP
      SAVE DMIN, P0, CP
      NAMELIST /NMIDYN/ ECC, DMIN, FLOSS
      NAMELIST /NMPICE/ P0, CP

      DATA OFIRST / .TRUE. /
      DATA ECC, DMIN, FLOSS / 2.D0, 2.0D-7, 17.D0 /
      DATA P0, CP / 2.D+5, 2.D+1 /

      IF (OFIRST) THEN
         OFIRST = .FALSE.
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMIDYN, END=99)
 99      WRITE(JFPAR, NMIDYN)
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMPICE, END=98)
 98      WRITE(JFPAR, NMPICE)

         C1 = 1.D0 + 1.D0 / ECC / ECC
         C2 = 4.D0 / ECC / ECC
         C3 = 2.D0 * (1.D0 - 1.D0 / ECC / ECC)
         C4 = 1.D0 / ECC / ECC
      END IF

      DO 100 IJ = IJTSTR, IJTEND+NXDIM+1
         PICE(IJ) = P0 * MICE(IJ) * EXP(- CP * (1.D0 - AICE(IJ)))
         DEL = SQRT(  C1 * (EXX(IJ) * EXX(IJ) + EYY(IJ) * EYY(IJ))
     &              + C2 * EXY(IJ) * EXY(IJ)
     &              + C3 * EXX(IJ) * EYY(IJ))
         DELTA(IJ) = MAX(DEL, DMIN)
         ZETA(IJ) = PICE(IJ) / DELTA(IJ) * 0.5D0
         ETA (IJ) = ZETA(IJ) * C4
         EMZ (IJ) = ETA(IJ) - ZETA(IJ)
         EPZ (IJ) = ETA(IJ) + ZETA(IJ)
  100 CONTINUE

      RETURN
      END
