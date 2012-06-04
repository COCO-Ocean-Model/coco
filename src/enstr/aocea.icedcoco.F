      SUBROUTINE OCSTUP(
     O               UB,     VB,     TB,     
     O               HB,   UBTB,   VBTB,      W,      R,
     O              AMV,    AHV,
     O               FT,   PTOP,
     I              DT1)
      use zocdim
      use zocgrd
      use zocfil
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  HISTORY
c     '02.10.10  H.Hasumi: from COCO3.4
c     '07.04.23  H.Hasumi
c     '07.09.25  H.Hasumi: arguments of CHEKIN
c     '07.10.04  H.Hasumi: 1-layer sea ice thermodynamics
c     '08.06.11  H.Hasumi: initial/final processing
c     '09.05.25  Y.Komuro: CMIP5 output code included
c     '09.10.06  Y.Komuro: FORSTO before PREDCI
c
c ---------------------------------------------------------------------

      REAL*8      TA(NXYZDM, NTDIM),     TB(NXYZDM, NTDIM)
      REAL*8      UA(NXYZDM),     UB(NXYZDM)
      REAL*8      VA(NXYZDM),     VB(NXYZDM)
      REAL*8      HA(NXYDIM),     HB(NXYDIM)
      REAL*8    UBTA(NXYDIM),   UBTB(NXYDIM)
      REAL*8    VBTA(NXYDIM),   VBTB(NXYDIM)
      REAL*8       W(NXYZDM),      R(NXYZDM)
      REAL*8     AMV(NXYZDM),    AHV(NXYZDM)

      REAL*8     AA(NXYDIM, 0:NIC),     AB(NXYDIM, 0:NIC)
      REAL*8    HIA(NXYDIM, 0:NIC),    HIB(NXYDIM, 0:NIC)
      REAL*8    UIA(NXYDIM),    UIB(NXYDIM)
      REAL*8    VIA(NXYDIM),    VIB(NXYDIM)
      REAL*8    TIA(NXYDIM, 0:NIC),    TIB(NXYDIM, 0:NIC)
      REAL*8    HSA(NXYDIM, 0:NIC),    HSB(NXYDIM, 0:NIC)
      REAL*8    TSI(NXYDIM, 0:NIC)

      REAL*8     FT(NXYDIM, NTDIM)
      REAL*8  SWABS(NXYDIM),     FS(NXYDIM)
      REAL*8   TAUX(NXYDIM),   TAUY(NXYDIM)
      REAL*8   PTOP(NXYDIM)

      REAL*8    UADV(NXYZDM),   VADV(NXYZDM),   WADV(NXYZDM, 9)
      REAL*8      GX(NXYZDM),     GY(NXYZDM)
      REAL*8      XX(NXYZDM),     YY(NXYZDM)
      REAL*8     GXX(NXYDIM),    GYY(NXYDIM)
      REAL*8     QAO(NXYDIM)
      REAL*8     QAI(NXYDIM, NIC),    QIO(NXYDIM, NIC)
      REAL*8     QII(NXYDIM, NIC)
      REAL*8     WEV(NXYDIM),    WSB(NXYDIM, NIC)
      REAL*8    PREC(NXYDIM),   SNOW(NXYDIM)
      REAL*8    ROFF(NXYDIM),   SOFF(NXYDIM)
      REAL*8  TAUAIX(NXYDIM), TAUAIY(NXYDIM)
      REAL*8  TAUAOX(NXYDIM), TAUAOY(NXYDIM)
      SAVE UADV, VADV, WADV, GX, GY, XX, YY, GXX, GYY,
     &     QAO, QAI, QIO, WEV, WSB, PREC, SNOW, ROFF, SOFF,
     &     TAUAIX, TAUAIY, TAUAOX, TAUAOY

      CHARACTER CTRNAM(NTDIM)*16, CFTNAM(NTDIM)*16
      SAVE CTRNAM, CFTNAM

#ifdef OPT_BODY
      REAL*8      TQ(NXYZDM, NTDIM)
      SAVE TQ
#endif

      REAL*8    SSFC(NXYDIM)

      REAL*8     AIG(NXYDIM),    HIG(NXYDIM),    HSG(NXYDIM)

      REAL*8     TT1,    TS1,    DT1,   TSS1
      INTEGER    NT1,   ITS1,  ITST1,  NTSS1
      LOGICAL OFLOUT(NFOMAX), OFLSTK(NFOMAX)

      INTEGER     IJ,      L
      INTEGER  IFPAR,  JFPAR

      CALL REWNML(IFPAR, JFPAR)
      WRITE(JFPAR, *) '*** OCSTUP ***'

      CTRNAM(1) = 'T'
      CTRNAM(2) = 'S'
      CFTNAM(1) = 'FT'
      CFTNAM(2) = 'FW'
      DO L = 3, NTDIM
         WRITE(CTRNAM(L), '(A6,I2.2)') 'TRACER', L
         WRITE(CFTNAM(L), '(A6,I2.2)') 'TRCFLX', L
      END DO

      DT = DT1

      CALL RDGEO
      CALL SVTSET
      CALL OVTSET(R, TB)
      CALL CHKSET

c *** Initialization for variables ***
#ifdef OPT_PARALLEL
      IF (MYRANK .LT. IJNODE) THEN
#endif
         CALL INISET(
     O               UADV,   VADV,   WADV,      R,
     I                 UB,     VB,     TB,
     I                 HB,   UBTB,   VBTB,
     I                  W,    AMV,    AHV)
#ifdef OPT_PARALLEL
      END IF
#endif

      DO 10 IJ = 1, NXYDIM
         GXX(IJ) = 0.D0
         GYY(IJ) = 0.D0
         PTOP(IJ) = 0.D0
   10 CONTINUE

      RETURN

c =====================================================================

      ENTRY OCEAN(
     M                   UA,     VA,     TA,
     M                   HA,   UBTA,   VBTA,
     M                   AB,    HIB,    UIB,    VIB,    TIB,    HSB,
     M                   UB,     VB,     TB,
     M                   HB,   UBTB,   VBTB,
     M                   AA,    HIA,    UIA,    VIA,    TIA,    HSA,
     M                    W,      R,    AMV,    AHV,
     M                   FT,  SWABS,     FS,
     M                 TAUX,   TAUY,   PTOP,    TSI,
     I                  NT1,    TT1,  ITST1,    TS1,   ITS1,
     I                NTSS1,   TSS1,
     I               OFLOUT, OFLSTK)

      IF (OINIT) THEN
         CALL PREDCI(
     M                   AB,    HIB,    UIB,    VIB,    TIB,    HSB,
     O                   FT,     FS,   TAUX,   TAUY,   PTOP,
     I                   AA,    HIA,    UIA,    VIA,    TIA,    HSA,
     I                   TB,     UB,     VB,     HB,     HA,
     I                  QAO,    QAI,    QII,    QIO,  SWABS,
     I                  WEV,    WSB,
     I                 PREC,   SNOW,   ROFF,   SOFF,
     I               TAUAIX, TAUAIY, TAUAOX, TAUAOY)
         CALL PREDCO(
     M                   HB,   UBTB,   VBTB,      W,      R,
     M                   UB,     VB,     TB,
     M                  AMV,    AHV,
     I                 TAUX,   TAUY,   PTOP,
     I                   FT,  SWABS,     FS,   SSFC,
#ifdef OPT_BODY
     I                   TQ,
#endif
     I                   UA,     VA,     TA,
     I                   HA,   UBTA,   VBTA,     AB,
     W                   GX,     GY,     XX,     YY,
     W                  GXX,    GYY,
     W                 UADV,   VADV,   WADV)
         RETURN
      END IF

      IF (OFINAL) THEN
#ifdef OPT_PARALLEL
         IF (MYRANK .LT. IJNODE) THEN
#endif
            CALL INISET(
     O               UADV,   VADV,   WADV,      R,
     I                 UB,     VB,     TB,
     I                 HB,   UBTB,   VBTB,
     I                  W,    AMV,    AHV)
#ifdef OPT_PARALLEL
         END IF
#endif
         CALL PREDCI(
     M                   AB,    HIB,    UIB,    VIB,    TIB,    HSB,
     O                   FT,     FS,   TAUX,   TAUY,   PTOP,
     I                   AA,    HIA,    UIA,    VIA,    TIA,    HSA,
     I                   TB,     UB,     VB,     HB,     HA,
     I                  QAO,    QAI,    QII,    QIO,  SWABS,
     I                  WEV,    WSB,
     I                 PREC,   SNOW,   ROFF,   SOFF,
     I               TAUAIX, TAUAIY, TAUAOX, TAUAOY)
         CALL PREDCO(
     M                   HB,   UBTB,   VBTB,      W,      R,
     M                   UB,     VB,     TB,
     M                  AMV,    AHV,
     I                 TAUX,   TAUY,   PTOP,
     I                   FT,  SWABS,     FS,   SSFC,
#ifdef OPT_BODY
     I                   TQ,
#endif
     I                   UA,     VA,     TA,
     I                   HA,   UBTA,   VBTA,     AB,
     W                   GX,     GY,     XX,     YY,
     W                  GXX,    GYY,
     W                 UADV,   VADV,   WADV)
         RETURN
      END IF

      TT    = TT1
      NT    = NT1
      TS    = TS1
      ITS   = ITS1
      ITST  = ITST1
      NTSS  = NTSS1
      TSS   = TSS1

      CALL CHKSTK(
     I            OFLSTK)

c *** Leap-frog ( A -> B ) ***

      IF ( ITST .EQ. 4 ) THEN
         CALL SFCFLX(
     O               QAO,    QAI,    QII,    QIO,  SWABS,    TSI,
     O               WEV,    WSB,
     O              PREC,   SNOW,   ROFF,   SOFF,
     O            TAUAIX, TAUAIY, TAUAOX, TAUAOY,
     O                FT,   PTOP,   SSFC,
     I                TB,     AB,    HIB,    TIB,    HSB)
#ifdef OPT_BODY
         CALL BDYFLX(
     O                TQ,
     I                TB)
#endif
         CALL PREDCI(
     M                   AB,    HIB,    UIB,    VIB,    TIB,    HSB,
     O                   FT,     FS,   TAUX,   TAUY,   PTOP,
     I                   AA,    HIA,    UIA,    VIA,    TIA,    HSA,
     I                   TB,     UB,     VB,     HB,     HA,
     I                  QAO,    QAI,    QII,    QIO,  SWABS,
     I                  WEV,    WSB,
     I                 PREC,   SNOW,   ROFF,   SOFF,
     I               TAUAIX, TAUAIY, TAUAOX, TAUAOY)
         CALL PREDCO(
     M                   HB,   UBTB,   VBTB,      W,      R,
     M                   UB,     VB,     TB,
     M                  AMV,    AHV,
     I                 TAUX,   TAUY,   PTOP,
     I                   FT,  SWABS,     FS,   SSFC,
#ifdef OPT_BODY
     I                   TQ,
#endif
     I                   UA,     VA,     TA,
     I                   HA,   UBTA,   VBTA,     AB,
     W                   GX,     GY,     XX,     YY,
     W                  GXX,    GYY,
     W                 UADV,   VADV,   WADV)

c *** Leap-frog ( B -> A ) ***

      ELSE IF ( ITST .EQ. 3 ) THEN
         CALL SFCFLX(
     O               QAO,    QAI,    QII,    QIO,  SWABS,    TSI,
     O               WEV,    WSB,
     O              PREC,   SNOW,   ROFF,   SOFF,
     O            TAUAIX, TAUAIY, TAUAOX, TAUAOY,
     O                FT,   PTOP,   SSFC,
     I                TA,     AA,    HIA,    TIA,    HSA)
#ifdef OPT_BODY
         CALL BDYFLX(
     O                TQ,
     I                TA)
#endif
         CALL PREDCI(
     M                   AA,    HIA,    UIA,    VIA,    TIA,    HSA,
     O                   FT,     FS,   TAUX,   TAUY,   PTOP,
     I                   AB,    HIB,    UIB,    VIB,    TIB,    HSB,
     I                   TA,     UA,     VA,     HA,     HB,
     I                  QAO,    QAI,    QII,    QIO,  SWABS,
     I                  WEV,    WSB,
     I                 PREC,   SNOW,   ROFF,   SOFF,
     I               TAUAIX, TAUAIY, TAUAOX, TAUAOY)
         CALL PREDCO(
     M                   HA,   UBTA,   VBTA,      W,      R,
     M                   UA,     VA,     TA,
     M                  AMV,    AHV,
     I                 TAUX,   TAUY,   PTOP,
     I                   FT,  SWABS,     FS,   SSFC,
#ifdef OPT_BODY
     I                   TQ,
#endif
     I                   UB,     VB,     TB,
     I                   HB,   UBTB,   VBTB,     AA,
     W                   GX,     GY,     XX,     YY,
     W                  GXX,    GYY,
     W                 UADV,   VADV,   WADV)

c *** Euler-backward ( B -> B ) ***

      ELSE IF ( ITST .EQ. 1 ) THEN
         CALL SFCFLX(
     O               QAO,    QAI,    QII,    QIO,  SWABS,    TSI,
     O               WEV,    WSB,
     O              PREC,   SNOW,   ROFF,   SOFF,
     O            TAUAIX, TAUAIY, TAUAOX, TAUAOY,
     O                FT,   PTOP,   SSFC,
     I                TB,     AB,    HIB,    TIB,    HSB)
#ifdef OPT_BODY
         CALL BDYFLX(
     O                TQ,
     I                TB)
#endif
         CALL FORSTI(
     O                   AA,    HIA,    UIA,    VIA,    TIA,    HSA,
     I                   AB,    HIB,    UIB,    VIB,    TIB,    HSB)
         CALL FORSTO(
     O                   UA,     VA,     TA,
     O                   HA,   UBTA,   VBTA,
     I                   UB,     VB,     TB,
     I                   HB,   UBTB,   VBTB)
         CALL PREDCI(
     M                   AA,    HIA,    UIA,    VIA,    TIA,    HSA,
     O                   FT,     FS,   TAUX,   TAUY,   PTOP,
     I                   AB,    HIB,    UIB,    VIB,    TIB,    HSB,
     I                   TA,     UA,     VA,     HA,     HB,
     I                  QAO,    QAI,    QII,    QIO,  SWABS,
     I                  WEV,    WSB,
     I                 PREC,   SNOW,   ROFF,   SOFF,
     I               TAUAIX, TAUAIY, TAUAOX, TAUAOY)
         CALL EXCNGI(
     O                   AA,    HIA,    UIA,    VIA,    TIA,    HSA,
     I                   AB,    HIB,    UIB,    VIB,    TIB,    HSB)
         IEULER = 1
         CALL PREDCO(
     M                   HA,   UBTA,   VBTA,      W,      R,
     M                   UA,     VA,     TA,
     M                  AMV,    AHV,
     I                 TAUX,   TAUY,   PTOP,
     I                   FT,  SWABS,     FS,   SSFC,
#ifdef OPT_BODY
     I                   TQ,
#endif
     I                   UB,     VB,     TB,
     I                   HB,   UBTB,   VBTB,     AB,
     W                   GX,     GY,     XX,     YY,
     W                  GXX,    GYY,
     W                 UADV,   VADV,   WADV)
         IEULER = 2
         CALL PREDCO(
     M                   HB,   UBTB,   VBTB,      W,      R,
     M                   UB,     VB,     TB,
     M                  AMV,    AHV,
     I                 TAUX,   TAUY,   PTOP,
     I                   FT,  SWABS,     FS,   SSFC,
#ifdef OPT_BODY
     I                   TQ,
#endif
     I                   UA,     VA,     TA,
     I                   HA,   UBTA,   VBTA,     AB,
     W                   GX,     GY,     XX,     YY,
     W                  GXX,    GYY,
     W                 UADV,   VADV,   WADV)

c *** Forward ( B -> B ) ***

      ELSE
         IEULER = 0
         CALL SFCFLX(
     O               QAO,    QAI,    QII,    QIO,  SWABS,    TSI,
     O               WEV,    WSB,
     O              PREC,   SNOW,   ROFF,   SOFF,
     O            TAUAIX, TAUAIY, TAUAOX, TAUAOY,
     O                FT,   PTOP,   SSFC,
     I                TB,     AB,    HIB,    TIB,    HSB)
#ifdef OPT_BODY
         CALL BDYFLX(
     O                TQ,
     I                TB)
#endif
         CALL FORSTI(
     O                   AA,    HIA,    UIA,    VIA,    TIA,    HSA,
     I                   AB,    HIB,    UIB,    VIB,    TIB,    HSB)
         CALL FORSTO(
     O                   UA,     VA,     TA,
     O                   HA,   UBTA,   VBTA,
     I                   UB,     VB,     TB,
     I                   HB,   UBTB,   VBTB)
         CALL PREDCI(
     M                   AA,    HIA,    UIA,    VIA,    TIA,    HSA,
     O                   FT,     FS,   TAUX,   TAUY,   PTOP,
     I                   AB,    HIB,    UIB,    VIB,    TIB,    HSB,
     I                   TA,     UA,     VA,     HA,     HB,
     I                  QAO,    QAI,    QII,    QIO,  SWABS,
     I                  WEV,    WSB,
     I                 PREC,   SNOW,   ROFF,   SOFF,
     I               TAUAIX, TAUAIY, TAUAOX, TAUAOY)
         CALL EXCNGI(
     O                   AA,    HIA,    UIA,    VIA,    TIA,    HSA,
     I                   AB,    HIB,    UIB,    VIB,    TIB,    HSB)
         CALL PREDCO(
     M                   HA,   UBTA,   VBTA,      W,      R,
     M                   UA,     VA,     TA,
     M                  AMV,    AHV,
     I                 TAUX,   TAUY,   PTOP,
     I                   FT,  SWABS,     FS,   SSFC,
#ifdef OPT_BODY
     I                   TQ,
#endif
     I                   UB,     VB,     TB,
     I                   HB,   UBTB,   VBTB,     AB,
     W                   GX,     GY,     XX,     YY,
     W                  GXX,    GYY,
     W                 UADV,   VADV,   WADV)
         CALL EXCNGO(
     M                   UA,     VA,     TA,
     M                   HA,   UBTA,   VBTA,
     M                   UB,     VB,     TB,
     M                   HB,   UBTB,   VBTB)
      END IF

c *** Output to file ***
      DO IJ = 1, NXYDIM
         AIG(IJ) = 0.0D0
         HIG(IJ) = 0.0D0
         HSG(IJ) = 0.0D0
      END DO

      IF ( ITST .EQ. 3 ) THEN
         CALL CHEKIN(    UA,    'U',
     &                   NX,     NY,     NZ, NXYZDM, 'OCN')
         CALL CHEKIN(    VA,    'V',
     &                   NX,     NY,     NZ, NXYZDM, 'OCN')
         DO L = 1, NTDIM
            CALL CHEKIN(     TA(1, L),      CTRNAM(L),
     &                     NX,     NY,     NZ, NXYZDM, 'OCN')
         END DO
         CALL CHEKIN(    HA,   'SH',
     &                   NX,     NY,      1, NXYDIM, 'SFC')
         CALL CHEKIN(  UBTA,  'UBT',
     &                   NX,     NY,      1, NXYDIM, 'SFC')
         CALL CHEKIN(  VBTA,  'VBT',
     &                   NX,     NY,      1, NXYDIM, 'SFC')
         CALL CHEKIN(    AA,   'AI',
     &                   NX,     NY,    NIC, NXYIDM, 'ICE')
         CALL CHEKIN(   HIA,   'HI',
     &                   NX,     NY,    NIC, NXYIDM, 'ICE')
         CALL CHEKIN(   UIA,   'UI',
     &                   NX,     NY,      1, NXYDIM, 'SFC')
         CALL CHEKIN(   VIA,   'VI',
     &                   NX,     NY,      1, NXYDIM, 'SFC')
         CALL CHEKIN(   TIA,   'TI',
     &                   NX,     NY,    NIC, NXYIDM, 'ICE')
         CALL CHEKIN(   HSA,   'HS',
     &                   NX,     NY,    NIC, NXYIDM, 'ICE')
         DO L = 1, NIC
            DO IJ = 1, NXYDIM
               AIG(IJ) = AIG(IJ) + AA(IJ, L)
               HIG(IJ) = HIG(IJ) + AA(IJ, L) * HIA(IJ, L)
               HSG(IJ) = HSG(IJ) + AA(IJ, L) * HSA(IJ, L)
            END DO
         END DO
      ELSE
         CALL CHEKIN(    UB,    'U',
     &                   NX,     NY,     NZ, NXYZDM, 'OCN')
         CALL CHEKIN(    VB,    'V',
     &                   NX,     NY,     NZ, NXYZDM, 'OCN')
         DO L = 1, NTDIM
            CALL CHEKIN(     TB(1, L),      CTRNAM(L),
     &                     NX,     NY,     NZ, NXYZDM, 'OCN')
         END DO
         CALL CHEKIN(    HB,   'SH',
     &                   NX,     NY,      1, NXYDIM, 'SFC')
         CALL CHEKIN(  UBTB,  'UBT',
     &                   NX,     NY,      1, NXYDIM, 'SFC')
         CALL CHEKIN(  VBTB,  'VBT',
     &                   NX,     NY,      1, NXYDIM, 'SFC')
         CALL CHEKIN(    AB,   'AI',
     &                   NX,     NY,    NIC, NXYIDM, 'ICE')
         CALL CHEKIN(   HIB,   'HI',
     &                   NX,     NY,    NIC, NXYIDM, 'ICE')
         CALL CHEKIN(   UIB,   'UI',
     &                   NX,     NY,      1, NXYDIM, 'SFC')
         CALL CHEKIN(   VIB,   'VI',
     &                   NX,     NY,      1, NXYDIM, 'SFC')
         CALL CHEKIN(   TIB,   'TI',
     &                   NX,     NY,    NIC, NXYIDM, 'ICE')
         CALL CHEKIN(   HSB,   'HS',
     &                   NX,     NY,    NIC, NXYIDM, 'ICE')
         DO L = 1, NIC
            DO IJ = 1, NXYDIM
               AIG(IJ) = AIG(IJ) + AB(IJ, L)
               HIG(IJ) = HIG(IJ) + AB(IJ, L) * HIB(IJ, L)
               HSG(IJ) = HSG(IJ) + AB(IJ, L) * HSB(IJ, L)
            END DO
         END DO
      END IF
         CALL CHEKIN(     W,    'W',
     &                   NX,     NY,     NZ, NXYZDM, 'OCN')
         CALL CHEKIN(   AMV,  'AMV',
     &                   NX,     NY,     NZ, NXYZDM, 'OCN')
         CALL CHEKIN(   AHV,  'AHV',
     &                   NX,     NY,     NZ, NXYZDM, 'OCN')
         DO L = 1, NTDIM
            CALL CHEKIN(      FT(1, L),      CFTNAM(L),
     &                      NX,     NY,      1, NXYDIM, 'SFC')
         END DO
         CALL CHEKIN( SWABS,'SWABS',
     &                   NX,     NY,      1, NXYDIM, 'SFC')
         CALL CHEKIN(    FS,    'FS',
     &                   NX,     NY,      1, NXYDIM, 'SFC')
         CALL CHEKIN(   TSI,   'TSI',
     &                   NX,     NY,    NIC, NXYIDM, 'ICE')

c        CHEKIN routine for FTX/Y/Z & FSX/Y/Z
         CALL CHKFTX

c        output section for CMIP5
c        AIG: total sea ice concentration, unit [ND]
c        HIG: sea ice thickness averaged over the entire grid cell
c             unit [cm]
c        HSG: mean snow thickness averaged over the entire grid cell
c             unit [cm]
         CALL CHEKIN(   AIG,  'AIG',
     &                   NX,     NY,      1, NXYDIM, 'SFC')
         CALL CHEKIN(   HIG,  'HIG',
     &                   NX,     NY,      1, NXYDIM, 'SFC')
         CALL CHEKIN(   HSG,  'HSG',
     &                   NX,     NY,      1, NXYDIM, 'SFC')

      CALL CHKOUT(
     I            OFLOUT)

      RETURN
      END
