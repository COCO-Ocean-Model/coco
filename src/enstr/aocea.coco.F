      SUBROUTINE OCSTUP(
     O               UB,     VB,     TB,     
     O               HB,   UBTB,   VBTB,      W,      R,
     O              AMV,    AHV,
     O               FT,
     I              DT1)
      use zocdim
      use zocgrd
      use zocfil
      use tflxt
      use tovtr
      use brdge
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  HISTORY
c     '00.03.17  H.Hasumi: for Ishizaki scheme
c     '00.05.30  H.Hasumi: parallelized COCO3
c     '00.12.07  H.Hasumi: combine par and non-par routines
c     '01.01.23  H.Hasumi: for hybrid vertical coordinate
c                          zonal filter
c     '01.02.02  H.Hasumi: divide the tracer routine
c     '01.02.04  H.Hasumi: combine body routines
c     '01.02.20  H.Hasumi: for incorporating BBL model
c     '01.12.03  H.Hasumi: subroutine FVSET is split into two
c     '01.12.07  H.Hasumi
c     '02.05.29  H.Nakano: tracer dimension
c     '07.04.23  H.Hasumi
c     '07.09.25  H.Hasumi: arguments of CHEKIN
c     '08.06.11  H.Hasumi: initial/final processing
c     '08.07.10  H.Hasumi: initial/final processing
c     '10.04.14  M.Kurogi: (COCO4.4 tripolar code by Dr. Suzuki)
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
      REAL*8      FT(NXYDIM, NTDIM)

      REAL*8    UADV(NXYZDM),   VADV(NXYZDM),   WADV(NXYZDM, 9)
      REAL*8      GX(NXYZDM),     GY(NXYZDM)
      REAL*8      XX(NXYZDM),     YY(NXYZDM)
      REAL*8     GXX(NXYDIM),    GYY(NXYDIM)
      REAL*8    TAUX(NXYDIM),   TAUY(NXYDIM)
      REAL*8    PTOP(NXYDIM)
      SAVE UADV, VADV, WADV, GX, GY, XX, YY, GXX, GYY, TAUX, TAUY, PTOP

      CHARACTER CTRNAM(NTDIM)*16, CFTNAM(NTDIM)*16
      SAVE CTRNAM, CFTNAM

#ifdef OPT_BODY
      REAL*8      TQ(NXYZDM, NTDIM)
      SAVE TQ
#endif

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
     I                  W,    AMV,    AHV,
     I               TAUX,   TAUY,
     I                 FT)
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
     M                UA,     VA,     TA,
     M                HA,   UBTA,   VBTA,
     M                UB,     VB,     TB,
     M                HB,   UBTB,   VBTB,
     M                 W,      R,    AMV,    AHV,
     M                FT,    
     I               NT1,    TT1,  ITST1,    TS1,   ITS1,
     I             NTSS1,   TSS1,
     G            OFLOUT, OFLSTK)

      IF (OINIT) THEN
         CALL PREDCO(
     M                HB,   UBTB,   VBTB,      W,      R,
     M                UB,     VB,     TB,
     M               AMV,    AHV,
     I              TAUX,   TAUY,     FT,   PTOP,
#ifdef OPT_BODY
     I                TQ,
#endif
     I                UA,     VA,     TA,
     I                HA,   UBTA,   VBTA,
     W                GX,     GY,     XX,     YY,
     W               GXX,    GYY,
     W              UADV,   VADV,   WADV)
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
     I                  W,    AMV,    AHV,
     I               TAUX,   TAUY,
     I                 FT)
#ifdef OPT_PARALLEL
         END IF
#endif
         CALL PREDCO(
     M                HB,   UBTB,   VBTB,      W,      R,
     M                UB,     VB,     TB,
     M               AMV,    AHV,
     I              TAUX,   TAUY,     FT,   PTOP,
#ifdef OPT_BODY
     I                TQ,
#endif
     I                UA,     VA,     TA,
     I                HA,   UBTA,   VBTA,
     W                GX,     GY,     XX,     YY,
     W               GXX,    GYY,
     W              UADV,   VADV,   WADV)
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
     O                FT,   TAUX,   TAUY,
     I                TB)
#ifdef OPT_BODY
         CALL BDYFLX(
     O                TQ,
     I                TB)
#endif
#ifdef OPT_PARALLEL
#ifdef OPT_TRIPOLE
      CALL SHIFT2(
     M              TAUX,   TAUY,
     D             NXDIM,  NYDIM,      1,
     I             -1.D0,      0,      0 )
#else
         CALL SHIFT2(
     M                 TAUX,   TAUY,
     D                NXDIM,  NYDIM,      1)
#endif
#else
         CALL STBCTU(
     M                 TAUX,   TAUY)
#endif
         CALL PREDCO(
     M                HB,   UBTB,   VBTB,      W,      R,
     M                UB,     VB,     TB,
     M               AMV,    AHV,
     I              TAUX,   TAUY,     FT,   PTOP,
#ifdef OPT_BODY
     I                TQ,
#endif
     I                UA,     VA,     TA,
     I                HA,   UBTA,   VBTA,
     W                GX,     GY,     XX,     YY,
     W               GXX,    GYY,
     W              UADV,   VADV,   WADV)

c *** Leap-frog ( B -> A ) ***

      ELSE IF ( ITST .EQ. 3 ) THEN
         CALL SFCFLX(
     O                FT,   TAUX,   TAUY,
     I                TA)
#ifdef OPT_BODY
         CALL BDYFLX(
     O                TQ,
     I                TA)
#endif
#ifdef OPT_PARALLEL
#ifdef OPT_TRIPOLE
      CALL SHIFT2(
     M              TAUX,   TAUY,
     D             NXDIM,  NYDIM,      1,
     I             -1.D0,      0,      0 )
#else
         CALL SHIFT2(
     M                 TAUX,   TAUY,
     D                NXDIM,  NYDIM,      1)
#endif
#else
         CALL STBCTU(
     M                 TAUX,   TAUY)
#endif
         CALL PREDCO(
     M                HA,   UBTA,   VBTA,      W,      R,
     M                UA,     VA,     TA,
     M               AMV,    AHV,
     I              TAUX,   TAUY,     FT,   PTOP,
#ifdef OPT_BODY
     I                TQ,
#endif
     I                UB,     VB,     TB,
     I                HB,   UBTB,   VBTB,
     W                GX,     GY,     XX,     YY,
     W               GXX,    GYY,
     W              UADV,   VADV,   WADV)

c *** Euler-backward ( B -> B ) ***

      ELSE IF ( ITST .EQ. 1 ) THEN
         CALL FORSTO(
     O                UA,     VA,     TA,
     O                HA,   UBTA,   VBTA,
     I                UB,     VB,     TB,
     I                HB,   UBTB,   VBTB)
         CALL SFCFLX(
     O                FT,   TAUX,   TAUY,
     I                TB)
#ifdef OPT_BODY
         CALL BDYFLX(
     O                TQ,
     I                TB)
#endif
#ifdef OPT_PARALLEL
#ifdef OPT_TRIPOLE
      CALL SHIFT2(
     M              TAUX,   TAUY,
     D             NXDIM,  NYDIM,      1,
     I             -1.D0,      0,      0 )
#else
         CALL SHIFT2(
     M                 TAUX,   TAUY,
     D                NXDIM,  NYDIM,      1)
#endif
#else
         CALL STBCTU(
     M                 TAUX,   TAUY)
#endif
         IEULER = 1
         CALL PREDCO(
     M                HA,   UBTA,   VBTA,      W,      R,
     M                UA,     VA,     TA,
     M               AMV,    AHV,
     I              TAUX,   TAUY,     FT,   PTOP,
#ifdef OPT_BODY
     I                TQ,
#endif
     I                UB,     VB,     TB,
     I                HB,   UBTB,   VBTB,
     W                GX,     GY,     XX,     YY,
     W               GXX,    GYY,
     W              UADV,   VADV,   WADV)
         IEULER = 2
         CALL PREDCO(
     M                HB,   UBTB,   VBTB,      W,      R,
     M                UB,     VB,     TB,
     M               AMV,    AHV,
     I              TAUX,   TAUY,     FT,   PTOP,
#ifdef OPT_BODY
     I                TQ,
#endif
     I                UA,     VA,     TA,
     I                HA,   UBTA,   VBTA,
     W                GX,     GY,     XX,     YY,
     W               GXX,    GYY,
     W              UADV,   VADV,   WADV)

c *** Forward ( B -> B ) ***

      ELSE
         IEULER = 0
         CALL FORSTO(
     O                UA,     VA,     TA,
     O                HA,   UBTA,   VBTA,
     I                UB,     VB,     TB,
     I                HB,   UBTB,   VBTB)
         CALL SFCFLX(
     O                FT,   TAUX,   TAUY,
     I                TB)
#ifdef OPT_BODY
         CALL BDYFLX(
     O                TQ,
     I                TB)
#endif
#ifdef OPT_PARALLEL
#ifdef OPT_TRIPOLE
      CALL SHIFT2(
     M              TAUX,   TAUY,
     D             NXDIM,  NYDIM,      1,
     I             -1.D0,      0,      0 )
#else
         CALL SHIFT2(
     M                 TAUX,   TAUY,
     D                NXDIM,  NYDIM,      1)
#endif
#else
         CALL STBCTU(
     M                 TAUX,   TAUY)
#endif
         CALL PREDCO(
     M                HA,   UBTA,   VBTA,      W,      R,
     M                UA,     VA,     TA,
     M               AMV,    AHV,
     I              TAUX,   TAUY,     FT,   PTOP,
#ifdef OPT_BODY
     I                TQ,
#endif
     I                UB,     VB,     TB,
     I                HB,   UBTB,   VBTB,
     W                GX,     GY,     XX,     YY,
     W               GXX,    GYY,
     W              UADV,   VADV,   WADV)
         CALL EXCNGO(
     M                UA,     VA,     TA,
     M                HA,   UBTA,   VBTA,
     M                UB,     VB,     TB,
     M                HB,   UBTB,   VBTB)

      END IF

c *** Output to file ***

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

c        CHEKIN routine for FTX/Y/Z & FSX/Y/Z
         CALL CHKFTX

      CALL CHKOUT(
     I            OFLOUT)

      RETURN
      END
