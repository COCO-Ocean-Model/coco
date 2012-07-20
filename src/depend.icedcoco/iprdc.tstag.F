      SUBROUTINE PREDCI(
     M                   AX,    HIX,    UIX,    VIX,    TIX,    HSX,
     O                   FT,     FS,   TAUX,   TAUY,   PTOP,
     I                   AY,    HIY,    UIY,    VIY,    TIY,    HSY,
     I                   TX,     UX,     VX,     HX,     HY,
     I                  QAO,    QAI,    QII,    QIO,  SWABS,
     I                  WEV,    WSB,
     I                 PREC,   SNOW,   ROFF,   SOFF,
     I               TAUAIX, TAUAIY, TAUAOX, TAUAOY)
      use zocdim
      use zocgrd
      use zocphy
      use ictrn
      use ifhea
      use ifwtr
      use ipadv
      use ipmmt
      use iprdg
      use iptmp
      use ucloc
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  HISTORY
c     '02.10.10  H.Hasumi: from MIROC3.1
c     '07.04.23  H.Hasumi
c     '07.09.25  H.Hasumi: multi-category sea ice
c     '07.09.28  H.Hasumi: 1-layer sea ice thermodynamics
c     '08.07.10  H.Hasumi: initial/final processing
c     '09.02.20  Y.Komuro: parallel forward time march
c     '09.05.25  Y.Komuro: CMIP5 output code included
c     '09.09.04  Y.Komuro: extra output code (FEX/FEY)
c     '09.09.26  Y.Komuro: bug fix
c     '10.04.14  M.Kurogi: staggered time stepping
c     '10.04.14  M.Kurogi: (COCO4.4 tripolar code by Dr. Suzuki)
c ---------------------------------------------------------------------

      REAL*8      TX(NXYDIM, NZDIM, NTDIM)
      REAL*8      UX(NXYDIM, NZDIM),     VX(NXYDIM, NZDIM)
      REAL*8      HX(NXYDIM),     HY(NXYDIM)
      REAL*8    PTOP(NXYDIM)

      REAL*8      AY(NXYDIM, 0:NIC),     AX(NXYDIM, 0:NIC)
      REAL*8     HIY(NXYDIM, 0:NIC),    HIX(NXYDIM, 0:NIC)
      REAL*8     UIY(NXYDIM),    UIX(NXYDIM)
      REAL*8     VIY(NXYDIM),    VIX(NXYDIM)
      REAL*8     TIY(NXYDIM, 0:NIC),    TIX(NXYDIM, 0:NIC)
      REAL*8     HSY(NXYDIM, 0:NIC),    HSX(NXYDIM, 0:NIC)

      REAL*8     QIO(NXYDIM, NIC),    QAI(NXYDIM, NIC)
      REAL*8     QII(NXYDIM, NIC)
      REAL*8     QAO(NXYDIM),  SWABS(NXYDIM)
      REAL*8     WEV(NXYDIM),    WSB(NXYDIM, NIC)
      REAL*8    PREC(NXYDIM),   SNOW(NXYDIM)
      REAL*8    ROFF(NXYDIM),   SOFF(NXYDIM)
      REAL*8  TAUAIX(NXYDIM), TAUAIY(NXYDIM)
      REAL*8  TAUAOX(NXYDIM), TAUAOY(NXYDIM)

      REAL*8      FT(NXYDIM, NTDIM),     FS(NXYDIM)
      REAL*8    TAUX(NXYDIM),   TAUY(NXYDIM)

      REAL*8     WAO(NXYDIM)
      REAL*8     WIO(NXYDIM, NIC),    WAS(NXYDIM, NIC)
      REAL*8     WIL(NXYDIM, NIC)
      REAL*8    SUBI(NXYDIM, NIC)
      REAL*8    EVAP(NXYDIM), ADJLAT(NXYDIM)
      REAL*8      AZ(NXYDIM, 0:NIC),    HIZ(NXYDIM, 0:NIC)
      REAL*8     HSZ(NXYDIM, 0:NIC),    TIZ(NXYDIM, 0:NIC)
      REAL*8    PICE(NXYDIM)

      REAL*8     EIX(NXYDIM, 0:NIC),    EIZ(NXYDIM, 0:NIC)

      REAL*8     FIX(NXYDIM, 0:NIC),    FIY(NXYDIM, 0:NIC)
      REAL*8     FSX(NXYDIM, 0:NIC),    FSY(NXYDIM, 0:NIC)
      REAL*8     FEX(NXYDIM, 0:NIC),    FEY(NXYDIM, 0:NIC)
      REAL*8   FTITD(NXYDIM)
      REAL*8  IGRFRA(NXYDIM), IGRCON(NXYDIM), IGRSNI(NXYDIM)
      REAL*8  INRLAT(NXYDIM)
      REAL*8  IMRSNO(NXYDIM), IMRISF(NXYDIM), IMRIBS(NXYDIM)

      LOGICAL OEOF

cc     for check
c      REAL*8  IMRTOT(NXYDIM), IGRTOT(NXYDIM), ITRTOT(NXYDIM)

      SAVE

      INTEGER     IJ,      L


      CALL CLCSTR('ICE')

#ifdef OPT_PARALLEL
      IF (      (MYRANK .GE. IJNODE)
     &    .AND. (.NOT. OINIT) .AND. (.NOT. OFINAL)) RETURN
#endif

      IF(OINIT) THEN
#ifdef OPT_TRIPOLE
         CALL RSTADD(PICE, OEOF, NXDIM, NYDIM, 1, 'PICE', 'SFC',
     &                                             1.D0,  0,  0 )
#else
         CALL RSTADD(PICE, OEOF, NXDIM, NYDIM, 1, 'PICE', 'SFC')
#endif
         IF(OEOF) THEN
            DO L = 0, NIC
            DO IJ = 1, NXYDIM
               AZ  (IJ, L) =  AX  (IJ, L)
               HIZ (IJ, L) =  HIX (IJ, L)
               EIZ (IJ, L) =  EIX (IJ, L)
               HSZ (IJ, L) =  HSX (IJ, L)
               TIZ (IJ, L) =  TIX (IJ, L) 
            END DO
            END DO
            CALL PRIDGE(
     O              PICE,
     M                AX,    HIX,    EIX,    HSX,    TIX,
     I                AZ,    HIZ,    EIZ,    HSZ,
     I               UIX,    VIX)
            DO L = 0, NIC
            DO IJ = 1, NXYDIM
               AX  (IJ, L) =  AZ  (IJ, L)
               HIX (IJ, L) =  HIZ (IJ, L)
               EIX (IJ, L) =  EIZ (IJ, L)
               HSX (IJ, L) =  HSZ (IJ, L)
               TIX (IJ, L) =  TIZ (IJ, L) 
            END DO
            END DO
         END IF
      END IF
      
      IF (OFINAL) THEN
         CALL FINADD(PICE, NXDIM, NYDIM, 1, 'PICE', 'SFC')
      END IF

      DO IJ = 1, NXYDIM
         EVAP(IJ) = WEV(IJ)
         FS  (IJ) = 0.D0
         TAUX(IJ) = 0.D0
         TAUY(IJ) = 0.D0
         AZ  (IJ, 0) = 1.D0
         HIZ (IJ, 0) = 0.D0
         HSZ (IJ, 0) = 0.D0
      END DO
      DO L = 1, NIC
         DO IJ = 1, NXYDIM
            SUBI(IJ, L) = 0.D0
            AZ  (IJ, L) = 0.D0
            HIZ (IJ, L) = 0.D0
            HSZ (IJ, L) = 0.D0
         END DO
      END DO
      DO L = 1, 2
         DO IJ = 1, NXYDIM
            FT(IJ, L) = 0.D0
         END DO
      END DO


      CALL CLCSTR('ICEDYN')
      CALL PMOMNT(
     M               UIX,    VIX,
     O              TAUX,   TAUY,
     I                AX,     AY,     AZ,
     I               HIX,    HIY,    HIZ,
     I               HSX,    HSY,    HSZ,
     I               UIY,    VIY,   PICE,
     I                UX,     VX,     HY,   PTOP,
     I            TAUAIX, TAUAIY, TAUAOX, TAUAOY)
#ifdef OPT_PARALLEL
#ifdef OPT_TRIPOLE
      CALL SHIFT2(
     M              TAUX,    TAUY, 
     D             NXDIM,   NYDIM,     1,
     I             -1.D0,      -1,    -1 )
#else
      CALL SHIFT2(
     M              TAUX,   TAUY,
     I             NXDIM,  NYDIM,      1)
#endif
#else
      CALL STBCTU(
     M               TAUX,   TAUY)
#endif
      CALL CLCEND('ICEDYN')




      CALL ICETMP(
     O               EIX,
     M               TIX,
     M               QAO,    QAI,    QIO,    QII,
     I                AX,    HIX)

      CALL FIHEAT(
     M                FT,
     O               WAO,    WIO,    WAS,    WIL,
     I                AX,     TX,     HX,
     I               QAO,    QAI,    QIO,    QII,  SWABS)
      CALL FWATER(
     M                AX,    HIX,    HSX,
     M              PREC,   SNOW,
     O              EVAP,   SUBI, ADJLAT,
     I               WEV,    WSB,   SOFF)
      CALL PTHERM(
     M                AX,    HIX,    HSX,
     M               EIX,    TIX,
     M              PREC,   SNOW,     FT,     FS,
     O             FTITD, IGRFRA, IGRCON, IGRSNI,
     O            INRLAT, IMRSNO, IMRISF, IMRIBS,
     I                TX,
     I               WIO,    WAO,    WAS,    WIL,
     I              EVAP,   SUBI,   ROFF, ADJLAT,
     I               QIO)
      CALL ICTRNS(
     M                AX,    HIX,    HSX,    EIX,    TIX,
     M                FT,     FS)

#ifdef OPT_PARALLEL
#ifdef OPT_TRIPOLE
      CALL SHIFT3(
     M                AX,    HIX,    HSX,
     D             NXDIM,  NYDIM,  NIC+1,
     I             1.0D0,      0,      0 )
      CALL SHIFT2(
     M               EIX,    TIX,
     D             NXDIM,  NYDIM,  NIC+1,
     I             1.0D0,      0,      0 )
#else
      CALL SHIFT3(
     M                AX,    HIX,    HSX,
     D             NXDIM,  NYDIM,  NIC+1)
      CALL SHIFT2(
     M               EIX,    TIX,
     D             NXDIM,  NYDIM,  NIC+1)
#endif
#else
      CALL STBCSI(
     M                AX,    HIX,    EIX,    TIX,    HSX)
#endif

      CALL PADVCT(
     M                AX,    HIX,    EIX,    HSX,    TIX,
     O                AZ,    HIZ,    EIZ,    HSZ,
     O               FIX,    FIY,    FSX,    FSY,
     O               FEX,    FEY,
     I               UIX,    VIX)
      CALL PRIDGE(
     O              PICE,
     M                AX,    HIX,    EIX,    HSX,    TIX,
     I                AZ,    HIZ,    EIZ,    HSZ,
     I               UIX,    VIX)
      CALL ICADJS(
     M                AX,    HIX,    HSX,    EIX)
      CALL ICTRNS(
     M                AX,    HIX,    HSX,    EIX,    TIX,
     M                FT,     FS)
#ifdef OPT_PARALLEL
#ifdef OPT_TRIPOLE
      CALL SHIFT3(
     M                AX,    HIX,    HSX,
     D             NXDIM,  NYDIM,  NIC+1,
     I             1.0D0,      0,      0 )
      CALL SHIFT2(
     M               EIX,    TIX,
     D             NXDIM,  NYDIM,  NIC+1,
     I             1.0D0,      0,      0 )
#else
      CALL SHIFT3(
     M                AX,    HIX,    HSX,
     D             NXDIM,  NYDIM,  NIC+1)
      CALL SHIFT2(
     M               EIX,    TIX,
     D             NXDIM,  NYDIM,  NIC+1)
#endif
#else
      CALL STBCSI(
     M                AX,    HIX,    EIX,    TIX,    HSX)
#endif


#ifdef OPT_PARALLEL
      IF (MYRANK .GE. IJNODE) RETURN
#endif

      DO L = 1, NIC
         DO IJ = 1, NXYDIM
            PTOP(IJ) = PTOP(IJ)
     &               + gravit * AY(IJ, L) *
     &                 (rhoi * HIY(IJ, L) + rhos * HSY(IJ, L))
         END DO
      END DO

      CALL CLCEND('ICE')

c     output section for CMIP5
c     FIX, FIY: eastward/northward sea ice transport
c               unit [cm^3(ice)/s]
      CALL CHEKIN(   FIX,  'FIX',
     &                NX,     NY,    NIC, NXYIDM, 'ICE')
      CALL CHEKIN(   FIY,  'FIY',
     &                NX,     NY,    NIC, NXYIDM, 'ICE')
c     FISX, FISY: eastward/northward snow transport
c                 unit [cm^3(snow)/s]
      CALL CHEKIN(   FSX, 'FISX',
     &                NX,     NY,    NIC, NXYIDM, 'ICE')
      CALL CHEKIN(   FSY, 'FISY',
     &                NX,     NY,    NIC, NXYIDM, 'ICE')
c     FTITD: heat flux into sea water due to sea ice thermodynamics
c            unit [erg/cm^2/s]
      CALL CHEKIN(  FTITD,  'FTITD',
     &                NX,     NY,      1, NXYDIM, 'SFC')
c     IGRFRA: frazil sea ice growth rate, unit [g/cm^2/s]
      CALL CHEKIN( IGRFRA, 'IGRFRA',
     &                NX,     NY,      1, NXYDIM, 'SFC')
c     IGRCON: congelation sea ice growth rate, unit [g/cm^2/s]
      CALL CHEKIN( IGRCON, 'IGRCON',
     &                NX,     NY,      1, NXYDIM, 'SFC')
c     IGRSNI: snow-ice formation rate, unit [g/cm^2/s]
      CALL CHEKIN( IGRSNI, 'IGRSNI',
     &                NX,     NY,      1, NXYDIM, 'SFC')
c     INRLAT: lateral sea ice net growth rate, unit [g/cm^2/s]
c      (COCO4.4 represents lateral melting process but not freezing,
c       thus this value will be zero or negative.)
      CALL CHEKIN( INRLAT, 'INRLAT',
     &                NX,     NY,      1, NXYDIM, 'SFC')
c     IMRSNO: snow melt rate, unit [g/cm^2/s]
      CALL CHEKIN( IMRSNO, 'IMRSNO',
     &                NX,     NY,      1, NXYDIM, 'SFC')
c     IMRISF: rate of melt at upper surface of sea ice, unit [g/cm^2/s]
      CALL CHEKIN( IMRISF, 'IMRISF',
     &                NX,     NY,      1, NXYDIM, 'SFC')
c     IMRIBS: rate of melt at sea ice base, unit [g/cm^2/s]
      CALL CHEKIN( IMRIBS, 'IMRIBS',
     &                NX,     NY,      1, NXYDIM, 'SFC')

c     extra output section
c     FEX, FEY: eastward/northward thermal transport by sea ice
c               unit [(erg/g)*(cm^3(ice)/s)]
      CALL CHEKIN(   FEX,  'FIEX',
     &                NX,     NY,    NIC, NXYIDM, 'ICE')
      CALL CHEKIN(   FEY,  'FIEY',
     &                NX,     NY,    NIC, NXYIDM, 'ICE')

cc     for check: not necessary for CMIP5 output
c      DO IJ = 1, NXYDIM
c         IGRTOT(IJ) = IGRFRA(IJ)+IGRCON(IJ)+IGRSNI(IJ)
c         IMRTOT(IJ) = IMRSNO(IJ)+IMRISF(IJ)+IMRIBS(IJ)-INRLAT(IJ)
c         ITRTOT(IJ) = IGRTOT(IJ)-IMRTOT(IJ)
c      END DO
c      CALL CHEKIN( IGRTOT, 'IGRTOT',
c     &                NX,     NY,      1, NXYDIM, 'SFC')
c      CALL CHEKIN( IMRTOT, 'IMRTOT',
c     &                NX,     NY,      1, NXYDIM, 'SFC')
c      CALL CHEKIN( ITRTOT, 'ITRTOT',
c     &                NX,     NY,      1, NXYDIM, 'SFC')

      RETURN
      END
