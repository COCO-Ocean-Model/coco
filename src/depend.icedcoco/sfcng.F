      SUBROUTINE SFCFLX(
     O               QAO,    QAI,    QII,    QIO,  SWABS,    TSI,
     O               WEV,    WSB,
     O              PREC,   SNOW,   ROFF,   SOFF,
     O            TAUAIX, TAUAIY, TAUAOX, TAUAOY,
     O                FT,   PTOP,   SSFC,
     I                 T,      A,     HI,     TI,    HSN)
      use zocdim
      use zocgrd
      use zocmsk
      use zocphy
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Estimate the surface forcing fluxes from the data given in the
c boundary condition files
c
c  HISTORY
c     '02.10.10  H.Hasumi: from COCO3.4 and MIROC3.1
c     '07.04.23  H.Hasumi
c     '07.09.25  H.Hasumi: multi-category sea ice
c     '07.09.26  H.Hasumi: sea ice (and snow) melting temperature
c                          is defined independently from TMELT
c     '07.09.28  H.Hasumi: 1-layer sea ice thermodynamics
c     '09.05.25  Y.Komuro: CMIP5 output code included
c     '09.09.26  Y.Komuro: bug fix
c     '10.04.14  M.Kurogi: (COCO4.4 tripolar code by Dr. Suzuki)
c
c
c ---------------------------------------------------------------------

      REAL*8 FACTM, FACTH, FACTW
      PARAMETER(FACTM = 1.0D+1, FACTH = 1.0D+3, FACTW = 1.0D+2)

      INTEGER NTYZ0
      PARAMETER(NTYZ0 = 3)

      REAL*8       T(NXYDIM, NZDIM, NTDIM)
      REAL*8       A(NXYDIM, 0:NIC)
      REAL*8      HI(NXYDIM, 0:NIC),    HSN(NXYDIM, 0:NIC)
      REAL*8      TI(NXYDIM, 0:NIC)

      REAL*8     WEV(NXYDIM),    WSB(NXYDIM, NIC)
      REAL*8    PREC(NXYDIM),   SNOW(NXYDIM)
      REAL*8    ROFF(NXYDIM),   SOFF(NXYDIM)
      REAL*8  TAUAIX(NXYDIM), TAUAIY(NXYDIM)
      REAL*8  TAUAOX(NXYDIM), TAUAOY(NXYDIM)
      REAL*8     QAI(NXYDIM, NIC),    QIO(NXYDIM, NIC)
      REAL*8     QII(NXYDIM, NIC)
      REAL*8     QAO(NXYDIM),  SWABS(NXYDIM)
      REAL*8     TSI(NXYDIM, 0:NIC)
      REAL*8      FT(NXYDIM, NTDIM)
      REAL*8    PTOP(NXYDIM)
      REAL*8    SSFC(NXYDIM)

      REAL*8    TAUX(NXYDIM),   TAUY(NXYDIM),   USFC(NXYDIM)
      REAL*8    TSFC(NXYDIM),   QSFC(NXYDIM),   WFLX(NXYDIM)
      REAL*8    SWNT(NXYDIM),   DWLW(NXYDIM),   PSFC(NXYDIM)
      REAL*8    GRTS(NXYDIM),   GRTB(NXYDIM)
      REAL*8   GRICE(NXYDIM),  GRSNW(NXYDIM),  GRICR(NXYDIM)
      REAL*8    GRZ0(NXYDIM, NTYZ0)
      REAL*8  GFLUXS(NXYDIM), TFLUXS(NXYDIM), QFLUXS(NXYDIM)
      REAL*8  WFLUXS(NXYDIM, 2)
      REAL*8  RFLXLU(NXYDIM), SFLXBL(NXYDIM)
      REAL*8   DGFDS(NXYDIM),  DTFDT(NXYDIM),  DTFDS(NXYDIM)
      REAL*8   DQFDS(NXYDIM)

      REAL*8     TMI
      SAVE TMI

      REAL*8   TISI(NXYDIM, 0:NIC)
      REAL*8     DIRDSN
      SAVE DIRDSN

*     [INTERNAL PARM] 
*     We need DFSNOW and DFICE for flux redistribution.
      INTEGER NRBND
      PARAMETER(NRBND = 3)

      REAL*8     DFSNOW                      !! diffusion coef. of snow
      REAL*8     SNWDMX                      !! maximum snow depth
      REAL*8     EPSSNW                      !! minimum snow
      REAL*8     SNWMAX                      !! ice forming snow
      REAL*8     ALBSNW( 2, NRBND )          !! snow albedo
      REAL*8     TALSNW( 2 )                 !! temp. for alb. change
      REAL*8     Z0SNW ( NTYZ0 )             !! roughness of snow
      REAL*8     SNWCRT                      !! snow amount for fraction=1
      REAL*8     SNWDEN                      !! snow density (kg/m**3)
*
      REAL*8     DFICE                       !! diffusion coef. of sea ice
      REAL*8     ALBICE( NRBND )             !! sea ice albedo
      REAL*8     Z0ICE ( NTYZ0 )             !! roughness of sea ice
      REAL*8     SICCRT                      !! ice amount for conc.=1
      REAL*8     SICDEN                      !! ice density (kg/m**3)
*
      SAVE  DFSNOW, SNWDMX, EPSSNW, SNWMAX,
     &      ALBSNW, TALSNW, Z0SNW ,
     &      SNWCRT, SNWDEN,  
     &       DFICE,  ALBICE, Z0ICE ,
     &      SICCRT, SICDEN
*
      NAMELIST  /NMSNOW/ DFSNOW, SNWDMX, EPSSNW, SNWMAX,
     &                   ALBSNW, TALSNW, Z0SNW ,
     &                   SNWCRT, SNWDEN  
      NAMELIST  /NMICE/  DFICE,  ALBICE, Z0ICE ,
     &                   SICCRT, SICDEN
*
      INTEGER     IJ,      L,      N,     NN
      INTEGER  IFPAR,  JFPAR
      LOGICAL OFIRST
      SAVE OFIRST

      REAL*8      SI
      NAMELIST /NMISLT/ SI
      SAVE SI

      DATA SI / 5.D0 /
      DATA OFIRST / .TRUE. /
*
      DATA       DFSNOW         / 0.4D0          /
      DATA       SNWDMX         / 5.D0           /
      DATA       EPSSNW         / 1.D-8          /        
      DATA       SNWMAX         / 1000.D0        /
      DATA       ALBSNW         / 0.75D0, 0.5D0,
     &                            0.75D0, 0.5D0, 
     &                            0.0D0 , 0.0D0      /
      DATA       TALSNW         / 258.15D0, 273.15D0 /
      DATA       Z0SNW          / 1.D-2, 1.D-3, 1.D-3 /
      DATA       SNWCRT         / 100.D0         /
      DATA       SNWDEN         / 400.D0         /
*
      DATA       DFICE          / 2.D0           /
      DATA       ALBICE         / 0.5D0, 0.5D0, 0.05D0 /
      DATA       Z0ICE          / 2.D-2, 2.D-3, 2.D-3 /
      DATA       SICCRT         /  300.D0        /
      DATA       SICDEN         / 1000.D0        /

      IF (OFIRST) THEN
         CALL REWNML(IFPAR, JFPAR)
         WRITE(JFPAR, *) '*** SFCFLX ***'
         READ(IFPAR, NMISLT, END=99)
 99      WRITE(JFPAR, NMISLT)
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMSNOW, END=98)
 98      WRITE(JFPAR, NMSNOW)
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMICE, END=97)
 97      WRITE(JFPAR, NMICE)
         TMI = dtds * SI         
         DIRDSN = DFICE / DFSNOW
         DO IJ = 1, NXYDIM
            TSFC(IJ) = 300.D0
            PSFC(IJ) = 1.D5
         END DO
         OFIRST = .FALSE.
      END IF

      DO L = 0, NIC
         DO IJ = 1, NXYDIM
            TISI(IJ, L) = 0.0D0
         END DO
      END DO

      DO L = 3, NTDIM
         N = 11 + (L - 3) * 2
         NN = N + 1
         CALL TMINTP(  TSFC,      N)
         CALL TMINTP(  QSFC,     NN)
         DO IJ = IJSTR, IJEND
            FT(IJ, L) = QSFC(IJ) * (TSFC(IJ) - T(IJ, KSTR, L)) *
     &                  DZ(IJ, KSTR) * AMSKT(IJ, KSTR)
         END DO
      END DO
      CALL TMINTP(  TAUX,      1)
      CALL TMINTP(  TAUY,      2)
      CALL TMINTP(  USFC,      3)
      CALL TMINTP(  TSFC,      4)
      CALL TMINTP(  QSFC,      5)
      CALL TMINTP(  WFLX,      6)
      CALL TMINTP(  SWNT,      7)
      CALL TMINTP(  DWLW,      8)
      CALL TMINTP(  PSFC,      9)
#ifdef OPT_SRST
      CALL TMINTP(  SSFC,     10)
#endif

#ifdef OPT_PARALLEL
#ifdef OPT_TRIPOLE
      CALL SHIFT2(
     M              TAUX,   TAUY,
     D             NXDIM,  NYDIM,      1,
     I             -1.D0,      0,      0 )
      CALL SHIFT1(
     M             PSFC,
     D             NXDIM,  NYDIM,      1,
     I             1.D0,      0,      0 )
#else
      CALL SHIFT3(
     M              TAUX,   TAUY,   PSFC,
     D             NXDIM,  NYDIM,      1)
#endif
#else
      CALL STBCTU(
     M              TAUX,   TAUY)
      CALL STBCSH(
     M              PSFC)
#endif

      DO IJ = IJSTR, IJEND
         TAUAIX(IJ) = (  TAUX(IJ) + TAUX(IJ+LE)
     &                 + TAUX(IJ+LN) + TAUX(IJ+LNE)) *
     &                0.25D0 * FACTM
         TAUAIY(IJ) = (  TAUY(IJ) + TAUY(IJ+LE)
     &                 + TAUY(IJ+LN) + TAUY(IJ+LNE)) *
     &                0.25D0 * FACTM
         TAUAOX(IJ) = (  TAUX(IJ) + TAUX(IJ+LE)
     &                 + TAUX(IJ+LN) + TAUX(IJ+LNE)) *
     &                0.25D0 * FACTM
         TAUAOY(IJ) = (  TAUY(IJ) + TAUY(IJ+LE)
     &                 + TAUY(IJ+LN) + TAUY(IJ+LNE)) *
     &                0.25D0 * FACTM
         SNOW(IJ) = 0.D0
         PREC(IJ) = - WFLX(IJ) * FACTW
         SOFF(IJ) = 0.D0
         ROFF(IJ) = 0.D0
         SWABS(IJ) = SWNT(IJ) * FACTH * A(IJ, 0)
      END DO

      DO L = 0, NIC

         IF (L .GT. 0) THEN
            DO IJ = IJSTR, IJEND
               GRTS (IJ) = TSI(IJ, L) + kelvin
               GRTB (IJ) = TI(IJ, L) + kelvin
               GRICE(IJ) = HI(IJ, L) * 5.D-3
               GRSNW(IJ) = HSN(IJ, L) * 1.D-2
               GRICR(IJ) = 1.D0
            END DO
         ELSE
            DO IJ = IJSTR, IJEND
               GRTS (IJ) = T(IJ, KSTR, 1) + kelvin
               GRTB (IJ) = T(IJ, KSTR, 1) + kelvin
               GRICE(IJ) = 0.D0
               GRSNW(IJ) = 0.D0
               GRICR(IJ) = 0.D0
            END DO
         END IF

         CALL OCNBCS_OMIP(
     O              GRZ0, GFLUXS,  DGFDS,
     I              GRTS,   GRTB,  GRICE,  GRSNW,  GRICR,
     I              USFC)
         CALL SFCFLX_OMIP(
     O            TFLUXS, QFLUXS,
     O             DTFDT,  DTFDS,  DQFDS,
     I              USFC,   TSFC,   QSFC,   PSFC,
     I              GRTS,   GRZ0)
         CALL OCNSLV_OMIP(
     M              GRTS, GFLUXS, TFLUXS, QFLUXS,
     O            WFLUXS, RFLXLU, SFLXBL,
     I             DTFDS,  DQFDS,  DGFDS,
     I              SWNT,   DWLW,
     I             GRICR,    TMI)

         IF (L .GT. 0) THEN
            DO IJ = IJSTR, IJEND
               QAI(IJ, L) = (GFLUXS(IJ) - SFLXBL(IJ)) * FACTH
               QII(IJ, L) = GFLUXS(IJ) * FACTH
               QIO(IJ, L) = (T(IJ, KSTR, 1) - TI(IJ, L)) *
     &                      cdi / HI(IJ, L) * 2.D0
               TSI(IJ, L) = MIN(GRTS(IJ), TMELT+TMI) - kelvin
               TISI(IJ, L) = TSI(IJ, L) +
     &                       (TI(IJ, L) - TSI(IJ, L)) *
     &                       DIRDSN * HSN(IJ, L) /
     &                       (0.5D0 * HI(IJ, L) + DIRDSN * HSN(IJ, L))
               WSB(IJ, L) = 0.D0 ! (sublimation speed) x (fraction)
            END DO
         ELSE
            DO IJ = IJSTR, IJEND
               QAO(IJ) = (GFLUXS(IJ) + SWNT(IJ)) * FACTH
               WEV(IJ) = 0.D0 ! (evaporation speed) x (fraction)
            END DO
         END IF
      END DO

c     output section for CMIP5
c     TISI: temperature at interface between sea ice and snow
c           unit [C]
      CALL CHEKIN(TISI,   'TISI',
     &            NX,     NY,    NIC, NXYIDM, 'ICE')

      DO IJ = 1, NXYDIM
         PTOP(IJ) = PSFC(IJ) * FACTM
      END DO

      RETURN
      END

c *********************************************************************

      SUBROUTINE OCNSLV_OMIP
     M         ( GDTS  , GFLUXS, TFLUXS, QFLUXS,
     O           WFLUXS, RFLXLU, SFLXBL,
     I           DTFDS , DQFDS , DGFDS ,
     I           RFLXSD, RFLXLD,
     I           GRICR , TMI) 
      use zocdim
      use zocphy
      IMPLICIT NONE
*   
*   [MODIFY] 
      REAL*8     GDTS  ( NXYDIM )                !! skin temperature
      REAL*8     GFLUXS( NXYDIM )                !! soil heat flux
      REAL*8     TFLUXS( NXYDIM )                !! flux of T
      REAL*8     QFLUXS( NXYDIM )                !! flux of q
*
*   [OUTPUT]
      REAL*8     WFLUXS( NXYDIM, 2 )             !! soil water flux
      REAL*8     RFLXLU( NXYDIM )                !! upward long wave
      REAL*8     SFLXBL( NXYDIM )                !! flux balance
*
*   [INPUT] 
      REAL*8     DTFDS ( NXYDIM )                !! dH/dTg
      REAL*8     DQFDS ( NXYDIM )                !! dE/dTg
      REAL*8     DGFDS ( NXYDIM )                !! dG/dTg
      REAL*8     RFLXSD( NXYDIM )                !! down. SW rad.
      REAL*8     RFLXLD( NXYDIM )                !! down. LW rad.
      REAL*8     GRICR ( NXYDIM )                !! snow/ice ratio
      REAL*8     TMI ! sea ice (and snow) melting temperature in C
*
*   [INTERNAL WORK] 
      REAL*8     ESUB, STG, DRFDS
      REAL*8     SFLUX, GSFLUX, DGSFDS
      REAL*8     GFLUXF
      REAL*8     SFLXBI, DSBDSI, DTI, EVAPI, GFLUXI
      REAL*8     FF, FI, DTX
      INTEGER    IJ
      INTEGER    IFPAR, JFPAR
*
*   [INTERNAL SAVE] 
      LOGICAL OFIRST
      LOGICAL OLWNET
      REAL*8    EMIS !! longwave emissivity (longwave co-albedo)
      REAL*8 FLWNET
      NAMELIST /NMLWEM/ EMIS, OLWNET
      SAVE OFIRST, EMIS, FLWNET
      DATA EMIS / 1.D0 /
      DATA OLWNET / .FALSE. /
      DATA       OFIRST / .TRUE.  /
*
      IF ( OFIRST ) THEN
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMLWEM, END=110)
 110     WRITE(JFPAR, NMLWEM)
         OFIRST = .FALSE.
         WRITE (JFPAR, *) ' @@@ OCNSLV: OCEAN HEAT BALANCE 98/06/30'
         IF (OLWNET) THEN
            FLWNET = 0.D0
         ELSE
            FLWNET = 1.D0
         END IF
      ENDIF
*
      ESUB = EL + EMELT
      DO 1000 IJ = IJSTR, IJEND
         STG            = EMIS*STB*GDTS( IJ )**4
         RFLXLU( IJ )   = (STG + ( 1.D0-EMIS )*RFLXLD( IJ )) * FLWNET
         DRFDS          = 4.D0*STG/GDTS( IJ ) * FLWNET
*
         SFLUX          = TFLUXS( IJ )
     &                  + RFLXLU( IJ ) - RFLXLD( IJ ) 
     &                  - RFLXSD( IJ )
         GSFLUX         = GFLUXS( IJ ) - SFLUX
         DGSFDS         = DGFDS ( IJ )
     &                  + DTFDS ( IJ ) + DRFDS
*
*        < ice free ocean surface >
         GFLUXF         = SFLUX + QFLUXS( IJ )*EL
*         
*        < ice/snow surface >
         SFLXBI         = GSFLUX - QFLUXS( IJ )*ESUB
         DSBDSI         = DGSFDS + DQFDS ( IJ )*ESUB
         DTI            = SFLXBI/DSBDSI
         DTI            = MIN( DTI, TMELT+TMI-GDTS( IJ ) )
         EVAPI          = QFLUXS( IJ ) + DQFDS( IJ )*DTI
         GFLUXI         = GFLUXS( IJ ) - DGFDS( IJ )*DTI
         SFLXBI         = SFLXBI       - DSBDSI     *DTI
*
*        < Ts change >
         FF             = 1.D0 - GRICR( IJ )
         FI             = GRICR( IJ )
         DTX            = FI * DTI
         GDTS  ( IJ )   = GDTS( IJ ) + DTX
*
*        < fluxs >
         TFLUXS( IJ )   = TFLUXS( IJ ) + DTFDS( IJ )*DTX
         RFLXLU( IJ )   = RFLXLU( IJ ) + DRFDS      *DTX
         QFLUXS( IJ )   = FF * QFLUXS( IJ ) + FI * EVAPI
         GFLUXS( IJ )   = FF * GFLUXF       + FI * GFLUXI 
         WFLUXS( IJ,1 ) = FF * QFLUXS( IJ )
         WFLUXS( IJ,2 ) = FI * EVAPI
         SFLXBL( IJ )   = FI * SFLXBI
 1000 CONTINUE 
*

      RETURN
      END

c *********************************************************************

      SUBROUTINE OSFBCS_OMIP
      use zocdim
      use zocphy
      IMPLICIT NONE

      INTEGER NTYZ0, NRBND
      PARAMETER(NTYZ0 = 3, NRBND = 3)

*   [OUTPUT] 
      REAL*8     GRZ0  ( NXYDIM, NTYZ0  )       !! surface roughness 
      REAL*8     FOGFLX( NXYDIM )               !! heat flux
      REAL*8     DGFDS ( NXYDIM )               !! dG/dTs
*   [INPUT] 
      REAL*8     GRTS  ( NXYDIM )               !! skin temperature
      REAL*8     GRTB  ( NXYDIM )               !! ice base temp.
      REAL*8     GRICE ( NXYDIM )               !! sea ice
      REAL*8     GRSNW ( NXYDIM )               !! snow smount
      REAL*8     GRICR ( NXYDIM )               !! ice fraction
      REAL*8     USFC  ( NXYDIM )               !! sfc wind
*   [INTERNAL WORK] 
      REAL*8     GRSNR ( NXYDIM )               !! snow fraction
*
*   [INTERNAL WORK] 
      INTEGER    IJ, L, M
      INTEGER    IFG
      REAL*8     TALSNX, ALBSNX, ALB0, DALB, TFACT, ALBX
      REAL*8     Z00, DZ0
      REAL*8     DFGT, DFGX
      INTEGER    IFPAR, JFPAR
*
*   [INTERNAL SAVE] 
      LOGICAL    OFIRST
      SAVE       OFIRST
*
*   [INTERNAL PARM] 
      INTEGER    NTYZZ
      PARAMETER  (NTYZZ =NTYZ0-1)
      REAL*8     DFSNOW                      !! diffusion coef. of snow
      REAL*8     SNWDMX                      !! maximum snow depth
      REAL*8     EPSSNW                      !! minimum snow
      REAL*8     SNWMAX                      !! ice forming snow
      REAL*8     ALBSNW( 2, NRBND )          !! snow albedo
      REAL*8     TALSNW( 2 )                 !! temp. for alb. change
      REAL*8     Z0SNW ( NTYZ0 )             !! roughness of snow
      REAL*8     SNWCRT                      !! snow amount for fraction=1
      REAL*8     SNWDEN                      !! snow density (kg/m**3)
*
      REAL*8     DFICE                       !! diffusion coef. of sea ice
      REAL*8     ALBICE( NRBND )             !! sea ice albedo
      REAL*8     Z0ICE ( NTYZ0 )             !! roughness of sea ice
      REAL*8     SICCRT                      !! ice amount for conc.=1
      REAL*8     SICDEN                      !! ice density (kg/m**3)
*
      REAL*8     Z0FCT                       !! heat z0/moumentum z0
      REAL*8     Z0MIN                       !! minimum z0
*
      REAL*8     DZOCN                       !! depth of ML Ocean
      REAL*8     DFOCN                       !! ocean dG/dTs
      REAL*8     ALBLO                       !! LW albedo (1-emis)
*
      NAMELIST  /NMSNOW/ DFSNOW, SNWDMX, EPSSNW, SNWMAX,
     &                   ALBSNW, TALSNW, Z0SNW ,
     &                   SNWCRT, SNWDEN  
      NAMELIST  /NMICE/  DFICE,  ALBICE, Z0ICE ,
     &                   SICCRT, SICDEN
      NAMELIST  /NMZ0/   Z0FCT,  Z0MIN
      NAMELIST  /NMOCN/  DZOCN,  DFOCN , ALBLO
*
      DATA       DFSNOW         / 0.4D0          /
      DATA       SNWDMX         / 5.D0           /
      DATA       EPSSNW         / 1.D-8          /        
      DATA       SNWMAX         / 1000.D0        /
      DATA       ALBSNW         / 0.75, 0.5,
     &                            0.75, 0.5, 
     &                            0.0 , 0.0      /
      DATA       TALSNW         / 258.15, 273.15 /
      DATA       Z0SNW          / 1.D-2, 1.D-3, 1.D-3 /
      DATA       SNWCRT         / 100.D0         /
      DATA       SNWDEN         / 400.D0         /
*
      DATA       DFICE          / 2.D0           /
      DATA       ALBICE         / 0.5, 0.5, 0.05 /
      DATA       Z0ICE          / 2.D-2, 2.D-3, 2.D-3 /
      DATA       SICCRT         /  300.D0        /
      DATA       SICDEN         / 1000.D0        /
*
      DATA       Z0FCT          / 0.1D0          /
      DATA       Z0MIN          / 1.D-6          /
*
      DATA       DZOCN          / 50.D0          /
      DATA       DFOCN          / 1.D10          /
      DATA       ALBLO          / 0.05           /
*
      DATA       OFIRST         / .TRUE.         /

c =====================================================================

      ENTRY      OCNBCS_OMIP
     O         ( GRZ0  , FOGFLX, DGFDS ,
     I           GRTS  , GRTB  , GRICE , GRSNW , GRICR ,
     I           USFC  )

      IF ( OFIRST ) THEN
         CALL REWNML(IFPAR, JFPAR)
         WRITE (JFPAR, *) ' @@@ OCNBCS: OCEAN SURFACE BC 98/07/29'
         OFIRST = .FALSE.
         READ(IFPAR, NMSNOW, END=130)
 130     WRITE(JFPAR, NMSNOW)
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMICE, END=140)
 140     WRITE(JFPAR, NMICE)
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMZ0, END=150)
 150     WRITE(JFPAR, NMZ0)
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMOCN, END=160)
 160     WRITE(JFPAR, NMOCN)
      ENDIF

      DO 600 IJ = IJSTR, IJEND
         IF (GRSNW(IJ) .GT. 0.D0) THEN
            GRSNR(IJ) = 1.D0
         ELSE
            GRSNR(IJ) = 0.D0
         END IF
  600 CONTINUE 

      CALL SEAZ0F_OMIP
     O         ( GRZ0(1,1) , GRZ0(1,2), GRZ0(1,3),
     I           USFC      )
      DO 2100 M = 1, NTYZ0
         DO 2100 IJ = IJSTR, IJEND
            IF ( GRICR( IJ ) .GT. 0.D0 ) THEN
               DZ0   = Z0ICE( M )   - GRZ0( IJ,M )
               Z00   = GRZ0( IJ,M ) + DZ0*GRICR( IJ )
               DZ0   = Z0SNW( M )   - Z00
               GRZ0( IJ,M ) = Z00   + DZ0*GRSNR( IJ )
            ENDIF
            GRZ0( IJ,M ) = MAX( GRZ0( IJ,M ), Z0MIN )
 2100 CONTINUE 

      DO 3100 IJ = IJSTR, IJEND
         IF (GRICR(IJ) .GT. 0.D0) THEN
            DFGT = DFICE / GRICE(IJ)
            DFGX = DFICE * DFSNOW
     &             / (DFICE * GRSNW(IJ) + DFSNOW * GRICE(IJ))
            DFGT = DFGT * (1.D0 - GRSNR(IJ)) + DFGX * GRSNR(IJ)
            FOGFLX(IJ) = DFGT * (GRTB(IJ) - GRTS(IJ))
            DGFDS (IJ) = DFGT
         ELSE
            FOGFLX(IJ) = 0.D0
            DGFDS (IJ) = DFOCN
         ENDIF
 3100 CONTINUE

      RETURN
      END

c *********************************************************************

      SUBROUTINE SEAZ0F_OMIP
     O         ( GRZ0M , GRZ0H , GRZ0E ,
     I           USFC  )
      use zocdim
      use zocphy
      IMPLICIT NONE
*
*      roughness of sea surface (Miller et al., 1992)
*
*   [OUTPUT] 
      REAL*8     GRZ0M  ( NXYDIM )        !! surface roughness (V)
      REAL*8     GRZ0H  ( NXYDIM )        !! surface roughness (T)
      REAL*8     GRZ0E  ( NXYDIM )        !! surface roughness (q)
*
*   [INPUT] 
      REAL*8     USFC   ( NXYDIM )        !! sfc wind speed
*
*   [INTERNAL WORK] 
      INTEGER    IJ
      REAL*8     USTAR
      INTEGER    IFPAR, JFPAR
*
*   [INTERNAL PARM]     
      REAL*8     Z0M0                     !! base
      REAL*8     Z0MR                     !! rough factor
      REAL*8     Z0MS                     !! smooth factor
      REAL*8     Z0H0                     !! base
      REAL*8     Z0HR                     !! rough factor
      REAL*8     Z0HS                     !! smooth factor
      REAL*8     Z0E0                     !! base
      REAL*8     Z0ER                     !! rough factor
      REAL*8     Z0ES                     !! smooth factor
      REAL*8     VISAIR                   !! kinematic viscosity 
      REAL*8     CM0                      !! bulk coef for USTAR
      REAL*8     USTRMN                   !! min(u*)
      REAL*8     Z0MMIN                   !! minimum
      REAL*8     Z0HMIN                   !! minimum
      REAL*8     Z0EMIN                   !! minimum
      LOGICAL    OFIRST
      NAMELIST  /NMSEAZ/
     &           Z0M0, Z0MR, Z0MS,
     &           Z0H0, Z0HR, Z0HS,
     &           Z0E0, Z0ER, Z0ES,
     &           VISAIR, CM0, USTRMN, 
     &           Z0MMIN, Z0HMIN, Z0EMIN
      DATA       Z0M0, Z0MR, Z0MS / 0.,     0.018, 0.11  /
      DATA       Z0H0, Z0HR, Z0HS / 1.4E-5, 0.   , 0.4   /
      DATA       Z0E0, Z0ER, Z0ES / 1.3E-4, 0.   , 0.62  /
      DATA       VISAIR           / 1.5E-5  / 
      DATA       CM0              / 1.D-3   /
      DATA       USTRMN           / 1.D-3   /
      DATA       Z0MMIN, Z0HMIN, Z0EMIN  / 3*1.E-5  /        
      DATA       OFIRST           / .TRUE. /
*
      IF ( OFIRST ) THEN
         OFIRST = .FALSE.
         CALL REWNML(IFPAR, JFPAR)
         WRITE(JFPAR, *) ' SEA SURFACE Z0 (Miller et al.) 98/06/19'
         READ(IFPAR, NMSEAZ, END=190)
  190    WRITE(JFPAR, NMSEAZ)
      ENDIF
*
      DO 1100 IJ = IJSTR, IJEND
         USTAR       = CM0 * USFC(IJ) * USFC(IJ)
         USTAR       = MAX(SQRT(USTAR), USTRMN)
*
         GRZ0M( IJ ) = Z0M0 + Z0MR*USTAR**2/GRAV + Z0MS*VISAIR/USTAR
         GRZ0H( IJ ) = Z0H0 + Z0HR*USTAR**2/GRAV + Z0HS*VISAIR/USTAR
         GRZ0E( IJ ) = Z0E0 + Z0ER*USTAR**2/GRAV + Z0ES*VISAIR/USTAR
*
         GRZ0M( IJ ) = MAX( GRZ0M( IJ ), Z0MMIN )            
         GRZ0H( IJ ) = MAX( GRZ0H( IJ ), Z0HMIN )            
         GRZ0E( IJ ) = MAX( GRZ0E( IJ ), Z0EMIN )            
 1100 CONTINUE
*
      RETURN
      END

c *********************************************************************

      SUBROUTINE SFCFLX_OMIP
     O         ( TFLUXS, QFLUXS,
     O           DTFDT , DTFDS , DQFDS , 
     I           USFC  , GDTA  , GDQA  , GDPS  ,
     I           GDTS  , GRZ0  )
      use zocdim
      use zocphy
      IMPLICIT NONE

      INTEGER NTYZ0
      PARAMETER(NTYZ0 = 3)
*
*   [MODIFY] 
      REAL*8     TFLUXS( NXYDIM )              !! flux of T
      REAL*8     QFLUXS( NXYDIM )              !! flux of q
*
*   [OUTPUT]  
      REAL*8     DUFDU ( NXYDIM )              !! -d(tau)/du
      REAL*8     DTFDT ( NXYDIM )              !! -dH/dTa
      REAL*8     DQFDQ ( NXYDIM )              !! -dLE/dqa
      REAL*8     DTFDS ( NXYDIM )              !! dH/dTs
      REAL*8     DQFDS ( NXYDIM )              !! dLE/dTs
*
*   [INPUT] 
      REAL*8     USFC  ( NXYDIM )              !! westerly u
      REAL*8     GDTA  ( NXYDIM )              !! temperature T
      REAL*8     GDQA  ( NXYDIM )              !! humidity q
      REAL*8     GDPS  ( NXYDIM )
      REAL*8     GDTS  ( NXYDIM )              !! surface skin temp
      REAL*8     GRZ0  ( NXYDIM, NTYZ0 )       !! surface roughness
*
*   [INTERNAL WORK] 
      INTEGER    IJ
      REAL*8     AKAPPA, EXI, RHO, QSAT, DQSAT
      INTEGER    IFPAR, JFPAR
      SAVE AKAPPA
*
*   [INTERNAL PARM] 
      LOGICAL    OFIRST

      REAL*8 ZA
      SAVE ZA
      NAMELIST /NMSAIR/ ZA
*
*   [INTERNAL FUNC] 
      DATA       OFIRST / .TRUE. /
      DATA ZA / 2.D0 /
*
!===== Define statement function 
#include        "zqsat.F"                /* saturate q */
!===== 
*
      IF ( OFIRST ) THEN
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMSAIR, END=110)
 110     WRITE(JFPAR, NMSAIR)
         OFIRST = .FALSE.
         WRITE(JFPAR, *) ' @@@ PSFCM: SURFACE FLUX 98/06/19'
         AKAPPA = RAIR/CP
      ENDIF
*
*
      CALL BLKCOF_OMIP
     O         ( DUFDU , DTFDS , DQFDQ , !! CMV, CHV, CEV            
     I           USFC  , GDTA  , GDQA  , GDTS  , GDPS  , GRZ0  )
*
      DO 6100 IJ = IJSTR, IJEND
         RHO = GDPS(IJ) / (RAIR + (RVAP - RAIR) * GDQA(IJ)) / GDTA(IJ)
         EXI = (GDPS(IJ) / (GDPS(IJ) - RHO * GRAV * ZA))**AKAPPA

         QSAT         = FQSAT ( GDTS( IJ ), GDPS(IJ) )
         DQSAT        = FDQSAT( GDTS( IJ ), QSAT )
*     
         DUFDU ( IJ ) =    RHO*DUFDU( IJ )
         DTFDS ( IJ ) = CP*RHO*DTFDS( IJ )
         DQFDQ ( IJ ) =    RHO*DQFDQ( IJ )
*
         DTFDT ( IJ ) =   DTFDS( IJ )*EXI
         DQFDS ( IJ ) =   DQFDQ( IJ )*DQSAT
*
         TFLUXS( IJ ) =   DTFDS( IJ )*( GDTS( IJ )-GDTA( IJ )*EXI )
         QFLUXS( IJ ) =   DQFDQ( IJ )*( QSAT      -GDQA( IJ )     )
 6100 CONTINUE
*
      RETURN
      END

c *********************************************************************

      SUBROUTINE BLKCOF_OMIP
     O         ( CMV   , CHV   , CEV   ,
     I           USFC  , GDTA  , GDQA  , GDTS  , GDPS  , GRZ0  )
      use zocdim
      use zocphy
      IMPLICIT NONE

      INTEGER NTYZ0
      PARAMETER(NTYZ0 = 3)
*
*   [OUTPUT] 
      REAL*8     CMV   ( NXYDIM )         !! bulk transfer coeff.: u
      REAL*8     CHV   ( NXYDIM )         !! bulk transfer coeff.: T
      REAL*8     CEV   ( NXYDIM )         !! bulk transfer coeff.: q
*
*   [INPUT] 
      REAL*8     USFC  ( NXYDIM )              !! u wind speed
      REAL*8     GDTA  ( NXYDIM )              !! temperature T
      REAL*8     GDTS  ( NXYDIM )              !! surface skin temp
      REAL*8     GDQA  ( NXYDIM )
      REAL*8     GDPS  ( NXYDIM )
      REAL*8     GRZ0  ( NXYDIM, NTYZ0 )       !! surface roughness
*
*   [INTERNAL WORK] 
      REAL*8     AKAPPA, EXI, RHO, VABS
      REAL*8 CL0, CL1
      INTEGER    IJ
      INTEGER    IFPAR, JFPAR
      SAVE AKAPPA
*
*   [INTERNAL PARM] 
      REAL*8     CMMIN, CHMIN, CEMIN
      REAL*8     CMMAX, CHMAX, CEMAX
      REAL*8     USMINM, USMINH, USMINE
      REAL*8     USMAXM, USMAXH, USMAXE
      LOGICAL    OFIRST
      REAL*8  ZA
      SAVE ZA
      NAMELIST /NMSAIR/ ZA
      NAMELIST  /NMSFCL/ CMMIN , CHMIN , CEMIN ,
     &                   CMMAX , CHMAX , CEMAX ,
     &                   USMINM, USMAXM,
     &                   USMINH, USMAXH,
     &                   USMINE, USMAXE
*
      DATA ZA / 2.D0 /
      DATA       CMMIN  /  1.E-5 /           !! min. bulk coef. of u
      DATA       CHMIN  /  1.E-5 /           !! min. bulk coef. of T
      DATA       CEMIN  /  1.E-5 /           !! min. bulk coef. of q
      DATA       CMMAX  /   1.   /           !! max. bulk coef. of u
      DATA       CHMAX  /   1.   /           !! max. bulk coef. of T
      DATA       CEMAX  /   1.   /           !! max. bulk coef. of q
      DATA       USMINM / 2.0    /           !! min. wind vel. for V
      DATA       USMINH / 2.0    /           !! min. wind vel. for T
      DATA       USMINE / 2.0    /           !! min. wind vel. for q
      DATA       USMAXM / 1000.  /           !! max wind vel. for V
      DATA       USMAXH / 1000.  /           !! max. wind vel. for heat
      DATA       USMAXE / 1000.  /           !! max wind vel. for q
      DATA       OFIRST / .TRUE. /
      
*
*   [ONCE] 
      IF ( OFIRST ) THEN
         CALL REWNML(IFPAR, JFPAR)
         WRITE(JFPAR, *) ' @@@ PSFCL: SURFACE BULK COEF. Kara' 
         READ(IFPAR, NMSAIR, END=110)
 110     WRITE(JFPAR, NMSAIR)
         OFIRST = .FALSE.
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMSFCL, END=190 )
  190    WRITE(JFPAR, NMSFCL)
         AKAPPA = RAIR / CP
      ENDIF

      DO IJ = IJSTR, IJEND
         RHO = GDPS(IJ) / (RAIR + (RVAP - RAIR) * GDQA(IJ)) / GDTA(IJ)
         EXI = (GDPS(IJ) / (GDPS(IJ) - RHO * GRAV * ZA))**AKAPPA

         VABS      = USFC(IJ)
         VABS      = MIN( MAX( VABS, USMINE ), USMAXE )
         CL0 = 1.D-3 * (   0.994D0 
     &                   + 0.061D0 * VABS
     &                   - 0.001D0 * VABS * VABS)
         CL1 = 1.D-3 * ( - 0.020D0 
     &                   + 0.691D0 * 1.D0 / VABS
     &                   - 0.817D0 * 1.D0 / VABS / VABS)
         CEV(IJ) = CL0 + CL1 * ( GDTS( IJ )-GDTA( IJ )*EXI ) 
         CHV(IJ) = 0.96D0 * CEV(IJ)
      ENDDO

      DO 5100 IJ = IJSTR, IJEND
         CMV( IJ ) = MAX( MIN( CMV( IJ ), CMMAX ), CMMIN )
         CHV( IJ ) = MAX( MIN( CHV( IJ ), CHMAX ), CHMIN )
         CEV( IJ ) = MAX( MIN( CEV( IJ ), CEMAX ), CEMIN )
 5100 CONTINUE

      DO 6100 IJ = IJSTR, IJEND
         VABS      = USFC(IJ)
         CMV( IJ ) = CMV( IJ )
     &               * MIN( MAX( VABS, USMINM ), USMAXM )
         CHV( IJ ) = CHV( IJ )
     &               * MIN( MAX( VABS, USMINH ), USMAXH )
         CEV( IJ ) = CEV( IJ )
     &               * MIN( MAX( VABS, USMINE ), USMAXE )
 6100 CONTINUE 

      RETURN
      END
#ifdef OPT_BODY
c *********************************************************************

      SUBROUTINE BDYFLX(
     O                 TQ,
     I                  T)
      use zocdim
      use zocgrd
      use zocmsk
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Estimate the surface forcing fluxes from the data given in the
c boundary condition files
c
c  HISTORY
c     '99.08.25  H.Hasumi
c     '02.06.02  H.Hasumi: tracer dimension
c
c ---------------------------------------------------------------------

      REAL*8      TQ(NXYDIM, NZDIM, NTDIM)
      REAL*8       T(NXYDIM, NZDIM, NTDIM)

      INTEGER     IJ,      K,      L,      N

      REAL*8  TBDY(NXYDIM, NZDIM, NTDIM), TDMB(NXYDIM, NZDIM, NTDIM)
c      COMMON /WORK/ TBDY, TDMB

c      DO L = 1, NTDIM
      DO L = 1, 2
         N = L
         CALL TMINTB(
     O                 TBDY(1, 1, L),
     I                    N)
      END DO
c      DO L = 1, NTDIM
c         N = L + NTDIM
      DO L = 1, 2
         N = L + 2
         CALL TMINTB(
     O                 TDMB(1, 1, L),
     I                    N)
      END DO

c      DO 120 L = 1, NTDIM
      DO 120 L = 1, 2
         DO 110 K = KSTR, KEND
            DO 100 IJ = IJTSTR, IJTEND
               TQ(IJ, K, L) = TDMB(IJ, K, L) *
     &                        (TBDY(IJ, K, L) - T(IJ, K, L)) *
     &                        AMSKT(IJ, K)
 100        CONTINUE
 110     CONTINUE
 120  CONTINUE

      RETURN
      END
#endif
