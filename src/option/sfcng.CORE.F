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
c     '08.06.30  Y.Komuro: from COCO4.2
c     '08.08.06  Y.Komuro: for COCO4.4; initial/final processing added
c     '08.08.26  Y.Komuro: diagnosing evapolation
c     '08.09.02  Y.Komuro: diagnosing wind stress
c     '08.09.08  Y.Komuro: bug(s) fix in BLKCOF_CORE
c     '08.09.16  Y.Komuro: latitude-dependent albedo
c     '08.09.19  Y.Komuro: bug fix in BLKCOF_CORE
c     '09.05.25  Y.Komuro: CMIP5 output code included
c     '09.09.26  Y.Komuro: bug fix
c     '10.04.14  M.kurogi
c     '10.04.14  M.Kurogi: for tripolar grid
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

      REAL*8    TAUX(NXYDIM),   TAUY(NXYDIM)
      REAL*8    USFC(NXYDIM),   VSFC(NXYDIM)
      REAL*8    TSFC(NXYDIM),   QSFC(NXYDIM),   PPLR(NXYDIM)
      REAL*8    SFLX(NXYDIM)
      REAL*8    SWNT(NXYDIM),   DWLW(NXYDIM),   PSFC(NXYDIM)
      REAL*8    GRTS(NXYDIM),   GRTB(NXYDIM)
      REAL*8   GRICE(NXYDIM),  GRSNW(NXYDIM),  GRICR(NXYDIM)
      REAL*8    GRZ0(NXYDIM, NTYZ0)
      REAL*8  GFLUXS(NXYDIM), TFLUXS(NXYDIM), QFLUXS(NXYDIM)
      REAL*8  WFLUXS(NXYDIM, 2)
      REAL*8  RFLXLU(NXYDIM), SFLXBL(NXYDIM)
      REAL*8   DGFDS(NXYDIM),  DTFDT(NXYDIM),  DTFDS(NXYDIM)
      REAL*8   DQFDS(NXYDIM),   SWDN(NXYDIM)
      REAL*8      FM(NXYDIM)

      REAL*8     TMI
      SAVE TMI

      REAL*8     DUFDU ( NXYDIM )              !! CMV
c      REAL*8     DTFDS ( NXYDIM )              !! CHV
      REAL*8     DQFDQ ( NXYDIM )              !! CEV
c      REAL*8     CMV(NXYDIM, 0:NIC), CHV(NXYDIM, 0:NIC)
c      REAL*8     CEV(NXYDIM, 0:NIC)
      
c      REAL*8 RALBSW(NXYDIM), ALBSW(NXYDIM, 0:NIC)
c      REAL*8 RWSB(NXYDIM, 0:NIC)
c      REAL*8    RQAI(NXYDIM, 0:NIC),   RQIO(NXYDIM, 0:NIC)
c      REAL*8    RQII(NXYDIM, 0:NIC)
c      REAL*8  FTATM(NXYDIM), SWNTWA(NXYDIM)
      REAL*8 RALBSW(NXYDIM, 0:NIC)
      REAL*8   TISI(NXYDIM, 0:NIC)
      REAL*8   WSBG(NXYDIM)

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
      INTEGER      I,      J
      INTEGER  IFPAR,  JFPAR
      LOGICAL OFIRST
      SAVE OFIRST

      REAL*8     SI
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

      DO IJ = 1, NXYDIM
c         FTATM(IJ) = 0.0D0
c         SWNTWA(IJ) = 0.0D0
         WSBG(IJ) = 0.0D0
      END DO
      DO L = 0, NIC
         DO IJ = 1, NXYDIM
c            RQAI(IJ, L) = 0.0D0
c            RQIO(IJ, L) = 0.0D0
c            RQII(IJ, L) = 0.0D0
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
      CALL TMINTP(  USFC,      1)
      CALL TMINTP(  VSFC,      2)
      CALL TMINTP(  TSFC,      3)
      CALL TMINTP(  QSFC,      4)
      CALL TMINTP(  PPLR,      5)
      CALL TMINTP(  SFLX,      6)
      CALL TMINTP(  SWNT,      7)
      CALL TMINTP(  DWLW,      8)
      CALL TMINTP(  PSFC,      9)
#ifdef OPT_SRST
      CALL TMINTP(  SSFC,     10)
#endif

!!!!!! kurogi 2009.11.27
#ifdef OPT_PARALLEL
#ifdef OPT_TRIPOLE
      CALL SHIFT1(
     M             PSFC,
     D             NXDIM,  NYDIM,      1,
     I             1.D0,      0,      0 )
#else
      CALL SHIFT1(
     M             PSFC,
     D             NXDIM,  NYDIM,      1)
#endif
#else
      CALL STBCSH(
     M              PSFC)
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      DO IJ = IJSTR, IJEND
         TAUAIX(IJ) = 0.0D0
         TAUAIY(IJ) = 0.0D0
         PREC(IJ) = - PPLR(IJ) * FACTW
         SNOW(IJ) = - SFLX(IJ) * FACTW
         SOFF(IJ) = 0.D0
         ROFF(IJ) = 0.D0
      END DO

      DO L = 0, NIC
         DO IJ = 1, NXYDIM
            FM(IJ) = A(IJ, L)
         END DO

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

         CALL OCNBCS_CORE(
     O              GFLUXS,  DGFDS,
     I              GRTS,   GRTB,  GRICE,  GRSNW,  GRICR )
         CALL SFCFLX_CORE(
     O            TFLUXS, QFLUXS,   TAUX,   TAUY,
     O             DTFDT,  DTFDS,  DQFDS,
     O             DUFDU,  DQFDQ,
     I              USFC,   VSFC,   TSFC,   QSFC,   PSFC,
     I              GRTS,      L )

c         DO IJ=1, NXYDIM
c            CMV(IJ, L) = DUFDU(IJ)
c            CHV(IJ, L) = DTFDS(IJ)
c            CEV(IJ, L) = DQFDQ(IJ)
c         END DO

         IF (L .GT. 0) THEN
            DO IJ = IJSTR, IJEND
               GRICE(IJ) = HI(IJ, L) * 1.D-2
            END DO
         END IF
         CALL OCNSLV_CORE(
     M              GRTS, GFLUXS, TFLUXS, QFLUXS,
     O            WFLUXS, RFLXLU, SFLXBL,   SWDN, RALBSW,
     I             DTFDS,  DQFDS,  DGFDS,
     I              SWNT,   DWLW,
     I             GRICR,  GRSNW,    TMI)

#ifdef OPT_PARALLEL
#ifdef OPT_TRIPOLE
      CALL SHIFT2(
     M              TAUX,   TAUY,
     D             NXDIM,  NYDIM,      1,
     I             -1.D0,      0,      0 )
      CALL SHIFT1(
     M                FM,
     D             NXDIM,  NYDIM,      1,
     I              1.D0,      0,      0 )
#else
      CALL SHIFT3(
     M              TAUX,   TAUY,    FM,
     D             NXDIM,  NYDIM,     1)
#endif
#else
      CALL STBCTU(
     M              TAUX,   TAUY)
      CALL STBCSH(
     M                FM)
#endif

         IF (L .GT. 0) THEN
            DO IJ = IJSTR, IJEND
               TAUAIX(IJ) = TAUAIX(IJ) 
     &         + ( TAUX(IJ) + TAUX(IJ+LE)
     &           + TAUX(IJ+LN) + TAUX(IJ+LNE) )
     &                     * 0.25D0 * FACTM * FM(IJ)
               TAUAIY(IJ) = TAUAIY(IJ)
     &         + ( TAUY(IJ) + TAUY(IJ+LE)
     &           + TAUY(IJ+LN) + TAUY(IJ+LNE) )
     &                     * 0.25D0 * FACTM * FM(IJ)
               QAI(IJ, L) = (GFLUXS(IJ) - SFLXBL(IJ)) * FACTH
               QII(IJ, L) = GFLUXS(IJ) * FACTH
               QIO(IJ, L) = (T(IJ, KSTR, 1) - TI(IJ, L)) *
     &                      cdi / HI(IJ, L) * 2.D0
c               RQAI(IJ, L) = QAI(IJ, L)
c               RQII(IJ, L) = QII(IJ, L)
c               RQIO(IJ, L) = QIO(IJ, L)
               SOFF(IJ) = 0.D0
               TSI(IJ, L) = MIN(GRTS(IJ), TMELT+TMI) - kelvin
               TISI(IJ, L) = TSI(IJ, L) +
     &                       (TI(IJ, L) - TSI(IJ, L)) *
     &                       DIRDSN * HSN(IJ, L) /
     &                       (0.5D0 * HI(IJ, L) + DIRDSN * HSN(IJ, L))
c               WSB(IJ, L) = 0.D0 ! (sublimation speed) x (fraction)
               WSB(IJ, L) = QFLUXS(IJ) * A(IJ, L) / DWATR * FACTW
               WSBG(IJ) = WSBG(IJ) +
     &                      QFLUXS(IJ) * A(IJ, L) / DWATR * FACTW
c               RWSB(IJ, L) = WSB(IJ, L)
c               ALBSW(IJ, L) = RALBSW(IJ)
c               FTATM(IJ) = FTATM(IJ) + QAI(IJ, L) * FM(IJ)
c               SWNTWA(IJ) = SWNTWA(IJ) + SWDN(IJ) * FM(IJ)
            END DO
         ELSE
            DO IJ = IJSTR, IJEND
               TAUAOX(IJ) = (  TAUX(IJ) + TAUX(IJ+LE)
     &                       + TAUX(IJ+LN) + TAUX(IJ+LNE)) *
     &                      0.25D0 * FACTM
               TAUAOY(IJ) = (  TAUY(IJ) + TAUY(IJ+LE)
     &                       + TAUY(IJ+LN) + TAUY(IJ+LNE)) *
     &                      0.25D0 * FACTM
               QAO(IJ) = (GFLUXS(IJ) + SWDN(IJ)) * FACTH
               SWABS(IJ) = SWDN(IJ) * FACTH * A(IJ, 0)
c               WEV(IJ) = 0.D0 ! (evaporation speed) x (fraction)
               WEV(IJ) = QFLUXS(IJ) * A(IJ, 0) / DWATR * FACTW
c               ALBSW(IJ, L) = RALBSW(IJ)
c               FTATM(IJ) = FTATM(IJ) + GFLUXS(IJ) * FACTH * A(IJ, 0)
c               SWNTWA(IJ) = SWNTWA(IJ) + SWDN(IJ) * A(IJ, 0)
            END DO
         END IF
      END DO

      DO IJ = 1, NXYDIM
         IF( A(IJ,0) .ne. 1.D0 ) THEN
            TAUAIX(IJ) = TAUAIX(IJ) / (1.0D0 - A(IJ, 0) )
            TAUAIY(IJ) = TAUAIY(IJ) / (1.0D0 - A(IJ, 0) )
         ENDIF
      ENDDO

#ifdef OPT_PARALLEL
#ifdef OPT_TRIPOLE
      CALL SHIFT2(TAUAOX,   TAUAOY,
     D             NXDIM,    NYDIM,    1,
     I             -1.D0,     -1,     -1 )

      CALL SHIFT2(TAUAIX,   TAUAIY,
     D             NXDIM,    NYDIM,    1,
     I             -1.D0,     -1,     -1 )
#else
      CALL SHIFT2(
     M              TAUAOX,   TAUAOY,
     D               NXDIM,    NYDIM,     1)
      CALL SHIFT2(
     M              TAUAIX,   TAUAIY,
     D               NXDIM,    NYDIM,     1)
#endif
#else
      CALL STBCTU(
     M              TAUAOX,   TAUAOY)
      CALL STBCTU(
     M              TAUAIX,   TAUAIY)
#endif

      CALL CHEKIN(WSBG,  'WSBG',
     &            NX,     NY,      1, NXYDIM, 'SFC')
      CALL CHEKIN(WEV,   'WEV',
     &            NX,     NY,      1, NXYDIM, 'SFC')
c      CALL CHEKIN(TAUAOX,'TAUAOX',
c     &            NX,     NY,      1, NXYDIM, 'SFC')
c      CALL CHEKIN(TAUAOY,'TAUAOY',
c     &            NX,     NY,      1, NXYDIM, 'SFC')
c      CALL CHEKIN(TAUAIX,'TAUAIX',
c     &            NX,     NY,      1, NXYDIM, 'SFC')
c      CALL CHEKIN(TAUAIY,'TAUAIY',
c     &            NX,     NY,      1, NXYDIM, 'SFC')
c      CALL CHEKIN(CMV,     'CMVO',
c     &            NX,     NY,      1, NXYDIM, 'SFC')
c      CALL CHEKIN(CHV,     'CHVO',
c     &            NX,     NY,      1, NXYDIM, 'SFC')
c      CALL CHEKIN(CEV,     'CEVO',
c     &            NX,     NY,      1, NXYDIM, 'SFC')
c      CALL CHEKIN(CMV,     'CMVI',
c     &            NX,     NY,    NIC, NXYIDM, 'ICE')
c      CALL CHEKIN(CHV,     'CHVI',
c     &            NX,     NY,    NIC, NXYIDM, 'ICE')
c      CALL CHEKIN(CEV,     'CEVI',
c     &            NX,     NY,    NIC, NXYIDM, 'ICE')
c      CALL CHEKIN(ALBSW, 'ALBSWO',
c     &            NX,     NY,      1, NXYDIM, 'SFC')
c      CALL CHEKIN(ALBSW, 'ALBSWI',
c     &            NX,     NY,    NIC, NXYIDM, 'ICE')
c      CALL CHEKIN(FTATM, 'FTATM',
c     &            NX,     NY,      1, NXYDIM, 'SFC')
c      CALL CHEKIN(SWNTWA, 'SWNT',
c     &            NX,     NY,      1, NXYDIM, 'SFC')
c      CALL CHEKIN(RQAI,    'QAI',
c     &            NX,     NY,    NIC, NXYIDM, 'ICE')
c      CALL CHEKIN(RQIO,    'QIO',
c     &            NX,     NY,    NIC, NXYIDM, 'ICE')
c      CALL CHEKIN(RQII,    'QII',
c     &            NX,     NY,    NIC, NXYIDM, 'ICE')

c     output section for CMIP5
c     TISI: temperature at interface between sea ice and snow
c           unit [C]
      CALL CHEKIN(TISI,   'TISI',
     &            NX,     NY,    NIC, NXYIDM, 'ICE')

      DO IJ = 1, NXYDIM
c         PTOP(IJ) = 0.D0
         PTOP(IJ) = PSFC(IJ) * FACTM
      END DO

      RETURN
      END

c *********************************************************************

      SUBROUTINE OCNSLV_CORE
     M         ( GDTS  , GFLUXS, TFLUXS, QFLUXS,
     O           WFLUXS, RFLXLU, SFLXBL, RFLXSD, RALBSW,
     I           DTFDS , DQFDS , DGFDS ,
     I           RFLXS , RFLXLD,
     I           GRICR , GRSNW ,    TMI) 
      use zocdim
      use zocgrd
      use zocfil
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
      REAL*8     RFLXSD( NXYDIM )                !! downward SW
      REAL*8     RALBSW( NXYDIM )                !! SW albedo
*
*   [INPUT] 
      REAL*8     DTFDS ( NXYDIM )                !! dH/dTg
      REAL*8     DQFDS ( NXYDIM )                !! dE/dTg
      REAL*8     DGFDS ( NXYDIM )                !! dG/dTg
      REAL*8     RFLXS ( NXYDIM )                !! SW rad. (net/down.)
      REAL*8     RFLXLD( NXYDIM )                !! down. LW rad.
      REAL*8     GRICR ( NXYDIM )                !! snow/ice ratio
      REAL*8     GRSNW ( NXYDIM )                !! snow thickness
      REAL*8     TMI ! sea ice (and snow) melting temperature in C
*
*   [INTERNAL WORK] 
      REAL*8     ESUB, STG, DRFDS
      REAL*8     SFLUX, GSFLUX, DGSFDS
      REAL*8     GFLUXF
      REAL*8     SFLXBI, DSBDSI, DTI, EVAPI, GFLUXI
      REAL*8     FF, FI, DTX
      REAL*8     EMIS, ALBSW, SNWALB, ICEALB, FBAREI, X
      REAL*8     OMEGA, SINIJ, CORT
      INTEGER    IJ
      INTEGER  IFPAR,  JFPAR
*
*   [INTERNAL SAVE] 
      LOGICAL OFIRST
      LOGICAL OLWNET   !! TRUE if longwave rad. is net flux
      LOGICAL OSWNET   !! TRUE if shortwave rad. is net flux
      LOGICAL OASOLD   !! TRUE if ocn. shortwave albedo depends on lat.
      REAL*8    EMISLO !! ocn. longwave emissivity (longwave co-albedo)
      REAL*8    EMISLI !! ice longwave emissivity (longwave co-albedo)
      REAL*8    ALBSWO !! ocn. shortwave albedo (if constant)
      REAL*8    ALBSWI !! bare ice shortwave albedo
      REAL*8    ALBSWS(2) !! snow-covered ice shortwave albedo
      REAL*8    ALSWST(2) !! snow-covered ice SW albedo, T dependency
      REAL*8    ALSDPT !! threshold thickness of snow-covered ice
      REAL*8    ASWO2D(NXYDIM) !! ocn. shortwave albedo distribution
      REAL*8 FLWNET, FSWALB, AICET1, DAICET, TSDPT
      NAMELIST /NMLWEM/ EMISLO, EMISLI, OLWNET
      NAMELIST /NMSWEM/ ALBSWO, ALBSWI, ALBSWS, ALSWST, ALSDPT,
     &                  OSWNET, OASOLD
      SAVE OFIRST, OLWNET, OSWNET, OASOLD, EMISLO, EMISLI, ALBSWO,ALBSWI
      SAVE ALBSWS, ALSWST, ALSDPT, FLWNET, FSWALB, AICET1, DAICET, TSDPT
      SAVE ASWO2D
      DATA EMISLO, EMISLI / 1.D0, 0.95D0 /
      DATA ALBSWO, ALBSWI / 0.066D0, 0.5D0 /
      DATA ALBSWS / 0.85D0, 0.65D0 /
      DATA ALSWST / -15.0D0, 0.0D0 /
      DATA ALSDPT / 1.0D0 /
      DATA OLWNET / .FALSE. /
      DATA OSWNET / .FALSE. /
      DATA OASOLD / .FALSE. /
      DATA       OFIRST / .TRUE.  /
*
      IF ( OFIRST ) THEN
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMLWEM, END=110)
 110     WRITE(JFPAR, NMLWEM)
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMSWEM, END=109)
 109     WRITE(JFPAR, NMSWEM)
         OFIRST = .FALSE.
         WRITE (JFPAR, *) ' @@@ OCNSLV: OCEAN HEAT BALANCE 98/06/30'
         IF (OLWNET) THEN
            FLWNET = 0.D0
         ELSE
            FLWNET = 1.D0
         END IF
         IF (OSWNET) THEN
            FSWALB = 0.D0
         ELSE
            FSWALB = 1.D0
         END IF
         AICET1 = ALSWST(1) + kelvin
         DAICET = ALSWST(2) - ALSWST(1)
         TSDPT = ALSDPT * 1.0D-2

c        '08.09.16: Large and Yeager (2008) latitude-dependant albedo
         IF (OASOLD) THEN
            WRITE (JFPAR, *)
     &       '== L-Y2008 latitude-dependant ocn. shortwave albedo =='
            OMEGA = 2.D0 * ( ATAN( 1.D0 )*4.D0 ) / 86400.D0
            DO IJ = IJSTR, IJEND
               CORT=(COR(IJ)+COR(IJ+LW)+COR(IJ+LSW)+COR(IJ+LS))*0.25D0
               SINIJ = MIN(MAX(CORT/(2.0D0 * OMEGA), -1.0D0), 1.0D0)
               ASWO2D(IJ) = 0.069D0 - 0.011 *
     &                       (1.0D0 - 2.0D0 * SINIJ * SINIJ)
c               WRITE(JFPAR, *) IJ, CORT, SINIJ, ASWO2D(IJ)
            END DO
         ELSE
            WRITE (JFPAR, *)
     &       '== Constant ocn. shortwave albedo,', ALBSWO, ' =='
            DO IJ = IJSTR, IJEND
               ASWO2D(IJ) = ALBSWO
            END DO
         END IF
      ENDIF
*
      ESUB = EL + EMELT
      DO 1000 IJ = IJSTR, IJEND
         EMIS = EMISLO * (1.0D0 - GRICR(IJ)) + EMISLI * GRICR(IJ)
         X = MIN( MAX( (GDTS(IJ) - AICET1) / DAICET, 0.0D0), 1.0D0)
         SNWALB = ALBSWS(1) * (1.0D0 - X) + ALBSWS(2) * X
         FBAREI = 0.5D0 + SIGN(0.5D0, TSDPT - GRSNW(IJ))
         ICEALB = SNWALB * (1.0D0 - FBAREI) + ALBSWI * FBAREI
         ALBSW = (ASWO2D(IJ) * (1.0D0 - GRICR(IJ)) + ICEALB * GRICR(IJ))
     &            * FSWALB
         RALBSW(IJ) = ALBSW
         STG            = EMIS*STB*GDTS( IJ )**4
         RFLXLU( IJ )   = (STG + ( 1.D0-EMIS )*RFLXLD( IJ )) * FLWNET
         DRFDS          = 4.D0*STG/GDTS( IJ ) * FLWNET
*
         RFLXSD( IJ )   = RFLXS( IJ ) * (1.0D0 - ALBSW)
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
      use zocfil
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
*   [INTERNAL WORK] 
      REAL*8     GRSNR ( NXYDIM )               !! snow fraction
*
*   [INTERNAL WORK] 
      INTEGER    IJ, L, M
      INTEGER  IFPAR,  JFPAR
      INTEGER    IFG
      REAL*8     TALSNX, ALBSNX, ALB0, DALB, TFACT, ALBX
      REAL*8     Z00, DZ0
      REAL*8     DFGT, DFGX
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
*
      DATA       Z0FCT          / 0.1D0          /
      DATA       Z0MIN          / 1.D-6          /
*
      DATA       DZOCN          / 50.D0          /
      DATA       DFOCN          / 1.D10          /
      DATA       ALBLO          / 0.05D0           /
*
      DATA       OFIRST         / .TRUE.         /

c =====================================================================

      ENTRY      OCNBCS_CORE
     O         ( FOGFLX, DGFDS ,
     I           GRTS  , GRTB  , GRICE , GRSNW , GRICR )

      IF ( OFIRST ) THEN
         CALL REWNML(IFPAR, JFPAR)
         WRITE (JFPAR, *) ' @@@ OCNBCS: OCEAN SURFACE BC 98/07/29'
         OFIRST = .FALSE.
         CALL REWNML(IFPAR, JFPAR)
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

      SUBROUTINE SFCFLX_CORE
     O         ( TFLUXS, QFLUXS, TAUX  ,  TAUY  ,
     O           DTFDT , DTFDS , DQFDS ,
     O           DUFDU,  DQFDQ,
     I           USFC  , VSFC  , GDTA  , GDQA  , GDPS  ,
     I           GDTS  , L )
      use zocdim
      use zocfil
      use zocphy
c      use zqsat
      IMPLICIT NONE

      INTEGER NTYZ0
      PARAMETER(NTYZ0 = 3)
*
*   [MODIFY] 
      REAL*8     TFLUXS( NXYDIM )              !! flux of T
      REAL*8     QFLUXS( NXYDIM )              !! flux of q
      REAL*8       TAUX( NXYDIM )
      REAL*8       TAUY( NXYDIM )
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
      REAL*8     VSFC  ( NXYDIM )              !! southerly v
      REAL*8     GDTA  ( NXYDIM )              !! temperature T
      REAL*8     GDQA  ( NXYDIM )              !! humidity q
      REAL*8     GDPS  ( NXYDIM )
      REAL*8     GDTS  ( NXYDIM )              !! surface skin temp
      INTEGER       L                          !! 1 = ice, 2 = ocean
*
*   [INTERNAL WORK] 
      INTEGER    IJ
      INTEGER  IFPAR,  JFPAR
      REAL*8     AKAPPA, EXI, RHO, QSAT, DQSAT, CQSAT, UABS
      SAVE AKAPPA
*
*   [INTERNAL PARM] 
      LOGICAL    OFIRST
      SAVE OFIRST

      REAL*8 ZA
      REAL*8     CMMIN, CHMIN, CEMIN
      REAL*8     CMMAX, CHMAX, CEMAX
      REAL*8     USMINM, USMINH, USMINE
      REAL*8     USMAXM, USMAXH, USMAXE
      SAVE ZA, CMMIN, CHMIN, CEMIN,
     &     CMMAX, CHMAX, CEMAX,
     &     USMINM, USMINH, USMINE,
     &     USMAXM, USMAXH, USMAXE
      NAMELIST /NMSAIR/ ZA
      NAMELIST  /NMSFCL/ CMMIN , CHMIN , CEMIN ,
     &                   CMMAX , CHMAX , CEMAX ,
     &                   USMINM, USMAXM,
     &                   USMINH, USMAXH,
     &                   USMINE, USMAXE
*
*   [INTERNAL FUNC] 
#include        "zqsat.F"                /* saturate q */
      DATA       OFIRST / .TRUE. /
      DATA       ZA     / 10.D0 /
      DATA       CMMIN  / 1.D-5 /           !! min. bulk coef. of u
      DATA       CHMIN  / 1.D-5 /           !! min. bulk coef. of T
      DATA       CEMIN  / 1.D-5 /           !! min. bulk coef. of q
      DATA       CMMAX  / 1.D0  /           !! max. bulk coef. of u
      DATA       CHMAX  / 1.D0  /           !! max. bulk coef. of T
      DATA       CEMAX  / 1.D0  /           !! max. bulk coef. of q
      DATA       USMINM / 0.5D0 /           !! min. wind vel. for V
      DATA       USMINH / 0.5D0 /           !! min. wind vel. for T
      DATA       USMINE / 0.5D0 /           !! min. wind vel. for q
      DATA       USMAXM / 1.0D3 /           !! max wind vel. for V
      DATA       USMAXH / 1.0D3 /           !! max. wind vel. for heat
      DATA       USMAXE / 1.0D3 /           !! max wind vel. for q
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
C     '08.06.05
C     Bulk coeff. are now diagnosed without distinction of ocean or ice 
      CALL BLKCOF_CORE
     O      ( DUFDU , DTFDS , DQFDQ , !! CMV, CHV, CEV            
     I        USFC  , VSFC  , GDTA  , GDQA  , GDTS  , GDPS ,    L )
      IF (L .EQ. 0) THEN   !! oceanic surface
c         CALL BLKCOF_CORE
c     O         ( DUFDU , DTFDS , DQFDQ , !! CMV, CHV, CEV            
c     I           USFC  , GDTA  , GDQA  , GDTS  , GDPS )
         CQSAT = 0.98D0
      ELSE   !! sea-ice surface
c         DO 6000 IJ = IJSTR, IJEND
c            UABS      = USFC(IJ)
c            DUFDU(IJ) = 1.63D-3
c     &               * MIN( MAX( UABS, USMINM ), USMAXM )
c            DTFDS(IJ) = 1.63D-3
c     &               * MIN( MAX( UABS, USMINH ), USMAXH )
c            DQFDQ(IJ) = 1.63D-3
c     &               * MIN( MAX( UABS, USMINE ), USMAXE )
c 6000    CONTINUE
         CQSAT = 1.0D0
      ENDIF
*
      RHO = 1.22D0
      DO 6100 IJ = IJSTR, IJEND
c         RHO = GDPS(IJ) / (RAIR + (RVAP - RAIR) * GDQA(IJ)) / GDTA(IJ)
         EXI = (GDPS(IJ) / (GDPS(IJ) - RHO * GRAV * ZA))**AKAPPA

         QSAT         = FQSAT ( GDTS( IJ ), GDPS(IJ) )
         DQSAT        = FDQSAT( GDTS( IJ ), QSAT )
         QSAT = CQSAT * QSAT
*     
         DUFDU ( IJ ) =    RHO*DUFDU( IJ )
         DTFDS ( IJ ) = CP*RHO*DTFDS( IJ )
         DQFDQ ( IJ ) =    RHO*DQFDQ( IJ )
*
         DTFDT ( IJ ) =   DTFDS( IJ )*EXI
         DQFDS ( IJ ) =   DQFDQ( IJ )*DQSAT
*
         TFLUXS( IJ ) =   DTFDS( IJ )*( GDTS( IJ )-GDTA( IJ )*EXI )
C     '08.08.29  Y.Komuro: QFLUXS must be >= 0
         QFLUXS( IJ ) =   DQFDQ( IJ )*
     &                    MAX( ( QSAT - GDQA( IJ ) ), 0.0D0 )
c         QFLUXS( IJ ) =   DQFDQ( IJ )*( QSAT      -GDQA( IJ )     )

         TAUX(IJ) = DUFDU(IJ) * USFC(IJ)
         TAUY(IJ) = DUFDU(IJ) * VSFC(IJ)
 6100 CONTINUE
*
      RETURN
      END

c *********************************************************************

      SUBROUTINE BLKCOF_CORE
     O         ( CMV   , CHV   , CEV   ,
     I           USFC  , VSFC  , GDTA  , GDQA  , GDTS  , GDPS,    L  )
      use zocdim
      use zocfil
      use zocphy
c      use zqsat
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
      REAL*8     USFC  ( NXYDIM )              !! 10 m westerly u
      REAL*8     VSFC  ( NXYDIM )              !! 10 m southerly v
      REAL*8     GDTA  ( NXYDIM )              !! 10 m air temperature
      REAL*8     GDTS  ( NXYDIM )              !! surface skin temp
      REAL*8     GDQA  ( NXYDIM )              !! 10 m specific hum.
      REAL*8     GDPS  ( NXYDIM )              !! sea level pressure
      INTEGER    L                        !! 1 for ice, 2 for ocean
*
*   [INTERNAL WORK] 
      REAL*8       PI,   EXI, AKAPPA
      REAL*8    THETA,THETAV,   RHO,  UABS,  QSAT, CMVSQR
      REAL*8    USTAR, TSTAR, QSTAR,  ZETA
      REAL*8        X, PSI_M, PSI_H,  UN10, FNCMVR,    XX,    X2
      REAL*8    ASTBL
c     bug fix 080422
      REAL*8    PHRCDK
c     bug fix 080908
      REAL*8    CDN10, CDN10R, CEN10, CHN10
      INTEGER      IJ,    N
      INTEGER  IFPAR,  JFPAR
      INTEGER  NITER
      SAVE  PI, EXI, AKAPPA
*
*   [INTERNAL PARM] 
      REAL*8     CMMIN, CHMIN, CEMIN
      REAL*8     CMMAX, CHMAX, CEMAX
      REAL*8     USMINM, USMINH, USMINE
      REAL*8     USMAXM, USMAXH, USMAXE
      LOGICAL    OFIRST
      REAL*8  ZA
      INTEGER NITERO, NITERI
      SAVE CMMIN, CHMIN, CEMIN, CMMAX, CHMAX, CEMAX
      SAVE USMINM, USMINH, USMINE, USMAXM, USMAXH, USMAXE
      SAVE OFIRST, ZA, NITERO, NITERI
      NAMELIST /NMSAIR/ ZA
      NAMELIST  /NMSFCL/ CMMIN , CHMIN , CEMIN ,
     &                   CMMAX , CHMAX , CEMAX ,
     &                   USMINM, USMAXM,
     &                   USMINH, USMAXH,
     &                   USMINE, USMAXE
      NAMELIST /NMSNIT/ NITERO, NITERI
*
      DATA       ZA     / 10.D0 /
      DATA       CMMIN  / 1.D-5 /           !! min. bulk coef. of u
      DATA       CHMIN  / 1.D-5 /           !! min. bulk coef. of T
      DATA       CEMIN  / 1.D-5 /           !! min. bulk coef. of q
      DATA       CMMAX  / 1.D0  /           !! max. bulk coef. of u
      DATA       CHMAX  / 1.D0  /           !! max. bulk coef. of T
      DATA       CEMAX  / 1.D0  /           !! max. bulk coef. of q
      DATA       USMINM / 0.5D0 /           !! min. wind vel. for V
      DATA       USMINH / 0.5D0 /           !! min. wind vel. for T
      DATA       USMINE / 0.5D0 /           !! min. wind vel. for q
      DATA       USMAXM / 1.0D3 /           !! max wind vel. for V
      DATA       USMAXH / 1.0D3 /           !! max. wind vel. for heat
      DATA       USMAXE / 1.0D3 /           !! max wind vel. for q
      DATA       NITERO / 2 /               !! N. of iteration for ocean
      DATA       NITERI / 5 /               !! N. of iteration for ice
      DATA       OFIRST / .TRUE. /
*
*   [INTERNAL FUNC] 
#include        "zqsat.F"                /* saturate q */
*
*   [ONCE] 
      IF ( OFIRST ) THEN
         CALL REWNML(IFPAR, JFPAR)
         WRITE(JFPAR, *)
     &        ' @@@ PSFCL: SURFACE BULK COEF. Large and Yeager' 
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMSAIR, END=110)
 110     WRITE(JFPAR, NMSAIR)
         OFIRST = .FALSE.
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMSFCL, END=190 )
  190    WRITE(JFPAR, NMSFCL)
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMSNIT, END=189 )
 189     WRITE(JFPAR, NMSNIT)
         PI = 4.0D0 * ATAN(1.0D0)
         AKAPPA = RAIR / CP
      ENDIF

      RHO = 1.22D0
      DO IJ = IJSTR, IJEND
c         RHO = GDPS(IJ) / (RAIR + (RVAP - RAIR) * GDQA(IJ)) / GDTA(IJ)
c         RHO = GDPS(IJ) / ( RAIR * GDTA(IJ)
c     &            * ( 1.0D0 + 0.608D0 * GDQA(IJ) ) )
         EXI = (GDPS(IJ) / (GDPS(IJ) - RHO * GRAV * ZA))**AKAPPA
         THETA = GDTA(IJ) * EXI
         THETAV = THETA * ( 1.0D0 + 0.608D0 * GDQA(IJ) )

         UABS = SQRT(USFC(IJ) * USFC(IJ) + VSFC(IJ) * VSFC(IJ))
         UABS = MIN( MAX( UABS, USMINE ), USMAXE )
c         QSAT = 0.98D0 * 640380.0D0 * exp( -5107.4D0 / GDTS(IJ)) / RHO
c        '08.09.08 bug fix:
c          Now the humidity factor (= 0.98) will be mulitiplied later,
c           since the value depend on surface condition (i.e., 0.98 or 1).
c         QSAT = 0.98D0 * FQSAT ( GDTS( IJ ), GDPS(IJ) )
         QSAT = FQSAT ( GDTS( IJ ), GDPS(IJ) )

c        First guess for bulk coefficients
C        '08.06.05
C        Bulk coeff. are now diagnosed without distinction of ocean or ice 
         IF (L .EQ. 0) THEN     !! oceanic surface
            QSAT = 0.98D0 * QSAT
            CDN10 =
     &         ( 2.70D0 / UABS + 0.142D0 + UABS / 13.09D0 ) * 1.0D-3
            CDN10R = SQRT( CDN10 )
            CEN10 = 34.6D-3 * CDN10R         
            ASTBL = 0.5D0 + SIGN( 0.5D0, THETA - GDTS(IJ) )
            CHN10 = ( 18.0D-3 * ASTBL + 32.7D-3 * (1.0D0 - ASTBL) )
     &                * CDN10R
            NITER = NITERO
         ELSE   !! sea-ice surface
            CDN10 = 1.63D-3
            CDN10R = SQRT( CDN10 )
            CEN10 = 1.63D-3
            CHN10 = 1.63D-3
            NITER = NITERI
         END IF
         CMV(IJ) = CDN10
         CEV(IJ) = CEN10
         CHV(IJ) = CHN10

c        Iteration for determing the coefficients
c        We assume the same height, 10 m above the sea surface,
c         for the atmospheric u, theta, and q.
c        In that case...
c         - zeta for u, theta, and q have the same value
c         - psi_h for u, theta, and q have the same value
c         - no need for shifting u, theta and q to 10 m height

         DO N=1, NITER
c           '08.09.19 bug fix: SQRT
            CMVSQR = SQRT( CMV(IJ) )
            USTAR = CMVSQR * UABS
            TSTAR = CHV(IJ) / CMVSQR * (THETA - GDTS(IJ))
            QSTAR = CEV(IJ) / CMVSQR * (GDQA(IJ) - QSAT)

            ZETA = CKARM * GRAV * ZA / ( USTAR * USTAR )
     &             * ( TSTAR / THETAV
     &               + QSTAR / ( GDQA(IJ) + 1.0D0/0.608D0 ) )

C --- from Dr. T.Suzuki's modification dated 080425 ---
C ---  originally from NCAR's code ? ---
            ZETA = SIGN( MIN( ABS(ZETA), 10.D0 ), ZETA)

            ASTBL = 0.5D0 + SIGN(0.5D0, ZETA)
            X2 = SQRT( 1.0D0 - 16.0D0 * MIN(ZETA, 0.0D0) )
            X = SQRT( X2 )
            PSI_M = ( -5.0D0 * ZETA ) * ASTBL 
     &              + ( 2.0D0 * LOG( ( 1.0D0 + X ) * 0.5D0 )
     &                  + LOG( (1.0D0 + X2 ) * 0.5D0 )
     &                  - 2.0D0 * ATAN(X) + PI * 0.5D0 )
     &                    * ( 1.0D0 - ASTBL )
            PSI_H = ( -5.0D0 * ZETA ) * ASTBL 
     &              + ( 2.0D0 * LOG( ( 1.0D0 + X2 ) * 0.5D0 ))
     &                  * ( 1.0D0 - ASTBL )

c           '08.09.08 bug fix:
c            The procedure for updating the transfar coefficients
c              was incorrect in the old code.
c            The coefficients must be diagnosed based on
c              neutral 10 m coefficients (i.e., CDN10 etc.).
c            --- the old code from here ---
c            XX = 1.0D0 / ( 1.0D0 - CMVSQR * PSI_M / CKARM )
c            UN10 = UABS * XX
c            CMV(IJ) = CMV(IJ) * XX * XX
cc           bug fix 080422
cc            XX = XX / ( 1.0D0 - CMVSQR * PSI_H / CKARM )
cc            CHV(IJ) = CHV(IJ) * XX
cc            CEV(IJ) = CEV(IJ) * XX
c            PHRCDK = PSI_H / (CMVSQR * CKARM)
c            CHV(IJ) = CHV(IJ) * XX / ( 1.0D0 - CHV(IJ) * PHRCDK )
c            CEV(IJ) = CEV(IJ) * XX / ( 1.0D0 - CEV(IJ) * PHRCDK )
c            CMVSQR = SQRT( CMV(IJ) )
c            --- the end of the old code ---
            IF (L .EQ. 0) THEN   !! re-calculate if oceanic surface
               UN10 = UABS / ( 1.0D0 - CDN10R * PSI_M / CKARM )
               CDN10 =
     &            ( 2.70D0 / UN10 + 0.142D0 + UN10 / 13.09D0 ) * 1.0D-3
               CDN10R = SQRT( CDN10 )
               CEN10 = 34.6D-3 * CDN10R         
               ASTBL = 0.5D0 + SIGN(0.5D0, ZETA)
               CHN10 = ( 18.0D-3 * ASTBL + 32.7D-3 * (1.0D0 - ASTBL) )
     &                   * CDN10R
            ENDIF
            
            XX = 1.0D0 / ( 1.0D0 - CDN10R * PSI_M / CKARM )
            CMV(IJ) = CDN10 * XX * XX
            PHRCDK = PSI_H / (CDN10R * CKARM)
            CHV(IJ) = CHN10 * XX / ( 1.0D0 - CHN10 * PHRCDK )
            CEV(IJ) = CEN10 * XX / ( 1.0D0 - CEN10 * PHRCDK )
         ENDDO
      ENDDO

      DO 5100 IJ = IJSTR, IJEND
         CMV( IJ ) = MAX( MIN( CMV( IJ ), CMMAX ), CMMIN )
         CHV( IJ ) = MAX( MIN( CHV( IJ ), CHMAX ), CHMIN )
         CEV( IJ ) = MAX( MIN( CEV( IJ ), CEMAX ), CEMIN )
 5100 CONTINUE

      DO 6100 IJ = IJSTR, IJEND
         UABS      = SQRT(USFC(IJ) * USFC(IJ) + VSFC(IJ) * VSFC(IJ))
         CMV( IJ ) = CMV( IJ )
     &               * MIN( MAX( UABS, USMINM ), USMAXM )
         CHV( IJ ) = CHV( IJ )
     &               * MIN( MAX( UABS, USMINH ), USMAXH )
         CEV( IJ ) = CEV( IJ )
     &               * MIN( MAX( UABS, USMINE ), USMAXE )
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
