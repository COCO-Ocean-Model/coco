module sfcng

! --- information -----------------------------------------------------
!
!  Estimate the surface forcing fluxes from the data given in the
! boundary condition files
!
!  HISTORY
!     '02.10.10  H.Hasumi: from COCO3.4 and MIROC3.1
!     '07.04.23  H.Hasumi
!     '07.09.25  H.Hasumi: multi-category sea ice
!     '07.09.26  H.Hasumi: sea ice (and snow) melting temperature
!                          is defined independently from TMELT
!     '07.09.28  H.Hasumi: 1-layer sea ice thermodynamics
!     '09.05.25  Y.Komuro: CMIP5 output code included
!     '09.09.26  Y.Komuro: bug fix
!     '10.04.14  M.Kurogi: (COCO4.4 tripolar code by Dr. Suzuki)
!     '12.09.04  H.Tatebe: Rewrite in F95 format
!     '21.05.26  Y.Komuro: snow aging & meltpond parametrization
!                (almost dummy, since albedo is implicit in this routine)
!
! ---------------------------------------------------------------------

  implicit none
#include "coco.h"
  private

  real(8),    parameter  ::  factm = 1.0d+1
  real(8),    parameter  ::  facth = 1.0d+3
  real(8),    parameter  ::  factw = 1.0d+2
  integer(4), parameter  ::  ntyz0 = 3,   nrbnd = 3
  integer(4), parameter  ::  ntyzz = ntyz0 - 1

  public  ::  sfcflx
#ifdef OPT_BODY
  public  ::  bdyflx
#endif

contains

  subroutine sfcflx(                                                  &
    &         qao,    qai,    qii,    qio,  swabs,    tsi,            &
    &         wev,    wsb,                                            &
    &        prec,   snow,   roff,   soff,                            &
    &      tauaix, tauaiy, tauaox, tauaoy,                            &
    &          ft,   ptop,   ssfc,                                    &
    &        dfdu,   dfbc,                                            &
    &           t,      a,     hi,     ti,    hsn,                    &
    &          as,    vmp,   frmp,                                    &
    &           u,      v )

    use zocdim,  only  :                                              &
    &     nxdim,   nydim,   nzdim,   ntdim,   nxydim,   nxyidm,       &
    &        nx,      ny,       nz,    nic,                           &
    &     ijstr,   ijend,    kstr,                                    &
    &        le,      ln,     lne
    use zocgrd,  only  :                                              &
    &        dz
    use zocmsk,  only  :                                              &
    &     amskt
    use zocphy,  only  :                                              &
    &    kelvin,    cdi,   dtds,  tmelt
    use utint
    use ufile
    use qckot
    use bshft

    implicit none

    real(8),    intent(inout)  ::    qao(nxydim)
    real(8),    intent(inout)  ::    qai(nxydim,nic),   qii(nxydim,nic)
    real(8),    intent(inout)  ::    qio(nxydim,nic)
    real(8),    intent(inout)  ::  swabs(nxydim)
    real(8),    intent(inout)  ::    tsi(nxydim,0:nic)
    real(8),    intent(inout)  ::    wev(nxydim),       wsb(nxydim,nic)
    real(8),    intent(inout)  ::   prec(nxydim),      snow(nxydim)
    real(8),    intent(inout)  ::   roff(nxydim),      soff(nxydim)
    real(8),    intent(inout)  :: tauaix(nxydim),    tauaiy(nxydim)
    real(8),    intent(inout)  :: tauaox(nxydim),    tauaoy(nxydim)
    real(8),    intent(inout)  ::     ft(nxydim,ntdim)
    real(8),    intent(inout)  ::   ptop(nxydim),      ssfc(nxydim)
    real(8),    intent(in)     ::      t(nxydim,nzdim,ntdim)
    real(8),    intent(in)     ::      a(nxydim,0:nic),  hi(nxydim,0:nic)
    real(8),    intent(in)     ::     ti(nxydim,0:nic), hsn(nxydim,0:nic)
    real(8),    intent(in)     ::     as(nxydim, 0:nic), vmp(nxydim, 0:nic)
    real(8),    intent(in)     ::   frmp(nxydim, 0:nic)
    real(8),    intent(in)     ::      u(nxydim,nzdim),   v(nxydim,nzdim) !! dummy
    real(8),    intent(out)    ::   dfdu(nxydim),   dfbc(nxydim)

!---- local variables
    real(8)        ::    taux(nxydim),   tauy(nxydim),   usfc(nxydim)
    real(8)        ::    tsfc(nxydim),   qsfc(nxydim),   wflx(nxydim)
    real(8)        ::    swnt(nxydim),   dwlw(nxydim),   psfc(nxydim)
    real(8)        ::    grts(nxydim),   grtb(nxydim)
    real(8)        ::   grice(nxydim),  grsnw(nxydim),  gricr(nxydim)
    real(8)        ::    grz0(nxydim,ntyz0)
    real(8)        ::  gfluxs(nxydim), tfluxs(nxydim), qfluxs(nxydim)
    real(8)        ::  wfluxs(nxydim,2)
    real(8)        ::  rflxlu(nxydim), sflxbl(nxydim)
    real(8)        ::   dgfds(nxydim),  dtfdt(nxydim),  dtfds(nxydim)
    real(8)        ::   dqfds(nxydim)
    real(8),  save ::     tmi,   dirdsn
    real(8)        ::    tisi(nxydim,0:nic)

!     [INTERNAL PARM] 
!     We need DFSNOW and DFICE for flux redistribution.

    real(8),  save ::     dfsnow                      !! diffusion coef. of snow
    real(8),  save ::     snwdmx                      !! maximum snow depth
    real(8),  save ::     epssnw                      !! minimum snow
    real(8),  save ::     snwmax                      !! ice forming snow
    real(8),  save ::     albsnw( 2, nrbnd )          !! snow albedo
    real(8),  save ::     talsnw( 2 )                 !! temp. for alb. change
    real(8),  save ::     z0snw ( ntyz0 )             !! roughness of snow
    real(8),  save ::     snwcrt                      !! snow amount for fraction=1
    real(8),  save ::     snwden                      !! snow density (kg/m**3)
                   
    real(8), save  ::     dfice                       !! diffusion coef. of sea ice
    real(8), save  ::     albice( nrbnd )             !! sea ice albedo
    real(8), save  ::     z0ice ( ntyz0 )             !! roughness of sea ice
    real(8), save  ::     siccrt                      !! ice amount for conc.=1
    real(8), save  ::     sicden                      !! ice density (kg/m**3)

    real(8), save  ::      si
    integer(4)     ::      ij,      l,      n,     nn
    integer(4)     ::   ifpar,  jfpar,  istat
    logical, save  ::  ofirst
    namelist  /nmsnow/  dfsnow, snwdmx, epssnw, snwmax,               &
    &                   albsnw, talsnw,  z0snw,                       &
    &                   snwcrt, snwden  
    namelist  /nmice/    dfice, albice,  z0ice, siccrt, sicden
    namelist /nmislt/ si
    data   dfsnow / 0.4d0          /
    data   snwdmx / 5.d0           /
    data   epssnw / 1.d-8          /        
    data   snwmax / 1000.d0        /
    data   albsnw /  0.75d0, 0.5d0,                                  &
    &                0.75d0, 0.5d0,                                  &
    &                0.0d0 , 0.0d0      /
    data   talsnw / 258.15d0, 273.15d0 /
    data   z0snw  / 1.d-2, 1.d-3, 1.d-3 /
    data   snwcrt / 100.d0         /
    data   snwden / 400.d0         /
           
    data   dfice  / 2.d0           /
    data   albice / 0.5d0, 0.5d0, 0.05d0 /
    data   z0ice  / 2.d-2, 2.d-3, 2.d-3 /
    data   siccrt /  300.d0        /
    data   sicden / 1000.d0        /
    data       si / 5.d0 /
    data   ofirst / .true. /

    if ( ofirst ) then

       call rewnml(ifpar, jfpar)
       write(jfpar, *) '*** SFCFLX ***'
       READ_NAMELIST( nmislt )
       READ_NAMELIST( nmsnow )
       READ_NAMELIST( nmice  )

       tmi = dtds * si         
       dirdsn = dfice / dfsnow
       do ij = 1, nxydim
          tsfc(ij) = 300.d0
          psfc(ij) = 1.d5
       end do
       ofirst = .false.

    end if

    tisi(1:nxydim,0:nic) = 0.d0

    do l = 3, ntdim

       n = 11 + (l - 3) * 2
       nn = n + 1
       call tmintp(  tsfc,      n)
       call tmintp(  qsfc,     nn)
       do ij = ijstr, ijend
          ft(ij, l) = qsfc(ij) * (tsfc(ij) - t(ij, kstr, l)) *        &
    &                  dz(ij, kstr) * amskt(ij, kstr)
       end do
    end do
      
    call tmintp(  taux,      1)
    call tmintp(  tauy,      2)
    call tmintp(  usfc,      3)
    call tmintp(  tsfc,      4)
    call tmintp(  qsfc,      5)
    call tmintp(  wflx,      6)
    call tmintp(  swnt,      7)
    call tmintp(  dwlw,      8)
    call tmintp(  psfc,      9)
#ifdef OPT_SRST
    call tmintp(  ssfc,     10)
#endif

!! 2021.05.31: Now dfdu and dfbc are dummy fluxes in OGCM.
  dfdu(:) = 0.0d0
  dfbc(:) = 0.0d0

#ifdef OPT_TRIPOLE
    call shift2(                                                      &
    &               taux,   tauy,                                     &
    &              nxdim,  nydim,      1,                             &
    &             -1.D0,      0,      0 )
    call shift1(                                                      &
    &               psfc,                                             &
    &              nxdim,  nydim,      1,                             &
    &               1.d0,      0,      0 )
#else
    call shift3(                                                      &
    &              taux,   tauy,   psfc,                              &
    &             nxdim,  nydim,      1)
#endif

    do ij = ijstr, ijend
       tauaix(ij) = (  taux(ij)    + taux(ij+le)                      &
    &                + taux(ij+ln) + taux(ij+lne)) *                  &
    &               0.25d0 * factm
       tauaiy(ij) = (  tauy(ij)    + tauy(ij+le)                      &
    &                + tauy(ij+ln) + tauy(ij+lne)) *                  &
    &               0.25d0 * factm
       tauaox(ij) = (  taux(ij)    + taux(ij+le)                      &
    &                + taux(ij+ln) + taux(ij+lne)) *                  & 
    &               0.25d0 * factm
       tauaoy(ij) = (  tauy(ij)    + tauy(ij+le)                      &
    &                + tauy(ij+ln) + tauy(ij+lne)) *                  &
    &               0.25d0 * factm
       snow(ij) = 0.d0
       prec(ij) = - wflx(ij) * factw
       soff(ij) = 0.d0
       roff(ij) = 0.d0
       swabs(ij) = swnt(ij) * facth * a(ij, 0)
    end do
    
    do l = 0, nic

       if (l > 0) then
          do ij = ijstr, ijend
             grts (ij) = tsi(ij, l) + kelvin
             grtb (ij) =  ti(ij, l) + kelvin
             grice(ij) =  hi(ij, l) * 5.d-3
             grsnw(ij) = hsn(ij, l) * 1.d-2
             gricr(ij) = 1.d0
          end do
       else
          do ij = ijstr, ijend
             grts (ij) = t(ij, kstr, 1) + kelvin
             grtb (ij) = t(ij, kstr, 1) + kelvin
             grice(ij) = 0.d0
             grsnw(ij) = 0.d0
             gricr(ij) = 0.d0
          end do
       end if

       call ocnbcs_omip(                                              &
    &               grz0, gfluxs,  dgfds,                             &
    &               grts,   grtb,  grice,  grsnw,  gricr,             &
    &               usfc)
       call sfcflx_omip(                                              &
    &            tfluxs, qfluxs,                                      &
    &             dtfdt,  dtfds,  dqfds,                              &
    &              usfc,   tsfc,   qsfc,   psfc,                      &
    &              grts,   grz0)
       call ocnslv_omip(                                              &
    &              grts, gfluxs, tfluxs, qfluxs,                      &
    &            wfluxs, rflxlu, sflxbl,                              &
    &             dtfds,  dqfds,  dgfds,                              &
    &              swnt,   dwlw,                                      &
    &             gricr,    tmi)

       if ( l > 0 ) then
          do ij = ijstr, ijend
             qai(ij, l) = (gfluxs(ij) - sflxbl(ij)) * facth
             qii(ij, l) = gfluxs(ij) * facth
             qio(ij, l) = (t(ij, kstr, 1) - ti(ij, l)) *             &
    &                      cdi / hi(ij, l) * 2.d0
             tsi(ij, l) = min(grts(ij), tmelt+tmi) - kelvin
             tisi(ij, l) = tsi(ij, l) +                              &
    &                       (ti(ij, l) - tsi(ij, l)) *               &
    &                       dirdsn * hsn(ij, l) /                    &
    &                       (0.5d0 * hi(ij, l) + dirdsn * hsn(ij, l))
             wsb(ij, l) = 0.d0 ! (sublimation speed) x (fraction)
          end do
       else
          do ij = ijstr, ijend
             qao(ij) = (gfluxs(ij) + swnt(ij)) * facth
             wev(ij) = 0.d0 ! (evaporation speed) x (fraction)
          end do
       end if
    end do
    
!     output section for CMIP5
!     TISI: temperature at interface between sea ice and snow
!           unit [C]
    call chekin(tisi,   'TISI',                             &
         &   'temperature at ice-snow interface', 'degC',   &
         &            nx,     ny,    nic, nxyidm, 'OCICET')

    do ij = 1, nxydim
       ptop(ij) = psfc(ij) * factm
    end do
    
  end subroutine sfcflx

! *********************************************************************

  subroutine ocnslv_omip(                                             &
    &       gdts,   gfluxs,   tfluxs,   qfluxs,                       &
    &     wfluxs,   rflxlu,   sflxbl,                                 &
    &      dtfds,    dqfds,    dgfds,                                 &
    &     rflxsd,   rflxld,                                           &
    &      gricr ,   tmi  ) 

    use zocdim,  only  :                                              &
    &    nxydim,   ntdim,     nic,                                    &
    &     ijstr,   ijend
    use zocphy,  only  :                                              &
    &    tmelt,    emelt,      el,    stb
    use ufile

    implicit none

    real(8),    intent(inout)  ::    gdts  ( nxydim )    !! skin temperature
    real(8),    intent(inout)  ::    gfluxs( nxydim )    !! soil heat flux
    real(8),    intent(inout)  ::    tfluxs( nxydim )    !! flux of T
    real(8),    intent(inout)  ::    qfluxs( nxydim )    !! flux of q
    real(8),    intent(inout)  ::    wfluxs( nxydim, 2 ) !! soil water flux
    real(8),    intent(inout)  ::    rflxlu( nxydim )    !! upward long wave
    real(8),    intent(inout)  ::    sflxbl( nxydim )    !! flux balance
    real(8),    intent(in)     ::    dtfds ( nxydim )    !! dH/dtg
    real(8),    intent(in)     ::    dqfds ( nxydim )    !! dE/dtg
    real(8),    intent(in)     ::    dgfds ( nxydim )    !! dG/dtg
    real(8),    intent(in)     ::    rflxsd( nxydim )    !! down. SW rad.
    real(8),    intent(in)     ::    rflxld( nxydim )    !! down. LW rad.
    real(8),    intent(in)     ::    gricr ( nxydim )    !! snow/ice ratio
    real(8),    intent(in)     ::    tmi ! sea ice (and snow) melting temperature in c

!   [INTERNAL WORK] 
    real(8)     ::   esub,    stg,  drfds
    real(8)     ::  sflux, gsflux, dgsfds
    real(8)     :: gfluxf
    real(8)     :: sflxbi, dsbdsi,    dti,  evapi,  gfluxi
    real(8)     ::     ff,     fi,    dtx
    integer(4)  ::     ij,  istat
    integer(4)  ::  ifpar, jfpar
!   [INTERNAL SAVE] 
    real(8),   save  ::    emis  !! longwave emissivity (longwave co-albedo)
    real(8),   save  ::  flwnet
    logical,   save  ::  ofirst,  olwnet
    namelist /nmlwem/ emis, olwnet
    data olwnet, ofirst / .false., .true. /
    data emis / 1.d0 /

    if ( ofirst ) then
       READ_NAMELIST( nmlwem )

       ofirst = .false.
       write (jfpar, *) ' @@@ OCNSLV: OCEAN HEAT BALANCE 98/06/30'
       if (olwnet) then
          flwnet = 0.d0
       else
          flwnet = 1.d0
       end if

    end if

    esub = el + emelt
    do ij = ijstr, ijend

       stg            = emis*stb*gdts( ij )**4
       rflxlu( ij )   = (stg + ( 1.d0-emis )*rflxld( ij )) * flwnet
       drfds          = 4.d0*stg/gdts( ij ) * flwnet

       sflux          = tfluxs( ij )                                  &
    &                 + rflxlu( ij ) - rflxld( ij )                   &
    &                 - rflxsd( ij )
       gsflux         = gfluxs( ij ) - sflux
       dgsfds         = dgfds ( ij )                                  &
    &                 + dtfds ( ij ) + drfds

!       < ice free ocean surface >
       gfluxf         = sflux + qfluxs( ij )*el

!       < ice/snow surface >
       sflxbi         = gsflux - qfluxs( ij )*esub
       dsbdsi         = dgsfds + dqfds ( ij )*esub
       dti            = sflxbi/dsbdsi
       dti            = min( dti, tmelt+tmi-gdts( ij ) )
       evapi          = qfluxs( ij ) + dqfds( ij )*dti
       gfluxi         = gfluxs( ij ) - dgfds( ij )*dti
       sflxbi         = sflxbi       - dsbdsi     *dti

!      < Ts change >
       ff             = 1.d0 - gricr( ij )
       fi             = gricr( ij )
       dtx            = fi * dti
       gdts  ( ij )   = gdts( ij ) + dtx

!      < fluxs >
       tfluxs( ij )   = tfluxs( ij ) + dtfds( ij )*dtx
       rflxlu( ij )   = rflxlu( ij ) + drfds      *dtx
       qfluxs( ij )   = ff * qfluxs( ij ) + fi * evapi
       gfluxs( ij )   = ff * gfluxf       + fi * gfluxi 
       wfluxs( ij,1 ) = ff * qfluxs( ij )
       wfluxs( ij,2 ) = fi * evapi
       sflxbl( ij )   = fi * sflxbi

    end do

  end subroutine ocnslv_omip


! =====================================================================

  subroutine  ocnbcs_omip(                                            &
    &      grz0  , fogflx, dgfds ,                                    &
    &      grts  , grtb  , grice , grsnw , gricr ,                    &
    &      usfc  )

    use zocdim,  only  :                                              &
    &    nxydim,   ntdim,     nic,                                    &
    &     ijstr,   ijend
    use ufile

    implicit none

    real(8),    intent(inout)  ::    grz0  ( nxydim, ntyz0  )    !! surface roughness 
    real(8),    intent(inout)  ::    fogflx( nxydim )            !! heat flux
    real(8),    intent(inout)  ::    dgfds ( nxydim )            !! dG/dTs
    real(8),    intent(in)     ::    grts  ( nxydim )            !! skin temperature
    real(8),    intent(in)     ::    grtb  ( nxydim )            !! ice base temp.
    real(8),    intent(in)     ::    grice ( nxydim )            !! sea ice
    real(8),    intent(in)     ::    grsnw ( nxydim )            !! snow smount
    real(8),    intent(in)     ::    gricr ( nxydim )            !! ice fraction
    real(8),    intent(in)     ::    usfc  ( nxydim )            !! sfc wind

!----- local variables
    real(8)            ::   grsnr(nxydim)      !! snow fraction
    real(8),     save  ::  dfsnow              !! diffusion coef. of snow
    real(8),     save  ::  snwdmx              !! maximum snow depth
    real(8),     save  ::  epssnw              !! minimum snow
    real(8),     save  ::  snwmax              !! ice forming snow
    real(8),     save  ::  albsnw( 2, nrbnd )  !! snow albedo
    real(8),     save  ::  talsnw( 2 )         !! temp. for alb. change
    real(8),     save  ::  z0snw ( ntyz0 )     !! roughness of snow
    real(8),     save  ::  snwcrt              !! snow amount for fraction=1
    real(8),     save  ::  snwden              !! snow density (kg/m**3)
    real(8),     save  ::   dfice              !! diffusion coef. of sea ice
    real(8),     save  ::  albice( nrbnd )     !! sea ice albedo
    real(8),     save  ::   z0ice( ntyz0 )     !! roughness of sea ice
    real(8),     save  ::  siccrt              !! ice amount for conc.=1
    real(8),     save  ::  sicden              !! ice density (kg/m**3)
    
    real(8),     save  ::   z0fct              !! heat z0/moumentum z0
    real(8),     save  ::   z0min              !! minimum z0
    
    real(8),     save  ::   dzocn              !! depth of ML ocean
    real(8),     save  ::   dfocn              !! ocean dg/dts
    real(8),     save  ::   alblo              !! LW albedo (1-emis)

    real(8)            ::    z00,    dz0
    real(8)            ::   dfgt,   dfgx
    integer(4)         ::     ij,      m,     ifg
    integer(4)         ::  istat,  ifpar,   jfpar
    logical,     save  :: ofirst
    data       dfsnow         / 0.4d0          /
    data       snwdmx         / 5.d0           /
    data       epssnw         / 1.d-8          /        
    data       snwmax         / 1000.d0        /
    data       albsnw         / 0.75, 0.5,                            &
    &                           0.75, 0.5,                            &
    &                           0.0 , 0.0      /
    data       talsnw         / 258.15, 273.15 /
    data       z0snw          / 1.d-2, 1.d-3, 1.d-3 /
    data       snwcrt         / 100.d0         /
    data       snwden         / 400.d0         /

    data       dfice          / 2.d0           /
    data       albice         / 0.5, 0.5, 0.05 /
    data       z0ice          / 2.d-2, 2.d-3, 2.d-3 /
    data       siccrt         /  300.d0        /
    data       sicden         / 1000.d0        /

    data       z0fct          / 0.1d0          /
    data       z0min          / 1.d-6          /

    data       dzocn          / 50.d0          /
    data       dfocn          / 1.d10          /
    data       alblo          / 0.05           /
    data  ofirst / .true. /
    namelist  /nmsnow/ dfsnow, snwdmx, epssnw, snwmax,                &
    &                  albsnw, talsnw, z0snw ,                        &
    &                  snwcrt, snwden  
    namelist  /nmice/  dfice,  albice, z0ice , siccrt, sicden
    namelist  /nmz0/   z0fct,  z0min
    namelist  /nmocn/  dzocn,  dfocn , alblo

    if ( ofirst ) then

       call rewnml(ifpar, jfpar)
       write (jfpar, *) ' @@@ OCNBCS: OCEAN SURFACE BC 98/07/29'
       ofirst = .false.
       READ_NAMELIST( nmsnow )
       READ_NAMELIST( nmice  )
       READ_NAMELIST( nmz0   )
       READ_NAMELIST( nmocn  )
    endif

    do ij = ijstr, ijend
       if ( grsnw(ij) > 0.d0 ) then
          grsnr(ij) = 1.d0
       else
          grsnr(ij) = 0.d0
       end if
    end do

    call seaz0f_omip( grz0(1,1) , grz0(1,2), grz0(1,3), usfc )

    do m = 1, ntyz0

       do ij = ijstr, ijend

          if ( gricr( ij ) > 0.d0 ) then
             dz0   = z0ice( m )   - grz0( ij,m )
             z00   = grz0( ij,m ) + dz0*gricr( ij )
             dz0   = z0snw( m )   - z00
             grz0( ij,m ) = z00   + dz0*grsnr( ij )
          end if
          grz0( ij,m ) = max( grz0( ij,m ), z0min )

       end do

       do ij = ijstr, ijend

          if ( gricr(ij) > 0.d0 ) then
             dfgt = dfice / grice(ij)
             dfgx = dfice * dfsnow                                    &
    &             / (dfice * grsnw(ij) + dfsnow * grice(ij))
             dfgt = dfgt * (1.d0 - grsnr(ij)) + dfgx * grsnr(ij)
             fogflx(ij) = dfgt * (grtb(ij) - grts(ij))
             dgfds (ij) = dfgt
          else
             fogflx(ij) = 0.d0
             dgfds (ij) = dfocn
          endif

       end do
       
    end do

  end subroutine ocnbcs_omip

! *********************************************************************

  subroutine seaz0f_omip(                                             &
    &    grz0m , grz0h , grz0e ,                                      &
    &    usfc   )

    use zocdim,  only  :                                              &
    &    nxydim,   ntdim,     nic,                                    &
    &     ijstr,   ijend
    use zocphy,  only  :                                              &
    &      grav
    use ufile

    implicit none
!
!      roughness of sea surface (Miller et al., 1992)
!
    real(8),    intent(inout)  ::  grz0m( nxydim )   !! surface roughness (V)
    real(8),    intent(inout)  ::  grz0h( nxydim )   !! surface roughness (T)
    real(8),    intent(inout)  ::  grz0e( nxydim )   !! surface roughness (q)
    real(8),    intent(in)     ::   usfc( nxydim )   !! sfc wind speed

    real(8),     save  ::     z0m0                   !! base
    real(8),     save  ::     z0mr                   !! rough factor
    real(8),     save  ::     z0ms                   !! smooth factor
    real(8),     save  ::     z0h0                   !! base
    real(8),     save  ::     z0hr                   !! rough factor
    real(8),     save  ::     z0hs                   !! smooth factor
    real(8),     save  ::     z0e0                   !! base
    real(8),     save  ::     z0er                   !! rough factor
    real(8),     save  ::     z0es                   !! smooth factor
    real(8),     save  ::     visair                 !! kinematic viscosity 
    real(8),     save  ::     cm0                    !! bulk coef for ustar
    real(8),     save  ::     ustrmn                 !! min(u*)
    real(8),     save  ::     z0mmin                 !! minimum
    real(8),     save  ::     z0hmin                 !! minimum
    real(8),     save  ::     z0emin                 !! minimum
    logical,     save  ::     ofirst
    real(8)            ::  ustar
    integer(4)         ::     ij
    integer(4)         ::  istat,  ifpar,   jfpar
    namelist  /nmseaz/                                                &
    &           z0m0, z0mr, z0ms,                                     &
    &           z0h0, z0hr, z0hs,                                     &
    &           z0e0, z0er, z0es,                                     &
    &           visair, cm0, ustrmn,                                  &
    &           z0mmin, z0hmin, z0emin
    data       z0m0, z0mr, z0ms / 0.,     0.018, 0.11  /
    data       z0h0, z0hr, z0hs / 1.4e-5, 0.   , 0.4   /
    data       z0e0, z0er, z0es / 1.3e-4, 0.   , 0.62  /
    data       visair           / 1.5e-5  / 
    data       cm0              / 1.d-3   /
    data       ustrmn           / 1.d-3   /
    data       z0mmin, z0hmin, z0emin  / 3*1.e-5  /        
    data       ofirst           / .true. /

    if ( ofirst ) then

       ofirst = .false.
       call rewnml(ifpar, jfpar)
       write(jfpar, *) ' SEA SURFACE Z0 (Miller et al.) 98/06/19'
       READ_NAMELIST( nmseaz )
    end if

    do ij = ijstr, ijend

       ustar       = cm0 * usfc(ij) * usfc(ij)
       ustar       = max(sqrt(ustar), ustrmn)

       grz0m( ij ) = z0m0 + z0mr*ustar**2/grav + z0ms*visair/ustar
       grz0h( ij ) = z0h0 + z0hr*ustar**2/grav + z0hs*visair/ustar
       grz0e( ij ) = z0e0 + z0er*ustar**2/grav + z0es*visair/ustar

       grz0m( ij ) = max( grz0m( ij ), z0mmin )            
       grz0h( ij ) = max( grz0h( ij ), z0hmin )            
       grz0e( ij ) = max( grz0e( ij ), z0emin )            

    end do

  end subroutine seaz0f_omip

! *********************************************************************

  subroutine sfcflx_omip(                                             &
    &           TFLUXS, QFLUXS,                                       &
    &           DTFDT , DTFDS , DQFDS ,                               &
    &           USFC  , GDTA  , GDQA  , GDPS  ,                       &
    &           GDTS  , GRZ0  ) 

    use zocdim,  only  :                                              &
    &    nxydim,   ntdim,     nic,                                    &
    &     ijstr,   ijend
    use zocphy,  only  :                                              &
    &        cp,    rvap,    rair,    grav
    use zqsat
    use ufile

    implicit none

    real(8),    intent(inout)  ::  tfluxs( nxydim )         !! flux of T
    real(8),    intent(inout)  ::  qfluxs( nxydim )         !! flux of q
    real(8),    intent(inout)  ::  dtfdt ( nxydim )         !! -dH/dta
    real(8),    intent(inout)  ::  dtfds ( nxydim )         !! dH/dts
    real(8),    intent(inout)  ::  dqfds ( nxydim )         !! dLE/dts
    real(8),    intent(in)     ::  usfc  ( nxydim )         !! westerly u
    real(8),    intent(in)     ::  gdta  ( nxydim )         !! temperature T
    real(8),    intent(in)     ::  gdqa  ( nxydim )         !! humidity q
    real(8),    intent(in)     ::  gdps  ( nxydim )
    real(8),    intent(in)     ::  gdts  ( nxydim )         !! surface skin temp
    real(8),    intent(in)     ::  grz0  ( nxydim, ntyz0 )  !! surface roughness

!----- local variables
    real(8)          ::   dufdu ( nxydim )      !! -d(tau)/du
    real(8)          ::   dqfdq ( nxydim )      !! -dLE/dqa
    real(8),    save ::  akappa,     za
    real(8)          ::     exi,    rho,   qsat,  dqsat
    integer(4)       ::      ij,  ifpar,  jfpar,  istat
    logical,    save ::  ofirst 
    data ofirst / .true. /
    data za / 2.d0 /
    namelist /nmsair/ za

    if ( ofirst ) then
       READ_NAMELIST( nmsair )
       ofirst = .false.
       write(jfpar, *) ' @@@ PSFCM: SURFACE FLUX 98/06/19'
       akappa = rair/cp
    end if

    call blkcof_omip                                  &
    &         ( dufdu , dtfds , dqfdq ,               &  !! cmv, chv, cev
    &           usfc  , gdta  , gdqa  , gdts  , gdps  , grz0  )

    do ij = ijstr, ijend

       rho = gdps(ij) / (rair + (rvap - rair) * gdqa(ij)) / gdta(ij)
       exi = (gdps(ij) / (gdps(ij) - rho * grav * za))**akappa

       qsat         = fqsat ( gdts( ij ), gdps(ij) )
       dqsat        = fdqsat( gdts( ij ), qsat )
     
       dufdu ( ij ) =    rho*dufdu( ij )
       dtfds ( ij ) = cp*rho*dtfds( ij )
       dqfdq ( ij ) =    rho*dqfdq( ij )

       dtfdt ( ij ) =   dtfds( ij )*exi
       dqfds ( ij ) =   dqfdq( ij )*dqsat

       tfluxs( ij ) =   dtfds( ij )*( gdts( ij )-gdta( ij )*exi )
       qfluxs( ij ) =   dqfdq( ij )*( qsat      -gdqa( ij )     )

    end do

  end subroutine sfcflx_omip

! *********************************************************************

  subroutine blkcof_omip                                              &
    &         ( cmv   , chv   , cev   ,                               &
    &           usfc  , gdta  , gdqa  , gdts  , gdps  , grz0  )

    use zocdim,  only  :                                              &
    &    nxydim,   ntdim,     nic,                                    &
    &     ijstr,   ijend
    use zocphy,  only  :                                              &
    &        cp,    rvap,    rair,    grav
    use ufile

    implicit none

    real(8),    intent(inout)  ::   cmv( nxydim )        !! bulk transfer coeff.: u
    real(8),    intent(inout)  ::   chv( nxydim )        !! bulk transfer coeff.: t
    real(8),    intent(inout)  ::   cev( nxydim )        !! bulk transfer coeff.: q
 
    real(8),    intent(in)     ::  usfc( nxydim )        !! u wind speed
    real(8),    intent(in)     ::  gdta( nxydim )        !! temperature t
    real(8),    intent(in)     ::  gdqa( nxydim )
    real(8),    intent(in)     ::  gdts( nxydim  )       !! surface skin temp
    real(8),    intent(in)     ::  gdps( nxydim )
    real(8),    intent(in)     ::  grz0( nxydim, ntyz0 ) !! surface roughness

!---- local variables
    real(8),      save  ::  akappa
    real(8),      save  ::   cmmin,  chmin,  cemin
    real(8),      save  ::   cmmax,  chmax,  cemax
    real(8),      save  ::  usminm, usminh, usmine
    real(8),      save  ::  usmaxm, usmaxh, usmaxe
    real(8),      save  ::      za
    real(8)             ::     cl0,    cl1
    real(8)             ::     exi,    rho,   vabs
    integer(4)          ::      ij
    integer(4)          ::   ifpar,  jfpar,  istat
    logical,      save  ::  ofirst
    namelist /nmsair/ za
    namelist  /nmsfcl/ cmmin , chmin , cemin ,                        &
    &                  cmmax , chmax , cemax ,                        &
    &                  usminm, usmaxm,                                &
    &                  usminh, usmaxh,                                &
    &                  usmine, usmaxe
    data      za / 2.d0 /
    data   cmmin /  1.e-5 /     !! min. bulk coef. of u
    data   chmin /  1.e-5 /     !! min. bulk coef. of t
    data   cemin /  1.e-5 /     !! min. bulk coef. of q
    data   cmmax /   1.   /     !! max. bulk coef. of u
    data   chmax /   1.   /     !! max. bulk coef. of t
    data   cemax /   1.   /     !! max. bulk coef. of q
    data  usminm / 2.0    /     !! min. wind vel. for v
    data  usminh / 2.0    /     !! min. wind vel. for t
    data  usmine / 2.0    /     !! min. wind vel. for q
    data  usmaxm / 1000.  /     !! max wind vel. for v
    data  usmaxh / 1000.  /     !! max. wind vel. for heat
    data  usmaxe / 1000.  /     !! max wind vel. for q
    data  ofirst / .true. /
      
    if ( ofirst ) then
       call rewnml(ifpar, jfpar)
       write(jfpar, *) ' @@@ PSFCL: SURFACE BULK COEF. Kara' 
       READ_NAMELIST( nmsair )
       ofirst = .false.

       READ_NAMELIST( nmsfcl )
       akappa = rair / cp
    end if

    do ij = ijstr, ijend

       rho = gdps(ij) / (rair + (rvap - rair) * gdqa(ij)) / gdta(ij)
       exi = (gdps(ij) / (gdps(ij) - rho * grav * za))**akappa

       vabs      = usfc(ij)
       vabs      = min( max( vabs, usmine ), usmaxe )

       cl0 = 1.d-3 * (   0.994d0                                      &
    &                  + 0.061d0 * vabs                               &
    &                  - 0.001d0 * vabs * vabs)
       cl1 = 1.d-3 * ( - 0.020d0                                      &
    &                  + 0.691d0 * 1.d0 / vabs                        &
    &                  - 0.817d0 * 1.d0 / vabs / vabs)
       cev(ij) = cl0 + cl1 * ( gdts( ij )-gdta( ij )*exi ) 
       chv(ij) = 0.96d0 * cev(ij)

    end do

    do ij = ijstr, ijend
       cmv( ij ) = max( min( cmv( ij ), cmmax ), cmmin )
       chv( ij ) = max( min( chv( ij ), chmax ), chmin )
       cev( ij ) = max( min( cev( ij ), cemax ), cemin )
    end do

    do ij = ijstr, ijend
       vabs      = usfc(ij)
       cmv( ij ) = cmv( ij ) * min( max( vabs, usminm ), usmaxm )
       chv( ij ) = chv( ij ) * min( max( vabs, usminh ), usmaxh )
       cev( ij ) = cev( ij ) * min( max( vabs, usmine ), usmaxe )
    end do
    
  end subroutine blkcof_omip

! --- information -----------------------------------------------------
!
!  Estimate the surface forcing fluxes from the data given in the
! boundary condition files
!
!  HISTORY
!     '99.08.25  H.Hasumi
!     '02.06.02  H.Hasumi: tracer dimension
!     '12.09.04  H.Tatebe: Rewrite in F95 format
! ---------------------------------------------------------------------

#ifdef OPT_BODY

  subroutine bdyflx(  tq,  t )
    
    use zocdim,   only  :                                             &
    &      nxydim,   nzdim,   ntdim,                                  &
    &      ijtstr,  ijtend,    kstr,   kend
    use zocmsk,   only  :                                             &
    &       amskt
    use utint

    implicit none

    real(8),    intent(inout)  ::  tq(nxydim, nzdim, ntdim)
    real(8),    intent(inout)  ::  t (nxydim, nzdim, ntdim)


    real(8),    save       ::    tbdy(nxydim, nzdim, ntdim) = 0.d0
    real(8),    save       ::    tdmb(nxydim, nzdim, ntdim) = 0.d0
    integer(4)   ::      ij,      k,      l,      n

    do l = 1, 2
       n = l
       call tmintb(  tbdy(1, 1, l),  n  )
    end do

    do l = 1, 2
       n = l + 2
       call tmintb(  tdmb(1, 1, l),  n  )
    end do

    do l = 1, 2
       do k = kstr, kend
          do ij = ijtstr, ijtend
             tq(ij, k, l) = tdmb(ij, k, l) *                          &
    &                      (tbdy(ij, k, l) - t(ij, k, l)) * amskt(ij, k)
          end do
       end do
    end do

  end subroutine bdyflx

#endif

end module sfcng


