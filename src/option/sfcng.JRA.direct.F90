module sfcng

! --- information -----------------------------------------------------
!
!  Estimate the surface forcing fluxes from the data given in the
! boundary condition files
!
!  HISTORY
!     '08.06.30  Y.Komuro: from COCO4.2
!     '08.08.06  Y.Komuro: for COCO4.4; initial/final processing added
!     '08.08.26  Y.Komuro: diagnosing evapolation
!     '08.09.02  Y.Komuro: diagnosing wind stress
!     '08.09.08  Y.Komuro: bug(s) fix in BLKCOF_CORE
!     '08.09.16  Y.Komuro: latitude-dependent albedo
!     '08.09.19  Y.Komuro: bug fix in BLKCOF_CORE
!     '09.05.25  Y.Komuro: CMIP5 output code included
!     '09.09.26  Y.Komuro: bug fix
!     '10.04.14  M.kurogi
!     '10.04.14  M.Kurogi: for tripolar grid
!     '12.10.11  Y.Komuro: for COCO5.0
!     '13.02.13  Y.Komuro: remove non-parallel code 
!     '21.05.26  Y.Komuro: snow aging & meltpond parametrization
!
! ---------------------------------------------------------------------

  use zocdim, only: &
    &     nx,     ny,  nxdim, nydim, &
    & nxydim, nxyidm,  nzdim,   kstr,   kend,  ntdim,    nic, &
    &  ijstr,  ijend, ijtstr, ijtend, &
    &     le,     lw,     ln,     ls,    lne,    lsw
  use zocgrd, only: &
    &     dz,    cor,  rangt
  use zocmsk, only: &
    &  amskt
  use zocphy, only: &
    &     cp,   grav,   rair,     el,  emelt,   rvap, &
    &  dwatr,    es0,    stb,  tmelt,  tqice,   epsv, &
    &  ckarm, kelvin,    cdi,   dtds,   rhoo,   rhos

  implicit none

  integer, parameter :: ntyz0 = 3, nrbnd = 3

! [namelist parameters] 
! namelist snow
  real(8), save :: dfsnow = 0.4d0              !! diffusion coef. of snow
  real(8), save :: snwdmx = 5.0d0              !! maximum snow depth
  real(8), save :: epssnw = 1.0d-8             !! minimum snow
  real(8), save :: snwmax = 1000.0d0           !! ice forming snow
  real(8), save :: albsnw( 2, nrbnd ) = &      !! snow albedo
    &                 reshape ( (/ 0.75d0, 0.5d0, &
    &                              0.75d0, 0.5d0,  &
    &                              0.0d0 , 0.0d0 /), &
    &                           shape = (/ 2, nrbnd /) )     
  real(8), save :: talsnw( 2 ) = &             !! temp. for alb. change
    &                 (/ 258.15d0, 273.15d0 /)                
  real(8), save :: z0snw ( ntyz0 ) = &         !! roughness of snow
    &                 (/ 1.d-2, 1.d-3, 1.d-3 /)
  real(8), save :: snwcrt = 100.0d0            !! snow amount for fraction=1
  real(8), save :: snwden = 400.0d0            !! snow density (kg/m**3)
! namelist nmice 
  real(8), save :: dfice = 2.0d0               !! diffusion coef. of sea ice
  real(8), save :: albice( nrbnd ) = &         !! sea ice albedo
    &                 (/ 0.5d0, 0.5d0, 0.05d0 /)
  real(8), save :: z0ice ( ntyz0 ) = &         !! roughness of sea ice
    &                 (/ 2.0d-2, 2.0d-3, 2.0d-3 /)
  real(8), save :: siccrt = 300.0d0            !! ice amount for conc.=1
  real(8), save :: sicden = 1000.0d0           !! ice density (kg/m**3)
! namelist nmswbr
  real(8), save :: rvis = 0.53d0   !! ratio of visible band in SW
  real(8), save :: rnir = 0.47d0   !! ratio of near-IR band in SW
  real(8), save :: rir = 0.0d0   !! ratio of IR (longwave) band in SW
! namelist nmsfrc
  logical, save :: oasfrc = .false. !! true if CCSM type snow fraction
  real(8), save :: alsdpt = 1.0d0
          !! threshold thickness of snow-covered ice in trad. scheme, cm
  real(8), save :: alspat = 2.0d0
                      !! snowpatch depth for CCSM type snow fraction, cm
! namelist nmsage
  logical, save :: osage = .false. !! true if snow aging param. is used
  real(8), save :: alssif( nrbnd ) = &  !! albedo of fresh snow on sea ice
    &                 (/ 0.75d0, 0.75d0, 0.0d0 /) 
  real(8), save :: alssio( nrbnd ) = &  !! albedo of old snow on sea ice
    &                 (/ 0.5d0, 0.5d0, 0.0d0 /) 
  real(8), save :: alfmax = 0.999d0 !! maximum value for albfct
  real(8), save :: snrfrs = 1.0d0 !! snowfall for refreshing snow surface [cm]
  real(8), save :: ftage = 5.0d3  !! aging factor, for agefr1
  real(8), save :: tauage = 2.0d6 !! time scale for aging
  real(8), save :: adirt0 = 0.3d0 !! dirt factor for agefr3, normal place
  real(8), save :: adirtc = 0.01d0 !! dirt factor for agefr3, clean place
  real(8), save :: adirts = 0.1d0 !! dirt factor for agefr3, coeff. for dscppm
  real(8), save :: adirtm = 1.0d0 !! dirt factor for agefr3, maximum
  real(8), save :: drsmax = 0.1d0 !! maximum ratio of dust to snow
  logical, save :: oadst = .false.
                      !! using dsdx/dsbx for aging insted of adirt0/adirtc
! namelist nmmpnd
  integer, save :: impnd = 0 !! 0: melt pond (MP) parametrization not used
                             !! 1: Holland et al. (2012) MP param.
                             !! 2: Hunke et al. (2013) MP param.
  real(8), save :: hminmp = 10.0d0 !! min. ice thickness for keeping MP [cm]
  real(8), save :: rtdpmp = 80.0d0 !! ratio of melt pond depth to frmp [cm/1]
  real(8), save :: rtmxmp = 0.9d0  !! max. ratio of MP depth to ice thickness
  !! The default dpscl in CICE is 1.0, but we set it to 0.1.
  !! (maybe due to slightly different implementation?)
  real(8), save :: dpscl = 0.1d0  !! permiability scale parameter [ND]
  real(8), save :: rmpcmn(0:2) = & !! minimum water catching rate of MP 
    &                 (/ 0.0d0, 0.15d0, 0.15d0 /) 
  real(8), save :: rmpcmx(0:2) = & !! maximum water catching rate of MP 
    &                 (/ 0.0d0, 0.7d0, 0.85d0 /) 
  real(8), save :: cmpfrz = 3.d-6 !! constant for melt pond freeze-up rate
  real(8), save :: tmpfrz = -2.0d0 !! ref. t for melt pond freeze-up [c]
  real(8), save :: albmpd( nrbnd ) = &  !! deep melt pond shortwave albedo
    &                 (/ 0.4d0, 0.1d0, 0.0d0 /)
  real(8), save :: almpdp(2) = &   !! MP sw albedo, depth dependency [cm]
    &                 (/ 0.5d0, 20.0d0 /)
  real(8), save :: frmpmn = 1.0d-14  !! empirical limiter for frmpx [ND]
  real(8), save :: vmpmin = 1.0d-12  !! empirical limiter for vmpx [cm]

! namelist nmislt
  real(8), save ::     si = 5.0d0              !! sea-ice salinity (psu)
! namelist nmlwem
  real(8), save :: emislo = 1.0d0   !! ocn. LW emissivity (LW co-albedo)
!  real(8), save :: emisli = 0.95d0  !! ice LW emissivity (LW co-albedo)
  logical, save :: olwnet = .false. !! true if longwave rad. is net flux
! namelist nmswem
  real(8), save :: albswo = 0.066d0 !! ocn. shortwave albedo (if constant)
!  real(8), save :: albswi = 0.5d0   !! bare ice shortwave albedo
!  real(8), save :: albsws(2) = &    !! snow-covered ice shortwave albedo
!    &                        (/ 0.85d0, 0.65d0 /)
!  real(8), save :: alswst(2) = &    !! snow-covered ice sw albedo, T dependency
!    &                        (/ -15.0d0, 0.0d0 /)
!  real(8), save :: alsdpt = 1.0d0   !! threshold thickness of snow-covered ice
  logical, save :: oswnet = .false. !! true if shortwave rad. is net flux
  logical, save :: oasold = .false. !! true if ocn. SW albedo depends on lat.
! namelist nmz0
  real(8), save :: z0fct = 0.1d0        !! heat z0/moumentum z0
  real(8), save :: z0min = 1.0d-6       !! minimum z0
! namelist nmocn
  real(8), save :: dzocn = 50.0d0       !! depth of ML ocean
  real(8), save :: dfocn = 1.0d10       !! ocean dG/dTs
  real(8), save :: alblo = 0.05d0       !! LW albedo (1-emis)
! namelist nmsfcl
  real(8), save :: za     = 10.0d0      !! height of atmospheric vars
  real(8), save :: cmmin  = 1.0d-5      !! min. bulk coef. of u
  real(8), save :: chmin  = 1.0d-5      !! min. bulk coef. of T
  real(8), save :: cemin  = 1.0d-5      !! min. bulk coef. of q
  real(8), save :: cmmax  = 1.0d0       !! max. bulk coef. of u
  real(8), save :: chmax  = 1.0d0       !! max. bulk coef. of T
  real(8), save :: cemax  = 1.0d0       !! max. bulk coef. of q
  real(8), save :: usminm = 0.5d0       !! min. wind vel. for v
  real(8), save :: usminh = 0.5d0       !! min. wind vel. for T
  real(8), save :: usmine = 0.5d0       !! min. wind vel. for q
  real(8), save :: usmaxm = 1.0d3       !! max wind vel. for v
  real(8), save :: usmaxh = 1.0d3       !! max. wind vel. for heat
  real(8), save :: usmaxe = 1.0d3       !! max wind vel. for q
! namelist nmsnit
  integer, save :: nitero = 2  !! n. of iteration for ocean surface
  integer, save :: niteri = 5  !! n. of iteration for ice surface
! namelist nmrwnd
  real(8), save :: alpha = 0.0d0  !! parameter for relative/absolute wind, 0.0: absolute -- 1.0: relative
! namelist nmtrdb
  logical, save :: trdbcp = .false.     !! trdbcp: True if PREC/SNOW/ROFF have traditional sign/unit (i.e., upward positive [m/s])
 
  namelist /nmsnow/  dfsnow, snwdmx, epssnw, snwmax, &
    &                albsnw, talsnw,  z0snw, &
    &                snwcrt, snwden  
  namelist /nmice/    dfice, albice,  z0ice, &
    &                siccrt, sicden
  namelist /nmswbr/    rvis,   rnir,    rir
  namelist /nmsfrc/  oasfrc, alsdpt, alspat
  namelist /nmsage/   osage, alssif, alssio, alfmax, snrfrs,  ftage, &
    &                tauage, adirt0, adirtc, adirts, adirtm, drsmax, &
    &                 oadst
  namelist /nmmpnd/   impnd, hminmp, rtdpmp, rtmxmp,  dpscl, &
    &                rmpcmn, rmpcmx, cmpfrz, tmpfrz, albmpd, almpdp, &
    &                frmpmn, vmpmin
  namelist /nmislt/      si
  namelist /nmlwem/  emislo, olwnet
  namelist /nmswem/  albswo, oswnet, oasold
  namelist /nmz0/     z0fct,  z0min
  namelist /nmocn/    dzocn,  dfocn,  alblo
  namelist /nmsair/      za
  namelist /nmsfcl/   cmmin,  chmin,  cemin, &
    &                 cmmax,  chmax,  cemax, &
    &                usminm, usmaxm, &
    &                usminh, usmaxh, &
    &                usmine, usmaxe
  namelist /nmsnit/  nitero, niteri
  namelist /nmrwnd/   alpha
  namelist /nmtrdb/  trdbcp

  private

  public :: sfcflx
#ifdef OPT_BODY
  public :: bdyflx
#endif

contains

subroutine sfcflx( &
  &                   qao,    qai,    qii,    qio,  swabs,    tsi, &
  &                   wev,    wsb, &
  &                  prec,   snow,   roff,   soff, &
  &                tauaix, tauaiy, tauaox, tauaoy, &
  &                    ft,   ptop,   ssfc, &
  &                  dfdu,   dfbc, &
  &                     t,      a,     hi,     ti,    hsn, &
  &                    as,    vmp,   frmp, &
  &                     u,      v )

  use qckot
  use ufile
  use utint
  use bshft
  use ucloc

  real(8), parameter :: factm = 1.0d+1, facth = 1.0d+3, factw = 1.0d+2, factmv = 1.0d-3

  real(8), intent(in)  ::       t(nxydim, nzdim, ntdim)
  real(8), intent(in)  ::       a(nxydim, 0:nic)
  real(8), intent(in)  ::      hi(nxydim, 0:nic),    hsn(nxydim, 0:nic)
  real(8), intent(in)  ::      ti(nxydim, 0:nic)
  real(8), intent(in)  ::      as(nxydim, 0:nic),    vmp(nxydim, 0:nic)
  real(8), intent(in)  ::    frmp(nxydim, 0:nic)
  real(8), intent(in)  ::       u(nxydim, nzdim),      v(nxydim, nzdim)

  real(8), intent(out) ::     wev(nxydim),    wsb(nxydim, nic)
  real(8), intent(out) ::    prec(nxydim),   snow(nxydim)
  real(8), intent(out) ::    roff(nxydim),   soff(nxydim)
  real(8), intent(out) ::  tauaix(nxydim), tauaiy(nxydim)
  real(8), intent(out) ::  tauaox(nxydim), tauaoy(nxydim)
  real(8), intent(out) ::     qai(nxydim, nic),    qio(nxydim, nic)
  real(8), intent(out) ::     qii(nxydim, nic)
  real(8), intent(out) ::     qao(nxydim),  swabs(nxydim)
  real(8), intent(out) ::     tsi(nxydim, 0:nic)
  real(8), intent(out) ::      ft(nxydim, ntdim)
  real(8), intent(out) ::    ptop(nxydim)
  real(8), intent(out) ::    ssfc(nxydim)
  real(8), intent(out) ::    dfdu(nxydim),   dfbc(nxydim)

  real(8) :: uo, vo
  real(8) ::    taux(nxydim)=0.d0,   tauy(nxydim)=0.d0
  real(8) ::    usfc(nxydim)=0.d0,   vsfc(nxydim)=0.d0
  real(8) ::     u10(nxydim)=0.d0,    v10(nxydim)=0.d0
  real(8) ::    tsfc(nxydim)=0.d0,   qsfc(nxydim)=0.d0,   pplr(nxydim)=0.d0
  real(8) ::    sflx(nxydim)=0.d0
  real(8) ::    swnt(nxydim)=0.d0,   dwlw(nxydim)=0.d0,   psfc(nxydim)=0.d0
  real(8) ::    grts(nxydim),   grtb(nxydim)
  real(8) ::   grice(nxydim),  grsnw(nxydim),  gricr(nxydim)
  real(8) ::   grasn(nxydim),  grvmp(nxydim), grfrmp(nxydim)
  real(8) ::    grz0(nxydim, ntyz0)
  real(8) ::  gfluxs(nxydim), tfluxs(nxydim)=0.d0, qfluxs(nxydim)=0.d0
  real(8) ::  wfluxs(nxydim, 2)
  real(8) ::  rflxlu(nxydim)=0.d0, sflxbl(nxydim)
  real(8) ::   dgfds(nxydim),  dtfdt(nxydim),  dtfds(nxydim)
  real(8) ::   dqfds(nxydim),   swdn(nxydim)
  real(8) ::      fm(nxydim)

  real(8) ::    swup(nxydim)=0.d0
  real(8) ::   swdnw(nxydim, 0:nic),           swupw(nxydim, 0:nic)
  real(8) ::   lwdnw(nxydim, 0:nic),           lwupw(nxydim, 0:nic)
  real(8) ::  swdnwg(nxydim), swupwg(nxydim)
  real(8) ::  swnetg(nxydim), lwnetg(nxydim)
  real(8) ::  lwdnwg(nxydim), lwupwg(nxydim)
  real(8) ::   senfx(nxydim, 0:nic),           latfx(nxydim, 0:nic)
  real(8) ::  senfxg(nxydim), latfxg(nxydim)
  real(8) ::    esub

  real(8), save ::     tmi,      factfw

  real(8) ::     dufdu ( nxydim )              !! CMV
! real(8) ::     dtfds ( nxydim )              !! CHV
  real(8) ::     dqfdq ( nxydim )              !! CEV
! real(8) ::     cmv(nxydim, 0:nic), chv(nxydim, 0:nic)
! real(8) ::     cev(nxydim, 0:nic)
      
 real(8) :: ralbsw(nxydim), albsw(nxydim, 0:nic)
! real(8) :: rwsb(nxydim, 0:nic)
! real(8) ::    rqai(nxydim, 0:nic),   rqio(nxydim, 0:nic)
! real(8) ::    rqii(nxydim, 0:nic)
! real(8) ::  ftatm(nxydim), swntwa(nxydim)
!  real(8) :: ralbsw(nxydim, 0:nic)
  real(8) ::   tisi(nxydim, 0:nic)
  real(8) ::   wsbg(nxydim), albswg(nxydim)
  
  real(8), save ::     dirdsn

  logical, save ::  ofirst = .true.

  integer ::     ij,      l,      n,     nn
  integer ::      i,      j
  integer ::  ifpar,  jfpar,  istat

  call clcstr('SFCFLX')
  
  if (ofirst) then
     call rewnml(ifpar, jfpar)
     write(jfpar, *) '*** sfcflx ***'
     call rewnml(ifpar, jfpar)
     read (ifpar, nmislt, iostat=istat)
     call cstnml(jfpar, 'sfcflx', 'nmislt', istat)
     write(jfpar, nmislt)
     call rewnml(ifpar, jfpar)
     read (ifpar, nmsnow, iostat=istat)
     call cstnml(jfpar, 'sfcflx', 'nmsnow', istat)
     write(jfpar, nmsnow)
     call rewnml(ifpar, jfpar)
     read (ifpar, nmice, iostat=istat)
     call cstnml(jfpar, 'sfcflx', 'nmice', istat)
     write(jfpar, nmice)
     call rewnml(ifpar, jfpar)
     read (ifpar, nmswbr, iostat=istat)
     call cstnml(jfpar, 'sfcflx', 'nmswbr', istat)
     write(jfpar, nmswbr)
     call rewnml(ifpar, jfpar)
     read (ifpar, nmsfrc, iostat=istat)
     call cstnml(jfpar, 'sfcflx', 'nmsfrc', istat)
     write(jfpar, nmsfrc)
     call rewnml(ifpar, jfpar)
     read (ifpar, nmsage, iostat=istat)
     call cstnml(jfpar, 'sfcflx', 'nmsage', istat)
     write(jfpar, nmsage)
     call rewnml(ifpar, jfpar)
     read (ifpar, nmmpnd, iostat=istat)
     call cstnml(jfpar, 'sfcflx', 'nmmpnd', istat)
     write(jfpar, nmmpnd)
     call rewnml(ifpar, jfpar)
     read (ifpar, nmlwem, iostat=istat)
     call cstnml(jfpar, 'sfcflx', 'nmlwem', istat)
     write(jfpar, nmlwem)
     call rewnml(ifpar, jfpar)
     read (ifpar, nmswem, iostat=istat)
     call cstnml(jfpar, 'sfcflx', 'nmswem', istat)
     write(jfpar, nmswem)
     call rewnml(ifpar, jfpar)
     read (ifpar, nmz0, iostat=istat)
     call cstnml(jfpar, 'sfcflx', 'nmz0', istat)
     write(jfpar, nmz0)
     call rewnml(ifpar, jfpar)
     read (ifpar, nmocn, iostat=istat)
     call cstnml(jfpar, 'sfcflx', 'nmocn', istat)
     write(jfpar, nmocn)
     call rewnml(ifpar, jfpar)
     read (ifpar, nmsair, iostat=istat)
     call cstnml(jfpar, 'sfcflx', 'nmsair', istat)
     write(jfpar, nmsair)
     call rewnml(ifpar, jfpar)
     read (ifpar, nmsfcl, iostat=istat)
     call cstnml(jfpar, 'sfcflx', 'nmsfcl', istat)
     write(jfpar, nmsfcl)
     call rewnml(ifpar, jfpar)
     read (ifpar, nmsnit, iostat=istat)
     call cstnml(jfpar, 'sfcflx', 'nmsnit', istat)
     write(jfpar, nmsnit)
     call rewnml(ifpar, jfpar)
     read (ifpar, nmrwnd, iostat=istat)
     call cstnml(jfpar, 'sfcflx', 'nmrwnd', istat)
     write(jfpar, nmrwnd)
     call rewnml(ifpar, jfpar)
     read (ifpar, nmtrdb, iostat=istat)
     call cstnml(jfpar, 'sfcflx', 'nmtrdb', istat)
     write(jfpar, nmtrdb)

     tmi = dtds * si         
     dirdsn = dfice / dfsnow
     do ij = 1, nxydim
        tsfc(ij) = 300.0d0
        psfc(ij) = 1.0d5
     end do
     if ( trdbcp ) then ! set conversion factor
        factfw = -1.0d0 * factw ! upward positive [m/s] --> downward positive [cm/s]
     else
        factfw = factmv * factw ! downward positive [kg/m^2/s] --> downward positive [cm/s]
     end if
     ofirst = .false.
  end if

  do ij = 1, nxydim
!     ftatm(ij) = 0.0d0
!     swntwa(ij) = 0.0d0
     wsbg(ij) = 0.0d0
     albswg(ij) = 0.0d0
  end do
  do l = 0, nic
     do ij = 1, nxydim
!        rqai(ij, l) = 0.0d0
!        rqio(ij, l) = 0.0d0
!        rqii(ij, l) = 0.0d0
        tisi(ij, l) = 0.0d0
     end do
  end do

  call clcstr('TMINTP')
  
  do l = 3, ntdim
     n = 12 + (l - 3) * 2
     nn = n + 1
     call tmintp(  tsfc,      n)
     call tmintp(  qsfc,     nn)
     do ij = ijstr, ijend
        ft(ij, l) = qsfc(ij) * (tsfc(ij) - t(ij, kstr, l)) * &
          &                    dz(ij, kstr) * amskt(ij, kstr)
     end do
  end do
  call tmintp_direct(   u10,      1)
  call tmintp_direct(   v10,      2)
  call tmintp_direct(  tsfc,      3)
  call tmintp_direct(  qsfc,      4)
  call tmintp_direct(  pplr,      5)
  call tmintp_direct(  sflx,      6)
  call tmintp_direct(  swnt,      7)
  call tmintp_direct(  dwlw,      8)
  call tmintp_direct(  psfc,      9)
  call tmintp_direct(  roff,     10)

  ssfc(:)=0.D0
#ifdef OPT_SRST
  call tmintp(  ssfc,     11)
#endif

  call clcend('TMINTP')
  
!! 2021.05.31: Now dfdu and dfbc are dummy fluxes in OGCM.
  dfdu(:) = 0.0d0
  dfbc(:) = 0.0d0

#ifdef OPT_TRIPOLE
  call shift2( &
    &            psfc,   roff,         &
    &           nxdim,  nydim,      1, &
    &            1.d0,      0,      0 )
#else
  call shift2( &
    &            psfc,   roff,         &
    &           nxdim,  nydim,      1 )
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do ij = ijstr, ijend
     tauaix(ij) = 0.0d0
     tauaiy(ij) = 0.0d0
     prec(ij) = pplr(ij) * factfw ! prec: downward is positive
     snow(ij) = sflx(ij) * factfw ! snow: downward is positive
     soff(ij) = 0.0d0
!     roff(ij) = 0.0d0
     roff(ij) = roff(ij) * factfw ! roff: downward is positive
  end do

  do ij = 1, nxydim
     usfc(ij) = cos( rangt(ij) ) * u10(ij) - sin( rangt(ij) ) * v10(ij)
     vsfc(ij) = sin( rangt(ij) ) * u10(ij) + cos( rangt(ij) ) * v10(ij)
  end do

  do ij = ijstr, ijend
     uo = 0.01d0 * 0.25d0 &
          & * ( u(ij,    kstr) + u(ij+lw,  kstr) &
          & +   u(ij+ls, kstr) + u(ij+lsw, kstr) )
     vo = 0.01d0 * 0.25d0 &
          & * ( v(ij,    kstr) + v(ij+lw,  kstr) &
          & +   v(ij+ls, kstr) + v(ij+lsw, kstr) )
     usfc(ij) = usfc(ij) - uo * a(ij, 0) * alpha
     vsfc(ij) = vsfc(ij) - vo * a(ij, 0) * alpha
  end do

  esub = el + emelt
  do l = 0, nic
     do ij = 1, nxydim
        fm(ij) = a(ij, l)
     end do

     if (l > 0) then
        do ij = ijstr, ijend
           grts (ij) = tsi(ij, l) + kelvin
           grtb (ij) = ti(ij, l) + kelvin
           grice(ij) = hi(ij, l) * 5.0d-3
           grsnw(ij) = hsn(ij, l) * 1.0d-2
           gricr(ij) = 1.0d0
           grasn(ij) = as(ij, l)
           grvmp(ij) = vmp(ij, l) * 1.0d-2
           grfrmp(ij) = frmp(ij, l)
        end do
     else
        do ij = ijstr, ijend
           grts (ij) = t(ij, kstr, 1) + kelvin
           grtb (ij) = t(ij, kstr, 1) + kelvin
           grice(ij) = 0.0d0
           grsnw(ij) = 0.0d0
           gricr(ij) = 0.0d0
           grasn(ij) = 0.0d0
           grvmp(ij) = 0.0d0
           grfrmp(ij) = 0.0d0
        end do
     end if

     call ocnbcs_core( &
       &               gfluxs,  dgfds, &
       &                 grts,   grtb,  grice,  grsnw,  gricr )
     call sfcflx_core( &
       &               tfluxs, qfluxs,   taux,   tauy, &
       &                dtfdt,  dtfds,  dqfds, &
       &                dufdu,  dqfdq, &
       &                 usfc,   vsfc,   tsfc,   qsfc,   psfc, &
       &                 grts,      l )
!     do ij=1, nxydim
!        cmv(ij, l) = dufdu(ij)
!        chv(ij, l) = dtfds(ij)
!        cev(ij, l) = dqfdq(ij)
!     end do

     if (l > 0) then
        do ij = ijstr, ijend
           grice(ij) = hi(ij, l) * 1.0d-2
        end do
     end if
     call ocnslv_core( &
       &                 grts, gfluxs, tfluxs, qfluxs, &
       &               wfluxs, rflxlu, sflxbl,   swdn, ralbsw, &
       &                 swup, &
       &                dtfds,  dqfds,  dgfds, &
       &                 swnt,   dwlw, &
       &                gricr,  grsnw,    tmi, &
       &                grasn,  grvmp, grfrmp )     

#ifdef OPT_TRIPOLE
     call shift2( &
       &            taux,   tauy, &
       &           nxdim,  nydim,      1, &
       &          -1.0d0,      0,      0 )
     call shift1( &
       &              fm, &
       &           nxdim,  nydim,      1, &
       &           1.0d0,      0,      0 )
#else
     call shift3( &
       &            taux,   tauy,    fm, &
       &           nxdim,  nydim,     1 )
#endif

     if (l > 0) then
        do ij = ijstr, ijend
           tauaix(ij) = tauaix(ij) &
             &        + ( taux(ij) + taux(ij+le) &
             &          + taux(ij+ln) + taux(ij+lne) ) &
             &                    * 0.25d0 * factm * fm(ij)
           tauaiy(ij) = tauaiy(ij) &
             &        + ( tauy(ij) + tauy(ij+le) &
             &          + tauy(ij+ln) + tauy(ij+lne) ) &
             &                    * 0.25d0 * factm * fm(ij)
           qai(ij, l) = (gfluxs(ij) - sflxbl(ij)) * facth
           qii(ij, l) = gfluxs(ij) * facth
           qio(ij, l) = (t(ij, kstr, 1) - ti(ij, l)) * &
             &            cdi / hi(ij, l) * 2.0d0
           soff(ij) = 0.0d0
           tsi(ij, l) = min(grts(ij), tmelt+tmi) - kelvin
           tisi(ij, l) = tsi(ij, l) + &
             &           (ti(ij, l) - tsi(ij, l)) * &
             &           dirdsn * hsn(ij, l) / &
             &           (0.5d0 * hi(ij, l) + dirdsn * hsn(ij, l))
           wsb(ij, l) = qfluxs(ij) * a(ij, l) / dwatr * factw
           wsbg(ij) = wsbg(ij) + &
             &        qfluxs(ij) * a(ij, l) / dwatr * factw
!           rqai(ij, l) = qai(ij, l)
!           rqii(ij, l) = qii(ij, l)
!           rqio(ij, l) = qio(ij, l)
!           rwsb(ij, l) = wsb(ij, l)
           albsw(ij, l) = ralbsw(ij)
           albswg(ij) = albswg(ij) + ralbsw(ij) * a(ij, l)
!           ftatm(ij) = ftatm(ij) + qai(ij, l) * fm(ij)
!           swntwa(ij) = swntwa(ij) + swdn(ij) * fm(ij)
        end do
        do ij = 1, nxydim
           latfx(ij, l) = qfluxs(ij) * esub
        end do
     else
        do ij = ijstr, ijend
           tauaox(ij) = (  taux(ij) + taux(ij+le) &
             &           + taux(ij+ln) + taux(ij+lne)) * &
             &          0.25d0 * factm
           tauaoy(ij) = (  tauy(ij) + tauy(ij+le) &
             &           + tauy(ij+ln) + tauy(ij+lne)) * &
             &          0.25d0 * factm
           qao(ij) = (gfluxs(ij) + swdn(ij)) * facth
           swabs(ij) = swdn(ij) * facth * a(ij, 0)
           wev(ij) = qfluxs(ij) * a(ij, 0) / dwatr * factw
           albsw(ij, l) = ralbsw(ij)
!           ftatm(ij) = ftatm(ij) + gfluxs(ij) * facth * a(ij, 0)
!           swntwa(ij) = swntwa(ij) + swdn(ij) * a(ij, 0)
        end do
        do ij = 1, nxydim
           latfx(ij, l) = qfluxs(ij) * el
        end do
     end if
     do ij = 1, nxydim
        swdnw(ij, l) = swnt(ij)
        swupw(ij, l) = swup(ij)
        lwdnw(ij, l) = dwlw(ij)
        lwupw(ij, l) = rflxlu(ij)
        senfx(ij, l) = tfluxs(ij)
     end do
  end do

  do ij = 1, nxydim
     if( a(ij,0) /= 1.d0 ) then
        tauaix(ij) = tauaix(ij) / (1.0d0 - a(ij, 0) )
        tauaiy(ij) = tauaiy(ij) / (1.0d0 - a(ij, 0) )
        albswg(ij) = albswg(ij) / (1.0d0 - a(ij, 0) )
     endif
  enddo

  swdnwg(:) = 0.d0
  swupwg(:) = 0.d0
  lwdnwg(:) = 0.d0
  lwupwg(:) = 0.d0
  senfxg(:) = 0.d0
  latfxg(:) = 0.d0
  do l = 0, nic
     swdnwg(:) = swdnwg(:) + swdnw(:, l) * facth * a(:, l) * amskt(:, kstr)
     swupwg(:) = swupwg(:) + swupw(:, l) * facth * a(:, l) * amskt(:, kstr)
     lwdnwg(:) = lwdnwg(:) + lwdnw(:, l) * facth * a(:, l) * amskt(:, kstr)
     lwupwg(:) = lwupwg(:) + lwupw(:, l) * facth * a(:, l) * amskt(:, kstr)
     senfxg(:) = senfxg(:) + senfx(:, l) * facth * a(:, l) * amskt(:, kstr)
     latfxg(:) = latfxg(:) + latfx(:, l) * facth * a(:, l) * amskt(:, kstr)
  end do
  swnetg(:) = swupwg(:) - swdnwg(:)
  lwnetg(:) = lwupwg(:) - lwdnwg(:)

  call chekin( swdnwg, 'SWDNWG', &
    &      'downward shortwave (ocn/ice top)', &
    &              'erg/cm^2/s', &
    &              nx,       ny,        1,   nxydim, 'OCSFCT')
  call chekin( swupwg, 'SWUPWG', &
    &      'upward shortwave (ocn/ice top)', &
    &              'erg/cm^2/s', &
    &              nx,       ny,        1,   nxydim, 'OCSFCT')
  call chekin( lwdnwg, 'LWDNWG', &
    &      'downard longwave (ocn/ice top)', &
    &              'erg/cm^2/s', &
    &              nx,       ny,        1,   nxydim, 'OCSFCT')
  call chekin( lwupwg, 'LWUPWG', &
    &      'upward longwave (ocn/ice top)', &
    &              'erg/cm^2/s', &
    &              nx,       ny,        1,   nxydim, 'OCSFCT')
  call chekin( lwnetg, 'LWNETG', &
    &      'net LW (ocn/ice top; up.:+)', &
    &              'erg/cm^2/s', &
    &              nx,       ny,        1,   nxydim, 'OCSFCT')
  call chekin( swnetg, 'SWNETG', &
    &      'net SW (ocn/ice top; up.:+)', &
    &              'erg/cm^2/s', &
    &              nx,       ny,        1,   nxydim, 'OCSFCT')
  call chekin( senfxg, 'SENFXG', &
    &      'sens. hflx (ocn/ice top; up.:+)', &
    &              'erg/cm^2/s', &
    &              nx,       ny,        1,   nxydim, 'OCSFCT')
  call chekin( latfxg, 'LATFXG', &
    &      'latent hflx (ocn/ice top; up.:+)', &
    &              'erg/cm^2/s', &
    &              nx,       ny,        1,   nxydim, 'OCSFCT')

#ifdef OPT_TRIPOLE
  call shift2( tauaox,   tauaoy, &
    &           nxdim,    nydim,    1, &
    &          -1.0d0,     -1,     -1 )

  call shift2( tauaix,   tauaiy, &
    &           nxdim,    nydim,    1, &
    &          -1.0d0,     -1,     -1 )
#else
  call shift2( &
    &          tauaox,   tauaoy, &
    &           nxdim,    nydim,     1 )
  call shift2( &
    &          tauaix,   tauaiy, &
    &           nxdim,    nydim,     1 )
#endif

  call chekin(   roff,  'ROFF', &
    &       'river runoff', 'cm/s', &
    &              nx,      ny,      1, nxydim, 'OCSFCT')
  call chekin(   usfc,  'USFC', &
    &       'zonal wind speed', 'm/s', &
    &              nx,      ny,      1, nxydim, 'OCSFCT')
  call chekin(   vsfc,  'VSFC', &
    &       'meridional wind speed', 'm/s', &
    &              nx,      ny,      1, nxydim, 'OCSFCT')

  call chekin(   wsbg,  'WSBG', &
    &       'sublimation from sea-ice surface', 'cm/s', &
    &              nx,      ny,      1, nxydim, 'OCSFCT')
  call chekin(    wev,   'WEV', &
    &           'evaporation from sea surface', 'cm/s', &
    &              nx,      ny,      1, nxydim, 'OCSFCT')
  call chekin(  albsw(1,0), 'ALBSWO', &
    &           'ocean surface albedo', 'ND', &
    &              nx,      ny,      1, nxydim, 'OCSFCT')
  call chekin(  albsw, 'ALBSWI', &
    &           'sea-ice surface albedo', 'ND', &
    &              nx,      ny,    nic, nxyidm, 'OCICET')
  call chekin( albswg, 'ALSWIG', &
    &           'sea-ice surface albedo', 'ND', &
    &              nx,      ny,      1, nxydim, 'OCSFCT')
!  call chekin( tauaox,'tauaox', &
!    &              nx,      ny,      1, nxydim, 'sfc')
!  call chekin( tauaoy,'tauaoy', &
!    &              nx,      ny,      1, nxydim, 'sfc')
!  call chekin( tauaix,'tauaix', &
!    &              nx,      ny,      1, nxydim, 'sfc')
!  call chekin( tauaiy,'tauaiy', &
!    &              nx,      ny,      1, nxydim, 'sfc')
!  call chekin(    cmv,  'cmvo', &
!    &              nx,      ny,      1, nxydim, 'sfc')
!  call chekin(    chv,  'chvo', &
!    &              nx,      ny,      1, nxydim, 'sfc')
!  call chekin(    cev,  'cevo', &
!    &              nx,      ny,      1, nxydim, 'sfc')
!  call chekin(    cmv,  'cmvi', &
!    &              nx,      ny,    nic, nxyidm, 'ice')
!  call chekin(    chv,  'chvi', &
!    &              nx,      ny,    nic, nxyidm, 'ice')
!  call chekin(    cev,  'cevi', &
!    &              nx,      ny,    nic, nxyidm, 'ice')
!  call chekin(  albsw,'albswo', &
!    &              nx,      ny,      1, nxydim, 'sfc')
!  call chekin(  albsw,'albswi', &
!    &              nx,      ny,    nic, nxyidm, 'ice')
!  call chekin(  ftatm, 'ftatm', &
!    &              nx,      ny,      1, nxydim, 'sfc')
!  call chekin( swntwa,  'swnt', &
!    &              nx,      ny,      1, nxydim, 'sfc')
!  call chekin(   rqai,   'qai', &
!    &              nx,      ny,    nic, nxyidm, 'ice')
!  call chekin(   rqio,   'qio', &
!    &              nx,      ny,    nic, nxyidm, 'ice')
!  call chekin(   rqii,   'qii', &
!    &              nx,      ny,    nic, nxyidm, 'ice')

! output section for CMIP5
! tisi: temperature at interface between sea ice and snow
!       unit [C]
  call chekin(tisi,   'TISI',                                       &
    &            'temperature at ice-snow interface', 'degC', &
    &            nx,     ny,    nic, nxyidm, 'OCICET')

  do ij = 1, nxydim
!     ptop(ij) = 0.d0
     ptop(ij) = psfc(ij) * factm
  end do

  call clcend('SFCFLX')
  
  return
end subroutine sfcflx

! *********************************************************************

subroutine ocnslv_core ( &
  &              gdts  , gfluxs, tfluxs, qfluxs, &
  &              wfluxs, rflxlu, sflxbl, rflxsd, ralbsw, &
  &              rflxsu, &
  &              dtfds , dqfds , dgfds , &
  &              rflxs , rflxld, &
  &              gricr , grsnw , tmi   , &
  &              grasn , grvmp , grfrmp )

  use ufile
   
  real(8), intent(inout) ::  gdts  ( nxydim )          !! skin temperature
  real(8), intent(inout) ::  gfluxs( nxydim )          !! soil heat flux
  real(8), intent(inout) ::  tfluxs( nxydim )          !! flux of t
  real(8), intent(inout) ::  qfluxs( nxydim )          !! flux of q

  real(8), intent(out)   ::  wfluxs( nxydim, 2 )       !! soil water flux
  real(8), intent(out)   ::  rflxlu( nxydim )          !! upward long wave
  real(8), intent(out)   ::  sflxbl( nxydim )          !! flux balance
  real(8), intent(out)   ::  rflxsd( nxydim )          !! downward SW
  real(8), intent(out)   ::  ralbsw( nxydim )          !! SW albedo
  
  real(8), intent(out)   ::  rflxsu( nxydim )          !! upward SW

  real(8), intent(in)    ::  dtfds ( nxydim )          !! dH/dTg
  real(8), intent(in)    ::  dqfds ( nxydim )          !! dE/dTg
  real(8), intent(in)    ::  dgfds ( nxydim )          !! dG/dTg
  real(8), intent(in)    ::  rflxs ( nxydim )          !! SW rad. (net/down.)
  real(8), intent(in)    ::  rflxld( nxydim )          !! down. LW rad.
  real(8), intent(in)    ::  gricr ( nxydim )          !! snow/ice ratio
  real(8), intent(in)    ::  grsnw ( nxydim )          !! snow thickness
  real(8), intent(in)    ::  tmi           !! sea ice melting temperature (C)
  real(8), intent(in)    ::  grasn ( nxydim )      !! snow age
  real(8), intent(in)    ::  grvmp ( nxydim )      !! melt pond volume
  real(8), intent(in)    ::  grfrmp( nxydim )      !! melt pond fraction

  real(8) ::  grsnr ( nxydim )      !! snow cover fraction
  real(8) ::  hsnow ( nxydim )      !! snow depth [m]
  real(8) ::  hmp ( nxydim )        !! melt pond depth [m]
  real(8) ::  rp ( nxydim )         !! retaind melt water ratio [ND]

  real(8) ::     esub, stg, drfds
  real(8) ::     sflux, gsflux, dgsfds
  real(8) ::     gfluxf
  real(8) ::     sflxbi, dsbdsi, dti, evapi, gfluxi
  real(8) ::     ff, fi, dtx
  real(8) ::     emis
  real(8) ::     x, albx, albsw, icealb
  real(8) ::     mpdalb, snwalb, smpalb, brialb
  real(8) ::     fbarei, fmpnd, fsnow, fsnwmp
  real(8) ::     hsneff, grsref, hslash, slsalb
  real(8) ::     omega, sinij, cort
  integer ::    ij, l
  integer ::  ifpar,  jfpar

  real(8), save :: aswo2d(nxydim) !! ocn. shortwave albedo distribution
  real(8), save :: tsdpt
  real(8), save :: flwnet, fswalb, fusemp, aicet1, daicet
  real(8), save :: emisli, alcice, alcsnw(2), alcsif, alcsio, alcmpd
  real(8), save :: dalmdp, falmdp
  real(8), save :: rsrro, rorros

  logical, save :: ofirst = .true.

  if ( ofirst ) then
     call rewnml(ifpar, jfpar)
     write (jfpar, *) ' @@@ OCNSLV: OCEAN HEAT BALANCE 98/06/30'
     ofirst = .false.

     if (olwnet) then
        flwnet = 0.d0
     else
        flwnet = 1.d0
     end if
     if (oswnet) then
        fswalb = 0.d0
     else
        fswalb = 1.d0
     end if
     if (oasfrc) then
        tsdpt = alspat * 1.0d-2
     else
        tsdpt = alsdpt * 1.0d-2
     end if
     if (impnd > 0) then
       fusemp = 1.0d0
     else
       fusemp = 0.0d0
     end if
     aicet1 = talsnw(1)
     daicet = talsnw(2) - talsnw(1)
     dalmdp = (almpdp(2) - almpdp(1)) * 1.0d-2
     falmdp = almpdp(1) * 1.0d-2
     alcice = rvis * albice(1) + rnir * albice(2) + rir * albice(3)
     do l = 1, 2
       alcsnw(l) = rvis * albsnw(l,1) + rnir * albsnw(l,2) &
         &       + rir * albsnw(l,3)
     end do
     alcsif = rvis * alssif(1) + rnir * alssif(2) + rir * alssif(3)
     alcsio = rvis * alssio(1) + rnir * alssio(2) + rir * alssio(3)
     alcmpd = rvis * albmpd(1) + rnir * albmpd(2) + rir * albmpd(3)
     emisli = 1.0d0 - albice(3)
     rsrro = rhos/rhoo
     rorros = rhoo/(rhoo-rhos)

!    '08.09.16: Large and Yeager (2008) latitude-dependant albedo
     if (oasold) then
        write (jfpar, *) &
          &   '== L-Y2008 latitude-dependant ocn. shortwave albedo =='
        omega = 2.0d0 * ( atan( 1.0d0 )*4.0d0 ) / 86400.d0
        do ij = ijstr, ijend
           cort=(cor(ij)+cor(ij+lw)+cor(ij+lsw)+cor(ij+ls))*0.25d0
           sinij = min(max(cort/(2.0d0 * omega), -1.0d0), 1.0d0)
           aswo2d(ij) = 0.069d0 - 0.011 * & 
             &          (1.0d0 - 2.0d0 * sinij * sinij)
        end do
     else
        write (jfpar, *) &
          &   '== constant ocn. shortwave albedo,', albswo, ' =='
        do ij = ijstr, ijend
           aswo2d(ij) = albswo
        end do
     end if
  endif

  do ij = ijstr, ijend
     if (oasfrc) then
        grsnr(ij) = grsnw(ij) / (tsdpt + grsnw(ij))
     else
        if (grsnw(ij) > tsdpt) then
           grsnr(ij) = 1.d0
        else
           grsnr(ij) = 0.d0
        end if
     end if
     if (grsnr(ij) > 0.0d0) then
        hsnow(ij) = grsnw(ij) / grsnr(ij)
     else
        hsnow(ij) = 0.0d0
     end if
     if (grfrmp(ij) > 0.0d0) then
        hmp(ij) = grvmp(ij) / grfrmp(ij)
     else
        hmp(ij) = 0.0d0
     end if
     rp(ij) = hmp(ij) + hsnow(ij) * rsrro
     if (rp(ij) > 0.0d0) then
        rp(ij) = hmp(ij) / rp(ij)
     else
        rp(ij) = 0.0d0
     end if
  end do

  esub = el + emelt
  brialb = alcice
  do ij = ijstr, ijend
     emis = emislo * (1.0d0 - gricr(ij)) + emisli * gricr(ij)
     albx = fusemp * min( max( &
       &            (hmp(ij) - falmdp) / dalmdp, 0.0d0), 1.0d0)
     mpdalb = alcice * (1.0d0 - albx) + alcmpd * albx
     if (osage) then
        snwalb = alcsif + grasn(ij) * ( alcsio - alcsif )
     else
        x = min( max( (gdts(ij) - aicet1) / daicet, 0.0d0), 1.0d0)
        snwalb = alcsnw(1) * (1.0d0 - x) + alcsnw(2) * x
     end if
     if (impnd == 2) then  !! Hunke MP param.
        fsnow  = (1.0d0 - grfrmp(ij)) * grsnr(ij)
        fmpnd  = grfrmp(ij) * (1.0d0 - grsnr(ij))
        fsnwmp = grfrmp(ij) * grsnr(ij)
        fbarei = 1.0d0 - (fsnow + fmpnd + fsnwmp)  !! = (1-grfrmp) * (1-grsnr)
        if (rp(ij) <= 0.15d0) then  !! all the MP water retained in snow
           smpalb = snwalb
        else
           if (grsnr(ij) > 0.0d0) then
              hsneff = hsnow(ij) - hmp(ij) * rorros
              grsref = max( 0.0d0, min( 1.0d0, &
                &      ( hsneff / (tsdpt + hsneff) ) / grsnr(ij) ) )
           else
              grsref = 0.0d0
           end if
           hslash = hmp(ij) + rsrro * hsnow(ij)
           albx = fmpnd * min( max( &
           &            (hslash - falmdp) / dalmdp, 0.0d0), 1.0d0)
           slsalb = alcice * (1.0d0 - albx) + alcmpd * albx
           smpalb = grsref * snwalb + (1.0d0 - grsref) * slsalb
        end if
     else
        fsnow  = grsnr(ij)
        fmpnd  = min(grfrmp(ij), 1.0d0-grsnr(ij))
        fbarei = 1.0d0 - (fsnow + fmpnd)
        fsnwmp = 0.0d0
        smpalb = 0.0d0
     end if
     icealb = snwalb * fsnow  &
       &    + mpdalb * fmpnd  &
       &    + smpalb * fsnwmp &
       &    + brialb * fbarei
     albsw = (aswo2d(ij) * (1.0d0 - gricr(ij)) + icealb * gricr(ij)) &
       &     * fswalb
     ralbsw(ij) = albsw
     stg            = emis*stb*gdts( ij )**4
     rflxlu( ij )   = (stg + ( 1.d0-emis )*rflxld( ij )) * flwnet
     drfds          = 4.d0*stg/gdts( ij ) * flwnet

     rflxsd( ij )   = rflxs( ij ) * (1.0d0 - albsw)
     rflxsu( ij )   = rflxs( ij ) * albsw
     sflux          = tfluxs( ij ) &
       &            + rflxlu( ij ) - rflxld( ij ) &
       &            - rflxsd( ij )
     gsflux         = gfluxs( ij ) - sflux
     dgsfds         = dgfds ( ij ) &
       &            + dtfds ( ij ) + drfds

!    < ice free ocean surface >
     gfluxf         = sflux + qfluxs( ij )*el
         
!    < ice/snow surface >
     sflxbi         = gsflux - qfluxs( ij )*esub
     dsbdsi         = dgsfds + dqfds ( ij )*esub
     dti            = sflxbi/dsbdsi
     dti            = min( dti, tmelt+tmi-gdts( ij ) )
     evapi          = qfluxs( ij ) + dqfds( ij )*dti
     gfluxi         = gfluxs( ij ) - dgfds( ij )*dti
     sflxbi         = sflxbi       - dsbdsi     *dti

!    < Ts change >
     ff             = 1.0d0 - gricr( ij )
     fi             = gricr( ij )
     dtx            = fi * dti
     gdts  ( ij )   = gdts( ij ) + dtx

!    < fluxs >
     tfluxs( ij )   = tfluxs( ij ) + dtfds( ij )*dtx
     rflxlu( ij )   = rflxlu( ij ) + drfds      *dtx
     qfluxs( ij )   = ff * qfluxs( ij ) + fi * evapi
     gfluxs( ij )   = ff * gfluxf       + fi * gfluxi 
     wfluxs( ij,1 ) = ff * qfluxs( ij )
     wfluxs( ij,2 ) = fi * evapi
     sflxbl( ij )   = fi * sflxbi

  end do

  return
end subroutine ocnslv_core

! *********************************************************************

subroutine ocnbcs_core ( &
  &                      fogflx, dgfds , &
  &                      grts  , grtb  , grice , grsnw , gricr )

  use ufile

  real(8), intent(out) :: fogflx( nxydim )      !! heat flux
  real(8), intent(out) :: dgfds ( nxydim )      !! dG/dTs
  
  real(8), intent(in)  :: grts  ( nxydim )      !! skin temperature
  real(8), intent(in)  :: grtb  ( nxydim )      !! ice base temp.
  real(8), intent(in)  :: grice ( nxydim )      !! sea ice
  real(8), intent(in)  :: grsnw ( nxydim )      !! snow smount
  real(8), intent(in)  :: gricr ( nxydim )      !! ice fraction

  real(8) :: grsnrf ( nxydim )      !! snow cover frac. for flux calc.

  integer ::     ij,     l,      m
  integer ::  ifpar,  jfpar
  integer ::    ifg
  real(8) :: talsnx, albsnx,   alb0,   dalb,  tfact
  real(8) ::    z00,    dz0
  real(8) ::   dfgt,   dfgx

  logical, save :: ofirst = .true.

  if ( ofirst ) then
     call rewnml(ifpar, jfpar)
     write (jfpar, *) ' @@@ OCNBCS: OCEAN SURFACE BC 98/07/29'
     ofirst = .false.
  endif

  do ij = ijstr, ijend
     if (grsnw(ij) .gt. 0.0d0) then
         grsnrf(ij) = 1.d0
     else
         grsnrf(ij) = 0.d0
     end if
  end do

  do ij = ijstr, ijend
     if (gricr(ij) > 0.0d0) then
        dfgt = dfice / grice(ij)
        dfgx = dfice * dfsnow &
          &    / (dfice * grsnw(ij) + dfsnow * grice(ij))
        dfgt = dfgt * (1.0d0 - grsnrf(ij)) + dfgx * grsnrf(ij)
        fogflx(ij) = dfgt * (grtb(ij) - grts(ij))
        dgfds (ij) = dfgt
     else
        fogflx(ij) = 0.0d0
        dgfds (ij) = dfocn
     endif
  end do

  return
end subroutine ocnbcs_core

! *********************************************************************

subroutine sfcflx_core ( &
  &                      tfluxs, qfluxs,   taux,    tauy, &
  &                       dtfdt,  dtfds,  dqfds, &
  &                       dufdu,  dqfdq, &
  &                        usfc,   vsfc,   gdta,    gdqa,   gdps, &
  &                        gdts,      l )

  use ufile
  use zqsat
  real(8), intent(out) :: tfluxs( nxydim )      !! flux of T
  real(8), intent(out) :: qfluxs( nxydim )      !! flux of Q
  real(8), intent(out) ::   taux( nxydim )      !! flux of U
  real(8), intent(out) ::   tauy( nxydim )      !! flux of V
  real(8), intent(out) :: dufdu ( nxydim )      !! -d(tau)/du
  real(8), intent(out) :: dtfdt ( nxydim )      !! -dH/dTa
  real(8), intent(out) :: dqfdq ( nxydim )      !! -dLE/dqa
  real(8), intent(out) :: dtfds ( nxydim )      !! dH/dTs
  real(8), intent(out) :: dqfds ( nxydim )      !! dLE/dTs
  
  real(8), intent(in)  :: usfc  ( nxydim )      !! westerly u
  real(8), intent(in)  :: vsfc  ( nxydim )      !! southerly v
  real(8), intent(in)  :: gdta  ( nxydim )      !! temperature T
  real(8), intent(in)  :: gdqa  ( nxydim )      !! humidity q
  real(8), intent(in)  :: gdps  ( nxydim )
  real(8), intent(in)  :: gdts  ( nxydim )      !! surface skin temp
  integer, intent(in)  :: l                     !! 1 = ice, 2 = ocean

  logical, save :: ofirst = .true.
  real(8), save :: akappa

  integer :: ij
  integer :: ifpar,  jfpar
  real(8) :: exi, rho, qsat, dqsat, cqsat, uabs

  if ( ofirst ) then
     call rewnml(ifpar, jfpar)
     write(jfpar, *) ' @@@ PSFCM: SURFACE FLUX 98/06/19'
     ofirst = .false.
     akappa = rair/cp
  endif

! '08.06.05
! Bulk coeff. are now diagnosed without distinction of ocean or ice 
  call blkcof_core ( &
    &                 dufdu,  dtfds,  dqfdq, &  !! CMV, CHV, CEV
    &                  usfc,   vsfc,   gdta,   gdqa,   gdts,   gdps,    l )
  if (l == 0) then   !! oceanic surface
     cqsat = 0.98d0
  else   !! sea-ice surface
     cqsat = 1.0d0
  endif

  rho = 1.22d0
  do ij = ijstr, ijend
!     rho = gdps(ij) / (rair + (rvap - rair) * gdqa(ij)) / gdta(ij)
     exi = (gdps(ij) / (gdps(ij) - rho * grav * za))**akappa

     qsat         = fqsat ( gdts( ij ), gdps(ij) )
     dqsat        = fdqsat( gdts( ij ), qsat )
     qsat = cqsat * qsat
     
     dufdu ( ij ) =    rho*dufdu( ij )
     dtfds ( ij ) = cp*rho*dtfds( ij )
     dqfdq ( ij ) =    rho*dqfdq( ij )

     dtfdt ( ij ) =   dtfds( ij )*exi
     dqfds ( ij ) =   dqfdq( ij )*dqsat

     tfluxs( ij ) =   dtfds( ij )*( gdts( ij )-gdta( ij )*exi )
!    '08.08.29  Y.Komuro: qfluxs must be >= 0
     qfluxs( ij ) =   dqfdq( ij )* &
     &                  max( ( qsat - gdqa( ij ) ), 0.0d0 )
!     qfluxs( ij ) =   dqfdq( ij )*( qsat      -gdqa( ij )     )

     taux(ij) = dufdu(ij) * usfc(ij)
     tauy(ij) = dufdu(ij) * vsfc(ij)
  end do

  return
end subroutine sfcflx_core

! *********************************************************************

subroutine blkcof_core ( &
  &                       cmv,    chv,    cev, &
  &                      usfc,   vsfc,   gdta,   gdqa,   gdts,   gdps,     l )

  use ufile
  use zqsat

  real(8), intent(out) :: cmv   ( nxydim )   !! bulk transfer coeff.: u
  real(8), intent(out) :: chv   ( nxydim )   !! bulk transfer coeff.: T
  real(8), intent(out) :: cev   ( nxydim )   !! bulk transfer coeff.: q

  real(8), intent(in)  :: usfc  ( nxydim )   !! 10 m westerly u
  real(8), intent(in)  :: vsfc  ( nxydim )   !! 10 m southerly v
  real(8), intent(in)  :: gdta  ( nxydim )   !! 10 m air temperature
  real(8), intent(in)  :: gdts  ( nxydim )   !! surface skin temp
  real(8), intent(in)  :: gdqa  ( nxydim )   !! 10 m specific hum.
  real(8), intent(in)  :: gdps  ( nxydim )   !! sea level pressure
  integer, intent(in)  :: l                  !! 1 for ice, 2 for ocean

  real(8), save ::       pi,   exi, akappa
  logical, save :: ofirst = .true.

  real(8) ::  theta, thetav,    rho,   uabs,   qsat,  cmvsqr
  real(8) ::  ustar,  tstar,  qstar,   zeta
  real(8) ::      x,  psi_m,  psi_h,   un10, fncmvr,      xx,     x2
  real(8) ::  astbl
! bug fix 080422
  real(8) :: phrcdk
! bug fix 080908
  real(8) ::  cdn10, cdn10r,  cen10,  chn10
  integer ::     ij,      n
  integer ::  ifpar,  jfpar
  integer ::  niter

  if ( ofirst ) then
     call rewnml(ifpar, jfpar)
     write(jfpar, *) &
       &      ' @@@ PSFCL: SURFACE BULK COEF. Large AND Yeager' 
     ofirst = .false.
     pi = 4.0d0 * atan(1.0d0)
     akappa = rair / cp
  endif

  rho = 1.22d0
  do ij = ijstr, ijend
!     rho = gdps(ij) / (rair + (rvap - rair) * gdqa(ij)) / gdta(ij)
!     rho = gdps(ij) / ( rair * gdta(ij) &
!     &              * ( 1.0d0 + 0.608d0 * gdqa(ij) ) )
     exi = (gdps(ij) / (gdps(ij) - rho * grav * za))**akappa
     theta = gdta(ij) * exi
     thetav = theta * ( 1.0d0 + 0.608d0 * gdqa(ij) )

     uabs = sqrt(usfc(ij) * usfc(ij) + vsfc(ij) * vsfc(ij))
     uabs = min( max( uabs, usmine ), usmaxe )
!    '08.09.08 bug fix:
!     Now the humidity factor (= 0.98) will be mulitiplied later,
!      since the value depend on surface condition (namely, 0.98 or 1).
!     qsat = 0.98d0 * fqsat ( gdts( ij ), gdps(ij) )
     qsat = fqsat ( gdts( ij ), gdps(ij) )

!    First guess for bulk coefficients
!    '08.06.05
!     Bulk coeff. are now diagnosed without distinction of ocean or ice 
     if (l == 0) then     !! oceanic surface
        qsat = 0.98d0 * qsat
        cdn10 = ( 2.70d0 / uabs + 0.142d0 + uabs / 13.09d0 ) * 1.0d-3
        cdn10r = sqrt( cdn10 )
        cen10 = 34.6d-3 * cdn10r         
        astbl = 0.5d0 + sign( 0.5d0, theta - gdts(ij) )
        chn10 = ( 18.0d-3 * astbl + 32.7d-3 * (1.0d0 - astbl) ) * cdn10r
        niter = nitero
     else   !! sea-ice surface
        cdn10 = 1.63d-3
        cdn10r = sqrt( cdn10 )
        cen10 = 1.63d-3
        chn10 = 1.63d-3
        niter = niteri
     end if
     cmv(ij) = cdn10
     cev(ij) = cen10
     chv(ij) = chn10

!    Iteration for determing the coefficients
!    The same height, 10 m above the sea surface, is assumed
!      for the atmospheri! u, theta, and q.
!    In that case...
!     - zeta for u, theta, and q have the same value
!     - psi_h for u, theta, and q have the same value
!     - no need for shifting u, theta and q to 10 m height

     do n=1, niter
!       '08.09.19 bug fix: sqrt
        cmvsqr = sqrt( cmv(ij) )
        ustar = cmvsqr * uabs
        tstar = chv(ij) / cmvsqr * (theta - gdts(ij))
        qstar = cev(ij) / cmvsqr * (gdqa(ij) - qsat)

        zeta = ckarm * grav * za / ( ustar * ustar ) &
          &    * ( tstar / thetav &
          &      + qstar / ( gdqa(ij) + 1.0d0/0.608d0 ) )

!       --- from Dr. T.Suzuki's modification dated 080425 ---
!       ---  originally from NCAR's code ? ---
        zeta = sign( min( abs(zeta), 10.d0 ), zeta)

        astbl = 0.5d0 + sign(0.5d0, zeta)
        x2 = sqrt( 1.0d0 - 16.0d0 * min(zeta, 0.0d0) )
        x = sqrt( x2 )
        psi_m = ( -5.0d0 * zeta ) * astbl &
     &          + ( 2.0d0 * log( ( 1.0d0 + x ) * 0.5d0 ) &
     &              + log( (1.0d0 + x2 ) * 0.5d0 ) &
     &              - 2.0d0 * atan(x) + pi * 0.5d0 ) &
     &                * ( 1.0d0 - astbl )
        psi_h = ( -5.0d0 * zeta ) * astbl &
     &          + ( 2.0d0 * log( ( 1.0d0 + x2 ) * 0.5d0 )) &
     &              * ( 1.0d0 - astbl )

!       '08.09.08 bug fix:
!        The procedure for updating the transfar coefficients
!          was incorrect in the old code.
!        The coefficients must be diagnosed based on
!          neutral 10 m coefficients (i.e., CDN10 etc.).
!        --- the old code from here ---
!        xx = 1.0d0 / ( 1.0d0 - cmvsqr * psi_m / ckarm )
!        un10 = uabs * xx
!        cmv(ij) = cmv(ij) * xx * xx
!!       bug fix 080422
!!        xx = xx / ( 1.0d0 - cmvsqr * psi_h / ckarm )
!!        chv(ij) = chv(ij) * xx
!!        cev(ij) = cev(ij) * xx
!        phrcdk = psi_h / (cmvsqr * ckarm)
!        chv(ij) = chv(ij) * xx / ( 1.0d0 - chv(ij) * phrcdk )
!        cev(ij) = cev(ij) * xx / ( 1.0d0 - cev(ij) * phrcdk )
!        cmvsqr = sqrt( cmv(ij) )
!        --- the end of the old code ---
        if (l == 0) then   !! re-calculate if oceanic surface
           un10 = uabs / ( 1.0d0 - cdn10r * psi_m / ckarm )
           cdn10 = ( 2.70d0 / un10 + 0.142d0 + un10 / 13.09d0 ) * 1.0d-3
           cdn10r = sqrt( cdn10 )
           cen10 = 34.6d-3 * cdn10r         
           astbl = 0.5d0 + sign(0.5d0, zeta)
           chn10 = ( 18.0d-3 * astbl + 32.7d-3 * (1.0d0 - astbl) ) * cdn10r
        endif
            
        xx = 1.0d0 / ( 1.0d0 - cdn10r * psi_m / ckarm )
        cmv(ij) = cdn10 * xx * xx
        phrcdk = psi_h / (cdn10r * ckarm)
        chv(ij) = chn10 * xx / ( 1.0d0 - chn10 * phrcdk )
        cev(ij) = cen10 * xx / ( 1.0d0 - cen10 * phrcdk )
     enddo
  enddo

  do ij = ijstr, ijend
     cmv( ij ) = max( min( cmv( ij ), cmmax ), cmmin )
     chv( ij ) = max( min( chv( ij ), chmax ), chmin )
     cev( ij ) = max( min( cev( ij ), cemax ), cemin )
  end do

  do ij = ijstr, ijend
     uabs      = sqrt(usfc(ij) * usfc(ij) + vsfc(ij) * vsfc(ij))
     cmv( ij ) = cmv( ij ) * min( max( uabs, usminm ), usmaxm )
     chv( ij ) = chv( ij ) * min( max( uabs, usminh ), usmaxh )
     cev( ij ) = cev( ij ) * min( max( uabs, usmine ), usmaxe )
  end do

  return
end subroutine blkcof_core

#ifdef OPT_BODY
! *********************************************************************

subroutine bdyflx( &
  &                    tq, &
  &                     t )

! --- information -----------------------------------------------------
!
!  Estimate the surface forcing fluxes from the data given in the
! boundary condition files
!
!  HISTORY
!     '99.08.25  H.Hasumi
!     '02.06.02  H.Hasumi: tracer dimension
!     '12.10.11  Y.Komuro: for COCO5.0
!
! ---------------------------------------------------------------------

  use utint

  real(8), intent(out) ::      tq(nxydim, nzdim, ntdim)
  real(8), intent(in)  ::       t(nxydim, nzdim, ntdim)

  integer ::     ij,      k,      l,      n
  
  real(8), save        ::    tbdy(nxydim, nzdim, ntdim) = 0.d0
  real(8), save        ::    tdmb(nxydim, nzdim, ntdim) = 0.d0

!  do l = 1, ntdim
  do l = 1, 2
     n = l
     call tmintb( &
       &            tbdy(1, 1, l), &
       &               n )
  end do
!  do l = 1, ntdim
!     n = l + ntdim
  do l = 1, 2
     n = l + 2
     call tmintb( &
       &            tdmb(1, 1, l), &
       &               n )
  end do

!  do l = 1, ntdim
  do l = 1, 2
     do k = kstr, kend
        do ij = ijtstr, ijtend
           tq(ij, k, l) = tdmb(ij, k, l) * &
             &            (tbdy(ij, k, l) - t(ij, k, l)) * &
             &            amskt(ij, k)
        end do
     end do
  end do

  return
end subroutine bdyflx
#endif

end module sfcng
