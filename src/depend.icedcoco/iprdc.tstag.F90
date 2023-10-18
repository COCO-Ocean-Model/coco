module iprdc

! --- information -----------------------------------------------------
!
!  HISTORY
!     '02.10.10  H.Hasumi: from MIROC3.1
!     '07.04.23  H.Hasumi
!     '07.09.25  H.Hasumi: multi-category sea ice
!     '07.09.28  H.Hasumi: 1-layer sea ice thermodynamics
!     '08.07.10  H.Hasumi: initial/final processing
!     '09.02.20  Y.Komuro: parallel forward time march
!     '09.05.25  Y.Komuro: CMIP5 output code included
!     '09.09.04  Y.Komuro: extra output code (FEX/FEY)
!     '09.09.26  Y.Komuro: bug fix
!     '10.04.14  M.Kurogi: staggered time stepping
!     '10.04.14  M.Kurogi: (COCO4.4 tripolar code by Dr. Suzuki)
!     '12.07.30  Y.Komuro: for COCO5.0
!     '13.02.12  Y.Komuro: remove non-parallel code 
!     '21.05.26  Y.Komuro: snow aging & meltpond parametrization 
!
! ---------------------------------------------------------------------

  use zocdim, only: &
    & nxydim,  nxdim,  nydim,  nzdim,  ntdim,    nic, nxyidm,   kstr, &
    &     nx,     ny, ijtstr, ijtend, &
    &  oinit, ofinal, &
    & myrank, ijnode
  use zocgrd, only: &
    &     ts
  use zocmsk, only: &
    &  amskt,  amskv
  use zocphy, only: &
    & gravit,   rhoi,   rhos

  implicit none

  integer, parameter :: nrbnd = 3  !! N. of radiation bands; VI, NIR, and IR

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

  namelist /nmsage/   osage, alssif, alssio, alfmax, snrfrs,  ftage, &
    &                tauage, adirt0, adirtc, adirts, adirtm, drsmax, &
    &                 oadst
  namelist /nmmpnd/   impnd, hminmp, rtdpmp, rtmxmp,  dpscl, &
    &                rmpcmn, rmpcmx, cmpfrz, tmpfrz, albmpd, almpdp, &
    &                frmpmn, vmpmin

  private

  public :: predci

contains

subroutine predci( &
  &                    ax,    hix,    uix,    vix,    tix,    hsx, &
  &                   asx,  frlvx,   vmpx,  frmpx,   dsdx,   dsbx, &
  &                    ft,     fs,   taux,   tauy,   ptop, &
  &                    ay,    hiy,    uiy,    viy,    tiy,    hsy, &
  &                   asy,  frlvy,   vmpy,  frmpy,   dsdy,   dsby, &
  &                    tx,     ux,     vx,     hx,     hy, &
  &                   qao,    qai,    qii,    qio,  swabs,    tsi, &
  &                   wev,    wsb, &
  &                  prec,   snow,   roff,   soff, &
  &                  dfdu,   dfbc, &
  &                tauaix, tauaiy, tauaox, tauaoy )

  use ictrn
  use ifhea
  use ifwtr
  use ipadv
  use ipmmt
  use iprdg
  use ipthm
  use iptmp
  use ucloc
  use brstt
  use qckag
  use qckot
  use ufile
  use bshft

  real(8), intent(inout) ::     ax(nxydim, 0:nic)
  real(8), intent(inout) ::    hix(nxydim, 0:nic)
  real(8), intent(inout) ::    uix(nxydim)
  real(8), intent(inout) ::    vix(nxydim)
  real(8), intent(inout) ::    tix(nxydim, 0:nic)
  real(8), intent(inout) ::    hsx(nxydim, 0:nic)
  
  real(8), intent(out)   ::     ft(nxydim, ntdim),     fs(nxydim)
  real(8), intent(out)   ::   taux(nxydim),   tauy(nxydim)
  real(8), intent(out)   ::   ptop(nxydim)

  real(8), intent(in)    ::     ay(nxydim, 0:nic)
  real(8), intent(in)    ::    hiy(nxydim, 0:nic)
  real(8), intent(in)    ::    uiy(nxydim)
  real(8), intent(in)    ::    viy(nxydim)
  real(8), intent(in)    ::    tiy(nxydim, 0:nic)
  real(8), intent(in)    ::    hsy(nxydim, 0:nic)
  real(8), intent(in)    ::     tx(nxydim, nzdim, ntdim)
  real(8), intent(in)    ::     ux(nxydim, nzdim),     vx(nxydim, nzdim)
  real(8), intent(in)    ::     hx(nxydim),     hy(nxydim)
  real(8), intent(inout) ::    qio(nxydim, nic),    qai(nxydim, nic)
  real(8), intent(inout) ::    qii(nxydim, nic)
  real(8), intent(inout) ::    qao(nxydim)
  real(8), intent(in)    ::  swabs(nxydim)
  real(8), intent(in)    ::    wev(nxydim)
  real(8), intent(inout) ::    wsb(nxydim, nic)
  real(8), intent(inout) ::   prec(nxydim),   snow(nxydim)
  real(8), intent(in)    ::   roff(nxydim),   soff(nxydim)
  real(8), intent(inout) :: tauaix(nxydim), tauaiy(nxydim)
  real(8), intent(inout) :: tauaox(nxydim), tauaoy(nxydim)

!---- arrays for sage.mp (snow aging & melt pond)
  real(8), intent(inout) ::    asx(nxydim, 0:nic)
  real(8), intent(inout) ::  frlvx(nxydim, 0:nic)
  real(8), intent(inout) ::   vmpx(nxydim, 0:nic)
  real(8), intent(inout) ::  frmpx(nxydim, 0:nic)
  real(8), intent(inout) ::   dsdx(nxydim, 0:nic)
  real(8), intent(inout) ::   dsbx(nxydim, 0:nic)
  real(8), intent(out)   ::    asy(nxydim, 0:nic)
  real(8), intent(out)   ::  frlvy(nxydim, 0:nic)
  real(8), intent(out)   ::   vmpy(nxydim, 0:nic)
  real(8), intent(out)   ::  frmpy(nxydim, 0:nic)
  real(8), intent(out)   ::   dsdy(nxydim, 0:nic)
  real(8), intent(out)   ::   dsby(nxydim, 0:nic)
  real(8), intent(in)    ::   dfdu(nxydim),   dfbc(nxydim)
  real(8), intent(in)    ::    tsi(nxydim, 0:nic)

  real(8), save ::    wao(nxydim)
  real(8), save ::    wio(nxydim, nic),    was(nxydim, nic)
  real(8), save ::    wil(nxydim, nic)
  real(8), save ::   subi(nxydim, nic)
  real(8), save ::   evap(nxydim)
  real(8), save :: wiadjs(nxydim), weadjs(nxydim)
  real(8), save ::     az(nxydim, 0:nic),    hiz(nxydim, 0:nic)
  real(8), save ::    hsz(nxydim, 0:nic),    tiz(nxydim, 0:nic)
  real(8), save ::   pice(nxydim)
  real(8), save ::  frmpz(nxydim, 0:nic)

  real(8), save ::    eix(nxydim, 0:nic),    eiz(nxydim, 0:nic)

  real(8), save ::    fix(nxydim, 0:nic),    fiy(nxydim, 0:nic)
  real(8), save ::    fsx(nxydim, 0:nic),    fsy(nxydim, 0:nic)
  real(8), save ::    fex(nxydim, 0:nic),    fey(nxydim, 0:nic)
  real(8), save ::  ftitd(nxydim)
  real(8), save :: igrfra(nxydim), igrcon(nxydim), igrsni(nxydim)
  real(8), save :: igrsfl(nxydim), inrlat(nxydim)
  real(8), save :: imrsno(nxydim), imrsmi(nxydim)
  real(8), save :: imrisf(nxydim), imribs(nxydim)
  real(8), save :: inrsbi(nxydim), inrsbs(nxydim)
  real(8), save :: imraji(nxydim), imrajs(nxydim)
  real(8), save :: impinc(nxydim, 0:nic), impfrz(nxydim, 0:nic)
  real(8), save :: improf(nxydim, 0:nic)
  real(8), save ::    fdd(nxydim),    fdb(nxydim)
  real(8), save :: sitfrc(nxydim), siuabs(nxydim)

  logical, save ::  oeof

!! for check
!  real(8) :: imrtot(nxydim), igrtot(nxydim), itrtot(nxydim)

  integer ::     ij,      k,      l
  integer ::  ifpar,  jfpar,  istat

  call clcstr('ICE')

  if (      (myrank .ge. ijnode) &
    & .and. (.not. oinit) .and. (.not. ofinal)) then
     return
  end if

  if (oinit) then
     call rewnml(ifpar, jfpar)
     read (ifpar, nmsage, iostat=istat)
     call cstnml(jfpar, 'predci', 'nmsage', istat)
     write(jfpar, nmsage)
     call rewnml(ifpar, jfpar)
     read (ifpar, nmmpnd, iostat=istat)
     call cstnml(jfpar, 'predci', 'nmmpnd', istat)
     write(jfpar, nmmpnd)

     do l = 0, nic
        do ij = 1, nxydim
           asx(ij, l) = 0.0d0 !! assume fresh snow
           frlvx(ij, l) = 1.0d0 !! assume level ice
           vmpx(ij, l) = 0.0d0
           frmpx(ij, l) = 0.0d0
           dsdx(ij, l) = 0.0d0   !! assume no dust
           dsbx(ij, l) = 0.0d0   !! assume no dust
        end do
     end do
#ifdef OPT_TRIPOLE
     call rstadd(pice, oeof, nxdim, nydim, 1, 'PICE', 'SFC', &
       &                                      1.d0,  0,  0 )
     if (osage) then
        call rstadd( asx, oeof, nxdim, nydim, nic+1, 'AS', 'ICE', &
          &                                      1.d0,  0,  0 )
     end if
     if (impnd > 0) then
        call rstadd(frlvx, oeof, nxdim, nydim, nic+1, 'FRLV', 'ICE', &
          &                                      1.d0,  0,  0 )
        call rstadd(vmpx, oeof, nxdim, nydim, nic+1, 'VMP', 'ICE', &
          &                                      1.d0,  0,  0 )
        call rstadd(frmpx, oeof, nxdim, nydim, nic+1, 'FRMP', 'ICE', &
          &                                      1.d0,  0,  0 )
     end if
     if (oadst) then
        call rstadd(dsdx, oeof, nxdim, nydim, nic+1, 'DSD', 'ICE', &
          &                                      1.d0,  0,  0 )
        call rstadd(dsbx, oeof, nxdim, nydim, nic+1, 'DSB', 'ICE', &
          &                                      1.d0,  0,  0 )
     end if
#else
     call rstadd(pice, oeof, nxdim, nydim, 1, 'PICE', 'SFC')
     if (osage) then
        call rstadd( asx, oeof, nxdim, nydim, nic+1, 'AS'  , 'ICE')
     end if
     if (impnd > 0) then
        call rstadd(frlvx, oeof, nxdim, nydim, nic+1, 'FRLV' , 'ICE')
        call rstadd(vmpx, oeof, nxdim, nydim, nic+1, 'VMP' , 'ICE')
        call rstadd(frmpx, oeof, nxdim, nydim, nic+1, 'FRMP' , 'ICE')
     end if
     if (oadst) then
        call rstadd(dsdx, oeof, nxdim, nydim, nic+1, 'DSD' , 'ICE')
        call rstadd(dsbx, oeof, nxdim, nydim, nic+1, 'DSB' , 'ICE')
     end if
#endif
     call idfrmp( &
       &          frlvx,   vmpx,  frmpx, &
       &         improf, &
       &             ax,    hix,    hsx)
     if (oeof) then
        do l = 0, nic
           do ij = 1, nxydim
              az  (ij, l) =  ax  (ij, l)
              hiz (ij, l) =  hix (ij, l)
              eiz (ij, l) =  eix (ij, l)
              hsz (ij, l) =  hsx (ij, l)
              tiz (ij, l) =  tix (ij, l) 
           end do
        end do
        call pridge( &
          &            pice, &
          &              ax,    hix,    eix,    hsx,    tix, &
          &             asx,  frlvx,   vmpx,  frmpx,   dsdx,   dsbx, & 
          &              az,    hiz,    eiz,    hsz, &
          &             uix,    vix )
        do l = 0, nic
           do ij = 1, nxydim
              ax  (ij, l) =  az  (ij, l)
              hix (ij, l) =  hiz (ij, l)
              eix (ij, l) =  eiz (ij, l)
              hsx (ij, l) =  hsz (ij, l)
              tix (ij, l) =  tiz (ij, l) 
           end do
        end do
        call idfrmp( &
        &          frlvx,   vmpx,  frmpx, &
        &         improf, &
        &             ax,    hix,    hsx) 
        do l = 0, nic
          do ij = 1, nxydim
             asy(ij, l) = asx(ij, l)
             frlvy(ij, l) = frlvx(ij, l)
             frmpy(ij, l) = frmpx(ij, l)
             vmpy(ij, l) = vmpx(ij, l)
             dsdy(ij, l) = dsdx(ij, l)
             dsby(ij, l) = dsbx(ij, l)
          end do
       end do         
     end if
  end if
      
  if (ofinal) then
     call finadd(pice, nxdim, nydim, 1, 'PICE', 'SFC')
     if (osage) then
        call finadd(asx,  nxdim, nydim, nic+1, 'AS'  , 'ICE')
     end if
     if (impnd > 0) then
        call finadd(frlvx, nxdim, nydim, nic+1, 'FRLV' , 'ICE')
        call finadd(vmpx, nxdim, nydim, nic+1, 'VMP' , 'ICE')
        call finadd(frmpx, nxdim, nydim, nic+1, 'FRMP' , 'ICE')
     end if
     if (oadst) then
        call finadd(dsdx, nxdim, nydim, nic+1, 'DSD' , 'ICE')
        call finadd(dsbx, nxdim, nydim, nic+1, 'DSB' , 'ICE')
     end if
  end if

  do ij = 1, nxydim
     evap(ij) = wev(ij)
     fs  (ij) = 0.d0
     taux(ij) = 0.d0
     tauy(ij) = 0.d0
     az  (ij, 0) = 1.d0
     hiz (ij, 0) = 0.d0
     hsz (ij, 0) = 0.d0
     fdd(ij) = 0.0d0
     fdb(ij) = 0.0d0
     wiadjs(ij) = 0.d0
     weadjs(ij) = 0.d0
     imraji(ij) = 0.d0
     imrajs(ij) = 0.d0
     igrfra(ij) = 0.d0
     igrcon(ij) = 0.d0
     igrsni(ij) = 0.d0
     igrsfl(ij) = 0.d0
     inrlat(ij) = 0.d0
     imrsno(ij) = 0.d0
     imrsmi(ij) = 0.d0
     imrisf(ij) = 0.d0
     imribs(ij) = 0.d0
     inrsbi(ij) = 0.d0
     inrsbs(ij) = 0.d0
  end do
  do l = 1, nic
     do ij = 1, nxydim
        subi(ij, l) = 0.d0
        az  (ij, l) = 0.d0
        hiz (ij, l) = 0.d0
        hsz (ij, l) = 0.d0
        improf(ij, l) = 0.0d0
      end do
  end do
  do l = 1, 2
     do ij = 1, nxydim
        ft(ij, l) = 0.d0
     end do
  end do

  call clcstr('ICEDYN')
  call pmomnt( &
    &            uix,    vix, &
    &           taux,   tauy, &
    &             ax,     ay,     az, &
    &            hix,    hiy,    hiz, &
    &            hsx,    hsy,    hsz, &
    &            uiy,    viy,   pice, &
    &             ux,     vx,     hy,   ptop, &
    &         tauaix, tauaiy, tauaox, tauaoy )
#ifdef OPT_TRIPOLE
  call shift2( &
    &           taux,    tauy, &
    &          nxdim,   nydim,     1, &
    &          -1.d0,      -1,    -1 )
#else
  call shift2( &
    &           taux,   tauy, &
    &          nxdim,  nydim,      1 )
#endif
  call clcend('ICEDYN')

  call ipsage( &
    &            asx, &
    &            hsx,    tsi,    snow,   dsdx,   dsbx)

  call icetmp( &
    &            eix, &
    &            tix, &
    &            qao,    qai,    qio,    qii, &
    &             ax,    hix )

  call fiheat( &
    &             ft, &
    &            wao,    wio,    was,    wil, &
    &             ax,     tx,     hx, &
    &            qao,    qai,    qio,    qii,  swabs )
  call fwater( &
    &             ax,    hix,    hsx,    eix,    tix, &
    &            asx,  frlvx,   vmpx,  frmpx,   dsdx,   dsbx, &
    &           prec,   snow,    fdd,    fdb, &
    &         inrsbi, inrsbs, imraji, imrajs, &
    &           subi, wiadjs, weadjs, &
    &            wev,    wsb,   soff )

  call cofpfw( &
    &           prec,   snow,   roff,   soff,   evap, &
    &             ax)

  call ptherm( &
    &             ax,    hix,    hsx,    eix,    tix, &
    &            asx,  frlvx,   vmpx,  frmpx,   dsdx,   dsbx, &
    &           prec,   snow,     ft,     fs,    fdd,    fdb, &
    &          ftitd, igrfra, igrcon, igrsni, &
    &         igrsfl, inrlat, &
    &         imrsno, imrsmi, imrisf, imribs, &
    &         impinc, impfrz, improf, &
    &             tx,    tsi, &
    &            wio,    wao,    was,    wil, &
    &           evap,   subi,   roff, wiadjs, weadjs, &
    &           dfdu,   dfbc, &
    &            qio )
  call idfrmp( &
    &          frlvx,   vmpx,  frmpx, &
    &         improf, &
    &             ax,    hix,    hsx)
  call ictrns( &
    &             ax,    hix,    hsx,    eix,    tix, &
    &            asx,  frlvx,   vmpx,  frmpx,   dsdx,   dsbx, &
    &             ft,     fs,    fdd,    fdb, &
    &         imraji, imrajs )
  call idfrmp( &
    &          frlvx,   vmpx,  frmpx, &
    &         improf, &
    &             ax,    hix,    hsx)

#ifdef OPT_TRIPOLE
  call shift3( &
    &             ax,    hix,    hsx, &
    &          nxdim,  nydim,  nic+1, &
    &          1.0d0,      0,      0 )
  call shift3( &
    &            eix,    tix,    asx, &
    &          nxdim,  nydim,  nic+1, &
    &          1.0d0,      0,      0 )
  call shift3( &
    &          frlvx,   vmpx,  frmpx, &
    &          nxdim,  nydim,  nic+1, &
    &          1.0d0,      0,      0 )
  call shift2( &
    &           dsdx,   dsbx, &
    &          nxdim,  nydim,  nic+1, &
    &          1.0d0,      0,      0 )
#else
  call shift3( &
    &             ax,    hix,    hsx, &
    &          nxdim,  nydim,  nic+1 )
  call shift3( &
    &            eix,    tix,    asx, &
    &          nxdim,  nydim,  nic+1)
  call shift3( &
    &          frlvx,   vmpx,  frmpx, &
    &          nxdim,  nydim,  nic+1)
  call shift2( &
    &           dsdx,   dsbx, &
    &          nxdim,  nydim,  nic+1)
#endif

  call padvct( &
    &             ax,    hix,    eix,    hsx,    tix, &
    &            asx,  frlvx,   vmpx,  frmpx,   dsdx,   dsbx, &
    &             az,    hiz,    eiz,    hsz, &
    &            fix,    fiy,    fsx,    fsy, &
    &            fex,    fey, &
    &            uix,    vix )
  call pridge( &
    &           pice, &
    &             ax,    hix,    eix,    hsx,    tix, &
    &            asx,  frlvx,   vmpx,  frmpx,   dsdx,   dsbx, &
    &             az,    hiz,    eiz,    hsz, &
    &            uix,    vix )
  call icadjs( &
    &             ax,    hix,    hsx,    eix, &
    &           vmpx,   dsdx,   dsbx )
  call ichflt( &
    &             ax,    hix,    hsx,    eix,    tix, &
    &            asx,  frlvx,   vmpx,  frmpx,   dsdx,   dsbx )  
  call ictrns( &
    &             ax,    hix,    hsx,    eix,    tix, &
    &            asx,  frlvx,   vmpx,  frmpx,   dsdx,   dsbx, &
    &             ft,     fs,    fdd,    fdb, &
    &         imraji, imrajs )
  call idfrmp( &
    &          frlvx,   vmpx,  frmpx, &
    &         improf, &
    &             ax,    hix,    hsx)
#ifdef OPT_TRIPOLE
  call shift3( &
    &             ax,    hix,    hsx, &
    &          nxdim,  nydim,  nic+1, &
    &          1.0d0,      0,      0 )
  call shift3( &
    &            eix,    tix,    asx, &
    &          nxdim,  nydim,  nic+1, &
    &          1.0d0,      0,      0 )
  call shift3( &
    &          frlvx,   vmpx,  frmpx, &
    &          nxdim,  nydim,  nic+1, &
    &          1.0d0,      0,      0 )
  call shift2( &
    &           dsdx,   dsbx, &
    &          nxdim,  nydim,  nic+1, &
    &          1.0d0,      0,      0 )
#else
  call shift3( &
    &             ax,    hix,    hsx, &
    &          nxdim,  nydim,  nic+1 )
  call shift3( &
    &            eix,    tix,    asx, &
    &          nxdim,  nydim,  nic+1)
  call shift3( &
    &          frlvx,   vmpx,  frmpx, &
    &          nxdim,  nydim,  nic+1)
  call shift2( &
    &           dsdx,   dsbx, &
    &          nxdim,  nydim,  nic+1)
#endif


  if (myrank .ge. ijnode) then
     return
  end if

  do l = 1, nic
     do ij = 1, nxydim
        ptop(ij) = ptop(ij) &
          &      + gravit * ay(ij, l) * &
          &        (rhoi * hiy(ij, l) + rhos * hsy(ij, l))
     end do
  end do

  call clcend('ICE')

  if (oinit .or. ofinal) then
     return
  end if

! output section for CMIP5
! FIX, FIY: eastward/northward sea ice transport
!           unit [cm^3(ice)/s]
  call chekin(   fix,  'FIX', &
       &            'eastward ice transport', 'cm^3(ice)/rad/s', &
       &             nx,     ny,    nic, nxyidm, 'OCICEX')
  call chekin(   fiy,  'FIY', &
     &            'northward ice transport', 'cm^3(ice)/rad/s', &
    &             nx,     ny,    nic, nxyidm, 'OCICEY')
! FISX, FISY: eastward/northward snow transport
!             unit [cm^3(snow)/s]
  call chekin(   fsx, 'FISX', &
     &            'eastward snow transport', 'cm^3(snow)/rad/s', &
    &             nx,     ny,    nic, nxyidm, 'OCICEX')
  call chekin(   fsy, 'FISY', &
     &            'northward snow transport', 'cm^3(snow)/rad/s', &
    &             nx,     ny,    nic, nxyidm, 'OCICEY')
! FTITD: heat flux into sea water due to sea ice thermodynamics
!        unit [erg/cm^2/s]
  call chekin(  ftitd,  'FTITD', &
     &     'heat flux into sea water due to sea ice thermodynamics', &
     &                 'erg/cm^2/s', &
    &              nx,     ny,      1, nxydim, 'OCSFCT')
! IGRFRA: frazil sea ice growth rate, unit [g/cm^2/s]
  call chekin( igrfra, 'IGRFRA', &
     &            'frazil sea ice growth rate', 'g/cm^2/s', &
    &              nx,     ny,      1, nxydim, 'OCSFCT')
! IGRCON: congelation sea ice growth rate, unit [g/cm^2/s]
  call chekin( igrcon, 'IGRCON', &
     &            'congelation sea ice growth rate', 'g/cm^2/s', &
    &              nx,     ny,      1, nxydim, 'OCSFCT')
! IGRSNI: snow-ice formation rate, unit [g/cm^2/s]
  call chekin( igrsni, 'IGRSNI', &
     &            'snow-ice formation rate', 'g/cm^2/s', &
    &              nx,     ny,      1, nxydim, 'OCSFCT')
! IGRSFL: snowfall rate on ice, unit [g/cm^2/s]
    call chekin( igrsfl, 'IGRSFL', &
    &            'snowfall rate on ice', 'g/cm^2/s', &
   &              nx,     ny,      1, nxydim, 'OCSFCT')
! INRLAT: lateral sea ice net growth rate, unit [g/cm^2/s]
!  (COCO4.4 represents lateral melting process but not freezing,
!   thus this value will be zero or negative.)
  call chekin( inrlat, 'INRLAT', &
     &            'lateral sea ice net growth rate', 'g/cm^2/s', &
    &              nx,     ny,      1, nxydim, 'OCSFCT')
! IMRSNO: snow melt rate, unit [g/cm^2/s]
  call chekin( imrsno, 'IMRSNO', &
     &            'snow melt rate', 'g/cm^2/s', &
    &              nx,     ny,      1, nxydim, 'OCSFCT')
! IMRSMI: snow decrease rate due to ice melt, unit [g/cm^2/s]
    call chekin( imrsmi, 'IMRSMI', &
    &            'snow decrease rate due to ice melt', 'g/cm^2/s', &
    &              nx,     ny,      1, nxydim, 'OCSFCT')
! IMRISF: rate of melt at upper surface of sea ice, unit [g/cm^2/s]
  call chekin( imrisf, 'IMRISF', &
     &            'rate of melt at upper surface of sea ice', &
     &                   'g/cm^2/s', &
    &              nx,     ny,      1, nxydim, 'OCSFCT')
! IMRIBS: rate of melt at sea ice base, unit [g/cm^2/s]
  call chekin( imribs, 'IMRIBS', &
     &            'rate of melt at sea ice base', 'g/cm^2/s', &
    &              nx,     ny,      1, nxydim, 'OCSFCT')
! INRSBI: net rate of ice sublimation, unit [g/cm^2/s]
  call chekin(inrsbi, 'INRSBI', &
    &         'rate of ice sublimation', 'g/cm^2/s', &
    &         nx, ny, 1, nxydim, 'OCSFCT')
! INRSBS: net rate of snow sublimation, unit [g/cm^2/s]
  CALL chekin(inrsbs, 'INRSBS', &
    &         'rate of snow sublimation', 'g/cm^2/s', &
    &         nx, ny, 1, nxydim, 'OCSFCT')
      
! IMRAJI: rate of ice melt in adjustment, unit [g/cm^2/s]
  CALL chekin(imraji, 'IMRAJI', &
    &         'rate of ice melt in adjustment', 'g/cm^2/s', &
    &         nx, ny, 1, nxydim, 'OCSFCT')
! IMRAJS: rate of snow melt in adjustment, unit [g/cm^2/s]
  CALL chekin(imrajs, 'IMRAJS', &
    &         'rate of snow melt in adjustment', 'g/cm^2/s', &
    &         nx, ny, 1, nxydim, 'OCSFCT')

! extra output section
! FEX, FEY: eastward/northward thermal transport by sea ice
!           unit [(erg/g)*(cm^3(ice)/s)]
  call chekin(   fex,  'FIEX', &
    &            'eastward ice heat transport', 'erg*cm^3/rad/g/s', &
    &             nx,     ny,    nic, nxyidm, 'OCICEX')
  call chekin(   fey,  'FIEY', &
    &            'northward ice heat transport', 'erg*cm^3/rad/g/s', &
    &             nx,     ny,    nic, nxyidm, 'OCICEY')

! output section for CMIP6 
  do ij = 1, nxydim
     if (ax(ij, 0) .lt. 1.0d0) then
        sitfrc(ij) = 1.0d0 * amskt(ij, kstr)
     else
        sitfrc(ij) = 0.0d0
     end if
     siuabs(ij) = sqrt(uix(ij)*uix(ij) + vix(ij)*vix(ij)) &
       &        * amskv(ij, kstr)
  end do

! ITFRAC: Fraction of time steps with sea ice (i.e., =1 if AI>0)
  call chekin( sitfrc, 'ITFRAC', &
    &          'fraction of time steps with sea ice', 'ND', &
    &              nx,     ny,      1, nxydim, 'OCSFCT')
! UIABS: Sea-ice speed
  call chekin( siuabs, 'UIABS', &
    &          'sea-ice speed', 'cm/s', &
    &              nx,     ny,      1, nxydim, 'OCSFCV')

!!     for check: not necessary for CMIP5 output
!  do ij = 1, nxydim
!     igrtot(ij) = igrfra(ij)+igrcon(ij)+igrsni(ij)
!     imrtot(ij) = imrsno(ij)+imrisf(ij)+imribs(ij)-inrlat(ij)
!     itrtot(ij) = igrtot(ij)-imrtot(ij)
!  end do
!  call chekin( igrtot, 'IGRTOT', &
!    &              nx,     ny,      1, nxydim, 'SFC')
!  call chekin( imrtot, 'IMRTOT', &
!    &              nx,     ny,      1, nxydim, 'SFC')
!  call chekin( itrtot, 'ITRTOT', &
!    &              nx,     ny,      1, nxydim, 'SFC')

! IMPINC: MP mass increase by snow/ice melting and rainfall [g/cm^2/s]
! IMPFRZ: MP mass decrease by MP freezing [g/cm^2/s]
! IMPROF: MP mass runoff by permiability, negative freeboard, and
!            lack of available space (too large fraction) [g/cm^2/s]
  call chekin( impinc, 'IMPINC', &
    &          'meltpond mass increase', 'g/cm^2/s', &
    &              nx,     ny,    nic, nxyidm, 'OCICET')
  call chekin( impfrz, 'IMPFRZ', &
    &          'meltpond mass change by freezing', 'g/cm^2/s', &
    &              nx,     ny,    nic, nxyidm, 'OCICET')
  call chekin( improf, 'IMPROF', &
    &          'meltpond mass runoff', 'g/cm^2/s', &
    &              nx,     ny,    nic, nxyidm, 'OCICET')

! ODFBC: dust fall rate, BC, [g/cm^2/s]
  call chekin(   dfbc, 'ODFBC', &
    &          'ocean dust fall rate, BC', 'g/cm^2/s', &
    &              nx,     ny,      1, nxydim, 'OCSFCT')
! ODFBC: dust fall rate, dust, [g/cm^2/s]
  call chekin(   dfdu, 'ODFDU', &
    &          'ocean dust fall rate, dust', 'g/cm^2/s', &
    &              nx,     ny,      1, nxydim, 'OCSFCT')
! FDD: dust flux into ocn: positive upward [g/cm^2/s]
  call chekin(    fdd,  'FDD', &
    &          'dust flux into ocn., upward positive', 'g/cm^2/s', &
    &              nx,     ny,      1, nxydim, 'OCSFCT')
! FDD: BC flux into ocn: positive upward [g/cm^2/s]
  call chekin(    fdb,  'FDB', &
    &          'BC flux into ocn., upward positive', 'g/cm^2/s', &
    &              nx,     ny,      1, nxydim, 'OCSFCT')

  return

end subroutine predci

end module iprdc
