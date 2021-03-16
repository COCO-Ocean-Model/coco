module aocea

! --- information -----------------------------------------------------
!
!  HISTORY
!     '02.10.10  H.Hasumi: from COCO3.4
!     '07.04.23  H.Hasumi
!     '07.09.25  H.Hasumi: arguments of CHEKIN
!     '07.10.04  H.Hasumi: 1-layer sea ice thermodynamics
!     '08.06.11  H.Hasumi: initial/final processing
!     '09.05.25  Y.Komuro: CMIP5 output code included
!     '09.10.06  Y.Komuro: FORSTO before PREDCI
!     '12.10.19  T.Suzuki: for COCO5.0 in F90
!     '13.02.12  Y.Komuro: remove non-parallel code 
!
! ---------------------------------------------------------------------

  use zocdim, only: &
    &     nx,     ny,     nz, &
    & nxydim, nxyzdm, nxyidm,  ntdim,    nic, &
    & myrank, ijnode, &
    &  oinit, ofinal
  use zocgrd, only: &
    &     dt, &
    &     tt,     ts,    tss, &
    &     nt,    its,   itst,   ntss, &
    & ieuler
  use zocfil, only: &
    & nfomax
      
  implicit none

  character(len=16), save :: ctrnam(ntdim), cftnam(ntdim)
  character(len=32), save :: ctrtit(ntdim), cfttit(ntdim)
  character(len=16), save :: ctruni(ntdim), cftuni(ntdim)
  real(8), save ::    uadv(nxyzdm),   vadv(nxyzdm),   wadv(nxyzdm, 9)
  real(8), save ::     gxx(nxydim),    gyy(nxydim)

  private

  public :: ocstup, ocean

contains

subroutine ocstup ( &
  &                     ub,     vb,     tb, &
  &                     hb,   ubtb,   vbtb,      w,      r, &
  &                    amv,    ahv, &
  &                     ft,   ptop, &
  &                    dt1 ) 

  use binst
  use brdge
  use qckot
  use tovtr
  use tslvt
  use ufile

  real(8), intent(out) ::     tb(nxyzdm, ntdim)
  real(8), intent(out) ::     ub(nxyzdm)
  real(8), intent(out) ::     vb(nxyzdm)
  real(8), intent(out) ::     hb(nxydim)
  real(8), intent(out) ::   ubtb(nxydim)
  real(8), intent(out) ::   vbtb(nxydim)
  real(8), intent(out) ::      w(nxyzdm),      r(nxyzdm)
  real(8), intent(out) ::    amv(nxyzdm),    ahv(nxyzdm)
  real(8), intent(in)  ::     ft(nxydim, ntdim)
  real(8), intent(out) ::   ptop(nxydim)
  real(8), intent(in)  ::    dt1

  integer ::     ij,      l
  integer ::  ifpar,  jfpar

  call rewnml(ifpar, jfpar)
  write(jfpar, *) '*** ocstup ***'

  ctrnam(1) = 'T'
  ctrnam(2) = 'S'
  cftnam(1) = 'FT'
  cftnam(2) = 'FW'
  do l = 3, ntdim
     write(ctrnam(l), '(a6,i2.2)') 'TRACER', l
     write(cftnam(l), '(a6,i2.2)') 'TRCFLX', l
  end do

  ctrtit(1) = 'ocean temperature'
  ctrtit(2) = 'ocean salinity'
  ctruni(1) = 'degC'
  ctruni(2) = 'psu'
  cfttit(1) = 'sea surface temperature flux'
  cfttit(2) = 'sea surface freshwater flux'
  cftuni(1) = 'K cm/s'
  cftuni(2) = 'cm/s'
  do l = 3, ntdim
     ctrtit(l) = ''
     ctruni(l) = ''
     cfttit(l) = ''
     cftuni(l) = ''
  end do

  DT = DT1

  call rdgeo
  call svtset
  call ovtset(r, tb)
  call chkset

! *** Initialization for variables ***
  if (myrank < ijnode) then
     call iniset( &
       &             uadv,   vadv,   wadv,      r, &
       &               ub,     vb,     tb, &
       &               hb,   ubtb,   vbtb, &
       &                w,    amv,    ahv )
  end if

  do ij = 1, nxydim
     gxx(ij) = 0.d0
     gyy(ij) = 0.d0
     ptop(ij) = 0.d0
  end do

  return
end subroutine ocstup

! =====================================================================

subroutine ocean ( &
  &                    ua,     va,     ta, &
  &                    ha,   ubta,   vbta, &
  &                    ab,    hib,    uib,    vib,    tib,    hsb, &
  &                    ub,     vb,     tb, &
  &                    hb,   ubtb,   vbtb, &
  &                    aa,    hia,    uia,    via,    tia,    hsa, &
  &                     w,      r,    amv,    ahv, &
  &                    ft,  swabs,     fs, &
  &                  taux,   tauy,   ptop,    tsi, &
  &                   nt1,    tt1,  itst1,    ts1,   its1, &
  &                 ntss1,   tss1, &
  &                oflout, oflstk )

  use aprdc
  use bfrch
  use binst
  use iprdc
  use qckot
  use sfcng
  use tflxt
  
  real(8), intent(inout) ::      ta(nxyzdm, ntdim),     tb(nxyzdm, ntdim)
  real(8), intent(inout) ::      ua(nxyzdm),     ub(nxyzdm)
  real(8), intent(inout) ::      va(nxyzdm),     vb(nxyzdm)
  real(8), intent(inout) ::      ha(nxydim),     hb(nxydim)
  real(8), intent(inout) ::    ubta(nxydim),   ubtb(nxydim)
  real(8), intent(inout) ::    vbta(nxydim),   vbtb(nxydim)
  real(8), intent(inout) ::       w(nxyzdm),      r(nxyzdm)
  real(8), intent(inout) ::     amv(nxyzdm),    ahv(nxyzdm)

  real(8), intent(inout) ::     aa(nxydim, 0:nic),     ab(nxydim, 0:nic)
  real(8), intent(inout) ::    hia(nxydim, 0:nic),    hib(nxydim, 0:nic)
  real(8), intent(inout) ::    uia(nxydim),    uib(nxydim)
  real(8), intent(inout) ::    via(nxydim),    vib(nxydim)
  real(8), intent(inout) ::    tia(nxydim, 0:nic),    tib(nxydim, 0:nic)
  real(8), intent(inout) ::    hsa(nxydim, 0:nic),    hsb(nxydim, 0:nic)
  real(8), intent(inout) ::    tsi(nxydim, 0:nic)

  real(8), intent(inout) ::     ft(nxydim, ntdim)
  real(8), intent(inout) ::  swabs(nxydim),     fs(nxydim)
  real(8), intent(inout) ::   taux(nxydim),   tauy(nxydim)
  real(8), intent(inout) ::   ptop(nxydim)

  real(8), intent(in)    ::    tt1,    ts1,   tss1
  integer, intent(in)    ::    nt1,   its1,  itst1,  ntss1
  logical, intent(in)    :: oflout(nfomax), oflstk(nfomax)

  real(8), save ::      gx(nxyzdm),     gy(nxyzdm)
  real(8), save ::      xx(nxyzdm),     yy(nxyzdm)
  real(8), save ::     qao(nxydim)
  real(8), save ::     qai(nxydim, nic),    qio(nxydim, nic)
  real(8), save ::     qii(nxydim, nic)
  real(8), save ::     wev(nxydim),    wsb(nxydim, nic)
  real(8), save ::    prec(nxydim),   snow(nxydim)
  real(8), save ::    roff(nxydim),   soff(nxydim)
  real(8), save ::  tauaix(nxydim), tauaiy(nxydim)
  real(8), save ::  tauaox(nxydim), tauaoy(nxydim)

#ifdef OPT_BODY
  real(8), save ::      tq(nxyzdm, ntdim)
#endif

  real(8) ::    ssfc(nxydim)
  real(8) ::     aig(nxydim),    hig(nxydim),    hsg(nxydim)

  integer ::     ij,      l

  if (oinit) then
     call predci ( &
       &               ab,    hib,    uib,    vib,    tib,    hsb, &
       &               ft,     fs,   taux,   tauy,   ptop, &
       &               aa,    hia,    uia,    via,    tia,    hsa, &
       &               tb,     ub,     vb,     hb,     ha, &
       &              qao,    qai,    qii,    qio,  swabs, &
       &              wev,    wsb, &
       &             prec,   snow,   roff,   soff, &
       &           tauaix, tauaiy, tauaox, tauaoy )
     call predco( &
       &               hb,   ubtb,   vbtb,      w,      r, &
       &               ub,     vb,     tb, &
       &              amv,    ahv, &
       &             taux,   tauy,   ptop, &
       &               ft,  swabs,     fs,   ssfc, &
#ifdef OPT_BODY
       &               tq, &
#endif
       &               ua,     va,     ta, &
       &               ha,   ubta,   vbta,     ab, &
       &               gx,     gy,     xx,     yy, &
       &              gxx,    gyy, &
       &             uadv,   vadv,   wadv )
     call putsig( &
       &               tb )
     return
  end if

  if (ofinal) then
     if (myrank < ijnode) then
        call iniset( &
          &            uadv,   vadv,   wadv,      r, &
          &              ub,     vb,     tb, &
          &              hb,   ubtb,   vbtb, &
          &               w,    amv,    ahv )
     end if
     call predci( &
       &               ab,    hib,    uib,    vib,    tib,    hsb, &
       &               ft,     fs,   taux,   tauy,   ptop, &
       &               aa,    hia,    uia,    via,    tia,    hsa, &
       &               tb,     ub,     vb,     hb,     ha, &
       &              qao,    qai,    qii,    qio,  swabs, &
       &              wev,    wsb, &
       &             prec,   snow,   roff,   soff, &
       &           tauaix, tauaiy, tauaox, tauaoy )
     call predco( &
       &               hb,   ubtb,   vbtb,      w,      r, &
       &               ub,     vb,     tb, &
       &              amv,    ahv, &
       &             taux,   tauy,   ptop, &
       &               ft,  swabs,     fs,   ssfc, &
#ifdef OPT_BODY
       &               tq, &
#endif
       &               ua,     va,     ta, &
       &               ha,   ubta,   vbta,     ab, &
       &               gx,     gy,     xx,     yy, &
       &              gxx,    gyy, &
       &             uadv,   vadv,   wadv )
     return
  end if

  tt    = tt1
  nt    = nt1
  ts    = ts1
  its   = its1
  itst  = itst1
  ntss  = ntss1
  tss   = tss1

  call chkstk( &
    &          oflstk )

! *** Leap-frog ( A -> B ) ***

  if ( itst == 4 ) then
     call sfcflx( &
       &             qao,    qai,    qii,    qio,  swabs,    tsi, &
       &             wev,    wsb, &
       &            prec,   snow,   roff,   soff, &
       &          tauaix, tauaiy, tauaox, tauaoy, &
       &              ft,   ptop,   ssfc, &
       &              tb,     ab,    hib,    tib,    hsb, &
       &              ub,     vb )
#ifdef OPT_BODY
     call bdyflx( &
       &              tq, &
       &              tb )
#endif
     call predci( &
       &              ab,    hib,    uib,    vib,    tib,    hsb, &
       &              ft,     fs,   taux,   tauy,   ptop, &
       &              aa,    hia,    uia,    via,    tia,    hsa, &
       &              tb,     ub,     vb,     hb,     ha, &
       &             qao,    qai,    qii,    qio,  swabs, &
       &             wev,    wsb, &
       &            prec,   snow,   roff,   soff, &
       &          tauaix, tauaiy, tauaox, tauaoy )
     call predco( &
       &              hb,   ubtb,   vbtb,      w,      r, &
       &              ub,     vb,     tb, &
       &             amv,    ahv, &
       &            taux,   tauy,   ptop, &
       &              ft,  swabs,     fs,   ssfc, &
#ifdef OPT_BODY
       &              tq, &
#endif
       &              ua,     va,     ta, &
       &              ha,   ubta,   vbta,     ab, &
       &              gx,     gy,     xx,     yy, &
       &             gxx,    gyy, &
       &            uadv,   vadv,   wadv )

! *** Leap-frog ( B -> A ) ***

  else if ( itst == 3 ) then
     call sfcflx( &
       &             qao,    qai,    qii,    qio,  swabs,    tsi, &
       &             wev,    wsb, &
       &            prec,   snow,   roff,   soff, &
       &          tauaix, tauaiy, tauaox, tauaoy, &
       &              ft,   ptop,   ssfc, &
       &              ta,     aa,    hia,    tia,    hsa, &
       &              ua,     va )
#ifdef OPT_BODY
     call bdyflx( &
       &              tq, &
       &              ta )
#endif
     call predci( &
       &              aa,    hia,    uia,    via,    tia,    hsa, &
       &              ft,     fs,   taux,   tauy,   ptop, &
       &              ab,    hib,    uib,    vib,    tib,    hsb, &
       &              ta,     ua,     va,     ha,     hb, &
       &             qao,    qai,    qii,    qio,  swabs, &
       &             wev,    wsb, &
       &            prec,   snow,   roff,   soff, &
       &          tauaix, tauaiy, tauaox, tauaoy )
     call predco( &
       &              ha,   ubta,   vbta,      w,      r, &
       &              ua,     va,     ta, &
       &             amv,    ahv, &
       &            taux,   tauy,   ptop, &
       &              ft,  swabs,     fs,   ssfc, &
#ifdef OPT_BODY
       &              tq, &
#endif
       &              ub,     vb,     tb, &
       &              hb,   ubtb,   vbtb,     aa, &
       &              gx,     gy,     xx,     yy, &
       &             gxx,    gyy, &
       &            uadv,   vadv,   wadv )

! *** Euler-backward ( B -> B ) ***

  else if ( itst == 1 ) then
     call sfcflx( &
       &             qao,    qai,    qii,    qio,  swabs,    tsi, &
       &             wev,    wsb, &
       &            prec,   snow,   roff,   soff, &
       &          tauaix, tauaiy, tauaox, tauaoy, &
       &              ft,   ptop,   ssfc, &
       &              tb,     ab,    hib,    tib,    hsb, &
       &              ub,     vb )
#ifdef OPT_BODY
     call bdyflx( &
       &              tq, &
       &              tb )
#endif
     call forsti( &
       &              aa,    hia,    uia,    via,    tia,    hsa, &
       &              ab,    hib,    uib,    vib,    tib,    hsb )
     call forsto( &
       &              ua,     va,     ta, &
       &              ha,   ubta,   vbta, &
       &              ub,     vb,     tb, &
       &              hb,   ubtb,   vbtb )
     call predci( &
       &              aa,    hia,    uia,    via,    tia,    hsa, &
       &              ft,     fs,   taux,   tauy,   ptop, &
       &              ab,    hib,    uib,    vib,    tib,    hsb, &
       &              ta,     ua,     va,     ha,     hb, &
       &             qao,    qai,    qii,    qio,  swabs, &
       &             wev,    wsb, &
       &            prec,   snow,   roff,   soff, &
       &          tauaix, tauaiy, tauaox, tauaoy )
     call excngi( &
       &              aa,    hia,    uia,    via,    tia,    hsa, &
       &              ab,    hib,    uib,    vib,    tib,    hsb )
     ieuler = 1
     call predco( &
       &              ha,   ubta,   vbta,      w,      r, &
       &              ua,     va,     ta, &
       &             amv,    ahv, &
       &            taux,   tauy,   ptop, &
       &              ft,  swabs,     fs,   ssfc, &
#ifdef OPT_BODY
       &              tq, &
#endif
       &              ub,     vb,     tb, &
       &              hb,   ubtb,   vbtb,     ab, &
       &              gx,     gy,     xx,     yy, &
       &             gxx,    gyy, &
       &            uadv,   vadv,   wadv )
     ieuler = 2
     call predco( &
       &              hb,   ubtb,   vbtb,      w,      r, &
       &              ub,     vb,     tb, &
       &             amv,    ahv, &
       &            taux,   tauy,   ptop, &
       &              ft,  swabs,     fs,   ssfc, &
#ifdef OPT_BODY
       &              tq, &
#endif
       &              ua,     va,     ta, &
       &              ha,   ubta,   vbta,     ab, &
       &              gx,     gy,     xx,     yy, &
       &             gxx,    gyy, &
       &            uadv,   vadv,   wadv )

! *** Forward ( B -> B ) ***

  else
     ieuler = 0
     call sfcflx( &
       &             qao,    qai,    qii,    qio,  swabs,    tsi, &
       &             wev,    wsb, &
       &            prec,   snow,   roff,   soff, &
       &          tauaix, tauaiy, tauaox, tauaoy, &
       &              ft,   ptop,   ssfc, &
       &              tb,     ab,    hib,    tib,    hsb, &
       &              ub,     vb )
#ifdef OPT_BODY
     call bdyflx( &
       &              tq, &
       &              tb )
#endif
     call forsti( &
       &              aa,    hia,    uia,    via,    tia,    hsa, &
       &              ab,    hib,    uib,    vib,    tib,    hsb )
     call forsto( &
       &              ua,     va,     ta, &
       &              ha,   ubta,   vbta, &
       &              ub,     vb,     tb, &
       &              hb,   ubtb,   vbtb )
     call predci( &
       &              aa,    hia,    uia,    via,    tia,    hsa, &
       &              ft,     fs,   taux,   tauy,   ptop, &
       &              ab,    hib,    uib,    vib,    tib,    hsb, &
       &              ta,     ua,     va,     ha,     hb, &
       &             qao,    qai,    qii,    qio,  swabs, &
       &             wev,    wsb, &
       &            prec,   snow,   roff,   soff, &
       &          tauaix, tauaiy, tauaox, tauaoy )
     call excngi( &
       &              aa,    hia,    uia,    via,    tia,    hsa, &
       &              ab,    hib,    uib,    vib,    tib,    hsb )
     call predco( &
       &              ha,   ubta,   vbta,      w,      r, &
       &              ua,     va,     ta, &
       &             amv,    ahv, &
       &            taux,   tauy,   ptop, &
       &              ft,  swabs,     fs,   ssfc, &
#ifdef OPT_BODY
       &              tq, &
#endif
       &              ub,     vb,     tb, &
       &              hb,   ubtb,   vbtb,     ab, &
       &              gx,     gy,     xx,     yy, &
       &             gxx,    gyy, &
       &            uadv,   vadv,   wadv )
     call excngo( &
       &              ua,     va,     ta, &
       &              ha,   ubta,   vbta, &
       &              ub,     vb,     tb, &
       &              hb,   ubtb,   vbtb )

  end if

! *** Output to file ***
  do ij = 1, nxydim
     aig(ij) = 0.0d0
     hig(ij) = 0.0d0
     hsg(ij) = 0.0d0
  end do

  if ( itst == 3 ) then
     call putsig( &
          &          ta )
     call chekin(    ua,    'U', &
          &      'ocean zonal velocity', 'cm/s', &
          &          nx,     ny,     nz, nxyzdm, 'OCLVTV' )
     call chekin(    va,    'V', &
          & 'ocean meridional velocity', 'cm/s', &
          &          nx,     ny,     nz, nxyzdm, 'OCLVTV' )
     do l = 1, ntdim
        call chekin(     ta(1, l),      ctrnam(l), &
             &          ctrtit(l),      ctruni(l), &
             &         nx,     ny,     nz, nxyzdm, 'OCLVTT' )
     end do
     call chekin(    ha,   'SH', &
          &        'sea surface height',   'cm', &
          &          nx,     ny,      1, nxydim, 'OCSFCT' )
     call chekin(  ubta,  'UBT', &
          &     'ocean zonal transport',         'cm^2/s', &
          &          nx,     ny,      1, nxydim, 'OCSFCV' )
     call chekin(  vbta,  'VBT', &
          &'ocean meridional transport',         'cm^2/s', &
          &          nx,     ny,      1, nxydim, 'OCSFCV' )

     call chekin(    aa,   'AI', &
          &         'ice concentration', 'N.D.', &
          &          nx,     ny,    nic, nxyidm, 'OCICET')
     call chekin(   hia,   'HI', &
          &             'ice thickness', 'cm', &
          &          nx,     ny,    nic, nxyidm, 'OCICET')
     call chekin(   uia,   'UI', &
          &        'ice zonal velocity', 'cm/s', &
          &          nx,     ny,      1, nxydim, 'OCSFCV')
     call chekin(   via,   'VI', &
          &   'ice meridional velocity', 'cm/s', &
          &          nx,     ny,      1, nxydim, 'OCSFCV')
     call chekin(   tia,   'TI', &
          &           'ice temperature', 'degC', &
          &          nx,     ny,    nic, nxyidm, 'OCICET')
     call chekin(   hsa,   'HS', &
          &            'snow thickness',   'cm', &
          &          nx,     ny,    nic, nxyidm, 'OCICET')
     do l = 1, nic
        do ij = 1, nxydim
           aig(ij) = aig(ij) + aa(ij, l)
           hig(ij) = hig(ij) + aa(ij, l) * hia(ij, l)
           hsg(ij) = hsg(ij) + aa(ij, l) * hsa(ij, l)
        end do
     end do
  else
     call putsig( &
          &          tb )
     call chekin(    ub,    'U', &
          &      'ocean zonal velocity', 'cm/s', &
          &          nx,     ny,     nz, nxyzdm, 'OCLVTV')
     call chekin(    vb,    'V', &
          & 'ocean meridional velocity', 'cm/s', &
          &          nx,     ny,     nz, nxyzdm, 'OCLVTV')
     do l = 1, ntdim
        call chekin(     tb(1, l),      ctrnam(l), &
             &          ctrtit(l),      ctruni(l), &
             &                 nx,     ny,     nz, nxyzdm, 'OCLVTT')
     end do
     call chekin(    hb,   'SH', &
          &        'sea surface height', 'cm', &
          &          nx,     ny,      1, nxydim, 'OCSFCT')
     call chekin(  ubtb,  'UBT', &
          &             'ocean zonal transport', 'cm^2/s', &
          &          nx,     ny,      1, nxydim, 'OCSFCV')
     call chekin(  vbtb,  'VBT', &
          &        'ocean meridional transport', 'cm^2/s', &
          &          nx,     ny,      1, nxydim, 'OCSFCV')
     call chekin(    ab,   'AI', &
          &         'ice concentration', 'N.D.', &
          &          nx,     ny,    nic, nxyidm, 'OCICET')
     call chekin(   hib,   'HI', &
          &             'ice thickness',   'cm', &
          &          nx,     ny,    nic, nxyidm, 'OCICET')
     call chekin(   uib,   'UI', &
          &        'ice zonal velocity', 'cm/s', &
          &          nx,     ny,      1, nxydim, 'OCSFCV')
     call chekin(   vib,   'VI', &
          &   'ice meridional velocity', 'cm/s', &
          &          nx,     ny,      1, nxydim, 'OCSFCV')
     call chekin(   tib,   'TI', &
          &           'ice temperature', 'degC', &
          &          nx,     ny,    nic, nxyidm, 'OCICET')
     call chekin(   hsb,   'HS', &
          &            'snow thickness',   'cm', &
          &          nx,     ny,    nic, nxyidm, 'OCICET')
     do l = 1, nic
        do ij = 1, nxydim
           aig(ij) = aig(ij) + ab(ij, l)
           hig(ij) = hig(ij) + ab(ij, l) * hib(ij, l)
           hsg(ij) = hsg(ij) + ab(ij, l) * hsb(ij, l)
        end do
     end do
  end if
  call chekin(     w,    'W', &
       &               'ocean vertical velocity', 'cm/s', &
       &          nx,     ny,     nz, nxyzdm, 'OCLVMT')
  call chekin(   amv,  'AMV', &
       &               'ocean vertical viscosity', 'cm^2/s', &
       &          nx,     ny,     nz, nxyzdm, 'OCLVMV')
  call chekin(   ahv,  'AHV', &
       &               'ocean vertical diffusivity', 'cm^2/s', &
       &          nx,     ny,     nz, nxyzdm, 'OCLVMT')
  do l = 1, ntdim
     call chekin(      ft(1, l),      cftnam(l), &
          &           cfttit(l),      cftuni(l), &
          &          nx,     ny,      1, nxydim, 'OCSFCT')
  end do
  call chekin( swabs,'SWABS', &
       &        'absorbed shortwave',      'erg/cm^2/s', &
       &          nx,     ny,      1, nxydim, 'OCSFCT')
  call chekin(    fs,   'FS', &
       &     'sea surface salt flux',       'psu cm/s', &
       &          nx,     ny,      1, nxydim, 'OCSFCT')
  call chekin(   tsi,  'TSI', &
       &       'sea ice surface temperature',   'degC', &
       &          nx,     ny,    nic, nxyidm, 'OCICET')

! chekin routine for ftx/y/z & fsx/y/z
  call chkftx

! output section for CMIP5
!  AIG: total sea ice concentration, unit [ND]
!  HIG: sea ice thickness averaged over the entire grid cell
!       unit [cm]
!  HSG: mean snow thickness averaged over the entire grid cell
!       unit [cm]
  call chekin(   aig,  'AIG', &
       &          'total sea ice concentration', 'N.D.', &
       &          nx,     ny,      1, nxydim, 'OCSFCT' )
  call chekin(   hig,  'HIG', &
       &    'sea ice thickness averaged over the entire grid cell', &
       &                'cm', &
       &          nx,     ny,      1, nxydim, 'OCSFCT' )
  call chekin(   hsg,  'HSG', &
       &  'mean snow thickness averaged over the entire grid cell', &
       &                'cm', &
       &          nx,     ny,      1, nxydim, 'OCSFCT' )

  call chkout( &
    &          oflout )

  return
end subroutine ocean

end module aocea
