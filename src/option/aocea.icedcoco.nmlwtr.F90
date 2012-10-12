module aocea

! --- information -----------------------------------------------------
!
!  HISTORY
!     '02.10.10  H.Hasumi: from COCO3.4
!     '07.04.23  H.Hasumi
!     '07.09.25  H.Hasumi: multi-category sea ice
!     '07.09.27  H.Hasumi: 1-layer sea ice thermodynamics
!     '08.06.11  H.Hasumi: initial/final processing
!     '08.07.10  H.Hasumi: initial/final processing
!     '08.08.28  Y.Komuro: normalization of water volume
!     '09.05.25  Y.Komuro: CMIP5 output code included
!     '09.10.06  Y.Komuro: FORSTO before PREDCI
!     '11.12.13  Y.Komuro: time scale for water normalization
!     '12.10.12  Y.Komuro: for COCO5.0
!
! ---------------------------------------------------------------------

  use zocdim, only: &
    &     nx,     ny,     nz, &
    &   istr,   iend,   jstr,   jend,   kstr, &
    &  nxdim, nxydim, nxyzdm, nxyidm,  ntdim,    nic, &
    & ijtstr, ijtend, &
#ifdef OPT_PARALLEL
    & myrank, ijnode,   ierr, &
#endif
    &  oinit, ofinal
  use zocgrd, only: &
    &     dx,     dy, &
    &    hxt,    hyt, &
    &     dt, &
    &     tt,     ts,    tss, &
    &     nt,    its,   itst,   ntss, &
    & ieuler
  use zocmsk, only: &
    &  amskt
  use zocfil, only: &
    & nfomax
  use zocphy, only: &
    &   rhoi,   rhos
      
  implicit none

  character(len=16), save :: ctrnam(ntdim), cftnam(ntdim)
  real(8), save ::    uadv(nxyzdm),   vadv(nxyzdm),   wadv(nxyzdm)
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
  real(8), intent(out) ::     ft(nxydim, ntdim)
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

  DT = DT1

  call rdgeo
  call svtset
  call ovtset(r, tb)
  call chkset

! *** Initialization for variables ***
#ifdef OPT_PARALLEL
  if (myrank < ijnode) then
#endif
     call iniset( &
       &             uadv,   vadv,   wadv,      r, &
       &               ub,     vb,     tb, &
       &               hb,   ubtb,   vbtb, &
       &                w,    amv,    ahv )
#ifdef OPT_PARALLEL
  end if
#endif

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
     return
  end if

  if (ofinal) then
#ifdef OPT_PARALLEL
     if (myrank < ijnode) then
#endif
        call iniset( &
          &            uadv,   vadv,   wadv,      r, &
          &              ub,     vb,     tb, &
          &              hb,   ubtb,   vbtb, &
          &               w,    amv,    ahv )
#ifdef OPT_PARALLEL
     end if
#endif
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
       &              tb,     ab,    hib,    tib,    hsb )
     call nmlwtr( &
       &            prec,    wev, &
       &              hb,     ab,    hib,    hsb )
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
       &              ta,     aa,    hia,    tia,    hsa )
     call nmlwtr( &
       &            prec,    wev, &
       &              ha,     aa,    hia,    hsa )
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
       &              tb,     ab,    hib,    tib,    hsb )
     call nmlwtr( &
       &            prec,    wev, &
       &              hb,     ab,    hib,    hsb )
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
       &              tb,     ab,    hib,    tib,    hsb )
     call nmlwtr( &
       &            prec,    wev, &
       &              hb,     ab,    hib,    hsb )
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
     call chekin(    ua,    'U', &
       &             nx,     ny,     nz, nxyzdm, 'OCN' )
     call chekin(    va,    'V', &
       &             nx,     ny,     nz, nxyzdm, 'OCN' )
     do l = 1, ntdim
        call chekin(     ta(1, l),      ctrnam(l), &
          &              nx,     ny,     nz, nxyzdm, 'OCN' )
     end do
     call chekin(    ha,   'SH', &
       &             nx,     ny,      1, nxydim, 'SFC' )
     call chekin(  ubta,  'UBT', &
       &             nx,     ny,      1, nxydim, 'SFC' )
     call chekin(  vbta,  'VBT', &
       &             nx,     ny,      1, nxydim, 'SFC' )
     call chekin(    aa,   'AI', &
       &             nx,     ny,    nic, nxyidm, 'ICE' )
     call chekin(   hia,   'HI', &
       &             nx,     ny,    nic, nxyidm, 'ICE' )
     call chekin(   uia,   'UI', &
       &             nx,     ny,      1, nxydim, 'SFC' )
     call chekin(   via,   'VI', &
       &             nx,     ny,      1, nxydim, 'SFC' )
     call chekin(   tia,   'TI', &
       &             nx,     ny,    nic, nxyidm, 'ICE' )
     call chekin(   hsa,   'HS', &
       &             nx,     ny,    nic, nxyidm, 'ICE' )
     do l = 1, nic
        do ij = 1, nxydim
           aig(ij) = aig(ij) + aa(ij, l)
           hig(ij) = hig(ij) + aa(ij, l) * hia(ij, l)
           hsg(ij) = hsg(ij) + aa(ij, l) * hsa(ij, l)
        end do
     end do
  else
     call chekin(    ub,    'U', &
       &             nx,     ny,     nz, nxyzdm, 'OCN' )
     call chekin(    vb,    'V', &
       &             nx,     ny,     nz, nxyzdm, 'OCN' )
     do l = 1, ntdim
        call chekin(     tb(1, l),      ctrnam(l), &
          &              nx,     ny,     nz, nxyzdm, 'OCN' )
     end do
     call chekin(    hb,   'SH', &
       &             nx,     ny,      1, nxydim, 'SFC' )
     call chekin(  ubtb,  'UBT', &
       &             nx,     ny,      1, nxydim, 'SFC' )
     call chekin(  vbtb,  'VBT', &
       &             nx,     ny,      1, nxydim, 'SFC' )
     call chekin(    ab,   'AI', &
       &             nx,     ny,    nic, nxyidm, 'ICE' )
     call chekin(   hib,   'HI', &
       &             nx,     ny,    nic, nxyidm, 'ICE' )
     call chekin(   uib,   'UI', &
       &             nx,     ny,      1, nxydim, 'SFC' )
     call chekin(   vib,   'VI', &
       &             nx,     ny,      1, nxydim, 'SFC' )
     call chekin(   tib,   'TI', &
       &             nx,     ny,    nic, nxyidm, 'ICE' )
     call chekin(   hsb,   'HS', &
       &             nx,     ny,    nic, nxyidm, 'ICE' )
     do l = 1, nic
        do ij = 1, nxydim
           aig(ij) = aig(ij) + ab(ij, l)
           hig(ij) = hig(ij) + ab(ij, l) * hib(ij, l)
           hsg(ij) = hsg(ij) + ab(ij, l) * hsb(ij, l)
        end do
     end do
  end if
  call chekin(     w,    'W', &
    &             nx,     ny,     nz, nxyzdm, 'OCN' )
  call chekin(   amv,  'AMV', &
    &             nx,     ny,     nz, nxyzdm, 'OCN' )
  call chekin(   ahv,  'AHV', &
    &             nx,     ny,     nz, nxyzdm, 'OCN' )
  do l = 1, ntdim
     call chekin(      ft(1, l),      cftnam(l), &
       &               nx,     ny,      1, nxydim, 'SFC' )
  end do
  call chekin( swabs,'SWABS', &
    &             nx,     ny,      1, nxydim, 'SFC' )
  call chekin(    fs,    'FS', &
    &             nx,     ny,      1, nxydim, 'SFC' )
  call chekin(   tsi,   'TSI', &
    &             nx,     ny,    nic, nxyidm, 'ICE' )

! chekin routine for ftx/y/z & fsx/y/z
  call chkftx

! output section for CMIP5
!  AIG: total sea ice concentration, unit [ND]
!  HIG: sea ice thickness averaged over the entire grid cell
!       unit [cm]
!  HSG: mean snow thickness averaged over the entire grid cell
!       unit [cm]
  call chekin(   aig,  'AIG', &
    &             nx,     ny,      1, nxydim, 'SFC' )
  call chekin(   hig,  'HIG', &
    &             nx,     ny,      1, nxydim, 'SFC' )
  call chekin(   hsg,  'HSG', &
    &             nx,     ny,      1, nxydim, 'SFC' )

  call chkout( &
    &          oflout )

  return
end subroutine ocean

! *********************************************************************

subroutine nmlwtr( &
  &                  prec,    wev, &
  &                    sh,      a,     hi,     hs )

! --- information -----------------------------------------------------
!
!  normalize water volume 
!
!  HISTORY
!     '08.08.28  Y.Komuro: for COCO4.4
!     '12.10.12  Y.Komuro: for COCO5.0
!
! ---------------------------------------------------------------------

  use qckot
  use ufile

#ifdef OPT_PARALLEL
#include "mpif.h"
#endif

  real(8), intent(inout) ::   prec(nxydim),    wev(nxydim)
  real(8), intent(in)    ::     sh(nxydim)
  real(8), intent(in)    ::      a(nxydim, 0:nic)
  real(8), intent(in)    ::     hi(nxydim, 0:nic)
  real(8), intent(in)    ::     hs(nxydim, 0:nic)

  real(8), save ::   garea(nxydim)
  real(8), save ::  rtardt

  logical, save :: ofirst = .true.

  real(8) ::  fwnmd(nxydim)
  real(8) :: vwtreq, vwteqt,  fwnml,  tarea

  integer ::     ij,      i,      j,      l
  integer ::  ifpar,  jfpar,  istat

  logical, save :: onmlwt = .false.
  real(8), save :: wdmp = 0.0d0

  namelist /nmnmlw/ onmlwt, wdmp

  if (ofirst) then
     call rewnml(ifpar, jfpar)
     read (ifpar, nmnmlw, iostat=istat)
     call cstnml(jfpar, 'nmlwtr', 'nmnmlw', istat)
     write(jfpar, nmnmlw)
     do ij = 1, nxydim
        garea(ij) = 0.0d0
        fwnmd(ij) = 0.0d0
     end do
     ofirst = .false.

     if (onmlwt) then
        write(jfpar, *) &
          &   '*** Normalization of water volume is applied. ***'
        vwteqt = 0.0d0
        tarea = 0.0d0
        do j=jstr, jend
           do i=istr, iend
              ij = nxdim*(j-1) + i
              garea(ij) = dx * dy(ij) * hxt(ij) * hyt(ij) &
                &       * amskt(ij, kstr)
              vwteqt = vwteqt + garea(ij)
           end do
        end do
#ifdef OPT_PARALLEL
        call mpi_allreduce( &
          &                vwteqt, tarea, 1, mpi_real8, &
          &               mpi_sum, mpi_comm_world, ierr)
#else
        tarea = vwteqt
#endif
        rtardt = 1.0d0 / tarea / (8.64d4 * wdmp)
     endif
  end if

  if (.not.onmlwt) then
     return
  endif

  vwteqt = 0.0d0
  do ij = ijtstr, ijtend
     vwtreq = sh(ij)
     do l = 1, nic
        vwtreq = vwtreq + a(ij, l) * &
          &      ( rhoi * hi(ij, l) + rhos * hs(ij, l) )
     end do
     vwteqt = vwteqt + vwtreq * garea(ij)
  end do
#ifdef OPT_PARALLEL
  call mpi_allreduce( &
    &                vwteqt, fwnml, 1, mpi_real8, &
    &               mpi_sum, mpi_comm_world, ierr)
#else
  fwnml = vwteqt
#endif
  fwnml = fwnml * rtardt

  do ij = ijtstr, ijtend
     prec(ij) = prec(ij) - min(fwnml, 0.0d0) * amskt(ij, kstr)
     wev(ij) = wev(ij) + max(fwnml, 0.0d0) * amskt(ij, kstr)
     fwnmd(ij) = fwnml
  end do

  call chekin( fwnmd, 'FWNML', &
    &             nx,     ny,      1, nxydim, 'SFC' )
      
  return
end subroutine nmlwtr

end module aocea
