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
!     '09.05.25  Y.Komuro: CMIP5 output code included
!     '09.10.06  Y.Komuro: FORSTO before PREDCI
!     '12.10.12  Y.Komuro: for COCO5.0
!     '13.02.13  Y.Komuro: remove non-parallel code
!     '21.05.26  Y.Komuro: snow aging & meltpond parametrization 
!
! ---------------------------------------------------------------------

  use zocdim, only: &
    &     nx,     ny,     nz, &
    & nxydim, nxyzdm, nxyidm,  nxdim,  ntdim,    nic, &
    &   istr,   iend,   jstr,   jend,   kstr, &
    & ijtstr, ijtend, &
    & myrank, inodes, jnodes, ijnode,  iroot,   ierr, &
    &  oinit, ofinal, mpi_comm_ogcm
  use zocgrd, only: &
    &     dx,     dy,    hxt,    hyt, &
    &     dt, &
    &     tt,     ts,    tss, &
    &     nt,    its,   itst,   ntss, &
    & ieuler
  use zocmsk, only: &
    &  amskt
  use zocfil, only: &
    & nfomax
      
  implicit none

  character(len=16), save :: ctrnam(ntdim), cftnam(ntdim)
  character(len=32), save :: ctrtit(ntdim), cfttit(ntdim)
  character(len=16), save :: ctruni(ntdim), cftuni(ntdim)
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
!$ use omp_lib

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

  integer ::  num_threads

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

!$omp parallel
!$ num_threads = omp_get_num_threads()
!$omp end parallel
!$ write(jfpar, *) "OMP_NUM_THREADS=", num_threads

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
  use qckag
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

!---- arrays for sage.mp (snow aging & melt pond)
  real(8), save ::    asa(nxydim, 0:nic),    asb(nxydim, 0:nic)
  real(8), save ::  frlva(nxydim, 0:nic),  frlvb(nxydim, 0:nic)
  real(8), save ::   vmpa(nxydim, 0:nic),   vmpb(nxydim, 0:nic)
  real(8), save ::  frmpa(nxydim, 0:nic),  frmpb(nxydim, 0:nic)
  real(8), save ::   dsda(nxydim, 0:nic),   dsdb(nxydim, 0:nic)
  real(8), save ::   dsba(nxydim, 0:nic),   dsbb(nxydim, 0:nic)
  real(8), save ::   dfdu(nxydim),   dfbc(nxydim)

#ifdef OPT_BODY
  real(8), save ::      tq(nxyzdm, ntdim)
#endif

  real(8) ::    ssfc(nxydim)
  real(8) ::     aig(nxydim),    hig(nxydim),    hsg(nxydim)
  real(8) ::     asg(nxydim),  frlvg(nxydim),   vmpg(nxydim)
  real(8) ::   frmpg(nxydim),   dsdg(nxydim),   dsbg(nxydim)

  integer ::     ij,      l

  if (oinit) then
     call predci ( &
       &               ab,    hib,    uib,    vib,    tib,    hsb, &
       &              asb,  frlvb,   vmpb,  frmpb,   dsdb,   dsbb, &
       &               ft,     fs,   taux,   tauy,   ptop, &
       &               aa,    hia,    uia,    via,    tia,    hsa, &
       &              asa,  frlva,   vmpa,  frmpa,   dsda,   dsba, &
       &               tb,     ub,     vb,     hb,     ha, &
       &              qao,    qai,    qii,    qio,  swabs,    tsi, &
       &              wev,    wsb, &
       &             prec,   snow,   roff,   soff, &
       &             dfdu,   dfbc, &
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
       &              asb,  frlvb,   vmpb,  frmpb,   dsdb,   dsbb, &
       &               ft,     fs,   taux,   tauy,   ptop, &
       &               aa,    hia,    uia,    via,    tia,    hsa, &
       &              asa,  frlva,   vmpa,  frmpa,   dsda,   dsba, &
       &               tb,     ub,     vb,     hb,     ha, &
       &              qao,    qai,    qii,    qio,  swabs,    tsi, &
       &              wev,    wsb, &
       &             prec,   snow,   roff,   soff, &
       &             dfdu,   dfbc, &
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
       &            dfdu,   dfbc, &
       &              tb,     ab,    hib,    tib,    hsb, &
       &             asb,   vmpb,  frmpb, &
       &              ub,     vb )
     call nmlper( &
       &             wev,   prec, &
       &             wsb,   snow,   roff,   soff)
#ifdef OPT_BODY
     call bdyflx( &
       &              tq, &
       &              tb )
#endif
     call predci( &
       &              ab,    hib,    uib,    vib,    tib,    hsb, &
       &             asb,  frlvb,   vmpb,  frmpb,   dsdb,   dsbb, &
       &              ft,     fs,   taux,   tauy,   ptop, &
       &              aa,    hia,    uia,    via,    tia,    hsa, &
       &             asa,  frlva,   vmpa,  frmpa,   dsda,   dsba, &
       &              tb,     ub,     vb,     hb,     ha, &
       &             qao,    qai,    qii,    qio,  swabs,    tsi, &
       &             wev,    wsb, &
       &            prec,   snow,   roff,   soff, &
       &            dfdu,   dfbc, &
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
       &            dfdu,   dfbc, &
       &              ta,     aa,    hia,    tia,    hsa, &
       &             asa,   vmpa,  frmpa, &
       &              ua,     va )
     call nmlper( &
       &             wev,   prec, &
       &             wsb,   snow,   roff,   soff)
#ifdef OPT_BODY
     call bdyflx( &
       &              tq, &
       &              ta )
#endif
     call predci( &
       &              aa,    hia,    uia,    via,    tia,    hsa, &
       &             asa,  frlva,   vmpa,  frmpa,   dsda,   dsba, &
       &              ft,     fs,   taux,   tauy,   ptop, &
       &              ab,    hib,    uib,    vib,    tib,    hsb, &
       &             asb,  frlvb,   vmpb,  frmpb,   dsdb,   dsbb, &
       &              ta,     ua,     va,     ha,     hb, &
       &             qao,    qai,    qii,    qio,  swabs,    tsi, &
       &             wev,    wsb, &
       &            prec,   snow,   roff,   soff, &
       &            dfdu,   dfbc, &
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
       &            dfdu,   dfbc, &
       &              tb,     ab,    hib,    tib,    hsb, &
       &             asb,   vmpb,  frmpb, &
       &              ub,     vb )
     call nmlper( &
       &             wev,   prec, &
       &             wsb,   snow,   roff,   soff)
#ifdef OPT_BODY
     call bdyflx( &
       &              tq, &
       &              tb )
#endif
     call forsti( &
       &              aa,    hia,    uia,    via,    tia,    hsa, &
       &              ab,    hib,    uib,    vib,    tib,    hsb )
     call forsta( &
       &             asa,  frlva,   vmpa,  frmpa,   dsda,   dsba, &
       &             asb,  frlvb,   vmpb,  frmpb,   dsdb,   dsbb )
     call forsto( &
       &              ua,     va,     ta, &
       &              ha,   ubta,   vbta, &
       &              ub,     vb,     tb, &
       &              hb,   ubtb,   vbtb )
     call predci( &
       &              aa,    hia,    uia,    via,    tia,    hsa, &
       &             asa,  frlva,   vmpa,  frmpa,   dsda,   dsba, &
       &              ft,     fs,   taux,   tauy,   ptop, &
       &              ab,    hib,    uib,    vib,    tib,    hsb, &
       &             asb,  frlvb,   vmpb,  frmpb,   dsdb,   dsbb, &
       &              ta,     ua,     va,     ha,     hb, &
       &             qao,    qai,    qii,    qio,  swabs,    tsi, &
       &             wev,    wsb, &
       &            prec,   snow,   roff,   soff, &
       &            dfdu,   dfbc, &
       &          tauaix, tauaiy, tauaox, tauaoy )
     call excngi( &
       &              aa,    hia,    uia,    via,    tia,    hsa, &
       &              ab,    hib,    uib,    vib,    tib,    hsb )
     call excnga( &
       &             asa,  frlva,   vmpa,  frmpa,   dsda,   dsba, &
       &             asb,  frlvb,   vmpb,  frmpb,   dsdb,   dsbb )
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
       &            dfdu,   dfbc, &
       &              tb,     ab,    hib,    tib,    hsb, &
       &             asb,   vmpb,  frmpb, &
       &              ub,     vb )
     call nmlper( &
       &             wev,   prec, &
       &             wsb,   snow,   roff,   soff)
#ifdef OPT_BODY
     call bdyflx( &
       &              tq, &
       &              tb )
#endif
     call forsti( &
       &              aa,    hia,    uia,    via,    tia,    hsa, &
       &              ab,    hib,    uib,    vib,    tib,    hsb )
     call forsta( &
       &             asa,  frlva,   vmpa,  frmpa,   dsda,   dsba, &
       &             asb,  frlvb,   vmpb,  frmpb,   dsdb,   dsbb )
     call forsto( &
       &              ua,     va,     ta, &
       &              ha,   ubta,   vbta, &
       &              ub,     vb,     tb, &
       &              hb,   ubtb,   vbtb )
     call predci( &
       &              aa,    hia,    uia,    via,    tia,    hsa, &
       &             asa,  frlva,   vmpa,  frmpa,   dsda,   dsba, &
       &              ft,     fs,   taux,   tauy,   ptop, &
       &              ab,    hib,    uib,    vib,    tib,    hsb, &
       &             asb,  frlvb,   vmpb,  frmpb,   dsdb,   dsbb, &
       &              ta,     ua,     va,     ha,     hb, &
       &             qao,    qai,    qii,    qio,  swabs,    tsi, &
       &             wev,    wsb, &
       &            prec,   snow,   roff,   soff, &
       &            dfdu,   dfbc, &
       &          tauaix, tauaiy, tauaox, tauaoy )
     call excngi( &
       &              aa,    hia,    uia,    via,    tia,    hsa, &
       &              ab,    hib,    uib,    vib,    tib,    hsb )
     call excnga( &
       &             asa,  frlva,   vmpa,  frmpa,   dsda,   dsba, &
       &             asb,  frlvb,   vmpb,  frmpb,   dsdb,   dsbb )
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

  call cofptu( &
    &               taux,   tauy)

! *** Output to file ***
  do ij = 1, nxydim
     aig(ij) = 0.0d0
     hig(ij) = 0.0d0
     hsg(ij) = 0.0d0
     asg(ij) = 0.0d0
     frlvg(ij) = 0.0d0
     vmpg(ij) = 0.0d0
     frmpg(ij) = 0.0d0
     dsdg(ij) = 0.0d0
     dsbg(ij) = 0.0d0
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
     call chekin(   asa,   'AS', &
          &                    'snow age', 'ND', &
          &          nx,     ny,    nic, nxyidm, 'OCICET')
     call chekin( frlva, 'FRLV', &
          &          'level ice fraction', 'ND', &
          &          nx,     ny,    nic, nxyidm, 'OCICET')
     call chekin(  vmpa,  'VMP', &
          &             'melt pond depth', 'cm', &
          &          nx,     ny,    nic, nxyidm, 'OCICET')
     call chekin( frmpa, 'FRMP', &
          &          'melt pond fraction', 'ND', &
          &          nx,     ny,    nic, nxyidm, 'OCICET')
     call chekin(  dsda, 'DSDU', &
          & 'concentration of dust, non-bc', 'g/cm^2', &
          &          nx,     ny,    nic, nxyidm, 'OCICET')
     call chekin(  dsba, 'DSBC', &
          & 'concentration of dust, bc', 'g/cm^2', &
          &          nx,     ny,    nic, nxyidm, 'OCICET')
     do l = 1, nic
        do ij = 1, nxydim
           aig(ij) = aig(ij) + aa(ij, l)
           hig(ij) = hig(ij) + aa(ij, l) * hia(ij, l)
           hsg(ij) = hsg(ij) + aa(ij, l) * hsa(ij, l)
           asg(ij) = asg(ij) + aa(ij, l) * asa(ij, l)
           frlvg(ij) = frlvg(ij) + aa(ij, l) * frlva(ij, l)
           vmpg(ij) = vmpg(ij) + aa(ij, l) * vmpa(ij, l)
           frmpg(ij) = frmpg(ij) + aa(ij, l) * frmpa(ij, l)
           dsdg(ij) = dsdg(ij) + aa(ij, l) * dsda(ij, l)
           dsbg(ij) = dsbg(ij) + aa(ij, l) * dsba(ij, l)
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
     call chekin(   asb,   'AS', &
          &                    'snow age', 'ND', &
          &          nx,     ny,    nic, nxyidm, 'OCICET')
     call chekin( frlvb, 'FRLV', &
          &          'level ice fraction', 'ND', &
          &          nx,     ny,    nic, nxyidm, 'OCICET')
     call chekin(  vmpb,  'VMP', &
          &             'melt pond depth', 'cm', &
          &          nx,     ny,    nic, nxyidm, 'OCICET')
     call chekin( frmpb, 'FRMP', &
          &          'melt pond fraction', 'ND', &
          &          nx,     ny,    nic, nxyidm, 'OCICET')
     call chekin(  dsdb, 'DSDU', &
          & 'concentration of dust, non-bc', 'g/cm^2', &
          &          nx,     ny,    nic, nxyidm, 'OCICET')
     call chekin(  dsbb, 'DSBC', &
          & 'concentration of dust, bc', 'g/cm^2', &
          &          nx,     ny,    nic, nxyidm, 'OCICET')
     do l = 1, nic
        do ij = 1, nxydim
           aig(ij) = aig(ij) + ab(ij, l)
           hig(ij) = hig(ij) + ab(ij, l) * hib(ij, l)
           hsg(ij) = hsg(ij) + ab(ij, l) * hsb(ij, l)
           asg(ij) = asg(ij) + ab(ij, l) * asb(ij, l)
           frlvg(ij) = frlvg(ij) + ab(ij, l) * frlvb(ij, l)
           vmpg(ij) = vmpg(ij) + ab(ij, l) * vmpb(ij, l)
           frmpg(ij) = frmpg(ij) + ab(ij, l) * frmpb(ij, l)
           dsdg(ij) = dsdg(ij) + ab(ij, l) * dsdb(ij, l)
           dsbg(ij) = dsbg(ij) + ab(ij, l) * dsbb(ij, l)
        end do
     end do
  end if
  do ij = 1, nxydim
     if (aig(ij) > 0.0d0) then
        frlvg(ij) = frlvg(ij) / aig(ij)
        frmpg(ij) = frmpg(ij) / aig(ij)
     end if
  end do
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
  call chekin(   asg,   'ASG', &
       &                    'snow age', 'ND', &
       &          nx,     ny,      1, nxydim, 'OCSFCT')
  call chekin( frlvg, 'FRLVG', &
       &          'level ice fraction', 'ND', &
       &          nx,     ny,      1, nxydim, 'OCSFCT')
  call chekin(  vmpg,  'VMPG', &
       &             'melt pond depth', 'cm', &
       &          nx,     ny,      1, nxydim, 'OCSFCT')
  call chekin( frmpg, 'FRMPG', &
       &          'melt pond fraction', 'ND', &
       &          nx,     ny,      1, nxydim, 'OCSFCT')
  call chekin(  dsdg, 'DSDUG', &
       & 'concentration of dust, non-bc', 'g/cm^2', &
       &          nx,     ny,      1, nxydim, 'OCSFCT')
  call chekin(  dsbg, 'DSBCG', &
       & 'concentration of dust, bc', 'g/cm^2', &
       &          nx,     ny,      1, nxydim, 'OCSFCT')

! surface flux output
  call chksfx

  call chkout( &
    &          oflout )

  return
end subroutine ocean

! *********************************************************************

subroutine nmlper( &
  &                  wev,   prec, &
  &                  wsb,   snow,   roff,   soff)


! --- information -----------------------------------------------------
!
!  normalize P-E+R
!
!  HISTORY
!     '21.04.07  Y.Komuro: for coco5.0
!
! ---------------------------------------------------------------------

  use qckag
  use qckot
  use ufile

#include "mpif.h"

  real(8), intent(inout) ::   prec(nxydim),    wev(nxydim)
  real(8), intent(in)    ::   roff(nxydim),   soff(nxydim)
  real(8), intent(in)    ::   snow(nxydim),    wsb(nxydim, nic)

  real(8), save ::  garea(nxydim)
  real(8), save ::  tarea, rtardt
  
  real(8) ::  fwnmd(nxydim)
  real(8) :: vwteqg(inodes*jnodes)
  real(8) :: vwtreq, vwteqt, fwnml

  integer ::    ij,      i,      j,      l
  integer :: ifpar,  jfpar,  istat
  logical, save :: ofirst = .true.

  logical, save :: onmper = .false.
  namelist /nmnper/ onmper

  if (ofirst) then
     call rewnml(ifpar, jfpar)
     read(ifpar, nmnper, iostat=istat)
     call cstnml(jfpar, 'nmlper', 'nmnper', istat)
     write(jfpar, nmnper)
     do ij = 1, nxydim
        garea(ij) = 0.0d0
        fwnmd(ij) = 0.0d0
     end do
     ofirst = .false.

     if (onmper) then
        write(jfpar, *) &
          &    '*** Normalization of P-E+R is applied. ***'
        vwteqt = 0.0d0
        tarea = 0.0d0
        do i = 1, inodes*jnodes
           vwteqg(i) = 0.0d0
        end do
        do j=jstr, jend
           do i=istr, iend
              ij = nxdim*(j-1) + i
              garea(ij) = dx * dy(ij) * hxt(ij) * hyt(ij) &
 &                        * amskt(ij, kstr)
              vwteqt = vwteqt + garea(ij)
           end do
        end do
!        write(jfpar, *) vwteqt
        call mpi_gather( &
          &  vwteqt, 1, mpi_real8, vwteqg(1), 1, mpi_real8, &
          &  iroot, mpi_comm_ogcm, ierr)
        if (myrank == iroot) then
           do i = 1, inodes*jnodes
              tarea = tarea + vwteqg(i)
           end do
        end if
        call mpi_bcast( &
          &  tarea, 1, mpi_real8, &
          &  iroot, mpi_comm_ogcm, ierr)
!        call mpi_allreduce( &
!          &   vwteqt, tarea, 1, mpi_real8, &
!          &   mpi_sum, mpi_comm_ogcm, ierr)
!        write(jfpar, *) tarea
        rtardt = 1.0d0 / tarea
     endif
  end if

  if (.not.onmper) then
     return
  endif

  vwteqt = 0.0d0
  fwnml = 0.0d0
  do i = 1, inodes*jnodes
     vwteqg(i) = 0.0d0
  end do
  do ij = ijtstr, ijtend
     vwtreq = prec(ij) + snow(ij) + roff(ij) + soff(ij) - wev(ij)
     do l = 1, nic
        vwtreq = vwtreq - wsb(ij, l)
     end do
     vwteqt = vwteqt + vwtreq * garea(ij)
  end do
!  write(jfpar, *) vwteqt
  call mpi_gather( &
    &  vwteqt, 1, mpi_real8, vwteqg(1), 1, mpi_real8, &
    &  iroot, mpi_comm_ogcm, ierr)
  if (myrank == iroot) then
     do i = 1, inodes*jnodes
        fwnml = fwnml + vwteqg(i)
     end do
  end if
  call mpi_bcast( &
    &  fwnml, 1, mpi_real8, &
    &  iroot, mpi_comm_ogcm, ierr)
!  call mpi_allreduce( &
!    &  vwteqt, fwnml, 1, mpi_real8, &
!    &  mpi_sum, mpi_comm_ogcm, ierr)
!  write(jfpar, *) fwnml
  fwnml = fwnml * rtardt

  do ij = ijtstr, ijtend
     prec(ij) = prec(ij) - min(fwnml, 0.0d0) * amskt(ij, kstr)
     wev(ij) = wev(ij) + max(fwnml, 0.0d0) * amskt(ij, kstr)
     fwnmd(ij) = fwnml
  end do

  call cofpnw( &
    &           fwnmd )
    
  call chekin( fwnmd, 'FWNML', &
    &         'Fw for normalizing surface height', 'cm/s', &
    &             nx,     ny,      1, nxydim, 'OCSFCT')      

  return
end subroutine nmlper

end module aocea
