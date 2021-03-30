module aocea

! --- information -----------------------------------------------------
!
!  HISTORY
!  HISTORY
!     '00.03.17  H.Hasumi: for Ishizaki scheme
!     '00.05.30  H.Hasumi: parallelized COCO3
!     '00.12.07  H.Hasumi: combine par and non-par routines
!     '01.01.23  H.Hasumi: for hybrid vertical coordinate
!                          zonal filter
!     '01.02.02  H.Hasumi: divide the tracer routine
!     '01.02.04  H.Hasumi: combine body routines
!     '01.02.20  H.Hasumi: for incorporating BBL model
!     '01.12.03  H.Hasumi: subroutine FVSET is split into two
!     '01.12.07  H.Hasumi
!     '02.05.29  H.Nakano: tracer dimension
!     '07.04.23  H.Hasumi
!     '07.09.25  H.Hasumi: arguments of CHEKIN
!     '08.06.11  H.Hasumi: initial/final processing
!     '08.07.10  H.Hasumi: initial/final processing
!     '10.04.14  M.Kurogi: (COCO4.4 tripolar code by Dr. Suzuki)
!     '13.01.07  T.Suzuki: for COCO5.0 in F90
! ---------------------------------------------------------------------

  use zocdim,  only  :  nxydim,  ntdim,  nxyzdm
  use zocfil,  only  :  nfomax

  implicit none

  private

  real(8),            save  ::  uadv(1:nxyzdm),  vadv(1:nxyzdm)
  real(8),            save  ::  wadv(1:nxyzdm, 9)
  real(8),            save  ::  taux(1:nxydim),  tauy(1:nxydim) 
  real(8),            save  ::  ptop(1:nxydim) 
  real(8),            save  ::   gxx(1:nxydim),   gyy(1:nxydim)
  real(8),            save  ::    xx(1:nxyzdm),    yy(1:nxyzdm)
  real(8),            save  ::    gx(1:nxyzdm),    gy(1:nxyzdm)
#ifdef OPT_BODY
  real(8),            save  ::    tq(1:nxyzdm,1:ntdim)
#endif

  character(len=16), save :: ctrnam(ntdim), cftnam(ntdim)
  character(len=32), save :: ctrtit(ntdim), cfttit(ntdim)
  character(len=16), save :: ctruni(ntdim), cftuni(ntdim)

  public  ::  ocstup
  public  ::  ocean

contains


  subroutine ocstup(     ub,     vb,     tb,                          &
      &                  hb,   ubtb,   vbtb,      w,      r,          &
      &                 amv,    ahv,                                  &
      &                  ft,                                          &
      &                 dt1 )

    use zocdim,  only  :  nxyzdm,  nxydim,  ntdim
    use zocnod,  only  :  myrank,  ijnode
    use zocgrd
    use zocfil,  only  :  nfomax
    use tslvt
    use tflxt
    use tovtr
    use brdge
    use binst
    use sfcng
    use ufile
    use qckot
    use bfrch
    use aprdc
    
    implicit none

    real(8),        intent(inout)  ::    ub(1:nxyzdm),  vb(1:nxyzdm)  
    real(8),        intent(inout)  ::    tb(1:nxyzdm,1:ntdim)
    real(8),        intent(inout)  ::    hb(1:nxydim)
    real(8),        intent(inout)  ::  ubtb(1:nxydim),  vbtb(1:nxydim)
    real(8),        intent(inout)  ::     w(1:nxyzdm),     r(1:nxyzdm)
    real(8),        intent(inout)  ::   amv(1:nxyzdm),   ahv(1:nxyzdm)
    real(8),        intent(inout)  ::    ft(1:nxydim,1:ntdim)
    real(8),        intent(in)     ::   dt1

    integer(4)   ::   ifpar,  jfpar
    integer(4)   ::      ij,      l

    call rewnml( ifpar, jfpar )
    write(jfpar, *) '*** OCSTUP ***'
    
    ctrnam(1) = 'T'
    ctrnam(2) = 'S'
    cftnam(1) = 'FT'
    cftnam(2) = 'FW'
    do l = 3, ntdim
       write(ctrnam(l), '(a6,i2.2)') 'TRACER', L
       write(cftnam(l), '(a6,i2.2)') 'TRCFLX', L
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

    dt = dt1

    call rdgeo
    call svtset
    call ovtset( r, tb )
    call chkset

! *** Initialization for variables ***
    if ( myrank < ijnode ) then
       call iniset(                                                 &
    &             uadv,   vadv,   wadv,      r,                     &
    &               ub,     vb,     tb,                             &
    &               hb,   ubtb,   vbtb,                             &
    &                w,    amv,    ahv,                             &
    &             taux,   tauy,                                     &
    &               ft  )
    end if

    do ij = 1, nxydim
       gxx(ij)  = 0.d0
       gyy(ij)  = 0.d0
       ptop(ij) = 0.d0
    end do

  end subroutine ocstup
    
! =====================================================================

  subroutine ocean(                                                   &
     &                ua,     va,     ta,                             &
     &                ha,   ubta,   vbta,                             &
     &                ub,     vb,     tb,                             &
     &                hb,   ubtb,   vbtb,                             &
     &                 w,      r,    amv,    ahv,                     &
     &                ft,                                             &
     &               nt1,    tt1,  itst1,    ts1,   its1,             &
     &             ntss1,   tss1,                                     &
     &            oflout, oflstk  )

    use zocdim,  only  :                                              &
     &       nxdim,   nydim,  nzdim,  ntdim,                          & 
     &          nx,      ny,     nz,                                  &
     &      nxyzdm,  nxydim,                                          &
     &       oinit,  ofinal
    use zocnod,  only  :  myrank,  ijnode
    use zocgrd
    use zocfil,  only  :  nfomax
    use tslvt
    use tflxt
    use tovtr
    use dvdif
    use brdge
    use binst
    use sfcng
    use ufile
    use qckag
    use qckot
    use bfrch
    use aprdc
    use bshft

    implicit none

    real(8),        intent(inout)  ::    ua(1:nxyzdm),    va(1:nxyzdm)  
    real(8),        intent(inout)  ::    ta(1:nxyzdm,1:ntdim)
    real(8),        intent(inout)  ::    ha(1:nxydim)
    real(8),        intent(inout)  ::  ubta(1:nxydim),  vbta(1:nxydim)
    real(8),        intent(inout)  ::    ub(1:nxyzdm),    vb(1:nxyzdm)  
    real(8),        intent(inout)  ::    tb(1:nxyzdm,1:ntdim)
    real(8),        intent(inout)  ::    hb(1:nxydim)
    real(8),        intent(inout)  ::  ubtb(1:nxydim),  vbtb(1:nxydim)
    real(8),        intent(inout)  ::     w(1:nxyzdm),     r(1:nxyzdm)
    real(8),        intent(inout)  ::   amv(1:nxyzdm),   ahv(1:nxyzdm)
    real(8),        intent(inout)  ::    ft(1:nxydim,1:ntdim)
    integer(4),     intent(in)     ::   nt1,    its1,    itst1,   ntss1
    real(8),        intent(in)     ::   tt1,     ts1,     tss1
    logical                        ::  oflout(1:nfomax)    
    logical                        ::  oflstk(1:nfomax)          

    integer(4)  ::  l

    
    if ( oinit ) then
       call predco(                                                   &
    &                hb,   ubtb,   vbtb,      w,      r,              &
    &                ub,     vb,     tb,                              &
    &               amv,    ahv,                                      &
    &              taux,   tauy,     ft,   ptop,                      &
#ifdef OPT_BODY
    &                tq,                                              &
#endif
    &                ua,     va,     ta,                              &
    &                ha,   ubta,   vbta,                              &
    &                gx,     gy,     xx,     yy,                      &
    &               gxx,    gyy,                                      &
    &              uadv,   vadv,   wadv  )
       call putsig(                                                   &
    &                tb )
       return
    end if

    if ( ofinal ) then
       if ( myrank < ijnode ) then
          call iniset(                                                &
    &               uadv,   vadv,   wadv,      r,                     &
    &                 ub,     vb,     tb,                             &
    &                 hb,   ubtb,   vbtb,                             &
    &                  w,    amv,    ahv,                             &
    &               taux,   tauy,                                     &
    &                 ft )
       end if
       call predco(                                                   &
    &                 hb,   ubtb,   vbtb,      w,      r,             &
    &                 ub,     vb,     tb,                             &
    &                amv,    ahv,                                     &
    &               taux,   tauy,     ft,   ptop,                     &
#ifdef OPT_BODY
    &                 tq,                                             &
#endif
    &                 ua,     va,     ta,                             &
    &                 ha,   ubta,   vbta,                             &
    &                 gx,     gy,     xx,     yy,                     &
    &                gxx,    gyy,                                     &
    &               uadv,   vadv,   wadv )
       return
    end if

    tt    = tt1
    nt    = nt1
    ts    = ts1
    its   = its1
    itst  = itst1
    ntss  = ntss1
    tss   = tss1

    call chkstk( oflstk )

! *** Leap-frog ( A -> B ) ***

    if ( itst == 4 ) then
       call sfcflx(  ft,   taux,   tauy,     tb  )
#ifdef OPT_BODY
       call bdyflx(  tq,     tb  )
#endif
#ifdef OPT_TRIPOLE
       call shift2(                                                   &
    &              taux,   tauy,                                      &
    &             nxdim,  nydim,      1,                              &
    &             -1.d0,      0,      0 )
#else
       call shift2(                                                   &
    &                 taux,   tauy,                                   &
    &                nxdim,  nydim,    1)
#endif
       call puttao(                                                   &
    &              taux,   tauy,                                      &
    &             1.0d0,  0.0d0  )
       call predco(                                                   &
    &                hb,   ubtb,   vbtb,      w,      r,              &
    &                ub,     vb,     tb,                              &
    &               amv,    ahv,                                      &
    &              taux,   tauy,     ft,   ptop,                      &
#ifdef OPT_BODY
    &                tq,                                              &
#endif     
    &                ua,     va,     ta,                              &
    &                ha,   ubta,   vbta,                              &
    &                gx,     gy,     xx,     yy,                      &
    &               gxx,    gyy,                                      &
    &              uadv,   vadv,   wadv  )

! *** Leap-frog ( B -> A ) ***

    else if ( itst == 3 ) then
       call sfcflx(  ft,   taux,   tauy,     ta )
#ifdef OPT_BODY
       call bdyflx(  tq,     ta  )  
#endif
#ifdef OPT_TRIPOLE
       call shift2(                                                   &
    &              taux,   tauy,                                      &
    &             nxdim,  nydim,      1,                              &
    &             -1.d0,      0,      0 )
#else
       call shift2(                                                   &
    &              taux,   tauy,                                      &
    &             nxdim,  nydim,      1)
#endif
       call puttao(                                                   &
    &              taux,   tauy,                                      &
    &             1.0d0,  0.0d0  )
       call predco(                                                   &
    &                ha,   ubta,   vbta,      w,      r,              &
    &                ua,     va,     ta,                              &
    &               amv,    ahv,                                      &
    &              taux,   tauy,     ft,   ptop,                      &
#ifdef OPT_BODY 
    &                tq,                                              &
#endif
    &                ub,     vb,     tb,                              &
    &                hb,   ubtb,   vbtb,                              &
    &                gx,     gy,     xx,     yy,                      &
    &               gxx,    gyy,                                      &
    &              uadv,   vadv,   wadv )

! *** Euler-backward ( B -> B ) ***

    else if ( itst == 1 ) then
       call forsto(                                                   &
    &                ua,     va,     ta,                              &
    &                ha,   ubta,   vbta,                              &
    &                ub,     vb,     tb,                              &
    &                hb,   ubtb,   vbtb )
       call sfcflx(  ft,   taux,   tauy,   tb  )
#ifdef OPT_BODY
       call bdyflx(  tq,     tb )
#endif
#ifdef OPT_TRIPOLE
       call shift2(                                                   &
    &              taux,   tauy,                                      &
    &             nxdim,  nydim,      1,                              &
    &             -1.d0,      0,      0 )
#else
       call shift2(                                                   &
    &              taux,   tauy,                                      &
    &             nxdim,  nydim,      1)
#endif
       call puttao(                                                   &
    &              taux,   tauy,                                      &
    &             1.0d0,  0.0d0  )
       ieuler = 1
       call predco(                                                   &
    &                ha,   ubta,   vbta,      w,      r,              &
    &                ua,     va,     ta,                              &
    &               amv,    ahv,                                      &
    &              taux,   tauy,     ft,   ptop,                      &
#ifdef OPT_BODY
    &                tq,                                              &
#endif 
    &                ub,     vb,     tb,                              &
    &                hb,   ubtb,   vbtb,                              &
    &                gx,     gy,     xx,     yy,                      &
    &               gxx,    gyy,                                      &
    &              uadv,   vadv,   wadv )
       ieuler = 2
       call predco(                                                   &
    &                hb,   ubtb,   vbtb,      w,      r,              &
    &                ub,     vb,     tb,                              &
    &               amv,    ahv,                                      &
    &              taux,   tauy,     ft,   ptop,                      &
#ifdef OPT_BODY
    &                tq,                                              &
#endif
    &                ua,     va,     ta,                              &
    &                ha,   ubta,   vbta,                              &
    &                gx,     gy,     xx,     yy,                      &
    &               gxx,    gyy,                                      &
    &              uadv,   vadv,   wadv)

! *** Forward ( B -> B ) ***
    else
       ieuler = 0
       call forsto(                                                   &
    &                ua,     va,     ta,                              &
    &                ha,   ubta,   vbta,                              &
    &                ub,     vb,     tb,                              &
    &                hb,   ubtb,   vbtb  )
       call sfcflx(  ft,   taux,   tauy,    tb )
#ifdef OPT_BODY
       call bdyflx(  tq,     tb )
#endif
#ifdef OPT_TRIPOLE
       call shift2(                                                   &
    &              taux,   tauy,                                      &
    &             nxdim,  nydim,      1,                              &
    &             -1.d0,      0,      0 )
#else
       call shift2(                                                   &
    &              taux,   tauy,                                      &
    &             nxdim,  nydim,      1)
#endif
       call puttao(                                                   &
    &              taux,   tauy,                                      &
    &             1.0d0,  0.0d0  )
       call predco(                                                   &
    &                ha,   ubta,   vbta,      w,      r,              &
    &                ua,     va,     ta,                              &
    &               amv,    ahv,                                      &
    &              taux,   tauy,     ft,   ptop,                      &
#ifdef OPT_BODY
    &                tq,                                              &
#endif
    &                ub,     vb,     tb,                              &
    &                hb,   ubtb,   vbtb,                              &
    &                gx,     gy,     xx,     yy,                      &
    &               gxx,    gyy,                                      &
    &              uadv,   vadv,   wadv  )
       call excngo(                                                   &
    &                ua,     va,     ta,                              &
    &                ha,   ubta,   vbta,                              &
    &                ub,     vb,     tb,                              &
    &                hb,   ubtb,   vbtb  )
    end if

    call cofptu( &
    &               taux,   tauy)

! *** Output to file ***
    if ( itst == 3 ) then
       call putsig(                                                   &
    &                  ta )
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
    else
       call putsig(                                                   &
    &                  tb )
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

    call chkftx

! surface flux output
    call chksfx

    call chkout( oflout )

  end subroutine ocean

end module aocea
