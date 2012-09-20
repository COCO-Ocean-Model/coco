module binst

! --- information -----------------------------------------------------
!
!  HISTORY
!     '99.08.17  H.Hasumi: from CCSR2-MASK
!     '01.02.16  H.Hasumi: for incorporating BBL model
!     '02.05.29  H.Nakano: tracer dimension
!     '02.06.02  H.Hasumi: combine parallel and nonparallel
!     '07.04.23  H.Hasumi
!     '09.02.23  Y.Komuro: call DDENST(diagnose R) instead of OVTURN
!                          bug fix (SHIFT UX/VX after STBBVT)
!     '10.04.14  M.kurogi
!     '10.04.14  M.Kurogi: (COCO4.4 tripolar code by Dr. Suzuki)
!
! ---------------------------------------------------------------------

  implicit none
  private

  public  ::  iniset  !  aocea

contains

  subroutine iniset(                                                  &
    &     uadv,   vadv,   wadv,      r,                               &
    &       ub,     vb,     tb,                                       &
    &       hb,   ubtb,   vbtb,                                       &
    &        w,    amv,    ahv,                                       &
    &     taux,   tauy,                                               &
    &       ft)

    use zocdim,  only  :                                              &
    &     nxdim,  nydim,  nxyzdm,  nxydim,  ntdim, nztdim
    use zocgrd
    use zocphy
    use zocmsk,  only  :  amskv
    use tovtr
    use dvlva
    
    implicit none

    real(8),    intent(inout)  ::   uadv(nxyzdm),  vadv(nxyzdm)
    real(8),    intent(inout)  ::   wadv(nxyzdm),     r(nxyzdm)
    real(8),    intent(inout)  ::     ub(nxyzdm),    vb(nxyzdm)
    real(8),    intent(in)     ::     tb(nxyzdm,ntdim)
    real(8),    intent(in)     ::     hb(nxydim)
    real(8),    intent(in)     ::   ubtb(nxydim),  vbtb(nxydim)
    real(8),    intent(in)     ::      w(nxyzdm)
    real(8),    intent(in)     ::    amv(nxyzdm),   ahv(nxyzdm)
    real(8),    intent(in)     ::   taux(nxydim),  tauy(nxydim)
    real(8),    intent(in)     ::     ft(nxydim,ntdim)

!---- local variables    
    real(8)       ::   uz(nxydim),    vz(nxydim)
    integer(4)    ::  ijk,   ij,    k

#ifdef OPT_TRIPOLE
    call shift1(    hb,                                               &
    &            nxdim,  nydim,      1,                               &
    &             1.d0,      0,      0 )
    call shift2(  ubtb,   vbtb,                                       &
    &            nxdim,  nydim,      1,                               &
    &            -1.d0,     -1,     -1 )
    call shift2(    ub,     vb,                                       &
    &            nxdim,  nydim,  nzdim,                               &
    &            -1.d0,     -1,     -1 ) 
    call shift1(     w,                                               &
    &            nxdim,  nydim,  nzdim,                               &
    &             1.d0,      0,      0 )
#else
    call shift3(                                                      &
    &               hb,   ubtb,   vbtb,                               &
    &            nxdim,  nydim,      1)
    call shift3(                                                      &
    &               ub,     vb,      w,                               &
    &            nxdim,  nydim,  nzdim)
#endif

    r(1:nxyzdm) = 0.d0

    call ddenst( r,  tb )
#ifdef OPT_BBL
    call stbbtr( tb )
#endif

#ifdef OPT_TRIPOLE
    call shift1(     r,                                               &
    &            nxdim,  nydim,  nzdim,                               &
    &             1.d0,      0,      0 )
    call shift1(    tb,                                               &
    &            nxdim,  nydim, nztdim,                               &
    &             1.d0,      0,      0 )
#else
    call shift1(                                                      &
    &                r,                                               &
    &            nxdim,  nydim,  nzdim)
    call shift1(                                                      &
    &               tb,                                               &
    &            nxdim,  nydim, nztdim)
#endif

    call stbctr( tb, r )
    
#ifdef OPT_BBL
    call rmmskv
    call admkvb
    do k = 1,  nzdim
       do ij = 1, nxydim
          ijk    =  ij + ( k - 1 ) * nxydim
          ub(ijk) = ub(ijk) * amskv(ij,k)
          vb(ijk) = vb(ijk) * amskv(ij,k)
       end do
    end do
    call rmmskt
    call admktb
    call admskv
#endif
    call velvad(                                                      &
    &                 uadv,   vadv,   wadv,                           &
    &                   ub,     vb,      w)
#ifdef OPT_BBL
    call stbbvt( ub, vb )
#ifdef OPT_TRIPOLE
    call shift2(    ub,     vb,                                       &
    &            nxdim,  nydim,  nzdim,                               &
    &            -1.d0,     -1,     -1 )
#else
    call shift2(                                                      &
    &               ub,     vb,
    &            nxdim,  nydim,  nzdim)
#endif
    call rmmskt
    call admkt1
    call rmmskv
    call admkv1
#endif

#ifdef OPT_TRIPOLE
    call shift1(   uadv,                                              &
    &             nxdim,  nydim,  nzdim,                              &
    &             -1.d0,      0,     -1 )
    call shift1(   vadv,                                              &
    &             nxdim,  nydim,  nzdim,                              &
    &             -1.d0,     -1,      0 )
    call shift1(   wadv,                                              &
    &             nxdim,  nydim,  nzdim,                              &
    &              1.d0,     -1,     -1 )
    call shift1(   amv,                                               &
    &            nxdim,  nydim,  nzdim,                               &
    &             1.d0,     -1,     -1 )
    call shift1(   ahv,                                               &
    &            nxdim,  nydim,  nzdim,                               &
    &             1.d0,      0,      0 )
    call shift2(                                                      &
    &             taux,   tauy,                                       &
    &            nxdim,  nydim,      1,                               &
    &            -1.d0,      0,      0 )
#else
    call shift3(                                                      &
    &             uadv,   vadv,   wadv,                               &
    &            nxdim,  nydim,  nzdim)
    call shift2(                                                      &
    &              amv,    ahv,                                       &
    &            nxdim,  nydim,  nzdim)
    call shift2(                                                      &
    &             taux,   tauy,                                       &
    &            nxdim,  nydim,      1)
#endif

  end subroutine iniset

end module binst

