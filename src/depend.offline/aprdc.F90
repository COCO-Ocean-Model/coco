module aprdc

! --- information -----------------------------------------------------
!
!  HISTORY
!     '02.10.10  H.Hasumi: from MIROC3.1-OMIP
!     '07.04.24  H.Hasumi
!     '14.08.05  M.Watanabe: Tripolar
!
! ---------------------------------------------------------------------

   use zocdim,  only  :   nxdim,    nydim,             &
      &                   nxyzdm,   nxydim,   ntdim,   &
      &                   nztdim,      nic
   use zocnod,  only  :   myrank,   ijnode

   implicit none

contains

   subroutine predco(                                                  &
      &                t,                                              &
      &               ft,                                              &
#ifdef OPT_BODY
      &               tq,                                              &
#endif
      &                u,      v,     ha,     hb,    ahv)

   use bchmk
   use bshft
   use bstbc
   use dwdns
   use tflxt
   use tslvt
   use ucloc
   use zocdim, only : nx, ny
   use qckot
   
    implicit none

    real(8), intent(inout) ::       t(nxyzdm, ntdim)
    real(8), intent(inout) ::      ft(nxydim, ntdim)
    real(8), intent(inout) ::       u(nxyzdm),      v(nxyzdm)
    real(8), intent(in)    ::      ha(nxydim),     hb(nxydim)
    real(8), intent(in)    ::     ahv(nxyzdm)

!---- local variables
    real(8) ::       w(nxyzdm),     hc(nxydim)
    real(8) ::   diffz(nxyzdm)
    real(8) ::     adt(nxyzdm, ntdim)
    real(8) ::       r(nxyzdm)
#ifdef OPT_BODY
    real(8) ::      tq(nxyzdm, ntdim)
#endif

    integer(4) ::    ijk,      n

    if (myrank >= ijnode) return
    
      call clcstr('HDIAG')
#ifdef OPT_BBL
      call rmmskt
      call admktb
#endif
      call wdenst(     w,    hc,     u,     v,    ha )
      call clcend('HDIAG')

      call clcstr('TRACE')
      call flxtrc( &
        &            adt,  diffz,             &
        &              t,     ha,             &
        &             hb,     hc,             &
        &              u,      v,      w,     &
        &            ahv)
#ifdef OPT_BBL
      call flxtrb(                            &
        &            adt,  diffz,             &
        &              t,                     &
        &              u,      v,      w,     &
        &            ahv)
      call rmmskt
      call admkt1
      call stbbgt(                            &
        &            adt,  diffz)
#endif
#ifdef OPT_BODY
      do n = 3, ntdim
         do ijk = 1, nxyzdm
            adt(ijk, n) = adt(ijk, n) + tq(ijk, n)
         end do
      end do
#endif
      call slvtrc(     t,   adt, diffz,    ha,    hb,    hc,    ft )
      call tundif(     t,    hc )
#ifdef OPT_BBL
      call stbbtr(                                         &
        &              t)
      call rmmskt
      call admktb
      call stbbt2(                                         &
        &              t)
      call rmmskt
      call admkt1
#endif

#ifdef OPT_TRIPOLE
      call shift1(                            &
        &              t,                     &
        &          nxdim,   nydim, nztdim,    &
        &           1.D0,       0,      0 )
#else
      call shift1(                            &
        &              t,                     &
        &          nxdim,   nydim, nztdim)
#endif
      call stbctr( &
        &              t,       r)
    call clcend('TRACE')

!---- for debug (tracer convervation)
!    call chekin(    hc,   'SH', &
!         &        'sea surface height',   'cm', &
!         &          nx,     ny,      1, nxydim, 'OCSFCT' )
!-----
    
  end subroutine predco

end module aprdc
