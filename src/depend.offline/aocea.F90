module aocea

! --- information -----------------------------------------------------
!
!  HISTORY
!     '02.10.10  H.Hasumi: from COCO3.4
!     '07.04.24  H.Hasumi
!     '14.08.05  M.Watanabe: Tripolar
!
! ---------------------------------------------------------------------     

   use zocdim, only: &
    &     nx,     ny,     nz, &
    & nxydim, nxyzdm, nxyidm,  nxdim,  nydim,  nzdim,  ntdim
   use zocgrd, only: &
    &     dx,     dy,    hxt,    hyt, &
    &     dt, &
    &     tt,     ts,    tss, &
    &     nt,    its,   itst,   ntss
   use zocfil, only: &
    & nfomax

   implicit none

   character(len=16), save :: ctrnam(ntdim), cftnam(ntdim)
   character(len=32), save :: ctrtit(ntdim), cfttit(ntdim)
   character(len=16), save :: ctruni(ntdim), cftuni(ntdim)

   private

   public :: ocstup, ocean

contains

subroutine ocstup( &
     &                t,      u,      v,     hb,    ahv, &
     &              tt1,    dt1)

   use brdge
   use bshft
   use bstbc
   use qckot
   use sfcng
   use tflxt
   use tslvt
   use ufile

   real(8), intent(out) ::       t(nxyzdm, ntdim)
   real(8), intent(out) ::       u(nxyzdm),      v(nxyzdm)
   real(8), intent(out) ::      hb(nxydim)
   real(8), intent(out) ::     ahv(nxyzdm)
   real(8), intent(in)  ::     tt1, dt1

   integer ::    ij,      l
   integer :: ifpar,  jfpar

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

   tt = tt1
   dt = dt1

   call rdgeo
   call svtset
   call chkset

   call ocnvar( &
   &                 t,      u,      v,     hb,    ahv)
#ifdef OPT_TRIPOLE
   call shift1(     hb, &
   &             nxdim,  nydim,      1, &
   &              1.d0,      0,      0 )
#else
   call shift1( &
   &                hb, &
   &             nxdim,  nydim,      1)
#endif

   return
end subroutine ocstup

! =====================================================================

subroutine ocean ( &
   &                    t,     ft, &
   &                    u,      v,     ha,     hb,    ahv, &
   &                  nt1,    tt1,  itst1,    ts1,   its1, &
   &                ntss1,   tss1, &
   &               oflout, oflstk)

   use aprdc
   use bshft
   use bstbc
   use qckot
   use sfcng

   real(8), intent(out) ::       t(nxyzdm, ntdim)
   real(8), intent(out) ::      ft(nxydim, ntdim)
   real(8), intent(out) ::       u(nxyzdm),      v(nxyzdm)
   real(8), intent(out) ::      ha(nxydim),     hb(nxydim)
   real(8), intent(out) ::     ahv(nxyzdm)
   real(8), intent(in)  ::     tt1
   real(8), intent(in)  ::    ts1,   tss1
   integer, intent(in)  ::    nt1,   its1,  itst1,  ntss1
   logical, intent(in)  :: oflout(nfomax), oflstk(nfomax)
   real(8)              ::       r(nxyzdm)
   integer ::    ij,      l
   
#ifdef OPT_BODY
   real(8),save ::      tq(nxyzdm, ntdim)
#endif

   tt    = tt1
   nt    = nt1
   ts    = ts1
   its   = its1
   itst  = itst1
   ntss  = ntss1
   tss   = tss1

   call chkstk( &
   &            oflstk)

   do ij = 1, nxydim
      ha(ij) = hb(ij)
   end do
   call ocnvar( &
   &                 t,      u,      v,     hb,    ahv)
#ifdef OPT_TRIPOLE
   call shift1(    hb, &
   &             nxdim,  nydim,      1, &
   &              1.d0,      0,      0 )
   call shift1(     t, &
   &             nxdim,  nydim, nztdim, &
   &              1.d0,      0,      0 )
   call shift2(     u,      v, &
   &             nxdim,  nydim,  nzdim, &
   &             -1.d0,     -1,     -1 )
   call shift1(   ahv, &
   &             nxdim,  nydim,  nzdim, &
   &              1.d0,      0,      0 )
#else
   call shift1( &
   &                hb, &
   &             nxdim,  nydim,      1)
   call shift1( &
   &                 t, &
   &             nxdim,  nydim,  nzdim*2)
   call shift3( &
   &                 u,      v,    ahv, &
   &             nxdim,  nydim,  nzdim)
#endif
   call stbctr(t,  r)
   call sfcflx( &
   &               ft, & ! only for passive tracer fluxes
   &                t)
#ifdef OPT_BODY
   call bdyflx( &
   &                tq, &
   &                 t)
#endif
   call predco( &
   &                    t, &
   &                   ft, &
#ifdef OPT_BODY
   &                   tq, &
#endif
   &                    u,      v,     ha,     hb,    ahv)

! *** Output to file ***
   do l = 3, ntdim
      call chekin(      t(1, l),      ctrnam(l), &
     &                ctrtit(l),      ctruni(l), &
     &                nx,     ny,     nz, nxyzdm, 'OCLVTT')
   end do
   do l = 3, ntdim
     call chekin(      ft(1, l),      cftnam(l), &
     &                cfttit(l),      cftuni(l), &
     &               nx,     ny,      1, nxydim, 'OCSFCT')
   end do

   call chkout( &
   &            oflout)

   return
end subroutine ocean

end module aocea
