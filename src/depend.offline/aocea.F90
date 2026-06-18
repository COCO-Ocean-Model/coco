module aocea

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

  subroutine ocstup(     t,      u,     v,    hb,   ahv,  &
     &                 tt1,    dt1 )

    use brdge
    use bshft
    use bstbc
    use qckot
    use sfcng
    use tflxt
    use tslvt
    use ufile

    implicit none
    
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
    
    call ocnvar(     t,     u,     v,    hb,   ahv )
#ifdef OPT_TRIPOLE
    call shift1(     hb,                    &
   &              nxdim,  nydim,      1,    &
   &               1.d0,      0,      0 )
#else
    call shift1(     hb,                    &
   &              nxdim,  nydim,      1)
#endif

  end subroutine ocstup

! =====================================================================

  subroutine ocean (    t,     ft,                         &
    &                   u,      v,     ha,     hb,    ahv, &
    &                 nt1,    tt1,  itst1,    ts1,   its1, &
    &               ntss1,   tss1,                         &
    &              oflout, oflstk)

    use aprdc
    use bshft
    use bstbc
    use qckot
    use sfcng

    implicit none

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

    call chkstk( oflstk )

    do ij = 1, nxydim
       ha(ij) = hb(ij)
    end do
    call ocnvar(     t,     u,     v,    hb,   ahv )
#ifdef OPT_TRIPOLE
    call shift1(     hb,                      &
    &             nxdim,  nydim,      1,      &
    &              1.d0,      0,      0 )     
    call shift1(      t,                      &
    &             nxdim,  nydim, nztdim,      &
    &              1.d0,      0,      0 )     
    call shift2(      u,      v,              &
    &             nxdim,  nydim,  nzdim,      &
    &             -1.d0,     -1,     -1 )     
    call shift1(    ahv,                      &
    &             nxdim,  nydim,  nzdim,      &
    &              1.d0,      0,      0 )
#else
    call shift1(     hb,                      &
    &             nxdim,  nydim,      1 )
    call shift1(      t,                      &
    &             nxdim,  nydim, nztdim )
    call shift3(      u,      v,    ahv,      &
    &             nxdim,  nydim,  nzdim)
#endif
    call stbctr(     t,     r )
    call sfcflx(    ft,     t )
#ifdef OPT_BODY
    call bdyflx(    tq,     t )
#endif
    call predco(     t,      &
    &               ft,      &
#ifdef OPT_BODY                  
    &               tq,      &
#endif
    &                u,     v,    ha,    hb,   ahv )

!    call south_bound_radiation(t(:,3:ntdim),v,0.d0)
!    call south_bound_fix(t(:,3:ntdim),0.d0)
    
! *** Output to file ***
    do l = 3, ntdim
       call chekin(      t(1, l),      ctrnam(l),  &
       &               ctrtit(l),      ctruni(l),  &
       &                nx,     ny,     nz, nxyzdm, 'OCLVTT')
   end do
   do l = 3, ntdim
      call chekin(     ft(1, l),       cftnam(l),  &
       &              cfttit(l),       cftuni(l),  &
       &               nx,     ny,      1, nxydim, 'OCSFCT')
   end do

   call chkout( oflout )

 end subroutine ocean

 subroutine south_bound_radiation(t_raw, v, t_bg)
   
   use zocnod, only: jrank
   use zocdim, only: istr, iend, jstr, kstr, kend
   use zocmsk, only: amskt, amskv
   use zocgrd, only: dzv, rym, hyu
   
   real(8), target, intent(inout) :: t_raw(nxyzdm, ntdim-2)
   real(8), intent(in)    :: t_bg
   real(8), intent(in)    :: v(nxyzdm)
   real(8), pointer :: t_ptr(:, :, :, :)
   real(8) :: vv, r
   integer :: i, j, k, n, ij, ijk

   t_ptr(1:nxdim, 1:nydim, 1:nzdim, 1:ntdim-2) => t_raw
   
   if (jrank == 0) then
      do n = 1, ntdim-2
      do k = kstr, kend
      do i = istr, iend
         ij = i + (jstr-1) * nxdim
         ijk = ij + (k-1) * nxydim
         vv = (v(ijk) * dzv(ij, k) * amskv(ij, k) + v(ijk-1) * dzv(ij-1, k) * amskv(ij-1, k)) * amskt(ij, k) * amskt(ij+nxdim, k)
         if (vv < 0.d0) then
            r = abs(vv) / (dzv(ij, k) * amskv(ij, k) + dzv(ij-1, k) * amskv(ij-1, k)) * ts * rym(ij) / (hyu(ij) + hyu(ij-1)) * 2.d0
            t_ptr(i, jstr+1, k, n) = (1.d0 - r) * t_ptr(i, jstr+1, k, n) + r * t_ptr(i, jstr+2, k, n)
         else
            t_ptr(i, jstr+1, k, n) = t_bg
         end if
         do j = 1, jstr
            t_ptr(i, j, k, n) = t_ptr(i, jstr+1, k, n)
         end do
      end do
      end do
      end do
   end if
   
 end subroutine south_bound_radiation

 subroutine south_bound_fix(t_raw, t_bg)
   
   use zocnod, only: jrank
   use zocdim, only: istr, iend, jstr, kstr, kend
   use zocmsk, only: amskt, amskv
   use zocgrd, only: dzv, rym, hyu
   
   real(8), target, intent(inout) :: t_raw(nxyzdm, ntdim-2)
   real(8), intent(in)    :: t_bg
   real(8), pointer :: t_ptr(:, :, :, :)
   integer :: i, j, k, n, ij, ijk
   integer :: ijstr_sbnd, ijend_sbnd

   t_ptr(1:nxdim, 1:nydim, 1:nzdim, 1:ntdim-2) => t_raw
   
   if (jrank == 0) then
      do n = 1, ntdim-2
      do k = kstr, kend
      do j = 1, jstr+1
         ijstr_sbnd = istr + (j - 1) * nxdim
         ijend_sbnd = iend + (j - 1) * nxdim
         t_ptr(istr:iend, j, k, n) = t_bg * amskt(ijstr_sbnd:ijend_sbnd, k) + t_ptr(istr:iend, j, k, n) * (1.d0 - amskt(ijstr_sbnd:ijend_sbnd, k))
      end do
      end do
      end do
   end if
   
 end subroutine south_bound_fix

end module aocea
