module tovtr

! --- information -----------------------------------------------------
!
!  convective adjustment for the unstable water column
!
!  history
!     '99.08.16  h.hasumi: from ccsr2-mask
!     '00.05.30  h.hasumi: change in the initialization of r
!     '00.12.14  h.hasumi: for hybrid vertical coordinate
!     '01.01.30  h.hasumi: for partial step bottom topography
!     '01.09.17  h.hasumi: gamma for bbl
!     '01.12.07  h.hasumi
!     '02.05.29  h.nakano: tracer dimension
!     '07.04.23  h.hasumi
!     '07.05.01  h.hasumi: mcdougall et al. (2003) eq. of state
!     '07.09.25  h.hasumi: arguments of chekin
!     '08.06.11  h.hasumi: initial/final processing
!     '08.07.10  h.hasumi: initial/final processing
!     '09.02.23  y.komuro: add ddenst (only diagnosing r)
!
! ---------------------------------------------------------------------

  use zocdim, only: &
    &      nx,     ny,     nz,     kz, &
    &  nxydim, nxyzdm,  nzdim,  ntdim, &
    &  ijtstr, ijtend,   kstr,   kend, &
    &   oinit, ofinal
  use zocgrd, only: &
    &      dz,     ds,   zbot
  use zocmsk, only: &
#ifdef OPT_BBL
    & amsktb, &
#endif
    & nbot

  implicit none
  private

  real(8), save :: c0(nzdim), c1(nzdim), c2(nzdim)
  real(8), save :: c3(nzdim), c4(nzdim), c5(nzdim), c6(nzdim)
  real(8), save :: d0(nzdim), d1(nzdim), d2(nzdim), d3(nzdim), d4(nzdim)
  real(8), save :: d5(nzdim), d6(nzdim), d7(nzdim), d8(nzdim), d9(nzdim)

  real(8), save :: gamma(nz) = 1.d0
  namelist /nmacct/ gamma

  public :: ovtset, ovturn, ddenst

contains 

subroutine ovtset( &
  &      r,      t )

  real(8), intent(out) ::       r(nxydim, nzdim)
  real(8), intent(in)  ::       t(nxydim, nzdim, ntdim)

  integer ::  ifpar,  jfpar,  istat
  integer ::     ij,      k
  real(8) ::     tl,     sl
  real(8) ::     p1,     p2

  call rewnml(ifpar, jfpar)
  write(jfpar, *) '*** OVTSET ***'
  read(ifpar, nmacct, iostat=istat)

  call secoef( &
     & c0(kstr), c1(kstr), c2(kstr), c3(kstr), &
     & c4(kstr), c5(kstr), c6(kstr), &
     & d0(kstr), d1(kstr), d2(kstr), d3(kstr),  d4(kstr), &
     & d5(kstr), d6(kstr), d7(kstr), d8(kstr),  d9(kstr) )
  write(jfpar, *) '*** the equation of state determined ***'

  do k = 1, nzdim
     do ij = 1, nxydim
        r(ij, k) = 0.d0
     end do
  end do
  do k = kstr, kend
     do ij = ijtstr, ijtend
        tl = t(ij, k, 1)
        sl = t(ij, k, 2)
        p1 = c0(k) &
           & + (c1(k) + (c2(k) + c3(k) * tl) * tl) * tl &
           & + (c4(k) + c5(k) * tl + c6(k) * sl) * sl
        p2 = d0(k) &
           & + (d1(k) + (d2(k) + (d3(k) + d4(k) * tl) * tl) * tl) * tl &
           & + (d5(k) + (d6(k) + d7(k) * tl * tl) * tl &
           &          + (d8(k) + d9(k) * tl * tl) * sqrt(sl)) * sl
        r(ij, k) = p1 / p2 - 1.d3
     end do
  end do

  return

end subroutine ovtset
! =====================================================================

subroutine ovturn( &
  &      r,      t,      h )

  real(8), intent(out)    ::       r(nxydim, nzdim)
  real(8), intent(inout)  ::       t(nxydim, nzdim, ntdim)
  real(8), intent(in)     ::       h(nxydim)

  real(8), save ::     zt(nxydim, nzdim)
  real(8) ::   conv(nxydim, nzdim)
  real(8) ::  dzsig(nxydim, nzdim)
  real(8) ::    ttl(nxydim, ntdim),     w2(nxydim)
  integer ::    lup(nxydim)
  real(8) ::     tu,     tl,     su,     sl
  real(8) ::     ru,     rl
  real(8) ::     p1,     p2
  integer ::     ij,      k,     kk,      n

! common /work/ ttl, w2, conv, dzsig, lup

  if (oinit .or. ofinal) then
     return
  end if

  do k = 1, nzdim
     do ij = 1, nxydim
        conv(ij, k) = 0.d0
     end do
  end do

  do k = kstr, kstr+kz-1
     do ij = 1, nxydim
        dzsig(ij, k) = ds(k) * (h(ij) + zbot) * gamma(k-kstr+1)
     end do
  end do
  do k = kstr+kz, kend
     do ij = 1, nxydim
        dzsig(ij, k) = dz(ij, k) * gamma(k-kstr+1)
     end do
  end do
#ifdef OPT_BBL
  do ij = ijtstr, ijtend
     k = nbot(ij)
     dzsig(ij, k) = dz(ij, k) * gamma(nz) * amsktb(ij) &
        &         + dzsig(ij, k) * (1.d0 - amsktb(ij))
  end do
#endif

  do ij = 1, nxydim
     zt(ij, kstr) = 0.d0
  end do
  do k = kstr, kend
     do ij = 1, nxydim
        zt(ij, k+1) = zt(ij, k) + dzsig(ij, k)
     end do
  end do

  do n = 1, ntdim
     do ij = ijtstr, ijtend
        ttl(ij, n) = t(ij, kstr, n) * dzsig(ij, kstr)
        lup(ij) = kstr
     end do
  end do

  do ij = ijtstr, ijtend
     do k = kstr+1, nbot(ij)
        tu = t(ij, k-1, 1)
        su = t(ij, k-1, 2)
        tl = t(ij, k, 1)
        sl = t(ij, k, 2)
        p1 = c0(k) &
           & + (c1(k) + (c2(k) + c3(k) * tu) * tu) * tu &
           & + (c4(k) + c5(k) * tu + c6(k) * su) * su
        p2 = d0(k) &
           & + (d1(k) + (d2(k) + (d3(k) + d4(k) * tu) * tu) * tu) * tu &
           & + (d5(k) + (d6(k) + d7(k) * tu * tu) * tu &
           &          + (d8(k) + d9(k) * tu * tu) * sqrt(su)) * su
        ru = p1 / p2
        p1 = c0(k) &
           & + (c1(k) + (c2(k) + c3(k) * tl) * tl) * tl &
           & + (c4(k) + c5(k) * tl + c6(k) * sl) * sl
        p2 = d0(k) &
           & + (d1(k) + (d2(k) + (d3(k) + d4(k) * tl) * tl) * tl) * tl &
           & + (d5(k) + (d6(k) + d7(k) * tl * tl) * tl &
           &          + (d8(k) + d9(k) * tl * tl) * sqrt(sl)) * sl
        rl = p1 / p2
        if (ru .gt. rl) then
           w2 (ij) = 1.d0 / (zt(ij, k+1) - zt(ij, lup(ij)))
           conv(ij, k) = 1.d0
           do n = 1, ntdim
              ttl(ij, n) = ttl(ij, n) &
                 &       + t(ij, k, n) * dzsig(ij, k)
              do kk = kstr, k
                 if (kk .ge. lup(ij)) then
                    t(ij, kk, n) = ttl(ij, n) * w2(ij)
                 end if
              end do
           end do
        else
           do n = 1, ntdim  
              ttl(ij, n) = t(ij, k, n) * dzsig(ij, k)
              lup(ij) = k
           end do
        end if
     end do
  end do

  do k = kstr, kend
     do ij = ijtstr, ijtend
        tl = t(ij, k, 1)
        sl = t(ij, k, 2)
        p1 = c0(k) &
           & + (c1(k) + (c2(k) + c3(k) * tl) * tl) * tl &
           & + (c4(k) + c5(k) * tl + c6(k) * sl) * sl
        p2 = d0(k) &
           & + (d1(k) + (d2(k) + (d3(k) + d4(k) * tl) * tl) * tl) * tl &
           & + (d5(k) + (d6(k) + d7(k) * tl * tl) * tl &
           &          + (d8(k) + d9(k) * tl * tl) * sqrt(sl)) * sl 
        r(ij, k) = p1 / p2 - 1.d3
     end do
  end do

  call chekin(  conv, 'CONV', &
     &            nx,     ny,     nz, nxyzdm, 'OCN')

  return

end subroutine ovturn
! =====================================================================

subroutine ddenst( &
  &      r,      t )

  real(8), intent(out) ::       r(nxydim, nzdim)
  real(8), intent(in)  ::       t(nxydim, nzdim, ntdim)

  integer ::     ij,      k
  real(8) ::     tl,     sl
  real(8) ::     p1,     p2

  do k = kstr, kend
     do ij = ijtstr, ijtend
        tl = t(ij, k, 1)
        sl = t(ij, k, 2)
        p1 = c0(k) &
           & + (c1(k) + (c2(k) + c3(k) * tl) * tl) * tl &
           & + (c4(k) + c5(k) * tl + c6(k) * sl) * sl
        p2 = d0(k) &
           & + (d1(k) + (d2(k) + (d3(k) + d4(k) * tl) * tl) * tl) * tl &
           & + (d5(k) + (d6(k) + d7(k) * tl * tl) * tl &
           &          + (d8(k) + d9(k) * tl * tl) * sqrt(sl)) * sl
        r(ij, k) = p1 / p2 - 1.d3
     end do
  end do

  return

end subroutine ddenst

end module tovtr
