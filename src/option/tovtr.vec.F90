module tovtr

! --- information -----------------------------------------------------
!
!  Convective adjustment for the unstable water column
!
!  HISTORY
!     '99.08.16  H.Hasumi: from CCSR2-MASK
!     '00.05.30  H.Hasumi: change in the initialization of R
!     '00.12.14  H.Hasumi: for hybrid vertical coordinate
!     '01.01.30  H.Hasumi: for partial step bottom topography
!     '01.09.17  H.Hasumi: GAMMA for BBL
!     '01.12.07  H.Hasumi
!     '02.05.29  H.Nakano: tracer dimension
!     '02.09.30  H.Hasumi: for vector machines
!     '07.04.23  H.Hasumi
!     '07.05.01  H.Hasumi: McDougall et al. (2003) eq. of state
!     '08.06.11  H.Hasumi: initial/final processing
!     '08.07.10  H.Hasumi: initial/final processing
!     '09.02.23  Y.Komuro: add DDENST (only diagnosing R)
!     '12.08.02  Y.Komuro: for COCO5.0
!
! ---------------------------------------------------------------------

  use zocdim, only: &
    &      nx,     ny,     nz,     kz, &
    &  nxydim, nxyzdm,  nzdim,  ntdim, &
    &  ijtstr, ijtend,   kstr,   kend, &
    &   oinit, ofinal
  use zocgrd, only: &
    &      dz,    dzm,     ds,   zbot
  use zocmsk, only: &
#ifdef OPT_BBL
    & amsktb, &
#endif
    &   amskt,  amftz,    nbot
  use zocphy, only: &
    & gravit,   rhoo

  implicit none
  private

  real(8), save :: c0(nzdim), c1(nzdim), c2(nzdim)
  real(8), save :: c3(nzdim), c4(nzdim), c5(nzdim), c6(nzdim)
  real(8), save :: d0(nzdim), d1(nzdim), d2(nzdim), d3(nzdim), d4(nzdim)
  real(8), save :: d5(nzdim), d6(nzdim), d7(nzdim), d8(nzdim), d9(nzdim)

  real(8), save :: c0s, c1s, c2s, c3s, c4s, c5s, c6s, &
    &              d0s, d1s, d2s, d3s, d4s, d5s, d6s, d7s, d8s, d9s
  real(8), save :: depth(nxydim, nzdim)

  real(8), save :: gamma(nz) = 1.d0
  real(8), save :: bcrit = 0.03d0
  integer, save :: kref = 1
  namelist /nmacct/ gamma
  namelist /nmmldt/ bcrit, kref

  public :: ovtset, ovturn, ddenst

contains 

subroutine ovtset( &
  &      r,      t )

  use xprst
  use ufile

  implicit none

  real(8), intent(out) ::       r(nxydim, nzdim)
  real(8), intent(in)  ::       t(nxydim, nzdim, ntdim)

  integer ::  ifpar,  jfpar,  istat
  integer ::     ij,      k
  real(8) ::     tl,     sl
  real(8) ::     p1,     p2

  call rewnml(ifpar, jfpar)
  write(jfpar, *) '*** OVTSET ***'
  read(ifpar, nmacct, iostat=istat)
  call cstnml(jfpar, 'ovtset', 'nmacct', istat)
  write(jfpar, nmacct)
  call rewnml(ifpar, jfpar)
  read(ifpar, nmmldt, iostat=istat)
  call cstnml(jfpar, 'ovtset', 'nmmldt', istat)
  write(jfpar, nmmldt)

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

  c0s=c0(kstr)
  c1s=c1(kstr)
  c2s=c2(kstr)
  c3s=c3(kstr)
  c4s=c4(kstr)
  c5s=c5(kstr)
  c6s=c6(kstr)
  d0s=d0(kstr)
  d1s=d1(kstr)
  d2s=d2(kstr)
  d3s=d3(kstr)
  d4s=d4(kstr)
  d5s=d5(kstr)
  d6s=d6(kstr)
  d7s=d7(kstr)
  d8s=d8(kstr)
  d9s=d9(kstr)

  do k = 1, nzdim
     do ij = 1, nxydim
        depth(ij, k) = 0.d0
     end do
  end do

  do k = kstr, kend
     do ij = 1, nxydim
        depth(ij, k) = depth(ij, k-1) + dz(ij, k)
     end do
  end do

  return

end subroutine ovtset
! =====================================================================

subroutine ovturn( &
  &      r,      t,      h )
  use qckot

  real(8), intent(out)    ::       r(nxydim, nzdim)
  real(8), intent(inout)  ::       t(nxydim, nzdim, ntdim)
  real(8), intent(in)     ::       h(nxydim)

  real(8), save ::     zt(nxydim, nzdim)
  real(8) ::   conv(nxydim, nzdim)
  real(8) :: cnvdep(nxydim)
  real(8) ::  dzsig(nxydim, nzdim)
  real(8) ::    ttl(nxydim, ntdim),     w2(nxydim)
  real(8) ::      n2(nxydim, nzdim),    mld(nxydim)
  real(8) :: dptmsig(nxydim, nzdim), dptsig(nxydim, nzdim)
  real(8) ::  dzmsig(nxydim, nzdim),   delb(nxydim, nzdim)
  integer ::    lup(nxydim)
  logical ::    lov(nxydim)
  real(8) ::     tu,     tl,     su,     sl,     tr,     sr
  real(8) ::     ru,     rl,     rr
  real(8) ::     p1,     p2
  integer ::     ij,      k,     kk,      n
  logical :: obtmld

! common /work/ ttl, w2, conv, dzsig, lup, lov

  if (oinit .or. ofinal) then
     return
  end if

  do k = 1, nzdim
     do ij = 1, nxydim
        dptmsig(ij, k) = 0.d0
        dptsig(ij, k) = 0.d0
        dzmsig(ij, k) = 0.d0
        delb(ij, k) = 0.d0
     end do
  end do
  do ij = 1, nxydim
     dzmsig(ij, kstr) = (h(ij) + zbot) * 0.5d0 * ds(kstr)
     dptmsig(ij, kstr) = dzmsig(ij, kstr)
     dptsig(ij, kstr) = (h(ij) + zbot) * ds(kstr)
  end do
  do k = kstr+1, kstr+kz-1
     do ij = 1, nxydim
        dzmsig(ij, k) = (h(ij) + zbot) * 0.5d0 * (ds(k-1) + ds(k))
        dptmsig(ij, k) = dptmsig(ij, k-1) + dzmsig(ij, k)
        dptsig(ij, k) = dptsig(ij, k-1) + (h(ij) + zbot) * ds(k)
     end do
  end do
  k = kstr+kz
  do ij = 1, nxydim
     dzmsig(ij, k) = (h(ij) + zbot) * 0.5d0 * ds(k-1) &
       &           + 0.5d0 * dz(ij, k)
     dptmsig(ij, k) = dptmsig(ij, k-1) + dzmsig(ij, k)
     dptsig(ij, k) = dptsig(ij, k-1) + dz(ij, k)
  end do
  do k = kstr+kz+1, kend
     do ij = 1, nxydim
        dzmsig(ij, k) = 0.5d0 * (dz(ij, k-1) + dz(ij, k))
        dptmsig(ij, k) = dptmsig(ij, k-1) + dzmsig(ij, k)
        dptsig(ij, k) = dptsig(ij, k-1) + dz(ij, k)
     end do
  end do
      
  do k = 1, nzdim
     do ij = 1, nxydim
        conv(ij, k) = 0.d0
        n2(ij, k) = 0.d0
     end do
  end do
  do ij = 1, nxydim
     mld(ij) = 0.d0
  end do

  do ij = 1, nxydim
     cnvdep(ij) = 0.d0
  end do

  do k = kstr, kstr+kz-1
     do ij = 1, nxydim
!        dzsig(ij, k) = ds(k) * (h(ij) + zbot) * gamma(k-kstr+1)
        dzsig(ij, k)=ds(k)*(h(ij)*amskt(ij,kstr)+zbot)*gamma(k-kstr+1)
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

  do k = kstr+1, kend
     do ij = ijtstr, ijtend
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
        ru = p1 / p2 * amskt(ij, k)
        p1 = c0(k) &
           & + (c1(k) + (c2(k) + c3(k) * tl) * tl) * tl &
           & + (c4(k) + c5(k) * tl + c6(k) * sl) * sl
        p2 = d0(k) &
           & + (d1(k) + (d2(k) + (d3(k) + d4(k) * tl) * tl) * tl) * tl &
           & + (d5(k) + (d6(k) + d7(k) * tl * tl) * tl &
           &          + (d8(k) + d9(k) * tl * tl) * sqrt(sl)) * sl
        rl = p1 / p2 * amskt(ij, k)
        if (ru > rl) then
           conv(ij, k) = 1.d0
           cnvdep(ij)=depth(ij,k)
           lov(ij) = .true.
           do n = 1, ntdim
              ttl(ij, n) = ttl(ij, n) &
                &        + t(ij, k, n) * dzsig(ij, k)
           end do
        else
           lup(ij) = k
           lov(ij) = .false.
           do n = 1, ntdim
              ttl(ij, n) = t(ij, k, n) * dzsig(ij, k)
           end do
        end if
        w2 (ij) = 1.d0 / (zt(ij, k+1) - zt(ij, lup(ij)))
     end do

     do kk = kstr, k
        do ij = ijtstr, ijtend
           if ((kk .ge. lup(ij)) .and. lov(ij)) then
              do n = 1, ntdim
                 t(ij, kk, n) = ttl(ij, n) * w2(ij)
              end do
           end if
        end do
     end do
  end do

! do k = kstr, kend
  k = kstr
  do ij = ijtstr, ijtend
     tl = t(ij, k, 1) * amskt(ij, k)
     sl = t(ij, k, 2) * amskt(ij, k)
     p1 = c0s &
     &  + (c1s + (c2s + c3s * tl) * tl) * tl &
     &  + (c4s + c5s * tl + c6s * sl) * sl
     p2 = d0s &
     &  + (d1s + (d2s + (d3s + d4s * tl) * tl) * tl) * tl &
     &  + (d5s + (d6s + d7s * tl * tl) * tl &
     &           + (d8s + d9s * tl * tl) * sqrt(sl)) * sl
     r(ij, k) = p1 / p2 - 1.d3
  enddo
! enddo

  do k = kstr+1, kend
     do ij = ijtstr, ijtend
        tl = t(ij, k, 1) * amskt(ij, k)
        sl = t(ij, k, 2) * amskt(ij, k)
        tu = t(ij, k-1, 1) * amskt(ij,k-1)
        su = t(ij, k-1, 2) * amskt(ij,k-1)
        tr = t(ij, kref+kstr-1, 1) * amskt(ij,kref+kstr-1)
        sr = t(ij, kref+kstr-1, 2) * amskt(ij,kref+kstr-1)
        p1 = c0(k) &
        &  + (c1(k) + (c2(k) + c3(k) * tl) * tl) * tl &
        &  + (c4(k) + c5(k) * tl + c6(k) * sl) * sl
        p2 = d0(k) &
        &  + (d1(k) + (d2(k) + (d3(k) + d4(k) * tl) * tl) * tl) * tl &
        &  + (d5(k) + (d6(k) + d7(k) * tl * tl) * tl &
        &           + (d8(k) + d9(k) * tl * tl) * sqrt(sl)) * sl
        rl = p1 / p2
        r(ij, k) = rl - 1.d3
        p1 = c0(k) &
        &  + (c1(k) + (c2(k) + c3(k) * tu) * tu) * tu &
        &  + (c4(k) + c5(k) * tu + c6(k) * su) * su
        p2 = d0(k) &
        &  + (d1(k) + (d2(k) + (d3(k) + d4(k) * tu) * tu) * tu) * tu &
        &  + (d5(k) + (d6(k) + d7(k) * tu * tu) * tu &
        &           + (d8(k) + d9(k) * tu * tu) * sqrt(su)) * su
        ru = p1 / p2
        p1 = c0(k) &
        &  + (c1(k) + (c2(k) + c3(k) * tr) * tr) * tr &
        &  + (c4(k) + c5(k) * tr + c6(k) * sr) * sr
        p2 = d0(k) &
        &  + (d1(k) + (d2(k) + (d3(k) + d4(k) * tr) * tr) * tr) * tr &
        &  + (d5(k) + (d6(k) + d7(k) * tr * tr) * tr &
        &           + (d8(k) + d9(k) * tr * tr) * sqrt(sr)) * sr
        rr = p1 / p2
        n2(ij, k) = - gravit / rhoo &
        &         * 1.0d-3 * (ru - rl) / dzm(ij, k) * amftz(ij, k)
        delb(ij, k) = - gravit * (rr - rl) / rl
     end do
  end do

  do ij = ijtstr, ijtend
     obtmld = .true.
     do k = kref+kstr, nbot(ij)
        if (delb(ij, k) > bcrit) then
           mld(ij) = dptmsig(ij, k-1) &
             &    + (bcrit-delb(ij,k-1)) / (delb(ij,k)-delb(ij,k-1)) &
             &    * dzmsig(ij, k)
           obtmld = .false.
           exit
        end if
     end do
     if (obtmld) then
        mld(ij) = dptsig(ij, nbot(ij))
     end if
  end do

  call chekin(    n2,     'N2', &
    &            'square of buoyancy frequency', '1/s^2', &
    &            nx,     ny,     nz, nxyzdm, 'OCLVMT')
  call chekin(   mld,    'MLD', &
    &            'mixed layer depth', 'cm', &
    &            nx,     ny,      1, nxydim, 'OCSFCT')
  call chekin(  conv,   'CONV', &
    &            'ocean conv. adj. index', 'N.D.', &
    &            nx,     ny,     nz, nxyzdm, 'OCLVMT')
  call chekin(cnvdep, 'CNVDEP', &
    &            'ocean conv. depth', 'cm', &
    &            nx,     ny,      1, nxydim, 'OCSFCT')

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
