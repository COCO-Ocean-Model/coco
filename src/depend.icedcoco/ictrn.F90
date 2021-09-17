module ictrn

! ---- information ----------------------------------------------------
!
!  Sea ice thickness category transfer
!
!  HISTORY
!     '03.07.30  H.Hasumi
!     '07.09.25  H.Hasumi: for COCO4
!     '08.06.11  H.Hasumi: initial/final processing
!     '08.07.10  H.Hasumi: initial/final processing
!     '12.07.19  Y.Komuro: for COCO5.0
!     '21.05.26  Y.Komuro: snow aging & meltpond parametrization 
!
! ---------------------------------------------------------------------

  use zocdim, only: &
    &     nx,     ny, nxydim,  ntdim,    nic,    kstr, ijtstr, ijtend, &
    &  oinit, ofinal
  use zocgrd, only: &
    &    hic,     ts
  use zocmsk, only: &
    &  amskt
  use zocphy, only: &
    &   rhoo,   rhoi,   rhos,   dtds,    cpi,    cpo,   hfus

  implicit none

  real(8), save ::    tmi

  real(8), save ::  amin = 1.0d-6,  amax = 1.d0,  si = 5.d0
  real(8), save :: hilmt = 1.0d4, hiref = 1.0d4
  integer, save ::  mic = nic
  logical, save :: ohiflt = .false.

  namelist /nmamin/ amin, amax, mic
  namelist /nmislt/ si
  namelist /nmhflt/ ohiflt, hilmt, hiref

  private

  public :: ictrns, icadjs, ichflt

contains

subroutine ictrns( &
  &                    ax,    hix,    hsx,    eix,    tix, &
  &                   asx,  frlvx,   vmpx,  frmpx,   dsdx,   dsbx, &
  &                    ft,     fs,    fdd,    fdb )
  use ufile
  use zocite

  real(8), intent(inout) ::     ax(nxydim, 0:nic)
  real(8), intent(inout) ::    hix(nxydim, 0:nic)
  real(8), intent(inout) ::    hsx(nxydim, 0:nic)
  real(8), intent(inout) ::    eix(nxydim, 0:nic),    tix(nxydim, 0:nic)
  real(8), intent(inout) ::    asx(nxydim, 0:nic),  frlvx(nxydim, 0:nic)
  real(8), intent(inout) ::   vmpx(nxydim, 0:nic),  frmpx(nxydim, 0:nic)
  real(8), intent(inout) ::   dsdx(nxydim, 0:nic),   dsbx(nxydim, 0:nic)
  real(8), intent(inout) ::     ft(nxydim, ntdim),     fs(nxydim)
  real(8), intent(inout) ::    fdd(nxydim),    fdb(nxydim)

  real(8) ::  axhix(nxydim, 0:nic),  axhsx(nxydim, 0:nic)
  real(8) ::  axeix(nxydim, 0:nic)
  real(8) ::  axasx(nxydim, 0:nic),  axvmp(nxydim, 0:nic)
  real(8) ::  axflv(nxydim, 0:nic),  axfmp(nxydim, 0:nic)
  real(8) ::  axdsd(nxydim, 0:nic),  axdsb(nxydim, 0:nic)
  real(8) ::     ci(nxydim)
!  common /work/ axhix, axhsx, axeix, ci

  real(8), save ::    rri,    rrs
  logical, save :: ofirst = .true.

  integer ::     ij,      k
  integer ::  ifpar,  jfpar,  istat


  if (oinit .or. ofinal) then
     return
  end if

  if (ofirst) then
     ofirst = .false.
     call rewnml(ifpar, jfpar)
     read (ifpar, nmamin, iostat=istat)
     call cstnml(jfpar, 'ictrns', 'nmamin', istat)
     write(jfpar, nmamin)
     call rewnml(ifpar, jfpar)
     read (ifpar, nmislt, iostat=istat)
     call cstnml(jfpar, 'ictrns', 'nmislt', istat)
     write(jfpar, nmislt)
     call rewnml(ifpar, jfpar)
     read (ifpar, nmhflt, iostat=istat)
     call cstnml(jfpar, 'ictrns', 'nmhflt', istat)
     write(jfpar, nmhflt)

     tmi = dtds * si
     rri = rhoo / rhoi
     rrs = rhoo / rhos
  end if

  do ij = ijtstr, ijtend
     axhix(ij, 0) = ax(ij, 0) * hix(ij, 0)
     axhsx(ij, 0) = ax(ij, 0) * hsx(ij, 0)
     axeix(ij, 0) = ax(ij, 0) * eix(ij, 0)
     axflv(ij, 0) = ax(ij, 0) * frlvx(ij, 0)
     axvmp(ij, 0) = ax(ij, 0) * vmpx(ij, 0)
     axfmp(ij, 0) = ax(ij, 0) * frmpx(ij, 0)
     axdsd(ij, 0) = ax(ij, 0) * dsdx(ij, 0)
     axdsb(ij, 0) = ax(ij, 0) * dsbx(ij, 0)
  end do

  do k = 1, nic-1
     do ij = ijtstr, ijtend
        if (hix(ij, k) .ge. hic(k+1)) then
!        if (      (ax(ij, k) .gt. 0.d0)
!          & .and. (hix(ij, k) .ge. hic(k+1))) then
           axhix(ij, k+1) = ax(ij, k+1) * hix(ij, k+1) &
             &            + ax(ij, k) * hix(ij, k)
           axhsx(ij, k+1) = ax(ij, k+1) * hsx(ij, k+1) &
             &            + ax(ij, k) * hsx(ij, k)
           axeix(ij, k+1) = ax(ij, k+1) * eix(ij, k+1) &
             &            + ax(ij, k) * eix(ij, k)
           axasx(ij, k+1) = ax(ij, k+1) * asx(ij, k+1) &
             &            + ax(ij, k) * asx(ij, k)
           axflv(ij, k+1) = ax(ij, k+1) * frlvx(ij, k+1) &
             &            + ax(ij, k) * frlvx(ij, k)
           axvmp(ij, k+1) = ax(ij, k+1) * vmpx(ij, k+1) &
             &            + ax(ij, k) * vmpx(ij, k)
           axfmp(ij, k+1) = ax(ij, k+1) * frmpx(ij, k+1) &
             &            + ax(ij, k) * frmpx(ij, k)
           axdsd(ij, k+1) = ax(ij, k+1) * dsdx(ij, k+1) &
             &            + ax(ij, k) * dsdx(ij, k)
           axdsb(ij, k+1) = ax(ij, k+1) * dsbx(ij, k+1) &
             &            + ax(ij, k) * dsbx(ij, k)             
           ax(ij, k+1) = ax(ij, k+1) + ax(ij, k)
           ax(ij, k) = 0.d0
           hix(ij, k) = hic(k)
           hsx(ij, k) = 0.d0
           eix(ij, k) = 0.d0
           tix(ij, k) = tmi
           asx(ij, k) = 0.d0
           frlvx(ij, k) = 1.d0
           vmpx(ij, k) = 0.d0
           frmpx(ij, k) = 0.d0
           dsdx(ij, k) = 0.d0
           dsbx(ij, k) = 0.d0
           hix(ij, k+1) = axhix(ij, k+1) / ax(ij, k+1)
           hsx(ij, k+1) = axhsx(ij, k+1) / ax(ij, k+1)
           eix(ij, k+1) = axeix(ij, k+1) / ax(ij, k+1)
           tix(ij, k+1) = ti(eix(ij, k+1)/hix(ij, k+1), si)
           asx(ij, k+1) = axasx(ij, k+1) / ax(ij, k+1)
           frlvx(ij, k+1) = axflv(ij, k+1) / ax(ij, k+1)
           vmpx(ij, k+1) = axvmp(ij, k+1) / ax(ij, k+1)
           frmpx(ij, k+1) = axfmp(ij, k+1) / ax(ij, k+1)
           dsdx(ij, k+1) = axdsd(ij, k+1) / ax(ij, k+1)
           dsbx(ij, k+1) = axdsb(ij, k+1) / ax(ij, k+1)
        end if
     end do
  end do
  do k = nic-1, 1, -1
     do ij = ijtstr, ijtend
        if (hix(ij, k+1) .lt. hic(k+1)) then
!        if (      (hix(ij, k+1) .lt. hic(k+1)) &
!          & .and. (ax(ij, k+1) .gt. 0.d0)) then
           axhix(ij, k) = ax(ij, k+1) * hix(ij, k+1) &
             &          + ax(ij, k) * hix(ij, k)
           axhsx(ij, k) = ax(ij, k+1) * hsx(ij, k+1) &
             &          + ax(ij, k) * hsx(ij, k)
           axeix(ij, k) = ax(ij, k+1) * eix(ij, k+1) &
             &          + ax(ij, k) * eix(ij, k)
           axasx(ij, k) = ax(ij, k+1) * asx(ij, k+1) &
             &          + ax(ij, k) * asx(ij, k)
           axflv(ij, k) = ax(ij, k+1) * frlvx(ij, k+1) &
             &          + ax(ij, k) * frlvx(ij, k)
           axvmp(ij, k) = ax(ij, k+1) * vmpx(ij, k+1) &
             &          + ax(ij, k) * vmpx(ij, k)
           axfmp(ij, k) = ax(ij, k+1) * frmpx(ij, k+1) &
             &          + ax(ij, k) * frmpx(ij, k)
           axdsd(ij, k) = ax(ij, k+1) * dsdx(ij, k+1) &
             &          + ax(ij, k) * dsdx(ij, k)
           axdsb(ij, k) = ax(ij, k+1) * dsbx(ij, k+1) &
             &          + ax(ij, k) * dsbx(ij, k)
           ax(ij, k) = ax(ij, k+1) + ax(ij, k)
           ax(ij, k+1) = 0.d0
           hix(ij, k) = axhix(ij, k) / ax(ij, k)
           hsx(ij, k) = axhsx(ij, k) / ax(ij, k)
           eix(ij, k) = axeix(ij, k) / ax(ij, k)
           tix(ij, k) = ti(eix(ij, k)/hix(ij, k), si)
           asx(ij, k) = axasx(ij, k) / ax(ij, k)
           frlvx(ij, k) = axflv(ij, k) / ax(ij, k)
           vmpx(ij, k) = axvmp(ij, k) / ax(ij, k)
           frmpx(ij, k) = axfmp(ij, k) / ax(ij, k)
           dsdx(ij, k) = axdsd(ij, k) / ax(ij, k)
           dsbx(ij, k) = axdsb(ij, k) / ax(ij, k)
           hix(ij, k+1) = hic(k+1)
           hsx(ij, k+1) = 0.d0
           eix(ij, k+1) = 0.d0
           tix(ij, k+1) = tmi
           asx(ij, k+1) = 0.d0
           frlvx(ij, k+1) = 1.d0
           vmpx(ij, k+1) = 0.d0
           frmpx(ij, k+1) = 0.d0
           dsdx(ij, k+1) = 0.d0
           dsbx(ij, k+1) = 0.d0
        end if
     end do
  end do
  do ij = ijtstr, ijtend
     if (hix(ij, 1) .lt. hic(1)) then
!     if (      (ax(ij, 1) .gt. 0.d0) &
!       & .and. (hix(ij, 1) .lt. hic(1))) then
        axhix(ij, 1) = ax(ij, 1) * hix(ij, 1)
        axhsx(ij, 1) = ax(ij, 1) * hsx(ij, 1)
        axeix(ij, 1) = ax(ij, 1) * eix(ij ,1)
        axvmp(ij, 1) = ax(ij, 1) * vmpx(ij, 1)
        axdsd(ij, 1) = ax(ij, 1) * dsdx(ij, 1)
        axdsb(ij, 1) = ax(ij, 1) * dsbx(ij, 1)
        hix(ij, 1) = hic(1)
        ax(ij, 1) = axhix(ij, 1) / hix(ij, 1)
        hsx(ij, 1) = axhsx(ij, 1) / ax(ij, 1)
        eix(ij, 1) = axeix(ij ,1) / ax(ij, 1)
        vmpx(ij, 1) = axvmp(ij, 1) / ax(ij, 1)
        dsdx(ij, 1) = axdsd(ij, 1) / ax(ij, 1)
        dsbx(ij, 1) = axdsb(ij, 1) / ax(ij, 1)
!       The following variables do not change:
!       ice temperature, snow age, 
!       level ice fraction, and melt pond fraction.
      end if
  end do

  do k = 1, nic
     do ij = ijtstr, ijtend
        if ((ax(ij, k) .gt. 0.d0) .and. (ax(ij, k) .lt. amin)) then
           axhix(ij, 0) = axhix(ij, 0) &
             &             + ax(ij, k) * hix(ij, k)
           axhsx(ij, 0) = axhsx(ij, 0) &
             &             + ax(ij, k) * hsx(ij, k)
           axeix(ij, 0) = axeix(ij, 0) &
             &             + ax(ij, k) * eix(ij, k)
           axvmp(ij, 0) = axvmp(ij, 0) &
             &          + ax(ij, k) * vmpx(ij, k)
           axdsd(ij, 0) = axdsd(ij, 0) &
             &          + ax(ij, k) * dsdx(ij, k)
           axdsb(ij, 0) = axdsb(ij, 0) &
             &          + ax(ij, k) * dsbx(ij, k)
           ax(ij, k) = 0.d0
           hix(ij, k) = hic(k)
           hsx(ij, k) = 0.d0
           eix(ij, k) = 0.d0
           tix(ij, k) = tmi
           asx(ij, k) = 0.d0
           frlvx(ij, k) = 1.d0
           vmpx(ij, k) = 0.d0
           frmpx(ij, k) = 0.d0
           dsdx(ij, k) = 0.d0
           dsbx(ij, k) = 0.d0
        end if
     end do
  end do

  do ij = ijtstr, ijtend
!     if (hix(ij, 0) .gt. 0.d0) then
     ft(ij, 2) = ft(ij, 2) &
       &       - (axhix(ij, 0) / rri + axhsx(ij, 0) / rrs) &
       &         / ts * amskt(ij, kstr)
     fs(ij) = fs(ij) &
       &    - axhix(ij, 0) / rri / ts * si * amskt(ij, kstr)
     ft(ij, 1) = ft(ij, 1) &
       &       - (  axhsx(ij, 0) / rrs * hfus &
       &         + axeix(ij, 0) / rri) / ts / cpo * &
       &       amskt(ij, kstr)
!     end if
!     melt pond is a virtual reservor
     fdd(ij) = fdd(ij) - axdsd(ij, 0) / ts * amskt(ij, kstr)
     fdb(ij) = fdb(ij) - axdsb(ij, 0) / ts * amskt(ij, kstr)
  end do

!  entry ic0set( &
!    &               ax,    hix,    hsx)

! *** initialize the category 0 ***
  do ij = ijtstr, ijtend
     ax(ij, 0) = 1.d0
     hix(ij, 0) = 0.d0
     hsx(ij, 0) = 0.d0
     eix(ij, 0) = 0.d0
     tix(ij, 0) = tmi
     frlvx(ij, 0) = 1.d0
     vmpx(ij, 0) = 0.d0
     frmpx(ij, 0) = 0.d0
     dsdx(ij, 0) = 0.d0
     dsbx(ij, 0) = 0.d0
  end do
  do k = 1, nic
     do ij = ijtstr, ijtend
        ax(ij, 0) = ax(ij, 0) - ax(ij, k)
     end do
  end do

! *** because of the finit precision, ax(ij, 0) could be a small negative
! *** value, which causes some problems
  do ij = ijtstr, ijtend
     ax(ij, 0) = max(0.d0, ax(ij ,0))
  end do

  return

end subroutine ictrns

! =====================================================================

subroutine icadjs( &
  &                    ax,    hix,    hsx,    eix, &
  &                  vmpx,   dsdx,   dsbx )

  real(8), intent(inout) ::     ax(nxydim, 0:nic)
  real(8), intent(inout) ::    hix(nxydim, 0:nic)
  real(8), intent(inout) ::    hsx(nxydim, 0:nic)
  real(8), intent(inout) ::    eix(nxydim, 0:nic)
  real(8), intent(inout) ::   vmpx(nxydim, 0:nic)
  real(8), intent(inout) ::   dsdx(nxydim, 0:nic),   dsbx(nxydim, 0:nic)

  real(8) ::  axhix(nxydim, 0:nic),  axhsx(nxydim, 0:nic)
  real(8) ::  axeix(nxydim, 0:nic)
  real(8) ::  axvmp(nxydim, 0:nic)
  real(8) ::  axdsd(nxydim, 0:nic),  axdsb(nxydim, 0:nic)
  real(8) ::     ci(nxydim)

  integer ::     ij,      k

  if (oinit .or. ofinal) then
     return
  end if

  if (amax .ge. 1.d0) then
     return
  end if

  do k = 1, mic
     do ij = ijtstr, ijtend
        axhix(ij, k) = ax(ij, k) * hix(ij, k)
        axhsx(ij, k) = ax(ij, k) * hsx(ij, k)
        axeix(ij, k) = ax(ij, k) * eix(ij, k)
        axvmp(ij, k) = ax(ij, k) * vmpx(ij, k)
        axdsd(ij, k) = ax(ij, k) * dsdx(ij, k)
        axdsb(ij, k) = ax(ij, k) * dsbx(ij, k)
     end do
  end do

  do ij = ijtstr, ijtend
     ci(ij) = 0.d0
     ax(ij, 0) = 1.d0
  end do
  do k = 1, mic
     do ij = ijtstr, ijtend
        ci(ij) = ci(ij) + ax(ij, k)
     end do
  end do

  do k = 1, mic
     do ij = ijtstr, ijtend
        if ((ci(ij) .gt. amax) .and. (ax(ij, k) .gt. 0.d0)) then
           ax(ij, k) = ax(ij, k) * amax / ci(ij)
           hix(ij, k) = axhix(ij, k) / ax(ij, k)
           hsx(ij, k) = axhsx(ij, k) / ax(ij, k)
           eix(ij, k) = axeix(ij, k) / ax(ij, k)
           vmpx(ij, k) = axvmp(ij, k) / ax(ij, k)
           dsdx(ij, k) = axdsd(ij, k) / ax(ij, k)
           dsbx(ij, k) = axdsb(ij, k) / ax(ij, k)
        end if
     end do
  end do

  do k = 1, nic
     do ij = ijtstr, ijtend
        ax(ij, 0) = ax(ij, 0) - ax(ij, k)
     end do
  end do

  return

end subroutine icadjs

! =====================================================================

subroutine ichflt( &
  &                    ax,    hix,    hsx,    eix,    tix, &
  &                   asx,  frlvx,   vmpx,  frmpx,   dsdx,   dsbx )
  use qckot
  use zocite

  real(8), intent(inout) ::     ax(nxydim, 0:nic)
  real(8), intent(inout) ::    hix(nxydim, 0:nic)
  real(8), intent(inout) ::    hsx(nxydim, 0:nic)
  real(8), intent(inout) ::    eix(nxydim, 0:nic)
  real(8), intent(inout) ::    tix(nxydim, 0:nic)
  real(8), intent(inout) ::    asx(nxydim, 0:nic)
  real(8), intent(inout) ::  frlvx(nxydim, 0:nic)
  real(8), intent(inout) ::   vmpx(nxydim, 0:nic)
  real(8), intent(inout) ::  frmpx(nxydim, 0:nic)
  real(8), intent(inout) ::   dsdx(nxydim, 0:nic)
  real(8), intent(inout) ::   dsbx(nxydim, 0:nic)

  real(8) ::  axhix(nxydim, 0:nic),  axhsx(nxydim, 0:nic)
  real(8) ::  axeix(nxydim, 0:nic)
  real(8) ::  axasx(nxydim, 0:nic),  axvmp(nxydim, 0:nic)
  real(8) ::  axflv(nxydim, 0:nic),  axfmp(nxydim, 0:nic)
  real(8) ::  axdsd(nxydim, 0:nic),  axdsb(nxydim, 0:nic)
  real(8) ::     ci(nxydim)
  real(8) :: daxhix(nxydim)

  real(8) ::    fax,    lax,  fdahi
  real(8) ::   fahi,   fahs,   faei
  real(8) ::   faas,  favmp,  fadsd, fadsb
  real(8) ::  faflv,  fafmp
  real(8) :: rdaxhi(nxydim)
  
  integer ::     ij,      k

  if (oinit .or. ofinal) then
     return
  end if

  if (.not.ohiflt) return

  do k = 1, nic
     do ij = ijtstr, ijtend
        axhix(ij, k) = ax(ij, k) * hix(ij, k)
        axhsx(ij, k) = ax(ij, k) * hsx(ij, k)
        axeix(ij, k) = ax(ij, k) * eix(ij, k)
        axvmp(ij, k) = ax(ij, k) * vmpx(ij, k)
        axdsd(ij, k) = ax(ij, k) * dsdx(ij, k)
        axdsb(ij, k) = ax(ij, k) * dsbx(ij, k)
     end do
  end do

  do ij = 1, nxydim
     rdaxhi(ij) = 0.0d0
  end do

  do ij = ijtstr, ijtend
     daxhix(ij) = ax(ij, nic) * max((hix(ij, nic) - hilmt), 0.0d0)
     rdaxhi(ij) = daxhix(ij)
     lax = max((ax(ij,0) - (1.0d0-amax)), 0.0d0)
     fax = min(lax, daxhix(ij) / hiref)
     ax(ij, nic) = ax(ij, nic) + fax
     daxhix(ij) = daxhix(ij) - fax * hiref
  end do

! The following variables change only with inter-category transfer:
!  snow age, level ice fraction, and melt pond fraction
  do k = 1, nic
     do ij = ijtstr, ijtend
        axasx(ij, k) = ax(ij, k) * asx(ij, k)
        axflv(ij, k) = ax(ij, k) * frlvx(ij, k)
        axfmp(ij, k) = ax(ij, k) * frmpx(ij, k)
     end do
  end do

  do k = nic-1, 1, -1
     do ij = ijtstr, ijtend
        if (ax(ij, k) > 0.0d0) then
           fax = min(ax(ij, k), daxhix(ij) / (hiref - hix(ij, k)))
           fdahi = fax * (hiref - hix(ij, k))
           daxhix(ij) = daxhix(ij) - fdahi
           fahi = fax * hix(ij, k)
           fahs = fax * hsx(ij, k)
           faei = fax * eix(ij, k)
           faas = fax * asx(ij, k)
           faflv = fax * frlvx(ij, k)
           favmp = fax * vmpx(ij, k)
           fafmp = fax * frmpx(ij, k)
           fadsd = fax * dsdx(ij, k)
           fadsb = fax * dsbx(ij, k)
           ax(ij, k) = ax(ij, k) - fax
           axhix(ij, k) = axhix(ij, k) - fahi
           axhsx(ij, k) = axhsx(ij, k) - fahs
           axeix(ij, k) = axeix(ij, k) - faei
           axasx(ij, k) = axasx(ij, k) - faas
           axflv(ij, k) = axflv(ij, k) - faflv
           axvmp(ij, k) = axvmp(ij, k) - favmp
           axfmp(ij, k) = axfmp(ij, k) - fafmp
           axdsd(ij, k) = axdsd(ij, k) - fadsd
           axdsb(ij, k) = axdsb(ij, k) - fadsb
           ax(ij, nic) = ax(ij, nic) + fax
           axhix(ij, nic) = axhix(ij, nic) + fahi
           axhsx(ij, nic) = axhsx(ij, nic) + fahs
           axeix(ij, nic) = axeix(ij, nic) + faei
           axasx(ij, nic) = axasx(ij, nic) + faas
           axflv(ij, nic) = axflv(ij, nic) + faflv
           axvmp(ij, nic) = axvmp(ij, nic) + favmp
           axfmp(ij, nic) = axfmp(ij, nic) + fafmp
           axdsd(ij, nic) = axdsd(ij, nic) + fadsd
           axdsb(ij, nic) = axdsb(ij, nic) + fadsb
        endif
     end do
  end do

  do ij = ijtstr, ijtend
     ax(ij, 0) = 1.d0
  end do

  do k = 1, nic
     do ij = ijtstr, ijtend
        ax(ij, 0) = ax(ij, 0) - ax(ij, k)
        if (ax(ij, k) .le. 0.d0) then
           ax(ij, k) = 0.d0
           hix(ij, k) = hic(k)
           hsx(ij, k) = 0.d0
           tix(ij, k) = tmi
           eix(ij, k) = 0.d0
           asx(ij, k) = 0.d0
           frlvx(ij, k) = 1.d0
           vmpx(ij, k) = 0.d0
           frmpx(ij, k) = 0.d0
           dsdx(ij, k) = 0.d0
           dsbx(ij, k) = 0.d0
        else
           hix(ij, k) = axhix(ij, k) / ax(ij, k)
           hsx(ij, k) = axhsx(ij, k) / ax(ij, k)
           eix(ij, k) = axeix(ij, k) / ax(ij, k)
           tix(ij, k) = ti(eix(ij, k)/hix(ij, k), si)
           asx(ij, k) = axasx(ij, k) / ax(ij, k)
           frlvx(ij, k) = axflv(ij, k) / ax(ij, k)
           vmpx(ij, k) = axvmp(ij, k) / ax(ij, k)
           frmpx(ij, k) = axfmp(ij, k) / ax(ij, k)
           dsdx(ij, k) = axdsd(ij, k) / ax(ij, k)
           dsbx(ij, k) = axdsb(ij, k) / ax(ij, k)
        endif
     end do
  end do

! Subtract the remaining daxhix from rdaxhi
  do ij = ijtstr, ijtend
     rdaxhi(ij) = rdaxhi(ij) - daxhix(ij)
  end do

  call chekin( rdaxhi, 'HLMAHI', &
    &           '', '', &
    &               nx,     ny,      1, nxydim, 'OCSFCT')

  return

end subroutine ichflt

end module ictrn
