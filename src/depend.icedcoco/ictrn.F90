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
!
! ---------------------------------------------------------------------

  use zocdim, only: &
    & nxydim,  ntdim,    nic,    kstr, ijtstr, ijtend, &
    &  oinit, ofinal
  use zocgrd, only: &
    &    hic,     ts
  use zocmsk, only: &
    &  amskt
  use zocphy, only: &
    &   rhoo,   rhoi,   rhos,   dtds,    cpi,    cpo,   hfus

  implicit none

  real(8), save ::  amin = 1.0d-6,  amax = 1.d0,  si = 5.d0
  integer, save ::  mic = nic

  namelist /nmamin/ amin, amax, mic
  namelist /nmislt/ si

  private

  public :: ictrns, icadjs

contains

subroutine ictrns( &
  &                    ax,    hix,    hsx,    eix,    tix, &
  &                    ft,     fs )
  use ufile
  use zocite

  real(8), intent(inout) ::     ax(nxydim, 0:nic)
  real(8), intent(inout) ::    hix(nxydim, 0:nic)
  real(8), intent(inout) ::    hsx(nxydim, 0:nic)
  real(8), intent(inout) ::    eix(nxydim, 0:nic),    tix(nxydim, 0:nic)
  real(8), intent(inout) ::     ft(nxydim, ntdim),     fs(nxydim)

  real(8) ::  axhix(nxydim, 0:nic),  axhsx(nxydim, 0:nic)
  real(8) ::  axeix(nxydim, 0:nic)
  real(8) ::     ci(nxydim)
!  common /work/ axhix, axhsx, axeix, ci

  real(8), save ::    rri,    rrs
  real(8), save ::    tmi
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

     tmi = dtds * si
     rri = rhoo / rhoi
     rrs = rhoo / rhos
  end if

  do ij = ijtstr, ijtend
     axhix(ij, 0) = ax(ij, 0) * hix(ij, 0)
     axhsx(ij, 0) = ax(ij, 0) * hsx(ij, 0)
     axeix(ij, 0) = ax(ij, 0) * eix(ij, 0)
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
           ax(ij, k+1) = ax(ij, k+1) + ax(ij, k)
           ax(ij, k) = 0.d0
           hix(ij, k) = hic(k)
           hsx(ij, k) = 0.d0
           eix(ij, k) = 0.d0
           tix(ij, k) = tmi
           hix(ij, k+1) = axhix(ij, k+1) / ax(ij, k+1)
           hsx(ij, k+1) = axhsx(ij, k+1) / ax(ij, k+1)
           eix(ij, k+1) = axeix(ij, k+1) / ax(ij, k+1)
           tix(ij, k+1) = ti(eix(ij, k+1)/hix(ij, k+1), si)
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
           ax(ij, k) = ax(ij, k+1) + ax(ij, k)
           ax(ij, k+1) = 0.d0
           hix(ij, k) = axhix(ij, k) / ax(ij, k)
           hsx(ij, k) = axhsx(ij, k) / ax(ij, k)
           eix(ij, k) = axeix(ij, k) / ax(ij, k)
           tix(ij, k) = ti(eix(ij, k)/hix(ij, k), si)
           hix(ij, k+1) = hic(k+1)
           hsx(ij, k+1) = 0.d0
           eix(ij, k+1) = 0.d0
           tix(ij, k+1) = tmi
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
        hix(ij, 1) = hic(1)
        ax(ij, 1) = axhix(ij, 1) / hix(ij, 1)
        hsx(ij, 1) = axhsx(ij, 1) / ax(ij, 1)
        eix(ij, 1) = axeix(ij ,1) / ax(ij, 1)
!       ice temperature does not change
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
           ax(ij, k) = 0.d0
           hix(ij, k) = hic(k)
           hsx(ij, k) = 0.d0
           eix(ij, k) = 0.d0
           tix(ij, k) = tmi
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
  &                    ax,    hix,    hsx,    eix )

  real(8), intent(inout) ::     ax(nxydim, 0:nic)
  real(8), intent(inout) ::    hix(nxydim, 0:nic)
  real(8), intent(inout) ::    hsx(nxydim, 0:nic)
  real(8), intent(inout) ::    eix(nxydim, 0:nic)

  real(8) ::  axhix(nxydim, 0:nic),  axhsx(nxydim, 0:nic)
  real(8) ::  axeix(nxydim, 0:nic)
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

end module ictrn
