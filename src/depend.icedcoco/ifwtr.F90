module ifwtr

! --- information -----------------------------------------------------
!
!  HISTORY
!     '02.10.10  H.Hasumi: from MIROC3.1
!     '07.04.23  H.Hasumi
!     '07.09.25  H.Hasumi: multi-category sea ice
!     '08.06.11  H.Hasumi: initial/final processing
!     '08.07.10  H.Hasumi: initial/final processing
!     '12.07.19  Y.Komuro: for COCO5.0
!
! ---------------------------------------------------------------------

  use zocdim, only: &
    & nxydim,    nic, ijtstr, ijtend, &
    &  oinit, ofinal
  use zocgrd, only: &
    &     ts
  use zocphy, only: &
    &   rhoo,   rhoi,   rhos

  implicit none

  private
  public :: fwater

contains

subroutine fwater( &
  &                    ax,    hix,    hsx, &
  &                  prec,   snow, &
  &                  evap,   subi, adjlat, &
  &                   wev,    wsb,   soff )

  real(8), intent(inout) ::     ax(nxydim, 0:nic),    hix(nxydim, 0:nic)
  real(8), intent(inout) ::    hsx(nxydim, 0:nic)
  real(8), intent(inout) ::   prec(nxydim),   snow(nxydim)
  real(8), intent(inout) ::    wsb(nxydim, nic)
  real(8), intent(out)   ::   evap(nxydim), adjlat(nxydim)
  real(8), intent(out)   ::   subi(nxydim, nic)
  real(8), intent(in)    ::    wev(nxydim),   soff(nxydim)
  
  real(8), save ::    rri,    rrs
  logical, save ::  ofirst = .true.

  real(8) ::   wdif
  real(8) ::    dhs,    dhi
  integer ::     ij,      k
  integer ::  ifpar,  jfpar

  real(8) ::     az(nxydim, 0:nic),    hiz(nxydim, 0:nic)
  real(8) ::    hsz(nxydim, 0:nic)
!      common /work/ az, hiz, hsz

  if (oinit .or. ofinal) then
     return
  end if

  if (ofirst) then
     ofirst = .false.
     rri = rhoo / rhoi
     rrs = rhoo / rhos
  end if

  do k = 0, nic
     do ij = 1, nxydim
        az (ij, k) = ax(ij, k)
        hiz(ij, k) = hix(ij, k)
        hsz(ij, k) = hsx(ij, k)
     end do
  end do

  do k = 1, nic
     do ij = ijtstr, ijtend
        if (az(ij, k) .gt. 0.d0) then
           if (hsz(ij, k) .gt. 0.d0) then
              dhs = ts * rrs * wsb(ij, k) / ax(ij, k)
              hsx(ij, k) = max(hsz(ij, k) - dhs, 0.d0)
              wsb(ij, k) = wsb(ij, k) &
                &        + ax(ij, k) * (hsx(ij, k) - hsz(ij, k)) &
                &          / rrs / ts
           end if
           dhi = ts * rri * wsb(ij, k) / ax(ij, k)
           hix(ij, k) = max(hiz(ij, k) - dhi, 0.d0)
           subi(ij, k) = ax(ij, k) * (hiz(ij, k) - hix(ij, k)) &
             &           / rri / ts
           wsb(ij, k) = wsb(ij, k) - subi(ij, k)
        end if
     end do
  end do

  do k = 1, nic
     do ij = ijtstr, ijtend
        evap(ij) = evap(ij) + wsb(ij, k)
     end do
  end do
  do ij = ijtstr, ijtend
     adjlat(ij) = evap(ij) - wev(ij)
  end do

  do ij = ijtstr, ijtend
!         prec(ij) = prec(ij) - snow(ij)
     snow(ij) = snow(ij) + soff(ij)
  end do

  return

end subroutine fwater

end module ifwtr
