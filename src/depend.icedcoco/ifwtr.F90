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
!     '21.05.26  Y.Komuro: snow aging & meltpond parametrization 
!
! ---------------------------------------------------------------------

  use zocdim, only: &
    & nxydim,    nic, ijtstr, ijtend, &
    &  oinit, ofinal
  use zocgrd, only: &
    &     ts,    hic
  use zocphy, only: &
    &   rhoo,   rhoi,   rhos,   hfus,   dtds

  implicit none

  real(8), save :: si = 5.d0

  namelist /nmislt/ si

  private
  public :: fwater

contains

subroutine fwater( &
  &                    ax,    hix,    hsx,    eix,    tix, &
  &                   asx,  frlvx,   vmpx,  frmpx,   dsdx,   dsbx, &
  &                  prec,   snow,    fdd,    fdb, &
  &                inrsbi, inrsbs, imraji, imrajs, &
  &                  evap,   subi, adjlat, wiadjs, &
  &                   wev,    wsb,   soff )

  use ufile
  use zocite

  real(8), intent(inout) ::     ax(nxydim, 0:nic),    hix(nxydim, 0:nic)
  real(8), intent(inout) ::    hsx(nxydim, 0:nic)
  real(8), intent(inout) ::    eix(nxydim, 0:nic),    tix(nxydim, 0:nic)
  real(8), intent(inout) ::    asx(nxydim, 0:nic),  frlvx(nxydim, 0:nic)
  real(8), intent(inout) ::   vmpx(nxydim, 0:nic),  frmpx(nxydim, 0:nic)
  real(8), intent(inout) ::   dsdx(nxydim, 0:nic),   dsbx(nxydim, 0:nic)
  real(8), intent(inout) ::   prec(nxydim),   snow(nxydim)
  real(8), intent(inout) ::    fdd(nxydim),    fdb(nxydim)
  real(8), intent(inout) ::    wsb(nxydim, nic)
  real(8), intent(inout) :: inrsbi(nxydim), inrsbs(nxydim)
  real(8), intent(inout) :: imraji(nxydim), imrajs(nxydim)
  real(8), intent(out)   ::   evap(nxydim), adjlat(nxydim), wiadjs(nxydim)
  real(8), intent(out)   ::   subi(nxydim, nic)
  real(8), intent(in)    ::    wev(nxydim),   soff(nxydim)
  
  real(8), save ::    rri,    rrs,    tmi
  logical, save ::  ofirst = .true.

  real(8) ::   wdif
  real(8) ::    dhs,    dhi,    dei,    dti, dhimax
  integer ::     ij,      k
  integer ::  ifpar,  jfpar,  istat

  real(8) ::     az(nxydim, 0:nic),    hiz(nxydim, 0:nic)
  real(8) ::    hsz(nxydim, 0:nic)
!      common /work/ az, hiz, hsz

  if (oinit .or. ofinal) then
     return
  end if

  if (ofirst) then
     ofirst = .false.
     call rewnml(ifpar, jfpar)
     read (ifpar, nmislt, iostat=istat)
     call cstnml(jfpar, 'fwater', 'nmislt', istat)
     write(jfpar, nmislt)
     rri = rhoo / rhoi
     rrs = rhoo / rhos
     tmi = dtds * si
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
              inrsbs(ij) = inrsbs(ij) &
                &        + ax(ij, k) * (hsx(ij, k) - hsz(ij, k)) &
                &          / rrs / ts
              if (hsx(ij, k) .le. 0.d0) then
                 fdd(ij) = fdd(ij) - ax(ij, k) * dsdx(ij, k) / ts
                 fdb(ij) = fdb(ij) - ax(ij, k) * dsbx(ij, k) / ts
                 dsdx(ij, k) = 0.d0
                 dsbx(ij, k) = 0.d0
              end if
           end if
           dhi = ts * rri * wsb(ij, k) / ax(ij, k)
           dhimax = eix(ij, k) / hfus
           if (dhi < dhimax) then 
              hix(ij, k) = hiz(ij, k) - dhi
              dei = (hix(ij, k) - hiz(ij, k)) * hfus
              eix(ij, k) = eix(ij, k) + dei
              tix(ij, k) = ti(eix(ij, k)/hix(ij, k), si)
              inrsbi(ij) = inrsbi(ij) &
                &        + ax(ij, k) * (hix(ij, k) - hiz(ij, k)) &
                &          / rri / ts
           else
!              write(0,*) '## dhi >= dhimax at ifwtr ##'
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
              wiadjs(ij) = wiadjs(ij) &
                &        - ax(ij, k) * (hiz(ij, k) - dhimax) &
                &          / rri / ts
              inrsbi(ij) = inrsbi(ij) &
                &        + ax(ij, k) * (-dhimax) &
                &          / rri / ts
              imraji(ij) = imraji(ij) &
                &        + ax(ij, k) * (hiz(ij, k) - dhimax) &
                &          / rri / ts
           end if
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
