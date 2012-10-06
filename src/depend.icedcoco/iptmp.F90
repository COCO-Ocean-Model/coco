module iptmp

! --- information -----------------------------------------------------
!
!  HISTORY
!     '07.09.30  H.Hasumi
!     '08.06.11  H.Hasumi: initial/final processing
!     '08.07.10  H.Hasumi: initial/final processing
!     '12.07.18  Y.Komuro: for COCO5.0
!
! ---------------------------------------------------------------------

  use zocdim, only: &
    &  nxydim,    nic, ijtstr, ijtend, &
    &  oinit, ofinal
  use zocgrd, only: &
    &     ts
  use zocphy, only: &
    &   rhoi,    cpi,   dtds,   hfus 

  implicit none
  private

  public :: icetmp

contains

subroutine icetmp( &
  &                   eix, &
  &                   tix, &
  &                   qao,    qai,    qio,    qii, &
  &                    ax,    hix )
  use ufile

  real(8), intent(out)   ::    eix(nxydim, 0:nic)
  real(8), intent(inout) ::    tix(nxydim, 0:nic)
  real(8), intent(inout) ::    qao(nxydim)
  real(8), intent(inout) ::    qai(nxydim, nic),    qio(nxydim, nic)
  real(8), intent(inout) ::    qii(nxydim, nic)
  real(8), intent(in)    ::     ax(nxydim, 0:nic),    hix(nxydim, 0:nic)

  integer ::     ij,      k
  integer ::  ifpar,  jfpar,  istat

  logical, save :: ofirst = .true.

  real(8), save :: si = 5.0d0

  namelist /nmislt/ si

!===== Define statement function 
#include "zocite.F90"
!===== 

  if (oinit .or. ofinal) then
     return
  end if

  if (ofirst) then
     ofirst = .false.
     call rewnml(ifpar, jfpar)
     read (ifpar, nmislt, iostat=istat)
     call cstnml(jfpar, 'icetmp', 'nmislt', istat)
     write(jfpar, nmislt)
  end if

  do k = 0, nic ! zero-set for k = 0 (hi = 0 for k = 0)
     do ij = 1, nxydim
        eix(ij, k) = ei(tix(ij, k), si) * hix(ij, k)
     end do
  end do

  do k = 1, nic
     do ij = ijtstr, ijtend
        eix(ij, k) = eix(ij, k) &
          &        + ts * (qii(ij, k) - qio(ij, k)) / rhoi
        tix(ij, k) = ti(eix(ij, k)/hix(ij, k), si)
     end do
  end do

  do ij = ijtstr, ijtend
     qao(ij) = qao(ij) * ax(ij, 0)
  end do
  do k = 1, nic
     do ij = ijtstr, ijtend
        qai(ij, k) = qai(ij, k) * ax(ij, k)
        qio(ij, k) = qio(ij, k) * ax(ij, k)
        qii(ij, k) = qii(ij, k) * ax(ij, k)
     end do
  end do

  return
end subroutine icetmp

end module iptmp
