module ifhea

! --- information -----------------------------------------------------
!
!  HISTORY
!     '02.10.10  H.Hasumi: from MIROC3.1
!     '07.04.23  H.Hasumi
!     '07.09.25  H.Hasumi: multicategory sea ice
!     '07.10.03  H.Hasumi: 1-layer sea ice thermodynamics
!     '08.06.11  H.Hasumi: initial/final processing
!     '08.07.10  H.Hasumi: initial/final processing
!     '12.07.19  Y.Komuro: for COCO5.0
!
! ---------------------------------------------------------------------

  use zocdim, only: &
    & nxydim,  nzdim,  ntdim,    nic,    kstr, ijtstr, ijtend, &
    &  oinit, ofinal
  use zocgrd, only: &
    &     ds,   zbot,     ts
  use zocphy, only: &
    &   rhoo,    cpo,   dtds 

  implicit none
#include "coco.h"
  real(8), save :: swcnv1(nxydim)

  private

  public :: fiheat, putswc

contains

subroutine fiheat( &
  &                    ft, &
  &                   wao,    wio,    was,    wil, &
  &                     a,      t,     sh, &
  &                   qao,    qai,    qio,    qii,  swabs )
  use ufile

  real(8), intent(inout) ::     ft(nxydim, ntdim)
  real(8), intent(out)   ::    wao(nxydim)
  real(8), intent(out)   ::    wio(nxydim, nic),    was(nxydim, nic)
  real(8), intent(out)   ::    wil(nxydim, nic)
  real(8), intent(in)    ::      a(nxydim, 0:nic)
  real(8), intent(in)    ::    qio(nxydim, nic),    qai(nxydim, nic)
  real(8), intent(in)    ::    qii(nxydim, nic)
  real(8), intent(in)    ::    qao(nxydim),  swabs(nxydim)
  real(8), intent(in)    ::      t(nxydim, nzdim, ntdim)
  real(8), intent(in)    ::     sh(nxydim)

  logical, save :: ofirst = .true.

  real(8) ::   wfrz(nxydim),    wib(nxydim),   wilm(nxydim)

  real(8) ::   tdev
  integer ::     ij,      k
  integer ::  ifpar,  jfpar,  istat

  real(8), save ::  fbtab = 0.5d0

  namelist /nmbtab/ fbtab

  if (oinit .or. ofinal) then
     return
  end if

  if (ofirst) then
     ofirst = .false.
     READ_NAMELIST( nmbtab )
     !$acc enter data create(wfrz, wib, wilm)
  end if
  
  !$acc kernels default(present)
  do ij = ijtstr, ijtend
     tdev    = t(ij, kstr, 1) - dtds * t(ij, kstr, 2)
     wfrz(ij) = - rhoo * cpo * tdev * &
       &         (sh(ij) + zbot) * ds(kstr) / ts
     ft (ij, 1) = tdev * (sh(ij) + zbot) * ds(kstr) / ts
     wao(ij) = a(ij, 0) * wfrz(ij) &
       &     + qao(ij) - swcnv1(ij) * swabs(ij)
  end do
  !$acc end kernels
  !$acc kernels default(present)
  do ij = ijtstr, ijtend
     if ((a(ij, 0) .lt. 1.d0) .and. (wao(ij) .lt. 0.d0)) then
        wib(ij) = fbtab * wao(ij) / (1.d0 - a(ij, 0))
        wilm(ij) = (1.d0 - fbtab) * wao(ij) / (1.d0 - a(ij, 0))
        wao(ij) = 0.d0
     else
        wib(ij) = 0.d0
        wilm(ij) = 0.d0
     end if
  end do
  !$acc end kernels
  !$acc kernels default(present)
  do k = 1, nic
     do ij = ijtstr, ijtend
        was(ij, k) = qai(ij, k) - qii(ij, k)
        wio(ij, k) = a(ij, k) * (wfrz(ij) + wib(ij)) &
          &        + qio(ij, k)
        wil(ij, k) = a(ij, k) * wilm(ij)
     end do
  end do
  !$acc end kernels
  return

end subroutine fiheat

! =====================================================================

subroutine putswc( &
  &                swconv )

  real(8), intent(in) :: swconv(nxydim)
  integer :: ij
  logical, save :: ofirst = .true.
  if (ofirst) then
     ofirst = .false.
     !$acc enter data create(swcnv1)
  end if
  
  !$acc kernels default(present)
  do ij = 1, nxydim
     swcnv1(ij) = swconv(ij)
  end do
  !$acc end kernels
  return
end subroutine putswc

end module ifhea
