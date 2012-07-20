module fbtav

! ---- information ----------------------------------------------------
!
!  HISTORY
!     '99.08.17  H.Hasumi
!     '07.04.23  H.Hasumi
!     '12.07.20  M.kurogi: for COCO5.0
!
! ---------------------------------------------------------------------

  implicit none
  private
  public btavst, btavad

  real(8), save ::   favr

contains
  subroutine btavst(                                                           &
   &                     ubtav,  vbtav)
  use zocdim, only : nxydim
  use zocgrd, only : ntss
  implicit none

  real(8), intent(out) ::  ubtav(nxydim),  vbtav(nxydim)
  integer :: ij

  do ij = 1, nxydim
     ubtav(ij) = 0.d0
     vbtav(ij) = 0.d0
  end do

  favr = 1.d0 / dble(ntss)

  return
  end subroutine btavst

! =====================================================================

  subroutine btavad(                                                           &
   &              ubtav,  vbtav,                                               &
   &                ubt,    vbt)
  use zocdim, only : nxydim
  implicit none
  real(8), intent(inout) ::  ubtav(nxydim),  vbtav(nxydim)
  real(8), intent(in)    ::    ubt(nxydim),    vbt(nxydim)
  integer :: ij

  do ij = 1, nxydim
     ubtav(ij) = ubtav(ij) + favr * ubt(ij)
     vbtav(ij) = vbtav(ij) + favr * vbt(ij)
  end do

  return
  end subroutine btavad
end module fbtav
