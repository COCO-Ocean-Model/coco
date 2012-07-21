module dvlta
! --- information -----------------------------------------------------
!
!  Velocity for tracer advection
!
!  HISTORY
!     '99.08.16  H.Hasumi
!     '01.01.30  H.Hasumi: for partial step bottom topography
!     '07.04.23  H.Hasumi
!     '08.07.10  H.Hasumi: initial/final processing
!     '09.02.23  Y.Komuro: bug fix by Dr. Kurogi
!                          (initialize UADV/VADV)
!     '12.07.21  M.kurogi: for COCO5.0
!
! ---------------------------------------------------------------------
private
public :: veltad

contains
  subroutine veltad(                                      &
   &               uadv,   vadv,                          &
   &                ubt,    vbt,      u,      v)
  use zocdim, only : nxydim, nzdim, kstr, kend, ijvstr, ijvend, oinit, ofinal
  use zocgrd, only : dzv, rdepv
  use zocmsk, only : amskv
  implicit none

  real(8), intent(out) ::  uadv(nxydim, nzdim),   vadv(nxydim, nzdim)
  real(8), intent(in)  ::     u(nxydim, nzdim),      v(nxydim, nzdim)
  real(8), intent(in)  ::   ubt(nxydim),    vbt(nxydim)

  real(8) ::  uavr(nxydim),   vavr(nxydim)
  real(8) ::  ubar(nxydim),   vbar(nxydim)
  integer ::  ij,      k

  if (oinit .or. ofinal) then
     return
  end if

  do k = 1, nzdim
     do ij = 1, nxydim
        uadv(ij, k) = 0.0d0
        vadv(ij, k) = 0.0d0
     end do
  end do

  do ij = ijvstr, ijvend
     uavr(ij) = u(ij, kstr) * dzv(ij, kstr) * amskv(ij, kstr)
     vavr(ij) = v(ij, kstr) * dzv(ij, kstr) * amskv(ij, kstr)
  end do
  do k = kstr+1, kend
     do ij = ijvstr, ijvend
            uavr(ij) = uavr(ij) + u(ij, k) * dzv(ij, k) * amskv(ij, k)
            vavr(ij) = vavr(ij) + v(ij, k) * dzv(ij, k) * amskv(ij, k)
     end do
  end do
  do ij = ijvstr, ijvend
     uavr(ij) = uavr(ij) * rdepv(ij)
     vavr(ij) = vavr(ij) * rdepv(ij)
     ubar(ij) = ubt(ij) * rdepv(ij)
     vbar(ij) = vbt(ij) * rdepv(ij)
  end do

  do k = kstr, kend
     do ij = ijvstr, ijvend
            uadv(ij, k) = (ubar(ij) + (u(ij, k) - uavr(ij))) * amskv(ij, k)
            vadv(ij, k) = (vbar(ij) + (v(ij, k) - vavr(ij))) * amskv(ij, k)
     end do
  end do

  return
  end   subroutine veltad
end module dvlta
