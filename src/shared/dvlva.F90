module dvlva
! --- information -----------------------------------------------------
!
!  Velocity for momentum advection
!
!  HISTORY
!     '99.08.17  H.Hasumi
!     '01.01.25  H.Hasumi: estimate of WADV is corrected
!     '01.01.30  H.Hasumi: for partial step bottom topography
!     '01.02.13  H.Hasumi: for incorporating BBL model
!     '01.05.02  H.Hasumi
!     '07.04.23  H.Hasumi
!     '08.07.10  H.Hasumi: initial/final processing
!     '09.02.23  Y.Komuro: bug fix by Dr. Kurogi
!                          (loop length extended)
!     '12.07.21  M.kurogi: for COCO5.0
!
! ---------------------------------------------------------------------
  private
  public :: velvad
 
contains

  subroutine velvad(                                                           &
   &               uadv,   vadv,   wadv,                                       &
   &                  u,      v,      w)
  use zocdim, only :                                                           &
   &  nxydim, nzdim, kstr, kend, kz, ijtstr, ijtend,                           &
   &  le, ln, lnw, lne, lw, ls, lsw, lse,                                      &
   &  oinit, ofinal

  use zocgrd, only : dy, dym, rym, dzv
  use zocmsk, only : amskt, amskv
  implicit none

  real(8), intent(out) :: uadv(nxydim, nzdim),   vadv(nxydim, nzdim)
  real(8), intent(out) :: wadv(nxydim, nzdim)
  real(8), intent(in)  ::    u(nxydim, nzdim),      v(nxydim, nzdim)
  real(8), intent(in)  ::    w(nxydim, nzdim)

  integer ::  ij,      k

  if (oinit .or. ofinal) then
     return
  end if

  do k = 1, nzdim
     do ij = 1, nxydim
        uadv(ij, k) = 0.d0
        vadv(ij, k) = 0.d0
        wadv(ij, k) = 0.d0
     end do
  end do

  do k = kstr, kstr+kz-1
     do ij = ijtstr, ijtend
          wadv(ij, k) = (  (w(ij+lne, k) + w(ij+ln, k)) * dy(ij+ln)            &
   &                     + (w(ij+le , k) + w(ij   , k)) * dy(ij)) *            &
!  &                    0.25d0 * rym(ij) * amfvz(ij, k)
   &                    0.25d0 * rym(ij) * amskv(ij, k)
          vadv(ij, k) = (  v(ij+le, k) + v(ij+lw, k)                           &
   &                     + 2.d0 * v(ij, k)                                     &
   &                     + v(ij+lse, k) + v(ij+lsw, k)                         &
   &                     + 2.d0 * v(ij+ls, k)) *                               &
   &                    0.125d0 * amskt(ij, k) * amskt(ij+le, k)
          uadv(ij, k) = (  (u(ij+lnw, k) + u(ij+ln, k)) * dy(ij+ln)            &
   &                     + (u(ij+lsw, k) + u(ij+ls, k)) * dy(ij)               &
   &                     + (u(ij+lw , k) + u(ij   , k)) * dym(ij) *            &
   &                       2.d0                                                &
   &                    ) * 0.125 * rym(ij) *                                  &
   &                    amskt(ij, k) * amskt(ij+ln, k)
     end do
  end do


  do k = kstr+kz, kend
     do ij = ijtstr, ijtend
          wadv(ij, k) = (  (w(ij+lne, k) + w(ij+ln, k)) * dy(ij+ln)            &
   &                     + (w(ij+le , k) + w(ij   , k)) * dy(ij)) *            &
!  &                    0.25d0 * rym(ij) * amfvz(ij, k)
   &                    0.25d0 * rym(ij) * amskv(ij, k)
          vadv(ij, k) = (  v(ij+le , k) * dzv(ij+le , k)                       &
   &                     + v(ij+lw , k) * dzv(ij+lw , k)                       &
   &                     + v(ij    , k) * dzv(ij    , k) * 2.d0                &
   &                     + v(ij+lse, k) * dzv(ij+lse, k)                       &
   &                     + v(ij+lsw, k) * dzv(ij+lsw, k)                       &
   &                     + v(ij+ls , k) * dzv(ij+ls , k) * 2.d0) *             &
   &                    0.125d0 * amskt(ij, k) * amskt(ij+le, k)
          uadv(ij, k) = (  (  u(ij+lnw, k) * dzv(ij+lnw, k)                    &
   &                        + u(ij+ln , k) * dzv(ij+ln , k)                    &
   &                        + u(ij+lw , k) * dzv(ij+lw , k)                    &
   &                        + u(ij    , k) * dzv(ij    , k)) *                 &
   &                       dy(ij+ln)                                           &
   &                     + (  u(ij    , k) * dzv(ij    , k)                    &
   &                        + u(ij+lw , k) * dzv(ij+lw , k)                    &
   &                        + u(ij+ls , k) * dzv(ij+ls , k)                    &
   &                        + u(ij+lsw, k) * dzv(ij+lsw, k)) *                 &
   &                       dy(ij)                                              &
   &                     ) * rym(ij) * 0.125d0 *                               &
   &                     amskt(ij, k) * amskt(ij+ln, k)
     end do
  end do

  return
  end subroutine velvad
end module dvlva
