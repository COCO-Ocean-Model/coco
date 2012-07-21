module dwdns
! --- information -----------------------------------------------------
!
!  Estimate the vertical velocity at T-points
!
!  HISTORY
!     '03.04.23  H.Hasumi: from COCO3.4
!     '07.04.23  H.Hasumi
!     '08.07.10  H.Hasumi: initial/final processing
!     '12.07.21  M.kurogi: for COCO5.0
!
! ---------------------------------------------------------------------
  private
  public :: wdenst

contains
  subroutine wdenst(                                                           &
   &                  w,                                                       &
   &                 ux,     vx,                                               &
   &                 hx,     hz)

  use zocdim, only :                                                           &
   &  nxdim, nxydim, nzdim, kstr, kend, kz, ijtstr, ijtend,                    &
   &  le, ln, lw, ls, lsw, oinit, ofinal

  use zocgrd, only : ts, zbot, dzv, hxu, hyu, rx, ry, rxt, ryt, ds
  use zocmsk, only : amskt
  implicit none

  real(8), intent(in)  ::  ux(nxydim, nzdim),     vx(nxydim, nzdim)
  real(8), intent(in)  ::  hx(nxydim),            hz(nxydim)
  real(8), intent(out) ::   w(nxydim, nzdim)

  real(8) ::  ftx(nxydim, nzdim),    fty(nxydim, nzdim)
  real(8) :: dhdt(nxydim), rhzbot(nxydim)

  integer ::   ij,      k
  integer :: ijlw,   ijls,  ijlsw

  if (oinit .or. ofinal) then
     return
  end if

  do k = 1, nzdim
     do ij = 1, nxydim
            ftx(ij, k) = 0.d0
            fty(ij, k) = 0.d0
     end do
  end do

  do ij = 1, nxydim
     dhdt(ij) = (hx(ij) - hz(ij)) / ts
     rhzbot(ij) = 1.d0 / (hz(ij) + zbot)
  end do

  do k = kstr, kend
     do ij = ijtstr-nxdim-1, ijtend+nxdim+nxdim+1
        ijlw = ij + lw
        ijls = ij + ls
        ijlsw = ij + lsw
        ftx(ij, k) = (  ux(ijlw, k) * dzv(ijlw, k) *                           &
   &                    hyu(ijlw)                                              &
   &                  + ux(ijlsw, k) * dzv(ijlsw, k) *                         &
   &                    hyu(ijlsw)) *                                          &
   &                 amskt(ij, k) * amskt(ijlw, k)
        fty(ij, k) = (  vx(ijls, k) * dzv(ijls, k) *                           &
   &                    hxu(ijls)                                              &
   &                  + vx(ijlsw, k) * dzv(ijlsw, k) *                         &
   &                    hxu(ijlsw)) *                                          &
   &                 amskt(ij, k) * amskt(ijls, k)
     end do
  end do

  do ij = ijtstr-nxdim-1, ijtend+nxdim+1
     w(ij, kend) = - (  (ftx(ij+le, kend) - ftx(ij, kend)) *                   &
   &                    rx                                                     &
   &                  + (fty(ij+ln, kend) - fty(ij, kend)) *                   &
   &                    ry(ij)) * rxt(ij) * ryt(ij) * 0.5d0
  end do
  do k = kend-1, kstr+kz, -1
     do ij = ijtstr-nxdim-1, ijtend+nxdim+1
        w(ij, k) = w(ij, k+1)                                                  &
   &             - (  (ftx(ij+le, k) - ftx(ij, k)) * rx                        &
   &                + (fty(ij+ln, k) - fty(ij, k)) * ry(ij)) *                 &
   &               rxt(ij) * ryt(ij) * 0.5d0
     end do
  end do
  k = kstr+kz-1
  do ij = ijtstr-nxdim-1, ijtend+nxdim+1
     w(ij, k) = w(ij, k+1) * rhzbot(ij)                                        &
   &          - (  (  (ftx(ij+le, k) - ftx(ij, k)) * rx                        &
   &                + (fty(ij+ln, k) - fty(ij, k)) * ry(ij)) *                 &
   &               rxt(ij) * ryt(ij) * 0.5d0                                   &
   &             + dhdt(ij) * ds(k)) * rhzbot(ij)
  end do
  do k = kstr+kz-2, kstr+1, -1
     do ij = ijtstr-nxdim-1, ijtend+nxdim+1
        w(ij, k) = w(ij, k+1)                                                  &
   &             - (  (  (ftx(ij+le, k) - ftx(ij, k)) * rx                     &
   &                   + (fty(ij+ln, k) - fty(ij, k)) * ry(ij)) *              &
   &                  rxt(ij) * ryt(ij) * 0.5d0                                &
   &                + dhdt(ij) * ds(k)) * rhzbot(ij)
     end do
  end do
  do ij = ijtstr-nxdim-1, ijtend+nxdim+1
     w(ij, kstr) = 0.d0
  end do

  return
  end subroutine wdenst
end module dwdns
