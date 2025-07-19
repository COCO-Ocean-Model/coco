module dwdns

! --- information -----------------------------------------------------
!
!  Estimate the vertical velocity at T-points
!
!  HISTORY
!     '03.04.23  H.Hasumi: from COCO3.4
!     '07.04.23  H.Hasumi
!
! ---------------------------------------------------------------------
   private
   public :: wdenst

contains
   subroutine wdenst(                                                           &
    &                  w,     hc,                                               &
    &                  u,      v,     ha)
   use zocdim, only :  &
    & nxdim, nxydim, nzdim, kstr, kend, kz, ijtstr, ijtend, &
    & le, ln, lw, ls, lsw, oinit, ofinal
   use zocgrd, only : ts, zbot, dzv, hxu, hyu, rx, ry, rxt, ryt, ds
   use zocmsk, only : amskt, amskv
   implicit none

   real(8),intent(inout) :: u(nxydim, nzdim),      v(nxydim, nzdim)
   real(8),intent(in)    ::       ha(nxydim)
   real(8),intent(out)   ::       hc(nxydim)
   real(8) ::       w(nxydim, nzdim)

   real(8) ::    ftx(nxydim, nzdim),    fty(nxydim, nzdim)
   real(8) ::    dhdt(nxydim), rhzbot(nxydim)

   integer ::     ij,      k
   integer ::   ijlw,   ijls,  ijlsw

   do k = 1, nzdim
      do ij = 1, nxydim
         ftx(ij, k) = 0.d0
         fty(ij, k) = 0.d0
      end do
   end do
   do ij = 1, nxydim
      dhdt(ij) = 0.d0
      hc(ij) = 0.d0
   end do

#ifdef OPT_BBL
   do k = kstr, kend
      do ij = 1, nxydim
         u(ij, k) = u(ij, k) * amskv(ij, k)
         v(ij, k) = v(ij, k) * amskv(ij, k)
      end do
   end do
#endif

   do k = kstr, kend
      do ij = ijtstr-nxdim-1, ijtend+nxdim+nxdim+1
         ijlw = ij + lw
         ijls = ij + ls
         ijlsw = ij + lsw
         ftx(ij, k) = (  u(ijlw, k) * dzv(ijlw, k) *                           &
    &                      hyu(ijlw)                                           &
    &                    + u(ijlsw, k) * dzv(ijlsw, k) *                       &
    &                      hyu(ijlsw)) *                                       &
    &                   amskt(ij, k) * amskt(ijlw, k)
            fty(ij, k) = (  v(ijls, k) * dzv(ijls, k) *                        &
    &                      hxu(ijls)                                           &
    &                    + v(ijlsw, k) * dzv(ijlsw, k) *                       &
    &                      hxu(ijlsw)) *                                       &
    &                   amskt(ij, k) * amskt(ijls, k)
      end do
   end do

   do k = kend, kstr, -1
      do ij = ijtstr-nxdim-1, ijtend+nxdim+1
         dhdt(ij) = dhdt(ij)                                                   &
    &               - (  (ftx(ij+le, k) - ftx(ij, k)) * rx                     &
    &                  + (fty(ij+ln, k) - fty(ij, k)) * ry(ij)) *              &
    &                 rxt(ij) * ryt(ij) * 0.5d0
      end do
   end do

   do ij = 1, nxydim
      hc(ij) = ha(ij) + dhdt(ij) * ts
      rhzbot(ij) = 1.d0 / (ha(ij) + zbot)
   end do

   do ij = ijtstr-nxdim-1, ijtend+nxdim+1
      w(ij, kend) = - (  (ftx(ij+le, kend) - ftx(ij, kend)) *                  &
    &                      rx                                                  &
    &                    + (fty(ij+ln, kend) - fty(ij, kend)) *                &
    &                      ry(ij)) * rxt(ij) * ryt(ij) * 0.5d0
   end do
   do k = kend-1, kstr+kz, -1
      do ij = ijtstr-nxdim-1, ijtend+nxdim+1
         w(ij, k) = w(ij, k+1)                                                 &
    &               - (  (ftx(ij+le, k) - ftx(ij, k)) * rx                     &
    &                  + (fty(ij+ln, k) - fty(ij, k)) * ry(ij)) *              &
    &                 rxt(ij) * ryt(ij) * 0.5d0
      end do
   end do
   k = kstr+kz-1
   do ij = ijtstr-nxdim-1, ijtend+nxdim+1
      w(ij, k) = w(ij, k+1) * rhzbot(ij)                                       &
    &            - (  (  (ftx(ij+le, k) - ftx(ij, k)) * rx                     &
    &                  + (fty(ij+ln, k) - fty(ij, k)) * ry(ij)) *              &
    &                 rxt(ij) * ryt(ij) * 0.5d0                                &
    &               + dhdt(ij) * ds(k)) * rhzbot(ij)
   end do
   do k = kstr+kz-2, kstr+1, -1
      do ij = ijtstr-nxdim-1, ijtend+nxdim+1
         w(ij, k) = w(ij, k+1)                                                 &
    &               - (  (  (ftx(ij+le, k) - ftx(ij, k)) * rx                  &
    &                     + (fty(ij+ln, k) - fty(ij, k)) * ry(ij)) *           &
    &                    rxt(ij) * ryt(ij) * 0.5d0                             &
    &                  + dhdt(ij) * ds(k)) * rhzbot(ij)
      end do
   end do
   do ij = ijtstr-nxdim-1, ijtend+nxdim+1
      w(ij, kstr) = 0.d0
   end do

   return
   end subroutine wdenst
end module dwdns
