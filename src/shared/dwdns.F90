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
  
  subroutine wdenst(                    &
   &                  w,                &
   &                 ux,     vx,        &
   &                 hx,     hz)

    use zocdim, only :                                            &
         &  nxdim, nxydim, nzdim, kstr, kend, kz, ijtstr, ijtend, &
         &  le, ln, lw, ls, lsw, oinit, ofinal
    use zocgrd, only : ts, zbot, dzv, hxu, hyu, rx, ry, rxt, ryt, ds
    use zocmsk, only : amskt, amskv

    implicit none

    real(8), intent(inout) ::   w(nxydim, nzdim)
    real(8), intent(inout) ::  ux(nxydim, nzdim)
    real(8), intent(inout) ::  vx(nxydim, nzdim)
    real(8), intent(inout) ::  hx(nxydim)
    real(8), intent(in)    ::  hz(nxydim)

    real(8) ::  ftx(nxydim, nzdim),    fty(nxydim, nzdim)
    real(8) :: dhdt(nxydim), rhzbot(nxydim)

    integer ::   ij,      k
    integer :: ijlw,   ijls,  ijlsw

#ifndef OPT_OFFLINE  
    if (oinit .or. ofinal) then
  !$acc enter data create(ftx,fty,dhdt,rhzbot)
       return
    end if
#endif  

  !$acc kernels default(present)
    do k = 1, nzdim
       do ij = 1, nxydim
          ftx(ij, k) = 0.d0
          fty(ij, k) = 0.d0
       end do
    end do

#ifndef OPT_OFFLINE
    do ij = 1, nxydim
       dhdt  (ij) = ( hx(ij) - hz(ij) ) / ts
       rhzbot(ij) = 1.d0 / (hz(ij) + zbot)
    end do
#else
    do ij = 1, nxydim
       dhdt(ij) = 0.d0
       hx  (ij) = 0.d0
    end do
#ifdef OPT_BBL
    do k = kstr, kend
       do ij = 1, nxydim
          ux(ij, k) = ux(ij, k) * amskv(ij, k)
          vx(ij, k) = vx(ij, k) * amskv(ij, k)
       end do
    end do
#endif
#endif

    do k = kstr, kend
       do ij = ijtstr-nxdim-1, ijtend+nxdim+nxdim+1
          ijlw = ij + lw
          ijls = ij + ls
          ijlsw = ij + lsw
          ftx(ij, k) = (  ux(ijlw, k) * dzv(ijlw, k) *       &
     &                    hyu(ijlw)                          &
     &                  + ux(ijlsw, k) * dzv(ijlsw, k) *     &
     &                    hyu(ijlsw)                         &
     &               ) * amskt(ij, k) * amskt(ijlw, k)
          fty(ij, k) = (  vx(ijls, k) * dzv(ijls, k) *       &
     &                    hxu(ijls)                          &
     &                  + vx(ijlsw, k) * dzv(ijlsw, k) *     &
     &                    hxu(ijlsw)                         &
     &               ) * amskt(ij, k) * amskt(ijls, k)
       end do
    end do

#ifdef OPT_OFFLINE    
    do k = kend, kstr, -1
       do ij = ijtstr-nxdim-1, ijtend+nxdim+1
          dhdt(ij) = dhdt(ij)                            &
     &             - (  (ftx(ij+le, k) - ftx(ij, k))     &
     &                * rx                               &
     &                + (fty(ij+ln, k) - fty(ij, k))     &
     &                * ry(ij) ) * rxt(ij) * ryt(ij) * 0.5d0
       end do
    end do
!---- correction of dhdt by fw
!    do ij = ijtstr-nxdim-1, ijtend+nxdim+1
!       dhdt(ij) = dhdt(ij) - fw(ij)
!    end do
    do ij = 1, nxydim
       hx(ij) = hz(ij) + dhdt(ij) * ts
       rhzbot(ij) = 1.d0 / (hz(ij) + zbot)
    end do
#endif    

    do ij = ijtstr-nxdim-1, ijtend+nxdim+1
       w(ij, kend) = - (  (ftx(ij+le, kend) - ftx(ij, kend))    &
     &                  * rx                                    &
     &                  + (fty(ij+ln, kend) - fty(ij, kend))    &
     &                  * ry(ij) ) * rxt(ij) * ryt(ij) * 0.5d0
    end do

    do k = kend-1, kstr+kz, -1
       do ij = ijtstr-nxdim-1, ijtend+nxdim+1
          w(ij, k) = w(ij, k+1)                                 &
     &              - (  (ftx(ij+le, k) - ftx(ij, k))           &
     &                 * rx                                     &
     &                 + (fty(ij+ln, k) - fty(ij, k))           &
     &                 * ry(ij) ) * rxt(ij) * ryt(ij) * 0.5d0
       end do
    end do

    k = kstr+kz-1
    do ij = ijtstr-nxdim-1, ijtend+nxdim+1
       w(ij, k) = w(ij, k+1) * rhzbot(ij)                       &
     &          - (  (  (ftx(ij+le, k) - ftx(ij, k))            &
     &                * rx                                      &
     &                + (fty(ij+ln, k) - fty(ij, k))            &
     &                * ry(ij) ) * rxt(ij) * ryt(ij) * 0.5d0    &
     &          + dhdt(ij) * ds(k) ) * rhzbot(ij)
    end do

    do k = kstr+kz-2, kstr+1, -1
       do ij = ijtstr-nxdim-1, ijtend+nxdim+1
          w(ij, k) = w(ij, k+1)                                 &
     &             - (  (  (ftx(ij+le, k) - ftx(ij, k))         &
     &                   * rx                                   &
     &                   + (fty(ij+ln, k) - fty(ij, k))         &
     &                   * ry(ij) ) * rxt(ij) * ryt(ij) * 0.5d0 &
     &              + dhdt(ij) * ds(k)) * rhzbot(ij)
       end do
    end do

    do ij = ijtstr-nxdim-1, ijtend+nxdim+1
       w(ij, kstr) = 0.d0
    end do
    
  !$acc end kernels
    return

  end subroutine wdenst
  
end module dwdns
