module bstbc
! --- information -----------------------------------------------------
!
!  Set boudary conditions
!
!  HISTORY
!     '02.10.10  H.Hasumi: from COCO3.4
!     '07.04.23  H.Hasumi
!     '07.09.25  H.Hasumi: multi-category sea ice
!     '07.09.27  H.Hasumi: 1-layer sea ice thermodynamics
!     '07.10.24  H.Hasumi: STBBT2 restored (for SOM)
!     '12.07.21  M.kurogi: for COCO5.0
!
!----- For BBL
!  BBL variables are defined at two levels: one at KEND and the other
! at the NBOT (or NBOTV). This routine keeps two values consistent.
!
!  HISTORY
!     '01.02.09  H.Hasumi
!     '12.07.21  M.kurogi: for COCO5.0
!
! ---------------------------------------------------------------------
  private
#ifdef OPT_BBL
  public ::   stbctr, stbbgv, stbbuv, stbbvt, stbbgt, stbbtr, stbbt2
#else
  public ::   stbctr
#endif

contains

  subroutine stbctr(                                                           &
   &                       t,      r)
  use zocdim, only : nxdim, nydim, nzdim, ntdim, kstr
  use zocmsk, only : nbot
  implicit none
#include "mpif.h"
  real(8), intent(inout) ::  t(nxdim, nydim, nzdim, ntdim)
  real(8), intent(in)    ::  r(nxdim, nydim, nzdim)
  integer ::     i,      j,      k,     ij,      n 

  do n = 1, ntdim
     do k = 1, kstr-1
        do j = 1, nydim
           do i = 1, nxdim
              t(i, j, k, n) = t(i, j, kstr, n)
           end do
        end do
     end do
  end do

  do n = 1, ntdim
     do j = 1, nydim
        do i = 1, nxdim
               ij = i + (j - 1) * nxdim
               t(i, j, nbot(ij)+1, n) = t(i, j, nbot(ij), n)
        end do
     end do
  end do

  return
  end subroutine stbctr

#ifdef OPT_BBL
  subroutine stbbgv(                                                           &
   &                    gx,     gy)
  use zocdim, only : nxydim, nzdim, kstr, kend, ijvstr, ijvend
  use zocmsk, only : amskv, amskvb, nbotv
  implicit none
  real(8), intent(inout) ::   gx(nxydim, nzdim),     gy(nxydim, nzdim)
  integer ::   ij,      k

  do ij = ijvstr, ijvend
     k = nbotv(ij)
     gx(ij, k) = gx(ij, k) * (1.d0 - amskvb(ij)) + gx(ij, kend) * amskvb(ij) 
     gy(ij, k) = gy(ij, k) * (1.d0 - amskvb(ij)) + gy(ij, kend) * amskvb(ij)
  end do

  do k = kstr, kend
     do ij = ijvstr, ijvend
        gx(ij, k) = gx(ij, k) * amskv(ij, k)
        gy(ij, k) = gy(ij, k) * amskv(ij, k)
     end do
  end do

  return
  end subroutine stbbgv

! *********************************************************************

  subroutine stbbuv(                                                           &
   &                       u,      v)
  use zocdim, only : nxydim, nzdim, kend, ijvstr, ijvend
  use zocmsk, only : amskvb, nbotv
  implicit none

  real(8), intent(inout) ::  u(nxydim, nzdim),      v(nxydim, nzdim)
  integer ::  ij,      k

  do ij = ijvstr, ijvend
     k = nbotv(ij)
     u(ij, kend) = u(ij, k) * amskvb(ij)
     v(ij, kend) = v(ij, k) * amskvb(ij)
  end do

  return
  end subroutine stbbuv

! *********************************************************************

  subroutine stbbvt(                                                           &
   &                      u,      v)
  use zocdim, only : nxydim, nzdim, kend, ijvstr, ijvend
  use zocmsk, only : amskvb, nbotv
  implicit none

  real(8), intent(inout) ::    u(nxydim, nzdim),      v(nxydim, nzdim)
  integer ::    ij,      k

  do ij = ijvstr, ijvend
     k = nbotv(ij)
     u(ij, k) = u(ij, kend) * amskvb(ij) + u(ij, k) * (1.d0 - amskvb(ij))
     v(ij, k) = v(ij, kend) * amskvb(ij) + v(ij, k) * (1.d0 - amskvb(ij))
  end do

  return
  end subroutine stbbvt

! *********************************************************************

  subroutine stbbgt(                                                           &
   &                     adt,  diffz)
  use zocdim, only : nxydim, nzdim, ntdim, kstr, kend, ijtstr, ijtend
  use zocmsk, only : amskt, amsktb, nbot
  implicit none

  real(8), intent(inout) ::    adt(nxydim, nzdim, ntdim)
  real(8), intent(inout) ::  diffz(nxydim, nzdim)
  integer ::  ij,      k,      n

  do n = 1, ntdim
     do ij = ijtstr, ijtend
        k = nbot(ij)
        adt(ij, k, n) = adt(ij, k, n) * (1.d0 - amsktb(ij))                    &
   &                  + adt(ij, kend, n) * amsktb(ij)
     end do

     do k = kstr, kend
        do ij = ijtstr, ijtend
           adt(ij, k, n) = adt(ij, k, n) * amskt(ij, k)
           diffz(ij, k) = diffz(ij, k) * amskt(ij, k)
        end do
     end do
  enddo

  return
  end subroutine stbbgt

! *********************************************************************

  subroutine stbbtr(                                                           &
   &                      t)
  use zocdim, only : nxydim, nzdim, ntdim, kend, ijtstr, ijtend
  use zocmsk, only : nbot
  implicit none

  real(8), intent(inout) ::   t(nxydim, nzdim, ntdim)
  integer ::   ij,      k,      n

  do n = 1, ntdim
     do ij = ijtstr, ijtend
        k = nbot(ij)
        t(ij, kend, n) = t(ij, k, n)
     end do
  enddo

  return
  end subroutine stbbtr

! *********************************************************************

  subroutine stbbt2(                                                           &
   &                      t)
  use zocdim, only : nxydim, nzdim, ntdim, kend, ijtstr, ijtend
  use zocmsk, only : amsktb, nbot
  implicit none

  real(8), intent(inout) ::   t(nxydim, nzdim, ntdim)      
  integer ::   ij,      k,      n

  do n = 1, ntdim
     do ij = ijtstr, ijtend
        k = nbot(ij)
        t(ij, k, n) = t(ij, kend, n) * amsktb(ij)                              &
   &                + t(ij, k, n) * (1.d0 - amsktb(ij))
     end do
  enddo

  return
  end subroutine stbbt2
#endif
end module bstbc
