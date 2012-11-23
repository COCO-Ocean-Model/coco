module bstbc

! --- information -----------------------------------------------------
!
!  Set boudary conditions
!  Note: For parallel processing, these routines should be replaced by
!        ones for inter-node data transfer.
!
!  HISTORY
!     '99.08.10  H.Hasumi
!     '02.05.29  H.Nakano: tracer dimension
!     '02.06.02  H.Hasumi: combine parallel and nonparallel
!     '07.04.23
!     '07.10.24  H.Hasumi: STBBT2 restored (for SOM)
!     '12.11.23  H.Tatebe: for COCO5.0 in F90
! ---------------------------------------------------------------------

  implicit none
  private

  public  ::  stbctr
#ifdef OPT_BBL
  public  ::  stbbgv,  stbbuv, stbbvt, stbbgt, stbbtr, stbbt2
#endif

contains

  subroutine stbctr( t,  r )

    use zocdim,  only  :                                              &
        &         nxdim,  nydim,  nzdim,  ntdim,                      &
        &          kstr
    use zocmsk,  only  :   nbot

    implicit none

#include "mpif.h"

    real(8),        intent(inout)  ::  t(1:nxdim,1:nydim,1:nzdim,1:ntdim)
    real(8),        intent(in)     ::  r(1:nxdim,1:nydim,1:nzdim)

    integer(4)  ::      i,      j,      k
    integer(4)  ::     ij,      n 

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

  end subroutine stbctr

#ifdef OPT_BBL

! --- information -----------------------------------------------------
!
!  BBL variables are defined at two levels: one at KEND and the other
! at the NBOT (or NBOTV). This routine keeps two values consistent.
!
!  HISTORY
!     '01.02.09  H.Hasumi
!
! ---------------------------------------------------------------------

  subroutine stbbgv( gx, gy )

    use zocdim,  only  :                                              &
        &        nxydim,  nzdim,  kstr,   kend,                       &
        &        ijvstr, ijvend
    use zocmsk,  only  :   nbotv,  amskv, amskvb
    
    implicit none

    real(8),        intent(inout)  ::  gx(1:nxydim,1:nzdim)
    real(8),        intent(inout)  ::  gy(1:nxydim,1:nzdim)

    integer(4)  ::   ij,  k

    do ij = ijvstr, ijvend

       k = nbotv(ij)
       gx(ij, k) = gx(ij, k   ) * (1.d0 - amskvb(ij))                 &
    &            + gx(ij, kend) * amskvb(ij)
       gy(ij, k) = gy(ij, k   ) * (1.d0 - amskvb(ij))                 &
    &            + gy(ij, kend) * amskvb(ij)

    end do

    do k = kstr, kend
       do ij = ijvstr, ijvend

          gx(ij, k) = gx(ij, k) * amskv(ij, k)
          gy(ij, k) = gy(ij, k) * amskv(ij, k)

       end do
    end do
    
  end subroutine stbbgv

! *********************************************************************

  subroutine stbbuv( u, v )


    use zocdim,  only  :                                              &
        &        nxydim,  nzdim,   kend,                              & 
        &        ijvstr, ijvend
    use zocmsk,  only  :   nbotv,  amskvb
    
    implicit none

    real(8),        intent(inout)  ::  u(1:nxydim,1:nzdim)
    real(8),        intent(inout)  ::  v(1:nxydim,1:nzdim)

    integer(4)  ::   ij,  k

    do ij = ijvstr, ijvend

       k = nbotv(ij)
       u(ij, kend) = u(ij, k) * amskvb(ij)
       v(ij, kend) = v(ij, k) * amskvb(ij)

    end do

  end subroutine stbbuv

! *********************************************************************

  subroutine stbbvt( u, v )

    use zocdim,  only  :                                              &
        &        nxydim,  nzdim,   kend,                              & 
        &        ijvstr, ijvend
    use zocmsk,  only  :   nbotv,  amskvb
    
    implicit none

    real(8),        intent(inout)  ::  u(1:nxydim,1:nzdim)
    real(8),        intent(inout)  ::  v(1:nxydim,1:nzdim)

    integer(4)  ::   ij,  k

    do ij = ijvstr, ijvend

       k = nbotv(ij)
       u(ij, k) = u(ij, kend) * amskvb(ij)                            &
    &           + u(ij, k   ) * (1.d0 - amskvb(ij))
       v(ij, k) = v(ij, kend) * amskvb(ij)                            &
    &           + v(ij, k   ) * (1.d0 - amskvb(ij))

    end do

  end subroutine stbbvt

! *********************************************************************

  subroutine stbbgt( adt,  diffz )

    use zocdim,  only  :                                              &
        &        nxydim,  nzdim,   ntdim,   kstr,   kend,             & 
        &        ijtstr, ijtend
    use zocmsk,  only  :   nbot,  amskt,  amsktb
    
    implicit none

    real(8),        intent(inout)  ::  adt (1:nxydim,1:nzdim,1:ntdim)
    real(8),        intent(inout)  :: diffz(1:nxydim,1:nzdim)

    integer(4)  ::   ij,  k,   n

    do n = 1, ntdim

       do ij = ijtstr, ijtend
          k = nbot(ij)
          adt(ij, k, n) = adt(ij, k   , n) * (1.d0 - amsktb(ij))      &
    &                   + adt(ij, kend, n) * amsktb(ij)
       end do

       do k = kstr, kend
          do ij = ijtstr, ijtend
             adt(ij, k, n) = adt(ij, k, n) * amskt(ij, k)
             diffz(ij, k) = diffz(ij, k) * amskt(ij, k)
          end do
       end do

    end do

  end subroutine stbbgt

! *********************************************************************

  subroutine stbbtr( t )

    use zocdim,  only  :                                              &
        &        nxydim,  nzdim,   ntdim,   kend,                     & 
        &        ijtstr, ijtend
    use zocmsk,  only  :   nbot
    
    implicit none

    real(8),        intent(inout)  ::  t(1:nxydim,1:nzdim,1:ntdim)

    integer(4)  ::   ij,  k,   n

    do n = 1, ntdim
       do ij = ijtstr, ijtend

          k = nbot(ij)
          t(ij, kend, n) = t(ij, k, n)

       end do
    end do

  end subroutine stbbtr

! *********************************************************************

  subroutine stbbt2( t )

    use zocdim,  only  :                                              &
        &        nxydim,  nzdim,   ntdim,   kend,                     & 
        &        ijtstr, ijtend
    use zocmsk,  only  :   nbot,  amsktb
    
    implicit none

    real(8),        intent(inout)  ::  t(1:nxydim,1:nzdim,1:ntdim)

    integer(4)  ::   ij,  k,   n

    do n = 1, ntdim
       do ij = ijtstr, ijtend

          k = nbot(ij)
          t(ij, k, n) = t(ij, kend, n) * amsktb(ij)                   &
    &                 + t(ij, k   , n) * (1.d0 - amsktb(ij))

       end do
    end do

  end subroutine stbbt2

#endif

end module bstbc
