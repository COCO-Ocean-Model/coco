module bgs2d
! --- information -----------------------------------------------------
!
!  Group communication (GATHER and SCATTER) routines for 2-dimensional
! variables.
!
!  HISTORY
!     '99.10.05  H.Hasumi
!     '07.04.23  H.Hasumi
!     '12.10.07  M.kurogi: for COCO5.0
!
! --------------------------------------------------------------------- 
 implicit none
 private
 public gs2dst, gather_2d, scatter_2d, scatter_2d_int, scatter_sfc

 integer, allocatable :: irstr(:), jrstr(:)
 integer, allocatable :: irend(:), jrend(:)
 integer, allocatable :: nsndc(:), nrcvc(:), ndisp(:)
 integer :: ngsnd,  nsrcv
 save irstr, jrstr, irend, jrend, nsndc, nrcvc, ndisp, &
   &     ngsnd, nsrcv

 real(8),  allocatable :: sndbuf(:, :), rcvbuf(:, :)
 integer, allocatable :: ircbuf(:, :)
 real(8),  allocatable :: dummy(:)
 integer, allocatable :: idumy(:)
 save sndbuf, rcvbuf, ircbuf, dummy, idumy

contains
 subroutine gs2dst
 use zocdim, only:  nx, ny, nxy, nxyg, nprocs,    &
  & inodes, jnodes, ijnode, igstr, jgstr, myrank

 implicit none
 integer ::  i, j, n

 allocate(irstr(0:nprocs-1))
 allocate(jrstr(0:nprocs-1))
 allocate(irend(0:nprocs-1))
 allocate(jrend(0:nprocs-1))
 allocate(nsndc(0:nprocs-1))
 allocate(nrcvc(0:nprocs-1))
 allocate(ndisp(0:nprocs-1))
 allocate(dummy(nxyg+nprocs))
 allocate(idumy(nxyg+nprocs))
 if (myrank .ge. ijnode) then
    allocate(sndbuf(1, 1))
    allocate(rcvbuf(1, 1))
    allocate(ircbuf(1, 1))
 else
    allocate(sndbuf(nx, ny))
    allocate(rcvbuf(nx, ny))
    allocate(ircbuf(nx, ny))
 end if

 sndbuf = 0.d0
 rcvbuf = 0.d0
 dummy  = 0.d0
 ircbuf = 0
 idumy  = 0

 do j = 0, jnodes-1
    do i = 0, inodes-1
       n = j * inodes + i
       irstr(n) = igstr + i * nx
       irend(n) = irstr(n) + nx - 1
       jrstr(n) = jgstr + j * ny
       jrend(n) = jrstr(n) + ny - 1
       nsndc(n) = nxy
       nrcvc(n) = nxy
    end do
 end do
 do n = ijnode, nprocs-1
    irstr(n) = 2
    irend(n) = 1
    jrstr(n) = 2
    jrend(n) = 1
    nsndc(n) = 1
    nrcvc(n) = 1
 end do
 ndisp(0) = 0
 do n = 1, nprocs-1
    ndisp(n) = ndisp(n-1) + nsndc(n-1)
 end do

 if (myrank .ge. ijnode) then
    ngsnd = 1
    nsrcv = 1
 else
    ngsnd = nxy
    nsrcv = nxy
 end if

 return 
 end subroutine gs2dst
! =====================================================================

 subroutine gather_2d(          &
     &                    qg,   &
     &                     q)
 use zocdim, only:  nx, ny, nxdim, nydim, nxgdim, nygdim, &
  &  nprocs, myrank, iroot, ijnode, istr, jstr, ierr, mpi_comm_ogcm


 implicit none
#include "mpif.h"
 real(8), intent (in ) ::     q(nxdim, nydim)
 real(8), intent (out) ::    qg(nxgdim, nygdim)
 integer ::  i, j, ijg, n


 if (myrank .lt. ijnode) then
    do j = 1, ny
       do i = 1, nx
          sndbuf(i, j) = q(i+istr-1, j+jstr-1)
       end do
    end do
 end if

 call mpi_gatherv(                                      &
  &                 sndbuf,  ngsnd,         mpi_real8,  &
  &                  dummy,  nrcvc,  ndisp, mpi_real8,  &
  &                  iroot, mpi_comm_ogcm, ierr)

 if (myrank .eq. iroot) then
    do n = 0, nprocs-1
       ijg = ndisp(n)
       do j = jrstr(n), jrend(n)
          do i = irstr(n), irend(n)
             ijg = ijg + 1
             qg(i, j) = dummy(ijg)
          end do
       end do
    end do
 end if

 return
 end subroutine gather_2d
! =====================================================================

 subroutine scatter_2d(         &
     &                      q,  &
     &                     qg)
 use zocdim, only:  nx, ny, nxdim, nydim, nxgdim, nygdim, &
  &  nprocs, myrank, iroot, ijnode, istr, jstr, ierr, mpi_comm_ogcm

 implicit none
#include "mpif.h"
 real(8), intent(out) ::      q(nxdim, nydim)
 real(8), intent(in ) ::    qg(nxgdim, nygdim)
 integer ::  i, j, ijg, n

 if (myrank .eq. iroot) then
    do n = 0, nprocs-1
       ijg = ndisp(n)
       do j = jrstr(n), jrend(n)
          do i = irstr(n), irend(n)
             ijg = ijg + 1
             dummy(ijg) = qg(i, j)
          end do
       end do
    end do
 end if

 call mpi_scatterv(                                     &
  &                   dummy,  nsndc,  ndisp, mpi_real8, &
  &                  rcvbuf,  nsrcv,         mpi_real8, &
  &                   iroot, mpi_comm_ogcm, ierr)

 if (myrank .lt. ijnode) then
    do j = 1, ny
       do i = 1, nx
          q(i+istr-1, j+jstr-1) = rcvbuf(i, j)
       end do
    end do
 end if

 return
 end subroutine scatter_2d
! =====================================================================

 subroutine scatter_2d_int(         &
     &                         iq,  &
     &                        iqg)
 use zocdim, only:  nx, ny, nxdim, nydim, nxgdim, nygdim, &
  &  nprocs, myrank, iroot, ijnode, istr, jstr, ierr, mpi_comm_ogcm

 implicit none
#include "mpif.h"
 integer, intent(out) ::   iq(nxdim, nydim)
 integer, intent(in ) ::  iqg(nxgdim, nygdim) 
 integer ::  i, j, ijg, n

 if (myrank .eq. iroot) then
    do n = 0, nprocs-1
       ijg = ndisp(n)
       do j = jrstr(n), jrend(n)
          do i = irstr(n), irend(n)
             ijg = ijg + 1
             idumy(ijg) = iqg(i, j)
          end do
       end do
    end do
 end if

 call mpi_scatterv(                                       &
  &                   idumy,  nsndc,  ndisp, mpi_integer, &
  &                  ircbuf,  nsrcv,         mpi_integer, &
  &                   iroot, mpi_comm_ogcm, ierr)

 if (myrank .lt. ijnode) then
    do j = 1, ny
       do i = 1, nx
          iq(i+istr-1, j+jstr-1) = ircbuf(i, j)
       end do
    end do
 end if

 return
 end subroutine scatter_2d_int
! =====================================================================

 subroutine scatter_sfc(            &
     &                         qs,  &
     &                        qsg)
 use zocdim, only:  nx, ny, igstr, igend, jgstr, jgend, &
  &  nprocs, myrank, iroot, ijnode, ierr, mpi_comm_ogcm

 implicit none
#include "mpif.h"
 real(8), intent(out) ::    qs(nx   , ny   )
 real(8), intent(in ) ::   qsg(igstr:igend, jgstr:jgend)
 integer ::  i, j, ijg, n

 if (myrank .eq. iroot) then
    do n = 0, nprocs-1
       ijg = ndisp(n)
       do j = jrstr(n), jrend(n)
          do i = irstr(n), irend(n)
             ijg = ijg + 1
             dummy(ijg) = qsg(i, j)
          end do
       end do
    end do
 end if

 call mpi_scatterv(                                     &
  &                   dummy,  nsndc,  ndisp, mpi_real8, &
  &                  rcvbuf,  nsrcv,         mpi_real8, &
  &                   iroot, mpi_comm_ogcm, ierr)

 if (myrank .lt. ijnode) then
    do j = 1, ny
       do i = 1, nx
          qs(i, j) = rcvbuf(i, j)
       end do
    end do
 end if

 return
 end subroutine scatter_sfc
end module bgs2d
