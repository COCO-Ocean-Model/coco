module bgsid
! --- information -----------------------------------------------------
!
!  Group communication (GATHER and SCATTER) routines for sea ice
! variables (multi-category thickness representation).
!
!  HISTORY
!     '07.09.25  H.Hasumi: from multi-category icedcoco3
!     '12.10.07  M.Kurogi: rewrite in F95 format
!
! ---------------------------------------------------------------------
 use zocdim
 implicit none
 private
 public gsidst, gather_id, scatter_id
 integer, allocatable :: irstr(:), jrstr(:)
 integer, allocatable :: irend(:), jrend(:)
 integer, allocatable :: nsndc(:), nrcvc(:), ndisp(:)
 integer :: ngsnd,  nsrcv
 save irstr, jrstr, irend, jrend, nsndc, nrcvc, ndisp, &
  &     ngsnd, nsrcv

 real(8),  allocatable :: sndbuf(:, :, :), rcvbuf(:, :, :)
 real(8),  allocatable :: dummy(:)
 save sndbuf, rcvbuf, dummy

contains
 subroutine gsidst
 use zocdim, only : nx, ny, nxy, nic, igstr, jgstr, &
  &  nprocs, myrank, ijnode
 implicit none
 integer :: i, j, k, n, ijkg

 allocate(irstr(0:nprocs-1))
 allocate(jrstr(0:nprocs-1))
 allocate(irend(0:nprocs-1))
 allocate(jrend(0:nprocs-1))
 allocate(nsndc(0:nprocs-1))
 allocate(nrcvc(0:nprocs-1))
 allocate(ndisp(0:nprocs-1))
 allocate(dummy(nxyg*nic+nprocs))

 if (myrank > ijnode) then
    allocate(sndbuf(1, 1, 1))
    allocate(rcvbuf(1, 1, 1))
 else
    allocate(sndbuf(nx, ny, nic))
    allocate(rcvbuf(nx, ny, nic))
 end if

 sndbuf = 0.d0
 rcvbuf = 0.d0
 dummy  = 0.d0

 do j = 0, jnodes-1
    do i = 0, inodes-1
       n = j * inodes + i
       irstr(n) = igstr + i * nx
       irend(n) = irstr(n) + nx - 1
       jrstr(n) = jgstr + j * ny
       jrend(n) = jrstr(n) + ny - 1
       nsndc(n) = nxy*nic
       nrcvc(n) = nxy*nic
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
    ngsnd = nxy*nic
    nsrcv = nxy*nic
 end if

 return
 end subroutine gsidst
! =====================================================================

  subroutine gather_id(         &
     &                    qg,  &
     &                     q)
 use zocdim, only : nic, ny, nic, nprocs, ijnode, myrank, iroot, istr, jstr, mpi_comm_ogcm
 implicit none
#include "mpif.h"
 real(8), intent(in ) ::    q(nxdim , nydim , 0:nic)
 real(8), intent(out) ::   qg(nxgdim, nygdim, 0:nic)

 integer :: i, j, k, n, ijkg

 if (myrank < ijnode) then
    do k = 1, nic
       do j = 1, ny
          do i = 1, nx
             sndbuf(i, j, k) = q(i+istr-1, j+jstr-1, k)
          end do
       end do
    end do
 end if

 call mpi_gatherv(                                    &
  &                 sndbuf,  ngsnd,         mpi_real8,  &
  &                  dummy,  nrcvc,  ndisp, mpi_real8,  &
  &                  iroot, mpi_comm_ogcm, ierr)

 if (myrank == iroot) then
    do n = 0, nprocs-1
       ijkg = ndisp(n)
       do k = 1, nic
          do j = jrstr(n), jrend(n)
             do i = irstr(n), irend(n)
                ijkg = ijkg + 1
                qg(i, j, k) = dummy(ijkg)
             end do
          end do
       end do
    end do
 end if

 return
 end subroutine gather_id
! =====================================================================

 subroutine scatter_id(         &
  &                      q,     &
  &                     qg)
 use zocdim, only : nx, ny, nic, nprocs, ijnode, myrank, iroot, istr, jstr, mpi_comm_ogcm
 implicit none
#include "mpif.h"
 real(8), intent(out) ::    q(nxdim , nydim , 0:nic)
 real(8), intent(in ) ::   qg(nxgdim, nygdim, 0:nic)

 integer :: i, j, k, n, ijkg

 if (myrank == iroot) then
    do n = 0, nprocs-1
       ijkg = ndisp(n)
       do k = 1, nic
          do j = jrstr(n), jrend(n)
             do i = irstr(n), irend(n)
                ijkg = ijkg + 1
                dummy(ijkg) = qg(i, j, k)
             end do
          end do
       end do
    end do
 end if

 call mpi_scatterv(                                     &
  &                   dummy,  nsndc,  ndisp, mpi_real8, &
  &                  rcvbuf,  nsrcv,         mpi_real8, &
  &                   iroot, mpi_comm_ogcm, ierr)

 if (myrank < ijnode) then
    do k = 1, nic
       do j = 1, ny
          do i = 1, nx
             q(i+istr-1, j+jstr-1, k) = rcvbuf(i, j, k)
          end do
       end do
    end do
 end if

 return
 end  subroutine scatter_id
end module bgsid
