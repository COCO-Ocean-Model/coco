module bgs3d
! --- information -----------------------------------------------------
!
!  Group communication (GATHER and SCATTER) routines for 3-dimensional
! variables.
!
!  HISTORY
!     '99.10.05  H.Hasumi
!     '07.04.23  H.Hasumi
!     '12.10.08  M.kurogi: for COCO5.0
!
! ---------------------------------------------------------------------
use zocdim
implicit none
private
public :: gs3dst, gather_3d, gather_chk_sig, scatter_3d, scatter_bdy, &
  &       gs3dst_sig

 integer, allocatable :: irstr(:), jrstr(:)
 integer, allocatable :: irend(:), jrend(:)
 integer, allocatable :: nsndc(:), nrcvc(:), ndisp(:)
 integer  ngsnd,  nsrcv
 save irstr, jrstr, irend, jrend, nsndc, nrcvc, ndisp, &
  &   ngsnd, nsrcv

 real(8),  allocatable :: sndbuf(:, :, :), rcvbuf(:, :, :)
 real(8),  allocatable :: dummy(:)
      save sndbuf, rcvbuf, dummy

 integer, allocatable, save :: nrcvcs(:), ndisps(:)
 real*8, allocatable, save :: sndbus(:, :, :)
 real*8, allocatable, save :: dummys(:)
 integer, save :: ngsnds


contains
 subroutine gs3dst
 use zocdim, only : nprocs, nx, ny, nxyzg, kstr, kend, &
  &   myrank, ijnode, inodes, jnodes
 implicit none
 integer  ::   i, j, n

 allocate(irstr(0:nprocs-1))
 allocate(jrstr(0:nprocs-1))
 allocate(irend(0:nprocs-1))
 allocate(jrend(0:nprocs-1))
 allocate(nsndc(0:nprocs-1))
 allocate(nrcvc(0:nprocs-1))
 allocate(ndisp(0:nprocs-1))
 allocate(dummy(nxyzg+nprocs))
 if (myrank .ge. ijnode) then
    allocate(sndbuf(1, 1, 1))
    allocate(rcvbuf(1, 1, 1))
 else
    allocate(sndbuf(nx, ny, kstr:kend))
    allocate(rcvbuf(nx, ny, kstr:kend))
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
       nsndc(n) = nxyz
       nrcvc(n) = nxyz
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
    ngsnd = nxyz
    nsrcv = nxyz
 end if

 return
 end subroutine gs3dst
! =====================================================================
 subroutine gs3dst_sig( &
   &                    nsigmx )
 use zocdim, only : nprocs, nx, ny, nxg, nyg, &
   &                myrank, inodes, jnodes, ijnode
 implicit none
 integer, intent(in) :: nsigmx
 integer  ::   i, j, n
    
 allocate(nrcvcs(0:nprocs-1))
 allocate(ndisps(0:nprocs-1))
 allocate(dummys(nxg*nyg*nsigmx+nprocs))
 if (myrank .ge. ijnode) then
    allocate(sndbus(1, 1, 1))
 else
    allocate(sndbus(nx, ny, nsigmx))
 end if

 sndbus = 0.d0
 dummys = 0.d0

 do j = 0, jnodes-1
    do i = 0, inodes-1
       n = j * inodes + i
       nrcvcs(n) = nx*ny*nsigmx
    end do
 end do
 do n = ijnode, nprocs-1
    nrcvcs(n) = 1
 end do
 ndisps(0) = 0
 do n = 1, nprocs-1
    ndisps(n) = ndisps(n-1) + nrcvcs(n-1)
 end do

 if (myrank .ge. ijnode) then
    ngsnds = 1
 else
    ngsnds = nx*ny*nsigmx
 end if

 return
 end subroutine gs3dst_sig
! =====================================================================

 subroutine gather_3d(          &
     &                    qg,   &
     &                     q)
 use zocdim, only : nprocs, nx, ny, istr, jstr, kstr, kend,      &
  &   nxdim, nydim, nzdim, nxgdim, nygdim, myrank, ijnode, iroot
 implicit none
#include "mpif.h"
 real(8), intent(in ) ::    q(nxdim , nydim , nzdim)
 real(8), intent(out) ::  qg(nxgdim, nygdim, nzdim)
 integer  ::   i, j, ijkg, k,  n

 if (myrank < ijnode) then
    do k = kstr, kend
       do j = 1, ny
          do i = 1, nx
             sndbuf(i, j, k) = q(i+istr-1, j+jstr-1, k)
          end do
       end do
    end do
 end if

 call mpi_gatherv(                                     &
  &                 sndbuf,  ngsnd,         mpi_real8, &
  &                  dummy,  nrcvc,  ndisp, mpi_real8, &
  &                  iroot, mpi_comm_world, ierr)

 if (myrank == iroot) then
    do n = 0, nprocs-1
       ijkg = ndisp(n)
       do k = kstr, kend
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
 end subroutine gather_3d
! =====================================================================

 subroutine gather_chk_sig(     &
     &                  qchksg,  &
     &                   qchks, nsigmx)
 use zocdim, only : nprocs, nx, ny, &
  &   igstr, igend, jgstr, jgend, myrank, ijnode, iroot
 implicit none
#include "mpif.h"
 real(8), intent(in ) ::  qchks(nx         , ny         , nsigmx)
 real(8), intent(out) :: qchksg(igstr:igend, jgstr:jgend, nsigmx)
 integer, intent(in) :: nsigmx
 integer  ::   i, j, ijkg, k,  n

 if (myrank < ijnode) then
    do k = 1, nsigmx
       do j = 1, ny
          do i = 1, nx
             sndbus(i, j, k) = qchks(i, j, k)
          end do
       end do
    end do
 end if

 call mpi_gatherv(                                      &
  &                 sndbus, ngsnds,         mpi_real8,  &
  &                 dummys, nrcvcs, ndisps, mpi_real8,  &
  &                  iroot, mpi_comm_world, ierr)

 if (myrank == iroot) then
    do n = 0, nprocs-1
       ijkg = ndisps(n)
       do k = 1, nsigmx
          do j = jrstr(n), jrend(n)
             do i = irstr(n), irend(n)
                ijkg = ijkg + 1
                qchksg(i, j, k) = dummys(ijkg)
             end do
          end do
       end do
    end do
 end if

 return
 end subroutine gather_chk_sig
! =====================================================================

 subroutine scatter_3d(         &
     &                      q,  &
     &                     qg)
 use zocdim, only : nprocs, nx, ny, istr, jstr, kstr, kend,      &
  &   nxdim, nydim, nzdim, nxgdim, nygdim, myrank, ijnode, iroot
 implicit none
#include "mpif.h"
 real(8), intent(out) ::    q(nxdim , nydim , nzdim)
 real(8), intent(in ) ::  qg(nxgdim, nygdim, nzdim)
 integer  ::   i, j, ijkg, k,  n

 if (myrank == iroot) then
    do n = 0, nprocs-1
       ijkg = ndisp(n)
       do k = kstr, kend
          do j = jrstr(n), jrend(n)
             do i = irstr(n), irend(n)
                ijkg = ijkg + 1
                dummy(ijkg) = qg(i, j, k)
             end do
          end do
       end do
    end do
 end if

 call mpi_scatterv(                                      &
  &                   dummy,  nsndc,  ndisp, mpi_real8,  &
  &                  rcvbuf,  nsrcv,         mpi_real8,  &
  &                   iroot, mpi_comm_world, ierr)

 if (myrank < ijnode) then
    do k = kstr, kend
       do j = 1, ny
          do i = 1, nx
             q(i+istr-1, j+jstr-1, k) = rcvbuf(i, j, k)
          end do
       end do
    end do
 end if

 return
 end subroutine scatter_3d
! =====================================================================

 subroutine scatter_bdy(        &
     &                    qchk, &
     &                   qchkg) 
 use zocdim, only : nprocs, nx, ny, kstr, kend,         &
  &   igstr, igend, jgstr, jgend, myrank, ijnode, iroot
 implicit none
#include "mpif.h"
 real(8), intent(out) ::  qchk(nx         , ny         , kstr:kend)
 real(8), intent(in ) :: qchkg(igstr:igend, jgstr:jgend, kstr:kend)
 integer  ::   i, j, ijkg, k,  n

      if (myrank == iroot) then
         do n = 0, nprocs-1
            ijkg = ndisp(n)
            do k = kstr, kend
               do j = jrstr(n), jrend(n)
                  do i = irstr(n), irend(n)
                     ijkg = ijkg + 1
                     dummy(ijkg) = qchkg(i, j, k)
                  end do
               end do
            end do
         end do
      end if

      call mpi_scatterv(                                    &
     &                   dummy,  nsndc,  ndisp, mpi_real8,  &
     &                  rcvbuf,  nsrcv,         mpi_real8,  &
     &                   iroot, mpi_comm_world, ierr)

      if (myrank < ijnode) then
         do k = kstr, kend
            do j = 1, ny
               do i = 1, nx
                  qchk(i, j, k) = rcvbuf(i, j, k)
               end do
            end do
         end do
      end if

      return
 end  subroutine scatter_bdy
end module bgs3d
