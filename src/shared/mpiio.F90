#ifdef OPT_IO_COCOMPI
! --- information -----------------------------------------------------
!
!  HISTORY
!     '15.04.07  M.Kurogi: MPI-IO
!
! ---------------------------------------------------------------------
module mpiio
  private
  public :: info_seq, mpi_write_header, mpi_read_chead,  &
       & mpi_read_sfc, mpi_read_bdy, mpi_read_root,      &
       & mpi_read_2d_intx, mpi_read_2d_dimx,             &
       & mpi_read_3d_dimx, mpi_read_2d, mpi_read_id,     &
       & mpi_read_3d, mpi_write_2d, mpi_write_id,        &
       & mpi_write_3d, reverse_real4, reverse_real8,     &
       & reverse_int4, mpi_read_direct,                  &
       & mpi_read_root_int, mpi_read_root_int_sgl,       &
       & mpi_iseof, mpi_read_root_char,                  &
       & info_seq8, reverse_int8
contains

  subroutine info_seq(fh, offset, nsize0)
  use zocdim
  implicit none
#include "mpif.h"
  integer :: nsize0, nsize
  integer :: fh
  integer   (kind = mpi_offset_kind) :: offset

  nsize=nsize0
  call reverse_int4(nsize) !swap endian

  call mpi_file_set_view(                       &
   &     fh, offset,                            &
   &     mpi_integer4, mpi_integer4,"native",   &
   &     mpi_info_null, ierr)
  if (myrank .eq. iroot) then
     call mpi_file_write(                       &
   &     fh, nsize, 1, mpi_integer4,            &
   &     mpi_status_ignore, ierr)
  end if
  offset = offset + 4

  end subroutine info_seq

  subroutine info_seq8(fh, offset, nsize0)
  use zocdim
  implicit none
#include "mpif.h"
  integer(8) :: nsize0, nsize
  integer :: fh
  integer   (kind = mpi_offset_kind) :: offset

  nsize=nsize0
  call reverse_int8(nsize) !swap endian

  call mpi_file_set_view(                       &
   &     fh, offset,                            &
   &     mpi_integer8, mpi_integer8,"native",   &
   &     mpi_info_null, ierr)
  if (myrank .eq. iroot) then
     call mpi_file_write(                       &
   &     fh, nsize, 1, mpi_integer8,            &
   &     mpi_status_ignore, ierr)
  end if
  offset = offset + 8

  end subroutine info_seq8

  subroutine mpi_write_header(chead, fh, offset)
  use zocdim
  implicit none
#include "mpif.h"

  character :: chead(64)*16
  integer :: fh
  integer (kind = mpi_offset_kind) :: offset
  integer :: nsize
  integer(8) :: nsize2
  nsize=64*16
  nsize2 = nsize

#ifdef OPT_IO_SEQUENTIAL
  call info_seq(fh, offset, nsize)
#elif defined(OPT_IO_SEQUENTIAL_H8)
  call info_seq8(fh, offset, nsize2)
#endif
!========= header
  call mpi_file_set_view(                      &
   &     fh, offset,                           &
   &     mpi_character,mpi_character,"native", &
   &     mpi_info_null,ierr)

  if (myrank .eq. iroot) then
     call mpi_file_write(                      &
     &     fh,chead,  nsize, mpi_character,    &
     &     mpi_status_ignore, ierr)
  end if
  offset = offset + nsize
#ifdef OPT_IO_SEQUENTIAL
  call info_seq(fh, offset, nsize)
#elif defined(OPT_IO_SEQUENTIAL_H8)
  call info_seq8(fh, offset, nsize2)
#endif

  return
  end subroutine mpi_write_header

!===============================================================
!     REVERSE_REAL4
!===============================================================

  subroutine reverse_real4(real4)
  implicit none

  real(4) :: real4, value
  integer(1) :: reverse(4), tmpval
  equivalence(value, reverse)

#ifdef OPT_IO_BYTESWAP
  reverse=0
  value = real4

  tmpval = reverse(1)
  reverse(1) = reverse(4)
  reverse(4) = tmpval
  tmpval = reverse(2)
  reverse(2) = reverse(3)
  reverse(3) = tmpval

  real4 = value
#endif
  return
  end subroutine reverse_real4

!===============================================================
!     REVERSE_REAL8
!===============================================================

  subroutine reverse_real8(real8)
  implicit none

  real(8) :: real8, value
  integer(1) :: reverse(8), tmpval
  equivalence(value, reverse)

#ifdef OPT_IO_BYTESWAP
  reverse=0
  value = real8

  tmpval = reverse(1)
  reverse(1) = reverse(8)
  reverse(8) = tmpval
  tmpval = reverse(2)
  reverse(2) = reverse(7)
  reverse(7) = tmpval
  tmpval = reverse(3)
  reverse(3) = reverse(6)
  reverse(6) = tmpval
  tmpval = reverse(4)
  reverse(4) = reverse(5)
  reverse(5) = tmpval

  real8 = value
#endif
  return
  end subroutine reverse_real8

!===============================================================
!     REVERSE_INT4
!===============================================================

  subroutine reverse_int4(int4)
  implicit none

  integer :: int4, value
  integer(1) :: reverse(4), tmpval
  equivalence(value, reverse)
#ifdef OPT_IO_BYTESWAP
  reverse=0
  value = int4

  tmpval = reverse(1)
  reverse(1) = reverse(4)
  reverse(4) = tmpval
  tmpval = reverse(2)
  reverse(2) = reverse(3)
  reverse(3) = tmpval

  int4 = value
#endif
  return
  end subroutine reverse_int4

!===============================================================
!     REVERSE_INT8
!===============================================================

  subroutine reverse_int8(int8)
  implicit none

  integer(8) :: int8, value
  integer(1) :: reverse(8), tmpval
  equivalence(value, reverse)

#ifdef OPT_IO_BYTESWAP
  reverse=0
  value = int8

  tmpval = reverse(1)
  reverse(1) = reverse(8)
  reverse(8) = tmpval
  tmpval = reverse(2)
  reverse(2) = reverse(7)
  reverse(7) = tmpval
  tmpval = reverse(3)
  reverse(3) = reverse(6)
  reverse(6) = tmpval
  tmpval = reverse(4)
  reverse(4) = reverse(5)
  reverse(5) = tmpval

  int8 = value
#endif
  return
  end subroutine reverse_int8

  subroutine mpi_read_chead(chead, fh, disp, icread)
  use zocdim
  implicit none
#include "mpif.h"
  character :: chead(64)*16
  integer :: fh
  integer (kind = mpi_offset_kind):: disp
  integer :: mpistat(mpi_status_size)
  integer :: icread
  integer :: ifpar, jfpar

#ifdef OPT_IO_SEQUENTIAL
  disp=disp+4
#elif defined(OPT_IO_SEQUENTIAL_H8)
  disp=disp+8
#endif
  call mpi_file_set_view(                            &
       &     fh, disp,                               &
       &     mpi_character, mpi_character,"native",  &
       &     mpi_info_null,ierr)

  if (myrank .eq. iroot) then
     call mpi_file_read(                            &
          &    fh, chead, 1024,                     &
          &    mpi_character, mpistat, ierr)

     call mpi_get_count(mpistat,mpi_character, icread,ierr)
  end if
  call mpi_bcast(chead, 1024, mpi_character,  &
       &                  iroot, mpi_comm_ogcm, ierr)
  call mpi_bcast(icread, 1, mpi_integer4,     &
       &                  iroot, mpi_comm_ogcm, ierr)
  disp=disp+ 1024
#ifdef OPT_IO_SEQUENTIAL
  disp=disp + 4
#elif defined(OPT_IO_SEQUENTIAL_H8)
  disp=disp + 8
#endif
  if(icread .ne. 1024) then
     disp=disp-1024
#ifdef OPT_IO_SEQUENTIAL
     disp=disp-8
#elif defined(OPT_IO_SEQUENTIAL_H8)
     disp=disp-16
#endif
  end if

  return
  end subroutine mpi_read_chead



  subroutine mpi_read_sfc(buf, fh, disp)
  use zocdim
  implicit none

#include "mpif.h"

  real(8) ::  buf(nx, ny)
  integer :: fh
  integer (kind = mpi_offset_kind):: disp
  integer ::  i, j, k
  integer :: ifile
  integer :: istart(2), igsize(2), isize(2)
  integer :: ifpar, jfpar

#ifdef OPT_IO_SEQUENTIAL
  disp=disp+4
#elif defined(OPT_IO_SEQUENTIAL_H8)
  disp=disp+8
#endif
  istart=(/irank*nx, jrank*ny/)
  igsize=(/nxg, nyg/)
  isize =(/nx , ny /)

  call mpi_type_create_subarray(     &
       &   2, igsize, isize, istart, &
       &   mpi_order_fortran,        &
       &   mpi_real8, ifile,  ierr)
  call mpi_type_commit(ifile, ierr)

  call mpi_file_set_view(           &
       &  fh, disp,                 &
       &  mpi_real8,ifile,"native", &
       &  mpi_info_null,ierr)

  call mpi_file_read_all(        &
       &  fh, buf, nx*ny,        &
       &   mpi_real8, mpi_status_ignore, ierr)


  do j=1,ny
     do i=1,nx
        call reverse_real8(buf(i,j))
     end do
  end do
  disp=disp+ nxg*nyg*8
#ifdef OPT_IO_SEQUENTIAL
  disp=disp+4
#elif defined(OPT_IO_SEQUENTIAL_H8)
  disp=disp+8
#endif
  return
  end subroutine mpi_read_sfc



  subroutine mpi_read_bdy(buf, fh, disp)
  use zocdim
  implicit none
#include "mpif.h"

  real(8) :: buf(nx, ny, nz)
  integer :: fh
  integer (kind = mpi_offset_kind):: disp
  integer :: i, j, k
  integer :: ifile
  integer :: istart(3), igsize(3), isize(3)
  integer :: ifpar, jfpar
  integer(8) :: int1, int2, int3, int4

#ifdef OPT_IO_SEQUENTIAL
  disp=disp+4
#elif defined(OPT_IO_SEQUENTIAL_H8)
  disp=disp+8
#endif
  istart=(/irank*nx, jrank*ny, 0/)
  igsize=(/nxg, nyg, nz/)
  isize =(/nx , ny , nz /)

  call mpi_type_create_subarray(       &
       &    3, igsize, isize, istart,  &
       &    mpi_order_fortran,         &
       &    mpi_real8, ifile, ierr)
  call mpi_type_commit(ifile, ierr)

  call mpi_file_set_view(             &
       &    fh, disp,                 &
       &    mpi_real8,ifile,"native", &
       &    mpi_info_null,ierr)

  call mpi_file_read_all(             &
       &    fh, buf, nx*ny*nz,        &
       &    mpi_real8, mpi_status_ignore, ierr)

  do k=1,nz
     do j=1,ny
        do i=1,nx
           call reverse_real8(buf(i,j,k))
        end do
     end do
  end do

#ifdef OPT_IO_SEQUENTIAL
  disp=disp+ nxg*nyg*nz*8
  disp=disp+4
#elif defined(OPT_IO_SEQUENTIAL_H8)
  int1=nxg
  int2=nyg
  int3=nz
  int4=8
  disp=disp+ int1*int2*int3*int4
  disp=disp+ 8
#else
  int1=nxg
  int2=nyg
  int3=nz
  int4=8
  disp=disp+ int1*int2*int3*int4
#endif

  return
  end subroutine mpi_read_bdy


  subroutine mpi_read_root(buf,nbuf, fh, disp)
  use zocdim
  implicit none
#include "mpif.h"

  real(8) :: buf(nbuf)
  integer :: nbuf, fh, i
  integer (kind = mpi_offset_kind):: disp
#ifdef OPT_IO_SEQUENTIAL
  disp=disp+4
#elif defined(OPT_IO_SEQUENTIAL_H8)
  disp=disp+8
#endif
  call mpi_file_set_view( fh, disp,    &
       &   mpi_real8,mpi_real8,"native", mpi_info_null,ierr)

  if (myrank .eq. iroot) then
     call mpi_file_read(fh, buf, nbuf, &
          &  mpi_real8, mpi_status_ignore, ierr)
     do i=1,nbuf
        call reverse_real8(buf(i))
     end do
  end if
#ifdef OPT_IO_SEQUENTIAL
  disp=disp+ nbuf*8 + 4
#elif defined(OPT_IO_SEQUENTIAL_H8)
  disp=disp+ nbuf*8 + 8
#else
  disp=disp+ nbuf*8
#endif
  return
  end subroutine mpi_read_root


  subroutine mpi_read_root_int(buf,nbuf, fh, disp)
    use zocdim
    implicit none
#include "mpif.h"

    integer :: buf(nbuf)
    integer :: nbuf, fh, i
    integer (kind = mpi_offset_kind):: disp
#ifdef OPT_IO_SEQUENTIAL
    disp=disp+4 
#elif defined(OPT_IO_SEQUENTIAL_H8)
    disp=disp+8
#endif
    call mpi_file_set_view( fh, disp,    &
         &   mpi_integer4,mpi_integer4,"native", mpi_info_null,ierr)

    if (myrank .eq. iroot) then
       call mpi_file_read(fh, buf, nbuf, &
            &  mpi_integer4, mpi_status_ignore, ierr)
       do i=1,nbuf
          call reverse_int4(buf(i))
       end do
    end if
#ifdef OPT_IO_SEQUENTIAL
    disp=disp+ nbuf*4 + 4
#elif defined(OPT_IO_SEQUENTIAL_H8)
    disp=disp+ nbuf*4 + 8
#else
    disp=disp+ nbuf*4
#endif
    return
  end subroutine mpi_read_root_int


  subroutine mpi_read_root_int_sgl(buf, fh, disp)
    use zocdim
    implicit none
#include "mpif.h"

    integer, parameter :: nbuf = 1
    integer :: buf
    integer :: fh, i
    integer (kind = mpi_offset_kind):: disp
    integer :: abuf(nbuf)
#ifdef OPT_IO_SEQUENTIAL
    disp=disp+4 
#elif defined(OPT_IO_SEQUENTIAL_H8)
    disp=disp+8
#endif
    call mpi_file_set_view( fh, disp,    &
         &   mpi_integer4,mpi_integer4,"native", mpi_info_null,ierr)

    if (myrank .eq. iroot) then
       call mpi_file_read(fh, abuf, nbuf, &
            &  mpi_integer4, mpi_status_ignore, ierr)
       call reverse_int4(abuf(1))
       buf=abuf(1)
    end if
#ifdef OPT_IO_SEQUENTIAL
    disp=disp+ nbuf*4 + 4
#elif defined(OPT_IO_SEQUENTIAL_H8)
    disp=disp+ nbuf*4 + 8
#else
    disp=disp+ nbuf*4
#endif
    return
  end subroutine mpi_read_root_int_sgl


  subroutine mpi_read_root_char(buf, nch, nelem, fh, disp)
    use zocdim
    implicit none
#include "mpif.h"

    integer :: nch, nelem, fh, i
    character(len=nch) :: buf(nelem)
    integer(kind=mpi_offset_kind) :: disp

#ifdef OPT_IO_SEQUENTIAL
    disp=disp+4 
#elif defined(OPT_IO_SEQUENTIAL_H8)
    disp=disp+8
#endif
    call mpi_file_set_view(                     &
      &  fh, disp,                              &
      &  mpi_character, mpi_character,"native", &
      &  mpi_info_null,ierr)

    if (myrank .eq. iroot) then
        call mpi_file_read(                          &
          &  fh, buf, nch*nelem,                     &
          &  mpi_character, mpi_status_ignore, ierr)
    end if

    disp=disp+nch*nelem
#ifdef OPT_IO_SEQUENTIAL
    disp=disp + 4
#elif defined(OPT_IO_SEQUENTIAL_H8)
    disp=disp + 8
#endif
    return
    end subroutine mpi_read_root_char


  subroutine mpi_read_2d_intx(buf, fh, disp)
  use zocdim
  implicit none

#include "mpif.h"

  integer :: buf(nxdim, nydim)
  integer :: tmp(nx, ny)
  integer :: fh
  integer (kind = mpi_offset_kind):: disp
  integer   i, j, k
  integer :: ifile
  integer :: istart(2), igsize(2), isize(2)
#ifdef OPT_IO_SEQUENTIAL
  disp=disp+4
#elif defined(OPT_IO_SEQUENTIAL_H8)
  disp=disp+8
#endif
  istart=(/irank*nx +igstr-1, jrank*ny+jgstr-1/)
  igsize=(/nxgdim, nygdim/)
  isize =(/nx    ,    ny /)

  call mpi_type_create_subarray(   &
       & 2, igsize, isize, istart, &
       & mpi_order_fortran,        &
       & mpi_integer4, ifile, ierr)
  call mpi_type_commit(ifile, ierr)

  call mpi_file_set_view(fh, disp,        &
       &     mpi_integer4,ifile,"native", &
       &     mpi_info_null,ierr)

  call mpi_file_read_all(fh, tmp, nx*ny,    &
       &  mpi_integer4, mpi_status_ignore, ierr)

  do j=1,ny
     do i=1,nx
        call reverse_int4(tmp(i,j))
        buf(i+istr-1,j+jstr-1)=tmp(i,j)
     end do
  end do
#ifdef OPT_IO_SEQUENTIAL
  disp=disp+ nxgdim*nygdim*4 + 4
#elif defined(OPT_IO_SEQUENTIAL_H8)
  disp=disp+ nxgdim*nygdim*4 + 8
#else
  disp=disp+ nxgdim*nygdim*4
#endif
  return
  end subroutine mpi_read_2d_intx


  subroutine mpi_read_2d_dimx(buf, fh, disp)
  use zocdim
  implicit none
#include "mpif.h"

  real(8) ::  buf(nxdim, nydim)
  real(8) ::  tmp(nx, ny)
  integer :: fh
  integer (kind = mpi_offset_kind):: disp
  integer :: i, j, k
  integer :: ifile
  integer :: istart(2), igsize(2), isize(2)

#ifdef OPT_IO_SEQUENTIAL
  disp=disp+4
#elif defined(OPT_IO_SEQUENTIAL_H8)
  disp=disp+8
#endif
  istart=(/irank*nx +igstr-1, jrank*ny+jgstr-1/)
  igsize=(/nxgdim, nygdim/)
  isize =(/nx    ,    ny /)

  call mpi_type_create_subarray(      &
       &    2, igsize, isize, istart, &
       &    mpi_order_fortran,        &
       &    mpi_real8, ifile, ierr)
  call mpi_type_commit(ifile, ierr)

  call mpi_file_set_view(fh, disp,    &
       &    mpi_real8,ifile,"native", &
       &    mpi_info_null,ierr)

  call mpi_file_read_all(fh, tmp, nx*ny,    &
       & mpi_real8, mpi_status_ignore, ierr)

  do j=1,ny
     do i=1,nx
        call reverse_real8(tmp(i,j))
        buf(i+istr-1,j+jstr-1)=tmp(i,j)
     end do
  end do

#ifdef OPT_IO_SEQUENTIAL
  disp=disp+ nxgdim*nygdim*8 + 4
#elif defined(OPT_IO_SEQUENTIAL_H8)
  disp=disp+ nxgdim*nygdim*8 + 8
#else
  disp=disp+ nxgdim*nygdim*8
#endif
  return
  end subroutine mpi_read_2d_dimx


  subroutine mpi_read_3d_dimx(buf, fh, disp)
  use zocdim
  implicit none

#include "mpif.h"

  real(8) ::  buf(nxdim, nydim, nzdim)
  real(8) ::  tmp(nx, ny, nzdim)
  integer :: fh
  integer (kind = mpi_offset_kind):: disp
  integer :: i, j, k
  integer :: ifile
  integer :: istart(3), igsize(3), isize(3)
  integer(8) :: int1, int2
#ifdef OPT_IO_SEQUENTIAL
  disp=disp+4
#elif defined(OPT_IO_SEQUENTIAL_H8)
  disp=disp+8
#endif
  istart=(/irank*nx +igstr-1, jrank*ny+jgstr-1, 0/)
  igsize=(/nxgdim, nygdim, nzdim/)
  isize =(/nx    ,    ny , nzdim/)

  call mpi_type_create_subarray(       &
       &    3, igsize, isize, istart,  &
       &    mpi_order_fortran,         &
       &    mpi_real8, ifile, ierr)
  call mpi_type_commit(ifile, ierr)

  call mpi_file_set_view(fh, disp,    &
       &    mpi_real8,ifile,"native", &
       &    mpi_info_null,ierr)

  call mpi_file_read_all(fh, tmp, nx*ny*nzdim, &
       &  mpi_real8, mpi_status_ignore, ierr)

  do k=1,nzdim
     do j=1,ny
        do i=1,nx
           call reverse_real8(tmp(i,j,k))
           buf(i+istr-1,j+jstr-1,k)=tmp(i,j,k)
        end do
     end do
  end do
#ifdef OPT_IO_SEQUENTIAL
  disp=disp+ nxyzgd*8 + 4
#elif defined(OPT_IO_SEQUENTIAL_H8)
  int1=nxyzgd
  int2=8
  disp=disp+ int1*int2
  disp=disp+8
#else
  int1=nxyzgd
  int2=8
  disp=disp+ int1*int2
#endif
  return
  end subroutine mpi_read_3d_dimx



  subroutine mpi_read_2d(buf, fh, disp)
  use zocdim
  implicit none
#include "mpif.h"

  real(8) ::  buf(nxdim, nydim)
  real(8) ::  tmp(nx, ny)
  integer :: fh
  integer (kind = mpi_offset_kind):: disp
  integer ::  i, j, k
  integer :: ifile
  integer :: istart(2), igsize(2), isize(2)
  integer :: ifpar,  jfpar

#ifdef OPT_IO_SEQUENTIAL
  disp=disp+4
#elif defined(OPT_IO_SEQUENTIAL_H8)
  disp=disp+8
#endif
  istart=(/irank*nx, jrank*ny/)
  igsize=(/nxg, nyg/)
  isize =(/nx , ny /)

  call mpi_type_create_subarray(       &
       &    2, igsize, isize, istart,  &
       &    mpi_order_fortran,         &
       &    mpi_real8, ifile, ierr)
  call mpi_type_commit(ifile, ierr)
  call mpi_file_set_view(fh, disp, mpi_real8,ifile,"native", mpi_info_null,ierr)
  call mpi_file_read_all(fh, tmp, nx*ny, mpi_real8, mpi_status_ignore, ierr)

  do j=1,ny
     do i=1,nx
        call reverse_real8(tmp(i,j))
     end do
  end do
  buf(istr:iend,jstr:jend)=tmp(:,:)

  disp=disp+ nxg*nyg*8
#ifdef OPT_IO_SEQUENTIAL
  disp=disp+4
#elif defined(OPT_IO_SEQUENTIAL_H8)
  disp=disp+8
#endif
  return
  end subroutine mpi_read_2d


  subroutine mpi_read_id(buf, fh, disp)
  use zocdim
  implicit none
#include "mpif.h"

  real(8) ::  tmp(nx, ny, nic)
  real(8) ::  buf(nxdim, nydim, 0:nic)
  integer :: fh
  integer (kind = mpi_offset_kind):: disp
  integer :: i, j, k
  integer :: ifile
  integer :: istart(3), igsize(3), isize(3)

#ifdef OPT_IO_SEQUENTIAL
  disp=disp+4
#elif defined(OPT_IO_SEQUENTIAL_H8)
  disp=disp+8
#endif
  istart=(/irank*nx, jrank*ny, 0/)
  igsize=(/nxg, nyg, nic/)
  isize =(/nx , ny , nic/)

  call mpi_type_create_subarray(3, igsize, isize, istart, &
       &   mpi_order_fortran, mpi_real8, ifile, ierr)
  call mpi_type_commit(ifile, ierr)
  call mpi_file_set_view(fh, disp, mpi_real8,ifile,"native",mpi_info_null,ierr)
  call mpi_file_read_all(fh, tmp, nx*ny*nic, mpi_real8, mpi_status_ignore, ierr)

  do k=1,nic
     do j=1,ny
        do i=1,nx
           call reverse_real8(tmp(i,j,k))
           buf(i+istr-1,j+jstr-1,k)=tmp(i,j,k)
        end do
     end do
  end do
  disp=disp+ nxg*nyg*nic*8
#ifdef OPT_IO_SEQUENTIAL
  disp=disp+4
#elif defined(OPT_IO_SEQUENTIAL_H8)
  disp=disp+8
#endif
  return
  end subroutine mpi_read_id


  subroutine mpi_read_3d(buf, fh, disp)
  use zocdim
  implicit none
#include "mpif.h"

  real(8) ::  buf(nxdim, nydim, nzdim)
  real(8) ::  tmp(nx, ny, nz)
  integer :: fh
  integer (kind = mpi_offset_kind):: disp
  integer :: i, j, k
  integer :: ifile
  integer :: istart(3), igsize(3), isize(3)
  integer(8) :: int1, int2, int3, int4

#ifdef OPT_IO_SEQUENTIAL
  disp=disp+4
#elif defined(OPT_IO_SEQUENTIAL_H8)
  disp=disp+8
#endif
  istart=(/irank*nx, jrank*ny, 0/)
  igsize=(/nxg, nyg, nz/)
  isize =(/nx , ny , nz/)

  call mpi_type_create_subarray(3,igsize, isize, istart, &
       &        mpi_order_fortran, mpi_real8, ifile, ierr)
  call mpi_type_commit(ifile, ierr)
  call mpi_file_set_view(fh, disp, mpi_real8,ifile,"native", mpi_info_null,ierr)
  call mpi_file_read_all(fh, tmp, nx*ny*nz, mpi_real8, mpi_status_ignore, ierr)

  do k=1,nz
     do j=1,ny
        do i=1,nx
           call reverse_real8(tmp(i,j,k))
           buf(i+istr-1, j+jstr-1, k+kstr-1)=tmp(i,j,k)
        end do
     end do
  end do

#ifdef OPT_IO_SEQUENTIAL
  disp=disp+ nxg*nyg*nz*8
  disp=disp+4
#elif defined(OPT_IO_SEQUENTIAL_H8)
  int1=nxg
  int2=nyg
  int3=nz
  int4=8
  disp=disp+ int1*int2*int3*int4
  disp=disp+8
#else
  int1=nxg
  int2=nyg
  int3=nz
  int4=8
  disp=disp+ int1*int2*int3*int4
#endif
  return
  end subroutine mpi_read_3d

  subroutine mpi_read_direct(chead, direct, fh, disp, icread)
  use zocdim
  use zocfil, only : nfstdo
  implicit none
#include "mpif.h"
  character, intent(out) :: chead(64)*16
  real(4), intent(out) :: direct(:,:)
  integer :: nx0,ny0
  integer :: fh
  integer (kind = mpi_offset_kind):: disp
  integer :: mpistat(mpi_status_size)
  integer :: i,j
  integer :: icread

  nx0=size(direct,1)
  ny0=size(direct,2)

  !=== chead ====
#ifdef OPT_IO_SEQUENTIAL_H8
  disp=disp+8
#else
  disp=disp+4
#endif
  call mpi_file_set_view(                            &
       &     fh, disp,                               &
       &     mpi_character, mpi_character,"native",  &
       &     mpi_info_null,ierr)

  if (myrank .eq. iroot) then
     call mpi_file_read(                            &
          &    fh, chead, 1024,                     &
          &    mpi_character, mpistat, ierr)

     call mpi_get_count(mpistat,mpi_character, icread,ierr)     
  end if
  call mpi_bcast(chead, 1024, mpi_character,  &
       &                  iroot, mpi_comm_ogcm, ierr)
  call mpi_bcast(icread, 1, mpi_integer4,     &
       &                  iroot, mpi_comm_ogcm, ierr)
  disp=disp+ 1024
#ifdef OPT_IO_SEQUENTIAL_H8
  disp=disp+8
#else
  disp=disp+4
#endif
  if(icread .ne. 1024) then
     write(nfstdo,*)'read error in mpi_read_direct'
     call mpi_abort(mpi_comm_ogcm, 1, ierr)
  end if

  !=== data ====
#ifdef OPT_IO_SEQUENTIAL_H8
  disp=disp+8
#else
  disp=disp+4
#endif
  call mpi_file_set_view( fh, disp,    &
       &   mpi_real4, mpi_real4,"native", mpi_info_null,ierr)

  if (myrank .eq. iroot) then
     call mpi_file_read(fh, direct, nx0*ny0, &
          &  mpi_real4, mpi_status_ignore, ierr)

     do j=1,ny0
        do i=1,nx0
           call reverse_real4(direct(i,j))
        end do
     end do
  end if
#ifdef OPT_IO_SEQUENTIAL_H8
  disp=disp+ nx0*ny0*4 + 8
#else
  disp=disp+ nx0*ny0*4 + 4
#endif

  call mpi_bcast(direct, nx0*ny0, mpi_real4, iroot, mpi_comm_ogcm, ierr)

  return
  end subroutine mpi_read_direct


  subroutine mpi_write_2d(buf, fh, disp)
  use zocdim
  implicit none
#include "mpif.h"

  real(8) :: buf(nxdim, nydim)
  real(8) :: tmp(nx, ny)
  integer :: fh
  integer (kind = mpi_offset_kind):: disp
  integer :: i, j, k, kdim
  integer :: ifile
  integer :: istart(2), igsize(2), isize(2)
  integer :: nsize
  integer(8) :: nsize2

  nsize=8*nxg*nyg
  nsize2=nsize
#ifdef OPT_IO_SEQUENTIAL
  call info_seq(fh, disp, nsize)
#elif defined(OPT_IO_SEQUENTIAL_H8)
  call info_seq8(fh, disp, nsize2)
#endif

  do j=1,ny
     do i=1,nx
        tmp(i,j)=buf(i+istr-1,j+jstr-1)
        call reverse_real8(tmp(i,j))
     end do
  end do

  istart=(/irank*nx, jrank*ny/)
  igsize=(/nxg, nyg/)
  isize =(/nx , ny /)

  call mpi_type_create_subarray(2, igsize, isize, istart, &
       &       mpi_order_fortran, mpi_real8, ifile, ierr)
  call mpi_type_commit(ifile, ierr)
  call mpi_file_set_view(fh, disp, mpi_real8, ifile, "native", mpi_info_null, ierr)
  call mpi_file_write_all(fh, tmp, nx*ny, mpi_real8, mpi_status_ignore, ierr)

  disp=disp+nsize
#ifdef OPT_IO_SEQUENTIAL
  call info_seq(fh, disp, nsize)
#elif defined(OPT_IO_SEQUENTIAL_H8)
  call info_seq8(fh, disp, nsize2)
#endif
  return
  end subroutine mpi_write_2d


  subroutine mpi_write_id(buf, fh, disp)
  use zocdim
  implicit none
#include "mpif.h"

  real(8) ::  buf(nxdim, nydim, 0:nic)
  real(8) ::  tmp(nx, ny, nic)
  integer :: fh
  integer (kind = mpi_offset_kind):: disp
  integer :: i, j, k, kdim
  integer :: ifile
  integer :: istart(3), igsize(3), isize(3)
  integer :: nsize
  integer(8) :: int1, int2, int3, int4,  nsize2

  int1=8
  int2=nxg
  int3=nyg
  int4=nic
  nsize2=int1*int2*int3*int4

  nsize= 8*nxg*nyg*nic
#ifdef OPT_IO_SEQUENTIAL
  call info_seq(fh, disp, nsize)
#elif defined(OPT_IO_SEQUENTIAL_H8)
  call info_seq8(fh, disp, nsize2)
#endif
  do k=1,nic
     do j=1,ny
        do i=1,nx
           tmp(i,j,k)=buf(i+istr-1,j+jstr-1,k)
           call reverse_real8(tmp(i,j,k))
        end do
     end do
  end do

  istart=(/irank*nx, jrank*ny, 0/)
  igsize=(/nxg, nyg, nic/)
  isize =(/nx , ny , nic/)

  call mpi_type_create_subarray(3, igsize, isize, &
       &    istart, mpi_order_fortran, mpi_real8, ifile, ierr)
  call mpi_type_commit(ifile, ierr)
  call mpi_file_set_view(fh, disp,mpi_real8,ifile,"native",mpi_info_null,ierr)
  call mpi_file_write_all(fh, tmp, nx*ny*nic, mpi_real8, mpi_status_ignore, ierr)
  disp=disp+nsize2
#ifdef OPT_IO_SEQUENTIAL
  call info_seq(fh, disp, nsize)
#elif defined(OPT_IO_SEQUENTIAL_H8)
  call info_seq8(fh, disp, nsize2)
#endif
  return
  end subroutine mpi_write_id


  subroutine mpi_write_3d(buf, fh, disp)
  use zocdim
  implicit none
#include "mpif.h"

  real(8) ::  buf(nxdim, nydim, nzdim)
  real(8) ::  tmp(nx, ny, nz)
  integer :: fh
  integer (kind = mpi_offset_kind):: disp
  integer :: i, j, k, kdim
  integer :: ifile
  integer :: istart(3), igsize(3), isize(3)
  integer :: nsize
  integer(8) :: int1, int2, int3, int4,  nsize2

  int1=8
  int2=nxg
  int3=nyg
  int4=nz
  nsize2=int1*int2*int3*int4

  nsize=8*nxg*nyg*nz
#ifdef OPT_IO_SEQUENTIAL
  call info_seq(fh, disp, nsize)
#elif defined(OPT_IO_SEQUENTIAL_H8)
  call info_seq8(fh, disp, nsize2)
#endif

  do k=1,nz
     do j=1,ny
        do i=1,nx
           tmp(i,j,k)=buf(i+istr-1,j+jstr-1,k+kstr-1)
           call reverse_real8(tmp(i,j,k))
        end do
     end do
  end do

  istart=(/irank*nx, jrank*ny, 0/)
  igsize=(/nxg, nyg, nz/)
  isize =(/nx , ny , nz/)

  call mpi_type_create_subarray( 3, igsize, isize, istart, &
       & mpi_order_fortran, mpi_real8, ifile, ierr)
  call mpi_type_commit(ifile, ierr)
  call mpi_file_set_view(fh, disp, mpi_real8,ifile,"native",mpi_info_null,ierr)
  call mpi_file_write_all(fh, tmp, nx*ny*nz, mpi_real8, mpi_status_ignore, ierr)
  disp=disp+nsize2
#ifdef OPT_IO_SEQUENTIAL
  call info_seq(fh, disp, nsize)
#elif defined(OPT_IO_SEQUENTIAL_H8)
  call info_seq8(fh, disp, nsize2)
#endif
  return
  end subroutine mpi_write_3d

  function mpi_iseof(fh, disp)
    use zocdim
    implicit none
#include "mpif.h"

    logical mpi_iseof
    integer fh
    integer (kind = mpi_offset_kind):: disp, siz

    call mpi_file_get_size(fh, siz, ierr)
    if (disp .ge. siz) then
       mpi_iseof = .true.
    else
       mpi_iseof = .false.
    end if
    return
  end function mpi_iseof

end module mpiio

#else
  subroutine mpi_io
  return
  end subroutine mpi_io
#endif
