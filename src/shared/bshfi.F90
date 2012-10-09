module bshfi
! --- information -----------------------------------------------------
!
!  One-on-one communication routines for a 2-dimensional integer
! variable.
!
!  HISTORY
!     '99.10.06  H.Hasumi
!     '07.04.23  H.Hasumi
!     '10.04.14  M.Kurogi: (COCO4.4 tripolar code by Dr. Suzuki)
!     '12.10.09  M.kurogi: for COCO5.0
!
! ---------------------------------------------------------------------
private
public :: shft1i, shftint, shftinv

contains
 subroutine shft1i(                  &
  &                      iq,         &
  &                    idim,   jdim)
 use zocdim, only : nxdim, nydim, istr, iend, jstr, jend, &
  & idown, iup, jdown, jup, icomm, jcomm
 implicit none
#include "mpif.h"

 integer, intent(in   ) ::  idim,   jdim
 integer, intent(inout) ::    iq(idim, jdim)

 integer :: isdbx1(icomm, nydim), isdbx2(icomm, nydim)
 integer :: isdby1(nxdim, jcomm), isdby2(nxdim, jcomm)
 integer :: irvbx1(icomm, nydim), irvbx2(icomm, nydim)
 integer :: irvby1(nxdim, jcomm), irvby2(nxdim, jcomm)
 integer ::      i,      j
 integer :: nbfdim

 if (idown .ne. mpi_proc_null) then
    do j = 1, nydim
       do i = 1, icomm
          isdbx1(i, j) = iq(i+istr-1, j)
       end do
    end do
 end if
 if (iup .ne. mpi_proc_null) then
    do j = 1, nydim
       do i = 1, icomm
          isdbx2(i, j) = iq(i+iend-icomm, j)
       end do
    end do
 end if
 nbfdim = nydim * icomm
 call shftxi(                    &
   &            irvbx1, irvbx2,  &
   &            isdbx1, isdbx2,  &
   &            nbfdim)
 if (idown .ne. mpi_proc_null) then
    do j = 1, nydim
       do i = 1, icomm
          iq(i, j) = irvbx1(i, j)
       end do
    end do
 end if
 if (iup .ne. mpi_proc_null) then
    do j = 1, nydim
       do i = 1, icomm
          iq(iend+i, j) = irvbx2(i, j)
       end do
    end do
 end if

 if (jdown .ne. mpi_proc_null) then
    do j = 1, jcomm
       do i = 1, nxdim
          isdby1(i, j) = iq(i, j+jstr-1)
       end do
    end do
 end if
 if (jup .ne. mpi_proc_null) then
    do j = 1, jcomm
       do i = 1, nxdim
          isdby2(i, j) = iq(i, j+jend-jcomm)
       end do
    end do
 end if
 nbfdim = nxdim * jcomm
 call shftyi(                    &
   &            irvby1, irvby2,  &
   &            isdby1, isdby2,  &
   &            nbfdim)
 if (jdown .ne. mpi_proc_null) then
    do j = 1, jcomm
       do i = 1, nxdim
          iq(i, j) = irvby1(i, j)
       end do
    end do
 end if
 if (jup .ne. mpi_proc_null) then
    do j = 1, jcomm
       do i = 1, nxdim
          iq(i, jend+j) = irvby2(i, j)
       end do
    end do
 end if

 return
 end subroutine shft1i

! *********************************************************************

 subroutine shftxi(             &
  &            irvbx1, irvbx2,  &
  &            isdbx1, isdbx2,  &
  &            nbfdim)
 use zocdim, only : nydim, idown, iup,  icomm, ierr
 implicit none
#include "mpif.h"

 integer, intent(in ) :: isdbx1(icomm, nydim), isdbx2(icomm, nydim)
 integer, intent(out) :: irvbx1(icomm, nydim), irvbx2(icomm, nydim)
 integer :: nbfdim

 integer :: isrqx1, isrqx2
 integer :: irrqx1, irrqx2
 integer :: istmpi(mpi_status_size)

 call mpi_isend(                                   &
  &               isdbx1, nbfdim, mpi_integer,     &
  &                idown,      1, mpi_comm_world,  &
  &               isrqx1,   ierr)
 call mpi_isend(                                   &
  &               isdbx2, nbfdim, mpi_integer,     &
  &                  iup,      2, mpi_comm_world,  &
  &               isrqx2,   ierr)
 call mpi_irecv(                                   &
  &               irvbx2, nbfdim, mpi_integer,     &
  &                  iup,      1, mpi_comm_world,  &
  &               irrqx1,   ierr)
 call mpi_irecv(                                   &
  &               irvbx1, nbfdim, mpi_integer,     &
  &                idown,      2, mpi_comm_world,  &
  &               irrqx2,   ierr)
 call mpi_wait(isrqx1, istmpi,   ierr)
 call mpi_wait(isrqx2, istmpi,   ierr)
 call mpi_wait(irrqx1, istmpi,   ierr)
 call mpi_wait(irrqx2, istmpi,   ierr)

 return
 end subroutine shftxi

! *********************************************************************

 subroutine shftyi(             &
  &            irvby1, irvby2,  &
  &            isdby1, isdby2,  &
  &            nbfdim)
 use zocdim, only : nxdim, jdown, jup,  jcomm, ierr
 implicit none
#include "mpif.h"

 integer, intent(in ) :: isdby1(nxdim, jcomm), isdby2(nxdim, jcomm)
 integer, intent(out) :: irvby1(nxdim, jcomm), irvby2(nxdim, jcomm)
 integer :: nbfdim

 integer :: isrqy1, isrqy2
 integer :: irrqy1, irrqy2
 integer :: istmpi(mpi_status_size)

 call mpi_isend(                                   &
  &               isdby1, nbfdim, mpi_integer,     &
  &                jdown,      3, mpi_comm_world,  &
  &               isrqy1,   ierr)
 call mpi_isend(                                   &
  &               isdby2, nbfdim, mpi_integer,     &
  &                  jup,      4, mpi_comm_world,  &
  &               isrqy2,   ierr)
 call mpi_irecv(                                   &
  &               irvby2, nbfdim, mpi_integer,     &
  &                  jup,      3, mpi_comm_world,  &
  &               irrqy1,   ierr)
 call mpi_irecv(                                   &
  &               irvby1, nbfdim, mpi_integer,     &
  &                jdown,      4, mpi_comm_world,  &
  &               irrqy2,   ierr)
 call mpi_wait(isrqy1, istmpi,   ierr)
 call mpi_wait(isrqy2, istmpi,   ierr)
 call mpi_wait(irrqy1, istmpi,   ierr)
 call mpi_wait(irrqy2, istmpi,   ierr)

 return
 end subroutine shftyi

#ifdef OPT_TRIPOLE
 subroutine shftint(                    &
  &                      iq,            &
  &                    idim,   jdim)
 use zocdim, only : nxdim, ny, jend, jupe, jupw, icomm, jcomm
 implicit none
#include "mpif.h"

 integer, intent(in   ) ::   idim,   jdim
 integer, intent(inout) ::     iq(idim, jdim)

 integer :: isdbx1(icomm, ny), isdbx2(icomm, ny)
 integer :: isdby1(nxdim, jcomm), isdby2(nxdim, jcomm)
 integer :: irvbx1(icomm, ny), irvbx2(icomm, ny)
 integer :: irvby1(nxdim, jcomm), irvby2(nxdim, jcomm)
 integer ::      i,      j
 integer :: nbfdim

 if (jupe .ne. mpi_proc_null) then
    do j = 1, jcomm
       do i = 1, nxdim
          isdby1(i, j) = iq(nxdim-i+1,jend-j+1)
       end do
    end do
 end if
 if (jupw .ne. mpi_proc_null) then
    do j = 1, jcomm
       do i = 1, nxdim
          isdby2(i, j) = iq(nxdim-i+1,jend-j+1)
       end do
    end do
 end if
 nbfdim = nxdim * jcomm
 call shftyin(                  &
  &            irvby1, irvby2,  &
  &            isdby1, isdby2,  &
  &            nbfdim)
 if (jupe .ne. mpi_proc_null) then
    do j = 1, jcomm
       do i = 1, nxdim
          iq(i, jend+j) = irvby1(i, j)
       end do
    end do
 end if
 if (jupw .ne. mpi_proc_null) then
    do j = 1, jcomm
       do i = 1, nxdim
          iq(i, jend+j) = irvby2(i, j)
       end do
    end do
 end if

 return
 end subroutine shftint


 subroutine shftinv(                  &
  &                      iq,          &
  &                    idim,   jdim)
 use zocdim, only : nxdim, ny, jend, jupe, jupw, icomm, jcomm
 implicit none
#include "mpif.h"

 integer, intent(in   ) ::   idim,   jdim
 integer, intent(inout) ::     iq(idim, jdim)

 integer :: isdbx1(icomm, ny), isdbx2(icomm, ny)
 integer :: isdby1(nxdim, jcomm), isdby2(nxdim, jcomm)
 integer :: irvbx1(icomm, ny), irvbx2(icomm, ny)
 integer :: irvby1(nxdim, jcomm), irvby2(nxdim, jcomm)
 integer ::      i,      j
 integer :: nbfdim

 if (jupe .ne. mpi_proc_null) then
    do j = 1, jcomm
       do i = 1, nxdim-1
          isdby1(i, j) = iq(nxdim-i,jend-j)
       end do
          isdby1(nxdim, j) = iq(1,jend-j)
    end do
 end if
 if (jupw .ne. mpi_proc_null) then
    do j = 1, jcomm
       do i = 1, nxdim-1
          isdby2(i, j) = iq(nxdim-i,jend-j)
       end do
          isdby2(nxdim, j) = iq(1,jend-j)
    end do
 end if
 nbfdim = nxdim * jcomm
 call shftyin(                  &
  &            irvby1, irvby2,  &
  &            isdby1, isdby2,  &
  &            nbfdim)
 if (jupe .ne. mpi_proc_null) then
    do j = 1, jcomm
       do i = 1, nxdim
          iq(i, jend+j) = irvby1(i, j)
       end do
    end do
 end if
 if (jupw .ne. mpi_proc_null) then
    do j = 1, jcomm
       do i = 1, nxdim
          iq(i, jend+j) = irvby2(i, j)
       end do
    end do
 end if

 return
 end subroutine shftinv
! *********************************************************************

 subroutine shftyin(            &
  &            irvby1, irvby2,  &
  &            isdby1, isdby2,  &
  &            nbfdim)
 use zocdim, only : nxdim, jcomm, jupe, jupw, ierr
 implicit none
#include "mpif.h"

 integer, intent(in ) :: isdby1(nxdim, jcomm), isdby2(nxdim, jcomm)
 integer, intent(out) :: irvby1(nxdim, jcomm), irvby2(nxdim, jcomm)
 integer :: nbfdim

 integer :: isrqy1, isrqy2
 integer :: irrqy1, irrqy2
 integer :: istmpi(mpi_status_size)

 call mpi_isend(                                    &
  &               isdby1, nbfdim, mpi_integer,      &
  &                 jupe,      3, mpi_comm_world,   &
  &               isrqy1,   ierr)
 call mpi_isend(                                    &
  &               isdby2, nbfdim, mpi_integer,      &
  &                 jupw,      4, mpi_comm_world,   &
  &               isrqy2,   ierr)
 call mpi_irecv(                                    &
  &               irvby2, nbfdim, mpi_integer,      &
  &                 jupw,      3, mpi_comm_world,   &
  &               irrqy1,   ierr)
 call mpi_irecv(                                    &
  &               irvby1, nbfdim, mpi_integer,      &
  &                 jupe,      4, mpi_comm_world,   &
  &               irrqy2,   ierr)
 call mpi_wait(isrqy1, istmpi,   ierr)
 call mpi_wait(isrqy2, istmpi,   ierr)
 call mpi_wait(irrqy1, istmpi,   ierr)
 call mpi_wait(irrqy2, istmpi,   ierr)

 return
 end subroutine shftyin
#endif
end module bshfi
