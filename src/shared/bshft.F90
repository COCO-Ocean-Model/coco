module bshft

! --- information -----------------------------------------------------
!
!  One-on-one communication routines for real*8 variables.
!
!  HISTORY
!     '99.10.06  H.Hasumi
!     '02.06.02  H.Hasumi
!     '07.04.23  H.Hasumi
!     '10.04.14  M.Kurogi: (COCO4.4 tripolar code by Dr. Suzuki)
!     '12.11.28  H.Tatebe: for COCO5.0 in F90
!     '13.02.12  Y.Komuro: bug fix (for non-tripole code)
!     '21.05.21  M.Kurogi: packed shift communication from MIROC5 (2013.05.23  Dr. Koji Ogochi)
! ---------------------------------------------------------------------

!  use zocdim,  only  :                                                &
!        nxdim,   nzdim,  nztdim,                                      &
!           nx,      ny,      nz,                                      &
!        icomm,   jcomm
  use zocdim
  use ufile, only : rewnml 
  implicit none
  
  private

!  [internal save]
  integer, parameter :: max_ksize0 = 10 * nzdim
  integer, parameter :: max_ksize  = 20 * nzdim * ntdim
  integer, parameter :: max_num_packed = 32

  real(8), save :: north_send(nxdim, jcomm, max_ksize)
  real(8), save :: south_send(nxdim, jcomm, max_ksize)
  real(8), save :: east_send (icomm, ny,    max_ksize)
  real(8), save :: west_send (icomm, ny,    max_ksize)

  real(8), save :: north_recv(nxdim, jcomm, max_ksize)
  real(8), save :: south_recv(nxdim, jcomm, max_ksize)
  real(8), save :: east_recv (icomm, ny,    max_ksize)
  real(8), save :: west_recv (icomm, ny,    max_ksize)

  logical, save :: pack_mode = .false.
  integer, save :: koffset(max_num_packed + 1)
  integer, save :: num_packed = 0
  logical, save :: is_tri_edge

#ifdef OPT_TRIPOLE
  real(8), save :: tri_send(nxdim, 0:jcomm, max_ksize)
  real(8), save :: tri_recv(nxdim, 0:jcomm, max_ksize)

  real(8), save :: tri_east_send(icomm, jcomm+1, max_ksize0)
  real(8), save :: tri_west_send(icomm, jcomm+1, max_ksize0)
  real(8), save :: tri_east_recv(icomm, jcomm+1, max_ksize0)
  real(8), save :: tri_west_recv(icomm, jcomm+1, max_ksize0)

  real(8), save :: facts(max_num_packed)
  integer, save :: ioffs(max_num_packed)
  integer, save :: joffs(max_num_packed)
#endif

  
  real(8)        ::   sdbfx1(1:icomm, 1:ny,    1:nztdim+nzdim)
  real(8)        ::   sdbfx2(1:icomm, 1:ny,    1:nztdim+nzdim)
  real(8)        ::   sdbfy1(1:nxdim, 1:jcomm, 1:nztdim+nzdim)
  real(8)        ::   sdbfy2(1:nxdim, 1:jcomm, 1:nztdim+nzdim)
  real(8)        ::   rvbfx1(1:icomm, 1:ny,    1:nztdim+nzdim)
  real(8)        ::   rvbfx2(1:icomm, 1:ny,    1:nztdim+nzdim)
  real(8)        ::   rvbfy1(1:nxdim, 1:jcomm, 1:nztdim+nzdim)
  real(8)        ::   rvbfy2(1:nxdim, 1:jcomm, 1:nztdim+nzdim)

#ifdef OPT_TRIPOLE
  real(8)        ::  sdbfx1n(1:icomm, 1:jcomm+1, 1:nztdim+nzdim)
  real(8)        ::  sdbfx2n(1:icomm, 1:jcomm+1, 1:nztdim+nzdim)
  real(8)        ::  rvbfx1n(1:icomm, 1:jcomm+1, 1:nztdim+nzdim)
  real(8)        ::  rvbfx2n(1:icomm, 1:jcomm+1, 1:nztdim+nzdim)
  
  real(8)        ::   sdbfn1(1:nxdim, 0:jcomm, 1:nztdim+nzdim)
  real(8)        ::   sdbfn2(1:nxdim, 1:jcomm, 1:nztdim+nzdim)
  real(8)        ::   rvbfn1(1:nxdim, 1:jcomm, 1:nztdim+nzdim)
  real(8)        ::   rvbfn2(1:nxdim, 0:jcomm, 1:nztdim+nzdim)
#endif

  integer(4)     ::        i,      j,      k,      n
  integer(4)     ::   nbfdim, nbfdm0,   istv
  integer(4)     :: ifpar, jfpar

  public  ::  shift1,  shift2,  shift3, shift_pack_begin, shift_pack_end

contains

!=======================================================================
  subroutine shift_pack_begin
    implicit none
#include "mpif.h"
    if (pack_mode) then
       call rewnml(ifpar, jfpar)
       write(jfpar,*)' ### shift_pack_begin: Illegal call'
       call mpi_abort(mpi_comm_ogcm, 1, ierr)
    end if
    pack_mode = .true.
    num_packed = 0
    koffset(1) = 0
    is_tri_edge = jrank .eq. jnodes -1

  end subroutine shift_pack_begin

!=======================================================================
  subroutine shift_pack_end
    integer :: nelems
    
    if (.not. pack_mode) return
    nelems = icomm * ny * koffset(num_packed + 1)
    call shiftx                         &
     &   ( west_recv, east_recv,        &
     &     west_send, east_send, nelems )

    do k = 1, koffset(num_packed + 1)
       do j = 1, jcomm
          do i = 1, icomm
             south_send(i,j,k) = west_recv(i,j,k)
             north_send(i,j,k) = west_recv(i,j-jcomm+ny,k)

             south_send(iend+i,j,k) = east_recv(i,j,k)
             north_send(iend+i,j,k) = east_recv(i,j-jcomm+ny,k)
          end do
       end do
    end do

    nelems = nxdim * jcomm * koffset(num_packed + 1)
    call shifty                          &
         &   ( south_recv, north_recv,   &
         &     south_send, north_send, nelems)

#ifdef OPT_TRIPOLE
    if (is_tri_edge) then
       do n = 1, num_packed
          if (ioffs(n) .eq. -1) then
             do k = koffset(n) + 1, koffset(n + 1)
                do j = 1, jcomm - joffs(n)
                   do i = 1, icomm - 1
                      tri_send(i,j+joffs(n),k) = facts(n) * east_recv(icomm - i, ny + 1 - j, k)
                   end do
                   do i = iend, nxdim - 1
                      tri_send(i,j+joffs(n),k) = facts(n) * west_recv(nxdim - i, ny + 1 - j, k)
                   end do
                end do
             end do
          else
             do k = koffset(n) + 1, koffset(n + 1)
                do j = 1, jcomm - joffs(n)
                   do i = 1, icomm
                      tri_send(i     ,j+joffs(n),k) = facts(n) * east_recv(icomm + 1 - i, ny + 1 - j, k)
                      tri_send(iend+i,j+joffs(n),k) = facts(n) * west_recv(icomm + 1 - i, ny + 1 - j, k)
                   end do
                end do
             end do
          endif
       end do

       nelems = nxdim * (jcomm + 1) * koffset(num_packed + 1)
       call shift_tri_edge      &
            &      ( tri_recv,  &
            &        tri_send, nelems )
    endif
#endif
    pack_mode = .false.
  end subroutine shift_pack_end

!=======================================================================
  subroutine shift_unpack(q1, id)
    implicit none
#include "mpif.h"
    real(8) :: q1(:,:,:)
    integer :: id, nelems, k0, kpacked
    
    if (id .lt. 1 .or. id .gt. num_packed) return

    k0 = koffset(id)
    kpacked = koffset(id+1) - k0

    do k = 1, kpacked
       do j = 1, ny
          do i = 1, icomm
             q1(i,      j+jstr-1, k) = west_recv(i, j, k0+k)
             q1(i+iend, j+jstr-1, k) = east_recv(i, j, k0+k)
          end do
       end do
    end do

    if (jdown .ne. mpi_proc_null) then
       do k = 1, kpacked
          do j = 1, jcomm
             do i = 1, nxdim
                q1(i, j, k) = south_recv(i, j, k0+k)
             end do
          end do
       end do
    end if

    if (jup .ne. mpi_proc_null) then
       do k = 1, kpacked
          do j = 1, jcomm
             do i = 1, nxdim
                q1(i, j+jend, k) = north_recv(i, j, k0+k)
             end do
          end do
       end do
    end if

#ifdef OPT_TRIPOLE
    if (is_tri_edge) then
       if (joffs(id) .eq. -1 .and. jupw .ne. mpi_proc_null) then
          if (inodes .eq. 1) then
             do k = 1, kpacked
                do i = nxdim / 2 + 1, nxdim
                   q1(i,jend,k) = tri_recv(i,0,k0+k)
                end do
             end do
          else
             do k = 1, kpacked
                do i = 1, nxdim
                   q1(i,jend,k) = tri_recv(i,0,k0+k)
                end do
             end do
          endif
       endif

       do k = 1, kpacked
          do j = 1, jcomm
             do i = 1, nxdim
                q1(i, jend+j, k) = tri_recv(i, j, k0+k)
             end do
          end do
       end do

       if (ioffs(id) .eq. -1) then
          if (kpacked .gt. max_ksize0) then
             call rewnml(ifpar, jfpar)
             write(jfpar,*)' ### packed_shift: exceed the limit of max_ksize0.'
             call mpi_abort(mpi_comm_ogcm, 1, ierr)
          endif

          do k = 1, kpacked
             do j = 1, jcomm + 1
                do i = 1, icomm
                   tri_west_send(i,j,k) = q1(i+istr-1    , j+jend-1, k)
                   tri_east_send(i,j,k) = q1(i+iend-icomm, j+jend-1, k)
                end do
             end do
          end do

          nelems = icomm * (jcomm + 1) * kpacked
          call shiftx( tri_west_recv, tri_east_recv, &
               &       tri_west_send, tri_east_send, nelems )

          do k = 1, kpacked
             do j = 1, jcomm + 1
                do i = 1, icomm
                   q1(i,     j+jend-1,k) = tri_west_recv(i,j,k)
                   q1(iend+i,j+jend-1,k) = tri_east_recv(i,j,k)
                end do
             end do
          end do
       endif
    endif
#endif
  end subroutine shift_unpack

!=======================================================================
  subroutine shift_pack(qq, ksize, fact2, ioff2, joff2)
    real(8) :: qq(:,:,:)
    integer :: ksize
    real(8) :: fact2
    integer :: ioff2, joff2
    integer i, j, k, k0

    if (num_packed .ge. max_num_packed) then
       call rewnml(ifpar, jfpar)
       write(jfpar,*)' ### packed_shift: exceed the limit of num_packed.'
       call mpi_abort(mpi_comm_ogcm, 1, ierr)
    endif

    num_packed = num_packed + 1
    k0 = koffset(num_packed)

    if (k0 + ksize .gt. max_ksize) then
       call rewnml(ifpar, jfpar)
       write(jfpar,*)' ### packed_shift: exceed the limit of max_ksize.'
       call mpi_abort(mpi_comm_ogcm, 1, ierr)
    endif

    do k = 1, ksize
       do j = 1, ny
          do i = 1, icomm
             west_send(i,j,k0+k) = qq(i+istr-1,    j+jstr-1,k)
             east_send(i,j,k0+k) = qq(i+iend-icomm,j+jstr-1,k)
          end do
       end do
    end do

    do k = 1, ksize
       do j = 1, jcomm
          do i = istr, iend
             south_send(i,j,k0+k) = qq(i,j+jstr-1,    k)
             north_send(i,j,k0+k) = qq(i,j+jend-jcomm,k)
          end do
       end do
    end do
    koffset(num_packed + 1) = k0 + ksize

#ifdef OPT_TRIPOLE
    if (is_tri_edge) then
       do k = 1, ksize
          do j = 1, jcomm - joff2
             do i = istr, iend
                tri_send(i+ioff2,j+joff2,k0+k) = fact2 * qq(nxdim + 1 - i, jend + 1 - j, k)
             end do
          end do
       end do
    endif
    facts(num_packed) = fact2
    ioffs(num_packed) = ioff2
    joffs(num_packed) = joff2
#endif
  end subroutine shift_pack

!=======================================================================

#ifdef OPT_TRIPOLE
  subroutine shift_tri_edge(recvbuf, sendbuf, nelems)
    implicit none
#include "mpif.h"

    real(8) :: recvbuf(*)
    real(8) :: sendbuf(*)
    integer :: nelems

    integer :: isrc, idest
    integer :: status(mpi_status_size)
    integer, parameter :: itag = 1234

    if (jupe .ne. mpi_proc_null) then
       isrc  = jupe
       idest = jupe
    else
       isrc  = jupw
       idest = jupw
    endif

    if (isrc .eq. myrank) then
       recvbuf(1:nelems) = sendbuf(1:nelems)
       return
    endif

    call mpi_sendrecv(sendbuf, nelems, mpi_double_precision,       &
         &                  idest, itag,                           &
         &                  recvbuf, nelems, mpi_double_precision, &
         &                  isrc, itag,                            &
         &                  mpi_comm_ogcm, status, ierr)

  end subroutine shift_tri_edge
#endif
!=======================================================================
  
  subroutine shift1(                                                  &
    &                 q1,                                             &
#ifndef OPT_TRIPOLE    
    &               idim,   jdim,   kdim )
#else
    &               idim,   jdim,   kdim,                             & 
    &               fact,   ioff,   joff )
#endif
    implicit none
    real(8),                  intent(inout)  ::    q1(1:idim,1:jdim,1:kdim)
    integer(4),               intent(in)     ::  idim,  jdim,  kdim
#ifdef OPT_TRIPOLE    
    real(8),                  intent(in)     ::  fact
    integer(4),               intent(in)     ::  ioff,  joff
#endif
    
    if (.not. pack_mode) then
       call instant_shift1(                                           &
    &                 q1,                                             &
#ifndef OPT_TRIPOLE    
    &               idim,   jdim,   kdim )
#else
    &               idim,   jdim,   kdim,                             & 
    &               fact,   ioff,   joff )
#endif

    end if
  end subroutine shift1

  subroutine shift2(                                                  &
    &                 q1,     q2,                                     &
#ifndef OPT_TRIPOLE
    &               idim,   jdim,   kdim )
#else
    &               idim,   jdim,   kdim,                             &
    &               fact,   ioff,   joff )
#endif
    implicit none
    real(8),                  intent(inout)  ::    q1(1:idim,1:jdim,1:kdim)
    real(8),                  intent(inout)  ::    q2(1:idim,1:jdim,1:kdim)
    integer(4),               intent(in)     ::  idim,  jdim,  kdim
#ifdef OPT_TRIPOLE
    real(8),                  intent(in)     ::  fact
    integer(4),               intent(in)     ::  ioff,  joff
#endif

    if (.not. pack_mode) then
       call instant_shift2(                                           &
    &                 q1,     q2,                                     &
#ifndef OPT_TRIPOLE
    &               idim,   jdim,   kdim )
#else
    &               idim,   jdim,   kdim,                             &
    &               fact,   ioff,   joff )
#endif
    end if
  end subroutine shift2

  subroutine shift3(                      &
    &                 q1,     q2,     q3, &
#ifndef OPT_TRIPOLE
    &               idim,   jdim,   kdim )
#else
    &               idim,   jdim,   kdim, &
    &               fact,   ioff,   joff )
#endif
    implicit none     
    real(8),                  intent(inout)  ::    q1(1:idim,1:jdim,1:kdim)
    real(8),                  intent(inout)  ::    q2(1:idim,1:jdim,1:kdim)
    real(8),                  intent(inout)  ::    q3(1:idim,1:jdim,1:kdim)
    integer(4),               intent(in)     ::  idim,  jdim,  kdim
#ifdef OPT_TRIPOLE
    real(8),                  intent(in)     ::  fact
    integer(4),               intent(in)     ::  ioff,  joff
#endif
    if (.not. pack_mode) then
       call instant_shift3(               &
    &                 q1,     q2,     q3, &
#ifndef OPT_TRIPOLE
    &               idim,   jdim,   kdim )
#else
    &               idim,   jdim,   kdim, &
    &               fact,   ioff,   joff )
#endif
    end if
  end subroutine shift3


!========================================================================================
    
  subroutine instant_shift1(                                          &
    &                 q1,                                             &
#ifndef OPT_TRIPOLE    
    &               idim,   jdim,   kdim )
#else
    &               idim,   jdim,   kdim,                             & 
    &               fact,   ioff,   joff )
#endif

    implicit none
     
#include "mpif.h"

    real(8),                  intent(inout)  ::    q1(1:idim,1:jdim,1:kdim)
    integer(4),               intent(in)     ::  idim,  jdim,  kdim
#ifdef OPT_TRIPOLE    
    real(8),                  intent(in)     ::  fact
    integer(4),               intent(in)     ::  ioff,  joff
#endif

    if (idown /= mpi_proc_null) then
       do k = 1, kdim
          do j = 1, ny
             do i = 1, icomm
                sdbfx1(i, j, k) = q1(i+istr-1, j+jstr-1, k)
             end do
          end do
       end do
    end if
    if (iup /= mpi_proc_null) then
       do k = 1, kdim
          do j = 1, ny
             do i = 1, icomm
                sdbfx2(i, j, k) = q1(i+iend-icomm, j+jstr-1, k)
             end do
          end do
       end do
    end if
    nbfdim = kdim * ny * icomm
    call shiftx(                                                      &
    &            rvbfx1, rvbfx2,                                      &
    &            sdbfx1, sdbfx2,                                      &
    &            nbfdim )
    if (idown /= mpi_proc_null) then
       do k = 1, kdim
          do j = 1, ny
             do i = 1, icomm
                q1(i, j+jstr-1, k) = rvbfx1(i, j, k)
             end do
          end do
       end do
    end if
    if (iup /= mpi_proc_null) then
       do k = 1, kdim
          do j = 1, ny
             do i = 1, icomm
                q1(iend+i, j+jstr-1, k) = rvbfx2(i, j, k)
             end do
          end do
       end do
    end if

    if (jdown /= mpi_proc_null) then
       do k = 1, kdim
          do j = 1, jcomm
             do i = 1, nxdim
                sdbfy1(i, j, k) = q1(i, j+jstr-1, k)
             end do
          end do
       end do
    end if
    if (jup /= mpi_proc_null) then
       do k = 1, kdim
          do j = 1, jcomm
             do i = 1, nxdim
                sdbfy2(i, j, k) = q1(i, j+jend-jcomm, k)
             end do
          end do
       end do
    end if
    nbfdim = kdim * nxdim * jcomm
    call shifty(                                                      &
    &            rvbfy1, rvbfy2,                                      &
    &            sdbfy1, sdbfy2,                                      &
    &            nbfdim  )
    if (jdown /= mpi_proc_null) then
       do k = 1, kdim
          do j = 1, jcomm
             do i = 1, nxdim
                q1(i, j, k) = rvbfy1(i, j, k)
             end do
          end do
       end do
    end if
    if (jup /= mpi_proc_null) then
       do k = 1, kdim
          do j = 1, jcomm
             do i = 1, nxdim
                q1(i, jend+j, k) = rvbfy2(i, j, k)
             end do
          end do
       end do
    end if
    
#ifdef OPT_TRIPOLE

    if ( joff == -1 ) then

       if (jupe /= mpi_proc_null) then
          do k = 1, kdim
             do j = 0, jcomm
                do i = 1+iabs(ioff), nxdim-iabs(ioff)
                   sdbfn1(i, j, k) = fact*q1(nxdim-i+1+ioff, jend-j, k)
                end do
             end do
          end do
       end if
       if (jupw /= mpi_proc_null) then
          do k = 1, kdim
             do j = 1, jcomm
                do i = 1+iabs(ioff), nxdim-iabs(ioff)
                   sdbfn2(i, j, k) = fact*q1(nxdim-i+1+ioff, jend-j, k)
                end do
             end do
          end do
       end if
       nbfdim = kdim * nxdim * jcomm
       nbfdm0 = kdim * nxdim * (jcomm+1)
       call shiftnv(                                                  &
    &            rvbfn1, rvbfn2,                                      &
    &            sdbfn1, sdbfn2,                                      &
    &            nbfdim, nbfdm0  )
       if (jupe /= mpi_proc_null) then
          do k = 1, kdim
             do j = 1, jcomm
                do i = 1, nxdim
                   q1(i, jend+j, k) = rvbfn1(i, j, k)
                end do
             end do
          end do
       end if

       istv = 1
       if( inodes == 1) istv = nxdim/2 + 1
       if (jupw /= mpi_proc_null) then
          do k = 1, kdim
             do j = 0, jcomm
                do i = istv, nxdim
                   q1(i, jend+j, k) = rvbfn2(i, j, k)
                end do
             end do
          end do
       end if
!-----------
    else
       if (jupe /= mpi_proc_null) then
          do k = 1, kdim
             do j = 1, jcomm
                do i = 1+iabs(ioff), nxdim-iabs(ioff)
                   sdbfy1(i, j, k) = fact*q1(nxdim-i+1+ioff,jend-j+1+joff,k)
                end do
             end do
          end do
       end if
       if (jupw /= mpi_proc_null) then
          do k = 1, kdim
             do j = 1, jcomm
                do i = 1+iabs(ioff), nxdim-iabs(ioff)
                   sdbfy2(i, j, k) = fact*q1(nxdim-i+1+ioff,jend-j+1+joff,k)
                end do
             end do
          end do
       end if
       nbfdim = kdim * nxdim * jcomm
       call shiftyn(                                                  &
    &            rvbfy1, rvbfy2,                                      &
    &            sdbfy1, sdbfy2,                                      &
    &            nbfdim  )
       if (jupe /= mpi_proc_null) then
          do k = 1, kdim
             do j = 1, jcomm
                do i = 1, nxdim
                   q1(i, jend+j, k) = rvbfy1(i, j, k)
                end do
             end do
          end do
       end if
       if (jupw /= mpi_proc_null) then
          do k = 1, kdim
             do j = 1, jcomm
                do i = 1, nxdim
                   q1(i, jend+j, k) = rvbfy2(i, j, k)
                end do
             end do
          end do
       end if
    end if

    if ( ioff /= 0 .and. jrank == jnodes-1 ) then

       if (idown /= mpi_proc_null) then
          do k = 1, kdim
             do j = 1, jcomm+1
                do i = 1, icomm
                   sdbfx1n(i, j, k) = q1(i+istr-1, j+jend-1, k)
                end do
             end do
          end do
       end if
       if (iup /= mpi_proc_null) then
          do k = 1, kdim
             do j = 1, jcomm+1
                do i = 1, icomm
                   sdbfx2n(i, j, k) = q1(i+iend-icomm, j+jend-1, k)
                end do
             end do
          end do
       end if
       nbfdim = kdim * (jcomm+1) * icomm
       call shiftx(                                                   &
    &            rvbfx1n, rvbfx2n,                                    &
    &            sdbfx1n, sdbfx2n,                                    &
    &            nbfdim  )
       if (idown /= mpi_proc_null) then
          do k = 1, kdim
             do j = 1, jcomm+1
                do i = 1, icomm
                   q1(i, j+jend-1, k) = rvbfx1n(i, j, k)
                end do
             end do
          end do
       end if
       if (iup /= mpi_proc_null) then
          do k = 1, kdim
             do j = 1, jcomm+1
                do i = 1, icomm
                   q1(iend+i, j+jend-1, k) = rvbfx2n(i, j, k)
                end do
             end do
          end do
       end if
    end if

#endif

  end subroutine instant_shift1

! =====================================================================

  subroutine instant_shift2(                                          &
    &                 q1,     q2,                                     &
#ifndef OPT_TRIPOLE
    &               idim,   jdim,   kdim )
#else
    &               idim,   jdim,   kdim,                             &
    &               fact,   ioff,   joff )
#endif

    implicit none

#include "mpif.h"

    real(8),                  intent(inout)  ::    q1(1:idim,1:jdim,1:kdim)
    real(8),                  intent(inout)  ::    q2(1:idim,1:jdim,1:kdim)
    integer(4),               intent(in)     ::  idim,  jdim,  kdim
#ifdef OPT_TRIPOLE
    real(8),                  intent(in)     ::  fact
    integer(4),               intent(in)     ::  ioff,  joff
#endif

    if (idown /= mpi_proc_null) then
       do k = 1, kdim
          do j = 1, ny
             do i = 1, icomm
                sdbfx1(i, j, k) = q1(i+istr-1, j+jstr-1, k)
                sdbfx1(i, j, k+kdim) = q2(i+istr-1, j+jstr-1, k)
             end do
          end do
       end do
    end if
    if (iup /= mpi_proc_null) then
       do k = 1, kdim
          do j = 1, ny
             do i = 1, icomm
                sdbfx2(i, j, k) = q1(i+iend-icomm, j+jstr-1, k)
                sdbfx2(i, j, k+kdim) = q2(i+iend-icomm, j+jstr-1, k)
             end do
          end do
       end do
    end if
    nbfdim = 2 * kdim * ny * icomm
    call shiftx(                                                      &
    &            rvbfx1, rvbfx2,                                      &
    &            sdbfx1, sdbfx2,                                      &
    &            nbfdim  )
    if (idown /= mpi_proc_null) then
       do k = 1, kdim
          do j = 1, ny
             do i = 1, icomm
                q1(i, j+jstr-1, k) = rvbfx1(i, j, k)
                q2(i, j+jstr-1, k) = rvbfx1(i, j, k+kdim)
             end do
          end do
       end do
    end if
    if (iup /= mpi_proc_null) then
       do k = 1, kdim
          do j = 1, ny
             do i = 1, icomm
                q1(iend+i, j+jstr-1, k) = rvbfx2(i, j, k)
                q2(iend+i, j+jstr-1, k) = rvbfx2(i, j, k+kdim)
             end do
          end do
       end do
    end if

    if (jdown /= mpi_proc_null) then
       do k = 1, kdim
          do j = 1, jcomm
             do i = 1, nxdim
                sdbfy1(i, j, k) = q1(i, j+jstr-1, k)
                sdbfy1(i, j, k+kdim) = q2(i, j+jstr-1, k)
             end do
          end do
       end do
    end if
    if (jup /= mpi_proc_null) then
       do k = 1, kdim
          do j = 1, jcomm
             do i = 1, nxdim
                sdbfy2(i, j, k) = q1(i, j+jend-jcomm, k)
                sdbfy2(i, j, k+kdim) = q2(i, j+jend-jcomm, k)
             end do
          end do
       end do
    end if
    nbfdim = 2 * kdim * nxdim * jcomm
    call shifty(                                                      &
    &            rvbfy1, rvbfy2,                                      &
    &            sdbfy1, sdbfy2,                                      &
    &            nbfdim  )
    if (jdown /= mpi_proc_null) then
       do k = 1, kdim
          do j = 1, jcomm
             do i = 1, nxdim
                q1(i, j, k) = rvbfy1(i, j, k)
                q2(i, j, k) = rvbfy1(i, j, k+kdim)
             end do
          end do
       end do
    end if
    if (jup /= mpi_proc_null) then
       do k = 1, kdim
          do j = 1, jcomm
             do i = 1, nxdim
                q1(i, jend+j, k) = rvbfy2(i, j, k)
                q2(i, jend+j, k) = rvbfy2(i, j, k+kdim)
             end do
          end do
       end do
    end if
#ifdef OPT_TRIPOLE
    if ( joff == -1 ) then

       if (jupe /= mpi_proc_null) then
          do k = 1, kdim
             do j = 0, jcomm
                do i = 1+iabs(ioff), nxdim-iabs(ioff)
                   sdbfn1(i, j, k) = fact*q1(nxdim-i+1+ioff, jend-j, k)
                   sdbfn1(i, j, k+kdim) = fact*q2(nxdim-i+1+ioff, jend-j, k)
                end do
             end do
          end do
       end if
       if (jupw /= mpi_proc_null) then
          do k = 1, kdim
             do j = 1, jcomm
                do i = 1+iabs(ioff), nxdim-iabs(ioff)
                   sdbfn2(i, j, k) = fact*q1(nxdim-i+1+ioff, jend-j, k)
                   sdbfn2(i, j, k+kdim) = fact*q2(nxdim-i+1+ioff, jend-j, k)
                end do
             end do
          end do
       end if
       nbfdim = 2 * kdim * nxdim * jcomm
       nbfdm0 = 2 * kdim * nxdim * (jcomm+1)
       call shiftnv(                                                  &
    &            rvbfn1, rvbfn2,                                      &
    &            sdbfn1, sdbfn2,                                      &
    &            nbfdim, nbfdm0  )
       if (jupe /= mpi_proc_null) then
          do k = 1, kdim
             do j = 1, jcomm
                do i = 1, nxdim
                   q1(i, jend+j, k) = rvbfn1(i, j, k)
                   q2(i, jend+j, k) = rvbfn1(i, j, k+kdim)
                end do
             end do
          end do
       end if

       istv=1
       if(inodes == 1) istv=nxdim/2+1
       if (jupw /= mpi_proc_null) then
          do k = 1, kdim
             do j = 0, jcomm
                do i = istv, nxdim
                   q1(i, jend+j, k) = rvbfn2(i, j, k)
                   q2(i, jend+j, k) = rvbfn2(i, j, k+kdim)
                end do
             end do
          end do
       end if

    else
       if (jupe /= mpi_proc_null) then
          do k = 1, kdim
             do j = 1, jcomm
                do i = 1+iabs(ioff), nxdim-iabs(ioff)
                   sdbfy1(i, j, k) = fact*q1(nxdim-i+1+ioff,jend-j+1+joff,k)
                   sdbfy1(i, j, k+kdim) = fact*q2(nxdim-i+1+ioff,jend-j+1+joff,k)
                end do
             end do
          end do
       end if
       if (jupw /= mpi_proc_null) then
          do k = 1, kdim
             do j = 1, jcomm
                do i = 1+iabs(ioff), nxdim-iabs(ioff)
                   sdbfy2(i, j, k) = fact*q1(nxdim-i+1+ioff,jend-j+1+joff,k)
                   sdbfy2(i, j, k+kdim) = fact*q2(nxdim-i+1+ioff,jend-j+1+joff,k)
                end do
             end do
          end do
       end if
       nbfdim = 2 * kdim * nxdim * jcomm
       call shiftyn(                                                  &
    &            rvbfy1, rvbfy2,                                      &
    &            sdbfy1, sdbfy2,                                      &
    &            nbfdim  )
       if (jupe /= mpi_proc_null) then
          do k = 1, kdim
             do j = 1, jcomm
                do i = 1, nxdim
                   q1(i, jend+j, k) = rvbfy1(i, j, k)
                   q2(i, jend+j, k) = rvbfy1(i, j, k+kdim)
                end do
             end do
          end do
       end if
       if (jupw /= mpi_proc_null) then
          do k = 1, kdim
             do j = 1, jcomm
                do i = 1, nxdim
                   q1(i, jend+j, k) = rvbfy2(i, j, k)
                   q2(i, jend+j, k) = rvbfy2(i, j, k+kdim)
                end do
             end do
          end do
       end if

    end if


    if ( ioff /= 0 .and. jrank == jnodes-1 ) then
       if (idown /= mpi_proc_null) then
          do k = 1, kdim
             do j = 1, jcomm+1
                do i = 1, icomm
                   sdbfx1n(i, j, k) = q1(i+istr-1, j+jend-1, k)
                   sdbfx1n(i, j, k+kdim) = q2(i+istr-1, j+jend-1, k)
                end do
             end do
          end do
       end if
       if (iup /= mpi_proc_null) then
          do k = 1, kdim
             do j = 1, jcomm+1
                do i = 1, icomm
                   sdbfx2n(i, j, k) = q1(i+iend-icomm, j+jend-1, k)
                   sdbfx2n(i, j, k+kdim) = q2(i+iend-icomm, j+jend-1, k)
                end do
             end do
          end do
       end if
       nbfdim = 2 * kdim * (jcomm+1) * icomm
       call shiftx(                                                   &
    &            rvbfx1n, rvbfx2n,                                    &
    &            sdbfx1n, sdbfx2n,                                    &
    &            nbfdim  )
       if (idown /= mpi_proc_null) then
          do k = 1, kdim
             do j = 1, jcomm+1
                do i = 1, icomm
                   q1(i, j+jend-1, k) = rvbfx1n(i, j, k)
                   q2(i, j+jend-1, k) = rvbfx1n(i, j, k+kdim)
                end do
             end do
          end do
       end if
       if (iup /= mpi_proc_null) then
          do k = 1, kdim
             do j = 1, jcomm+1
                do i = 1, icomm
                   q1(iend+i, j+jend-1, k) = rvbfx2n(i, j, k)
                   q2(iend+i, j+jend-1, k) = rvbfx2n(i, j, k+kdim)
                end do
             end do
          end do
       end if
    end if

#endif
 
  end subroutine instant_shift2

! =====================================================================

  subroutine instant_shift3(              &
    &                 q1,     q2,     q3, &
#ifndef OPT_TRIPOLE
    &               idim,   jdim,   kdim )
#else
    &               idim,   jdim,   kdim, &
    &               fact,   ioff,   joff )
#endif

    implicit none
     
#include "mpif.h"

    real(8),                  intent(inout)  ::    q1(1:idim,1:jdim,1:kdim)
    real(8),                  intent(inout)  ::    q2(1:idim,1:jdim,1:kdim)
    real(8),                  intent(inout)  ::    q3(1:idim,1:jdim,1:kdim)
    integer(4),               intent(in)     ::  idim,  jdim,  kdim
#ifdef OPT_TRIPOLE
    real(8),                  intent(in)     ::  fact
    integer(4),               intent(in)     ::  ioff,  joff
#endif

    if (idown /= mpi_proc_null) then
       do k = 1, kdim
          do j = 1, ny
             do i = 1, icomm
                sdbfx1(i, j, k) = q1(i+istr-1, j+jstr-1, k)
                sdbfx1(i, j, k+kdim) = q2(i+istr-1, j+jstr-1, k)
                sdbfx1(i, j, k+kdim*2) = q3(i+istr-1, j+jstr-1, k)
             end do
          end do
       end do
    end if
    if (iup /= mpi_proc_null) then
       do k = 1, kdim
          do j = 1, ny
             do i = 1, icomm
                sdbfx2(i, j, k) = q1(i+iend-icomm, j+jstr-1, k)
                sdbfx2(i, j, k+kdim) = q2(i+iend-icomm, j+jstr-1, k)
                sdbfx2(i, j, k+kdim*2) = q3(i+iend-icomm, j+jstr-1, k)
             end do
          end do
       end do
    end if
    nbfdim = 3 * kdim * ny * icomm
    call shiftx(                                                      &
    &            rvbfx1, rvbfx2,                                      &
    &            sdbfx1, sdbfx2,                                      &
    &            nbfdim  ) 
    if (idown /= mpi_proc_null) then
       do k = 1, kdim
          do j = 1, ny
             do i = 1, icomm
                q1(i, j+jstr-1, k) = rvbfx1(i, j, k)
                q2(i, j+jstr-1, k) = rvbfx1(i, j, k+kdim)
                q3(i, j+jstr-1, k) = rvbfx1(i, j, k+kdim*2)
             end do
          end do
       end do
    end if
    if (iup /= mpi_proc_null) then
       do k = 1, kdim
          do j = 1, ny
             do i = 1, icomm
                q1(iend+i, j+jstr-1, k) = rvbfx2(i, j, k)
                q2(iend+i, j+jstr-1, k) = rvbfx2(i, j, k+kdim)
                q3(iend+i, j+jstr-1, k) = rvbfx2(i, j, k+kdim*2)
             end do
          end do
       end do
    end if

    if (jdown /= mpi_proc_null) then
       do k = 1, kdim
          do j = 1, jcomm
             do i = 1, nxdim
                sdbfy1(i, j, k) = q1(i, j+jstr-1, k)
                sdbfy1(i, j, k+kdim) = q2(i, j+jstr-1, k)
                sdbfy1(i, j, k+kdim*2) = q3(i, j+jstr-1, k)
             end do
          end do
       end do
    end if
    if (jup /= mpi_proc_null) then
       do k = 1, kdim
          do j = 1, jcomm
             do i = 1, nxdim
                sdbfy2(i, j, k) = q1(i, j+jend-jcomm, k)
                sdbfy2(i, j, k+kdim) = q2(i, j+jend-jcomm, k)
                sdbfy2(i, j, k+kdim*2) = q3(i, j+jend-jcomm, k)
             end do
          end do
       end do
    end if
    nbfdim = 3 * kdim * nxdim * jcomm
    call shifty(                                                      &
    &            rvbfy1, rvbfy2,                                      &
    &            sdbfy1, sdbfy2,                                      &
    &            nbfdim  )
    if (jdown /= mpi_proc_null) then
       do k = 1, kdim
          do j = 1, jcomm
             do i = 1, nxdim
                q1(i, j, k) = rvbfy1(i, j, k)
                q2(i, j, k) = rvbfy1(i, j, k+kdim)
                q3(i, j, k) = rvbfy1(i, j, k+kdim*2)
             end do
          end do
       end do
    end if
    if (jup /= mpi_proc_null) then
       do k = 1, kdim
          do j = 1, jcomm
             do i = 1, nxdim
                q1(i, jend+j, k) = rvbfy2(i, j, k)
                q2(i, jend+j, k) = rvbfy2(i, j, k+kdim)
                q3(i, jend+j, k) = rvbfy2(i, j, k+kdim*2)
             end do
          end do
       end do
    end if
#ifdef OPT_TRIPOLE
    if ( joff == -1 ) then

       if (jupe /= mpi_proc_null) then
          do k = 1, kdim
             do j = 0, jcomm
                do i = 1+iabs(ioff), nxdim-iabs(ioff)
                   sdbfn1(i, j, k) = fact*q1(nxdim-i+1+ioff, jend-j, k)
                   sdbfn1(i, j, k+kdim) = fact*q2(nxdim-i+1+ioff, jend-j, k)
                   sdbfn1(i, j, k+2*kdim) = fact*q3(nxdim-i+1+ioff, jend-j, k)
                end do
             end do
          end do
       end if
       if (jupw /= mpi_proc_null) then
          do k = 1, kdim
             do j = 1, jcomm
                do i = 1+iabs(ioff), nxdim-iabs(ioff)
                   sdbfn2(i, j, k) = fact*q1(nxdim-i+1+ioff, jend-j, k)
                   sdbfn2(i, j, k+kdim) = fact*q2(nxdim-i+1+ioff, jend-j, k)
                   sdbfn2(i, j, k+2*kdim) = fact*q3(nxdim-i+1+ioff, jend-j, k)
                end do
             end do
          end do
       end if
       nbfdim = 3 * kdim * nxdim * jcomm
       nbfdm0 = 3 * kdim * nxdim * (jcomm+1)
       call shiftnv(                                                  &
    &            rvbfn1, rvbfn2,                                      &
    &            sdbfn1, sdbfn2,                                      &
    &            nbfdim, nbfdm0  )
       if (jupe /= mpi_proc_null) then
          do k = 1, kdim
             do j = 1, jcomm
                do i = 1, nxdim
                   q1(i, jend+j, k) = rvbfn1(i, j, k)
                   q2(i, jend+j, k) = rvbfn1(i, j, k+kdim)
                   q3(i, jend+j, k) = rvbfn1(i, j, k+2*kdim)
                end do
             end do
          end do
       end if

       istv=1
       if(inodes == 1) istv=nxdim/2+1
       if (jupw /= mpi_proc_null) then
          do k = 1, kdim
             do j = 0, jcomm
                do i = istv, nxdim
                   q1(i, jend+j, k) = rvbfn2(i, j, k)
                   q2(i, jend+j, k) = rvbfn2(i, j, k+kdim)
                   q3(i, jend+j, k) = rvbfn2(i, j, k+2*kdim)
                end do
             end do
          end do
       end if

    else
       if (jupe /= mpi_proc_null) then
          do k = 1, kdim
             do j = 1, jcomm
                do i = 1+iabs(ioff), nxdim-iabs(ioff)
                   sdbfy1(i,j,k) = fact*q1(nxdim-i+1+ioff,jend-j+1+joff,k)
                   sdbfy1(i,j,k+kdim) = fact*q2(nxdim-i+1+ioff,jend-j+1+joff,k)
                   sdbfy1(i,j,k+2*kdim) = fact*q3(nxdim-i+1+ioff,jend-j+1+joff,k)
                end do
             end do
          end do
       end if
       if (jupw /= mpi_proc_null) then
          do k = 1, kdim
             do j = 1, jcomm
                do i = 1+iabs(ioff), nxdim-iabs(ioff)
                   sdbfy2(i,j,k) = fact*q1(nxdim-i+1+ioff,jend-j+1+joff,k)
                   sdbfy2(i,j,k+kdim) = fact*q2(nxdim-i+1+ioff,jend-j+1+joff,k)
                   sdbfy2(i,j,k+2*kdim) = fact*q3(nxdim-i+1+ioff,jend-j+1+joff,k)
                end do
             end do
          end do
       end if
       nbfdim = 3 * kdim * nxdim * jcomm
       call shiftyn(                                                  &
    &            rvbfy1, rvbfy2,                                      &
    &            sdbfy1, sdbfy2,                                      &
    &            nbfdim  )
       if (jupe /= mpi_proc_null) then
          do k = 1, kdim
             do j = 1, jcomm
                do i = 1, nxdim
                   q1(i, jend+j, k) = rvbfy1(i, j, k)
                   q2(i, jend+j, k) = rvbfy1(i, j, k+kdim)
                   q3(i, jend+j, k) = rvbfy1(i, j, k+2*kdim)
                end do
             end do
          end do
       end if
       if (jupw /= mpi_proc_null) then
          do k = 1, kdim
             do j = 1, jcomm
                do i = 1, nxdim
                   q1(i, jend+j, k) = rvbfy2(i, j, k)
                   q2(i, jend+j, k) = rvbfy2(i, j, k+kdim)
                   q3(i, jend+j, k) = rvbfy2(i, j, k+2*kdim)
                end do
             end do
          end do
       end if
    end if


    if ( ioff /=0 .and. jrank == jnodes-1 ) then
       if (idown /= mpi_proc_null) then
          do k = 1, kdim
             do j = 1, jcomm+1
                do i = 1, icomm
                   sdbfx1n(i, j, k) = q1(i+istr-1, j+jend-1, k)
                   sdbfx1n(i, j, k+kdim) = q2(i+istr-1, j+jend-1, k)
                   sdbfx1n(i, j, k+kdim*2) = q3(i+istr-1, j+jend-1, k)
                end do
             end do
          end do
       end if
       if (iup /= mpi_proc_null) then
          do k = 1, kdim
             do j = 1, jcomm+1
                do i = 1, icomm
                   sdbfx2n(i, j, k) = q1(i+iend-icomm, j+jend-1, k)
                   sdbfx2n(i, j, k+kdim) = q2(i+iend-icomm, j+jend-1, k)
                   sdbfx2n(i, j, k+kdim*2) = q3(i+iend-icomm, j+jend-1, k)
                end do
             end do
          end do
       end if
       nbfdim = 3 * kdim * (jcomm+1) * icomm
       call shiftx(                                                   &
    &            rvbfx1n, rvbfx2n,                                    &
    &            sdbfx1n, sdbfx2n,                                    &
    &            nbfdim  )
       if (idown /= mpi_proc_null) then
          do k = 1, kdim
             do j = 1, jcomm+1
                do i = 1, icomm
                   q1(i, j+jend-1, k) = rvbfx1n(i, j, k)
                   q2(i, j+jend-1, k) = rvbfx1n(i, j, k+kdim)
                   q3(i, j+jend-1, k) = rvbfx1n(i, j, k+kdim*2)
                end do
             end do
          end do
       end if
       if (iup /= mpi_proc_null) then
          do k = 1, kdim
             do j = 1, jcomm+1
                do i = 1, icomm
                   q1(iend+i, j+jend-1, k) = rvbfx2n(i, j, k)
                   q2(iend+i, j+jend-1, k) = rvbfx2n(i, j, k+kdim)
                   q3(iend+i, j+jend-1, k) = rvbfx2n(i, j, k+kdim*2)
                end do
             end do
          end do
       end if
    end if

#endif

  end subroutine instant_shift3

! *********************************************************************

  subroutine shiftx(                                                  &
    &            rvbfx1l, rvbfx2l,                                    &
    &            sdbfx1l, sdbfx2l,                                    &
    &            nbfdim  )

    implicit none

#include "mpif.h"

    real(8),      intent(inout) ::  rvbfx1l(:,:,:)
    real(8),      intent(inout) ::  rvbfx2l(:,:,:)
    real(8),      intent(in)    ::  sdbfx1l(:,:,:)
    real(8),      intent(in)    ::  sdbfx2l(:,:,:)
    integer(4),   intent(inout) ::  nbfdim

!---- internal work
    integer(4)  ::  isrqx1, isrqx2
    integer(4)  ::  irrqx1, irrqx2
    integer(4)  ::  istmpi(mpi_status_size)

    call mpi_isend(                                                  &
    &              sdbfx1l, nbfdim, mpi_real8,                       &
    &                idown,      1, mpi_comm_ogcm,                  &
    &               isrqx1,   ierr)
    call mpi_isend(                                                  &
    &              sdbfx2l, nbfdim, mpi_real8,                       &
    &                  iup,      2, mpi_comm_ogcm,                  &
    &               isrqx2,   ierr) 
    call mpi_irecv(                                                  &
    &              rvbfx2l, nbfdim, mpi_real8,                       &
    &                  iup,      1, mpi_comm_ogcm,                  &
    &               irrqx1,   ierr)
    call mpi_irecv(                                                  &
    &              rvbfx1l, nbfdim, mpi_real8,                       &
    &                idown,      2, mpi_comm_ogcm,                  &
    &               irrqx2,   ierr)

    call mpi_wait( isrqx1, istmpi, ierr )
    call mpi_wait( isrqx2, istmpi, ierr )
    call mpi_wait( irrqx1, istmpi, ierr )
    call mpi_wait( irrqx2, istmpi, ierr )

  end subroutine shiftx

! *********************************************************************

  subroutine shifty(                                                  &
    &            rvbfy1l, rvbfy2l,                                    &
    &            sdbfy1l, sdbfy2l,                                    &
    &            nbfdim )

    implicit none

#include "mpif.h"

    real(8),      intent(inout) ::  rvbfy1l(:,:,:)
    real(8),      intent(inout) ::  rvbfy2l(:,:,:)
    real(8),      intent(in)    ::  sdbfy1l(:,:,:)
    real(8),      intent(in)    ::  sdbfy2l(:,:,:)
    integer(4),   intent(inout) ::  nbfdim

!---- internal work
    integer(4)  ::  isrqy1, isrqy2
    integer(4)  ::  irrqy1, irrqy2
    integer(4)  ::  istmpi(mpi_status_size)

    call mpi_isend(                                                   &
    &              sdbfy1l, nbfdim, mpi_real8,                        &
    &                jdown,      3, mpi_comm_ogcm,                   &
    &               isrqy1,   ierr)
    call mpi_isend(                                                   &
    &              sdbfy2l, nbfdim, mpi_real8,                        &
    &                  jup,      4, mpi_comm_ogcm,                   &
    &               isrqy2,   ierr)
    call mpi_irecv(                                                   &
    &              rvbfy2l, nbfdim, mpi_real8,                        &
    &                  jup,      3, mpi_comm_ogcm,                   &
    &               irrqy1,   ierr)
    call mpi_irecv(                                                   &
    &              rvbfy1l, nbfdim, mpi_real8,                        &
    &                jdown,      4, mpi_comm_ogcm,                   &
    &               irrqy2,   ierr)

    call mpi_wait( isrqy1, istmpi, ierr )
    call mpi_wait( isrqy2, istmpi, ierr )
    call mpi_wait( irrqy1, istmpi, ierr )
    call mpi_wait( irrqy2, istmpi, ierr )

  end subroutine shifty

#ifdef OPT_TRIPOLE

  subroutine shiftyn(                                                 &
    &            rvbfy1l, rvbfy2l,                                    &
    &            sdbfy1l, sdbfy2l,                                    &
    &            nbfdim )

    implicit none

#include "mpif.h"

    real(8),      intent(inout) ::  rvbfy1l(:,:,:)
    real(8),      intent(inout) ::  rvbfy2l(:,:,:)
    real(8),      intent(in)    ::  sdbfy1l(:,:,:)
    real(8),      intent(in)    ::  sdbfy2l(:,:,:)
    integer(4),   intent(inout) ::  nbfdim

!---- internal work
    integer(4)  ::  isrqy1, isrqy2
    integer(4)  ::  irrqy1, irrqy2
    integer(4)  ::  istmpi(mpi_status_size)

    call mpi_isend(                                                   &
    &             sdbfy1l, nbfdim, mpi_real8,                         &
    &                jupe,       3, mpi_comm_ogcm,                   &
    &               isrqy1,   ierr)
    call mpi_isend(                                                   &
    &              sdbfy2l, nbfdim, mpi_real8,                        &
    &                 jupw,      4, mpi_comm_ogcm,                   &
    &               isrqy2,   ierr) 
    call mpi_irecv(                                                   &
    &              rvbfy2l, nbfdim, mpi_real8,                        &
    &                 jupw,      3, mpi_comm_ogcm,                   &
    &               irrqy1,   ierr)
    call mpi_irecv(                                                   &
    &              rvbfy1l, nbfdim, mpi_real8,                        &
    &                 jupe,      4, mpi_comm_ogcm,                   &
    &               irrqy2,   ierr)

    call mpi_wait( isrqy1, istmpi, ierr )
    call mpi_wait( isrqy2, istmpi, ierr )
    call mpi_wait( irrqy1, istmpi, ierr )
    call mpi_wait( irrqy2, istmpi, ierr )

  end subroutine shiftyn

  subroutine shiftnv(                                                 &
    &            rvbfn1l, rvbfn2l,                                    &
    &            sdbfn1l, sdbfn2l,                                    &
    &             nbfdim, nbfdm0   )

    implicit none

#include "mpif.h"

    real(8),      intent(inout) ::  rvbfn1l(:,:,:)
    real(8),      intent(inout) ::  rvbfn2l(:,:,:)
    real(8),      intent(in)    ::  sdbfn1l(:,:,:)
    real(8),      intent(in)    ::  sdbfn2l(:,:,:)
    integer(4),   intent(inout) ::  nbfdim,   nbfdm0

!---- internal work
    integer(4)  ::  isrqy1, isrqy2
    integer(4)  ::  irrqy1, irrqy2
    integer(4)  ::  istmpi(mpi_status_size)

    call mpi_isend(                                                   &
    &              sdbfn1l, nbfdm0, mpi_real8,                        &
    &                jupe,       3, mpi_comm_ogcm,                   &
    &               isrqy1,   ierr)
    call mpi_isend(                                                   &
    &              sdbfn2l, nbfdim, mpi_real8,                        &
    &                 jupw,      4, mpi_comm_ogcm,                   &
    &               isrqy2,   ierr) 
    call mpi_irecv(                                                   &
    &              rvbfn2l, nbfdm0, mpi_real8,                        &
    &                 jupw,      3, mpi_comm_ogcm,                   &
    &               irrqy1,   ierr)
    call mpi_irecv(                                                   &
    &              rvbfn1l, nbfdim, mpi_real8,                        &
    &                 jupe,      4, mpi_comm_ogcm,                   &
    &               irrqy2,   ierr)

    call mpi_wait( isrqy1, istmpi, ierr )
    call mpi_wait( isrqy2, istmpi, ierr )
    call mpi_wait( irrqy1, istmpi, ierr )
    call mpi_wait( irrqy2, istmpi, ierr )

  end subroutine shiftnv

#endif

end module bshft
