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
!     '21.05.21  M.Kurogi: packed shift communication from MIROC6 (2013.05.23  Dr. Koji Ogochi)
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

  real(8)        :: qb(nxdim, nydim, nztdim+nzdim)
  
  real(8), allocatable :: sdbffy(:,:,:), rvbffy(:,:,:)
#endif
  
  integer(4)     ::        i,      j,      k,      n
  integer(4)     ::   nbfdim, nbfdm0,   istv
  integer(4)     :: ifpar, jfpar

  logical, save :: shift_gpu=.true.,shift_rdgeo=.false.

  public  ::  shift1,  shift2,  shift3, shift_pack_begin, shift_pack_end, shift_unpack, shift_gpu, shift_rdgeo
#ifdef OPT_TRIPOLE
  public  ::  shiftf1
#endif

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
    call shifts                         &
     &   (    nelems,    nelems,        &
     &     west_recv, east_recv,        &
     &     west_send, east_send,        &
               idown,       iup  )

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
    call shifts                          &
         &   (     nelems,     nelems,   &
         &     south_recv, north_recv,   &
         &     south_send, north_send,   &
         &          jdown,        jup  )

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
    real(8) :: q1(nxdim,nydim,*)
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
          call shifts(        nelems,        nelems, &
               &       tri_west_recv, tri_east_recv, &
               &       tri_west_send, tri_east_send, &
               &               idown,           iup  )

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
    real(8)                                  ::  fact
    integer(4)                               ::  ioff,  joff

    if (oinit) shift_gpu=.false.
    if (shift_rdgeo) then
       shift_gpu=.false.
    end if
    
    if (pack_mode) then
       call shift_pack(q1, kdim, fact, ioff, joff)
    else
       call instant_shift(                                            &
    &                 q1,                                             &
#ifndef OPT_TRIPOLE
    &               idim,   jdim,   kdim )
#else
    &               idim,   jdim,   kdim,                             &
    &               fact,   ioff,   joff )
#endif
    end if
    shift_gpu=.true.
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
    real(8)                                  ::  fact
    integer(4)                               ::  ioff,  joff

    if (oinit) shift_gpu=.false.
    if (shift_rdgeo) then
       shift_gpu=.false.
    end if
    
    if (pack_mode) then
       call shift_pack(q1, kdim, fact, ioff, joff)
       call shift_pack(q2, kdim, fact, ioff, joff)
    else
       !call shift_pack_begin
       !call shift_pack(q1, kdim, fact, ioff, joff)
       !call shift_pack(q2, kdim, fact, ioff, joff)
       !call shift_pack_end
       !call shift_unpack(q1, 1)
       !call shift_unpack(q2, 2)

       !$acc data create(qb) if(shift_gpu)
       !$acc kernels default(present) if(shift_gpu)
       qb(:,:,     1:kdim  ) = q1(:,:,1:kdim)
       qb(:,:,kdim+1:2*kdim) = q2(:,:,1:kdim)
       !$acc end kernels
       call instant_shift(                                            &
    &                 qb,                                             &
#ifndef OPT_TRIPOLE
    &               idim,   jdim, 2*kdim )
#else
    &               idim,   jdim, 2*kdim,                             &
    &               fact,   ioff,   joff )
#endif
       !$acc kernels default(present) if(shift_gpu)
       q1(:,:,1:kdim)=qb(:,:,     1:kdim  )
       q2(:,:,1:kdim)=qb(:,:,kdim+1:2*kdim)
       !$acc end kernels
       !$acc end data
    end if
    shift_gpu=.true.
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
    real(8)                                  ::  fact
    integer(4)                               ::  ioff,  joff

    if (oinit) shift_gpu=.false.
    if (shift_rdgeo) then
       shift_gpu=.false.
    end if
    
    if (pack_mode) then
       call shift_pack(q1, kdim, fact, ioff, joff)
       call shift_pack(q2, kdim, fact, ioff, joff)
       call shift_pack(q3, kdim, fact, ioff, joff)
    else
       !call shift_pack_begin
       !call shift_pack(q1, kdim, fact, ioff, joff)
       !call shift_pack(q2, kdim, fact, ioff, joff)
       !call shift_pack(q3, kdim, fact, ioff, joff)
       !call shift_pack_end
       !call shift_unpack(q1, 1)
       !call shift_unpack(q2, 2)
       !call shift_unpack(q3, 3)

       !$acc data create(qb) if(shift_gpu)
       !$acc kernels default(present) if(shift_gpu)
       qb(:,:,       1:kdim  ) = q1(:,:,1:kdim)
       qb(:,:,  kdim+1:2*kdim) = q2(:,:,1:kdim)
       qb(:,:,2*kdim+1:3*kdim) = q3(:,:,1:kdim)
       !$acc end kernels
       call instant_shift(                                            &
    &                 qb,                                             &
#ifndef OPT_TRIPOLE
    &               idim,   jdim,   3*kdim )
#else
    &               idim,   jdim,   3*kdim,                           &
    &               fact,   ioff,   joff )
#endif
       !$acc kernels default(present) if(shift_gpu)
       q1(:,:,1:kdim)=qb(:,:,       1:kdim  )
       q2(:,:,1:kdim)=qb(:,:,  kdim+1:2*kdim)
       q3(:,:,1:kdim)=qb(:,:,2*kdim+1:3*kdim)
       !$acc end kernels
       !$acc end data
    end if
    shift_gpu=.true.
  end subroutine shift3


!========================================================================================
  subroutine instant_shift(                                           &
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
    logical, save :: ofirst= .true.

    if(ofirst) then
       ofirst = .false.
       !$acc enter data create(sdbfx1,sdbfx2,rvbfx1,rvbfx2, sdbfy1,sdbfy2,rvbfy1,rvbfy2)
#ifdef OPT_TRIPOLE
       !$acc enter data create(sdbfn1,sdbfn2,rvbfn1,rvbfn2, sdbfx1n,sdbfx2n,rvbfx1n,rvbfx2n)
#endif
    end if

    if (idown /= mpi_proc_null) then
       !$acc kernels default(present) if(shift_gpu)
       do k = 1, kdim
          do j = 1, ny
             do i = 1, icomm
                sdbfx1(i, j, k) = q1(i+istr-1, j+jstr-1, k)
             end do
          end do
       end do
       !$acc end kernels
    end if
    if (iup /= mpi_proc_null) then
       !$acc kernels default(present) if(shift_gpu)
       do k = 1, kdim
          do j = 1, ny
             do i = 1, icomm
                sdbfx2(i, j, k) = q1(i+iend-icomm, j+jstr-1, k)
             end do
          end do
       end do
       !$acc end kernels
    end if
    
    nbfdim = kdim * ny * icomm
    call shifts(nbfdim, nbfdim, rvbfx1, rvbfx2, sdbfx1, sdbfx2, idown, iup)
    
    if (idown /= mpi_proc_null) then
       !$acc kernels default(present) if(shift_gpu)
       do k = 1, kdim
          do j = 1, ny
             do i = 1, icomm
                q1(i, j+jstr-1, k) = rvbfx1(i, j, k)
             end do
          end do
       end do
       !$acc end kernels
    end if
    if (iup /= mpi_proc_null) then
       !$acc kernels default(present) if(shift_gpu)
       do k = 1, kdim
          do j = 1, ny
             do i = 1, icomm
                q1(iend+i, j+jstr-1, k) = rvbfx2(i, j, k)
             end do
          end do
       end do
       !$acc end kernels
    end if

    if (jdown /= mpi_proc_null) then
       !$acc kernels default(present) if(shift_gpu)
       do k = 1, kdim
          do j = 1, jcomm
             do i = 1, nxdim
                sdbfy1(i, j, k) = q1(i, j+jstr-1, k)
             end do
          end do
       end do
       !$acc end kernels
    end if
    if (jup /= mpi_proc_null) then
       !$acc kernels default(present) if(shift_gpu)
       do k = 1, kdim
          do j = 1, jcomm
             do i = 1, nxdim
                sdbfy2(i, j, k) = q1(i, j+jend-jcomm, k)
             end do
          end do
       end do
       !$acc end kernels
    end if
    
    nbfdim = kdim * nxdim * jcomm
    call shifts(nbfdim, nbfdim, rvbfy1, rvbfy2, sdbfy1, sdbfy2, jdown, jup)
    
    if (jdown /= mpi_proc_null) then
       !$acc kernels default(present) if(shift_gpu)
       do k = 1, kdim
          do j = 1, jcomm
             do i = 1, nxdim
                q1(i, j, k) = rvbfy1(i, j, k)
             end do
          end do
       end do
       !$acc end kernels
    end if
    if (jup /= mpi_proc_null) then
       !$acc kernels default(present) if(shift_gpu)
       do k = 1, kdim
          do j = 1, jcomm
             do i = 1, nxdim
                q1(i, jend+j, k) = rvbfy2(i, j, k)
             end do
          end do
       end do
       !$acc end kernels
    end if
    
#ifdef OPT_TRIPOLE

    if ( joff == -1 ) then

       if (jupe /= mpi_proc_null) then
          !$acc kernels default(present) if(shift_gpu)
          do k = 1, kdim
             do j = 0, jcomm
                do i = 1+iabs(ioff), nxdim-iabs(ioff)
                   sdbfn1(i, j, k) = fact*q1(nxdim-i+1+ioff, jend-j, k)
                end do
             end do
          end do
          !$acc end kernels
       end if
       if (jupw /= mpi_proc_null) then
          !$acc kernels default(present) if(shift_gpu)
          do k = 1, kdim
             do j = 1, jcomm
                do i = 1+iabs(ioff), nxdim-iabs(ioff)
                   sdbfn2(i, j, k) = fact*q1(nxdim-i+1+ioff, jend-j, k)
                end do
             end do
          end do
          !$acc end kernels
       end if
       
       nbfdim = kdim * nxdim * jcomm
       nbfdm0 = kdim * nxdim * (jcomm+1)
       call shifts(nbfdim, nbfdm0, rvbfn1, rvbfn2, sdbfn1, sdbfn2, jupe, jupw)
       
       if (jupe /= mpi_proc_null) then
          !$acc kernels default(present) if(shift_gpu)
          do k = 1, kdim
             do j = 1, jcomm
                do i = 1, nxdim
                   q1(i, jend+j, k) = rvbfn1(i, j, k)
                end do
             end do
          end do
          !$acc end kernels
       end if

       istv = 1
       if( inodes == 1) istv = nxdim/2 + 1
       if (jupw /= mpi_proc_null) then
          !$acc kernels default(present) if(shift_gpu)
          do k = 1, kdim
             do j = 0, jcomm
                do i = istv, nxdim
                   q1(i, jend+j, k) = rvbfn2(i, j, k)
                end do
             end do
          end do
          !$acc end kernels
       end if
!-----------
    else
       if (jupe /= mpi_proc_null) then
          !$acc kernels default(present) if(shift_gpu)
          do k = 1, kdim
             do j = 1, jcomm
                do i = 1+iabs(ioff), nxdim-iabs(ioff)
                   sdbfy1(i, j, k) = fact*q1(nxdim-i+1+ioff,jend-j+1+joff,k)
                end do
             end do
          end do
          !$acc end kernels
       end if
       if (jupw /= mpi_proc_null) then
          !$acc kernels default(present) if(shift_gpu)
          do k = 1, kdim
             do j = 1, jcomm
                do i = 1+iabs(ioff), nxdim-iabs(ioff)
                   sdbfy2(i, j, k) = fact*q1(nxdim-i+1+ioff,jend-j+1+joff,k)
                end do
             end do
          end do
          !$acc end kernels
       end if
       
       nbfdim = kdim * nxdim * jcomm
       call shifts(nbfdim, nbfdim, rvbfy1, rvbfy2, sdbfy1, sdbfy2, jupe, jupw)
       
       if (jupe /= mpi_proc_null) then
          !$acc kernels default(present) if(shift_gpu)
          do k = 1, kdim
             do j = 1, jcomm
                do i = 1, nxdim
                   q1(i, jend+j, k) = rvbfy1(i, j, k)
                end do
             end do
          end do
          !$acc end kernels
       end if
       if (jupw /= mpi_proc_null) then
          !$acc kernels default(present) if(shift_gpu)
          do k = 1, kdim
             do j = 1, jcomm
                do i = 1, nxdim
                   q1(i, jend+j, k) = rvbfy2(i, j, k)
                end do
             end do
          end do
          !$acc end kernels
       end if
    end if

    if ( ioff /= 0 .and. jrank == jnodes-1 ) then

       if (idown /= mpi_proc_null) then
          !$acc kernels default(present) if(shift_gpu)
          do k = 1, kdim
             do j = 1, jcomm+1
                do i = 1, icomm
                   sdbfx1n(i, j, k) = q1(i+istr-1, j+jend-1, k)
                end do
             end do
          end do
          !$acc end kernels
       end if
       if (iup /= mpi_proc_null) then
          !$acc kernels default(present) if(shift_gpu)
          do k = 1, kdim
             do j = 1, jcomm+1
                do i = 1, icomm
                   sdbfx2n(i, j, k) = q1(i+iend-icomm, j+jend-1, k)
                end do
             end do
          end do
          !$acc end kernels
       end if
       
       nbfdim = kdim * (jcomm+1) * icomm
       call shifts(nbfdim, nbfdim, rvbfx1n, rvbfx2n, sdbfx1n, sdbfx2n, idown, iup)
       
       if (idown /= mpi_proc_null) then
          !$acc kernels default(present) if(shift_gpu)
          do k = 1, kdim
             do j = 1, jcomm+1
                do i = 1, icomm
                   q1(i, j+jend-1, k) = rvbfx1n(i, j, k)
                end do
             end do
          end do
          !$acc end kernels
       end if
       if (iup /= mpi_proc_null) then
          !$acc kernels default(present) if(shift_gpu)
          do k = 1, kdim
             do j = 1, jcomm+1
                do i = 1, icomm
                   q1(iend+i, j+jend-1, k) = rvbfx2n(i, j, k)
                end do
             end do
          end do
          !$acc end kernels
       end if
    end if
#endif
  end subroutine instant_shift
  
  subroutine shifts(nbfdim, nbfdim0, rbf1, rbf2, sbf1, sbf2, n_down, n_up)
    implicit none
#include "mpif.h"
    real(8),      intent(inout) ::  rbf1(:,:,:)
    real(8),      intent(inout) ::  rbf2(:,:,:)
    real(8),      intent(in)    ::  sbf1(:,:,:)
    real(8),      intent(in)    ::  sbf2(:,:,:)
    integer(4),   intent(in)    ::  nbfdim, nbfdim0, n_down, n_up
!---- internal work
    integer(4)  ::  is1, is2, ir1, ir2
    integer(4)  ::  istmpi(mpi_status_size)

    !$acc host_data use_device(sbf1) if(shift_gpu)
    call mpi_isend(sbf1, nbfdim0, mpi_real8, n_down, 1, mpi_comm_world, is1, ierr)
    !$acc end host_data
    
    !$acc host_data use_device(sbf2) if(shift_gpu)
    call mpi_isend(sbf2, nbfdim,  mpi_real8,   n_up, 2, mpi_comm_world, is2, ierr) 
    !$acc end host_data
    
    !$acc host_data use_device(rbf2) if(shift_gpu)
    call mpi_irecv(rbf2, nbfdim0, mpi_real8,   n_up, 1, mpi_comm_world, ir1, ierr)
    !$acc end host_data
    
    !$acc host_data use_device(rbf1) if(shift_gpu)
    call mpi_irecv(rbf1, nbfdim,  mpi_real8, n_down, 2, mpi_comm_world, ir2, ierr)
    !$acc end host_data
    
    call mpi_wait( is1, istmpi, ierr )
    call mpi_wait( is2, istmpi, ierr )
    call mpi_wait( ir1, istmpi, ierr )
    call mpi_wait( ir2, istmpi, ierr )
  end subroutine shifts


#ifdef OPT_TRIPOLE
!========================================================================================
    
  subroutine shiftf1(                                                &
   &                 q1,                                             &
   &               idim,   jdim,   kdim )

  implicit none
    
#include "mpif.h"

  real(8),                  intent(inout)  ::    q1(1:idim,1:jdim,1:kdim)
  integer(4),               intent(in)     ::  idim,  jdim,  kdim

  allocate(sdbffy(nxdim, 1, kdim), rvbffy(nxdim, 1, kdim))
  if (jupfy .ne. mpi_proc_null) then
     do k = 1, kdim
        do i = 1, nxdim
           sdbffy(i, 1, k) = q1(i, jend+1, k)
        end do
     end do
  end if
  nbfdim = kdim * nxdim * 1
  call shiftfy( &
    &         rvbffy, &
    &         sdbffy, &
    &         nbfdim, kdim)
  if (jdownfy .ne. mpi_proc_null) then
     do k = 1, kdim
        do i = 1, nxdim
           q1(i, jstr, k) = rvbffy(i, 1, k)
        end do
     end do
  end if
  deallocate(sdbffy, rvbffy)

  return
  end subroutine shiftf1

! *********************************************************************
  subroutine shiftfy(                                                 &
    &            rvbffy,                                              &
    &            sdbffy,                                              &
    &            nbfdim, nzsdim )

  implicit none

#include "mpif.h"

  integer(4),   intent(in)    ::  nbfdim, nzsdim
  real(8),      intent(inout) ::  rvbffy(nxdim, 1, nzsdim)
  real(8),      intent(in)    ::  sdbffy(nxdim, 1, nzsdim)

!---- internal work
  integer(4)  ::  isrqfy1
  integer(4)  ::  irrqfy1
  integer(4)  ::  istmpi(mpi_status_size)

  call mpi_isend( &
    &               sdbffy, nbfdim, mpi_real8, &
    &                jupfy,      5, mpi_comm_world, &
    &              isrqfy1,   ierr)
  call mpi_irecv( &
    &               rvbffy, nbfdim, mpi_real8, &
    &              jdownfy,      5, mpi_comm_world, &
    &              irrqfy1,   ierr)

  call mpi_wait(isrqfy1, istmpi,   ierr)
  call mpi_wait(irrqfy1, istmpi,   ierr)

  return
  end subroutine shiftfy
#endif

end module bshft
