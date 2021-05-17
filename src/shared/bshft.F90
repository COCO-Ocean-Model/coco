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
! ---------------------------------------------------------------------

!  use zocdim,  only  :                                                &
!        nxdim,   nzdim,  nztdim,                                      &
!           nx,      ny,      nz,                                      &
!        icomm,   jcomm
  use zocdim
  use ufile, only : rewnml 
  implicit none
  integer, parameter :: max_num_packed = 32
  
  logical, save :: pack_mode = .false.
  integer, save :: koffset(max_num_packed +1)
  integer, save :: num_packed = 0
  logical, save :: is_tri_edge
  
  private

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
  integer(4)     :: ifpar, jfpar, ierr
  
  public  ::  shift1,  shift2,  shift3, shift_pack_begin

contains

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

