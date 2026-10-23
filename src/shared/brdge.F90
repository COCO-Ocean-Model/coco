module brdge
! --- information -----------------------------------------------------
!
!  Reading grid information and masking array
!
!  HISTORY
!     '03.04.21  H.Hasumi: from COCO3.4
!     '07.09.25  H.Hasumi: multi-category sea ice
!     '07.11.06  T.Suzuki: for ES
!     '07.12.20  H.Hasumi: MASK file format is modified
!     '10.04.14  M.Kurogi: (COCO4.4 tripolar code by Dr. Suzuki)
!     '12.07.20  M.kurogi: for COCO5.0
!     '15.04.07  M.Kurogi: MPI-IO
! ---------------------------------------------------------------------
implicit none
private 
public rdgeo

contains
 subroutine rdgeo
 use zocdim, only :                                                            &
  & nxdim, nydim, nzdim,nic, nxydim, nxygdm, nxyzgd, nxygdm,                   &
  & jstr, jend, ijstr, ijend,  kstr, kend,                                     & 
  & ijtstr, ijtend, ijvstr, ijvend,                                            &
#ifdef OPT_TRIPOLE
  & jupe, jupw,                                                                &
#endif
  & iroot, ierr,  jup, jdown, myrank, mpi_comm_ogcm


#ifdef OPT_BBL
 use zocmsk, only :                                                            &
  & amskt, amftx, amfty, amftz,  amskv,  amfvx,  amfvy,  amfvz,  amskb,  nbot, &
  & amsktb, amskt0, amskt1, amskvb, amskv0, amskv1, nbotv
#else
 use zocmsk, only :                                                            &
  &  amskt, amftx, amfty, amftz,  amskv,  amfvx,  amfvy,  amfvz,  amskb,  nbot 
#endif


 use zocgrd, only :                                                            &
  & dx, dy, dym, dz, dzm, dzv, dz0, ds, dsm, hic, rea, zbot, cor, dept, rdepv, &
  & hxt, hxu, hyt, hyu, hxyt, hxyu, hyxt, hyxu,                                &
#ifdef OPT_EXMASK
  & glont, glatt, rangt,                                                       &
#endif
  & rx, ry, rym, rxt, rxu, ryt, ryu, rs, rsm

 use zocfil, only : ncf
 use ufile
 use bshfi
 use bshft
#ifdef OPT_IO_COCOMPI
 use mpiio
#else
 use bgs2d
 use bgs3d
#endif

#ifdef OPT_OFFLINE 
 use bchmk
#endif

 implicit none
#include "mpif.h"
#include "coco.h"
#ifndef OPT_IO_COCOMPI
 integer,allocatable :: g2di(:)
 real(8),allocatable :: g2d(:), g3d(:)
#endif
  integer ::  n
  integer :: ij,      k
  integer ::  i,      j,    ijs,   ije,   ijn,  ijw
  integer ::  ifpar,  jfpar
  integer :: nfmask
  integer :: istat

#ifdef OPT_IO_COCOMPI
  integer :: mpi_fh
  integer (kind = mpi_offset_kind) :: disp
#endif
  character(ncf) :: cfmask
  data cfmask / 'not-specified' /
  namelist /nmmask/ cfmask

  call rewnml(ifpar, jfpar)
  write(jfpar, *) '*** rdgeo ***'
  READ_NAMELIST( nmmask )

#ifdef OPT_IO_COCOMPI
  call mpi_filopn(mpi_fh, cfmask, 'READ')
  disp=0
#else
  if (myrank .eq. iroot) then
     call filopn(nfmask, cfmask, 'read')
     allocate(g2d(nxygdm))
     allocate(g3d(nxyzgd))
     allocate(g2di(nxygdm))
  else
     allocate(g2d(1))
     allocate(g3d(1))
     allocate(g2di(1))
  end if
#endif

  shift_rdgeo=.true.

  call read_mask_scalar(dx)
  call read_mask_2d(dy,  fact=1.d0, ioff= 0, joff= 0, vinit=1.d0, fill_nsb=.true.)
  call read_mask_2d(dym, fact=1.d0, ioff=-1, joff=-1, vinit=1.d0, fill_nsb=.true.)
  call read_mask_3d(dz,  fact=1.d0, ioff= 0, joff= 0, vinit=1.d0, fill_nsb=.true.)
  call read_mask_3d(dzm, fact=1.d0, ioff= 0, joff= 0, vinit=1.d0, fill_nsb=.true.)
  call read_mask_3d(dzv, fact=1.d0, ioff=-1, joff=-1, vinit=1.d0, fill_nsb=.true.)
  call read_mask_1d(dz0, nzdim)
  call read_mask_1d(ds,  nzdim)
  call read_mask_1d(dsm, nzdim)
  call read_mask_1d(hic, nic+1)
  call read_mask_scalar(rea)
  call read_mask_scalar(zbot)

  call read_mask_2d(  cor, fact= 1.d0, ioff=-1, joff=-1, vinit=1.d0, fill_nsb=.true.)
  call read_mask_2d( dept, fact= 1.d0, ioff= 0, joff= 0, vinit=1.d0, fill_nsb=.true.)
  call read_mask_2d(rdepv, fact= 1.d0, ioff=-1, joff=-1, vinit=1.d0, fill_nsb=.true.)
  call read_mask_2d(  hxt, fact= 1.d0, ioff= 0, joff= 0, vinit=1.d0, fill_nsb=.true.)
  call read_mask_2d(  hxu, fact= 1.d0, ioff=-1, joff=-1, vinit=1.d0, fill_nsb=.true.)
  call read_mask_2d(  hyt, fact= 1.d0, ioff= 0, joff= 0, vinit=1.d0, fill_nsb=.true.)
  call read_mask_2d(  hyu, fact= 1.d0, ioff=-1, joff=-1, vinit=1.d0, fill_nsb=.true.)
  call read_mask_2d( hxyt, fact=-1.d0, ioff= 0, joff= 0, vinit=1.d0, fill_nsb=.true.)
  call read_mask_2d( hxyu, fact=-1.d0, ioff=-1, joff=-1, vinit=1.d0, fill_nsb=.true.)
  call read_mask_2d( hyxt, fact=-1.d0, ioff= 0, joff= 0, vinit=1.d0, fill_nsb=.true.)
  call read_mask_2d( hyxu, fact=-1.d0, ioff=-1, joff=-1, vinit=1.d0, fill_nsb=.true.)

  call read_mask_2i(nbot, vinit=kstr-1)
  call shft1i(      nbot, nxdim, nydim)
#ifdef OPT_TRIPOLE
  call shftint(     nbot, nxdim, nydim)
#endif

#ifdef OPT_BBL
  call read_mask_2i(nbotv, vinit=kstr-1)
  call shft1i(      nbotv, nxdim, nydim)
#ifdef OPT_TRIPOLE
  call shftinv(     nbotv, nxdim, nydim)
  call shft1i(      nbotv, nxdim, nydim)
#endif
#endif
  
  call read_mask_3d(amskt, fact=1.d0, ioff= 0, joff= 0, vinit=0.d0)
  call read_mask_3d(amftx, fact=1.d0, ioff= 1, joff= 0, vinit=0.d0)
  call read_mask_3d(amfty,                              vinit=0.d0)
  call read_mask_3d(amftz, fact=1.d0, ioff= 0, joff= 0, vinit=0.d0)

#ifdef OPT_BBL
  call read_mask_2d(amsktb, fact=1.d0, ioff= 0, joff= 0, vinit=0.d0)
  call read_mask_2d(amskt0, fact=1.d0, ioff= 0, joff= 0, vinit=0.d0)
  call read_mask_2d(amskt1, fact=1.d0, ioff= 0, joff= 0, vinit=0.d0)
#endif

  call read_mask_3d(amskv, fact=1.d0, ioff=-1, joff=-1, vinit=0.d0)
  call read_mask_3d(amfvx, fact=1.d0, ioff= 0, joff=-1, vinit=0.d0)
  call read_mask_3d(amfvy, fact=1.d0, ioff=-1, joff= 0, vinit=0.d0)
  call read_mask_3d(amfvz, fact=1.d0, ioff=-1, joff=-1, vinit=0.d0)

#ifdef OPT_BBL
  call read_mask_2d(amskvb, fact=1.d0, ioff=-1, joff=-1, vinit=0.d0)
  call read_mask_2d(amskv0, fact=1.d0, ioff=-1, joff=-1, vinit=0.d0)
  call read_mask_2d(amskv1, fact=1.d0, ioff=-1, joff=-1, vinit=0.d0)
#endif
  call read_mask_3d(amskb,  fact=1.d0, ioff=-1, joff=-1, vinit=0.d0)

#ifdef OPT_EXMASK
  call read_mask_2d(glont, fact=1.d0, ioff=0, joff=0, fill_nsb=.true.)
  call read_mask_2d(glatt, fact=1.d0, ioff=0, joff=0, fill_nsb=.true.)
  call read_mask_2d(rangt, fact=1.d0, ioff=0, joff=0, fill_nsb=.true.)
#endif

#ifdef OPT_IO_COCOMPI
  call mpi_filcls(mpi_fh)
#else
  if (myrank .eq. iroot) then
     call filcls(nfmask)
  endif
  deallocate(g2d)
  deallocate(g3d)
  deallocate(g2di)
#endif
  
#ifdef OPT_TRIPOLE
  if (jupe .ne. mpi_proc_null.or.jupw .ne. mpi_proc_null) then
     do k=1, nzdim
        do j = jend+1, nydim
           do i = 1, nxdim
              ijs = (j - 2) * nxdim + i
              ij  = (j - 1) * nxdim + i
              amfty(ij,k) = amskt(ij,k)*amskt(ijs,k) * min(dz(ij,k),dz(ijs,k))
           end do
        end do
     end do
  end if
#ifdef OPT_BBL
  if (jupe .ne. mpi_proc_null.or.jupw .ne. mpi_proc_null) then
     do j = jend+1, nydim
        do i = 1, nxdim
           ijs = (j - 2) * nxdim + i
           ij  = (j - 1) * nxdim + i
           amfty(ij,kend) =  amsktb(ij)*amsktb(ijs) * min(dz(ij,kend),dz(ijs,kend))
        end do
     end do
  endif
#endif
#endif
  call shift1(amfty, nxdim,  nydim,  nzdim, 1.d0,  0,  1)
  
  do ij = ijstr, ijend
     ijtstr = ij
     if (amskt(ij, kstr) .eq. 1.d0) exit
  end do
  do ij = ijend, ijstr, -1
     ijtend = ij
     if (amskt(ij, kstr) .eq. 1.d0) exit
  end do

  do ij = ijstr, ijend
     ijvstr = ij
     if (amskv(ij, kstr) .eq. 1.d0) exit
  end do

  do ij = ijend, ijstr, -1
     ijvend = ij
     if (amskv(ij, kstr) .eq. 1.d0) exit
  end do

  rx = 1.d0 / dx
  do ij = 1, nxydim
     ry( ij) = 1.d0 / dy( ij)
     rym(ij) = 1.d0 / dym(ij)
     rxt(ij) = 1.d0 / hxt(ij)
     rxu(ij) = 1.d0 / hxu(ij)
     ryt(ij) = 1.d0 / hyt(ij)
     ryu(ij) = 1.d0 / hyu(ij)
  end do
  do k = 1, nzdim
     rs( k) = 1.d0 / ds( k)
     rsm(k) = 1.d0 / dsm(k)
  end do

#ifdef OPT_OFFLINE  
#ifdef OPT_BBL
  call rmmskv
  call admkvb
#endif
#endif
!$acc enter data copyin(amskt, amftx, amfty, amftz,  amskv,  amfvx,  amfvy,  amfvz,  amskb,  nbot)
!$acc enter data copyin(dy, dym, dz, dzm, dzv, dz0, ds, dsm, hic, cor, dept, rdepv)
!$acc enter data copyin(hxt, hxu, hyt, hyu, hxyt, hxyu, hyxt, hyxu)
!$acc enter data copyin(ry, rym, rxt, rxu, ryt, ryu, rs, rsm)
#ifdef OPT_EXMASK
!$acc enter data copyin(glont, glatt, rangt)
#endif
  shift_rdgeo=.false.
  return
contains
  subroutine read_mask_scalar(d0)
    implicit none
    real(8), intent(out) :: d0
    real(8) :: d1(1)
#ifdef OPT_IO_COCOMPI
    call mpi_read_root(d1,1,mpi_fh, disp)
    d0=d1(1)
#else
    if (myrank .eq. iroot) then
       read(nfmask)d0
    end if
#endif
    call mpi_bcast(d0, 1, mpi_real8, iroot, mpi_comm_ogcm, ierr)
  end subroutine read_mask_scalar

  subroutine read_mask_1d(d1, ndim)
    implicit none
    integer, intent(in)  :: ndim
    real(8), intent(out) :: d1(ndim)
#ifdef OPT_IO_COCOMPI
    call mpi_read_root(d1, ndim, mpi_fh, disp)
#else
    if (myrank .eq. iroot) then
       read(nfmask)d1
    end if
#endif
    call mpi_bcast(d1, ndim, mpi_real8, iroot, mpi_comm_ogcm, ierr)
  end subroutine read_mask_1d
  
  subroutine read_mask_2i(i2, vinit)
    implicit none
    integer, intent(out) :: i2(nxdim, nydim)
    integer, intent(in), optional :: vinit
    if( present(vinit) ) then
       i2(:,:)=vinit
    end if
#ifdef OPT_IO_COCOMPI
    call mpi_read_2d_intx( i2, mpi_fh, disp)
#else
  if (myrank .eq. iroot) then
     read(nfmask) g2di
  end if
  call scatter_2d_int(i2, g2di)
#endif
  end subroutine read_mask_2i
  
  subroutine read_mask_2d(d2, fact, ioff, joff, vinit, fill_nsb)
    implicit none
    real(8), intent(out) :: d2(nxdim, nydim)
    integer, intent(in), optional :: ioff, joff
    real(8), intent(in), optional :: fact
    real(8), intent(in), optional :: vinit
    logical, intent(in), optional :: fill_nsb
    if( present(vinit) ) then
       d2(:,:)=vinit
    end if
#ifdef OPT_IO_COCOMPI
    call mpi_read_2d_dimx(d2, mpi_fh, disp)
#else
    if (myrank .eq. iroot) then
       read(nfmask) g2d
    end if
    call scatter_2d(d2, g2d)
#endif
    if( present(fact) ) then
       call shift1(d2, nxdim, nydim, 1, fact, ioff, joff)
    end if

    if( present(fill_nsb) ) then
       if(.not. fill_nsb) return
#ifndef OPT_TRIPOLE
       if (jup .eq. mpi_proc_null) then
          do j = jend+1, nydim
             do i=1, nxdim
                d2(i,j)=d2(i,jend)
             end do
          end do
       end if
#endif
       if (jdown .eq. mpi_proc_null) then
          do j = 1, jstr-1
             do i = 1, nxdim
                d2(i,j)=d2(i,jstr)
             end do
          end do
       end if
    end if
  end subroutine read_mask_2d

  subroutine read_mask_3d(d3, fact, ioff, joff, vinit, fill_nsb)
    implicit none
    real(8), intent(out) :: d3(nxdim, nydim, nzdim)
    integer, intent(in), optional :: ioff, joff
    real(8), intent(in), optional :: fact
    real(8), intent(in), optional :: vinit
    logical, intent(in), optional :: fill_nsb
    
    if( present(vinit) ) then
       d3(:,:,:)=vinit
    end if
#ifdef OPT_IO_COCOMPI
    call mpi_read_3d_dimx(d3, mpi_fh, disp)
#else
    if (myrank .eq. iroot) then
       read(nfmask) g3d
    end if
    call scatter_3d(d3, g3d)
#endif
    if( present(fact) ) then
       call shift1(d3, nxdim, nydim, nzdim, fact, ioff, joff)
    end if

    if( present(fill_nsb) ) then
       if(.not. fill_nsb) return
#ifndef OPT_TRIPOLE
       if (jup .eq. mpi_proc_null) then
          do k = 1, nzdim
             do j = jend+1, nydim
                do i=1, nxdim
                   d3(i,j,k)=d3(i,jend,k)
                end do
             end do
          end do
       end if
#endif
       if (jdown .eq. mpi_proc_null) then
          do k = 1, nzdim
             do j = 1, jstr-1
                do i = 1, nxdim
                   d3(i,j,k)=d3(i,jstr,k)
                end do
             end do
          end do
       end if
    end if
  end subroutine read_mask_3d
  
  end subroutine rdgeo
end module brdge
