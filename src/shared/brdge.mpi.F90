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
 use mpiio
#ifdef OPT_OFFLINE 
 use bchmk
#endif

 implicit none
#include "mpif.h"

  integer ::  n
  integer :: ij,      k
  integer ::  i,      j,    ijs,   ije,   ijn,  ijw
  integer ::  ifpar,  jfpar
  integer :: istat
  integer :: mpi_fh
  integer (kind = mpi_offset_kind) :: disp
  real(8) :: dim1(1)
  character(ncf) :: cfmask
  data cfmask / 'not-specified' /
  namelist /nmmask/ cfmask

  call rewnml(ifpar, jfpar)
  write(jfpar, *) '*** rdgeo ***'
  read(ifpar, nmmask, iostat=istat)
  call cstnml(jfpar, 'rdgeo', 'nmmask', istat)

  call mpi_filopn(mpi_fh, cfmask, 'READ')

  shift_rdgeo=.true.
  
  amskt = 0.d0
  amftx = 0.d0
  amfty = 0.d0
  amftz = 0.d0
  amskv = 0.d0
  amfvx = 0.d0
  amfvy = 0.d0
  amfvz = 0.d0
  amskb = 0.d0
  dy = 1.d0
  dym = 1.d0
  dz  = 1.d0
  dzm = 1.d0
  dzv = 1.d0
  cor = 1.d0
  hxt = 1.d0
  hxu = 1.d0
  hyt = 1.d0
  hyu = 1.d0
  hxyt = 1.d0
  hxyu = 1.d0
  hyxt = 1.d0
  hyxu = 1.d0
  dept = 1.d0
  rdepv = 1.d0
  nbot = kstr - 1
#ifdef OPT_BBL
  nbotv = kstr - 1
  amsktb = 0.d0
  amskt0 = 0.d0
  amskt1 = 0.d0
  amskvb = 0.d0
  amskv0 = 0.d0
  amskv1 = 0.d0
#endif

  disp=0
  call mpi_read_root(dim1,1,mpi_fh, disp)
  dx=dim1(1)
  call mpi_bcast(dx, 1, mpi_real8, iroot, mpi_comm_ogcm, ierr)

  call mpi_read_2d_dimx(dy , mpi_fh, disp)
  call mpi_read_2d_dimx(dym, mpi_fh, disp)
  call mpi_read_3d_dimx(dz , mpi_fh, disp)
  call mpi_read_3d_dimx(dzm, mpi_fh, disp)
  call mpi_read_3d_dimx(dzv, mpi_fh, disp)

  call mpi_read_root(dz0, nzdim, mpi_fh, disp)
  call mpi_bcast(dz0, nzdim, mpi_real8, iroot, mpi_comm_ogcm, ierr)

  call mpi_read_root(ds, nzdim, mpi_fh, disp)
  call mpi_bcast(ds, nzdim, mpi_real8, iroot, mpi_comm_ogcm, ierr)

  call mpi_read_root(dsm, nzdim, mpi_fh, disp)
  call mpi_bcast(dsm, nzdim, mpi_real8, iroot, mpi_comm_ogcm, ierr)

  call mpi_read_root(hic, nic+1, mpi_fh, disp)
  call mpi_bcast(hic, nic+1, mpi_real8, iroot, mpi_comm_ogcm, ierr)

  call mpi_read_root(dim1, 1, mpi_fh, disp)
  rea=dim1(1)
  call mpi_bcast(rea, 1, mpi_real8, iroot, mpi_comm_ogcm, ierr)

  call mpi_read_root(dim1, 1, mpi_fh, disp)
  zbot=dim1(1)
  call mpi_bcast(zbot, 1, mpi_real8, iroot, mpi_comm_ogcm, ierr)

  call mpi_read_2d_dimx(  cor, mpi_fh, disp)
  call mpi_read_2d_dimx( dept, mpi_fh, disp)
  call mpi_read_2d_dimx(rdepv, mpi_fh, disp)
  call mpi_read_2d_dimx(  hxt, mpi_fh, disp)
  call mpi_read_2d_dimx(  hxu, mpi_fh, disp)
  call mpi_read_2d_dimx(  hyt, mpi_fh, disp)
  call mpi_read_2d_dimx(  hyu, mpi_fh, disp)
  call mpi_read_2d_dimx( hxyt, mpi_fh, disp)
  call mpi_read_2d_dimx( hxyu, mpi_fh, disp)
  call mpi_read_2d_dimx( hyxt, mpi_fh, disp)
  call mpi_read_2d_dimx( hyxu, mpi_fh, disp)

  call shift2(  dy, dept,        nxdim, nydim,     1,  1.d0,  0,  0)
  call shift3( dym,  cor, rdepv, nxdim, nydim,     1,  1.d0, -1, -1)
  call shift2(  dz,  dzm,        nxdim, nydim, nzdim,  1.d0,  0,  0)
  call shift1( dzv,              nxdim, nydim, nzdim,  1.d0, -1, -1)
  call shift2( hxt,  hyt,        nxdim, nydim,     1,  1.d0,  0,  0)
  call shift2( hxu,  hyu,        nxdim, nydim,     1,  1.d0, -1, -1)
  call shift2(hxyt, hyxt,        nxdim, nydim,     1, -1.d0,  0,  0)
  call shift2(hxyu, hyxu,        nxdim, nydim,     1, -1.d0, -1, -1)

#ifndef OPT_TRIPOLE
  if (jup .eq. mpi_proc_null) then
     do j = jend+1, nydim
        do i = 1, nxdim
           ije = (jend - 1) * nxdim + i
           ij  = (j - 1) * nxdim + i
           dy(ij) = dy(ije)
           dym(ij) = dym(ije)
           cor(ij) = cor(ije)
           dept(ij) = dept(ije)
           rdepv(ij) = rdepv(ije)
           hxt(ij) = hxt(ije)
           hxu(ij) = hxu(ije)
           hyt(ij) = hyt(ije)
           hyu(ij) = hyu(ije)
           hxyt(ij) = hxyt(ije)
           hxyu(ij) = hxyu(ije)
           hyxt(ij) = hyxt(ije)
           hyxu(ij) = hyxu(ije)
        end do
     end do
     do k = 1, nzdim
        do j = jend+1, nydim
           do i = 1, nxdim
              ije = (jend - 1) * nxdim + i
              ij  = (j - 1) * nxdim + i
              dz(ij, k) = dz(ije, k)
              dzm(ij, k) = dzm(ije, k)
              dzv(ij, k) = dzv(ije, k)
           end do
        end do
     end do
  end if
#endif
  if (jdown .eq. mpi_proc_null) then
     do j = 1, jstr-1
        do i = 1, nxdim
           ijs = (jstr - 1) * nxdim + i
           ij  = (j - 1) * nxdim + i
           dy(ij) = dy(ijs)
           dym(ij) = dym(ijs)
           cor(ij) = cor(ijs)
           dept(ij) = dept(ijs)
           rdepv(ij) = rdepv(ijs)
           hxt(ij) = hxt(ijs)
           hxu(ij) = hxu(ijs)
           hyt(ij) = hyt(ijs)
           hyu(ij) = hyu(ijs)
           hxyt(ij) = hxyt(ijs)
           hxyu(ij) = hxyu(ijs)
           hyxt(ij) = hyxt(ijs)
           hyxu(ij) = hyxu(ijs)
        end do
     end do
     do k = 1, nzdim
        do j = 1, jstr-1
           do i = 1, nxdim
              ijs = (jstr - 1) * nxdim + i
              ij  = (j - 1) * nxdim + i
              dz(ij, k) = dz(ijs, k)
              dzm(ij, k) = dzm(ijs, k)
              dzv(ij, k) = dzv(ijs, k)
           end do
        end do
     end do
  end if

  call mpi_read_2d_intx( nbot, mpi_fh, disp)
#ifdef OPT_BBL
  call mpi_read_2d_intx(nbotv, mpi_fh, disp)
#endif

  call mpi_read_3d_dimx(amskt, mpi_fh, disp)
  call mpi_read_3d_dimx(amftx, mpi_fh, disp)
  call mpi_read_3d_dimx(amfty, mpi_fh, disp)
  call mpi_read_3d_dimx(amftz, mpi_fh, disp)

#ifdef OPT_BBL
  call mpi_read_2d_dimx(amsktb, mpi_fh, disp)
  call mpi_read_2d_dimx(amskt0, mpi_fh, disp)
  call mpi_read_2d_dimx(amskt1, mpi_fh, disp)
#endif

  call mpi_read_3d_dimx(amskv, mpi_fh, disp)
  call mpi_read_3d_dimx(amfvx, mpi_fh, disp)
  call mpi_read_3d_dimx(amfvy, mpi_fh, disp)
  call mpi_read_3d_dimx(amfvz, mpi_fh, disp)

#ifdef OPT_BBL
  call mpi_read_2d_dimx(amskvb, mpi_fh, disp)
  call mpi_read_2d_dimx(amskv0, mpi_fh, disp)
  call mpi_read_2d_dimx(amskv1, mpi_fh, disp)
#endif
  call mpi_read_3d_dimx(amskb, mpi_fh, disp)

  call shift2(amskt, amftz, nxdim, nydim, nzdim, 1.d0,  0,  0)
  call shift2(amskv, amfvz, nxdim, nydim, nzdim, 1.d0, -1, -1)
  call shift1(amskb,        nxdim, nydim, nzdim, 1.d0, -1, -1)
  call shft1i( nbot, nxdim, nydim)

#ifdef OPT_TRIPOLE
  call shftint(nbot, nxdim, nydim)

  if (jupe .ne. mpi_proc_null.or.jupw .ne. mpi_proc_null) then
  do j = jend+1, nydim
     do i = 1, nxdim
        ijs = (j - 2) * nxdim + i
        ij  = (j - 1) * nxdim + i
        do k=1,nzdim
            amfty(ij,k) = amskt(ij,k)*amskt(ijs,k)                             &
   &                    * min(dz(ij,k),dz(ijs,k))
         enddo
     enddo
  enddo
  endif
#endif
  call shift1(amftx, nxdim,  nydim,  nzdim, 1.d0,  1,  0)
  call shift1(amfty, nxdim,  nydim,  nzdim, 1.d0,  0,  1)
  call shift1(amfvx, nxdim,  nydim,  nzdim, 1.d0,  0, -1)
  call shift1(amfvy, nxdim,  nydim,  nzdim, 1.d0, -1,  0)

#ifdef OPT_BBL
  call shift3(amsktb, amskt0, amskt1, nxdim, nydim, 1, 1.d0,  0,  0)
  call shift3(amskvb, amskv0, amskv1, nxdim, nydim, 1, 1.d0, -1, -1)
  call shft1i( nbotv, nxdim, nydim)

#ifdef OPT_TRIPOLE
  call shftinv(nbotv, nxdim, nydim)
  call shft1i( nbotv, nxdim, nydim)

  if (jupe .ne. mpi_proc_null.or.jupw .ne. mpi_proc_null) then
  do j = jend+1, nydim
     do i = 1, nxdim
        ijs = (j - 2) * nxdim + i
        ij  = (j - 1) * nxdim + i
        amfty(ij,kend) =  amsktb(ij)*amsktb(ijs)                               &
   &                    * min(dz(ij,kend),dz(ijs,kend))
     enddo
  enddo
  endif
  call shift1(amftx, nxdim,  nydim,  nzdim, 1.d0,  1,  0)
  call shift1(amfty, nxdim,  nydim,  nzdim, 1.d0,  0,  1)
  call shift1(amfvx, nxdim,  nydim,  nzdim, 1.d0,  0, -1)
  call shift1(amfvy, nxdim,  nydim,  nzdim, 1.d0, -1,  0)
#endif
#endif

#ifdef OPT_EXMASK
  call mpi_read_2d_dimx(glont, mpi_fh, disp)
  call mpi_read_2d_dimx(glatt, mpi_fh, disp)
  call mpi_read_2d_dimx(rangt, mpi_fh, disp)

  call shift3(glont, glatt, rangt, nxdim, nydim, 1, 1.d0, 0, 0)

#ifndef OPT_TRIPOLE
  if (jup .eq. mpi_proc_null) then
     do j = jend+1, nydim
        do i = 1, nxdim
           ije = (jend - 1) * nxdim + i
           ij  = (j - 1) * nxdim + i
           glont(ij) = glont(ije)
           glatt(ij) = glatt(ije)
           rangt(ij) = rangt(ije)
        end do
     end do
  end if
#endif

  if (jdown .eq. mpi_proc_null) then
     do j = 1, jstr-1
        do i = 1, nxdim
           ijs = (jstr - 1) * nxdim + i
           ij  = (j - 1) * nxdim + i
           glont(ij) = glont(ijs)
           glatt(ij) = glatt(ijs)
           rangt(ij) = rangt(ijs)
        end do
     end do
  end if
#endif


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

#ifdef OPT_OFFLINE
#ifdef OPT_BBL
  call rmmskv
  call admkvb
#endif
#endif

  call mpi_filcls(mpi_fh)

  rx = 1.d0 / dx
  do ij = 1, nxydim
     ry(ij) = 1.d0 / dy(ij)
     rym(ij) = 1.d0 / dym(ij)
     rxt(ij) = 1.d0 / hxt(ij)
     rxu(ij) = 1.d0 / hxu(ij)
     ryt(ij) = 1.d0 / hyt(ij)
     ryu(ij) = 1.d0 / hyu(ij)
  end do
  do k = 1, nzdim
     rs(k) = 1.d0 / ds(k)
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
  end subroutine rdgeo
end module brdge
