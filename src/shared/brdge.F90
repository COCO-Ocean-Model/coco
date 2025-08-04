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
!
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
 use bgs2d
 use bgs3d
 use bshfi
 use bshft
#ifdef OPT_OFFLINE 
 use bchmk
#endif
 implicit none
#include "mpif.h"


  integer,allocatable :: g2di(:)
  real(8),allocatable :: g2d(:), g3d(:)
  integer ::  n
  integer :: ij,      k
  integer ::  i,      j,    ijs,   ije,   ijn,  ijw
  integer ::  ifpar,  jfpar
  integer :: nfmask
  integer :: istat

  character(ncf) :: cfmask
  data cfmask / 'not-specified' /
  namelist /nmmask/ cfmask

  call rewnml(ifpar, jfpar)
  write(jfpar, *) '*** rdgeo ***'
  read(ifpar, nmmask, iostat=istat)
  call cstnml(jfpar, 'rdgeo', 'nmmask', istat)

  if (myrank .eq. iroot) then
    call filopn(nfmask, cfmask, 'read')
  endif

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


  if (myrank .eq. iroot) then
     allocate(g2d(nxygdm))
     allocate(g3d(nxyzgd))
     allocate(g2di(nxygdm))
  else
     allocate(g2d(1))
     allocate(g3d(1))
     allocate(g2di(1))
  end if

  if (myrank .eq. iroot) then
     read(nfmask) dx
  end if
  call mpi_bcast(dx, 1, mpi_real8, iroot, mpi_comm_ogcm, ierr)

  if (myrank .eq. iroot) then
     read(nfmask) g2d
  end if
  call scatter_2d(dy, g2d)

  if (myrank .eq. iroot) then
     read(nfmask) g2d
  end if
  call scatter_2d(dym, g2d)

  if (myrank .eq. iroot) then
     read(nfmask) g3d
  end if
  call scatter_3d(dz, g3d)

  if (myrank .eq. iroot) then
     read(nfmask) g3d
  end if
  call scatter_3d(dzm, g3d)

  if (myrank .eq. iroot) then
     read(nfmask) g3d
  end if
  call scatter_3d(dzv, g3d)

  if (myrank .eq. iroot) then
     read(nfmask) dz0
  end if
  call mpi_bcast(dz0, nzdim, mpi_real8, iroot, mpi_comm_ogcm, ierr)

  if (myrank .eq. iroot) then
     read(nfmask) ds
  end if
  call mpi_bcast(ds, nzdim, mpi_real8, iroot, mpi_comm_ogcm, ierr)

  if (myrank .eq. iroot) then
     read(nfmask) dsm
  end if
  call mpi_bcast(dsm, nzdim, mpi_real8, iroot, mpi_comm_ogcm, ierr)

  if (myrank .eq. iroot) then
     read(nfmask) hic
  end if
  call mpi_bcast(hic, nic+1, mpi_real8, iroot, mpi_comm_ogcm, ierr)

  if (myrank .eq. iroot) then
     read(nfmask) rea
  end if
  call mpi_bcast(rea, 1, mpi_real8, iroot, mpi_comm_ogcm, ierr)

  if (myrank .eq. iroot) then
     read(nfmask) zbot
  end if
  call mpi_bcast(zbot, 1, mpi_real8, iroot, mpi_comm_ogcm, ierr)

  if (myrank .eq. iroot) then
     read(nfmask) g2d
  end if
  call scatter_2d(cor, g2d)

  if (myrank .eq. iroot) then
     read(nfmask) g2d
  end if
  call scatter_2d(dept, g2d)

  if (myrank .eq. iroot) then
     read(nfmask) g2d
  end if
  call scatter_2d(rdepv, g2d)

  if (myrank .eq. iroot) then
     read(nfmask) g2d
  end if
  call scatter_2d(hxt, g2d)

  if (myrank .eq. iroot) then
     read(nfmask) g2d
  end if
  call scatter_2d(hxu, g2d)

  if (myrank .eq. iroot) then
     read(nfmask) g2d
  end if
  call scatter_2d(hyt, g2d)

  if (myrank .eq. iroot) then
     read(nfmask) g2d
  end if
  call scatter_2d(hyu, g2d)

  if (myrank .eq. iroot) then
     read(nfmask) g2d
  end if
  call scatter_2d(hxyt, g2d)

  if (myrank .eq. iroot) then
     read(nfmask) g2d
  end if
  call scatter_2d(hxyu, g2d)

  if (myrank .eq. iroot) then
     read(nfmask) g2d
  end if
  call scatter_2d(hyxt, g2d)

  if (myrank .eq. iroot) then
     read(nfmask) g2d
  end if
  call scatter_2d(hyxu, g2d)

#ifdef OPT_TRIPOLE
  call shift2(    dy,     dept,                                                &
   &             nxdim,  nydim,      1,                                        &
   &              1.d0,      0,      0 )
  call shift3(   dym,    cor,  rdepv,                                          &
   &             nxdim,  nydim,      1,                                        &
   &              1.d0,     -1,     -1 )
  call shift2(    dz,    dzm,                                                  &
   &             nxdim,  nydim,  nzdim,                                        &
   &              1.d0,      0,      0 )
  call shift1(   dzv,                                                          &
   &             nxdim,  nydim,  nzdim,                                        &
   &              1.d0,     -1,     -1 )
  call shift2(   hxt,    hyt,                                                  &
   &             nxdim,  nydim,      1,                                        &
   &              1.d0,      0,      0 )
  call shift2(   hxu,    hyu,                                                  &
   &             nxdim,  nydim,      1,                                        &
   &              1.d0,     -1,     -1 )
  call shift2(  hxyt,   hyxt,                                                  &
   &             nxdim,  nydim,      1,                                        &
   &             -1.d0,      0,      0 )
  call shift2(  hxyu,   hyxu,                                                  &
   &             nxdim,  nydim,      1,                                        &
   &             -1.d0,     -1,     -1 )
#else
  call shift2(                                                                 &
   &                dy,    dym,                                                &
   &             nxdim,  nydim,      1)
  call shift3(                                                                 &
   &                dz,    dzm,    dzv,                                        &
   &             nxdim,  nydim,  nzdim)
  call shift3(                                                                 &
   &               cor,   dept,  rdepv,                                        &
   &             nxdim,  nydim,      1)
  call shift3(                                                                 &
   &               hxt,    hxu,    hyt,                                        &
   &             nxdim,  nydim,      1)
  call shift3(                                                                 &
   &               hyu,   hxyt,   hxyu,                                        &
   &             nxdim,  nydim,      1)
  call shift2(                                                                 &
   &              hyxt,   hyxu,                                                &
   &             nxdim,  nydim,      1)

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

  if (myrank .eq. iroot) then
     read(nfmask) g2di
  end if
  call scatter_2d_int(nbot, g2di)

#ifdef OPT_BBL
  if (myrank .eq. iroot) then
     read(nfmask) g2di
  end if
  call scatter_2d_int(nbotv, g2di)
#endif

  if (myrank .eq. iroot) then
     read(nfmask) g3d
  end if
  call scatter_3d(amskt, g3d)

  if (myrank .eq. iroot) then
     read(nfmask) g3d
  end if
  call scatter_3d(amftx, g3d)

  if (myrank .eq. iroot) then
     read(nfmask) g3d
  end if
  call scatter_3d(amfty, g3d)

  if (myrank .eq. iroot) then
     read(nfmask) g3d
  end if
  call scatter_3d(amftz, g3d)

#ifdef OPT_BBL
  if (myrank .eq. iroot) then
     read(nfmask) g2d
  end if
  call scatter_2d(amsktb, g2d)

  if (myrank .eq. iroot) then
     read(nfmask) g2d
  end if
  call scatter_2d(amskt0, g2d)

  if (myrank .eq. iroot) then
     read(nfmask) g2d
  end if
  call scatter_2d(amskt1, g2d)
#endif

  if (myrank .eq. iroot) then
     read(nfmask) g3d
  end if
  call scatter_3d(amskv, g3d)

  if (myrank .eq. iroot) then
     read(nfmask) g3d
  end if
  call scatter_3d(amfvx, g3d)

  if (myrank .eq. iroot) then
     read(nfmask) g3d
  end if
  call scatter_3d(amfvy, g3d)

  if (myrank .eq. iroot) then
     read(nfmask) g3d
  end if
  call scatter_3d(amfvz, g3d)

#ifdef OPT_BBL
  if (myrank .eq. iroot) then
     read(nfmask) g2d
  end if
  call scatter_2d(amskvb, g2d)

  if (myrank .eq. iroot) then
     read(nfmask) g2d
  end if
  call scatter_2d(amskv0, g2d)

  if (myrank .eq. iroot) then
     read(nfmask) g2d
  end if
  call scatter_2d(amskv1, g2d)
#endif

  if (myrank .eq. iroot) then
     read(nfmask) g3d
  endif
  call scatter_3d(amskb, g3d)

#ifdef OPT_TRIPOLE
  call shift2(                                                                 &
   &             amskt,  amftz,                                                &
   &             nxdim,  nydim,  nzdim,                                        &
   &              1.d0,      0,      0 )
  call shift2(                                                                 &
   &             amskv,  amfvz,                                                &
   &             nxdim,  nydim,  nzdim,                                        &
   &              1.d0,     -1,     -1 )
  call shift1(                                                                 &
   &             amskb,                                                        &
   &             nxdim,  nydim,  nzdim,                                        &
   &              1.d0,     -1,     -1 )

  call shft1i(                                                                 &
   &              nbot,                                                        &
   &             nxdim,  nydim)
  call shftint(                                                                &
   &              nbot,                                                        &
   &             nxdim,  nydim)

  if (jupe .ne. mpi_proc_null.or.jupw .ne. mpi_proc_null) then
!      do j = jend+1, nydim
!         do i = 2, nxdim
!            ijw = (j - 1) * nxdim + i  - 1
!            ij  = (j - 1) * nxdim + i
!            do k=1,nzdim
!               amftx(ij,k) = amskt(ij,k)*amskt(ijw,k)
!     &                     * min(dz(ij,k),dz(ijw,k))
!             enddo
!         enddo
!      enddo
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
!      do j = jend+1, nydim-1
!         do i = 2, nxdim
!            ijn = j * nxdim + i
!            ijw = (j - 1) * nxdim + i  - 1
!            ij  = (j - 1) * nxdim + i
!            do k= 1, kstr+kz-1
!                amfvx(ij,k) = amskt(ij,k)*amskt(ijn,k)
!            enddo
!            do k=kstr+kz,nzdim
!                amfvx(ij,k) = amskt(ij,k)*amskt(ijn,k)
!     &                      * min(dzv(ij,k),dzv(ijw,k))
!             enddo
!         enddo
!      enddo
!      do j = jend+1, nydim
!         do i = 1, nxdim-1
!            ijs = (j - 2) * nxdim + i
!            ije = (j - 1) * nxdim + i  + 1
!            ij  = (j - 1) * nxdim + i
!            do k= 1, kstr+kz-1
!                amfvy(ij,k) = amskt(ij,k)*amskt(ije,k)
!            enddo
!            do k= kstr+kz,nzdim
!                amfvy(ij,k) = amskt(ij,k)*amskt(ije,k)
!     &                      * min(dzv(ij,k),dzv(ijs,k))
!            enddo
!         enddo
!      enddo
  endif

  call shift1(                                                                 &
   &             amftx,                                                        &
   &             nxdim,  nydim,  nzdim,                                        &
   &             1.d0 ,      1,      0)
  call shift1(                                                                 &
   &             amfty,                                                        &
   &             nxdim,  nydim,  nzdim,                                        &
   &             1.d0 ,      0,      1)
  call shift1(                                                                 &
   &             amfvx,                                                        &
   &             nxdim,  nydim,  nzdim,                                        &
   &             1.d0 ,      0,     -1)
  call shift1(                                                                 &
   &             amfvy,                                                        &
   &             nxdim,  nydim,  nzdim,                                        &
   &             1.d0 ,     -1,      0)
#else
  call shift3(                                                                 &
   &             amskt,  amskv,  amskb,                                        &
   &             nxdim,  nydim,  nzdim)
  call shift3(                                                                 &
   &             amftx,  amfty,  amftz,                                        &
   &             nxdim,  nydim,  nzdim)
  call shift3(                                                                 &
   &             amfvx,  amfvy,  amfvz,                                        &
   &             nxdim,  nydim,  nzdim)
  call shft1i(                                                                 &
   &              nbot,                                                        &
   &             nxdim,  nydim)
#endif

#ifdef OPT_BBL
#ifdef OPT_TRIPOLE
  call shift3(                                                                 &
   &            amsktb, amskt0, amskt1,                                        &
   &             nxdim,  nydim,      1,                                        &
   &              1.d0,      0,      0 )
  call shift3(                                                                 &
   &            amskvb, amskv0, amskv1,                                        &
   &             nxdim,  nydim,      1,                                        &
   &              1.d0,     -1,     -1 )
  call shft1i(                                                                 &
   &             nbotv,                                                        &
   &             nxdim,  nydim)
  call shftinv(                                                                &
   &             nbotv,                                                        &
   &             nxdim,  nydim)
  call shft1i(                                                                 &
   &             nbotv,                                                        &
   &             nxdim,  nydim)

  if (jupe .ne. mpi_proc_null.or.jupw .ne. mpi_proc_null) then
!      do j = jend+1, nydim
!         do i = 2, nxdim
!            ijs = (j - 2) * nxdim + i
!            ij  = (j - 1) * nxdim + i
!            amftx(ij,kend) =  amsktb(ij)*amsktb(ijw)
!     &                      * min(dz(ij,kend),dz(ijw,kend))
!         enddo
!      enddo
  do j = jend+1, nydim
     do i = 1, nxdim
        ijs = (j - 2) * nxdim + i
        ij  = (j - 1) * nxdim + i
        amfty(ij,kend) =  amsktb(ij)*amsktb(ijs)                               &
   &                    * min(dz(ij,kend),dz(ijs,kend))
     enddo
  enddo
!      do j = jend+1, nydim-1
!         do i = 2, nxdim
!            ijn = j * nxdim + i
!            ijw = (j - 1) * nxdim + i  - 1
!            ij  = (j - 1) * nxdim + i
!            amfvx(ij,kend) =  amsktb(ij)*amsktb(ijn)
!     &                      * min(dzv(ij,kend),dzv(ijw,kend))
!         enddo
!      enddo
!      do j = jend+1, nydim
!         do i = 1, nxdim-1
!            ijs = (j - 2) * nxdim + i
!            ije = (j - 1) * nxdim + i  + 1
!            ij  = (j - 1) * nxdim + i
!            amfvy(ij,kend) =  amsktb(ij)*amsktb(ije)
!     &                      * min(dzv(ij,kend),dzv(ijs,kend))
!         enddo
!      enddo
  endif
  call shift1(                                                                 &
   &             amftx,                                                        &
   &             nxdim,  nydim,  nzdim,                                        &
   &             1.d0 ,      1,      0)
  call shift1(                                                                 &
   &             amfty,                                                        &
   &             nxdim,  nydim,  nzdim,                                        &
   &             1.d0 ,      0,      1)
  call shift1(                                                                 &
   &             amfvx,                                                        &
   &             nxdim,  nydim,  nzdim,                                        &
   &             1.d0 ,      0,     -1)
  call shift1(                                                                 &
   &             amfvy,                                                        &
   &             nxdim,  nydim,  nzdim,                                        &
   &             1.d0 ,     -1,      0)
#else
  call shft1i(                                                                 &
   &             nbotv,                                                        &
   &             nxdim,  nydim)
  call shift3(                                                                 &
   &            amsktb, amskt0, amskt1,                                        &
   &             nxdim,  nydim,      1)
  call shift3(                                                                 &
   &            amskvb, amskv0, amskv1,                                        &
   &             nxdim,  nydim,      1)
#endif
#endif

#ifdef OPT_EXMASK
  if(myrank .eq. iroot) THEN
     read(nfmask) g2d
  end if
  call scatter_2d(glont, g2d)

  if(myrank .eq. iroot) THEN
     read(nfmask) g2d
  end if
  call scatter_2d(glatt, g2d)

  if(myrank .eq. iroot) THEN
     read(nfmask) g2d
  end if
  call scatter_2d(rangt, g2d)

#ifdef OPT_TRIPOLE
  call shift3(                                                                 &
   &           glont,  glatt,  rangt,                                          &
   &           nxdim,  nydim,      1,                                          &
   &            1.d0,      0,      0 )
#else
  call shift3(                                                                 &
   &           glont,  glatt,  rangt,                                          &
   &           nxdim,  nydim,      1)
  
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

  if (myrank .eq. iroot) then
     call filcls(nfmask)
  endif
  deallocate(g2d)
  deallocate(g3d)
  deallocate(g2di)

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

end subroutine rdgeo

end module brdge
