module ctnuv

! --- information -----------------------------------------------------
!
!  Estimate the pressure gradient and Coriolis terms of the equation
! of motion.
!
!  HISTORY
!     '03.04.22  H.Hasumi: from COCO3.4
!     '07.04.23  H.Hasumi
!     '07.05.01  H.Hasumi: for McDougall et al. (2003) eq. of state
!     '08.06.11  H.Hasumi: initial/final processing
!     '08.07.10  H.Hasumi: initial/final processing
!     '09.01.15  T.Suzuki: bug fix vertical integration thanks to Hasumi
!     '10.04.14  M.Kurogi: staggered time stepping
!     '12.06.14  Y.Komuro: for COCO5.0
!
! ---------------------------------------------------------------------

  use zocdim, only: &
    &  nxdim,  nydim, nxydim,  nzdim,  ntdim, &
    &   kstr,   kend,     kz, &
    & ijvstr, ijvend, &
    &     le,     lw,     ln,     ls,    lne, &
    &  oinit, ofinal
  use zocgrd, only: &
    &     dy,     dz,    dzv,    dz0,     ds, &
    &     rx,    rym,   zbot,    cor, &
    &     hxu,   hyu,   hxyu,   hyxu,    rxu,    ryu 
  use zocmsk, only: &
#ifdef OPT_BBL
    & amskvb, nbotv, &
#endif
    &  amskv,   nbot
  use zocphy, only: &
    & gravit,   rhoo
  use brstt

  implicit none
  private

  real(8), save :: px(nxydim, nzdim),    pxm(nxydim, nzdim)
  real(8), save :: py(nxydim, nzdim),    pym(nxydim, nzdim)

  public :: tnduvd
#ifdef OPT_BBL
  public :: tnduvb
#endif

contains 

subroutine tnduvd( &
  &    gxx,    gyy,     gx,     gy, &
  &     xx,     yy,     uy,     vy, &
  &     hy,      r,     ux,     vx, & 
  &   uadv,   vadv )

  real(8), intent(out)   ::    gxx(nxydim),           gyy(nxydim)
  real(8), intent(inout) ::     gx(nxydim, nzdim),     gy(nxydim, nzdim)
  real(8), intent(inout) ::     xx(nxydim, nzdim),     yy(nxydim, nzdim)
  real(8), intent(in)    ::     uy(nxydim, nzdim),     vy(nxydim, nzdim)
  real(8), intent(in)    ::      r(nxydim, nzdim),     hy(nxydim)
  real(8), intent(in)    ::     ux(nxydim, nzdim),     vx(nxydim, nzdim)
  real(8), intent(in)    ::   uadv(nxydim, nzdim),   vadv(nxydim, nzdim)

  real(8), save ::     cof

  real(8) ::      uu(nxydim)       ,     uv(nxydim)
  real(8) ::      vu(nxydim)       ,     vv(nxydim)
  real(8) ::   hvbot(nxydim)
  real(8) ::   dzsig(nxydim, nzdim)

  integer ::     ij,      k
  logical, save :: ofirst = .true.

! common /work/ px, pxm,  py, pym, uu, uv, vu, vv, &
!   &              hvbot, dzsig
! common /work2/ px, pxm,  py, pym

  real(8), parameter :: cx1 = 23.d0/12.d0, cx2 = -16.d0/12.d0, &
    &                   cx3 = 5.d0/12.d0
  real(8), save ::    xx1(nxydim),    yy1(nxydim)
  real(8), save ::    xx2(nxydim),    yy2(nxydim)
  real(8), save ::    xx3(nxydim),    yy3(nxydim)
  integer, save ::  ncall = 0
  logical, save ::   oeof

  if (oinit) then
     do ij=1,nxydim
        xx2(ij)=0.d0
        xx3(ij)=0.d0
        yy2(ij)=0.d0
        yy3(ij)=0.d0
     end do

#ifdef OPT_TRIPOLE
     call rstadd(xx2, oeof, nxdim, nydim, 1, 'XX2', 'SFC', &
        &                                 -1.d0,  -1,  -1)
     call rstadd(xx3, oeof, nxdim, nydim, 1, 'XX3', 'SFC', &
        &                                 -1.d0,  -1,  -1)
     call rstadd(yy2, oeof, nxdim, nydim, 1, 'YY2', 'SFC', &
        &                                 -1.d0,  -1,  -1)
     call rstadd(yy3, oeof, nxdim, nydim, 1, 'YY3', 'SFC', &
        &                                 -1.d0,  -1,  -1)
#else
     call rstadd(xx2, oeof, nxdim, nydim, 1, 'XX2', 'SFC')
     call rstadd(xx3, oeof, nxdim, nydim, 1, 'XX3', 'SFC')
     call rstadd(yy2, oeof, nxdim, nydim, 1, 'YY2', 'SFC')
     call rstadd(yy3, oeof, nxdim, nydim, 1, 'YY3', 'SFC')
#endif
     return
  end if

  if (ofinal) then
     call finadd(xx2, nxdim, nydim, 1, 'XX2', 'SFC')
     call finadd(xx3, nxdim, nydim, 1, 'XX3', 'SFC')
     call finadd(yy2, nxdim, nydim, 1, 'YY2', 'SFC')
     call finadd(yy3, nxdim, nydim, 1, 'YY3', 'SFC')
     return
  end if


  if (ofirst) then
     ofirst = .false.
     cof = -0.25d0 * gravit / rhoo * 1.d-3
  end if

#ifdef ADF_
!$acc data copyin(dy, ds) 
!$acc data copyin(cor) 
!$acc data copyin(ryu, rym, rxu) 
!$acc data copyin(hyu, hxu)  
!$acc data copyin(amskv) 
!$acc data copyin(dzv) 
!$acc data copyin(hyxu, hxyu) 
!---  arguments
!$acc data copy(gx, gy, xx, yy) 
!$acc data copy(gxx, gyy) 
!$acc data copyin(uy, vy) 
!$acc data copyin(r, hy) 
!$acc data copyin(ux, vx) 
!$acc data copyin(uadv, vadv) 
!$acc data create(hvbot, dzsig)  
!$acc data create(uu, vu, uv, vv)   
!$acc data copy(pxm, pym, px, py)   
#endif

#ifdef OMP_
!$omp parallel do private(ij)
#elif  ACC_
!$acc data copyin(dy, ds, rym) 
!$acc data copyin(hy) 
!$acc data copy(hvbot, dzsig)  
!$acc kernels 
#endif
  do ij = ijvstr, ijvend+nxdim
     hvbot(ij) = (  (hy(ij) + hy(ij+le)) * dy(ij) &
        &         + (hy(ij+ln) + hy(ij+lne)) * dy(ij+ln)) * &
        &        rym(ij) * 0.25d0 &
        &      + zbot
  end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels 
#endif

#ifdef OMP_
!$omp parallel do private(ij, k)
#elif  ACC_
!$acc kernels 
#endif
  do k = kstr, kstr+kz-1
     do ij = ijvstr, ijvend+nxdim
        dzsig(ij, k) = ds(k) * hvbot(ij)
     end do
  end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copyin(dy, ds, rym)
!$acc end data ! copyin(hy)
!$acc end data ! copy(hvbot, dzsig)
#endif

#ifdef OMP_
!$omp parallel do private(ij)
#elif  ACC_
!$acc data copy(gx, gy) 
!$acc data copyin(vx, ux, r) 
!$acc data copy(gxx, gyy) 
!$acc data copy(xx, yy) 
!$acc data copyin(ryu, rxu, rym, cor, amskv) 
!$acc data copy(py, pxm, pym, px)   
!--- local varibles
!$acc data copyin(dzsig)   
!$acc kernels 
#endif
  do ij = ijvstr, ijvend
     pxm(ij, kstr) = (  r(ij+lne, kstr) + r(ij+le, kstr) &
        &             - r(ij+ln , kstr) - r(ij   , kstr)) * &
        &            cof * rx * rxu(ij) * dzsig(ij, kstr) * &
        &            amskv(ij, kstr)
     px (ij, kstr) = pxm(ij, kstr)
     gx (ij, kstr) = gx(ij, kstr) &
        &          + (cor(ij) * vx(ij, kstr) + px(ij, kstr)) * &
        &            amskv(ij, kstr)
     xx (ij, kstr) = xx(ij, kstr) &
        &          + px(ij, kstr) * amskv(ij, kstr) 
     gxx(ij)       = xx(ij, kstr) * dzsig(ij, kstr) * &
        &            amskv(ij, kstr)
     pym(ij, kstr) = (  r(ij+lne, kstr) + r(ij+ln, kstr) &
        &             - r(ij+le , kstr) - r(ij   , kstr)) * &
        &            cof * rym(ij) * ryu(ij) * dzsig(ij, kstr) * &
        &            amskv(ij, kstr)
     py (ij, kstr) = pym(ij, kstr)
     gy (ij, kstr) = gy(ij, kstr) &
        &          - (cor(ij) * ux(ij, kstr) - py(ij, kstr)) * &
        &            amskv(ij, kstr)
     yy (ij, kstr) = yy(ij, kstr) &
        &          + py(ij, kstr) * amskv(ij, kstr)
     gyy(ij)       = yy(ij, kstr) * dzsig(ij, kstr) * &
        &            amskv(ij, kstr)
  end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copy(gx, gy)
!$acc end data ! copyin(vx, ux, r)
!$acc end data ! copy(gxx, gyy)
!$acc end data ! copy(xx, yy)
!$acc end data ! copyin(ryu, rxu, rym, cor, amskv)
!$acc end data ! copy(py, pxm, pym, px)
!--- local varibles
!$acc end data ! copyin(dzsig)
#endif

#ifdef OMP_
!$omp parallel
#elif  ACC_
!$acc data copy(xx, yy, gx, gy, gxx, gyy) 
!$acc data copyin(r, ux, vx) 
!$acc data copy(pxm, pym, px, py)   
!$acc data copyin(amskv, cor, rym, rxu, ryu) 
!--- local varibles
!$acc data copyin(dzsig)   
!$acc kernels 
!$acc loop seq
#endif
  do k = kstr+1, kstr+kz-1
#ifdef OMP_
!$omp do private(ij)
#elif  ACC_
!$acc loop independent
#endif
     do ij = ijvstr, ijvend
        pxm(ij, k) = (  r(ij+lne, k) + r(ij+le, k) &
           &          - r(ij+ln , k) - r(ij   , k)) * &
           &         cof * rx * rxu(ij) * dzsig(ij, k) * &
           &         amskv(ij, k)
        px (ij, k) = (px(ij, k-1) + pxm(ij, k-1) + pxm(ij, k)) * &
           &         amskv(ij, k)
        gx (ij, k) = gx(ij, k) &
           &       + (cor(ij) * vx(ij, k) + px(ij, k)) * &
           &         amskv(ij, k)
        xx (ij, k) = xx(ij, k) &
           &       + px(ij, k) * amskv(ij, k)
        gxx(ij)    = gxx(ij) &
           &       + xx(ij, k) * dzsig(ij, k) * amskv(ij, k)
        pym(ij, k) = (  r(ij+lne, k) + r(ij+ln, k) &
           &          - r(ij+le , k) - r(ij   , k)) * &
           &         cof * rym(ij) * ryu(ij) * dzsig(ij, k) * &
           &         amskv(ij, k)
        py (ij, k) = (py(ij, k-1) + pym(ij, k-1) + pym(ij, k)) * &
           &         amskv(ij, k)
        gy (ij, k) = gy(ij, k) &
           &       - (cor(ij) * ux(ij, k) - py(ij, k)) * &
           &         amskv(ij, k)
        yy (ij, k) = yy(ij, k) &
           &       + py(ij, k) * amskv(ij, k)
        gyy(ij)    = gyy(ij) + yy(ij, k) * dzsig(ij, k) * &
           &         amskv(ij, k)
     end do
  end do
#ifdef OMP_
!$omp end parallel
#elif  ACC_
!$acc end kernels
!$acc end data ! copy(xx, yy, gx, gy, gxx, gyy)
!$acc end data ! copyin(r, ux, vx)
!$acc end data ! copy(pxm, pym, px, py)
!$acc end data ! copyin(amskv, cor, rym, rxu, ryu)
!--- local varibles
!$acc end data ! copyin(dzsig)
#endif

#ifdef OMP_
!$omp parallel
#elif  ACC_
!$acc data copyin(dzv, amskv, cor, rym, rxu, ryu) 
!$acc data copy(py, pxm, px, pym)  
!$acc data copyin(ux, r, vx) 
!$acc data copy(yy, xx, gx, gyy, gxx, gy) 
!$acc kernels 
!$acc loop seq
#endif
  do k = kstr+kz, kend
#ifdef OMP_
!$omp do private(ij)
#elif  ACC_
!$acc loop independent
#endif
     do ij = ijvstr, ijvend
        pxm(ij, k) = (  r(ij+lne, k) - r(ij+ln, k) &
           &          + r(ij+le , k) - r(ij   , k)) * &
           &         cof * rx * rxu(ij) * dzv(ij, k) * &
           &         amskv(ij, k)
        px (ij, k) = (px(ij, k-1) + pxm(ij, k-1) + pxm(ij, k)) * &
           &         amskv(ij, k)
        gx (ij, k) = gx(ij, k) &
           &       + (cor(ij) * vx(ij, k) + px(ij, k)) * &
           &         amskv(ij, k)
        xx (ij, k) = xx(ij, k) &
           &       + px(ij, k) * amskv(ij, k)
        gxx(ij)    = gxx(ij) + xx(ij, k) * dzv(ij, k) * amskv(ij, k)
        pym(ij, k) = (  r(ij+lne, k) + r(ij+ln, k) &
           &          - r(ij+le , k) - r(ij   , k)) * &
           &         cof * rym(ij) * ryu(ij) * dzv(ij, k) * &
           &         amskv(ij, k)
        py (ij, k) = (py(ij, k-1) + pym(ij, k-1) + pym(ij, k)) * &
           &         amskv(ij, k)
        gy (ij, k) = gy(ij, k) &
           &       - (cor(ij) * ux(ij, k) - py(ij, k)) * &
           &         amskv(ij, k)
        yy (ij, k) = yy(ij, k) &
           &       + py(ij, k) * amskv(ij, k)
        gyy(ij)    = gyy(ij) + yy(ij, k) * dzv(ij, k) * amskv(ij, k)
     end do
  end do
#ifdef OMP_
!$omp end parallel
#elif  ACC_
!$acc end kernels
!$acc end data ! copyin(dzv, amskv, cor, rym, rxu, ryu)
!$acc end data ! copy(py, pxm, px, pym)
!$acc end data ! copyin(ux, r, vx)
!$acc end data ! copy(yy, xx, gx, gyy, gxx, gy)
#endif

#ifdef OMP_
!$omp parallel do private(ij)
#elif  ACC_
!$acc data copy(xx1, yy1) 
!$acc kernels 
#endif
  do ij=1,nxydim
     xx1(ij)=0.d0
     yy1(ij)=0.d0
  end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copy(xx1, yy1)
#endif

#ifdef OMP_
!$omp parallel
#elif  ACC_
!$acc data copyin(uadv, vadv, uy, vy) 
!$acc data copy(yy1, xx1) 
!--- local varibales
!$acc data copyin(dzsig)    
!$acc data copy(uu, vu, uv, vv)   
!$acc data copyin(ryu, hxu, rym, hyu, rxu) 
!$acc kernels 
!$acc loop seq
#endif
  do k = kstr, kstr+kz-1
#ifdef OMP_
!$omp do private(ij)
#endif
     do ij = ijvstr, ijvend+nxdim
        uu(ij) = uadv(ij, k) * 0.5d0 * &
           &    (  uy(ij, k) * hyu(ij) &
           &     + uy(ij+lw, k) * hyu(ij+lw))
        vu(ij) = vadv(ij, k) * 0.5d0 * &
           &    (  uy(ij, k) * hxu(ij) &
           &     + uy(ij+ls, k) * hxu(ij+ls))
        uv(ij) = uadv(ij, k) * 0.5d0 * &
           &    (  vy(ij, k) * hyu(ij) &
           &     + vy(ij+lw, k) * hyu(ij+lw))
        vv(ij) = vadv(ij, k) * 0.5d0 * &
           &    (  vy(ij, k) * hxu(ij) &
           &     + vy(ij+ls, k) * hxu(ij+ls))
     end do
#ifdef OMP_
!$omp do private(ij)
#endif
     do ij = ijvstr, ijvend
!x        gxx(ij) = gxx(ij)
        xx1(ij) = xx1(ij) &
           &    - (  (uu(ij+le) - uu(ij)) * rx &
           &       + (vu(ij+ln) - vu(ij)) * rym(ij)) * &
           &      rxu(ij) * ryu(ij) * dzsig(ij,k)
!x        gyy(ij) = gyy(ij)
        yy1(ij) = yy1(ij) &
           &    - (  (uv(ij+le) - uv(ij)) * rx &
           &       + (vv(ij+ln) - vv(ij)) * rym(ij)) * &
           &      rxu(ij) * ryu(ij) * dzsig(ij,k)
     end do
  end do
#ifdef OMP_
!$omp end parallel
#elif  ACC_
!$acc end kernels
!$acc end data ! copyin(uadv, vadv, uy, vy)
!$acc end data ! copy(yy1, xx1)
!--- local varibales
!$acc end data ! copyin(dzsig)
!$acc end data ! copy(uu, vu, uv, vv)
!$acc end data ! copyin(ryu, hxu, rym, hyu, rxu)
#endif
! -----------------------------------------------------------
!      do ij = 1, nxydim
!         uu(ij) = 0.d0
!         uv(ij) = 0.d0
!         vu(ij) = 0.d0
!         vv(ij) = 0.d0
!      end do
!
!      do k = kstr, kstr+kz-1
!         do ij = ijvstr, ijvend+nxdim
!            uu(ij) = uu(ij) &
!        &      + uadv(ij, k) * 0.5d0 * &
!        &        (  uy(ij, k) * hyu(ij) &
!        &         + uy(ij+lw, k) * hyu(ij+lw)) * ds(k)
!            vu(ij) = vu(ij) &
!        &      + vadv(ij, k) * 0.5d0 * &
!        &        (  uy(ij, k) * hxu(ij) &
!        &         + uy(ij+ls, k) * hxu(ij+ls)) * ds(k) &
!            uv(ij) = uv(ij) &
!        &      + uadv(ij, k) * 0.5d0 * &
!        &        (  vy(ij, k) * hyu(ij) &
!        &         + vy(ij+lw, k) * hyu(ij+lw)) * ds(k)
!            vv(ij) = vv(ij) &
!        &      + vadv(ij, k) * 0.5d0 * &
!        &        (  vy(ij, k) * hxu(ij) &
!        &         + vy(ij+ls, k) * hxu(ij+ls)) * ds(k)
!         end do
!      end do
!
!      do ij = ijvstr, ijvend
!         gxx(ij) = gxx(ij) &
!        &    - (  (uu(ij+le) - uu(ij)) * rx &
!        &       + (vu(ij+ln) - vu(ij)) * rym(ij)) * &
!        &      rxu(ij) * ryu(ij) * hvbot(ij)
!         gyy(ij) = gyy(ij) &
!        &    - (  (uv(ij+le) - uv(ij)) * rx &
!        &       + (vv(ij+ln) - vv(ij)) * rym(ij)) * &
!        &      rxu(ij) * ryu(ij) * hvbot(ij)
!      enddo
! ------------------------------------------------------------
#ifdef OMP_
!$omp parallel do private(ij)
#elif  ACC_
!--- local
!$acc data copy(uu, vu, uv, vv)  
!$acc kernels 
#endif
  do ij = 1, nxydim
     uu(ij) = 0.d0
     uv(ij) = 0.d0
     vu(ij) = 0.d0
     vv(ij) = 0.d0
  end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!--- local
!$acc end data ! copy(uu, vu, uv, vv)
#endif
!
#ifdef OMP_
!$omp parallel
#elif  ACC_
!$acc data copyin(uadv, vadv, uy, vy)  
!$acc data copyin(hyu, hxu)  
!--- local
!$acc data copy(uu, vu, uv, vv)   
!$acc kernels 
!$acc loop seq
#endif
  do k = kstr+kz, kend
#ifdef OMP_
!$omp do private(ij)
#endif
     do ij = ijvstr, ijvend+nxdim
        uu(ij) = uu(ij) &
           &   + uadv(ij, k) * 0.5d0 * &
           &     (  uy(ij, k) * hyu(ij) &
           &      + uy(ij+lw, k) * hyu(ij+lw))
        vu(ij) = vu(ij) &
           &   + vadv(ij, k) * 0.5d0 * &
           &     (  uy(ij, k) * hxu(ij) &
           &      + uy(ij+ls, k) * hxu(ij+ls))
        uv(ij) = uv(ij) &
           &   + uadv(ij, k) * 0.5d0 * &
           &     (  vy(ij, k) * hyu(ij) &
           &      + vy(ij+lw, k) * hyu(ij+lw))
        vv(ij) = vv(ij) &
           &   + vadv(ij, k) * 0.5d0 * &
           &     (  vy(ij, k) * hxu(ij) &
           &      + vy(ij+ls, k) * hxu(ij+ls))
     end do
  end do
#ifdef OMP_
!$omp end parallel
#elif  ACC_
!$acc end kernels
!--- local
!$acc end data ! copy(uu, vu, uv, vv)
!$acc end data ! copyin(uadv, vadv, uy, vy) 
!$acc end data ! copyin(hyu, hxu) 
#endif

#ifdef OMP_
!$omp parallel do private(ij)
#elif  ACC_
!$acc data copyin(uu, vu, uv, vv)   
!$acc data copy(xx1, yy1) 
!$acc data copyin(rym, rxu, ryu) 
!$acc kernels 
#endif
  do ij = ijvstr, ijvend
!x     gxx(ij) = gxx(ij)
     xx1(ij) = xx1(ij) &
        &    - (  (uu(ij+le) - uu(ij)) * rx &
        &       + (vu(ij+ln) - vu(ij)) * rym(ij)) * &
        &      rxu(ij) * ryu(ij)
!x     gyy(ij) = gyy(ij)
     yy1(ij) = yy1(ij) &
        &    - (  (uv(ij+le) - uv(ij)) * rx &
        &       + (vv(ij+ln) - vv(ij)) * rym(ij)) * &
        &      rxu(ij) * ryu(ij)
  end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copyin(uu, vu, uv, vv)
!$acc end data ! copy(xx1, yy1)
!$acc end data ! copyin(rym, rxu, ryu)
#endif

#ifdef OMP_
!$omp parallel do private(ij)
#elif  ACC_
!$acc data copyin(uy, vy) 
!$acc data copyin(amskv) 
!$acc data copyin(dzsig)   
!$acc data copy(uu, uv, vv)   
!$acc kernels 
#endif
  do ij = ijvstr, ijvend
     uu(ij) = uy(ij, kstr) * uy(ij, kstr) * &
        &     dzsig(ij, kstr) * amskv(ij, kstr)
     uv(ij) = uy(ij, kstr) * vy(ij, kstr) * &
        &     dzsig(ij, kstr) * amskv(ij, kstr)
     vv(ij) = vy(ij, kstr) * vy(ij, kstr) * &
        &     dzsig(ij, kstr) * amskv(ij, kstr)
  end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copyin(uy, vy)
!$acc end data ! copyin(amskv)
!$acc end data ! copyin(dzsig)
!$acc end data ! copy(uu, uv, vv)
#endif

#ifdef OMP_
!$omp parallel
#elif  ACC_
!$acc data copyin(uy, vy) 
!$acc data copyin(amskv) 
!$acc data copyin(dzsig)  
!$acc data copy(uu, uv, vv)  
!$acc kernels 
!$acc loop seq
#endif
  do k = kstr+1, kstr+kz-1
#ifdef OMP_
!$omp do private(ij)
#endif
     do ij = ijvstr, ijvend
        uu(ij) = uu(ij) &
           &   + uy(ij, k) * uy(ij, k) * &
           &     dzsig(ij, k) * amskv(ij, k)
        uv(ij) = uv(ij) &
           &   + uy(ij, k) * vy(ij, k) * &
           &     dzsig(ij, k) * amskv(ij, k)
        vv(ij) = vv(ij) &
           &   + vy(ij, k) * vy(ij, k) * &
           &     dzsig(ij, k) * amskv(ij, k)
     end do
  end do
#ifdef OMP_
!$omp end parallel
#elif  ACC_
!$acc end kernels
!$acc end data ! copyin(uy, vy)
!$acc end data ! copyin(amskv)
!$acc end data ! copyin(dzsig)
!$acc end data ! copy(uu, uv, vv)
#endif

#ifdef OMP_
!$omp parallel
#elif  ACC_
!$acc data copy(uu, uv, vv)   
!$acc data copyin(uy, vy) 
!$acc data copyin(amskv, dzv) 
!$acc kernels 
!$acc loop seq
#endif
  do k = kstr+kz, kend
#ifdef OMP_
!$omp do private(ij)
#endif
     do ij = ijvstr, ijvend
        uu(ij) = uu(ij) &
           &   + uy(ij, k) * uy(ij, k) * dzv(ij, k) * amskv(ij, k)
        uv(ij) = uv(ij) &
           &   + uy(ij, k) * vy(ij, k) * dzv(ij, k) * amskv(ij, k)
        vv(ij) = vv(ij) &
           &   + vy(ij, k) * vy(ij, k) * dzv(ij, k) * amskv(ij, k)
     end do
  end do
#ifdef OMP_
!$omp end parallel
#elif  ACC_
!$acc end kernels
!$acc end data ! copyin(uy, vy)
!$acc end data ! copyin(amskv, dzv)
!$acc end data ! copy(uu, uv, vv)
#endif

#ifdef OMP_
!$omp parallel do private(ij)
#elif  ACC_
!$acc data copyin(vv, uv, uu)   
!$acc data copy(xx1, yy1) 
!$acc data copyin(hyxu, hxyu, amskv) 
!$acc kernels 
#endif
  do ij = ijvstr, ijvend
!x     gxx(ij) = gxx(ij) + vv(ij) * hyxu(ij) - uv(ij) * hxyu(ij)
!x     gyy(ij) = gyy(ij) + uu(ij) * hxyu(ij) - uv(ij) * hyxu(ij)
     xx1(ij) = xx1(ij) + vv(ij) * hyxu(ij) - uv(ij) * hxyu(ij)
     yy1(ij) = yy1(ij) + uu(ij) * hxyu(ij) - uv(ij) * hyxu(ij)
     xx1(ij) = xx1(ij) * amskv(ij,kstr)
     yy1(ij) = yy1(ij) * amskv(ij,kstr)
  end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copyin(vv, uv, uu)
!$acc end data ! copy(xx1, yy1)
!$acc end data ! copyin(hyxu, hxyu, amskv)
#endif


!cccccccc adams-bashforth scheme
  ncall=ncall +1

  if( (ncall .ge. 3) .or. (.not. oeof)) then
     ncall=3
#ifdef OMP_
!$omp parallel do private(ij)
#elif  ACC_
!$acc data copy(gxx, gyy) 
!$acc data copyin(xx3, yy1, yy3, yy2, xx2, xx1) 
!$acc kernels 
#endif
     do ij = ijvstr, ijvend
        gxx(ij) =  gxx(ij) &
           &    +  cx1 *xx1(ij)  + cx2* xx2(ij) + cx3*xx3(ij) 

        gyy(ij) =  gyy(ij) &
           &    +  cx1 *yy1(ij)  + cx2* yy2(ij) + cx3*yy3(ij) 
     end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copyin(xx3, yy1, yy3, yy2, xx2, xx1)
!$acc end data ! copy(gxx, gyy)
#endif
  else
     if(ncall .eq. 1) then ! forward
#ifdef OMP_
!$omp parallel do private(ij)
#elif  ACC_
!$acc data copy(gxx, gyy) 
!$acc data copyin(xx1, yy1) 
!$acc kernels 
#endif
        do ij = ijvstr, ijvend
           gxx(ij) =  gxx(ij) + xx1(ij)
           gyy(ij) =  gyy(ij) + yy1(ij)
        end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copyin(xx1, yy1)
!$acc end data ! copy(gxx, gyy)
#endif
     else if(ncall .eq. 2) then !2nd order ab       
#ifdef OMP_
!$omp parallel do private(ij)
#elif  ACC_
!$acc data copy(gxx, gyy) 
!$acc data copyin(yy2, yy1, xx1, xx2) 
!$acc kernels 
#endif
        do ij = ijvstr, ijvend
           gxx(ij) =  gxx(ij) &
              &    +  1.5d0 *xx1(ij)  -0.5d0* xx2(ij)

           gyy(ij) =  gyy(ij) &
              &    +  1.5d0 *yy1(ij)  -0.5d0* yy2(ij)
        end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copyin(yy2, yy1, xx1, xx2)
!$acc end data ! copy(gxx, gyy)
#endif
     end if
  end if


#ifdef OMP_
!$omp parallel do private(ij)
#elif  ACC_
!$acc data copyin(yy1, xx1) 
!$acc data copy(xx2, yy2) 
!$acc data copy(xx3, yy3) 
!$acc kernels 
#endif
  do ij = ijvstr, ijvend
     xx3(ij)=xx2(ij) !n-1 > n-2
     yy3(ij)=yy2(ij)

     xx2(ij)=xx1(ij) !n> n-1
     yy2(ij)=yy1(ij)
  end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copyin(yy1, xx1)
!$acc end data ! copy(xx2, yy2)
!$acc end data ! copy(xx3, yy3)
#endif

#ifdef ADF_
!$acc end data ! copyin(dy, ds)
!$acc end data ! copyin(cor)
!$acc end data ! copyin(ryu, rym, rxu)
!$acc end data ! copyin(hyu, hxu) 
!$acc end data ! copyin(amskv)
!$acc end data ! copyin(dzv)
!$acc end data 
!---  arguments
!$acc end data ! copy(gx, gy, xx, yy)
!$acc end data ! copy(gxx, gyy)
!$acc end data ! copyin(uy, vy)
!$acc end data ! copyin(r, hy)
!$acc end data ! copyin(ux, vx)
!$acc end data ! copyin(uadv, vadv)
!$acc end data ! copy(hvbot, dzsig)  
!$acc end data ! copy(pxm, pym, px, py)   
!$acc end data ! copy(uu, vu, uv, vv)   
#endif

  return
end subroutine tnduvd

#ifdef OPT_BBL
! *********************************************************************

subroutine tnduvb( &
  &    gxx,    gyy,     gx,     gy, &
  &     xx,     yy,     uy,     vy,     ty, &
  &     ux,     vx )

! --- information -----------------------------------------------------
!
!  pressure gradient and coriolis terms for the bbl momentum eqs. see
! nakano and suginohara (2002, jpo).
!
! ---------------------------------------------------------------------

  use xprst

  real(8), intent(inout) ::     gxx(nxydim)       ,    gyy(nxydim)
  real(8), intent(inout) ::      gx(nxydim, nzdim),     gy(nxydim, nzdim)
  real(8), intent(inout) ::      xx(nxydim, nzdim),     yy(nxydim, nzdim)
  real(8), intent(in)    ::      ux(nxydim, nzdim),     vx(nxydim, nzdim)
  real(8), intent(in)    ::      uy(nxydim, nzdim),     vy(nxydim, nzdim)
  real(8), intent(in)    ::      ty(nxydim, nzdim, ntdim)     

  real(8), save ::  c0b(nxydim), c1b(nxydim), c2b(nxydim), c3b(nxydim)
  real(8), save ::  c4b(nxydim), c5b(nxydim), c6b(nxydim)
  real(8), save ::  d0b(nxydim), d1b(nxydim), d2b(nxydim), d3b(nxydim)
  real(8), save ::  d4b(nxydim), d5b(nxydim), d6b(nxydim), d7b(nxydim)
  real(8), save ::  d8b(nxydim), d9b(nxydim)

  real(8), save ::  dept0(nxydim)
  real(8), save ::    cof
  logical, save :: ofirst = .true.

  real(8) ::      tl,     sl,   rtmp
  real(8) ::       p,     pe,     pn,    pne
  integer ::      ij,      k,     kk
  integer ::    ijle,   ijln,   ijne

  real(8), parameter :: cx1 = 23.d0/12.d0, cx2 = -16.d0/12.d0, &
    &                   cx3 = 5.d0/12.d0
  real(8), save ::    xx1(nxydim),    yy1(nxydim)
  real(8), save ::    xx2(nxydim),    yy2(nxydim)
  real(8), save ::    xx3(nxydim),    yy3(nxydim)
  integer, save ::  ncall = 0
  logical, save ::   oeof

!===== Define statement function 
  real(8) ::    rbbl,     tb,     sb
  real(8) ::  c0, c1, c2, c3, c4, c5, c6
  real(8) ::  d0, d1, d2, d3, d4, d5, d6, d7, d8, d9

  rbbl(tb, sb, c0, c1, c2, c3, c4, c5, c6, &
     &     d0, d1, d2, d3, d4, d5, d6, d7, d8, d9) = &
     &    (c0 + (c1 + (c2 + c3 * tb) * tb) * tb &
     &        + (c4 + c5 * tb + c6 * sb) * sb) &
     &  / (d0 + (d1 + (d2 + (d3 + d4 * tb) * tb) * tb) * tb &
     &        + (d5 + (d6 + d7 * tb * tb) * tb &
     &              + (d8 + d9 * tb * tb) * sqrt(sb)) * sb) &
     &  - 1.d3
!===== 

  if (oinit) then
     do ij=1,nxydim
        xx1(ij)=0.d0
        xx2(ij)=0.d0
        xx3(ij)=0.d0
        yy1(ij)=0.d0
        yy2(ij)=0.d0
        yy3(ij)=0.d0
     end do
#ifdef OPT_TRIPOLE
     call rstadd(xx2, oeof, nxdim, nydim, 1, 'XX2', 'SFC', &
        &                                 -1.d0,  -1,  -1)
     call rstadd(xx3, oeof, nxdim, nydim, 1, 'XX3', 'SFC', &
        &                                 -1.d0,  -1,  -1)
     call rstadd(yy2, oeof, nxdim, nydim, 1, 'YY2', 'SFC', &
        &                                 -1.d0,  -1,  -1)
     call rstadd(yy3, oeof, nxdim, nydim, 1, 'YY3', 'SFC', &
        &                                 -1.d0,  -1,  -1)
#else
     call rstadd(xx2, oeof, nxdim, nydim, 1, 'XX2', 'SFC')
     call rstadd(xx3, oeof, nxdim, nydim, 1, 'XX3', 'SFC')
     call rstadd(yy2, oeof, nxdim, nydim, 1, 'YY2', 'SFC')
     call rstadd(yy3, oeof, nxdim, nydim, 1, 'YY3', 'SFC')
#endif
     return
  end if

  if (ofinal) then
     call finadd(xx2, nxdim, nydim, 1, 'XX2', 'SFC')
     call finadd(xx3, nxdim, nydim, 1, 'XX3', 'SFC')
     call finadd(yy2, nxdim, nydim, 1, 'YY2', 'SFC')
     call finadd(yy3, nxdim, nydim, 1, 'YY3', 'SFC')
     return
  end if


  if (ofirst) then
     ofirst = .false.
     call secfbb( &
        &    c0b,    c1b,    c2b,    c3b, &
        &    c4b,    c5b,    c6b, &
        &    d0b,    d1b,    d2b,    d3b,    d4b, &
        &    d5b,    d6b,    d7b,    d8b,    d9b )

     cof = -0.5d0 * gravit / rhoo * 1.d-3
     do ij = 1, nxydim
        dept0(ij) = 0.d0
        do k = kstr, nbot(ij)-1
           dept0(ij) = dept0(ij) + dz0(k)
        end do
     end do
  end if

#ifdef ADF_
!$acc data copyin(rxu, rym, nbotv) 
!$acc data copyin(nbot, ryu, cor, amskvb, dzv, dz0) 
!$acc data copyin(amskvb, hxyu, dzv, hyxu) 
!$acc data copy(gxx, gyy) 
!$acc data copy(gx, gy) 
!$acc data copy(xx, yy) 
!$acc data copyin(ux, vx) 
!$acc data copyin(uy, vy) 
!$acc data copyin(ty) 
!$acc data copyin(c0b, c1b, c2b, c3b) 
!$acc data copyin(c4b, c5b, c6b) 
!$acc data copyin(d0b, d1b, d2b, d3b) 
!$acc data copyin(d4b, d5b, d6b, d7b) 
!$acc data copyin(d8b, d9b) 
!$acc data copyin(dept0) 
#endif

#ifdef OMP_
!$omp parallel do private(ij, k, ijle, ijln, ijne, tl, sl, rtmp, p, kk, pe, pn, pne)
#elif  ACC_
!$acc data copyin(vx, ux, ty) 
!$acc data copy(yy, xx, gyy, gx, gxx, gy) 
!$acc data copyin(px, py)   
!$acc data copyin(rxu, rym, nbotv) 
!$acc data copyin(nbot, ryu, cor, amskvb, dzv, dz0) 
!$acc kernels 
#endif
  do ij = ijvstr, ijvend

     k = max(nbotv(ij)-1, 1)
     ijle = ij + le
     ijln = ij + ln
     ijne = ij + lne

     tl = ty(ij, k, 1)
     sl = ty(ij, k, 2)
     rtmp = rbbl(tl, sl, &
        &        c0b(ij), c1b(ij), c2b(ij), c3b(ij), &
        &        c4b(ij), c5b(ij), c6b(ij), &
        &        d0b(ij), d1b(ij), d2b(ij), d3b(ij), d4b(ij), &
        &        d5b(ij), d6b(ij), d7b(ij), d8b(ij), d9b(ij))
     p = rtmp * dz0(k) * 0.5d0
     do kk = k+1, nbot(ij)-1
        tl = ty(ij, kk, 1)
        sl = ty(ij, kk, 2)
        rtmp = rbbl(tl, sl, &
           &        c0b(ij), c1b(ij), c2b(ij), c3b(ij), &
           &        c4b(ij), c5b(ij), c6b(ij), &
           &        d0b(ij), d1b(ij), d2b(ij), d3b(ij), d4b(ij), &
           &        d5b(ij), d6b(ij), d7b(ij), d8b(ij), d9b(ij))
        p = p + rtmp * dz0(kk)
     end do
     tl = ty(ij, kend, 1)
     sl = ty(ij, kend, 2)
     rtmp = rbbl(tl, sl, &
        &        c0b(ij), c1b(ij), c2b(ij), c3b(ij), &
        &        c4b(ij), c5b(ij), c6b(ij), &
        &        d0b(ij), d1b(ij), d2b(ij), d3b(ij), d4b(ij), &
        &        d5b(ij), d6b(ij), d7b(ij), d8b(ij), d9b(ij))
     p = p + rtmp * dz0(kend) * 0.5d0

     tl = ty(ijle, k, 1)
     sl = ty(ijle, k, 2)
     rtmp = rbbl(tl, sl, &
        &        c0b(ij), c1b(ij), c2b(ij), c3b(ij), &
        &        c4b(ij), c5b(ij), c6b(ij), &
        &        d0b(ij), d1b(ij), d2b(ij), d3b(ij), d4b(ij), &
        &        d5b(ij), d6b(ij), d7b(ij), d8b(ij), d9b(ij))
     pe = rtmp * dz0(k) * 0.5d0
     do kk = k+1, nbot(ijle)-1
        tl = ty(ijle, kk, 1)
        sl = ty(ijle, kk, 2)
        rtmp = rbbl(tl, sl, &
           &        c0b(ij), c1b(ij), c2b(ij), c3b(ij), &
           &        c4b(ij), c5b(ij), c6b(ij), &
           &        d0b(ij), d1b(ij), d2b(ij), d3b(ij), d4b(ij), &
           &        d5b(ij), d6b(ij), d7b(ij), d8b(ij), d9b(ij))
        pe = pe + rtmp * dz0(kk)
     end do
     tl = ty(ijle, kend, 1)
     sl = ty(ijle, kend, 2)
     rtmp = rbbl(tl, sl, &
        &        c0b(ij), c1b(ij), c2b(ij), c3b(ij), &
        &        c4b(ij), c5b(ij), c6b(ij), &
        &        d0b(ij), d1b(ij), d2b(ij), d3b(ij), d4b(ij), &
        &        d5b(ij), d6b(ij), d7b(ij), d8b(ij), d9b(ij))
     pe = pe + rtmp * dz0(kend) * 0.5d0
     
     tl = ty(ijln, k, 1)
     sl = ty(ijln, k, 2)
     rtmp = rbbl(tl, sl, &
        &        c0b(ij), c1b(ij), c2b(ij), c3b(ij), &
        &        c4b(ij), c5b(ij), c6b(ij), &
        &        d0b(ij), d1b(ij), d2b(ij), d3b(ij), d4b(ij), &
        &        d5b(ij), d6b(ij), d7b(ij), d8b(ij), d9b(ij))
     pn = rtmp * dz0(k) * 0.5d0
     do kk = k+1, nbot(ijln)-1
        tl = ty(ijln, kk, 1)
        sl = ty(ijln, kk, 2)
        rtmp = rbbl(tl, sl, &
           &        c0b(ij), c1b(ij), c2b(ij), c3b(ij), &
           &        c4b(ij), c5b(ij), c6b(ij), &
           &        d0b(ij), d1b(ij), d2b(ij), d3b(ij), d4b(ij), &
           &        d5b(ij), d6b(ij), d7b(ij), d8b(ij), d9b(ij))
        pn = pn + rtmp * dz0(kk)
     end do
     tl = ty(ijln, kend, 1)
     sl = ty(ijln, kend, 2)
     rtmp = rbbl(tl, sl, &
        &        c0b(ij), c1b(ij), c2b(ij), c3b(ij), &
        &        c4b(ij), c5b(ij), c6b(ij), &
        &        d0b(ij), d1b(ij), d2b(ij), d3b(ij), d4b(ij), &
        &        d5b(ij), d6b(ij), d7b(ij), d8b(ij), d9b(ij))
     pn = pn + rtmp * dz0(kend) * 0.5d0
     
     tl = ty(ijne, k, 1)
     sl = ty(ijne, k, 2)
     rtmp = rbbl(tl, sl, &
        &        c0b(ij), c1b(ij), c2b(ij), c3b(ij), &
        &        c4b(ij), c5b(ij), c6b(ij), &
        &        d0b(ij), d1b(ij), d2b(ij), d3b(ij), d4b(ij), &
        &        d5b(ij), d6b(ij), d7b(ij), d8b(ij), d9b(ij))
     pne = rtmp * dz0(k) * 0.5d0
     do kk = k+1, nbot(ijne)-1
        tl = ty(ijne, kk, 1)
        sl = ty(ijne, kk, 2)
        rtmp = rbbl(tl, sl, &
           &        c0b(ij), c1b(ij), c2b(ij), c3b(ij), &
           &        c4b(ij), c5b(ij), c6b(ij), &
           &        d0b(ij), d1b(ij), d2b(ij), d3b(ij), d4b(ij), &
           &        d5b(ij), d6b(ij), d7b(ij), d8b(ij), d9b(ij))
        pne = pne + rtmp * dz0(kk)
     end do
     tl = ty(ijne, kend, 1)
     sl = ty(ijne, kend, 2)
     rtmp = rbbl(tl, sl, &
        &        c0b(ij), c1b(ij), c2b(ij), c3b(ij), &
        &        c4b(ij), c5b(ij), c6b(ij), &
        &        d0b(ij), d1b(ij), d2b(ij), d3b(ij), d4b(ij), &
        &        d5b(ij), d6b(ij), d7b(ij), d8b(ij), d9b(ij))
     pne = pne + rtmp * dz0(kend) * 0.5d0

     tl = 0.25d0 * (  ty(ij  , kend, 1) + ty(ijle, kend, 1) &
        &           + ty(ijln, kend, 1) + ty(ijne, kend, 1))
     sl = 0.25d0 * (  ty(ij  , kend, 2) + ty(ijle, kend, 2) &
        &           + ty(ijln, kend, 2) + ty(ijne, kend, 2))
     rtmp = rbbl(tl, sl, &
        &        c0b(ij), c1b(ij), c2b(ij), c3b(ij), &
        &        c4b(ij), c5b(ij), c6b(ij), &
        &        d0b(ij), d1b(ij), d2b(ij), d3b(ij), d4b(ij), &
        &        d5b(ij), d6b(ij), d7b(ij), d8b(ij), d9b(ij))

     px(ij, kend) = px(ij, k) &
        &         + (  pne + pe - pn - p &
        &            - rtmp * (  dept0(ijne) + dept0(ijle) &
        &                      - dept0(ijln) - dept0(ij  ))) * &
        &           cof * rx * rxu(ij) * amskvb(ij)
     gx(ij, kend) = gx(ij, kend) &
        &         + (cor(ij) * vx(ij, kend) + px(ij, kend)) * &
        &           amskvb(ij)
     xx(ij, kend) = xx(ij, kend) + px(ij, kend) * amskvb(ij)
     gxx(ij)      = gxx(ij) &
        &         + xx(ij, kend) * dzv(ij, kend) * amskvb(ij)

     py(ij, kend) = py(ij, k) &
        &         + (  pne + pn - pe - p &
        &            - rtmp * (  dept0(ijne) + dept0(ijln) &
        &                      - dept0(ijle) - dept0(ij  ))) * &
        &           cof * rym(ij) * ryu(ij) * amskvb(ij)
     gy(ij, kend) = gy(ij, kend) &
        &         - (cor(ij) * ux(ij, kend) - py(ij, kend)) * &
        &           amskvb(ij)
     yy(ij, kend) = yy(ij, kend) + py(ij, kend) * amskvb(ij)
     gyy(ij)      = gyy(ij) &
        &         + yy(ij, kend) * dzv(ij, kend) * amskvb(ij)

  end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copyin(vx, ux, ty)
!$acc end data ! copy(yy, xx, gyy, gx, gxx, gy)
!$acc end data ! copyin(px, py)
!$acc end data ! copyin(rxu, rym, nbotv)
!$acc end data ! copyin(nbot, ryu, cor, amskvb, dzv, dz0)
#endif

#ifdef OMP_
!$omp parallel do private(ij)
#elif  ACC_
!$acc data copy(yy1, xx1) 
!$acc data copyin(amskvb, hxyu, dzv, hyxu) 
!$acc data copyin(uy, vy) 
!$acc kernels 
#endif
  do ij = ijvstr, ijvend
!x     gxx(ij) = gxx(ij) +
     xx1(ij) = &
        &    (  vy(ij, kend) * vy(ij, kend) * hyxu(ij) &
        &     - uy(ij, kend) * vy(ij, kend) * hxyu(ij)) * &
        &    dzv(ij, kend) * amskvb(ij)
!x     gyy(ij) = gyy(ij) +
     yy1(ij) = &
        &    (  uy(ij, kend) * uy(ij, kend) * hxyu(ij) &
        &     - uy(ij, kend) * vy(ij, kend) * hyxu(ij)) * &
        &    dzv(ij, kend) * amskvb(ij)
  end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copy(yy1, xx1)
!$acc end data ! copyin(amskvb, hxyu, dzv, hyxu)
!$acc end data ! copyin(uy, vy)
#endif

!cccccccc adams-bashforth scheme
  ncall=ncall +1

  if( (ncall .ge. 3) .or. (.not. oeof)) then
     ncall=3
#ifdef OMP_
!$omp parallel do private(ij)
#elif  ACC_
!$acc data copy(gyy, gxx) 
!$acc data copyin(xx3, yy1, yy3, yy2, xx2, xx1) 
!$acc kernels 
#endif
     do ij = ijvstr, ijvend
        gxx(ij) =  gxx(ij) &
           &    +  cx1 *xx1(ij)  + cx2* xx2(ij) + cx3*xx3(ij) 

        gyy(ij) =  gyy(ij) &
           &    +  cx1 *yy1(ij)  + cx2* yy2(ij) + cx3*yy3(ij) 
     end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copyin(xx3, yy1, yy3, yy2, xx2, xx1)
!$acc end data ! copy(gyy, gxx)
#endif
  else
     if(ncall .eq. 1) then ! forward
#ifdef OMP_
!$omp parallel do private(ij)
#elif  ACC_
!$acc data copyin(yy1, xx1) 
!$acc data copy(gxx, gyy) 
!$acc kernels 
#endif
        do ij = ijvstr, ijvend
           gxx(ij) =  gxx(ij) + xx1(ij)
           gyy(ij) =  gyy(ij) + yy1(ij)
        end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copyin(yy1, xx1)
!$acc end data ! copy(gxx, gyy)
#endif
     else if(ncall .eq. 2) then !2nd order ab       
#ifdef OMP_
!$omp parallel do private(ij)
#elif  ACC_
!$acc data copyin(yy2, yy1, xx1, xx2) 
!$acc data copy(gxx, gyy) 
!$acc kernels 
#endif
        do ij = ijvstr, ijvend
           gxx(ij) =  gxx(ij) &
              &    +  1.5d0 *xx1(ij)  -0.5d0* xx2(ij)

           gyy(ij) =  gyy(ij) &
              &    +  1.5d0 *yy1(ij)  -0.5d0* yy2(ij)
        end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copyin(yy2, yy1, xx1, xx2)
!$acc end data ! copy(gxx, gyy)
#endif
     end if
  end if


#ifdef OMP_
!$omp parallel do private(ij)
#elif  ACC_
!$acc data copy(xx3, yy3) 
!$acc data copy(yy2, xx2) 
!$acc data copyin(yy1, xx1) 
!$acc kernels 
#endif
  do ij = ijvstr, ijvend
     xx3(ij)=xx2(ij) !n-1 > n-2
     yy3(ij)=yy2(ij)

     xx2(ij)=xx1(ij) !n> n-1
     yy2(ij)=yy1(ij)
  end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copy(xx3, yy3)
!$acc end data ! copy(yy2, xx2)
!$acc end data ! copyin(yy1, xx1)
#endif

#ifdef ADF_
!$acc end data ! copyin(rxu, rym, nbotv)
!$acc end data ! copyin(nbot, ryu, cor, amskvb, dzv, dz0)
!$acc end data ! copyin(amskvb, hxyu, dzv, hyxu)
!$acc end data ! copy(gxx, gyy)
!$acc end data ! copy(gx, gy)
!$acc end data ! copy(xx, yy)
!$acc end data ! copyin(ux, vx)
!$acc end data ! copyin(uy, vy)
!$acc end data ! copyin(ty)
!$acc end data ! copyin(c0b, c1b, c2b, c3b)
!$acc end data ! copyin(c4b, c5b, c6b)
!$acc end data ! copyin(d0b, d1b, d2b, d3b)
!$acc end data ! copyin(d4b, d5b, d6b, d7b)
!$acc end data ! copyin(d8b, d9b)
!$acc end data ! copyin(dept0)
#endif

  return

end subroutine tnduvb
! *********************************************************************
! do not use internal function now, instead we use statement function.
!function rbbl(tb, sb, c0, c1, c2, c3, c4, c5, c6, &
!  &           d0, d1, d2, d3, d4, d5, d6, d7, d8, d9)
!
!  real(8) :: rbbl
!  real(8), intent(in) :: tb, sb, c0, c1, c2, c3, c4, c5, c6, &
!    &                    d0, d1, d2, d3, d4, d5, d6, d7, d8, d9
!
!  rbbl = &
!     &   (c0 + (c1 + (c2 + c3 * tb) * tb) * tb &
!     &       + (c4 + c5 * tb + c6 * sb) * sb) &
!     & / (d0 + (d1 + (d2 + (d3 + d4 * tb) * tb) * tb) * tb &
!     &       + (d5 + (d6 + d7 * tb * tb) * tb &
!     &             + (d8 + d9 * tb * tb) * sqrt(sb)) * sb) &
!     & - 1.d3
!
!  return
!
!end function rbbl
#endif

end module ctnuv
