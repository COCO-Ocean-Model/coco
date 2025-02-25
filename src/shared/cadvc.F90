module cadvc

! --- information -----------------------------------------------------
!
!  Advection and metric terms of the equation of motion.
!
!  HISTORY
!     '03.04.22  H.Hasumi: from COCO3.4
!     '07.04.23  H.Hasumi
!     '08.06.11  H.Hasumi: initial/final processing
!     '08.07.10  H.Hasumi: initial/final processing
!     '12.09.02  H.Tatebe: for COCO5.0 in F90
! ---------------------------------------------------------------------

  implicit none

  private

  real(8),   parameter,   private  ::  c1 =  23.d0 / 12.d0
  real(8),   parameter,   private  ::  c2 = -16.d0 / 12.d0
  real(8),   parameter,   private  ::  c3 =   5.d0 / 12.d0

  logical,                   save  ::  ofirst
  data ofirst / .true. /

  public  ::  advvel
#ifdef OPT_BBL
  public  ::  advvlb
#endif

contains

  subroutine advvel(                                                  &
         &      gx,     gy,     xx,     yy,                           &
         &      uy,     ux,     vy,     vx,                           &
         &      hx,     hy,                                           &
         &    uadv,   vadv,   wadv )

    use zocdim,  only :                                               &
         &  nxydim,  nxdim,  nydim,  nzdim,                           &
         &    kstr,   kend,     kz,                                   &
         &   ijvstr, ijvend,                                          &
         &      le,     lw,     ln,     ls,    lne,                   &
         &   oinit,  ofinal
    use zocgrd,  only :                                               &
         &    zbot,                                                   &  
         &      dy,    dzv,                                           &  
         &      rx,    rym,     rs,                                   & 
         &     hxu,    hyu,   hxyu,   hyxu,                           &
         &     rxu,    ryu   
    use zocmsk,  only :  amskv,  amfvx,  amfvy
    use brstt
    
    implicit none

    real(8),   intent(inout)  ::     gx(nxydim,nzdim),    gy(nxydim,nzdim)
    real(8),   intent(inout)  ::     xx(nxydim,nzdim),    yy(nxydim,nzdim)
    real(8),   intent(in)     ::     uy(nxydim,nzdim),    ux(nxydim,nzdim)
    real(8),   intent(in)     ::     vy(nxydim,nzdim),    vx(nxydim,nzdim)
    real(8),   intent(in)     ::     hx(nxydim),          hy(nxydim)
    real(8),   intent(in)     ::   uadv(nxydim,nzdim),  vadv(nxydim,nzdim)
    real(8),   intent(in)     ::   wadv(nxydim,nzdim)

!---- local variables
    real(8),     save  ::    xx1(nxydim,nzdim),  yy1(nxydim,nzdim) 
    real(8),     save  ::    xx2(nxydim,nzdim),  yy2(nxydim,nzdim) 
    real(8),     save  ::    xx3(nxydim,nzdim),  yy3(nxydim,nzdim) 
    integer(4),  save  ::  ncall
    logical,     save  ::  oeof   
    data ncall / 0 /

    real(8)            ::    fux(nxydim,nzdim),  fvx(nxydim,nzdim) 
    real(8)            ::    fuy(nxydim,nzdim),  fvy(nxydim,nzdim) 
    real(8)            ::    fuz(nxydim,nzdim),  fvz(nxydim,nzdim) 
    real(8)            ::     rz(nxydim,nzdim),  rzm(nxydim,nzdim) 
    real(8)            ::  hvbot(nxydim),     hvbotx(nxydim) 
    real(8)            ::    div(nxydim,kstr:kstr+kz-1)
    integer(4)         ::     ij,      k

    if ( oinit ) then
       do k = 1, nzdim
          do ij = 1, nxydim
             xx1(ij,k) = 0.d0
             xx2(ij,k) = 0.d0
             xx3(ij,k) = 0.d0
             yy1(ij,k) = 0.d0
             yy2(ij,k) = 0.d0
             yy3(ij,k) = 0.d0
          end do
       end do
#ifdef OPT_TRIPOLE
       call rstadd(xx2, oeof, nxdim, nydim, nzdim, 'XX2', 'OCN',      &
    &                                            -1.d0,  -1,  -1)
       call rstadd(xx3, oeof, nxdim, nydim, nzdim, 'XX3', 'OCN',      &
    &                                            -1.d0,  -1,  -1)
       call rstadd(yy2, oeof, nxdim, nydim, nzdim, 'YY2', 'OCN',      &
    &                                            -1.d0,  -1,  -1)
       call rstadd(yy3, oeof, nxdim, nydim, nzdim, 'YY3', 'OCN',      &
    &                                            -1.d0,  -1,  -1)
#else
       call rstadd(xx2, oeof, nxdim, nydim, nzdim, 'XX2', 'OCN')
       call rstadd(xx3, oeof, nxdim, nydim, nzdim, 'XX3', 'OCN')
       call rstadd(yy2, oeof, nxdim, nydim, nzdim, 'YY2', 'OCN')
       call rstadd(yy3, oeof, nxdim, nydim, nzdim, 'YY3', 'OCN')
#endif
       return
    end if

    if ( ofinal ) then
       call finadd(xx2, nxdim, nydim, nzdim, 'XX2', 'OCN')
       call finadd(xx3, nxdim, nydim, nzdim, 'XX3', 'OCN')
       call finadd(yy2, nxdim, nydim, nzdim, 'YY2', 'OCN')
       call finadd(yy3, nxdim, nydim, nzdim, 'YY3', 'OCN')
       return
    end if

#ifdef ADF_
!$acc data copy(gx, gy, xx, yy) 
!$acc data copyin(uy, ux, vy, vx, uadv, vadv, wadv) 
!$acc data copyin(dy) 
!$acc data copyin(rym) 
!$acc data copyin(rs) 
!$acc data copyin(dzv) 
!$acc data copyin(rxu, ryu)  
!$acc data copyin(amskv)  
!$acc data copyin(hxu) 
!$acc data copyin(hyu)     
!$acc data copyin(hxyu, hyxu) 
!$acc data copy(xx1, yy1, xx2, yy2, xx3, yy3)  
!$acc data create(fux, fvx, fuy, fvy, fuz, fvz) 
!$acc data create(hvbot, hvbotx) 
!$acc data create(rz) 
!$acc data create(div) 
#endif

#ifdef OMP_
!$omp parallel do private(k, ij)
#elif  ACC_
!$acc data copy(fux, fvx, fuy, fvy, fuz, fvz) 
!$acc kernels 
#endif
    do k = 1, nzdim
       do ij = 1, nxydim
          fux(ij, k) = 0.d0
          fvx(ij, k) = 0.d0
          fuy(ij, k) = 0.d0
          fvy(ij, k) = 0.d0
          fuz(ij, k) = 0.d0
          fvz(ij, k) = 0.d0
       end do
    end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copy(fux, fvx, fuy, fvy, fuz, fvz)
#endif

#ifdef OMP_
!$omp parallel do private(ij)
#elif  ACC_
!$acc data copyin(hy, hx) 
!$acc data copyin(rym, dy) 
!$acc data copy(hvbot, hvbotx) 
!$acc kernels 
#endif
    do ij = ijvstr, ijvend
       hvbot(ij) = (  (hy(ij)    + hy(ij+le) ) * dy(ij)               &
    &               + (hy(ij+ln) + hy(ij+lne)) * dy(ij+ln)) *         &
    &              rym(ij) * 0.25d0                                   &
    &            + zbot
       hvbot(ij) = 1.d0 / hvbot(ij)
       hvbotx(ij) = (  (hx(ij)    + hx(ij+le) ) * dy(ij)              &
    &                + (hx(ij+ln) + hx(ij+lne)) * dy(ij+ln)) *        &
    &               rym(ij) * 0.25d0                                  &
    &             + zbot
       hvbotx(ij) = 1.d0 / hvbotx(ij)
    end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copyin(hy, hx)
!$acc end data ! copyin(rym, dy)
!$acc end data ! copy(hvbot, hvbotx)
#endif

#ifdef OMP_
!$omp parallel do private(k, ij)
#elif  ACC_
!$acc data copyin(hvbotx) 
!$acc data copyin(rs) 
!$acc data copy(rz) 
!$acc kernels 
#endif
    do k = kstr, kstr+kz-1
       do ij = ijvstr, ijvend
          rz (ij, k) = 1.d0 * rs (k) * hvbotx(ij)
       end do
    end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copy(rz)
!$acc end data ! copyin(rs)
!$acc end data ! copyin(hvbotx)
#endif

#ifdef OMP_
!$omp parallel do private(k, ij)
#elif  ACC_
!$acc data copy(rz) 
!$acc data copyin(dzv) 
!$acc kernels 
#endif
    do k = kstr+kz, kend
       do ij = ijvstr, ijvend
          rz (ij, k) = 1.d0 / dzv(ij, k)
       end do
    end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copyin(dzv)
!$acc end data ! copy(rz)
#endif

#ifdef OMP_
!$omp parallel do private(k, ij)
#elif  ACC_
!$acc data copyin(uadv, vadv, wadv)  
!$acc data copyin(hxu, hyu, rym, rs, rxu, ryu, amskv)  
!$acc data copy(div) 
!$acc kernels 
#endif
    do k = kstr, kstr+kz-2
       do ij = ijvstr, ijvend
          div(ij, k) = (  (  (  uadv(ij+le, k) *                      &
    &                          (hyu(ij+le) + hyu(ij))                 &
    &                        - uadv(ij, k) *                          &
    &                          (hyu(ij) + hyu(ij+lw))) *              &
    &                       0.5d0 * rx                                &
    &                     + (  vadv(ij+ln, k) *                       &
    &                          (hxu(ij+ln) + hxu(ij))                 &
    &                        - vadv(ij, k) *                          &
    &                          (hxu(ij) + hxu(ij+ls))) *              &
    &                       0.5d0 * rym(ij)) * rxu(ij) * ryu(ij)      &
    &                  + (wadv(ij, k) - wadv(ij, k+1)) * rs(k)) *     &
    &                 amskv(ij, k)
       end do
    end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copyin(uadv, vadv, wadv) 
!$acc end data ! copyin(hxu, hyu, rym, rs, rxu, ryu, amskv) 
!$acc end data ! copy(div)
#endif
    k = kstr+kz-1
#ifdef OMP_
!$omp parallel do private(ij)
#elif  ACC_
!$acc data copyin(amskv) 
!$acc data copyin(hxu, hyu, rxu, ryu, rs, rym) 
!$acc data copyin(hvbot) 
!$acc data copy(div) 
!$acc data copyin(uadv, vadv, wadv) 
!$acc kernels 
#endif
    do ij = ijvstr, ijvend
       div(ij, k) = (  (  (  uadv(ij+le, k) *                         &
    &                       (hyu(ij+le) + hyu(ij))                    &
    &                     - uadv(ij, k) *                             &
    &                       (hyu(ij) + hyu(ij+lw))) *                 &
    &                    0.5d0 * rx                                   &
    &                  + (  vadv(ij+ln, k) *                          &
    &                       (hxu(ij+ln) + hxu(ij))                    &
    &                     - vadv(ij, k) *                             &
    &                       (hxu(ij) + hxu(ij+ls))) *                 &
    &                    0.5d0 * rym(ij)) * rxu(ij) * ryu(ij)         &
    &               + (  wadv(ij, k)                                  &
    &                  - wadv(ij, k+1) * hvbot(ij)) * rs(k)) *        &
    &              amskv(ij, k)
    end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copyin(amskv)
!$acc end data ! copyin(hxu, hyu, rxu, ryu, rs, rym)
!$acc end data ! copyin(hvbot)
!$acc end data ! copy(div)
!$acc end data ! copyin(uadv, vadv, wadv)
#endif

#ifdef OMP_
!$omp parallel do private(ij)
#elif  ACC_
!$acc data copy(fuz, fvz) 
!$acc kernels 
#endif
    do ij = ijvstr, ijvend
       fuz(ij, kstr) = 0.d0
       fvz(ij, kstr) = 0.d0
    end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copy(fuz, fvz)
#endif

#ifdef OMP_
!$omp parallel do private(k, ij)
#elif  ACC_
!$acc data copyin(wadv, vy, uy) 
!$acc data copy(fuz, fvz) 
!$acc kernels 
#endif
    do k = kstr+1, kend
       do ij = ijvstr, ijvend
          fuz(ij, k) =                                                &
    &                - wadv(ij, k) * 0.5d0 *                          &
    &                  (uy(ij, k-1) + uy(ij, k))
          fvz(ij, k) =                                                &
    &                - wadv(ij, k) * 0.5d0 *                          &
    &                  (vy(ij, k-1) + vy(ij, k))
       end do
    end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copy(fuz, fvz)
!$acc end data ! copyin(wadv, vy, uy)
#endif

#ifdef OMP_
!$omp parallel do private(ij, k)
#elif  ACC_
!$acc data copyin(uy, vy, vadv) 
!$acc data copyin(hxu) 
!$acc data copy(fuy, fvy) 
!$acc kernels 
#endif
    do k = kstr, kend
       do ij = ijvstr, ijvend+nxdim
          fuy(ij, k) =                                                &
    &                - vadv(ij, k) * 0.5d0 *                          &
    &                  (  uy(ij, k) * hxu(ij)                         &
    &                   + uy(ij+ls, k) * hxu(ij+ls))
          fvy(ij, k) =                                                &
    &                - vadv(ij, k) * 0.5d0 *                          &
    &                  (  vy(ij, k) * hxu(ij)                         &
    &                   + vy(ij+ls, k) * hxu(ij+ls))
       end do
    end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copyin(uy, vy, vadv)
!$acc end data ! copyin(hxu)
!$acc end data ! copy(fuy, fvy)
#endif

#ifdef OMP_
!$omp parallel do private(k, ij)
#elif  ACC_
!$acc data copyin(uadv, vy, uy)     
!$acc data copyin(hyu)     
!$acc data copy(fux, fvx) 
!$acc kernels 
#endif
    do k = kstr, kend
       do ij = ijvstr, ijvend+1
          fux(ij, k) =                                                &
    &                - uadv(ij, k) * 0.5d0 *                          &
    &                  (  uy(ij, k) * hyu(ij)                         &
    &                   + uy(ij+lw, k) * hyu(ij+lw))
          fvx(ij, k) =                                                &
    &                - uadv(ij, k) * 0.5d0 *                          &
    &                  (  vy(ij, k) * hyu(ij)                         &
    &                   + vy(ij+lw, k) * hyu(ij+lw))
       end do
    end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copyin(uadv, vy, uy)    
!$acc end data ! copyin(hyu)    
!$acc end data ! copy(fux, fvx)
#endif

#ifdef OMP_
!$omp parallel do private(k, ij)
#elif  ACC_
!$acc data copyin(vy, uy) 
!$acc data copyin(fvz, fux, fuy, fvy, fvx, fuz, div) 
!$acc data copyin(rxu, rs, rym, hyxu, amskv, ryu, hxyu) 
!$acc data copy(xx1, yy1) 
!$acc kernels 
#endif
    do k = kstr, kstr+kz-2
       do ij = ijvstr, ijvend
          xx1(ij, k) = (                                              &
    &                 + (  (fux(ij+le, k) - fux(ij, k)) * rx          &
    &                    + (fuy(ij+ln, k) - fuy(ij, k)) * rym(ij)) *  &
    &                   rxu(ij) * ryu(ij)                             & 
    &                 + (fuz(ij, k) - fuz(ij, k+1)) * rs(k)           &
    &                 + uy(ij, k) * div(ij, k)                        &
    &                 + vy(ij, k) * vy(ij, k) * hyxu(ij)              &
    &                 - uy(ij, k) * vy(ij, k) * hxyu(ij)) *           &
    &                amskv(ij, k)
          yy1(ij, k) = (                                              &
    &                 + (  (fvx(ij+le, k) - fvx(ij, k)) * rx          &
    &                    + (fvy(ij+ln, k) - fvy(ij, k)) * rym(ij)) *  &
    &                   rxu(ij) * ryu(ij)                             &
    &                 + (fvz(ij, k) - fvz(ij, k+1)) * rs(k)           &
    &                 + vy(ij, k) * div(ij, k)                        &
    &                 + uy(ij, k) * uy(ij, k) * hxyu(ij)              &
    &                 - uy(ij, k) * vy(ij, k) * hyxu(ij)) *           &
    &                amskv(ij, k)
       end do
    end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copyin(vy, uy)
!$acc end data ! copyin(fvz, fux, fuy, fvy, fvx, fuz, div)
!$acc end data ! copyin(rxu, rs, rym, hyxu, amskv, ryu, hxyu)
!$acc end data ! copy(xx1, yy1)
#endif

    k = kstr+kz-1
#ifdef OMP_
!$omp parallel do private(ij)
#elif  ACC_
!$acc data copy(xx1, yy1) 
!$acc data copyin(uy, vy) 
!$acc data copyin(fvy, fuz, fvx, fux, fuy, fvz, div, hvbot) 
!$acc data copyin(hxyu, ryu, amskv, rym, rs, rxu, hyxu) 
!$acc kernels 
#endif
    do ij = ijvstr, ijvend
       xx1(ij, k) = (                                                 &
    &                + (  (fux(ij+le, k) - fux(ij, k)) * rx           &
    &                   + (fuy(ij+ln, k) - fuy(ij, k)) * rym(ij)) *   &
    &                  rxu(ij) * ryu(ij)                              &
    &                + (  fuz(ij, k)                                  &
    &                   - fuz(ij, k+1) * hvbot(ij)) * rs(k)           &
    &                + uy(ij, k) * div(ij, k)                         &
    &                + vy(ij, k) * vy(ij, k) * hyxu(ij)               &
    &                - uy(ij, k) * vy(ij, k) * hxyu(ij)) *            &
    &               amskv(ij, k)
       yy1(ij, k) = (                                                 &
    &                + (  (fvx(ij+le, k) - fvx(ij, k)) * rx           &
    &                   + (fvy(ij+ln, k) - fvy(ij, k)) * rym(ij)) *   &
    &                  rxu(ij) * ryu(ij)                              &
    &                + (  fvz(ij, k)                                  &
    &                   - fvz(ij, k+1) * hvbot(ij)) * rs(k)           &
    &                + vy(ij, k) * div(ij, k)                         &
    &                + uy(ij, k) * uy(ij, k) * hxyu(ij)               &
    &                - uy(ij, k) * vy(ij, k) * hyxu(ij)) *            &
    &               amskv(ij, k)
    end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copy(xx1, yy1)
!$acc end data ! copyin(uy, vy)
!$acc end data ! copyin(fvy, fuz, fvx, fux, fuy, fvz, div, hvbot)
!$acc end data ! copyin(hxyu, ryu, amskv, rym, rs, rxu, hyxu)
#endif

#ifdef OMP_
!$omp parallel do private(k, ij)
#elif  ACC_
!$acc data copy(xx1, yy1) 
!$acc data copyin(vy, uy) 
!$acc data copyin(fvz, fux, fuy, fvy, fvx, fuz, rz) 
!$acc data copyin(ryu, amskv, rxu, hyxu, rym, hxyu) 
!$acc kernels 
#endif
    do k = kstr+kz, kend
       do ij = ijvstr, ijvend
          xx1(ij, k) = (                                              &
    &                 + (  (  (  fux(ij+le, k)                        &
    &                          - fux(ij   , k)) * rx                  &
    &                       + (  fuy(ij+ln, k)                        &
    &                          - fuy(ij   , k)) * rym(ij)) *          &
    &                      rxu(ij) * ryu(ij)                          &
    &                    + fuz(ij, k) - fuz(ij, k+1)) * rz(ij, k)     &
    &                 + vy(ij, k) * vy(ij, k) * hyxu(ij)              &
    &                 - uy(ij, k) * vy(ij, k) * hxyu(ij)) *           &
    &                amskv(ij, k)
          yy1(ij, k) = (                                              &
    &                 + (  (  (  fvx(ij+le, k)                        &
    &                          - fvx(ij   , k)) * rx                  &
    &                       + (  fvy(ij+ln, k)                        &
    &                          - fvy(ij   , k)) * rym(ij)) *          &
    &                      rxu(ij) * ryu(ij)                          &
    &                    + fvz(ij, k) - fvz(ij, k+1)) * rz(ij, k)     &
    &                 + uy(ij, k) * uy(ij, k) * hxyu(ij)              &
    &                 - uy(ij, k) * vy(ij, k) * hyxu(ij)) *           &
    &                amskv(ij, k)
       end do
    end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copy(xx1, yy1)
!$acc end data ! copyin(vy, uy)
!$acc end data ! copyin(fvz, fux, fuy, fvy, fvx, fuz, rz)
!$acc end data ! copyin(ryu, amskv, rxu, hyxu, rym, hxyu)
#endif

!---- Adams-Bashforth scheme by M. Kurogi
    ncall = ncall + 1
    if ( ( ncall >= 3 ) .or. ( .not. oeof ) ) then
       ncall = 3
#ifdef OMP_
!$omp parallel do private(k, ij)
#elif  ACC_
!$acc data copy(gx, gy) 
!$acc data copyin(xx1, xx2, xx3) 
!$acc data copyin(yy1, yy2, yy3) 
!$acc kernels 
#endif
       do k = kstr, kend
          do ij = ijvstr, ijvend
             gx(ij, k) = gx(ij, k)                                    &
    &             + c1 * xx1(ij,k) + c2 * xx2(ij,k) + c3 * xx3(ij, k) 
             gy(ij, k) =  gy(ij, k)                                   &
    &             + c1 * yy1(ij,k) + c2 * yy2(ij,k) + c3 * yy3(ij, k) 
          end do
       end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copy(gx, gy)
!$acc end data ! copyin(yy1, yy2, yy3)
!$acc end data ! copyin(xx1, xx2, xx3)
#endif
    else 
       if ( ncall == 1 ) then ! forward
#ifdef OMP_
!$omp parallel do private(k, ij)
#elif  ACC_
!$acc data copyin(xx1, yy1) 
!$acc data copy(gx, gy) 
!$acc kernels 
#endif
          do k = kstr, kend
             do ij = ijvstr, ijvend
                gx(ij, k) = gx(ij, k) + xx1(ij, k)
                gy(ij, k) = gy(ij, k) + yy1(ij, k)
             end do
          end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copyin(xx1, yy1)
!$acc end data ! copy(gx, gy)
#endif
       else if ( ncall == 2 ) then !2nd order AB
#ifdef OMP_
!$omp parallel do private(k, ij)
#elif  ACC_
!$acc data copyin(xx1, xx2, yy1, yy2) 
!$acc data copy(gx, gy) 
!$acc kernels 
#endif
          do k = kstr, kend
             do ij = ijvstr, ijvend
                gx(ij, k) = gx(ij, k)                                 &
    &                 + 1.5d0 * xx1(ij,k) - 0.5d0 * xx2(ij,k)
                gy(ij, k) =  gy(ij, k)                                &
    &                 + 1.5d0 * yy1(ij,k) - 0.5d0 * yy2(ij,k)
             end do
          end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copyin(xx1, xx2, yy1, yy2)
!$acc end data ! copy(gx, gy)
#endif
       end if
    end if

#ifdef OMP_
!$omp parallel do private(k, ij)
#elif  ACC_
!$acc data copyin(xx1, yy1) 
!$acc data copy(xx2, yy2, xx3, yy3) 
!$acc kernels 
#endif
    do k = kstr, kend
       do ij = ijvstr, ijvend
         
          xx3(ij, k) = xx2(ij, k) ! n-1 > n-2
          yy3(ij, k) = yy2(ij, k)
         
          xx2(ij, k) = xx1(ij, k) ! n   > n-1
          yy2(ij, k) = yy1(ij, k)

       end do
    end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copyin(xx1, yy1)
!$acc end data ! copy(xx2, yy2, xx3, yy3)
#endif

#ifdef ADF_
!$acc end data ! copy(gx, gy, xx, yy)
!$acc end data ! copyin(uy, ux, vy, vx, uadv, vadv, wadv)
!$acc end data ! copyin(dy)
!$acc end data ! copyin(rym)
!$acc end data ! copyin(rs)
!$acc end data ! copyin(dzv)
!$acc end data ! copyin(rxu, ryu) 
!$acc end data ! copyin(amskv) 
!$acc end data ! copyin(hxu)
!$acc end data ! copyin(hyu)    
!$acc end data ! copyin(hxyu, hyxu)
!$acc end data ! copy(fux, fvx, fuy, fvy, fuz, fvz) 
!$acc end data ! copy(hvbot, hvbotx) 
!$acc end data ! copy(rz) 
!$acc end data ! copy(div) 
!$acc end data ! copy(xx1, yy1, xx2, yy2, xx3, yy3)  
#endif
   
  end subroutine advvel

#ifdef OPT_BBL

! --- information -----------------------------------------------------
!
!  Advection term of the BBL momentum eqs.
!
! ---------------------------------------------------------------------
  subroutine advvlb(                                                  &
         &      gx,     gy,     xx,     yy,                           &
         &      uy,     ux,     vy,     vx,                           &
         &    uadv,   vadv,   wadv)

    use zocdim,  only :                                               &
         &  nxydim,  nxdim,  nydim,  nzdim,                           &
         &    kstr,   kend,     kz,                                   &
         &  ijvstr, ijvend,                                           &
         &      le,     lw,     ln,     ls,                           &
         &   oinit,  ofinal
    use zocgrd,  only :                                               &
         &     dzv,                                                   &  
         &      rx,    rym,                                           &
         &     hxu,    hyu,   hxyu,   hyxu,                           &
         &     rxu,    ryu   
    use zocmsk,  only :  amskvb,  nbotv
    use brstt
    
    implicit none

    real(8),   intent(inout)  ::     gx(nxydim,nzdim),    gy(nxydim,nzdim)
    real(8),   intent(inout)  ::     xx(nxydim,nzdim),    yy(nxydim,nzdim)
    real(8),   intent(in)     ::     uy(nxydim,nzdim),    ux(nxydim,nzdim)
    real(8),   intent(in)     ::     vy(nxydim,nzdim),    vx(nxydim,nzdim)
    real(8),   intent(in)     ::   uadv(nxydim,nzdim),  vadv(nxydim,nzdim)
    real(8),   intent(in)     ::   wadv(nxydim,nzdim)

!---- local variables
    real(8),     save  ::    xx1(nxydim),  yy1(nxydim) 
    real(8),     save  ::    xx2(nxydim),  yy2(nxydim) 
    real(8),     save  ::    xx3(nxydim),  yy3(nxydim) 
    integer(4),  save  ::  ncall
    logical,     save  ::  oeof   
    data ncall / 0 /

    real(8)            ::    fux(nxydim),  fvx(nxydim)
    real(8)            ::    fuy(nxydim),  fvy(nxydim) 
    real(8)            ::    fuz(nxydim),  fvz(nxydim)
    real(8),    save   ::     rz(nxydim)
    integer(4)         ::     ij,     kup

    if ( oinit ) then
       do ij = 1, nxydim
          xx1(ij) = 0.d0
          xx2(ij) = 0.d0
          xx3(ij) = 0.d0
          yy1(ij) = 0.d0
          yy2(ij) = 0.d0
          yy3(ij) = 0.d0
       end do
#ifdef OPT_TRIPOLE
       call rstadd(xx2, oeof, nxdim, nydim, 1, 'XX2', 'SFC',          &
    &                                        -1.d0,  -1,  -1)
       call rstadd(xx3, oeof, nxdim, nydim, 1, 'XX3', 'SFC',          &
    &                                        -1.d0,  -1,  -1)
       call rstadd(yy2, oeof, nxdim, nydim, 1, 'YY2', 'SFC',          &
    &                                        -1.d0,  -1,  -1)
       call rstadd(yy3, oeof, nxdim, nydim, 1, 'YY3', 'SFC',          &
    &                                        -1.d0,  -1,  -1)
#else
       call rstadd(xx2, oeof, nxdim, nydim, 1, 'XX2', 'SFC')
       call rstadd(xx3, oeof, nxdim, nydim, 1, 'XX3', 'SFC')
       call rstadd(yy2, oeof, nxdim, nydim, 1, 'YY2', 'SFC')
       call rstadd(yy3, oeof, nxdim, nydim, 1, 'YY3', 'SFC')
#endif
       return
    end if

    if ( ofinal ) then
       call finadd(xx2, nxdim, nydim, 1, 'XX2', 'SFC')
       call finadd(xx3, nxdim, nydim, 1, 'XX3', 'SFC')
       call finadd(yy2, nxdim, nydim, 1, 'YY2', 'SFC')
       call finadd(yy3, nxdim, nydim, 1, 'YY3', 'SFC')
       return
    end if

    if ( ofirst ) then
       ofirst = .false.
       do ij = 1, nxydim
          rz(ij) = 1.d0 / dzv(ij, kend)
       end do
    end if

#ifdef ADF_
!$acc data create(fux, fuy, fuz, fvy, fvx, fvz) 
!$acc data copyin(uadv, vadv, wadv) 
!$acc data copyin(vy, uy)  
!$acc data copy(gx, gy) 
!$acc data copyin(nbotv) 
!$acc data copyin(hxu)  
!$acc data copyin(hyu) 
!$acc data copyin(hxyu, rym, hyxu, rxu, amskvb, ryu) 

!$acc data copyin(rz) 
!$acc data copy(xx3, yy1, yy3, yy2, xx2, xx1) 
#endif

#ifdef OMP_
!$omp parallel do private(ij)
#elif  ACC_
!$acc data copy(fux, fuy, fvy, fvx) 
!$acc kernels 
#endif
    do ij = 1, nxydim
       fux(ij) = 0.d0
       fvx(ij) = 0.d0
       fuy(ij) = 0.d0
       fvy(ij) = 0.d0
    end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copy(fux, fuy, fvy, fvx)
#endif

#ifdef OMP_
!$omp parallel do private(ij, kup)
#elif  ACC_
!$acc data copyin(wadv, vy, uy) 
!$acc data copyin(nbotv) 
!$acc data copy(fvz, fuz) 
!$acc kernels 
#endif
    do ij = ijvstr, ijvend
       kup = max(nbotv(ij)-1, 1)                                      
       fuz(ij) = - wadv(ij, kend) * 0.5d0 *                           &
    &             (uy(ij, kup) + uy(ij, kend))
       fvz(ij) = - wadv(ij, kend) * 0.5d0 *                           &
    &             (vy(ij, kup) + vy(ij, kend))
    end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copyin(wadv, vy, uy)
!$acc end data ! copyin(nbotv)
!$acc end data ! copy(fvz, fuz)
#endif

#ifdef OMP_
!$omp parallel do private(ij)
#elif  ACC_
!$acc data copy(fuy, fvy) 
!$acc data copyin(vadv, vy, uy)  
!$acc data copyin(hxu)  
!$acc kernels 
#endif
    do ij = ijvstr, ijvend+nxdim
       fuy(ij) = - vadv(ij, kend) * 0.5d0 *                           &
    &             (  uy(ij, kend) * hxu(ij)                           &
    &              + uy(ij+ls, kend) * hxu(ij+ls))
       fvy(ij) = - vadv(ij, kend) * 0.5d0 *                           &
    &             (  vy(ij, kend) * hxu(ij)                           &
    &              + vy(ij+ls, kend) * hxu(ij+ls))
    end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copy(fuy, fvy)
!$acc end data ! copyin(vadv, vy, uy) 
!$acc end data ! copyin(hxu) 
#endif

#ifdef OMP_
!$omp parallel do private(ij)
#elif  ACC_
!$acc data copy(fux, fvx) 
!$acc data copyin(uadv, vy, uy) 
!$acc data copyin(hyu) 
!$acc kernels 
#endif
    do ij = ijvstr, ijvend+1
       fux(ij) = - uadv(ij, kend) * 0.5d0 *                           &
    &             (  uy(ij, kend) * hyu(ij)                           &
    &              + uy(ij+lw, kend) * hyu(ij+lw))
       fvx(ij) = - uadv(ij, kend) * 0.5d0 *                           &
    &             (  vy(ij, kend) * hyu(ij)                           &
    &              + vy(ij+lw, kend) * hyu(ij+lw))
    end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copy(fux, fvx)
!$acc end data ! copyin(uadv, vy, uy)
!$acc end data ! copyin(hyu)
#endif

#ifdef OMP_
!$omp parallel do private(ij)
#elif  ACC_
!$acc data copy(xx1, yy1) 
!$acc data copyin(fvy, fuz, fvx, fux, fuy, fvz, rz) 
!$acc data copyin(uy, vy) 
!$acc data copyin(hxyu, rym, hyxu, rxu, amskvb, ryu) 
!$acc kernels 
#endif
    do ij = ijvstr, ijvend
       xx1(ij) = (                                                    &
    &                 + (  (  (fux(ij+le) - fux(ij)) * rx             &
    &                       + (fuy(ij+ln) - fuy(ij)) * rym(ij)) *     &
    &                      rxu(ij) * ryu(ij)                          &
    &                    + fuz(ij)) * rz(ij)                          &
    &                 + vy(ij, kend) * vy(ij, kend) * hyxu(ij)        &
    &                 - uy(ij, kend) * vy(ij, kend) * hxyu(ij)) *     &
    &                amskvb(ij) 
       yy1(ij) = (                                                    &
    &                 + (  (  (fvx(ij+le) - fvx(ij)) * rx             &
    &                       + (fvy(ij+ln) - fvy(ij)) * rym(ij)) *     &
    &                      rxu(ij) * ryu(ij)                          &
    &                    + fvz(ij)) * rz(ij)                          &
    &                 + uy(ij, kend) * uy(ij, kend) * hxyu(ij)        &
    &                 - uy(ij, kend) * vy(ij, kend) * hyxu(ij)) *     &
    &                amskvb(ij)
    end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copy(xx1, yy1)
!$acc end data ! copyin(fvy, fuz, fvx, fux, fuy, fvz, rz)
!$acc end data ! copyin(uy, vy)
!$acc end data ! copyin(hxyu, rym, hyxu, rxu, amskvb, ryu)
#endif

!---- Adams-Bashforth scheme by M. Kurogi
    ncall = ncall + 1

    if ( ( ncall >= 3 ) .or. ( .not. oeof ) ) then
       ncall = 3
#ifdef OMP_
!$omp parallel do private(ij)
#elif  ACC_
!$acc data copy(gx, gy) 
!$acc data copyin(xx3, yy1, yy3, yy2, xx2, xx1) 
!$acc kernels 
#endif
       do ij = ijvstr, ijvend
          gx(ij, kend) =  gx(ij, kend)                                &
    &        + c1 * xx1(ij) + c2 * xx2(ij) + c3 * xx3(ij) 
          gy(ij, kend) =  gy(ij, kend)                                &
    &        + c1 * yy1(ij) + c2 * yy2(ij) + c3 * yy3(ij) 
       end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copyin(xx3, yy1, yy3, yy2, xx2, xx1) 
!$acc end data ! copy(gx, gy)
#endif
    else
       if ( ncall == 1 ) then ! forward
#ifdef OMP_
!$omp parallel do private(ij)
#elif  ACC_
!$acc data copy(gx, gy) 
!$acc data copyin(yy1, xx1) 
!$acc kernels 
#endif
          do ij = ijvstr, ijvend
             gx(ij, kend) =  gx(ij, kend) + xx1(ij)
             gy(ij, kend) =  gy(ij, kend) + yy1(ij)
          end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copy(gx, gy)
!$acc end data ! copyin(yy1, xx1)
#endif
       else if (ncall == 2 ) then !2nd order AB
#ifdef OMP_
!$omp parallel do private(ij)
#elif  ACC_
!$acc data copy(gx, gy) 
!$acc data copyin(yy1, yy2, xx2, xx1) 
!$acc kernels 
#endif
          do ij = ijvstr, ijvend
             gx(ij, kend) = gx(ij, kend) + 1.5d0 * xx1(ij) - 0.5d0 * xx2(ij)
             gy(ij, kend) = gy(ij, kend) + 1.5D0 * yy1(ij) - 0.5d0 * yy2(ij)
          end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copyin(yy1, yy2, xx2, xx1)
!$acc end data ! copy(gx, gy)
#endif
       end if
    end if

#ifdef OMP_
!$omp parallel do private(ij)
#elif  ACC_
!$acc data copy(xx3, yy3, xx2, yy2) 
!$acc data copyin(xx1, yy1) 
!$acc kernels 
#endif
    do ij = ijvstr, ijvend
    
       xx3(ij) = xx2(ij) ! n-1 > n-2
       yy3(ij) = yy2(ij)

       xx2(ij) = xx1(ij) ! n   > n-1
       yy2(ij) = yy1(ij)

    end do
#ifdef OMP_
!$omp end parallel do
#elif  ACC_
!$acc end kernels
!$acc end data ! copy(xx3, yy3, xx2, yy2)
!$acc end data ! copyin(xx1, yy1)
#endif

#ifdef ADF_
!$acc end data ! create(fux, fuy, fuz, fvy, fvx, fvz)
!$acc end data ! copyin(uadv, vadv, wadv)
!$acc end data ! copyin(vy, uy) 
!$acc end data ! copy(gx, gy)
!$acc end data ! copyin(nbotv)
!$acc end data ! copyin(hxu) 
!$acc end data ! copyin(hyu)

!$acc end data ! copyin(rz)
!$acc end data ! copy(xx3, yy1, yy3, yy2, xx2, xx1)
!$acc end data ! copyin(hxyu, rym, hyxu, rxu, amskvb, ryu)
#endif

  end subroutine advvlb

#endif

end module cadvc

