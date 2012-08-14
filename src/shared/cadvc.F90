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
!     '12.08.14  H.Tatebe: for COCO5.0 in F90
! ---------------------------------------------------------------------

  implicit none

  private
  public  ::  advvel, advvlb

  logical,     save  ::  ofirst
  data ofirst / .true. /

contains

  subroutine advvel(                                                  &
         &      gx,     gy,     xx,     yy,                           &
         &      uy,     ux,     vy,     vx,                           &
         &      hx,     hy,                                           &
         &    uadv,   vadv,   wadv)

    use zocdim,  only :                                               &
         &  nxydim,  nxdim,  nzdim,                                   &
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
    
    implicit none

    real(8),   intent(inout)  ::     gx(nxydim,nzdim),    gy(nxydim,nzdim)
    real(8),   intent(inout)  ::     xx(nxydim,nzdim),    yy(nxydim,nzdim)
    real(8),   intent(in)     ::     uy(nxydim,nzdim),    ux(nxydim,nzdim)
    real(8),   intent(in)     ::     vy(nxydim,nzdim),    vx(nxydim,nzdim)
    real(8),   intent(in)     ::     hx(nxydim),          hy(nxydim)
    real(8),   intent(in)     ::   uadv(nxydim,nzdim),  vadv(nxydim,nzdim)
    real(8),   intent(in)     ::   wadv(nxydim,nzdim)

!---- local variables
    real(8)            ::    fux(nxydim,nzdim),  fvx(nxydim,nzdim) 
    real(8)            ::    fuy(nxydim,nzdim),  fvy(nxydim,nzdim) 
    real(8)            ::    fuz(nxydim,nzdim),  fvz(nxydim,nzdim) 
    real(8)            ::     rz(nxydim,nzdim),  rzm(nxydim,nzdim) 
    real(8)            ::  hvbot(nxydim),     hvbotx(nxydim) 
    real(8)            ::    div(nxydim,kstr:kstr+kz-1)
    integer(4)         ::     ij,      k

    if ( oinit .or. ofinal ) then
       return
    end if

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

    do k = kstr, kstr+kz-1
       do ij = ijvstr, ijvend
          rz (ij, k) = 1.d0 * rs (k) * hvbotx(ij)
       end do
    end do
    do k = kstr+kz, kend
       do ij = ijvstr, ijvend
          rz (ij, k) = 1.d0 / dzv(ij, k)
       end do
    end do
    
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
    k = kstr+kz-1
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

    do ij = ijvstr, ijvend
       fuz(ij, kstr) = 0.d0
       fvz(ij, kstr) = 0.d0
    end do

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

    do k = kstr, kstr+kz-2
       do ij = ijvstr, ijvend
          gx(ij, k) = (  gx(ij, k)                                    &
    &                 + (  (fux(ij+le, k) - fux(ij, k)) * rx          &
    &                    + (fuy(ij+ln, k) - fuy(ij, k)) * rym(ij)) *  &
    &                   rxu(ij) * ryu(ij)                             & 
    &                 + (fuz(ij, k) - fuz(ij, k+1)) * rs(k)           &
    &                 + uy(ij, k) * div(ij, k)                        &
    &                 + vy(ij, k) * vy(ij, k) * hyxu(ij)              &
    &                 - uy(ij, k) * vy(ij, k) * hxyu(ij)) *           &
    &                amskv(ij, k)
          gy(ij, k) = (  gy(ij, k)                                    &
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

    k = kstr+kz-1
    do ij = ijvstr, ijvend
       gx(ij, k) = (  gx(ij, k)                                       &
    &                + (  (fux(ij+le, k) - fux(ij, k)) * rx           &
    &                   + (fuy(ij+ln, k) - fuy(ij, k)) * rym(ij)) *   &
    &                  rxu(ij) * ryu(ij)                              &
    &                + (  fuz(ij, k)                                  &
    &                   - fuz(ij, k+1) * hvbot(ij)) * rs(k)           &
    &                + uy(ij, k) * div(ij, k)                         &
    &                + vy(ij, k) * vy(ij, k) * hyxu(ij)               &
    &                - uy(ij, k) * vy(ij, k) * hxyu(ij)) *            &
    &               amskv(ij, k)
       gy(ij, k) = (  gy(ij, k)                                       &
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

    do k = kstr+kz, kend
       do ij = ijvstr, ijvend
          gx(ij, k) = (  gx(ij, k)                                    &
    &                 + (  (  (  fux(ij+le, k)                        &
    &                          - fux(ij   , k)) * rx                  &
    &                       + (  fuy(ij+ln, k)                        &
    &                          - fuy(ij   , k)) * rym(ij)) *          &
    &                      rxu(ij) * ryu(ij)                          &
    &                    + fuz(ij, k) - fuz(ij, k+1)) * rz(ij, k)     &
    &                 + vy(ij, k) * vy(ij, k) * hyxu(ij)              &
    &                 - uy(ij, k) * vy(ij, k) * hxyu(ij)) *           &
    &                amskv(ij, k)
         gy(ij, k) = (  gy(ij, k)                                     &
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
         &  nxydim,  nxdim,  nzdim,                                   &
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
    
    implicit none

    real(8),   intent(inout)  ::     gx(nxydim,nzdim),    gy(nxydim,nzdim)
    real(8),   intent(inout)  ::     xx(nxydim,nzdim),    yy(nxydim,nzdim)
    real(8),   intent(in)     ::     uy(nxydim,nzdim),    ux(nxydim,nzdim)
    real(8),   intent(in)     ::     vy(nxydim,nzdim),    vx(nxydim,nzdim)
    real(8),   intent(in)     ::   uadv(nxydim,nzdim),  vadv(nxydim,nzdim)
    real(8),   intent(in)     ::   wadv(nxydim,nzdim)

!---- local variables
    real(8)            ::    fux(nxydim),  fvx(nxydim)
    real(8)            ::    fuy(nxydim),  fvy(nxydim) 
    real(8)            ::    fuz(nxydim),  fvz(nxydim)
    real(8),    save   ::     rz(nxydim)
    integer(4)         ::     ij,     kup

    if ( oinit .or. ofinal ) then
       return
    end if

    if ( ofirst ) then
       ofirst = .false.
       do ij = 1, nxydim
          rz   (ij) = 1.d0 / dzv(ij, kend)
       end do
    end if

    do ij = 1, nxydim
       fux(ij) = 0.d0
       fvx(ij) = 0.d0
       fuy(ij) = 0.d0
       fvy(ij) = 0.d0
    end do

    do ij = ijvstr, ijvend
       kup = max(nbotv(ij)-1, 1)                                      
       fuz(ij) = - wadv(ij, kend) * 0.5d0 *                           &
    &             (uy(ij, kup) + uy(ij, kend))
       fvz(ij) = - wadv(ij, kend) * 0.5d0 *                           &
    &             (vy(ij, kup) + vy(ij, kend))
    end do

    do ij = ijvstr, ijvend+nxdim
       fuy(ij) = - vadv(ij, kend) * 0.5d0 *                           &
    &             (  uy(ij, kend) * hxu(ij)                           &
    &              + uy(ij+ls, kend) * hxu(ij+ls))
       fvy(ij) = - vadv(ij, kend) * 0.5d0 *                           &
    &             (  vy(ij, kend) * hxu(ij)                           &
    &              + vy(ij+ls, kend) * hxu(ij+ls))
    end do

    do ij = ijvstr, ijvend+1
       fux(ij) = - uadv(ij, kend) * 0.5d0 *                           &
    &             (  uy(ij, kend) * hyu(ij)                           &
    &              + uy(ij+lw, kend) * hyu(ij+lw))
       fvx(ij) = - uadv(ij, kend) * 0.5d0 *                           &
    &             (  vy(ij, kend) * hyu(ij)                           &
    &              + vy(ij+lw, kend) * hyu(ij+lw))
    end do

    do ij = ijvstr, ijvend
       gx(ij, kend) = (  gx(ij, kend)                                 &
    &                 + (  (  (fux(ij+le) - fux(ij)) * rx             &
    &                       + (fuy(ij+ln) - fuy(ij)) * rym(ij)) *     &
    &                      rxu(ij) * ryu(ij)                          &
    &                    + fuz(ij)) * rz(ij)                          &
    &                 + vy(ij, kend) * vy(ij, kend) * hyxu(ij)        &
    &                 - uy(ij, kend) * vy(ij, kend) * hxyu(ij)) *     &
    &                amskvb(ij) 
       gy(ij, kend) = (  gy(ij, kend)                                 &
    &                 + (  (  (fvx(ij+le) - fvx(ij)) * rx             &
    &                       + (fvy(ij+ln) - fvy(ij)) * rym(ij)) *     &
    &                      rxu(ij) * ryu(ij)                          &
    &                    + fvz(ij)) * rz(ij)                          &
    &                 + uy(ij, kend) * uy(ij, kend) * hxyu(ij)        &
    &                 - uy(ij, kend) * vy(ij, kend) * hyxu(ij)) *     &
    &                amskvb(ij)
    end do

  end subroutine advvlb

#endif

end module cadvc

