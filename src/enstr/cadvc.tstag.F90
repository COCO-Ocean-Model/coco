module cadvc

! --- information -----------------------------------------------------
!
!  Advection and metric terms of the equation of motion.
!
!  HISTORY
!     '03.05.13  H.Hasumi: from COCO3.4
!     '06.04.19  T.Suzuki:
!     '07.04.23  H.Hasumi
!     '08.06.11  H.Hasumi: initial/final processing
!     '08.07.10  H.Hasumi: initial/final processing
!     '10.04.14  M.Kurogi: staggered time stepping
!     '10.04.14  M.Kurogi: (COCO4.4 tripolar code by Dr. Suzuki)
!     '12.10.23. T.Suzuki: for COCO5.0
!
! ---------------------------------------------------------------------
  implicit none

  private
  public  ::  advvel
#ifdef OPT_BBL
  public  ::  advvlb
#endif

contains

  subroutine advvel(                                                  &
         &            fux,    fuy,   fune,   fuse,                    &
         &            fvx,    fvy,   fvne,   fvse,                    &
         &             gx,     gy,     xx,     yy,                    &
         &             uy,     ux,     vy,     vx,                    &
         &             hx,     hy,                                    &
         &           uadv,   vadv,   wadv)
    use zocdim,  only :                                               &
         &  nxydim,  nxdim,  nydim,   nzdim,                          &
         &    kstr,   kend,     kz,                                   &
         &   ijvstr, ijvend,                                          &
         &      le,     lw,     ln,     ls,    lne,                   &
         &     lsw,    lnw,    lse,                                   &
         &   oinit,  ofinal
    use zocgrd,  only :                                               &
         &    zbot,                                                   &  
         &      dx,     dy,    dzv,                                   &  
         &      rx,    rym,     rs,                                   &
         &     hxt,    hyt,                                           & 
         &     hxu,    hyu,   hxyu,   hyxu,                           &
         &     rxu,    ryu   
#ifdef OPT_BBL
    use zocmsk,  only :  amskv,  amfvx,  amfvy,   nbotv,  amskvb
#else
    use zocmsk,  only :  amskv,  amfvx,  amfvy
#endif
    use bchmk
    use brstt

    real(8), parameter :: c1=23.d0/12.d0, c2=-16.d0/12.d0, c3=5.d0/12.d0
    real(8), save ::      xx1(nxydim, nzdim),     yy1(nxydim, nzdim)
    real(8), save ::      xx2(nxydim, nzdim),     yy2(nxydim, nzdim)
    real(8), save ::      xx3(nxydim, nzdim),     yy3(nxydim, nzdim)
    integer, save :: ncall = 0
    logical, save :: oeof
    real(8),intent(in) ::      uy(nxydim, nzdim),     ux(nxydim, nzdim)
    real(8),intent(in) ::      vy(nxydim, nzdim),     vx(nxydim, nzdim)
    real(8),intent(in) ::      hy(nxydim),            hx(nxydim)
    real(8),intent(in) ::    uadv(nxydim, nzdim),   vadv(nxydim, nzdim)
    real(8),intent(in) ::    wadv(nxydim, nzdim, 9)
    real(8),intent(inout) ::      gx(nxydim, nzdim),     gy(nxydim, nzdim)
    real(8),intent(inout) ::      xx(nxydim, nzdim),     yy(nxydim, nzdim)
    real(8),intent(out) ::     fux(nxydim, nzdim),    fvx(nxydim, nzdim)
    real(8),intent(out) ::     fuy(nxydim, nzdim),    fvy(nxydim, nzdim)
    real(8),intent(out) ::    fune(nxydim, nzdim),   fvne(nxydim, nzdim)
    real(8),intent(out) ::    fuse(nxydim, nzdim),   fvse(nxydim, nzdim)

    real(8),save ::     cxn(nxydim, nzdim),    cxs(nxydim, nzdim)
    real(8),save ::     cye(nxydim, nzdim),    cyw(nxydim, nzdim)
    real(8),save ::     cne(nxydim, nzdim),    cse(nxydim, nzdim)

    real(8) ::     fuz(nxydim, nzdim),    fvz(nxydim, nzdim)
    real(8) ::    fuzu(nxydim, nzdim),   fvzu(nxydim, nzdim)
    real(8) ::    fuzd(nxydim, nzdim),   fvzd(nxydim, nzdim)
    real(8) ::      rz(nxydim, nzdim),    rzm(nxydim, nzdim)
 
    real(8) ::     div(nxydim, kstr:kstr+kz-1)
    real(8) ::   hvbot(nxydim), hvbotx(nxydim)
    integer     ij,      k

    logical, save :: ofirst = .true.

    if (oinit) then
       do k=1,nzdim
       do ij=1,nxydim
          xx1(ij,k)=0.d0
          xx2(ij,k)=0.d0
          xx3(ij,k)=0.d0
          yy1(ij,k)=0.d0
          yy2(ij,k)=0.d0
          yy3(ij,k)=0.d0
       end do
       end do
#ifdef OPT_TRIPOLE
      call rstadd(xx2, oeof, nxdim, nydim, nzdim, 'XX2', 'OCN', &
     &                                           -1.d0,  -1,  -1)
      call rstadd(xx3, oeof, nxdim, nydim, nzdim, 'XX3', 'OCN', &
     &                                           -1.d0,  -1,  -1)
      call rstadd(yy2, oeof, nxdim, nydim, nzdim, 'YY2', 'OCN', &
     &                                           -1.d0,  -1,  -1)
      call rstadd(yy3, oeof, nxdim, nydim, nzdim, 'YY3', 'OCN', &
     &                                           -1.d0,  -1,  -1)
#else
      call rstadd(xx2, oeof, nxdim, nydim, nzdim, 'XX2', 'OCN')
      call rstadd(xx3, oeof, nxdim, nydim, nzdim, 'XX3', 'OCN')
      call rstadd(yy2, oeof, nxdim, nydim, nzdim, 'YY2', 'OCN')
      call rstadd(yy3, oeof, nxdim, nydim, nzdim, 'YY3', 'OCN')
#endif
         return
      end if

      if (ofinal) then
         call finadd(xx2, nxdim, nydim, nzdim, 'XX2', 'OCN')
         call finadd(xx3, nxdim, nydim, nzdim, 'XX3', 'OCN')
         call finadd(yy2, nxdim, nydim, nzdim, 'YY2', 'OCN')
         call finadd(yy3, nxdim, nydim, nzdim, 'YY3', 'OCN')
         return
      end if


    if (ofirst) then
       ofirst = .false.
#ifdef OPT_BBL
    call rmmskv
#endif
    do k = 1, nzdim
       do ij = 1, nxydim
          cxn(ij, k) = 0.d0
          cxs(ij, k) = 0.d0
          cye(ij, k) = 0.d0
          cyw(ij, k) = 0.d0
          cne(ij, k) = 0.d0
          cse(ij, k) = 0.d0
       enddo
    enddo
    do k = kstr, kend
       do ij = ijvstr-nxdim-1, ijvend+nxdim+1
          cxn(ij, k) = amskv(ij, k) * amskv(ij+lw, k) *          &
    &                (3.d0 - amskv(ij+ls, k) - amskv(ij+lsw, k)  &
    &                      + amskv(ij+ls, k) * amskv(ij+lsw, k))
          cxs(ij, k) = amskv(ij+ls, k) * amskv(ij+lsw, k) *      &
    &                (3.d0 - amskv(ij, k) - amskv(ij+lw, k)      &
    &                      + amskv(ij, k) * amskv(ij+lw, k))
          cye(ij, k) = amskv(ij, k) * amskv(ij+ls, k) *          &
    &                (3.d0 - amskv(ij+lw, k) - amskv(ij+lsw, k)  &
    &                      + amskv(ij+lw, k) * amskv(ij+lsw, k))
          cyw(ij, k) = amskv(ij+lw, k) * amskv(ij+lsw, k) *      &
    &                (3.d0 - amskv(ij, k) - amskv(ij+ls, k)      &
    &                      + amskv(ij, k) * amskv(ij+ls, k))
          cne(ij, k) = amskv(ij, k) * amskv(ij+lsw, k) *         &
    &                (3.d0 - amskv(ij+lw, k) - amskv(ij+ls, k))
          cse(ij, k) = amskv(ij+lw, k) * amskv(ij+ls, k) *       &
    &                (3.d0 - amskv(ij, k) - amskv(ij+lsw, k))
       enddo
    enddo

!#ifdef OPT_BBL
!         call admkv1
!#endif
    end if

    do k = 1, nzdim
       do ij = 1, nxydim
          fux (ij, k) = 0.d0
          fuy (ij, k) = 0.d0
          fune(ij, k) = 0.d0
          fuse(ij, k) = 0.d0
          fvx (ij, k) = 0.d0
          fvy (ij, k) = 0.d0
          fvne(ij, k) = 0.d0
          fvse(ij, k) = 0.d0
          fuz (ij, k) = 0.d0
          fuzu(ij, k) = 0.d0
          fuzd(ij, k) = 0.d0
          fvz (ij, k) = 0.d0
          fvzu(ij, k) = 0.d0
          fvzd(ij, k) = 0.d0
       enddo
    enddo

      do ij = ijvstr, ijvend
         hvbot(ij) = (  (hy(ij)    + hy(ij+le) ) * dy(ij) &
     &                + (hy(ij+ln) + hy(ij+lne)) * dy(ij+ln)) * &
     &               rym(ij) * 0.25d0 &
     &             + zbot
         hvbot(ij) = 1.d0 / hvbot(ij)
         hvbotx(ij) = (  (hx(ij)    + hx(ij+le) ) * dy(ij) &
     &                + (hx(ij+ln) + hx(ij+lne)) * dy(ij+ln)) * &
     &               rym(ij) * 0.25d0 &
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

      do k = kstr, kstr+kz-1
         do ij = ijvstr, ijvend+nxdim+1
            fvy(ij, k) = (  cye(ij, k) * vadv(ij, k) * hxt(ij) &
     &                    + cyw(ij+le, k) * vadv(ij+le, k) * hxt(ij+le)) &
     &                   / 6.d0
            fvx(ij, k) = (  cxn(ij, k) * uadv(ij, k) * hyt(ij) &
     &                    + cxs(ij+ln, k) * uadv(ij+ln, k) * hyt(ij+ln)) &
     &                   / 6.d0
            fvne(ij, k) = cne(ij, k) / 6.d0 * &
     &                    (  uadv(ij, k) * dy(ij) * hyt(ij) &
     &                     + vadv(ij, k) * dx * hxt(ij))
            fvse(ij, k) = cse(ij, k) / 6.d0 * &
     &                    (  uadv(ij, k) * dy(ij) * hyt(ij) &
     &                     - vadv(ij, k) * dx * hxt(ij))
         end do
      end do

      do k = kstr+1, kstr+kz-1
         do ij = ijvstr, ijvend
            fvz(ij, k) = wadv(ij, k, 1)
            fvzu(ij, k) = wadv(ij, k, 2) &
     &                  + wadv(ij, k, 3) &
     &                  + wadv(ij, k, 4) &
     &                  + wadv(ij, k, 5) &
     &                  + wadv(ij, k, 6) &
     &                  + wadv(ij, k, 7) &
     &                  + wadv(ij, k, 8) &
     &                  + wadv(ij, k, 9)
            fvzd(ij, k) = wadv(ij+ln , k+1, 6) &
     &                  + wadv(ij+lne, k+1, 7) &
     &                  + wadv(ij+le , k+1, 8) &
     &                  + wadv(ij+lse, k+1, 9) &
     &                  + wadv(ij+ls , k+1, 2) &
     &                  + wadv(ij+lsw, k+1, 3) &
     &                  + wadv(ij+lw , k+1, 4) &
     &                  + wadv(ij+lnw, k+1, 5) 
         end do 
      end do
      do ij = ijvstr, ijvend
         fvz(ij, kstr+kz) = wadv(ij, kstr+kz, 1) * hvbot(ij)
         fvzd(ij, kstr) = wadv(ij+ln , kstr+1, 6) &
     &                  + wadv(ij+lne, kstr+1, 7) &
     &                  + wadv(ij+le , kstr+1, 8) &
     &                  + wadv(ij+lse, kstr+1, 9) &
     &                  + wadv(ij+ls , kstr+1, 2) &
     &                  + wadv(ij+lsw, kstr+1, 3) &
     &                  + wadv(ij+lw , kstr+1, 4) &
     &                  + wadv(ij+lnw, kstr+1, 5) 
         fvzd(ij, kstr+kz-1) = (  wadv(ij+ln , kstr+kz, 6) &
     &                          + wadv(ij+lne, kstr+kz, 7) &
     &                          + wadv(ij+le , kstr+kz, 8) &
     &                          + wadv(ij+lse, kstr+kz, 9) &
     &                          + wadv(ij+ls , kstr+kz, 2) &
     &                          + wadv(ij+lsw, kstr+kz, 3) &
     &                          + wadv(ij+lw , kstr+kz, 4) &
     &                          + wadv(ij+lnw, kstr+kz, 5)) * hvbot(ij)
      end do

      do k = kstr, kstr+kz-1
         do ij = ijvstr, ijvend
            div(ij, k) = (  (  (fvx(ij+le, k) - fvx(ij, k)) * rx &
     &                       + (fvy(ij+ln, k) - fvy(ij, k)) * rym(ij)) * &
     &                      rxu(ij) * ryu(ij) &
     &                    + (  fvne(ij+lne, k) - fvne(ij, k) &
     &                       + fvse(ij+le, k) - fvse(ij+ln, k)) * &
     &                      rx * rym(ij) * rxu(ij) * ryu(ij) &
     &                    + (fvz(ij, k) - fvz(ij, k+1)) * rs(k) &
     &                    + (fvzu(ij, k) - fvzd(ij, k)) * rs(k)) * &
     &                   amskv(ij, k)
         end do
      end do

    do ij = ijvstr, ijvend
        fuzd(ij, kstr) = - 0.5d0 *             &
    &   (  wadv(ij+ln , kstr+1, 6) *           &
    &      (uy(ij+ln , kstr+1) + uy(ij, kstr)) &
    &    + wadv(ij+lne, kstr+1, 7) *           &
    &      (uy(ij+lne, kstr+1) + uy(ij, kstr)) &
    &    + wadv(ij+le , kstr+1, 8) *           &
    &      (uy(ij+le , kstr+1) + uy(ij, kstr)) &
    &    + wadv(ij+lse, kstr+1, 9) *           &
    &      (uy(ij+lse, kstr+1) + uy(ij, kstr)) &
    &    + wadv(ij+ls , kstr+1, 2) *           &
    &      (uy(ij+ls , kstr+1) + uy(ij, kstr)) &
    &    + wadv(ij+lsw, kstr+1, 3) *           &
    &      (uy(ij+lsw, kstr+1) + uy(ij, kstr)) &
    &    + wadv(ij+lw , kstr+1, 4) *           &
    &      (uy(ij+lw , kstr+1) + uy(ij, kstr)) &
    &    + wadv(ij+lnw, kstr+1, 5) *           &
    &      (uy(ij+lnw, kstr+1) + uy(ij, kstr)))

        fvzd(ij, kstr) = - 0.5d0 *             &
    &   (  wadv(ij+ln , kstr+1, 6) *           &
    &      (vy(ij+ln , kstr+1) + vy(ij, kstr)) &
    &    + wadv(ij+lne, kstr+1, 7) *           &
    &      (vy(ij+lne, kstr+1) + vy(ij, kstr)) & 
    &    + wadv(ij+le , kstr+1, 8) *           &
    &      (vy(ij+le , kstr+1) + vy(ij, kstr)) &
    &    + wadv(ij+lse, kstr+1, 9) *           &
    &      (vy(ij+lse, kstr+1) + vy(ij, kstr)) &
    &    + wadv(ij+ls , kstr+1, 2) *           &
    &      (vy(ij+ls , kstr+1) + vy(ij, kstr)) &
    &    + wadv(ij+lsw, kstr+1, 3) *           &
    &      (vy(ij+lsw, kstr+1) + vy(ij, kstr)) &
    &    + wadv(ij+lw , kstr+1, 4) *           &
    &      (vy(ij+lw , kstr+1) + vy(ij, kstr)) &
    &    + wadv(ij+lnw, kstr+1, 5) *           &
    &      (vy(ij+lnw, kstr+1) + vy(ij, kstr)))
    enddo 

    do k = kstr+1, kend
       do ij = ijvstr, ijvend
           fuz(ij, k) = &
    &                 - wadv(ij, k, 1) * 0.5d0 * &
    &                   (uy(ij, k-1) + uy(ij, k))
           fvz(ij, k) = &
    &                 - wadv(ij, k, 1) * 0.5d0 * &
    &                   (vy(ij, k-1) + vy(ij, k))

           fuzu(ij, k) = - 0.5d0 * &
    &      (  wadv(ij, k, 2) * (uy(ij, k) + uy(ij+ln , k-1)) &
    &       + wadv(ij, k, 3) * (uy(ij, k) + uy(ij+lne, k-1)) &
    &       + wadv(ij, k, 4) * (uy(ij, k) + uy(ij+le , k-1)) &
    &       + wadv(ij, k, 5) * (uy(ij, k) + uy(ij+lse, k-1)) &
    &       + wadv(ij, k, 6) * (uy(ij, k) + uy(ij+ls , k-1)) &
    &       + wadv(ij, k, 7) * (uy(ij, k) + uy(ij+lsw, k-1)) &
    &       + wadv(ij, k, 8) * (uy(ij, k) + uy(ij+lw , k-1)) &
    &       + wadv(ij, k, 9) * (uy(ij, k) + uy(ij+lnw, k-1)))
           fuzd(ij, k) = - 0.5d0 * &
    &      (  wadv(ij+ln , k+1, 6) * (uy(ij+ln , k+1) + uy(ij, k)) &
    &       + wadv(ij+lne, k+1, 7) * (uy(ij+lne, k+1) + uy(ij, k)) &
    &       + wadv(ij+le , k+1, 8) * (uy(ij+le , k+1) + uy(ij, k)) &
    &       + wadv(ij+lse, k+1, 9) * (uy(ij+lse, k+1) + uy(ij, k)) &
    &       + wadv(ij+ls , k+1, 2) * (uy(ij+ls , k+1) + uy(ij, k)) &
    &       + wadv(ij+lsw, k+1, 3) * (uy(ij+lsw, k+1) + uy(ij, k)) &
    &       + wadv(ij+lw , k+1, 4) * (uy(ij+lw , k+1) + uy(ij, k)) &
    &       + wadv(ij+lnw, k+1, 5) * (uy(ij+lnw, k+1) + uy(ij, k)))

           fvzu(ij, k) = - 0.5d0 * &
    &      (  wadv(ij, k, 2) * (vy(ij, k) + vy(ij+ln , k-1)) &
    &       + wadv(ij, k, 3) * (vy(ij, k) + vy(ij+lne, k-1)) &
    &       + wadv(ij, k, 4) * (vy(ij, k) + vy(ij+le , k-1)) &
    &       + wadv(ij, k, 5) * (vy(ij, k) + vy(ij+lse, k-1)) &
    &       + wadv(ij, k, 6) * (vy(ij, k) + vy(ij+ls , k-1)) &
    &       + wadv(ij, k, 7) * (vy(ij, k) + vy(ij+lsw, k-1)) &
    &       + wadv(ij, k, 8) * (vy(ij, k) + vy(ij+lw , k-1)) &
    &       + wadv(ij, k, 9) * (vy(ij, k) + vy(ij+lnw, k-1)))
           fvzd(ij, k) = - 0.5d0 * &
    &      (  wadv(ij+ln , k+1, 6) * (vy(ij+ln , k+1) + vy(ij, k)) &
    &       + wadv(ij+lne, k+1, 7) * (vy(ij+lne, k+1) + vy(ij, k)) &
    &       + wadv(ij+le , k+1, 8) * (vy(ij+le , k+1) + vy(ij, k)) &
    &       + wadv(ij+lse, k+1, 9) * (vy(ij+lse, k+1) + vy(ij, k)) &
    &       + wadv(ij+ls , k+1, 2) * (vy(ij+ls , k+1) + vy(ij, k)) &
    &       + wadv(ij+lsw, k+1, 3) * (vy(ij+lsw, k+1) + vy(ij, k)) &
    &       + wadv(ij+lw , k+1, 4) * (vy(ij+lw , k+1) + vy(ij, k)) &
    &       + wadv(ij+lnw, k+1, 5) * (vy(ij+lnw, k+1) + vy(ij, k)))
       enddo
    enddo
#ifdef OPT_BBL
    do ij = ijvstr, ijvend
       fuzd(ij,kend-1) = 0.d0
       fvzd(ij,kend-1) = 0.d0
       k = nbotv(ij) - amskvb(ij)
           fuzd(ij, k) = fuzd(ij, k) - 0.5d0 * &
    &      (  wadv(ij+ln , kend, 6) * (uy(ij+ln , kend) + uy(ij, k)) &
    &       + wadv(ij+lne, kend, 7) * (uy(ij+lne, kend) + uy(ij, k)) &
    &       + wadv(ij+le , kend, 8) * (uy(ij+le , kend) + uy(ij, k)) &
    &       + wadv(ij+lse, kend, 9) * (uy(ij+lse, kend) + uy(ij, k)) &
    &       + wadv(ij+ls , kend, 2) * (uy(ij+ls , kend) + uy(ij, k)) &
    &       + wadv(ij+lsw, kend, 3) * (uy(ij+lsw, kend) + uy(ij, k)) &
    &       + wadv(ij+lw , kend, 4) * (uy(ij+lw , kend) + uy(ij, k)) &
    &       + wadv(ij+lnw, kend, 5) * (uy(ij+lnw, kend) + uy(ij, k)))
           fvzd(ij, k) = fvzd(ij, k) - 0.5d0 * &
    &      (  wadv(ij+ln , kend, 6) * (vy(ij+ln , kend) + vy(ij, k)) &
    &       + wadv(ij+lne, kend, 7) * (vy(ij+lne, kend) + vy(ij, k)) &
    &       + wadv(ij+le , kend, 8) * (vy(ij+le , kend) + vy(ij, k)) &
    &       + wadv(ij+lse, kend, 9) * (vy(ij+lse, kend) + vy(ij, k)) &
    &       + wadv(ij+ls , kend, 2) * (vy(ij+ls , kend) + vy(ij, k)) &
    &       + wadv(ij+lsw, kend, 3) * (vy(ij+lsw, kend) + vy(ij, k)) &
    &       + wadv(ij+lw , kend, 4) * (vy(ij+lw , kend) + vy(ij, k)) &
    &       + wadv(ij+lnw, kend, 5) * (vy(ij+lnw, kend) + vy(ij, k)))
    enddo
#endif

    do k = kstr, kend
       do  ij = ijvstr, ijvend+nxdim
           fuy(ij, k) = &
    &                 - (  cye(ij, k) * vadv(ij, k) * hxt(ij)           &
    &                    + cyw(ij+le, k) * vadv(ij+le, k) * hxt(ij+le)) &
    &                   / 6.d0 * (uy(ij, k) + uy(ij+ls, k)) * 0.5d0
           fvy(ij, k) = &
    &                 - (  cye(ij, k) * vadv(ij, k) * hxt(ij)           &
    &                    + cyw(ij+le, k) * vadv(ij+le, k) * hxt(ij+le)) &
    &                   / 6.d0 * (vy(ij, k) + vy(ij+ls, k)) * 0.5d0
       enddo
    enddo

    do k = kstr, kend
       do ij = ijvstr, ijvend+1
         fux(ij, k) = &
    &                 - (  cxn(ij, k) * uadv(ij, k) * hyt(ij)           &
    &                    + cxs(ij+ln, k) * uadv(ij+ln, k) * hyt(ij+ln)) &
    &                   / 6.d0 * (uy(ij, k) + uy(ij+lw, k)) * 0.5d0
         fvx(ij, k) = &
    &                 - (  cxn(ij, k) * uadv(ij, k) * hyt(ij)           &
    &                    + cxs(ij+ln, k) * uadv(ij+ln, k) * hyt(ij+ln)) &
    &                   / 6.d0 * (vy(ij, k) + vy(ij+lw, k)) * 0.5d0
       enddo
    enddo

    do k = kstr, kend
       do ij = ijvstr, ijvend+nxdim+1
           fune(ij, k) = &
    &                  - cne(ij, k) / 6.d0 *                 &
    &                    (  uadv(ij, k) * dy(ij) * hyt(ij)   &
    &                     + vadv(ij, k) * dx * hxt(ij)) *    &
    &                    (uy(ij, k) + uy(ij+lsw, k)) * 0.5d0
           fvne(ij, k) = &
    &                  - cne(ij, k) / 6.d0 *                 &
    &                    (  uadv(ij, k) * dy(ij) * hyt(ij)   &
    &                     + vadv(ij, k) * dx * hxt(ij)) *    &
    &                    (vy(ij, k) + vy(ij+lsw, k)) * 0.5d0
       enddo
    enddo

    do k = kstr, kend
       do ij = ijvstr, ijvend+nxdim
           fuse(ij, k) = &
    &                  - cse(ij, k) / 6.d0 *                 &
    &                    (  uadv(ij, k) * dy(ij) * hyt(ij)   &
    &                     - vadv(ij, k) * dx * hxt(ij)) *    &
    &                    (uy(ij+lw, k) + uy(ij+ls, k)) * 0.5d0
           fvse(ij, k) = &
    &                  - cse(ij, k) / 6.d0 *                 &
    &                    (  uadv(ij, k) * dy(ij) * hyt(ij)   &
    &                     - vadv(ij, k) * dx * hxt(ij)) *    &
    &                    (vy(ij+lw, k) + vy(ij+ls, k)) * 0.5d0
       enddo
    enddo

    do k = kstr, kstr+kz-2
       do ij = ijvstr, ijvend
!            gx(ij, k) = (  gx(ij, k)                                   &
            xx1(ij, k) = ( &
    &                   + (  (fux(ij+le, k) - fux(ij, k)) * rx         &
    &                      + (fuy(ij+ln, k) - fuy(ij, k)) * rym(ij)) * &
    &                     rxu(ij) * ryu(ij)                            &
    &                   + (  fune(ij+lne, k) - fune(ij   , k)          &
    &                      + fuse(ij+le , k) - fuse(ij+ln, k)) *       &
    &                     rx * rym(ij) * rxu(ij) * ryu(ij)             &
    &                   + (fuz(ij, k) - fuz(ij, k+1)) * rs(k)          &
    &                   + (fuzu(ij, k) - fuzd(ij, k)) * rs(k)          &
    &                   + uy(ij, k) * div(ij, k)                       &
    &                   + vy(ij, k) * vy(ij, k) * hyxu(ij)             &
    &                   - uy(ij, k) * vy(ij, k) * hxyu(ij)) *          &
    &                  amskv(ij, k)
!           gy(ij, k) = (  gy(ij, k)                                    &
           yy1(ij, k) = ( &
    &                   + (  (fvx(ij+le, k) - fvx(ij, k)) * rx         &
    &                      + (fvy(ij+ln, k) - fvy(ij, k)) * rym(ij)) * &
    &                     rxu(ij) * ryu(ij)                            &
    &                   + (  fvne(ij+lne, k) - fvne(ij   , k)          &
    &                      + fvse(ij+le , k) - fvse(ij+ln, k)) *       &
    &                     rx * rym(ij) * rxu(ij) * ryu(ij)             &
    &                   + (fvz(ij, k) - fvz(ij, k+1)) * rs(k)          &
    &                   + (fvzu(ij, k) - fvzd(ij, k)) * rs(k)          &
    &                   + vy(ij, k) * div(ij, k)                       &
    &                   + uy(ij, k) * uy(ij, k) * hxyu(ij)             &
    &                   - uy(ij, k) * vy(ij, k) * hyxu(ij)) *          &
    &                  amskv(ij, k)
       enddo
    enddo

    k = kstr+kz-1
    do ij = ijvstr, ijvend
!           gx(ij, k) = (  gx(ij, k)                                   &
            xx1(ij, k) = ( &
    &                  + (  (fux(ij+le, k) - fux(ij, k)) * rx         &
    &                     + (fuy(ij+ln, k) - fuy(ij, k)) * rym(ij)) * &
    &                    rxu(ij) * ryu(ij)                            &
    &                  + (  fune(ij+lne, k) - fune(ij   , k)          &
    &                     + fuse(ij+le , k) - fuse(ij+ln, k)) *       &
    &                    rx * rym(ij) * rxu(ij) * ryu(ij)             &
    &                  + (  fuz(ij, k)                                &
    &                     - fuz(ij, k+1) * hvbot(ij)) * rs(k)         &
    &                  + (  fuzu(ij, k)                               &
    &                     - fuzd(ij, k) * hvbot(ij)) * rs(k)          &
    &                  + uy(ij, k) * div(ij, k)                       &
    &                  + vy(ij, k) * vy(ij, k) * hyxu(ij)             &
    &                  - uy(ij, k) * vy(ij, k) * hxyu(ij)) *          &
    &                 amskv(ij, k)
!           gy(ij, k) = (  gy(ij, k)                                   &
            yy1(ij, k) = ( &
    &                  + (  (fvx(ij+le, k) - fvx(ij, k)) * rx         &
    &                     + (fvy(ij+ln, k) - fvy(ij, k)) * rym(ij)) * &
    &                    rxu(ij) * ryu(ij)                            &
    &                  + (  fvne(ij+lne, k) - fvne(ij   , k)          &
    &                     + fvse(ij+le , k) - fvse(ij+ln, k)) *       &
    &                    rx * rym(ij) * rxu(ij) * ryu(ij)             &
    &                  + (  fvz(ij, k)                                &
    &                     - fvz(ij, k+1) * hvbot(ij)) * rs(k)         &
    &                  + (  fvzu(ij, k)                               &
    &                     - fvzd(ij, k) * hvbot(ij)) * rs(k)          &
    &                  + vy(ij, k) * div(ij, k)                       &
    &                  + uy(ij, k) * uy(ij, k) * hxyu(ij)             &
    &                  - uy(ij, k) * vy(ij, k) * hyxu(ij)) *          &
    &                 amskv(ij, k)
    enddo

    do k = kstr+kz, kend
       do ij = ijvstr, ijvend
!            gx(ij, k) = (  gx(ij, k)                                &
            xx1(ij, k) = ( &            
    &                   + (  (  (  fux(ij+le, k)                    &
    &                            - fux(ij, k)) * rx                 &
    &                         + (  fuy(ij+ln, k)                    &
    &                            - fuy(ij, k)) * rym(ij)) *         &
    &                        rxu(ij) * ryu(ij)                      &
    &                      + (  fune(ij+lne, k) - fune(ij   , k)    &
    &                         + fuse(ij+le , k) - fuse(ij+ln, k)) * &
    &                        rx * rym(ij) * rxu(ij) * ryu(ij)       &
    &                      + fuz(ij, k) - fuz(ij, k+1)              &
    &                      + fuzu(ij, k) - fuzd(ij, k)) * rz(ij, k) &
    &                   + vy(ij, k) * vy(ij, k) * hyxu(ij)          &
    &                   - uy(ij, k) * vy(ij, k) * hxyu(ij)) *       &
    &                  amskv(ij, k)
!           gy(ij, k) = (  gy(ij, k)                                 &
            yy1(ij, k) = ( &
    &                   + (   (  (  fvx(ij+le, k)                   &
    &                             - fvx(ij, k)) * rx                &
    &                          + (  fvy(ij+ln, k)                   &
    &                             - fvy(ij, k)) * rym(ij)) *        &
    &                         rxu(ij) * ryu(ij)                     &
    &                      + (  fvne(ij+lne, k) - fvne(ij   , k)    &
    &                         + fvse(ij+le , k) - fvse(ij+ln, k)) * &
    &                        rx * rym(ij) * rxu(ij) * ryu(ij)       &
    &                      + fvz(ij, k) - fvz(ij, k+1)              &
    &                      + fvzu(ij, k) - fvzd(ij, k)) * rz(ij, k) &
    &                   + uy(ij, k) * uy(ij, k) * hxyu(ij)          &
    &                   - uy(ij, k) * vy(ij, k) * hyxu(ij)) *       &
    &                  amskv(ij, k)
       enddo
    enddo

! ---- Adams-Bashforth scheme
      ncall=ncall +1
      if( (ncall .ge. 3) .or. (.not. oeof)) then
        ncall=3
        do k = kstr, kend
        do ij = ijvstr, ijvend
            gx(ij, k) =  gx(ij, k) &
     & +  c1 *xx1(ij, k)  + c2* xx2(ij, k) + c3*xx3(ij, k) 

            gy(ij, k) =  gy(ij, k) &
     & +  c1 *yy1(ij, k)  + c2* yy2(ij, k) + c3*yy3(ij, k) 
        end do
        end do
      else 
         if(ncall .eq. 1) then ! forward
            do k = kstr, kend
            do ij = ijvstr, ijvend
            gx(ij, k) =  gx(ij, k)  +  xx1(ij, k)
            gy(ij, k) =  gy(ij, k)  +  yy1(ij, k)
            end do
            end do
         else if(ncall .eq. 2) then !2nd order ab
            do k = kstr, kend
            do ij = ijvstr, ijvend
               gx(ij, k) =  gx(ij, k) & 
     &       +  1.5d0 *xx1(ij, k)  - 0.5d0* xx2(ij, k)

               gy(ij, k) =  gy(ij, k) & 
     &       +  1.5d0 *yy1(ij, k)  - 0.5d0* yy2(ij, k)
            end do
            end do
         end if
      end if


      do k = kstr, kend
      do ij = ijvstr, ijvend
          xx3(ij, k)=xx2(ij, k) !n-1 > n-2
          yy3(ij, k)=yy2(ij, k)

          xx2(ij, k)=xx1(ij, k) !n> n-1
          yy2(ij, k)=yy1(ij, k)
      end do
      end do

      return

   end subroutine advvel
#ifdef OPT_BBL
! *********************************************************************
  subroutine advvlb(                                                  &
         &            fux,    fuy,   fune,   fuse,                    &
         &            fvx,    fvy,   fvne,   fvse,                    &
         &             gx,     gy,     xx,     yy,                    &
         &             uy,     ux,     vy,     vx,                    &
         &           uadv,   vadv,   wadv)
    use zocdim,  only :                                               &
         &  nxydim,  nxdim,  nydim,  nzdim,                           &
         &    kstr,   kend,     kz,                                   &
         &  ijvstr, ijvend,                                           &
         &      le,     lw,     ln,     ls,                           &
         &     lsw,    lnw,    lse,    lne,                           &
         &   oinit,  ofinal
    use zocgrd,  only :                                               &
         &      dx,     dy,                                           &  
         &     dzv,    dzm,                                           &  
         &      rx,    rym,                                           &
         &     hxt,    hyt,                                           &
         &     hxu,    hyu,   hxyu,   hyxu,                           &
         &     rxu,    ryu   
    use zocmsk,  only :  amskvb,  nbotv

    use brstt

    implicit none

    real(8), parameter :: c1=23.d0/12.d0, c2=-16.d0/12.d0, c3=5.d0/12.d0
    real(8), save ::      xx1(nxydim),     yy1(nxydim)
    real(8), save ::      xx2(nxydim),     yy2(nxydim)
    real(8), save ::      xx3(nxydim),     yy3(nxydim)
    integer, save :: ncall = 0
    logical, save :: oeof

    real(8),intent(in) ::      uy(nxydim, nzdim),     ux(nxydim, nzdim)
    real(8),intent(in) ::      vy(nxydim, nzdim),     vx(nxydim, nzdim)
    real(8),intent(in) ::    uadv(nxydim, nzdim),   vadv(nxydim, nzdim)
    real(8),intent(in) ::    wadv(nxydim, nzdim, 9)
    real(8),intent(inout) ::      gx(nxydim, nzdim),     gy(nxydim, nzdim)
    real(8),intent(inout) ::      xx(nxydim, nzdim),     yy(nxydim, nzdim)
    real(8),intent(out) ::     fux(nxydim, nzdim),    fvx(nxydim, nzdim)
    real(8),intent(out) ::     fuy(nxydim, nzdim),    fvy(nxydim, nzdim)
    real(8),intent(out) ::    fune(nxydim, nzdim),   fvne(nxydim, nzdim)
    real(8),intent(out) ::    fuse(nxydim, nzdim),   fvse(nxydim, nzdim)

    real(8),save ::     cxn(nxydim),    cxs(nxydim)
    real(8),save ::     cye(nxydim),    cyw(nxydim)
    real(8),save ::     cne(nxydim),    cse(nxydim)
    real(8),save ::      rz(nxydim),    rzm(nxydim)

    real(8) ::     fuz(nxydim, nzdim),    fvz(nxydim, nzdim)
    real(8) ::    fuzu(nxydim, nzdim),   fvzu(nxydim, nzdim)
    real(8) ::    fuzd(nxydim, nzdim),   fvzd(nxydim, nzdim)

    integer ::    ij,      k,    kup

    logical, save :: ofirst = .true.

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
     &                                        -1.d0,  -1,  -1)
      call rstadd(xx3, oeof, nxdim, nydim, 1, 'XX3', 'SFC', &
     &                                        -1.d0,  -1,  -1)
      call rstadd(yy2, oeof, nxdim, nydim, 1, 'YY2', 'SFC', &
     &                                        -1.d0,  -1,  -1)
      call rstadd(yy3, oeof, nxdim, nydim, 1, 'YY3', 'SFC', &
     &                                        -1.d0,  -1,  -1)
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

    do ij = 1, nxydim
         rz   (ij) = 1.d0 / dzv(ij, kend)
         rzm  (ij) = 1.d0 / dzm(ij, kend)
    end do

    do  ij = 1, nxydim
         cxn(ij) = 0.d0
         cxs(ij) = 0.d0
         cye(ij) = 0.d0
         cyw(ij) = 0.d0
         cne(ij) = 0.d0
         cse(ij) = 0.d0
    enddo
    do ij = ijvstr-nxdim-1, ijvend+nxdim+1
        cxn(ij) = amskvb(ij) * amskvb(ij+lw) *           &
    &             (3.d0 - amskvb(ij+ls) - amskvb(ij+lsw) &
    &                   + amskvb(ij+ls) * amskvb(ij+lsw))
        cxs(ij) = amskvb(ij+ls) * amskvb(ij+lsw) *   &
    &             (3.d0 - amskvb(ij) - amskvb(ij+lw) &
    &                   + amskvb(ij) * amskvb(ij+lw))
        cye(ij) = amskvb(ij) * amskvb(ij+ls) *           &
    &             (3.d0 - amskvb(ij+lw) - amskvb(ij+lsw) &
    &                   + amskvb(ij+lw) * amskvb(ij+lsw))
        cyw(ij) = amskvb(ij+lw) * amskvb(ij+lsw) *   &
    &             (3.d0 - amskvb(ij) - amskvb(ij+ls) &
    &                   + amskvb(ij) * amskvb(ij+ls))
        cne(ij) = amskvb(ij) * amskvb(ij+lsw) *  &
    &             (3.d0 - amskvb(ij+lw) - amskvb(ij+ls)) 
        cse(ij) = amskvb(ij+lw) * amskvb(ij+ls) * &
    &             (3.d0 - amskvb(ij) - amskvb(ij+lsw))
    enddo
    end if
    do ij = ijvstr, ijvend
       kup = max(nbotv(ij)-1, 1)
       fuz(ij, kend) = &
    &                 - wadv(ij, kend, 1) * 0.5d0 * &
    &                   (uy(ij, kup) + uy(ij, kend))
       fvz(ij, kend) = &
    &                 - wadv(ij, kend, 1) * 0.5d0 * &
    &                   (vy(ij, kup) + vy(ij, kend))

       k = nbotv(ij+ln ) - amskvb(ij+ln )
       fuz(ij, kend) = fuz(ij,kend) - 0.5d0 * &
    &      ( wadv(ij, kend, 2) * (uy(ij, kend) + uy(ij+ln , k)))
       k = nbotv(ij+lne) - amskvb(ij+lne)
       fuz(ij, kend) = fuz(ij,kend) - 0.5d0 * &
    &      ( wadv(ij, kend, 3) * (uy(ij, kend) + uy(ij+lne, k)))
       k = nbotv(ij+le ) - amskvb(ij+le )
       fuz(ij, kend) = fuz(ij,kend) - 0.5d0 * &
    &      ( wadv(ij, kend, 4) * (uy(ij, kend) + uy(ij+le , k)))
       k = nbotv(ij+lse) - amskvb(ij+lse)
       fuz(ij, kend) = fuz(ij,kend) - 0.5d0 * &
    &      ( wadv(ij, kend, 5) * (uy(ij, kend) + uy(ij+lse, k)))
       k = nbotv(ij+ls ) - amskvb(ij+ls )
       fuz(ij, kend) = fuz(ij,kend) - 0.5d0 * &
    &      ( wadv(ij, kend, 6) * (uy(ij, kend) + uy(ij+ls , k)))
       k = nbotv(ij+lsw) - amskvb(ij+lsw)
       fuz(ij, kend) = fuz(ij,kend) - 0.5d0 * &
    &      ( wadv(ij, kend, 7) * (uy(ij, kend) + uy(ij+lsw, k)))
       k = nbotv(ij+lw ) - amskvb(ij+lw )
       fuz(ij, kend) = fuz(ij,kend) - 0.5d0 * &
    &      ( wadv(ij, kend, 8) * (uy(ij, kend) + uy(ij+lw , k)))
       k = nbotv(ij+lnw ) - amskvb(ij+lnw )
       fuz(ij, kend) = fuz(ij,kend) - 0.5d0 * &
    &      ( wadv(ij, kend, 9) * (uy(ij, kend) + uy(ij+lnw, k)))

       k = nbotv(ij+ln ) - amskvb(ij+ln )
       fvz(ij, kend) = fvz(ij,kend) - 0.5d0 * &
    &      ( wadv(ij, kend, 2) * (vy(ij, kend) + vy(ij+ln , k)))
       k = nbotv(ij+lne) - amskvb(ij+lne)
       fvz(ij, kend) = fvz(ij,kend) - 0.5d0 * &
    &      ( wadv(ij, kend, 3) * (vy(ij, kend) + vy(ij+lne, k)))
       k = nbotv(ij+le ) - amskvb(ij+le )
       fvz(ij, kend) = fvz(ij,kend) - 0.5d0 * &
    &      ( wadv(ij, kend, 4) * (vy(ij, kend) + vy(ij+le , k)))
       k = nbotv(ij+lse) - amskvb(ij+lse)
       fvz(ij, kend) = fvz(ij,kend) - 0.5d0 * &
    &      ( wadv(ij, kend, 5) * (vy(ij, kend) + vy(ij+lse, k)))
       k = nbotv(ij+ls ) - amskvb(ij+ls )
       fvz(ij, kend) = fvz(ij,kend) - 0.5d0 * &
    &      ( wadv(ij, kend, 6) * (vy(ij, kend) + vy(ij+ls , k)))
       k = nbotv(ij+lsw) - amskvb(ij+lsw)
       fvz(ij, kend) = fvz(ij,kend) - 0.5d0 * &
    &      ( wadv(ij, kend, 7) * (vy(ij, kend) + vy(ij+lsw, k)))
       k = nbotv(ij+lw ) - amskvb(ij+lw )
       fvz(ij, kend) = fvz(ij,kend) - 0.5d0 * &
    &      ( wadv(ij, kend, 8) * (vy(ij, kend) + vy(ij+lw , k)))
       k = nbotv(ij+lnw ) - amskvb(ij+lnw )
       fvz(ij, kend) = fvz(ij,kend) - 0.5d0 * &
    &      ( wadv(ij, kend, 9) * (vy(ij, kend) + vy(ij+lnw, k)))
    end do

    do ij = ijvstr, ijvend+nxdim
        fuy(ij, kend) = &
    &                 - (  cye(ij) * vadv(ij, kend) * hxt(ij)           &
    &                    + cyw(ij+le) * vadv(ij+le, kend) * hxt(ij+le)) &
    &                    / 6.d0 *                                       &
    &                   (uy(ij, kend) + uy(ij+ls, kend)) * 0.5d0
        fvy(ij, kend) = &
    &                 - (  cye(ij) * vadv(ij, kend) * hxt(ij)           & 
    &                    + cyw(ij+le) * vadv(ij+le, kend) * hxt(ij+le)) &
    &                   / 6.d0 *                                        &
    &                   (vy(ij, kend) + vy(ij+ls, kend)) * 0.5d0
    end do

    do ij = ijvstr, ijvend+1
        fux(ij, kend) = &
    &                 - (  cxn(ij) * uadv(ij, kend) * hyt(ij)           &
    &                    + cxs(ij+ln) * uadv(ij+ln, kend) * hyt(ij+ln)) &
    &                   / 6.d0 *                                        &
    &                   (uy(ij, kend) + uy(ij+lw, kend)) * 0.5d0
        fvx(ij, kend) = & 
    &                 - (  cxn(ij) * uadv(ij, kend) * hyt(ij)           &
    &                    + cxs(ij+ln) * uadv(ij+ln, kend) * hyt(ij+ln)) &
    &                   / 6.d0 *                                        &
    &                   (vy(ij, kend) + vy(ij+lw, kend)) * 0.5d0
    end do

    do ij = ijvstr, ijvend+nxdim+1
       fune(ij, kend) = &
    &                  - cne(ij) / 6.d0 *                     &
    &                    (  uadv(ij, kend) * dy(ij) * hyt(ij) &
    &                     + vadv(ij, kend) * dx * hxt(ij)) *  &
    &                    (uy(ij, kend) + uy(ij+lsw, kend)) * 0.5d0
       fvne(ij, kend) = &
    &                  - cne(ij) / 6.d0 *                     &
    &                    (  uadv(ij, kend) * dy(ij) * hyt(ij) &
    &                     + vadv(ij, kend) * dx * hxt(ij)) *  &
    &                    (vy(ij, kend) + vy(ij+lsw, kend)) * 0.5d0
    end do

    do ij = ijvstr, ijvend+nxdim
        fuse(ij, kend) = &
    &                  - cse(ij) / 6.d0 *                     &
    &                    (  uadv(ij, kend) * dy(ij) * hyt(ij) &
    &                     - vadv(ij, kend) * dx * hxt(ij)) *  &
    &                    (uy(ij+lw, kend) + uy(ij+ls, kend)) * 0.5d0
        fvse(ij, kend) = &
    &                  - cse(ij) / 6.d0 *                     &
    &                    (  uadv(ij, kend) * dy(ij) * hyt(ij) &
    &                     - vadv(ij, kend) * dx * hxt(ij)) *  &
    &                    (vy(ij+lw, kend) + vy(ij+ls, kend)) * 0.5d0
    end do
    do ij = ijvstr, ijvend
!       gx(ij, kend) = (  gx(ij, kend)                               &
        xx1(ij) = ( &
    &                   + (  (  (  fux(ij+le, kend)                 &
    &                            - fux(ij, kend)) * rx              &
    &                         + (  fuy(ij+ln, kend)                 &
    &                            - fuy(ij, kend)) * rym(ij)) *      &
    &                        rxu(ij) * ryu(ij)                      &
    &                      + (  fune(ij+lne, kend)                  &
    &                         - fune(ij    , kend)                  &
    &                         + fuse(ij+le , kend)                  &
    &                         - fuse(ij+ln , kend)) *               &
    &                        rx * rym(ij) * rxu(ij) * ryu(ij)       &
    &                      + fuz(ij, kend)) * rz(ij)                &
    &                   + vy(ij, kend) * vy(ij, kend) * hyxu(ij)    &
    &                   - uy(ij, kend) * vy(ij, kend) * hxyu(ij)) * &
    &                  amskvb(ij)
!       gy(ij, kend) = (  gy(ij, kend)                               &
        yy1(ij) = ( &
    &                   + (   (  (  fvx(ij+le, kend)                &
    &                             - fvx(ij, kend)) * rx             &
    &                          + (  fvy(ij+ln, kend)                &
    &                             - fvy(ij, kend)) * rym(ij)) *     &
    &                         rxu(ij) * ryu(ij)                     &
    &                      + (  fvne(ij+lne, kend)                  &
    &                         - fvne(ij    , kend)                  &
    &                         + fvse(ij+le , kend)                  &
    &                         - fvse(ij+ln , kend)) *               &
    &                        rx * rym(ij) * rxu(ij) * ryu(ij)       &
    &                      + fvz(ij, kend)) * rz(ij)                &
    &                   + uy(ij, kend) * uy(ij, kend) * hxyu(ij)    &
    &                   - uy(ij, kend) * vy(ij, kend) * hyxu(ij)) * &
    &                  amskvb(ij)
    end do

! --- Adams-Bashforth scheme

      ncall=ncall +1

      if( (ncall .ge. 3) .or. (.not. oeof)) then
      ncall=3
        do ij = ijvstr, ijvend
            gx(ij, kend) =  gx(ij, kend) &
     & +  c1 *xx1(ij)  + c2* xx2(ij) + c3*xx3(ij) 

            gy(ij, kend) =  gy(ij, kend) &
     & +  c1 *yy1(ij)  + c2* yy2(ij) + c3*yy3(ij) 
        end do
      else
        if(ncall .eq. 1) then ! forward
           do ij = ijvstr, ijvend
            gx(ij, kend) =  gx(ij, kend) + xx1(ij)
            gy(ij, kend) =  gy(ij, kend) + yy1(ij)
           end do
        else if(ncall .eq. 2) then !2nd order ab       
          do ij = ijvstr, ijvend
            gx(ij, kend) =  gx(ij, kend) & 
     &    +  1.5d0 *xx1(ij)  -0.5d0* xx2(ij)

            gy(ij, kend) =  gy(ij, kend) &
     &    +  1.5d0 *yy1(ij)  -0.5d0* yy2(ij)
          end do
        end if
      end if




      do ij = ijvstr, ijvend
          xx3(ij)=xx2(ij) !n-1 > n-2
          yy3(ij)=yy2(ij)

          xx2(ij)=xx1(ij) !n> n-1
          yy2(ij)=yy1(ij)
      end do

      return
  end subroutine advvlb
#endif
end module cadvc
