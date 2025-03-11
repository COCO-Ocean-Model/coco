module aprdc

! --- information -----------------------------------------------------
!
!     '02.10.10  H.Hasumi: from MIROC3.1-OMIP
!     '07.04.23  H.Hasumi
!     '08.06.11  H.Hasumi: initial/final processing
!     '08.07.10  H.Hasumi: initial/final processing
!     '09.01.20  T.Suzuki: bug fix (SHIFT UX/VX after STBBVT)
!     '09.10.06  Y.Komuro: bug fix (skip barotro. part when INIT/FINAL)
!     '10.04.14  M.Kurogi: staggered time stepping
!     '10.04.14  M.Kurogi: (COCO4.4 tripolar code by Dr. Suzuki)
!     '12.09.04  Y.Komuro: (change surface water flux by T. Suzuki)
!     '12.10.23  T.Suzuki: for COCO5.0 in F90
!
! ---------------------------------------------------------------------

  use zocdim,  only  :   nxyzdm,   nxydim,   ntdim
  implicit none
  private

  real(8)        ::    tmp(nxyzdm, ntdim)
  real(8)        ::     hz(nxydim)
  real(8)        ::   htmp(nxydim),  ubtmp(nxydim),  vbtmp(nxydim)
  real(8), save  ::  ubtav(nxydim),  vbtav(nxydim)

  real(8)        ::    fux(nxyzdm),   fvx(nxyzdm)
  real(8)        ::    fuy(nxyzdm),   fvy(nxyzdm)
  real(8)        ::   fune(nxyzdm),  fvne(nxyzdm)
  real(8)        ::   fuse(nxyzdm),  fvse(nxyzdm)

  real(8), save  :: ubtav2(nxydim), vbtav2(nxydim),    hav(nxydim)
  real(8)        ::   fact
  integer        ::     nb

  public  ::  predco

contains

  subroutine predco(                                                  &
    &        hx,   ubtx,   vbtx,      w,      r,                      &
    &        ux,     vx,     tx,                                      &
    &       amv,    ahv,                                              &
    &      taux,   tauy,   ptop,                                      & 
    &        ft,  swabs,     fs,   ssfc,                              &
#ifdef OPT_BODY
    &        tq,                                                      &
#endif
    &        uy,     vy,     ty,                                      &
    &        hy,   ubty,   vbty,     ax,                              &
    &        gx,     gy,     xx,     yy,                              &
    &       gxx,    gyy,                                              &
    &      uadv,   vadv,   wadv)

    use zocdim
    use zocgrd
    use cadvc
    use csrvl
    use ctnuv
    use cvisc
    use cvlrd
    use dvdif
    use fshlw
    use fshdf
    use tflxt
    use tovtr
    use tslvt
    use ucloc
    use fbtav
    use dvlta
    use dvlva
    use dwdns
    use bstbc
    use brstt
    use bchmk
    use qckot
    use bshft

    implicit none

    real(8),    intent(inout)  ::     hx(nxydim) 
    real(8),    intent(inout)  ::   ubtx(nxydim),   vbtx(nxydim) 
    real(8),    intent(inout)  ::      w(nxyzdm),      r(nxyzdm)
    real(8),    intent(inout)  ::     ux(nxyzdm),     vx(nxyzdm)
    real(8),    intent(inout)  ::     tx(nxyzdm, ntdim)
    real(8),    intent(inout)  ::    amv(nxyzdm),    ahv(nxyzdm)
    real(8),    intent(in)     ::   taux(nxydim),   tauy(nxydim)
    real(8),    intent(in)     ::   ptop(nxydim)
    real(8),    intent(inout)  ::     ft(nxydim, ntdim)
    real(8),    intent(in)     ::  swabs(nxydim),     fs(nxydim)
    real(8),    intent(in)     ::   ssfc(nxydim)
#ifdef OPT_BODY
    real(8),    intent(in)     ::     tq(nxyzdm, ntdim)
#endif
    real(8),    intent(in)     ::     uy(nxyzdm),     vy(nxyzdm)
    real(8),    intent(inout)  ::     ty(nxyzdm, ntdim)
    real(8),    intent(inout)  ::     hy(nxydim) 
    real(8),    intent(inout)  ::   ubty(nxydim),   vbty(nxydim) 
    real(8),    intent(in)     ::     ax(nxydim, 0:nic)
    real(8),    intent(inout)  ::     gx(nxyzdm),     gy(nxyzdm)
    real(8),    intent(inout)  ::     xx(nxyzdm),     yy(nxyzdm)
    real(8),    intent(inout)  ::    gxx(nxydim),    gyy(nxydim)
    real(8),    intent(inout)  ::   uadv(nxyzdm),   vadv(nxyzdm)
    real(8),    intent(inout)  ::   wadv(nxyzdm, 9)

!---- local variables
    real(8)     ::    hxb(nxydim)
    real(8), save ::   h1(nxydim) 
    integer(4)  ::  itsplt,     ij,    ijk,      n
    logical     ::    oeof

    if (       ( myrank >= ijnode )                                   &
    &    .and. (.not. oinit) .and. (.not. ofinal)) return

    if(oinit) then
    !$acc enter data create(   tmp)
    !$acc enter data create(    hz)
    !$acc enter data create(  htmp,  ubtmp,  vbtmp)
    !$acc enter data create( ubtav,  vbtav)
    !$acc enter data create(   fux,   fvx)
    !$acc enter data create(   fuy,   fvy)
    !$acc enter data create(  fune,  fvne)
    !$acc enter data create(  fuse,  fvse)
    !$acc enter data create( ubtav2, vbtav2,    hav)
    !$acc enter data create( hxb, h1)
       
#ifdef OPT_TRIPOLE
    call rstadd(h1, oeof, nxdim, nydim, 1, 'H1', 'SFC',     &
    &                                        1.d0,  0,  0 )
#else
    call rstadd(h1, oeof, nxdim, nydim, 1, 'H1', 'SFC')
#endif
    !$acc update device(h1)
    if(oeof) then
      !$acc kernels default(present)
      do ij = 1, nxydim
         h1(ij) = hx(ij)
      end do
      !$acc end kernels
    end if  
    end if

    if (ofinal) then
       !$acc update self(h1)
       call finadd(h1, nxdim, nydim, 1, 'H1', 'SFC')
    end if

! *** vertical diffusivity and viscosity ***
    call clcstr('COEFF')
    call vdiff (   amv,    ahv,                                       &
    &                uy,     vy,      r,   taux,   tauy,              &
    &                ty,     hy  )
#ifdef OPT_BBL
    call vdiffb(   amv,    ahv  )
#endif
    call shift_pack_begin
    call shift1(amv, nxdim, nydim, nzdim, 1.d0, -1, -1)
    call shift1(ahv, nxdim, nydim, nzdim, 1.d0,  0,  0)
    call shift_pack_end
    call shift_unpack(amv, 1)
    call shift_unpack(ahv, 2)
    call clcend('COEFF')

! *** baroclinic flow ***
    call clcstr('BRCLI')
#ifdef OPT_BBL
    call rmmskv
#endif
    call srcvel(                                                      &
    &                   gx,     gy,     xx,     yy,                   &
    &                   ux,     vx  )
    call vscvel(                                                      &
    &                   gx,     gy,     xx,     yy,                   &
    &                   uy,     ux,     vy,     vx,                   &
    &                   h1,     hy,                                   &
    &                  amv,   taux,   tauy  )
    call advvel(                                                      &
    &                  fux,    fuy,   fune,   fuse,                   &
    &                  fvx,    fvy,   fvne,   fvse,                   &
    &                   gx,     gy,     xx,     yy,                   &
    &                   uy,     ux,     vy,     vx,                   &
    &                   h1,     hy,                                   &
    &                 uadv,   vadv,   wadv  )
#ifdef OPT_BBL
    call srcvlb(                                                      &
    &                   gx,     gy,     xx,     yy,                   &
    &                   ux,     vx)
    call vscvlb(                                                      &
    &                   gx,     gy,     xx,     yy,                   &
    &                   uy,     ux,     vy,     vx,                   &
    &                  amv)
    call advvlb(                                                      &
    &                  fux,    fuy,   fune,   fuse,                   &
    &                  fvx,    fvy,   fvne,   fvse,                   &
    &                   gx,     gy,     xx,     yy,                   &
    &                   uy,     ux,     vy,     vx,                   &
    &                 uadv,   vadv,   wadv)
#endif
    call tnduvd(                                                      &
    &                  gxx,    gyy,                                   &
    &                   gx,     gy,                                   &
    &                   xx,     yy,                                   &
    &                   uy,     vy,     hy,      r,                   &
    &                   ux,     vx,                                   &
    &                  fux,    fuy,   fune,   fuse,                   &
    &                  fvx,    fvy,   fvne,   fvse)
#ifdef OPT_BBL
    call tnduvb(                                                      &
    &                  gxx,    gyy,                                   &
    &                   gx,     gy,                                   &
    &                   xx,     yy,                                   &
    &                   uy,     vy,     ty,                           &
    &                   ux,     vx,                                   &
    &                  fux,    fuy,   fune,   fuse,                   &
    &                  fvx,    fvy,   fvne,   fvse)
    call admkv1
    call stbbgv(  gx,  gy  )
#endif
    call velrds(                                                      &
    &                   ux,     vx,                                   &
    &                   gx,     gy,    amv,                           &
    &                   h1 )
#ifdef OPT_BBL
    call stbbuv(       ux,     vx  )
#endif
    call clcend('BRCLI')
! *** barotropic flow and surface elevation ***
    call clcstr('BRTRO')
    if (oinit .or. ofinal) then
       call modgxy(   gxx,    gyy,                                    &
    &                ubtx,   vbtx) 
    else
       if (ieuler == 2) then
          !$acc kernels default(present)
          do ij = 1, nxydim
             hz   (ij) = hx   (ij)
             htmp (ij) = hx   (ij)
             ubtmp(ij) = ubtx (ij)
             vbtmp(ij) = vbtx (ij)
             hx   (ij) = hy   (ij)
             ubtx (ij) = ubty (ij)
             vbtx (ij) = vbty (ij)
             hy   (ij) = htmp (ij)
             ubty (ij) = ubtmp(ij)
             vbty (ij) = vbtmp(ij)
          end do
          !$acc end kernels
       else
          !$acc kernels default(present)
          do ij = 1, nxydim
             hz(ij) = hx(ij)
          end do
          !$acc end kernels

          call modgxy(  gxx,    gyy,                                  &
    &                  ubtx,   vbtx  )

          call shift_pack_begin
          call shift2(    gxx, gyy, nxdim, nydim, 1, -1.d0, -1, -1)
          call shift1(ft(:,2),      nxdim, nydim, 1,  1.d0,  0,  0)
          call shift_pack_end
          call shift_unpack(    gxx, 1)
          call shift_unpack(    gyy, 2)
          call shift_unpack(ft(:,2), 3)
          
          nb = ntss * 2
          !$acc kernels default(present)
          ubtav (:) = 0.d0
          vbtav (:) = 0.d0
          ubtav2(:) = 0.d0
          vbtav2(:) = 0.d0
          hav   (:) = hx(:) / dble(nb+1)
          !$acc end kernels
          
          do itsplt = 1, nb
             !$acc kernels default(present)
             htmp (:) = hx  (:)
             ubtmp(:) = ubtx(:)
             vbtmp(:) = vbtx(:)
             !$acc end kernels
             call shalow(                                             &
    &                    htmp,  ubtmp,  vbtmp,                        &
    &                      hx,   ubtx,   vbtx,                        &
    &                     gxx,    gyy,   ptop,  ft(1,2))

             call shift_pack_begin
             call shift2(ubtmp, vbtmp, nxdim, nydim, 1, -1.d0, -1, -1)
             call shift1( htmp,        nxdim, nydim, 1,  1.d0,  0,  0)
             call shift_pack_end
             call shift_unpack(ubtmp, 1)
             call shift_unpack(vbtmp, 2)
             call shift_unpack( htmp, 3)
             
             fact = 2.d0 * dble(nb-itsplt+1) / dble(nb * (nb+1))
             !$acc kernels default(present)
             ubtav (:) = ubtav (:) + ubtmp(:) * fact
             vbtav (:) = vbtav (:) + vbtmp(:) * fact
             ubtav2(:) = ubtav2(:) + ubtmp(:) / dble(nb)
             vbtav2(:) = vbtav2(:) + vbtmp(:) / dble(nb)
             !$acc end kernels
             call shalow(                                             &
    &                      hx,   ubtx,   vbtx,                        &
    &                    htmp,  ubtmp,  vbtmp,                        &
    &                     gxx,    gyy,   ptop,  ft(1,2))


             call shift_pack_begin
             call shift2(ubtx, vbtx, nxdim, nydim, 1, -1.d0, -1, -1)
             call shift1(  hx,       nxdim, nydim, 1,  1.d0,  0,  0)
             call shift_pack_end
             call shift_unpack(ubtx, 1)
             call shift_unpack(vbtx, 2)
             call shift_unpack(  hx, 3)

             !$acc kernels default(present)
             hav(:) = hav(:) + hx(:) / dble(nb+1)
             !$acc end kernels
         end do

         call shift_pack_begin
         call shift2( ubtav,  vbtav, nxdim,  nydim, 1, -1.d0, -1, -1)
         call shift2(ubtav2, vbtav2, nxdim,  nydim, 1, -1.d0, -1, -1)
         call shift_pack_end
         call shift_unpack( ubtav, 1)
         call shift_unpack( vbtav, 2)
         call shift_unpack(ubtav2, 3)
         call shift_unpack(vbtav2, 4)

       end if
       !$acc kernels default(present)
       hx  (:) = hav   (:)
       hxb (:) = hx    (:)
       ubtx(:) = ubtav2(:)
       vbtx(:) = vbtav2(:)
       h1  (:) = 0.5d0 * (hz(:) + hx(:))
       !$acc end kernels
    end if
    call clcend('BRTRO')
! *** velocity for tracer advection ***

    call clcstr('TDIAG')
#ifdef OPT_BBL
    call rmmskv
    call admkvb
#endif
    call veltad(                                                      &
    &                 UADV,   VADV,                                   &
    &                 UBTAV,  VBTAV,     UX,     VX)
!                     UBTAV,  VBTAV,     UY,     VY)

#ifdef OPT_TRIPOLE
    if (.not.(oinit .or. ofinal)) then
#endif
       call shift2(uadv, vadv, nxdim, nydim, nzdim, -1.d0, -1, -1)
#ifdef OPT_TRIPOLE
    end if
#endif
    
#ifdef OPT_BBL
    call rmmskt
    call admktb
#endif
    call wdenst(                                                      &
    &                 wadv,                                           &
    &                 uadv,   vadv,                                   &
    &                   hx,     hz)
    call clcend('TDIAG')
! *** tracer ***
    call clcstr('TRACE')
    if (oinit .or. ofinal) then
       call flxtrc(                                                   &
    &                  tmp,     xx,                                   &
    &                   tx,     hx,                                   &
    &                   ty,     hz,                                   &
    &                 uadv,   vadv,   wadv,                           &
    &                  ahv)
#ifdef OPT_BBL
       call flxtrb(                                                   &
    &                  tmp,     xx,                                   &
    &                   tx,                                           &
    &                   ty,                                           &
    &                 uadv,   vadv,   wadv,                           &
    &                  ahv)
#endif
       call slvtrc(                                                   &
    &                   tx,     hx,                                   &
    &                  tmp,     xx,                                   &
    &                   ft,  swabs,     fs,     hz,   ssfc,           &
    &                   ax)
       call shdiff(     tx,     hx)
       call ovturn(      r,     tx,     hx)
    else
       if (ieuler == 2) then
          !$acc kernels default(present)
          do n = 1, ntdim
             do ijk = 1, nxyzdm
                gx(ijk)    = tx(ijk, n)
                tx(ijk, n) = ty(ijk, n)
                ty(ijk, n) = gx(ijk)
             end do
          end do
          !$acc end kernels
       else
          call flxtrc(                                                &
    &                   tmp,     xx,                                  &
    &                    tx,     hx,                                  &
    &                    ty,     hz,                                  &
    &                  uadv,   vadv,   wadv,                          &
    &                   ahv )
#ifdef OPT_BBL
          call flxtrb(                                                &
    &                  tmp,     xx,                                   &
    &                   tx,                                           &
    &                   ty,                                           &
    &                 uadv,   vadv,   wadv,                           &
    &                  ahv)
          call rmmskt
          call admkt1
          call stbbgt(  tmp,     xx)
#endif
#ifdef OPT_BODY
          !$acc kernels default(present)
          do n = 1, ntdim
             do ijk = 1, nxyzdm
                tmp(ijk, n) = tmp(ijk, n) + tq(ijk, n)
             end do
          end do
          !$acc end kernels
#endif
          call slvtrc(                                                &
    &                   tx,     hx,                                   &
    &                  tmp,     xx,                                   &
    &                   ft,  swabs,     fs,     hz,   ssfc,           &
    &                   ax)

          call shift_pack_begin
          call shift1(tx, nxdim, nydim, nztdim, 1.d0, 0, 0)
          call shift1(hx, nxdim, nydim,      1, 1.d0, 0, 0)
          call shift_pack_end
          call shift_unpack(tx, 1)
          call shift_unpack(hx, 2)
    
          call shdiff(   tx,     hx)
          call clcstr('TUNDIF')
          call tundif(   tx,     hx)
          call clcend('TUNDIF')
          call ovturn(    r,     tx,     hx)
#ifdef OPT_BBL
          call stbbtr(   tx   )
#endif

          call shift_pack_begin
          call shift1( r, nxdim, nydim,  nzdim,  1.d0, 0, 0)
          call shift1(tx, nxdim, nydim, nztdim,  1.d0, 0, 0)
          call shift1(hx, nxdim, nydim,      1,  1.d0, 0, 0)
          call shift_pack_end
          call shift_unpack( r, 1)
          call shift_unpack(tx, 2)
          call shift_unpack(hx, 3)

          call stbctr(    tx,      r)
       end if
    end if
    call clcend('TRACE')
! *** velocity for momentum advection ***

    call clcstr('VDIAG')
    if (.not. (oinit .or. ofinal)) then
       !$acc kernels default(present)
       do ijk = 1, nxyzdm
          uadv(ijk) = ux(ijk)
          vadv(ijk) = vx(ijk)
       end do
       !$acc end kernels
    end if
    call veltad(                                                      &
    &               ux,     vx,                                       &
    &             ubtav, vbtav,   uadv,   vadv)
!    &             ubtx,   vbtx,   uadv,   vadv)

    call shift2(ux, vx, nxdim, nydim, nzdim, -1.d0, -1, -1)

#ifdef OPT_BBL
    call rmmskt
    call admktb
#endif
    call wdenst(                                                      &
    &              w,                                                 &
    &             ux,     vx,                                         &
    &            hxb,     hz)
!    &             hx,     hz)

#ifdef OPT_BBL
    call admskv
#endif
    call velvad(                                                      &
    &             uadv,   vadv,   wadv,                               &
    &               ux,     vx,      w ) 
#ifdef OPT_BBL
    call velvab( &
    &             wadv,                                               &
    &                w) 
    call stbbvt(    ux,     vx    )

    call shift2(ux, vx, nxdim, nydim, nzdim, -1.d0, -1, -1)

    call rmmskv
    call admkv1
#endif

    call shift_pack_begin
    call shift1(uadv,      nxdim, nydim, nzdim, -1.d0,  0,  0)
    call shift1(vadv,      nxdim, nydim, nzdim, -1.d0,  0,  0)
    call shift1(wadv(:,1), nxdim, nydim, nzdim,  1.d0, -1, -1)
    call shift1(wadv(:,2), nxdim, nydim, nzdim,  1.d0, -1, -1)
    call shift1(wadv(:,3), nxdim, nydim, nzdim,  1.d0, -1, -1)
    call shift1(wadv(:,4), nxdim, nydim, nzdim,  1.d0, -1, -1)
    call shift1(wadv(:,5), nxdim, nydim, nzdim,  1.d0, -1, -1)
    call shift1(wadv(:,6), nxdim, nydim, nzdim,  1.d0, -1, -1)
    call shift1(wadv(:,7), nxdim, nydim, nzdim,  1.d0, -1, -1)
    call shift1(wadv(:,8), nxdim, nydim, nzdim,  1.d0, -1, -1)
    call shift1(wadv(:,9), nxdim, nydim, nzdim,  1.d0, -1, -1)
    call shift_pack_end
    call shift_unpack(     uadv,  1)
    call shift_unpack(     vadv,  2)
    call shift_unpack(wadv(:,1),  3)
    call shift_unpack(wadv(:,2),  4)
    call shift_unpack(wadv(:,3),  5)
    call shift_unpack(wadv(:,4),  6)
    call shift_unpack(wadv(:,5),  7)
    call shift_unpack(wadv(:,6),  8)
    call shift_unpack(wadv(:,7),  9)
    call shift_unpack(wadv(:,8), 10)
    call shift_unpack(wadv(:,9), 11)

#ifdef OPT_TRIPOLE
    call excngw(  wadv  )
#endif
    call clcend('VDIAG')

      if ( (.not. oinit) .and. (.not. ofinal)) then
         call chekin(  ubtav,  'UBTAV',                                 &
              &      'ocean zonal transport', 'cm^2/s', &    
              &           nx,     ny,      1, nxydim, 'OCSFCV')
         call chekin(  vbtav,  'VBTAV',                                 &
              & 'ocean meridional transport', 'cm^2/s', &    
              &           nx,     ny,      1, nxydim, 'OCSFCV')
      end if
  end subroutine predco

end module aprdc

