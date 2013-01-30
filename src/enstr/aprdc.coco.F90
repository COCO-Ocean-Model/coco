module aprdc

! --- information -----------------------------------------------------
!
!  HISTORY
!     '00.03.17  H.Hasumi: for Ishizaki scheme
!     '00.05.30  H.Hasumi: parallelized COCO3
!     '00.06.15  H.Hasumi: raise efficiency for parallelization
!     '00.12.07  H.Hasumi: combine par and non-par routines
!     '00.12.14  H.Hasumi: for consistency with the surface height
!                          diffusion
!     '01.01.23  H.Hasumi: for hybrid vertical coordinate
!                          zonal filter
!     '01.02.04  H.Hasumi: combine body routines
!     '01.02.20  H.Hasumi: for incorporating BBL model
!     '01.05.08  H.Hasumi
!     '01.12.03  H.Hasumi: subroutine FLXVEL is split into two
!     '02.05.29  H.Nakano: tracer dimension
!     '07.04.23  H.Hasumi
!     '08.06.11  H.Hasumi: initial/final processing
!     '08.07.10  H.Hasumi: initial/final processing
!     '09.01.20  T.Suzuki: bug fix (SHIFT UX/VX after STBBVT)
!     '10.04.14  M.Kurogi: (COCO4.4 tripolar code by Dr. Suzuki)
!     '12.09.04  Y.Komuro: (change surface water flux by T. Suzuki)
!     '13.01.07  T.Suzuki: for COCO5.0 in F90
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

  public  ::  predco

contains

  subroutine predco(                                                  &
    &        hx,   ubtx,   vbtx,      w,      r,                      &
    &        ux,     vx,     tx,                                      &
    &       amv,    ahv,                                              &
    &      taux,   tauy,     ft,   ptop,                              & 
#ifdef OPT_BODY
    &        tq,                                                      &
#endif
    &        uy,     vy,     ty,                                      &
    &        hy,   ubty,   vbty,                                      &
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
    real(8),    intent(in)     ::     ft(nxydim, ntdim)
#ifdef OPT_BODY
    real(8),    intent(in)     ::     tq(nxyzdm, ntdim)
#endif
    real(8),    intent(in)     ::     uy(nxyzdm),     vy(nxyzdm)
    real(8),    intent(inout)  ::     ty(nxyzdm, ntdim)
    real(8),    intent(inout)  ::     hy(nxydim) 
    real(8),    intent(inout)  ::   ubty(nxydim),   vbty(nxydim) 
    real(8),    intent(inout)  ::     gx(nxyzdm),     gy(nxyzdm)
    real(8),    intent(inout)  ::     xx(nxyzdm),     yy(nxyzdm)
    real(8),    intent(inout)  ::    gxx(nxydim),    gyy(nxydim)
    real(8),    intent(inout)  ::   uadv(nxyzdm),   vadv(nxyzdm)
    real(8),    intent(inout)  ::   wadv(nxyzdm, 9)

!---- local variables
    integer(4)  ::  itsplt,     ij,    ijk,      n
    logical     ::    oeof

    if (       ( myrank >= ijnode )                                   &
    &    .and. (.not. oinit) .and. (.not. ofinal)) return

! *** vertical diffusivity and viscosity ***
    
    call clcstr('COEFF')
    call vdiff (   amv,    ahv,                                       &
    &                uy,     vy,      r,   taux,   tauy,              &
    &                ty,     hy  )
#ifdef OPT_BBL
    call vdiffb(   amv,    ahv  )
#endif
#ifdef OPT_TRIPOLE
    call shift1(    amv,                                              &
    &             nxdim,  nydim,  nzdim,                              &
    &              1.D0,     -1,     -1   )
    call shift1(    ahv,                                              &
    &             nxdim,  nydim,  nzdim,                              &
    &              1.d0,      0,      0   )
#else
    call shift2(    amv,    ahv,                                      &
    &             nxdim,  nydim,  nzdim)
#endif
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
    &                   hx,     hy,                                   &
    &                  amv,   taux,   tauy  )
    call advvel(                                                      &
    &                  fux,    fuy,   fune,   fuse,                   &
    &                  fvx,    fvy,   fvne,   fvse,                   &
    &                   gx,     gy,     xx,     yy,                   &
    &                   uy,     ux,     vy,     vx,                   &
    &                   hx,     hy,                                   &
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
    &                   hx )
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
       else
          do ij = 1, nxydim
             hz(ij) = hx(ij)
          end do
          call modgxy(  gxx,    gyy,                                  &
    &                  ubtx,   vbtx  )
#ifdef OPT_TRIPOLE
          call shift2(   gxx,    gyy,                                 &
    &                  nxdim,  nydim,      1,                         &
    &                  -1.D0,     -1,     -1 )
#else
          call shift2(                                                &
    &                  gxx,    gyy,                                   &
    &                nxdim,  nydim,      1)
#endif
          call btavst( ubtav,  vbtav )
          do itsplt = 1, ntss

             htmp (1:nxydim) = hx  (1:nxydim)
             ubtmp(1:nxydim) = ubtx(1:nxydim)
             vbtmp(1:nxydim) = vbtx(1:nxydim)

             call shalow(                                             &
    &                    htmp,  ubtmp,  vbtmp,                        &
    &                      hx,   ubtx,   vbtx,                        &
    &                     gxx,    gyy,   ptop,  ft(1,2))
             call btavad(                                             &
    &                   ubtav,  vbtav,                                &
    &                   ubtmp,  vbtmp)
             call shalow(                                             &
    &                      hx,   ubtx,   vbtx,                        &
    &                    htmp,  ubtmp,  vbtmp,                        &
    &                     gxx,    gyy,   ptop,  ft(1,2))
#ifdef OPT_TRIPOLE
             call shift2(  ubtx,   vbtx,                              &
    &                     nxdim,  nydim,      1,                      &
    &                     -1.d0,     -1,     -1 )
             call shift1(    hx,                                      &
    &                     nxdim,  nydim,      1,                      &
    &                      1.d0,      0,      0 )
#else
             call shift3(                                             &
    &                       hx,   ubtx,   vbtx,                       &
    &                    nxdim,  nydim,      1)
#endif

          end do
#ifdef OPT_TRIPOLE
          call shift2( ubtav,  vbtav,                                 &
    &                  nxdim,  nydim,      1,                         &
    &                  -1.d0,     -1,     -1 )
#else
          call shift2(                                                &
    &                 ubtav,  vbtav,                                  &
    &                 nxdim,  nydim,      1)
#endif
       end if
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
                     UBTAV,  VBTAV,     UY,     VY)
#ifdef OPT_TRIPOLE
    if (.not.(oinit .or. ofinal)) then
       call shift2(   uadv,    vadv,                                  &
    &                nxdim,   nydim,  nzdim,                          &
    &                -1.d0,      -1,     -1 )
    end if
#else
    call shift2(                                                      &
    &                uadv,   vadv,                                    &
    &               nxdim,  nydim,  nzdim)
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
    &                   ft,     hz )
       call shdiff(     tx,     hx)
       call ovturn(      r,     tx,     hx)
    else
       if (ieuler == 2) then
          do n = 1, ntdim
             do ijk = 1, nxyzdm
                gx(ijk)    = tx(ijk, n)
                tx(ijk, n) = ty(ijk, n)
                ty(ijk, n) = gx(ijk)
             end do
          end do
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
          do n = 1, ntdim
             do ijk = 1, nxyzdm
                tmp(ijk, n) = tmp(ijk, n) + tq(ijk, n)
             end do
          end do
#endif
          call slvtrc(                                                &
    &                   tx,     hx,                                   &
    &                  tmp,     xx,                                   &
    &                   ft,     hz)
#ifdef OPT_TRIPOLE
          call shift1(    tx,                                         &
    &                  nxdim,   nydim, nztdim,                        &
    &                   1.d0,       0,      0 )
          call shift1(    hx,                                         &
    &                  nxdim,   nydim,      1,                        &
    &                   1.d0,       0,      0 )
#else
          call shift1(                                                &
    &                   tx,                                           &
    &                nxdim,  nydim, nztdim)
          call shift1(                                                &
    &                   hx,                                           &
    &                nxdim,  nydim,      1)
#endif
          call shdiff(   tx,     hx)
!          call clcstr('TUNDIF')
!          call tundif(   tx,     hx)
!          call clcend('TUNDIF')
          call ovturn(    r,     tx,     hx)
#ifdef OPT_BBL
          call stbbtr(   tx   )
#endif
#ifdef OPT_TRIPOLE
          call shift1(   r,                                           &
    &                nxdim,   nydim,  nzdim,                          &
    &                 1.d0,       0,      0 )
          call shift1(  tx,                                           &
    &                nxdim,   nydim, nztdim,                          &
    &                 1.d0,       0,      0 )
          call shift1(  hx,                                           &
    &                nxdim,   nydim,      1,                          &
    &                 1.d0,       0,      0 )
#else
          call shift1(                                                &
    &                     r,                                          &
    &                 nxdim,  nydim,  nzdim)
          call shift1(                                                &
    &                    tx,                                          &
    &                 nxdim,  nydim, nztdim)
          call shift1(                                                &
    &                    hx,                                          &
    &                 nxdim,  nydim,      1)
#endif
          call stbctr(    tx,      r)
       end if
    end if
    call clcend('TRACE')
    
! *** velocity for momentum advection ***

    call clcstr('VDIAG')
    if (.not. (oinit .or. ofinal)) then
       do ijk = 1, nxyzdm
          uadv(ijk) = ux(ijk)
          vadv(ijk) = vx(ijk)
       end do
    end if
    call veltad(                                                      &
    &               ux,     vx,                                       &
    &             ubtx,   vbtx,   uadv,   vadv)
#ifdef OPT_TRIPOLE
    call shift2(    ux,     vx,                                       &
    &            nxdim,  nydim,  nzdim,                               &
    &            -1.d0,     -1,     -1 )
#else
    call shift2(                                                      &
    &               ux,     vx,                                       &
    &            nxdim,  nydim,  nzdim)
#endif
#ifdef OPT_BBL
    call rmmskt
    call admktb
#endif
    call wdenst(                                                      &
    &              w,                                                 &
    &             ux,     vx,                                         &
    &             hx,     hz)
#ifdef OPT_BBL
    call admskv
#endif
    call velvad(                                                      &
    &             uadv,   vadv,   wadv,                               &
    &               ux,     vx,      w) 
#ifdef OPT_BBL
    call velvab( &
    &             wadv,                                               &
    &                w) 
    call stbbvt(    ux,     vx    )
#ifdef OPT_TRIPOLE
    call shift2(    ux,     vx,                                       &
    &            nxdim,  nydim,  nzdim,                               &
    &            -1.d0,     -1,     -1 )
#else
    call shift2(                                                      &
    &               ux,     vx,                                       &
    &            nxdim,  nydim,  nzdim)
#endif
    call rmmskv
    call admkv1
#endif
#ifdef OPT_TRIPOLE
    call shift2(                                                      &
    &             uadv,   vadv,                                       &
    &            nxdim,  nydim,  nzdim,                               &
    &            -1.d0,      0,      0)
    call shift3(                                                      &
    &            wadv(1,1), wadv(1,2), wadv(1,3),                     &
    &            nxdim,  nydim,  nzdim,                               &
    &             1.d0,     -1,     -1)
    call shift3(                                                      &
    &            wadv(1,4), wadv(1,5), wadv(1,6),                     &
    &            nxdim,  nydim,  nzdim,                               &
    &             1.d0,     -1,     -1)
    call shift3(                                                      &
    &            wadv(1,7), wadv(1,8), wadv(1,9),                     &
    &            nxdim,  nydim,  nzdim,                               &
    &             1.d0,     -1,     -1)
    call excngw(  wadv  )
#else
    call shift2(                                                      &
    &             uadv,   vadv,                                       &
    &            nxdim,  nydim,  nzdim)
    call shift3(                                                      &
    &            wadv(1,1), wadv(1,2), wadv(1,3),                     &
    &            nxdim,  nydim,  nzdim)
    call shift3(                                                      &
    &            wadv(1,4), wadv(1,5), wadv(1,6),                     &
    &            nxdim,  nydim,  nzdim)
    call shift3(                                                      &
    &            wadv(1,7), wadv(1,8), wadv(1,9),                     &
    &            nxdim,  nydim,  nzdim)
#endif
    call clcend('VDIAG')

  end subroutine predco

end module aprdc

