module aprdc

! --- information -----------------------------------------------------
!
!  HISTORY
!     '00.05.30  H.Hasumi: parallelized COCO3
!     '00.06.15  H.Hasumi: raise efficiency for parallelization
!     '00.12.07  H.Hasumi: combine par and non-par routines
!     '00.12.12  H.Hasumi: for hybrid vertical coordinate
!     '01.01.19  H.Hasumi: add sea surface pressure as a BC variable
!     '01.01.23  H.Hasumi: zonal filter
!     '01.02.02  H.Hasumi: divide the tracer routine
!     '01.02.04  H.Hasumi: combine body routines
!     '01.02.08  H.Hasumi: incorporate BBL model
!     '01.12.03  H.Hasumi: subroutine VELADV is split into two
!     '02.05.29  H.Nakano: tracer dimension
!     '02.06.02  H.Hasumi: trivial modifications to the above
!     '07.04.23  H.Hasumi:
!     '08.06.11  H.Hasumi: initial/final processing
!     '08.07.10  H.Hasumi: initial/final processing
!     '08.11.17  Y.Komuro: bug fix (skip barotro. part when INIT/FINAL)
!     '09.01.20  T.Suzuki: bug fix (SHIFT UX/VX after STBBVT)
!     '10.04.14  M.Kurogi: staggered time stepping
!     '10.04.14  M.Kurogi: (COCO4.4 tripolar code by Dr. Suzuki)
!     '12.01.30  Y.Komuro: (change surface water flux by T. Suzuki)
!     '12.11.23  H.Tatebe: for COCO5.0 in F90
! ---------------------------------------------------------------------

  implicit none

  private
  public  :: predco

contains
  
  subroutine predco(      hx,   ubtx,   vbtx,      w,      r,         &
     &                    ux,     vx,     tx,                         &
     &                   amv,    ahv,                                 &
     &                  taux,   tauy,     ft,   ptop,                 &
#ifdef OPT_BODY
     &                    tq,                                         &
#endif 
     &                    uy,     vy,     ty,                         &
     &                    hy,   ubty,   vbty,                         &
     &                    gx,     gy,     xx,     yy,                 &
     &                   gxx,    gyy,                                 &
     &                  uadv,   vadv,   wadv )

    use zocdim,  only  :                                              &
            nxdim,  nydim,  ntdim,  nxydim,  nxyzdm,  nztdim,         &
               nx,     ny,     nz,                                    &
            oinit, ofinal
    use zocnod,  only  :   myrank,  ijnode
    use zocgrd
    use cadvc
    use csrvl
    use ctnuv
    use cvisc
    use cvlrd
    use dvdif
    use fshlw
    use fshdf
    use tslvt
    use tflxt
    use tovtr
    use fbtav
    use dvlta
    use dvlva
    use dwdns
    use ucloc
    use brstt
    use bstbc
    use qckot
    use bchmk
    use bshft
    
    implicit none

    real(8),        intent(inout)  ::    hx(1:nxydim)
    real(8),        intent(inout)  ::  ubtx(1:nxydim),  vbtx(1:nxydim)
    real(8),        intent(inout)  ::     w(1:nxyzdm),     r(1:nxyzdm)
    real(8),        intent(inout)  ::    ux(1:nxyzdm),    vx(1:nxyzdm)
    real(8),        intent(inout)  ::    tx(1:nxyzdm,1:ntdim)
    real(8),        intent(inout)  ::   amv(1:nxyzdm),   ahv(1:nxyzdm)
    real(8),        intent(in)     ::  taux(1:nxydim),  tauy(1:nxydim) 
    real(8),        intent(in)     ::    ft(1:nxydim,1:ntdim)
    real(8),        intent(in)     ::  ptop(1:nxydim) 
#ifdef OPT_BODY
    real(8),        intent(in)     ::    tq(1:nxyzdm,1:ntdim)
#endif
    real(8),        intent(in)     ::    uy(1:nxyzdm),    vy(1:nxyzdm)
    real(8),        intent(inout)  ::    ty(1:nxyzdm,1:ntdim)
    real(8),        intent(inout)  ::    hy(1:nxydim)
    real(8),        intent(inout)  ::  ubty(1:nxydim),  vbty(1:nxydim)
    real(8),        intent(inout)  ::    gx(1:nxyzdm),    gy(1:nxyzdm)
    real(8),        intent(inout)  ::    xx(1:nxyzdm),    yy(1:nxyzdm)
    real(8),        intent(inout)  ::   gxx(1:nxydim),   gyy(1:nxydim)
    real(8),        intent(inout)  ::  uadv(1:nxyzdm),  vadv(1:nxyzdm)
    real(8),        intent(inout)  ::  wadv(1:nxyzdm)

!---- local variables
    real(8)               ::    tmp(1:nxyzdm,1:ntdim)
    real(8)               ::     hz(1:nxydim),   htmp(1:nxydim)
    real(8)               ::  ubtmp(1:nxydim),  vbtmp(1:nxydim)
    real(8),        save  ::  ubtav(1:nxydim),  vbtav(1:nxydim)
    integer(4)            ::  itsplt,     ij,    ijk,      n

    if ( ( myrank >= ijnode ) .and.                                  &
    &    ( .not. oinit      ) .and.                                  &
    &    ( .not. ofinal     )        ) return

! *** vertical diffusivity and viscosity ***

    call clcstr('COEFF')
    call vdiff (   amv,    ahv,                                       &
    &               uy,     vy,      r,   taux,   tauy,               &
    &               ty,     hy  )
#ifdef OPT_BBL
    call vdiffb(   amv,    ahv  )
#endif

#ifdef OPT_TRIPOLE
    call shift1(   amv,                                               &
    &            nxdim,  nydim,  nzdim,                               &
    &             1.d0,     -1,     -1 )
    call shift1(   ahv,                                               &
    &            nxdim,  nydim,  nzdim,                               &
    &             1.d0,      0,      0 )
#else
    call shift2(   amv,    ahv,  nxdim,  nydim,  nzdim   )
#endif

    call clcend('COEFF')

! *** baroclinic flow ***
    
    call clcstr('BRCLI')
#ifdef OPT_BBL
    call rmmskv
#endif
    call srcvel(    gx,     gy,     xx,     yy,                       &
    &               ux,     vx )
    call vscvel(    gx,     gy,     xx,     yy,                       &
    &               uy,     ux,     vy,     vx,                       &
    &               hx,     hy,                                       &
    &              amv,   taux,   tauy )
    call advvel(                                                      &
    &               gx,     gy,     xx,     yy,                       &    
    &               uy,     ux,     vy,     vx,                       &
    &               hx,     hy,                                       &
    &             uadv,   vadv,   wadv )
#ifdef OPT_BBL
    call srcvlb(    gx,     gy,     xx,     yy,                       &
    &               ux,     vx )
    call vscvlb(    gx,     gy,     xx,     yy,                       &
    &               uy,     ux,     vy,     vx,                       &
    &              amv )
    call advvlb(    gx,     gy,     xx,     yy,                       &
    &               uy,     ux,     vy,     vx,                       &
    &             uadv,   vadv,   wadv )
#endif
    call tnduvd(   gxx,    gyy,                                       &
    &               gx,     gy,                                       &
    &               xx,     yy,                                       &
    &               uy,     vy,     hy,      r,                       &
    &               ux,     vx,   uadv,   vadv )
#ifdef OPT_BBL
    call tnduvb(   gxx,    gyy,                                       &
    &               gx,     gy,                                       &
    &               xx,     yy,                                       &
    &               uy,     vy,     ty,                               &
    &               ux,     vx )
    call admkv1
    call stbbgv(    gx,     gy )
#endif
    call velrds(    ux,     vx,                                       &
    &               gx,     gy,    amv,                               &
    &               hx )
#ifdef OPT_BBL
    call stbbuv(   ux,     vx )
#endif
    call clcend('BRCLI')

! *** barotropic flow and surface elevation ***

    call clcstr('BRTRO')
    if ( oinit .or. ofinal ) then
         call modgxy(   gxx,    gyy,   ubtx,   vbtx  )
    else
       if ( ieuler == 2 ) then
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
          call modgxy(   gxx,    gyy,    ubtx,   vbtx  )

#ifdef OPT_TRIPOLE
          call shift2(   gxx,    gyy,                                 &
    &                  nxdim,  nydim,      1,                         &
    &                  -1.d0,     -1,     -1  )
#else
          call shift2(   gxx,    gyy,   nxdim,  nydim,      1 )
#endif
          call btavst( ubtav,  vbtav )
          do itsplt = 1, ntss

             do ij = 1, nxydim
                htmp (ij) = hx  (ij)
                ubtmp(ij) = ubtx(ij)
                vbtmp(ij) = vbtx(ij)
             end do
             call shalow(   htmp,  ubtmp,  vbtmp,                     &
    &                        hx,   ubtx,   vbtx,                      &
    &                       gxx,    gyy,   ptop,  ft(1,2) )
             call btavad( ubtav,  vbtav,  ubtmp,  vbtmp )
             call shalow(    hx,   ubtx,   vbtx,                      &
    &                      htmp,  ubtmp,  vbtmp,                      &
    &                       gxx,    gyy,   ptop,  ft(1,2) )
#ifdef OPT_TRIPOLE
             call shift2(  ubtx,   vbtx,                              &
    &                     nxdim,  nydim,      1,                      &
    &                     -1.d0,     -1,     -1 )
             call shift1(    hx,                                      &
    &                     nxdim,  nydim,      1,                      &
    &                      1.d0,      0,      0 )
#else
             call shift3(   hx,   ubtx,   vbtx,                       &
    &                    nxdim,  nydim,      1 )
#endif
          end do
#ifdef OPT_TRIPOLE
          call shift2( ubtav,  vbtav,                                 &
    &                  nxdim,  nydim,      1,                         &
    &                  -1.D0,     -1,     -1  )
#else
          call shift2( ubtav,  vbtav,                                 &
    &                  nxdim,  nydim,      1  )
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
    call veltad(   uadv,   vadv,                                      &
    &             ubtav,  vbtav,     uy,     vy  )
#ifdef OPT_TRIPOLE
    if ( .not. ( oinit .or. ofinal ) ) then
       call shift2(   uadv,    vadv,                                  &
    &                nxdim,   nydim,  nzdim,                          &
    &                -1.d0,      -1,     -1 )
    end if
#else
    call shift2(   uadv,   vadv,                                      &
    &             nxdim,  nydim,  nzdim )
#endif

#ifdef OPT_BBL
    call rmmskt
    call admktb
#endif
    call wdenst(  wadv,                                              &
    &             uadv,   vadv,                                      &
    &               hx,     hz )
    call clcend('TDIAG')

! *** tracer ***

    call clcstr('TRACE')
    if ( oinit .or. ofinal ) then
       call flxtrc(     tmp,     xx,                                 &
    &                    tx,     hx,                                 &
    &                    ty,     hz,                                 &
    &                  uadv,   vadv,   wadv,                         &
    &                   ahv )
#ifdef OPT_BBL
       call flxtrb(     tmp,     xx,                                 &
    &                    tx,                                         &
    &                    ty,                                         &
    &                  uadv,   vadv,   wadv,                         &
    &                   ahv )
#endif
       call slvtrc(      tx,     hx,                                 &
    &                   tmp,     xx,                                 &
    &                    ft,     hz )
       call shdiff(      tx,     hx ) 
       call ovturn(       r,     tx,    hx  )
    else
       if ( ieuler == 2 ) then
          do n = 1, ntdim
             do ijk = 1, nxyzdm
                gx(ijk)    = tx(ijk, n)
                tx(ijk, n) = ty(ijk, n)
                ty(ijk, n) = gx(ijk)
             end do
          end do
       else
          call flxtrc(   tmp,    xx,                                  &
    &                    tx,     hx,                                  &
    &                    ty,     hz,                                  &
    &                  uadv,   vadv,   wadv,                          &
    &                   ahv )
#ifdef OPT_BBL
          call flxtrb(   tmp,     xx,                                 &
    &                    tx,                                          &
    &                    ty,                                          &
    &                  uadv,   vadv,   wadv,                          &
    &                   ahv )
          call rmmskt
          call admkt1
          call stbbgt(   tmp,     xx )
#endif
#ifdef OPT_BODY
          do n = 1, ntdim
             do ijk = 1, nxyzdm
                tmp(ijk, n) = tmp(ijk, n) + tq(ijk, n)
             end do
          end do
#endif
          call slvtrc(    tx,     hx,                                 &
    &                    tmp,     xx,                                 &
    &                     ft,     hz  )
#ifdef OPT_TRIPOLE
          call shift1(     tx,                                        &
    &                   nxdim,   nydim, nztdim,                       &
    &                    1.d0,       0,      0 )
          call shift1(     hx,                                        &
    &                   nxdim,   nydim,      1,                       &
    &                    1.d0,       0,      0 )
#else
          call shift1(     tx,                                        &
    &                   nxdim,  nydim, nztdim )
          call shift1(     hx,                                        &
    &                   nxdim,  nydim,      1 )
#endif
          call shdiff(     tx,     hx )
          call ovturn(      r,     tx,    hx )
#ifdef OPT_BBL
          call stbbtr(     tx )
#endif
#ifdef OPT_TRIPOLE
          call shift1(    r,                                          &
    &                 nxdim,   nydim,  nzdim,                         &
    &                  1.d0,       0,      0 )
          call shift1(   tx,                                          &
    &                 nxdim,   nydim, nztdim,                         &
    &                  1.d0,       0,      0 )
          call shift1(   hx,                                          &
    &                 nxdim,   nydim,      1,                         &
    &                  1.D0,       0,      0 )
#else
          call shift1(    r,                                          &
    &                 nxdim,  nydim,  nzdim)
          call shift1(   tx,                                          &
    &                 nxdim,  nydim, nztdim)
          call shift1(   hx,                                          &
    &                 nxdim,  nydim,      1)
#endif
          call stbctr(   tx,      r)
       end if
    end if
    call clcend('TRACE')
      
! *** velocity for momentum advection ***

    call clcstr('VDIAG')
    if ( .not. ( oinit .or. ofinal ) ) then
       do ijk = 1, nxyzdm
          uadv(ijk) = ux(ijk)
          vadv(ijk) = vx(ijk)
       end do
    end if
    call veltad(    ux,     vx,                                       &
    &             ubtx,   vbtx,   uadv,   vadv )
#ifdef OPT_TRIPOLE
    call shift2(    ux,     vx,                                       &
    &            nxdim,  nydim,  nzdim,                               &
    &            -1.d0,     -1,     -1 )
#else
    call shift2(    ux,     vx,                                       &
    &            nxdim,  nydim,  nzdim )
#endif
#ifdef OPT_BBL
    call rmmskt
    call admktb
#endif
    call wdenst(    w,                                                &
    &              ux,     vx,                                        &
    &              hx,     hz )
#ifdef OPT_BBL
    call admskv
#endif
    call velvad( uadv,   vadv,   wadv,                                &
    &              ux,     vx,      w )
#ifdef OPT_BBL
    call stbbvt(   ux,     vx )
#ifdef OPT_TRIPOLE
    call shift2(    ux,     vx,                                       &
    &            nxdim,  nydim,  nzdim,                               &
    &            -1.d0,     -1,     -1 )
#else
    call shift2(    ux,     vx,                                       &
    &            nxdim,  nydim,  nzdim )
#endif
#else
    call stbcuv(    ux,     vx )
#endif
#ifdef OPT_TRIPOLE
    call shift1(  uadv,                                               &
    &            nxdim,  nydim,  nzdim,                               &
    &            -1.d0,      0,     -1 )
    call shift1(  vadv,                                               &
    &            nxdim,  nydim,  nzdim,                               &
    &            -1.d0,     -1,      0 )
    call shift1(  wadv,                                               &
    &            nxdim,  nydim,  nzdim,                               &
    &             1.d0,     -1,     -1 )
#else
    call shift3(  uadv,   vadv,   wadv,                               &
    &            nxdim,  nydim,  nzdim )
#endif
    call clcend('VDIAG')

  end subroutine predco
  
end module aprdc
