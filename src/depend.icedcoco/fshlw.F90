module fshlw

! ---- information ----------------------------------------------------
!
!  HISTORY
!     '03.04.23  H.Hasumi: from COCO3.4
!     '06.01.10  H.Hasumi: bug fix (calculation of SYY)
!     '07.04.23  H.Hasumi: Formulation changed
!     '07.07.30  H.Hasumi: Loop boundary for S?? calculation reduced
!     '08.06.11  H.Hasumi: initial/final processing
!     '08.07.10  H.Hasumi: initial/final processing
!     '12.01.30  Y.Komuro: (change surface water flux by T. Suzuki)
!     '12.06.29  H.Tatebe: for COCO5.0 in F90
! ---------------------------------------------------------------------

  use zocdim,  only :  nxydim, nzdim, ntdim, nxdim, nydim, nx, ny, istr, jstr, kstr
#ifdef OPT_WIDE_BT_SHFT
  use bt_shft, only :  nxydim_w
#endif
  implicit none

  private
  public  ::  modgxy,  shalow,  bvfreq,  ttsp,  ubtwof,  vbtwof, lnovis, otide, brtro_filt

  real(8),     save  ::    fux(nxydim),     fuy(nxydim)
  real(8),     save  ::    fvx(nxydim),     fvy(nxydim)
  real(8),     save  ::    fhx(nxydim),     fhy(nxydim)
  real(8),     save  ::     gu(nxydim),      gv(nxydim)
  real(8),     save  ::    sxx(nxydim),     syy(nxydim)
  real(8),     save  ::    sxy(nxydim),     syx(nxydim)
  real(8),     save  ::     gh(nxydim)

#ifdef OPT_WIDE_BT_SHFT
  real(8),     save  ::  amskt_w(nxydim_w), amskv_w(nxydim_w)
  real(8),     save  ::     ry_w(nxydim_w),   rxt_w(nxydim_w),   ryt_w(nxydim_w)
  real(8),     save  ::    rym_w(nxydim_w),   rxu_w(nxydim_w),   ryu_w(nxydim_w)
  real(8),     save  ::    hxu_w(nxydim_w),   hyu_w(nxydim_w)
  real(8),     save  ::    cor_w(nxydim_w), rdepv_w(nxydim_w)
  real(8),     save  ::     gh_w(nxydim_w),   fhx_w(nxydim_w),   fhy_w(nxydim_w)
  real(8),     save  ::     gu_w(nxydim_w),    gv_w(nxydim_w)
#endif

  real(8),     save  :: amhmod(nxydim)
  real(8),     save  ::   accb,   acc
  data accb, acc / -1.d0, 1.d0 /

  logical,     save  ::  lnovis = .false.
  character(len=16)  ::  chead(1:64)

! for tide
  logical, save :: otide = .false.
  real(8), save :: ubtwof(nxydim), vbtwof(nxydim)
  real(8), save  ::   bvf(nxydim)
  real(8), save  ::   rzm(nxydim, nzdim)
  real(8), save  ::   cfb
  real(8), save, dimension(nzdim) :: c0, c1, c2, c3, c4, c5, c6
  real(8), save, dimension(nzdim) :: d0, d1, d2, d3, d4, d5, d6, d7, d8, d9
  real(8), save  ::  beta = 0.d0, alpha = 1.d0
  real(8), save  :: bvfmax = 1.d30
  real(8), save  ::  rghn(nxydim)
  real(8), save  ::  ttsp

contains
#ifdef OPT_WIDE_BT_SHFT
  subroutine brtro_filt(gxx, gyy, fw, ptop,       &
             &          ubtx, vbtx, hx,           &
             &          ubtav, vbtav)
    use zocgrd, only : ntss
    use bshft

    use bt_shft, only : &
         bt_shift1, bt_shift2, bt_shift3, bt_shift_pack_begin, bt_shift_pack_end, bt_shift_unpack, &
         istr_w, jstr_w, iend_w, jend_w, nxdim_w, nydim_w, ncomm

    real(8), intent(inout) ::   gxx(nxdim, nydim),   gyy(nxdim, nydim), fw(nxdim, nydim)
    real(8), intent(in)    ::  ptop(nxdim, nydim)
    real(8), intent(inout) ::  ubtx(nxdim, nydim),  vbtx(nxdim, nydim), hx(nxdim, nydim)
    real(8), intent(out)   :: ubtav(nxdim, nydim), vbtav(nxdim, nydim)

    real(8) :: ubtav2(nxdim, nydim), vbtav2(nxdim, nydim),  hav(nxdim, nydim)
    
    real(8) ::   gxx_w(nxdim_w, nydim_w)
    real(8) ::   gyy_w(nxdim_w, nydim_w)
    real(8) ::    fw_w(nxdim_w, nydim_w)
    real(8) ::  ptop_w(nxdim_w, nydim_w)

    real(8) ::     h_w(nxdim_w, nydim_w)
    real(8) ::   ubt_w(nxdim_w, nydim_w)
    real(8) ::   vbt_w(nxdim_w, nydim_w)
    real(8) ::  htmp_w(nxdim_w, nydim_w)
    real(8) :: ubtmp_w(nxdim_w, nydim_w)
    real(8) :: vbtmp_w(nxdim_w, nydim_w)

    real(8) ::   fact
    integer ::   i, j, nb, itsplt
    logical,     save  ::  ofirst = .true.

    if (ofirst) then
       !$acc enter data create(ubtav2, vbtav2, hav)
       !$acc enter data create(gxx_w, gyy_w, fw_w, ptop_w, h_w, ubt_w, vbt_w)
       !$acc enter data create(htmp_w, ubtmp_w, vbtmp_w)
    end if
    
    !$acc kernels default(present)
    do j =1, nydim_w
       do i =1, nxdim_w
          gxx_w  (i,j)=0.d0
          gyy_w  (i,j)=0.d0
          ptop_w (i,j)=0.d0
          fw_w   (i,j)=0.d0
          h_w    (i,j)=0.d0
          ubt_w  (i,j)=0.d0
          vbt_w  (i,j)=0.d0
       end do
    end do
    !$acc end kernels
    
    call modgxy_w(gxx, gyy, ubtx, vbtx)

    !$acc kernels default(present)
    do j =1, ny
       do i =1, nx
          gxx_w (i+istr_w-1, j+jstr_w-1) = gxx (i+istr-1, j+jstr-1)
          gyy_w (i+istr_w-1, j+jstr_w-1) = gyy (i+istr-1, j+jstr-1)
          fw_w  (i+istr_w-1, j+jstr_w-1) = fw  (i+istr-1, j+jstr-1)
          ptop_w(i+istr_w-1, j+jstr_w-1) = ptop(i+istr-1, j+jstr-1)
          h_w   (i+istr_w-1, j+jstr_w-1) = hx  (i+istr-1, j+jstr-1)
          ubt_w (i+istr_w-1, j+jstr_w-1) = ubtx(i+istr-1, j+jstr-1)
          vbt_w (i+istr_w-1, j+jstr_w-1) = vbtx(i+istr-1, j+jstr-1)
       end do
    end do
    !$acc end kernels
    
    call bt_shift_pack_begin
    call bt_shift2( gxx_w,  gyy_w,      nxdim_w, nydim_w, 1, -1.d0, -1, -1)
    call bt_shift3(  fw_w, ptop_w, h_w, nxdim_w, nydim_w, 1,  1.d0,  0,  0)
    call bt_shift2( ubt_w,  vbt_w,      nxdim_w, nydim_w, 1, -1.d0, -1, -1)
    call bt_shift_pack_end
    call bt_shift_unpack( gxx_w, 1)
    call bt_shift_unpack( gyy_w, 2)
    call bt_shift_unpack(  fw_w, 3)
    call bt_shift_unpack(ptop_w, 4)
    call bt_shift_unpack(   h_w, 5)
    call bt_shift_unpack( ubt_w, 6)
    call bt_shift_unpack( vbt_w, 7)
    
    nb = ntss * 2
    !$acc kernels default(present)
    do j = 1, nydim
       do i = 1, nxdim
          ubtav (i,j) = 0.d0
          vbtav (i,j) = 0.d0
          ubtav2(i,j) = 0.d0
          vbtav2(i,j) = 0.d0
          hav   (i,j) = hx(i,j) / dble(nb+1)
       end do
    end do
    !$acc end kernels
          
    do itsplt = 1, nb
       !$acc kernels default(present)
       do j = 1, nydim_w
          do i = 1, nxdim_w
             htmp_w (i,j) = h_w  (i,j)
             ubtmp_w(i,j) = ubt_w(i,j)
             vbtmp_w(i,j) = vbt_w(i,j)
          end do
       end do
       !$acc end kernels
       call shalow_w(                                           &
    &              htmp_w,  ubtmp_w,  vbtmp_w,                  &
    &                 h_w,    ubt_w,    vbt_w,                  &
    &               gxx_w,    gyy_w,   ptop_w,  fw_w)

       call shalow_w(                                           &
    &                 h_w,    ubt_w,    vbt_w,                  &
    &              htmp_w,  ubtmp_w,  vbtmp_w,                  &
    &               gxx_w,    gyy_w,   ptop_w,  fw_w)

       if (mod(itsplt, ncomm/2) == 0) then
          call bt_shift_pack_begin
          call bt_shift2(ubt_w, vbt_w, nxdim_w, nydim_w, 1, -1.d0, -1, -1)
          call bt_shift1(  h_w,        nxdim_w, nydim_w, 1,  1.d0,  0,  0)
          call bt_shift_pack_end
          call bt_shift_unpack(ubt_w, 1)
          call bt_shift_unpack(vbt_w, 2)
          call bt_shift_unpack(  h_w, 3)
       end if
       fact = 2.d0 * dble(nb-itsplt+1) / dble(nb * (nb+1))
       !$acc kernels default(present)
       do j = 1, ny
          do i = 1, nx
             ubtav (i+istr-1, j+jstr-1) = ubtav (i+istr-1, j+jstr-1) + ubtmp_w(i+istr_w-1, j+jstr_w-1) * fact
             vbtav (i+istr-1, j+jstr-1) = vbtav (i+istr-1, j+jstr-1) + vbtmp_w(i+istr_w-1, j+jstr_w-1) * fact
             ubtav2(i+istr-1, j+jstr-1) = ubtav2(i+istr-1, j+jstr-1) + ubtmp_w(i+istr_w-1, j+jstr_w-1) / dble(nb)
             vbtav2(i+istr-1, j+jstr-1) = vbtav2(i+istr-1, j+jstr-1) + vbtmp_w(i+istr_w-1, j+jstr_w-1) / dble(nb)
             hav   (i+istr-1, j+jstr-1) = hav   (i+istr-1, j+jstr-1) + h_w    (i+istr_w-1, j+jstr_w-1) / dble(nb+1)
          end do
       end do
       !$acc end kernels
    end do

    call shift_pack_begin
    call shift2( ubtav,  vbtav, nxdim,  nydim, 1, -1.d0, -1, -1)
    call shift2(ubtav2, vbtav2, nxdim,  nydim, 1, -1.d0, -1, -1)
    call shift1(   hav,         nxdim,  nydim, 1,  1.d0,  0,  0)
    call shift_pack_end
    call shift_unpack( ubtav, 1)
    call shift_unpack( vbtav, 2)
    call shift_unpack(ubtav2, 3)
    call shift_unpack(vbtav2, 4)
    call shift_unpack(   hav, 5)

    !$acc kernels default(present)
    do j = 1, nydim
       do i = 1, nxdim
          hx  (i,j) = hav   (i,j)
          ubtx(i,j) = ubtav2(i,j)
          vbtx(i,j) = vbtav2(i,j)
       end do
    end do
    !$acc end kernels
  end subroutine brtro_filt

  subroutine modgxy_w(                                                &
         &      gxx,    gyy,                                          &
         &      ubtx,   vbtx  )

    use zocphy,   only :  gravit, rhoo
    use zocgrd,   only :  ry, rxt, ryt, rym, rxu, ryu, hxu, hyu, rdepv, cor
    use bt_shft, only : bt_shift1, bt_shift2, bt_shift3, nxdim_w, nydim_w, nxydim_w, istr_w, jstr_w, ncomm
    use zocmsk,  only :  amskt,  amskv
    use ufile
    implicit none
    real(8),   intent(in)  ::    gxx(nxydim),   gyy(nxydim)
    real(8),   intent(in)  ::   ubtx(nxydim),  vbtx(nxydim)
    logical,     save  ::  ofirst = .true.
    integer(4)  ::  ifpar, jfpar, istat, ij, ij_w, i, j

    namelist /nmaccb/ accb
    namelist /nmaccv/ acc
    namelist /nmnovis/ lnovis
    namelist /nmtide/  otide
    
    if ( ofirst ) then
       ofirst = .false.
       call rewnml( ifpar, jfpar )
       read( ifpar, nmaccb, iostat = istat ) 
       call cstnml( jfpar, 'modgxy_w', 'nmaccb', istat )
       write( jfpar, nmaccb )
       call rewnml( ifpar, jfpar )
       read( ifpar, nmaccv, iostat = istat ) 
       call cstnml( jfpar, 'modgxy_w', 'nmaccv', istat )
       write( jfpar, nmaccv )
       call rewnml( ifpar, jfpar )
       read( ifpar, nmnovis, iostat = istat )
       call cstnml( jfpar, 'modgxy_w', 'nmnovis', istat )
       write( jfpar, nmnovis )
       call rewnml( ifpar, jfpar )
       read( ifpar, nmtide, iostat = istat )
       call cstnml( jfpar, 'modgxy_w', 'nmtide', istat )
       write( jfpar, nmtide )

       if (ncomm > min(nx,ny)) then
          write(jfpar,*) 'Error. ncomm should be <= min(nx,ny)'
          stop
       end if
       if( (.not. lnovis) .or. otide) then
          write(jfpar,*) 'Error. OPT_WIDE_BT_SHFT is currently applicable ' &
                       //'for lnovis=.ture.  and otide=.false. '
          stop
       end if

       if ( accb <= 0.d0 ) then
          accb = acc
       end if
       !$acc enter data create(amskt_w, amskv_w)
       !$acc enter data create(ry_w, rxt_w, ryt_w)
       !$acc enter data create(rym_w, rxu_w, ryu_w)
       !$acc enter data create(hxu_w, hyu_w)
       !$acc enter data create(cor_w, rdepv_w)
       !$acc enter data create(gh_w, fhx_w, fhy_w)
       !$acc enter data create(gu_w, gv_w)
       
       !$acc kernels default(present)
       do ij = 1, nxydim_w
          amskt_w(ij) = 0.d0
          amskv_w(ij) = 0.d0
          ry_w   (ij) = 0.d0
          rxt_w  (ij) = 0.d0
          ryt_w  (ij) = 0.d0
          rym_w  (ij) = 0.d0
          rxu_w  (ij) = 0.d0
          ryu_w  (ij) = 0.d0
          hxu_w  (ij) = 0.d0
          hyu_w  (ij) = 0.d0
          cor_w  (ij) = 1.d0
          rdepv_w(ij) = 1.d0
          gh_w   (ij) = 0.d0
       end do
       do j = 1, ny
          do i = 1, nx
             ij   = (i+ istr  -1)  + (j+jstr  -1 -1) *nxdim
             ij_w = (i+ istr_w-1)  + (j+jstr_w-1 -1) *nxdim_w
             amskt_w(ij_w) = amskt(ij, kstr)
             amskv_w(ij_w) = amskv(ij, kstr)
             ry_w   (ij_w) = ry   (ij)
             rxt_w  (ij_w) = rxt  (ij)
             ryt_w  (ij_w) = ryt  (ij)
             rym_w  (ij_w) = rym  (ij)
             rxu_w  (ij_w) = rxu  (ij)
             ryu_w  (ij_w) = ryu  (ij)
             hxu_w  (ij_w) = hxu  (ij)
             hyu_w  (ij_w) = hyu  (ij)
             cor_w  (ij_w) = cor  (ij)
             rdepv_w(ij_w) = rdepv(ij)
             gh_w   (ij_w) = gravit / rdepv(ij)
          end do
       end do
       !$acc end kernels
       call bt_shift1(amskt_w,               nxdim_w,  nydim_w, 1,  1.d0,  0,  0)
       call bt_shift1(amskv_w,               nxdim_w,  nydim_w, 1,  1.d0, -1, -1)
       call bt_shift3(  ry_w, rxt_w, ryt_w,  nxdim_w,  nydim_w, 1,  1.d0,  0,  0)
       call bt_shift3( rym_w, rxu_w, ryu_w,  nxdim_w,  nydim_w, 1,  1.d0, -1, -1)
       call bt_shift2( hxu_w, hyu_w,         nxdim_w,  nydim_w, 1,  1.d0, -1, -1)
       call bt_shift3(cor_w, rdepv_w, gh_w,  nxdim_w,  nydim_w, 1,  1.d0, -1, -1)
    end if
  end subroutine modgxy_w

  subroutine shalow_w(                                                &
         &       hx,   ubtx,   vbtx,                                  &
         &       hy,   ubty,   vbty,                                  &
         &      gxx,    gyy,   ptop,    fw )

    use zocphy,   only :  rhoo
    use zocgrd,   only :  tss, rx 
    use bt_shft, only :  nxydim_w, nxdim_w
    implicit none
 
    real(8),   intent(inout)  ::    hx(nxydim_w),  ubtx(nxydim_w),  vbtx(nxydim_w)
    real(8),   intent(in)     ::    hy(nxydim_w),  ubty(nxydim_w),  vbty(nxydim_w)
    real(8),   intent(in)     ::   gxx(nxydim_w),   gyy(nxydim_w)
    real(8),   intent(in)     ::  ptop(nxydim_w),    fw(nxydim_w)

    integer, parameter :: le = 1, lw = -1, ln = nxdim_w, ls = -nxdim_w
    integer, parameter :: lne = nxdim_w+1, lsw = -nxdim_w-1

    real(8)     ::     cf
    integer(4)  ::     ij

    !$acc kernels default(present) async
    do ij = 1, nxydim_w
       fhx_w(ij) = 0.d0
       fhy_w(ij) = 0.d0
    end do

    do ij = nxdim_w+2, nxydim_w
       fhx_w(ij) = - (  ubty(ij+lw)  * hyu_w(ij+lw)                   &
    &                 + ubty(ij+lsw) * hyu_w(ij+lsw)) * 0.5d0
    end do

    do ij = nxdim_w+2, nxydim_w
       fhy_w(ij) = - (  vbty(ij+ls)  * hxu_w(ij+ls)                   &
    &                 + vbty(ij+lsw) * hxu_w(ij+lsw)) * 0.5d0
    end do

    do ij = 1, nxydim_w - ln
       hx(ij) = hx(ij)                                                &
    &         + tss * (  (fhx_w(ij+le) - fhx_w(ij)) * rx              &
    &                  + (fhy_w(ij+ln) - fhy_w(ij)) * ry_w(ij)) *     &
    &           rxt_w(ij) * ryt_w(ij) * amskt_w(ij)                   &
    &         - tss * fw(ij) * amskt_w(ij)
    end do
    !$acc end kernels

    !$acc kernels default(present) async
    do ij = 1, nxydim_w - lne
       gu_w(ij) = gxx(ij) + cor_w(ij) * vbtx(ij)                     &
    &         - (  ( (hy(ij+lne)+hy(ij+le))                          &
    &               -(hy(ij+ln )+hy(ij   )) ) * gh_w(ij)             &
    &            + ( (ptop(ij+lne)+ptop(ij+le))                      &
    &               -(ptop(ij+ln )+ptop(ij   )) )                    &
    &              / rdepv_w(ij) / rhoo                              &
    &           ) * 0.5d0 * rx * rxu_w(ij)
       gv_w(ij) = gyy(ij) - cor_w(ij) * ubtx(ij)                     &
    &         - (  ( (hy(ij+lne)+hy(ij+ln))                          &
    &               -(hy(ij+le )+hy(ij   )) ) * gh_w(ij)             &
    &            + ( (ptop(ij+lne)+ptop(ij+ln))                      &
    &               -(ptop(ij+le )+ptop(ij   )) )                    &
    &              / rdepv_w(ij) / rhoo                              &
    &           ) * 0.5d0 * rym_w(ij) * ryu_w(ij)
    end do

    do ij = 1, nxydim_w - lne
       cf = cor_w(ij) * tss / accb * 0.5d0
       ubtx(ij) = ubtx(ij)                                           &
    &           + tss / accb / (1.d0 + cf * cf) *                    &
    &             (gu_w(ij) + cf * gv_w(ij)) * amskv_w(ij)
       vbtx(ij) = vbtx(ij)                                           &
    &           + tss / accb / (1.d0 + cf * cf) *                    &
    &             (gv_w(ij) - cf * gu_w(ij)) * amskv_w(ij)
    end do
    !$acc end kernels
    !$acc wait
  end subroutine shalow_w

#else
  subroutine brtro_filt(gxx, gyy, fw, ptop,       &
             &          ubtx, vbtx, hx,           &
             &          ubtav, vbtav)
    use zocgrd, only : ntss
    use bshft
    real(8), intent(inout) ::   gxx(nxydim),   gyy(nxydim), fw(nxydim)
    real(8), intent(in)    ::   ptop(nxydim)
    real(8), intent(inout) ::  ubtx(nxydim),  vbtx(nxydim), hx(nxydim)
    real(8), intent(out)   :: ubtav(nxydim), vbtav(nxydim)

    real(8) ::   htmp(nxydim),  ubtmp(nxydim), vbtmp(nxydim)
    real(8) :: ubtav2(nxydim), vbtav2(nxydim),  hav(nxydim)
    real(8) ::   fact
    integer ::   ij, nb, itsplt
    logical     :: ofirst=.true.

    if (ofirst) then
       !$acc enter data create(htmp, ubtmp, vbtmp)
       !$acc enter data create(ubtav2, vbtav2, hav)
    end if

    call modgxy(gxx, gyy, ubtx, vbtx)

    call shift_pack_begin
    call shift2(gxx, gyy, nxdim, nydim, 1, -1.d0, -1, -1)
    call shift1( fw,      nxdim, nydim, 1,  1.d0,  0,  0)
    call shift_pack_end
    call shift_unpack( gxx, 1)
    call shift_unpack( gyy, 2)
    call shift_unpack(  fw, 3)
    
    nb = ntss * 2
    !$acc kernels default(present)
    do ij = 1, nxydim
       ubtav (ij) = 0.d0
       vbtav (ij) = 0.d0
       ubtav2(ij) = 0.d0
       vbtav2(ij) = 0.d0
       hav   (ij) = hx(ij) / dble(nb+1)
    end do
    !$acc end kernels
          
    do itsplt = 1, nb
       !$acc kernels default(present)
       do ij = 1, nxydim
          htmp (ij) = hx  (ij)
          ubtmp(ij) = ubtx(ij)
          vbtmp(ij) = vbtx(ij)
       end do
       !$acc end kernels
       call shalow(                                             &
    &              htmp,  ubtmp,  vbtmp,                        &
    &                hx,   ubtx,   vbtx,                        &
    &               gxx,    gyy,   ptop,  fw)

       if( (.not. lnovis) .or. otide) then
          call shift_pack_begin
          call shift2(ubtmp, vbtmp, nxdim, nydim, 1, -1.d0, -1, -1)
          call shift1( htmp,        nxdim, nydim, 1,  1.d0,  0,  0)
          call shift_pack_end
          call shift_unpack(ubtmp, 1)
          call shift_unpack(vbtmp, 2)
          call shift_unpack( htmp, 3)
       end if

       call shalow(                                             &
    &                hx,   ubtx,   vbtx,                        &
    &              htmp,  ubtmp,  vbtmp,                        &
    &               gxx,    gyy,   ptop,  fw)

       call shift_pack_begin
       call shift2(ubtx, vbtx, nxdim, nydim, 1, -1.d0, -1, -1)
       call shift1(  hx,       nxdim, nydim, 1,  1.d0,  0,  0)
       call shift_pack_end
       call shift_unpack(ubtx, 1)
       call shift_unpack(vbtx, 2)
       call shift_unpack(  hx, 3)

       fact = 2.d0 * dble(nb-itsplt+1) / dble(nb * (nb+1))
       !$acc kernels default(present)
       do ij = 1, nxydim
          ubtav (ij) = ubtav (ij) + ubtmp(ij) * fact
          vbtav (ij) = vbtav (ij) + vbtmp(ij) * fact
          ubtav2(ij) = ubtav2(ij) + ubtmp(ij) / dble(nb)
          vbtav2(ij) = vbtav2(ij) + vbtmp(ij) / dble(nb)
          hav   (ij) = hav   (ij) + hx   (ij) / dble(nb+1)
       end do
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

    !$acc kernels default(present)
    do ij = 1, nxydim
       hx  (ij) = hav   (ij)
       ubtx(ij) = ubtav2(ij)
       vbtx(ij) = vbtav2(ij)
    end do
    !$acc end kernels
  end subroutine brtro_filt
#endif

  subroutine modgxy(                                                  &
         &      gxx,    gyy,                                          &
         &      ubtx,   vbtx  )

    use zocdim,  only :                                               &
         &     nxg,    nyg, nxgdim, nygdim,                           &
         &  nxydim,  nxdim,  nydim,  nzdim,                           &
         &    kstr,   kend,     kz,                                   &
         &   ijstr,  ijend, ijvstr, ijvend,                           &
         &   igstr,  jgstr,                                           &
         &      le,     lw,     ln,     ls,                           &
         &     lnw,    lne,    lse,    lsw,                           &
         &   oinit,  ofinal
    use zocphy,   only :  gravit, rhoo
    use zocgrd,   only :                                              &
         &   rdepv,                                                   &  
         &      dy,    hxt,    rea,                                   &  
         &      rx,     ry,    rym,                                   & 
         &     hxu,    hyu,   hxyu,   hyxu,                           &
         &     rxu,    ryu,                                           &
         &     dzm
    use zocmsk,  only :  amskv,  amfvx,  amfvy, nbot
    use zocnod,  only :  iroot,  myrank
    use zocfil,  only :  ncf
    use ufile
#ifdef OPT_IO_COCOMPI
    use mpiio
#else
    use bgs2d
#endif
    use bshft
    use xprst

    implicit none
#include "mpif.h"
    
    real(8),   intent(inout)  ::    gxx(nxydim),   gyy(nxydim)
    real(8),   intent(in)     ::   ubtx(nxydim),  vbtx(nxydim)

!----- local variables
    integer(4)  ::     ij,      i,      j
    integer(4)  ::   ijle,   ijln,   ijlw,   ijls
    integer(4)  ::   ijnw,   ijse,   ijsw
    integer(4)  ::  ifpar,  jfpar,   istat

    namelist /nmaccb/ accb
    namelist /nmaccv/ acc

    real(8),     save  ::  amh
    integer(4)         ::  iam,    nfamh
    character(len=ncf) ::  cfamh
    namelist /nmvish/  amh
    namelist /nmcvis/  iam, cfamh
    namelist /nmnovis/ lnovis
    data amh, iam, cfamh / 0.d0, 0, 'not-specified' /

#ifdef OPT_IO_COCOMPI
    integer :: mpi_fh
    integer :: icread
    integer (kind = mpi_offset_kind) :: disp
#else
    real(8), allocatable :: buf2(:, :),  g2d(:, :)
#endif

!---- for tide
    integer :: nfrghn, ksfr = 1
    character(len=ncf) :: cfrghn = 'not-specified'
    real(8) :: lscale = 1.d4, pi, rghmax = 1.d30
    real(8) :: rght(nxydim)
    namelist /nmtide/  otide
    namelist /nmtsal/   beta
    namelist /nmtbdy/  alpha
    namelist /nmbfjl/ cfrghn, lscale, ksfr, rghmax, bvfmax
    
    logical,     save  ::  ofirst = .true.

    if ( oinit .or. ofinal ) then
       return
    end if

    if ( ofirst ) then

       ofirst = .false.
       call rewnml( ifpar, jfpar )
       read( ifpar, nmaccb, iostat = istat ) 
       call cstnml( jfpar, 'modgxy', 'nmaccb', istat )
       write( jfpar, nmaccb )
       call rewnml( ifpar, jfpar )
       read( ifpar, nmaccv, iostat = istat ) 
       call cstnml( jfpar, 'modgxy', 'nmaccv', istat )
       write( jfpar, nmaccv )
       call rewnml( ifpar, jfpar )
       read( ifpar, nmvish, iostat = istat ) 
       call cstnml( jfpar, 'modgxy', 'nmvish', istat )
       write( jfpar, nmvish )
       call rewnml( ifpar, jfpar )
       read( ifpar, nmcvis, iostat = istat ) 
       call cstnml( jfpar, 'modgxy', 'nmcvis', istat )
       write( jfpar, nmcvis )
       call rewnml( ifpar, jfpar )
       read( ifpar, nmnovis, iostat = istat ) 
       call cstnml( jfpar, 'modgxy', 'nmnovis', istat )
       write( jfpar, nmnovis )
       call rewnml( ifpar, jfpar )
       read( ifpar, nmtide, iostat = istat ) 
       call cstnml( jfpar, 'modgxy', 'nmtide', istat )
       write( jfpar, nmtide )
       if ( otide ) then
          call rewnml( ifpar, jfpar )
          read( ifpar, nmtsal, iostat = istat ) 
          call cstnml( jfpar, 'modgxy', 'nmtsal', istat )
          write( jfpar, nmtsal )
          call rewnml( ifpar, jfpar )
          read( ifpar, nmtbdy, iostat = istat ) 
          call cstnml( jfpar, 'modgxy', 'nmtbdy', istat )
          write( jfpar, nmtbdy )
          call rewnml( ifpar, jfpar )
          read( ifpar, nmbfjl, iostat = istat ) 
          call cstnml( jfpar, 'modgxy', 'nmbfjl', istat )
          write( jfpar, nmbfjl )
       end if
       
       if ( accb <= 0.d0 ) then
          accb = acc
       end if
       
       do ij = 1, nxydim
          gh(ij) = gravit / rdepv(ij)
       end do

       if ( otide ) then
          do ij = 1, nxydim
             gh(ij) = gh(ij) * (1.d0 - beta)
          end do
          call secoef( &
               & c0(kstr), c1(kstr), c2(kstr), c3(kstr), &
               & c4(kstr), c5(kstr), c6(kstr), &
               & d0(kstr), d1(kstr), d2(kstr), d3(kstr),  d4(kstr), &
               & d5(kstr), d6(kstr), d7(kstr), d8(kstr),  d9(kstr) )

          rzm(:,:) = 1.d0 / dzm(:,:)
#ifdef OPT_IO_COCOMPI
          call mpi_filopn(mpi_fh, cfrghn, 'READ')
          disp=0
          call mpi_read_2d(rght, mpi_fh  , disp)
          call mpi_filcls(mpi_fh)
#else
          allocate ( buf2(1:nxg,1:nyg), g2d(1:nxgdim,1:nygdim) )
          buf2(1:nxg,1:nyg) = 0.d0
          g2d (1:nxgdim,1:nygdim) = 0.d0
          if ( myrank == iroot ) then
             call filopn( nfrghn, cfrghn, 'READ' )
             rewind( nfrghn )
             read( nfrghn ) chead
             read( nfrghn ) buf2
             do j = 1, nyg
                do i = 1, nxg
                   g2d(igstr+i-1, jgstr+j-1) = buf2(i, j)
                end do
             end do
             call filcls( nfrghn )
          end if
          call scatter_2d( rght, g2d )
          deallocate ( buf2, g2d )
#endif
          !$acc enter data copyin(rght)
          !$acc enter data create(rghn)

          call shift1(rght, nxdim, nydim, 1, 1.d0, 0, 0)

          ksfr = ksfr + kstr-1
          pi = 4.d0 * atan(1.d0)
          lscale = pi / lscale
          !$acc kernels default(present)
          do ij = ijvstr, ijvend
             rghn(ij) = 0.25d0* ( rght(ij   )+ rght(ij+le ) &
                  &              +rght(ij+ln)+ rght(ij+lne))
             rghn(ij) = min( rghmax, rghn(ij) )
             if ( min(nbot(ij), nbot(ij+le), &
                  &   nbot(ij+ln), nbot(ij+lne)) .ge. ksfr ) then
                rghn(ij) = lscale * rghn(ij) * 1.d2 * 0.25d0 * amskv(ij, kstr)
             else
                rghn(ij) = 0.d0
             end if
          end do
          !$acc end kernels
          cfb = gravit / rhoo * 1.d-3

          call shift1(rghn, nxdim, nydim, 1, 1.d0, 0, 0)
       end if
       
       if ( lnovis ) then
!---- fshlw.novis setting
          amhmod(1:nxydim) = 0.d0
          write(jfpar,*) ' AMH is zero in fshlw'
       else
       
       if ( iam < 0 ) then
!---- spatially constant
          write(jfpar,*) 'spatially constant AMH'
          amhmod(1:nxydim) = amh
       end if
       if ( iam == 0 ) then
!---- zonal resolution dependent (cvisc.clat setting)
          write(jfpar,*) 'dx dependent AMH'
          do ij = 1, nxydim
             amhmod(ij) = amh * hxt(ij) / rea
          end do
       end if
       if ( iam == 1 ) then
!---- given by input file
          write(jfpar,*) 'AMH is given by a file'          
#ifdef OPT_IO_COCOMPI
          call mpi_filopn(mpi_fh, cfamh, 'READ')
          disp=0
          call mpi_read_chead(chead, mpi_fh, disp, icread)
          call mpi_read_2d(amhmod, mpi_fh  , disp)
          call mpi_filcls(mpi_fh)
#else
          allocate ( buf2(1:nxg,1:nyg), g2d(1:nxgdim,1:nygdim) )
          buf2(1:nxg,1:nyg) = 0.d0
          g2d (1:nxgdim,1:nygdim) = 0.d0
          if ( myrank == iroot ) then
             call filopn( nfamh, cfamh, 'READ' )
!---- for MIROC
!               CALL IFLOPN(
!     O                     NFAMH,    IERR,
!     I                     CFAMH,  'READ', 'UNFORMATTED') 
             rewind( nfamh )
             read( nfamh ) chead
             read( nfamh ) buf2
             do j = 1, nyg
                do i = 1, nxg
                   g2d(igstr+i-1, jgstr+j-1) = buf2(i, j)
                end do
             end do
             call filcls( nfamh )
!---- for MIROC
!               CLOSE(NFAMH)
          end if
          call scatter_2d( amhmod, g2d )
          deallocate ( buf2, g2d )
#endif
          call shift1(amhmod, nxdim, nydim, 1, 1.d0, 0, 0)
       end if ! iam
       end if ! lnovis
       !$acc enter data copyin(gh, amhmod)
       !$acc enter data create(fux,     fuy)
       !$acc enter data create(fvx,     fvy)
       !$acc enter data create(fhx,     fhy)
       !$acc enter data create( gu,      gv)
       !$acc enter data create(sxx,     syy)
       !$acc enter data create(sxy,     syx)
       
       !$acc enter data create(ubtwof, vbtwof)
       !$acc enter data create(bvf)
       !$acc enter data copyin(rzm)
       
       !$acc kernels
       ubtwof(:)=0.d0
       vbtwof(:)=0.d0
       !$acc end kernels
    end if
    
    if (lnovis) return
!$acc kernels default(present)
!$omp parallel do privete( ij )
    do ij = 1, nxydim
       fux(ij) = 0.d0
       fuy(ij) = 0.d0
       fvx(ij) = 0.d0
       fvy(ij) = 0.d0
    end do
!$omp end parallel do
       
!$omp parallel do &
!$omp private( ij, ijls, ijlw, ijln, ijle, ijnw, ijse, ijsw )
    do ij = ijstr-nxdim-1, ijend+nxdim+1
       ijls = ij + ls
       ijlw = ij + lw
       ijln = ij + ln
       ijle = ij + le
       ijnw = ij + lnw
       ijse = ij + lse
       ijsw = ij + lsw
       sxx(ij) = (ubtx(ij) - ubtx(ijlw)) * amhmod(ij) *               &
    &             rx / (hxu(ij) + hxu(ijlw)) * 2.d0                   &
    &          + (vbtx(ij) + vbtx(ijlw)) * amhmod(ij) *               &
    &            (hxyu(ij) + hxyu(ijlw)) * 0.25d0                     &
    &          - (  vbtx(ijln) + vbtx(ijnw)                           & 
    &             - vbtx(ijls) - vbtx(ijsw)) * amhmod(ij)             &
    &            / (hyu(ij) + hyu(ijlw)) / (dy(ij) + dy(ijln))        &
    &          - (ubtx(ij) + ubtx(ijlw)) * amhmod(ij) *               &
    &            (hyxu(ij) + hyxu(ijlw)) * 0.25d0
       sxy(ij) = (ubtx(ij) - ubtx(ijls)) * amhmod(ij) *               &
    &            ry(ij) / (hyu(ij) + hyu(ijls)) * 2.d0                &
    &          - (ubtx(ij) + ubtx(ijls)) * amhmod(ij) *               &
    &            (hxyu(ij) + hxyu(ijls)) * 0.25d0                     &
    &          + (  vbtx(ijle) + vbtx(ijse)                           &
    &             - vbtx(ijlw) - vbtx(ijsw)) * amhmod(ij)             &
    &            / (hxu(ij) + hxu(ijls)) * rx * 0.5d0                 &
    &          - (vbtx(ij) + vbtx(ijls)) * amhmod(ij) *               &
    &           (hyxu(ij) + hyxu(ijls)) * 0.25d0
       syy(ij) = (vbtx(ij) - vbtx(ijls)) * amhmod(ij) *               &
    &            ry(ij) / (hyu(ij) + hyu(ijls)) * 2.d0                &
    &          + (ubtx(ij) + ubtx(ijls)) * amhmod(ij) *               &
    &            (hyxu(ij) + hyxu(ijls)) * 0.25d0                     &
    &          - (  ubtx(ijle) + ubtx(ijse)                           &
    &             - ubtx(ijlw) - ubtx(ijsw)) * amhmod(ij)             &
    &            / (hxu(ij) + hxu(ijls)) * rx * 0.5d0                 &
    &          - (vbtx(ij) + vbtx(ijls)) * amhmod(ij) *               & 
    &            (hxyu(ij) + hxyu(ijls)) * 0.25d0
       syx(ij) = (vbtx(ij) - vbtx(ijlw)) * amhmod(ij) *               &
    &            rx / (hxu(ij) + hxu(ijlw)) * 2.d0                    &
    &          - (vbtx(ij) + vbtx(ijlw)) * amhmod(ij) *               &
    &            (hyxu(ij) + hyxu(ijlw)) * 0.25d0                     &
    &          + (  ubtx(ijln) + ubtx(ijnw)                           &
    &             - ubtx(ijls) - ubtx(ijsw)) * amhmod(ij)             &
    &            / (hyu(ij) + hyu(ijlw)) / (dy(ij) + dy(ijln))        &
    &          - (ubtx(ij) + ubtx(ijlw)) * amhmod(ij) *               &
    &            (hxyu(ij) + hxyu(ijlw)) * 0.25d0
    end do
!$omp end parallel do
    
!$omp parallel do &
!$omp private( ij, ijlw )
    do ij = ijvstr-nxdim-1, ijvend+1
       ijlw = ij + lw
       fux(ij) = sxx(ij) * (hyu(ij) + hyu(ijlw)) *                   &
    &                      (hyu(ij) + hyu(ijlw)) * 0.25d0
       fvx(ij) = syx(ij) * (hyu(ij) + hyu(ijlw)) *                   &
    &                      (hyu(ij) + hyu(ijlw)) * 0.25d0
    end do
!$omp end parallel do
!$omp parallel do &
!$omp private( ij, ijls )
    do ij = ijvstr-nxdim-1, ijvend+nxdim
       ijls = ij + ls
       fuy(ij) = sxy(ij) * (hxu(ij) + hxu(ijls)) *                   &
    &                      (hxu(ij) + hxu(ijls)) * 0.25d0
       fvy(ij) = syy(ij) * (hxu(ij) + hxu(ijls)) *                   &
    &                      (hxu(ij) + hxu(ijls)) * 0.25d0
    end do
!$omp end parallel do

!$omp parallel do private( ij )
    do ij = ijvstr-nxdim-1, ijvend+nxdim+1
       gxx(ij) = gxx(ij)                                             &
    &          - (  (fux(ij+le) - fux(ij)) * rx * ryu(ij)            & 
    &             + (fuy(ij+ln) - fuy(ij)) * rym(ij) * rxu(ij)) *    &
    &            rxu(ij) * ryu(ij) * amskv(ij, kstr)
       gyy(ij) = gyy(ij)                                             &
    &          - (  (fvx(ij+le) - fvx(ij)) * rx * ryu(ij)            &
    &             + (fvy(ij+ln) - fvy(ij)) * rym(ij) * rxu(ij)) *    &
    &            rxu(ij) * ryu(ij) * amskv(ij, kstr)
    end do
!$omp end parallel do
!$acc end kernels
  end subroutine modgxy

! =====================================================================
  subroutine shalow(                                                  &
         &       hx,   ubtx,   vbtx,                                  &
         &       hy,   ubty,   vbty,                                  &
         &      gxx,    gyy,   ptop,    fw )
    use zocdim,  only :                                               &
         &     nxg,    nyg, nxgdim, nygdim,                           &
         &  nxydim,  nxdim,  nzdim,                                   &
         &    kstr,   kend,     kz,                                   &
         &   ijstr,  ijend, ijtstr, ijtend, ijvstr, ijvend,           &
         &      le,     lw,     ln,     ls,                           &
         &     lnw,    lne,    lse,    lsw
    use zocphy,   only :  rhoo
    use zocgrd,   only :                                              &
         &     cor,  rdepv,    tss,                                   &  
         &      dy,    hxt,    rea,                                   &  
         &      rx,     ry,    rym,                                   & 
         &     hxu,    hyu,   hxyu,   hyxu,                           &
         &     rxu,    ryu,    rxt,    ryt
    use zocmsk,  only :  amskt,  amskv
    use zocnod,  only :  iroot,  myrank
    use zocfil,  only :  ncf
    use ucloc

    implicit none
 
    real(8),   intent(inout)  ::    hx(nxydim),  ubtx(nxydim),  vbtx(nxydim)
    real(8),   intent(in)     ::    hy(nxydim),  ubty(nxydim),  vbty(nxydim)
    real(8),   intent(in)     ::   gxx(nxydim),   gyy(nxydim)
    real(8),   intent(in)     ::  ptop(nxydim),    fw(nxydim)

!----- for tide
    real(8), save             ::    abv(nxydim)
    real(8), save             ::  tidep(nxydim)
    
!----- local variables
    real(8)     ::     cf
    integer(4)  ::     ij,      i,      j
    integer(4)  ::   ijle,   ijln,   ijlw,   ijls
    integer(4)  ::   ijnw,   ijse,   ijsw
    integer(4)  ::  ifpar,  jfpar,   istat
    logical     :: ofirst=.true.

    if (ofirst) then
       ofirst=.false.
       !$acc enter data create(abv, tidep)
       !$acc kernels default(present)
       abv(:)=0.d0
       tidep(:)=0.d0
       !$acc end kernels
    end if
    
!$acc kernels default(present)
!$omp parallel do private( ij )
    do ij = 1, nxydim
       fux(ij) = 0.d0
       fuy(ij) = 0.d0
       fvx(ij) = 0.d0
       fvy(ij) = 0.d0
       fhx(ij) = 0.d0
       fhy(ij) = 0.d0
    end do
!$omp end parallel do
    
!$omp parallel do private( ij )
    do ij = ijtstr-nxdim-1, ijtend+nxdim+2
       fhx(ij) =  - (  ubty(ij+lw) * hyu(ij+lw)                       &
    &                + ubty(ij+lsw) * hyu(ij+lsw)) * 0.5d0
    end do
!$omp end parallel do
!$omp parallel do private( ij )
    do ij = ijtstr-nxdim-1, ijtend+nxdim+nxdim+1
       fhy(ij) = - (  vbty(ij+ls) * hxu(ij+ls)                        &
    &               + vbty(ij+lsw) * hxu(ij+lsw)) * 0.5d0
    end do
!$omp end parallel do
 
!$omp parallel do private( ij )
    do ij = ijtstr-nxdim-1, ijtend+nxdim+1
       hx(ij) = hx(ij)                                                &
    &         + tss * (  (fhx(ij+le) - fhx(ij)) * rx                  &
    &                  + (fhy(ij+ln) - fhy(ij)) * ry(ij)) *           &
    &           rxt(ij) * ryt(ij) * amskt(ij, kstr)                   &
    &         - tss * fw(ij) * amskt(ij, kstr)
    end do
!$omp end parallel do
!$acc end kernels
    if ( otide ) then
       !$acc kernels default(present)
       do ij = ijstr-nxdim-1, ijend+nxdim+1
          abv(ij) = ( bvf(ij   ) + bvf(ij+le ) &
               &    + bvf(ij+ln) + bvf(ij+lne) ) * rghn(ij)
       end do
       !$acc end kernels
       call clcstr('TIDE')
       call tide(tidep)
       call clcend('TIDE')
    end if

    if ( .not. lnovis ) then
!$acc kernels default(present)
!$omp parallel do &
!$omp private( ij, ijls, ijlw, ijln, ijle, ijnw, ijse, ijsw )
    do ij = ijstr-nxdim-1, ijend+nxdim+1
       ijls = ij + ls
       ijlw = ij + lw
       ijln = ij + ln
       ijle = ij + le
       ijnw = ij + lnw
       ijse = ij + lse
       ijsw = ij + lsw
       sxx(ij) = (ubtx(ij) - ubtx(ijlw)) * amhmod(ij) *               &
    &             rx / (hxu(ij) + hxu(ijlw)) * 2.d0                   &
    &          + (vbtx(ij) + vbtx(ijlw)) * amhmod(ij) *               &
    &            (hxyu(ij) + hxyu(ijlw)) * 0.25d0                     &
    &          - (  vbtx(ijln) + vbtx(ijnw)                           &
    &             - vbtx(ijls) - vbtx(ijsw)) * amhmod(ij)             &
    &            / (hyu(ij) + hyu(ijlw)) / (dy(ij) + dy(ijln))        &
    &          - (ubtx(ij) + ubtx(ijlw)) * amhmod(ij) *               &
    &            (hyxu(ij) + hyxu(ijlw)) * 0.25d0                     
       sxy(ij) = (ubtx(ij) - ubtx(ijls)) * amhmod(ij) *               &
    &            ry(ij) / (hyu(ij) + hyu(ijls)) * 2.d0                &
    &          - (ubtx(ij) + ubtx(ijls)) * amhmod(ij) *               & 
    &            (hxyu(ij) + hxyu(ijls)) * 0.25d0                     &
    &          + (  vbtx(ijle) + vbtx(ijse)                           &
    &             - vbtx(ijlw) - vbtx(ijsw)) * amhmod(ij)             &
    &            / (hxu(ij) + hxu(ijls)) * rx * 0.5d0                 &
    &          - (vbtx(ij) + vbtx(ijls)) * amhmod(ij) *               &
    &            (hyxu(ij) + hyxu(ijls)) * 0.25d0                     
       syy(ij) = (vbtx(ij) - vbtx(ijls)) * amhmod(ij) *               &
    &            ry(ij) / (hyu(ij) + hyu(ijls)) * 2.d0                &
    &          + (ubtx(ij) + ubtx(ijls)) * amhmod(ij) *               &
    &            (hyxu(ij) + hyxu(ijls)) * 0.25d0                     &
    &          - (  ubtx(ijle) + ubtx(ijse)                           &
    &             - ubtx(ijlw) - ubtx(ijsw)) * amhmod(ij)             &
    &            / (hxu(ij) + hxu(ijls)) * rx * 0.5d0                 &
    &          - (vbtx(ij) + vbtx(ijls)) * amhmod(ij) *               &
    &            (hxyu(ij) + hxyu(ijls)) * 0.25d0                     
       syx(ij) = (vbtx(ij) - vbtx(ijlw)) * amhmod(ij) *               &
    &            rx / (hxu(ij) + hxu(ijlw)) * 2.d0                    &
    &          - (vbtx(ij) + vbtx(ijlw)) * amhmod(ij) *               &
    &            (hyxu(ij) + hyxu(ijlw)) * 0.25d0                     &
    &          + (  ubtx(ijln) + ubtx(ijnw)                           &
    &             - ubtx(ijls) - ubtx(ijsw)) * amhmod(ij)             &
    &            / (hyu(ij) + hyu(ijlw)) / (dy(ij) + dy(ijln))        &
    &          - (ubtx(ij) + ubtx(ijlw)) * amhmod(ij) *               &
    &            (hxyu(ij) + hxyu(ijlw)) * 0.25d0
    end do
!$omp end parallel do
 
!$omp parallel do private( ij, ijlw )
    do ij = ijvstr-nxdim-1, ijvend+1
       ijlw = ij + lw
       fux(ij) = sxx(ij) * (hyu(ij) + hyu(ijlw)) *                    &
    &                      (hyu(ij) + hyu(ijlw)) * 0.25d0
       fvx(ij) = syx(ij) * (hyu(ij) + hyu(ijlw)) *                    &
    &                      (hyu(ij) + hyu(ijlw)) * 0.25d0
    end do
!$omp end parallel do
 
!$omp parallel do private( ij, ijls )
    do ij = ijvstr-nxdim-1, ijvend+nxdim
       ijls = ij + ls
       fuy(ij) = sxy(ij) * (hxu(ij) + hxu(ijls)) *                    & 
    &                      (hxu(ij) + hxu(ijls)) * 0.25d0
       fvy(ij) = syy(ij) * (hxu(ij) + hxu(ijls)) *                    &
    &                      (hxu(ij) + hxu(ijls)) * 0.25d0
    end do
!$omp end parallel do
!$acc end kernels
    end if ! lnovis

!$acc kernels default(present)
!$omp parallel do private( ij )
    do ij = ijstr-nxdim-1, ijend
       gu(ij) = gxx(ij) + cor(ij) * vbtx(ij)                         &
    &         + (  (fux(ij+le) - fux(ij)) * rx * ryu(ij)             &
    &            + (fuy(ij+ln) - fuy(ij)) * rym(ij) * rxu(ij)) *     &
    &           rxu(ij) * ryu(ij)                                    &
    &         - (  ( (hy(ij+lne)+hy(ij+le))                          &
    &               -(hy(ij+ln )+hy(ij   )) ) * gh(ij)               &
    &            + ( (ptop(ij+lne)+ptop(ij+le))                      &
    &               -(ptop(ij+ln )+ptop(ij   )) )                    &
    &              / rdepv(ij) / rhoo                                &
    &            + ( (tidep(ij+lne)+tidep(ij+le))                    &
    &              - (tidep(ij+ln )+tidep(ij   )) )*alpha/ rdepv(ij) &
    &           ) * 0.5d0 * rx * rxu(ij)                             &
    &         - abv(ij) * (ubtx(ij) - ubtwof(ij)) * rdepv(ij)
       gv(ij) = gyy(ij) - cor(ij) * ubtx(ij)                         &
    &         + (  (fvx(ij+le) - fvx(ij)) * rx * ryu(ij)             &
    &            + (fvy(ij+ln) - fvy(ij)) * rym(ij) * rxu(ij)) *     &
    &           rxu(ij) * ryu(ij)                                    &
    &         - (  ( (hy(ij+lne)+hy(ij+ln))                          &
    &               -(hy(ij+le )+hy(ij   )) ) * gh(ij)               &
    &            + ( (ptop(ij+lne)+ptop(ij+ln))                      &
    &               -(ptop(ij+le )+ptop(ij   )) )                    &
    &              / rdepv(ij) / rhoo                                &
    &            + ( (tidep(ij+lne)+tidep(ij+ln))                    &
    &              - (tidep(ij+le )+tidep(ij   )) )*alpha/ rdepv(ij) &
    &           ) * 0.5d0 * rym(ij) * ryu(ij)                        &
    &         - abv(ij) * (vbtx(ij) - vbtwof(ij)) * rdepv(ij)
    end do
!$omp end parallel do
 
!$omp parallel do private( ij, cf )
    do ij = ijstr-nxdim-1, ijend
       cf = cor(ij) * tss / accb * 0.5d0
       ubtx(ij) = ubtx(ij)                                           &
    &           + tss / accb / (1.d0 + cf * cf) *                    &
    &             (gu(ij) + cf * gv(ij)) * amskv(ij, kstr)
       vbtx(ij) = vbtx(ij)                                           &
    &           + tss / accb / (1.d0 + cf * cf) *                    &
    &             (gv(ij) - cf * gu(ij)) * amskv(ij, kstr)
    end do
!$omp end parallel do
!$acc end kernels
  end subroutine shalow

  subroutine bvfreq(t)
    use zocdim, only : ijtend, ijtstr, kstr, nxdim, nydim
    use zocmsk, only : nbot
    use bshft
    real(8), intent(in) :: t(nxydim, nzdim, ntdim)
    real(8) :: p1, p2, tu, tl, su, sl, ru, rl
    integer :: ij, k, ku
    
    do ij = 1, nxydim
       bvf(ij) = 0.d0
    end do

    do ij = ijtstr, ijtend
       k = nbot(ij)
       if(k .lt. kstr) cycle
       ku = k-1
       tu = t(ij, ku, 1)
       su = t(ij, ku, 2)
       tl = t(ij, k, 1)
       sl = t(ij, k, 2)
       p1 = c0(k) &
       &      + (c1(k) + (c2(k) + c3(k) * tu) * tu) * tu &
       &      + (c4(k) + c5(k) * tu + c6(k) * su) * su
       p2 = d0(k) &
       &      + (d1(k) + (d2(k) + (d3(k) + d4(k) * tu) * tu) * tu) * tu &
       &      + (d5(k) + (d6(k) + d7(k) * tu * tu) * tu &
       &               + (d8(k) + d9(k) * tu * tu) * sqrt(su)) * su
       ru = p1 / p2
       p1 = c0(k) &
       &      + (c1(k) + (c2(k) + c3(k) * tl) * tl) * tl &
       &      + (c4(k) + c5(k) * tl + c6(k) * sl) * sl
       p2 = d0(k) &
       &      + (d1(k) + (d2(k) + (d3(k) + d4(k) * tl) * tl) * tl) * tl &
       &      + (d5(k) + (d6(k) + d7(k) * tl * tl) * tl &
       &               + (d8(k) + d9(k) * tl * tl) * sqrt(sl)) * sl
       rl = p1 / p2
       bvf(ij) = min(bvfmax, max(0.d0, cfb * (rl-ru) * rzm(ij, k)))
       bvf(ij) = sqrt(bvf(ij))
    end do

    call shift1(bvf, nxdim, nydim, 1, 1.d0, 0, 0)

  end subroutine bvfreq

  SUBROUTINE TIDE(TIDEP)

    use zocdim
#ifdef OPT_EXMASK
    use zocgrd,  only :  glont,   glatt
#endif
    use zocnod,  only :  iroot,  myrank
    use ufile
    use ucaln
    
    IMPLICIT NONE
#include "mpif.h"
      real*8 gravc ! gravitational constant [cm^3 s^-2 g^-1]
      parameter(gravc = 6.674D-8)
      real*8 masss, massm ! mass of Sun and Moon [g]
      parameter(masss = 1.989D33, massm = 7.348D25)
      real*8 au ! Astronomical Unit [cm]
      parameter(au = 1.495978707D+13)
      real*8 rearth ! radius of Earth [cm]
      parameter(rearth = 6.371D8)

      REAL*8   TIDEP(NXYDIM)
#ifndef OPT_EXMASK
      REAL*8   GLONT(NXYDIM),  GLATT(NXYDIM) ! dummy
#endif
      
      INTEGER  IFPAR,  JFPAR
      INTEGER      I,     IJ
      REAL*8    ANGM,   ANGS
      REAL*8      JD,     JC
      REAL*8      XX,     YY,     ZZ
      REAL*8    ZETA,    ZED,  THETA
      REAL*8     ALP,    DEL
      REAL*8     GST
      REAL*8       R(6)
      REAL*8     LON(2),    LAT(2),    DIS(2)
      INTEGER  NTARG

      LOGICAL OFIRST
      REAL*8      PI,    PI2,    RAD
      REAL*8      GS,     GM
      REAL*8     SID,    SIH ! seconds in day/hour
      REAL*8  JD2000 ! Julian day for YYYY/MM/DD/hh/mm/ss = 2000/01/01/12/00/00
      REAL*8   T2000 ! Total second in model for 2000/01/01/12/00/00
      INTEGER  J2000(6)
      INTEGER   NCTR
      SAVE OFIRST, PI, PI2, RAD, GS, GM, SID, SIH, J2000, T2000, NCTR

      INTEGER NRFILE,  NRECL,  KSIZE
      CHARACTER CPLEPH*256
      SAVE NRFILE, NRECL, KSIZE, CPLEPH
      NAMELIST /NMJPLD/ NRFILE, NRECL, KSIZE, CPLEPH

      REAL*8     ANG
      REAL*8    LON1,   LAT1,   LON2,   LAT2
      ANG(LON1, LAT1, LON2, LAT2) = &
     &     ACOS(  COS(LAT1) * COS(LON1) * COS(LAT2) * COS(LON2) &
     &     + COS(LAT1) * SIN(LON1) * COS(LAT2) * SIN(LON2) &
     &     + SIN(LAT1) * SIN(LAT2))

      DATA NRFILE / 12 /
      DATA NRECL / 4 /
      DATA KSIZE / 2036 /
      DATA CPLEPH / 'JPLEPH' /

      DATA OFIRST / .TRUE. /
      DATA J2000 / 2000, 1, 1, 12, 0, 0 /
      DATA JD2000 / 2451545.D0 /
      DATA NCTR / 3 / ! For Earth

      IF (OFIRST) THEN
         OFIRST = .FALSE.
         CALL REWNML(IFPAR, JFPAR)
         READ(IFPAR, NMJPLD, END=99)
 99      WRITE(JFPAR, NMJPLD)
         PI = ATAN(1.D0) * 4.D0
         PI2 = 2.D0 * PI
         RAD = PI / 180.D0
         SID = 8.64D4
         SIH = 3.6D3
         GM = - 0.75D0 * gravc * massm * rearth**2 / au**3
         GS = - 0.75D0 * gravc * masss * rearth**2 / au**3
         CALL CYH2SS(T2000, J2000)
      END IF

      JD = (TTSP - T2000) / SID + JD2000
      JC = (TTSP - T2000) / SID / 36525.D0
      GST = 18.697374558D0 & ! Greenwich mean sidereal time
     &     + DMOD(1.00273790935D0 * (JD - JD2000) * 24.D0, 24.D0)
      GST = DMOD(GST * 15.D0 * RAD, PI2)

      ZETA = (  2306.2181D0 * JC &
     &     + 0.30188D0 * JC * JC &
     &     + 0.017998D0 * JC * JC * JC) / SIH * RAD
      ZED = (  2306.2181D0 * JC &
     &     + 1.09468D0 * JC * JC &
     &     + 0.018203D0 * JC * JC * JC) / SIH * RAD
      THETA = (  2004.3109D0 * JC &
     &     - 0.42665D0 * JC * JC &
     &     - 0.041833D0 * JC * JC * JC) / SIH * RAD

      DO I = 1, 2 ! I=1: Moon; I=2: Sun
         NTARG = 9 + I

         IF (MYRANK .EQ. IROOT) THEN
            CALL PLEPH(JD, NTARG, NCTR, R, &
     &                 NRFILE, NRECL, KSIZE, CPLEPH)
         END IF
         CALL MPI_BCAST(R, 6, MPI_REAL8, IROOT, MPI_COMM_WORLD, IERR)

         DIS(I) = SQRT(R(1)*R(1) + R(2)*R(2) + R(3)*R(3))

         XX = (  COS(ZETA) * COS(ZED) * COS(THETA) &
     &         - SIN(ZETA) * SIN(ZED)) * R(1) &
     &      + (- SIN(ZETA) * COS(ZED) * COS(THETA) &
     &         - COS(ZETA) * SIN(ZED)) * R(2) &
     &      - COS(ZED) * SIN(THETA) * R(3)
         YY = (  COS(ZETA) * SIN(ZED) * COS(THETA) &
     &         + SIN(ZETA) * COS(ZED)) * R(1) &
     &      + (- SIN(ZETA) * SIN(ZED) * COS(THETA) &
     &         + COS(ZETA) * COS(ZED)) * R(2) &
     &      - SIN(ZED) * SIN(THETA) * R(3)
         ZZ = COS(ZETA) * SIN(THETA) * R(1) &
     &      - SIN(ZETA) * SIN(THETA) * R(2) &
     &      + COS(THETA) * R(3)

         IF(XX == 0.D0 .AND. YY == 0.D0)THEN
            ALP = 0.D0
         ELSE
            ALP = ATAN2(YY, XX)
         END IF
         LON(I) = DMOD(ALP - GST, PI2)
         LAT(I) = ASIN(ZZ / DIS(I))
      END DO

      DO IJ = IJSTR-NXDIM-1, IJEND+NXDIM+1
         ANGM = ANG(GLONT(IJ), GLATT(IJ), LON(1), LAT(1))
         ANGS = ANG(GLONT(IJ), GLATT(IJ), LON(2), LAT(2))
         TIDEP(IJ) = GM * COS(2.D0 * ANGM) / (DIS(1)**3) &
     &             + GS * COS(2.D0 * ANGS) / (DIS(2)**3)
      END DO

      RETURN
    END SUBROUTINE TIDE

! *********************************************************************
! The following part was taken from testeph.f distributed by JPL
! *********************************************************************
!++++++++++++++++++++++++++

      SUBROUTINE PLEPH ( ET, NTARG, NCENT, RRD , &
           &                 NRFILE, NRECL, KSIZE, NAMFIL)

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IMPLICIT INTEGER (I-N)

      INTEGER NMAX
      PARAMETER (NMAX = 1000)

      INTEGER NRFILE, NRECL, KSIZE
      CHARACTER NAMFIL*256

      DIMENSION RRD(6),ET2Z(2),ET2(2),PV(6,13)
      DIMENSION PVST(6,11),PNUT(4)
      DIMENSION SS(3),CVAL(NMAX),PVSUN(6),ZIPS(2)
      DATA ZIPS/2*0.d0/

      LOGICAL BSAVE,KM,BARY
      LOGICAL FIRST
      DATA FIRST/.TRUE./

      INTEGER LIST(12),IPT(39),DENUM

      COMMON/EPHHDR/CVAL,SS,AU,EMRAT,DENUM,NCON,IPT

      COMMON/STCOMX/KM,BARY,PVSUN

      integer ifpar, jfpar

!C     INITIALIZE ET2 FOR 'STATE' AND SET UP COMPONENT COUNT
!C
      ET2(1)=ET
      ET2(2)=0.D0
      GO TO 11

!C     ENTRY POINT 'DPLEPH' FOR DOUBLY-DIMENSIONED TIME ARGUMENT 
!C          (SEE THE DISCUSSION IN THE SUBROUTINE STATE)

      ENTRY DPLEPH(ET2Z,NTARG,NCENT,RRD)

      ET2(1)=ET2Z(1)
      ET2(2)=ET2Z(2)

  11  IF(FIRST) CALL STATE(ZIPS,LIST,PVST,PNUT, &
     &                     NRFILE,NRECL,KSIZE,NAMFIL)
      FIRST=.FALSE.

  96  IF(NTARG .EQ. NCENT) RETURN

      DO I=1,12
        LIST(I)=0
      ENDDO

!C     CHECK FOR NUTATION CALL

      IF(NTARG.NE.14) GO TO 97
        IF(IPT(35).GT.0) THEN
          LIST(11)=2
          CALL STATE(ET2,LIST,PVST,PNUT, &
     &               NRFILE,NRECL,KSIZE,NAMFIL)
          DO I=1,4
            RRD(I)=PNUT(I)
          ENDDO
          RRD(5) = 0.d0
          RRD(6) = 0.d0
      RETURN
          RETURN
        ELSE
          DO I=1,4
            RRD(I)=0.d0
          ENDDO
          WRITE(6,297)
  297     FORMAT(' *****  NO NUTATIONS ON THE EPHEMERIS FILE  *****')
          STOP
        ENDIF

!C     CHECK FOR LIBRATIONS

  97  CONTINUE
      DO I=1,6
        RRD(I)=0.d0
      ENDDO

      IF(NTARG.NE.15) GO TO 98
        IF(IPT(38).GT.0) THEN
          LIST(12)=2
          CALL STATE(ET2,LIST,PVST,PNUT, &
     &               NRFILE,NRECL,KSIZE,NAMFIL)
          DO I=1,6
            RRD(I)=PVST(I,11)
          ENDDO
      RETURN
          RETURN
        ELSE
          WRITE(6,298)
  298     FORMAT(' *****  NO LIBRATIONS ON THE EPHEMERIS FILE  *****')
          STOP
        ENDIF

!C       FORCE BARYCENTRIC OUTPUT BY 'STATE'

  98  BSAVE=BARY
      BARY=.TRUE.

!C       SET UP PROPER ENTRIES IN 'LIST' ARRAY FOR STATE CALL

      DO I=1,2
        K=NTARG
        IF(I .EQ. 2) K=NCENT
        IF(K .LE. 10) LIST(K)=2
        IF(K .EQ. 10) LIST(3)=2
        IF(K .EQ. 3) LIST(10)=2
        IF(K .EQ. 13) LIST(3)=2
      ENDDO

!C       MAKE CALL TO STATE

      CALL STATE(ET2,LIST,PVST,PNUT, &
     &           NRFILE,NRECL,KSIZE,NAMFIL)

      DO I=1,10
        DO J = 1,6
          PV(J,I) = PVST(J,I)
        ENDDO
      ENDDO

      IF(NTARG .EQ. 11 .OR. NCENT .EQ. 11) THEN
      DO I=1,6
        PV(I,11)=PVSUN(I)
      ENDDO
      ENDIF

      IF(NTARG .EQ. 12 .OR. NCENT .EQ. 12) THEN
        DO I=1,6
          PV(I,12)=0.D0
        ENDDO
      ENDIF

      IF(NTARG .EQ. 13 .OR. NCENT .EQ. 13) THEN
        DO I=1,6
          PV(I,13) = PVST(I,3)
        ENDDO
      ENDIF

      IF(NTARG*NCENT .EQ. 30 .AND. NTARG+NCENT .EQ. 13) THEN
        DO I=1,6
          PV(I,3)=0.D0
        ENDDO
        GO TO 99
      ENDIF

      IF(LIST(3) .EQ. 2) THEN
        DO I=1,6
          PV(I,3)=PVST(I,3)-PVST(I,10)/(1.D0+EMRAT)
        ENDDO
      ENDIF

      IF(LIST(10) .EQ. 2) THEN
        DO I=1,6
          PV(I,10) = PV(I,3)+PVST(I,10)
        ENDDO
      ENDIF

  99  DO I=1,6
        RRD(I)=PV(I,NTARG)-PV(I,NCENT)
      ENDDO

      BARY=BSAVE
      RETURN
    END SUBROUTINE PLEPH

!C+++++++++++++++++++++++++++++++++

      SUBROUTINE INTERP(BUF,T,NCF,NCM,NA,IFL,PV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

      SAVE
      DOUBLE PRECISION BUF(NCF,NCM,*),T(2),PV(NCM,*),PC(18),VC(18)
!C
      DATA NP/2/
      DATA NV/3/
      DATA TWOT/0.D0/
      DATA PC(1),PC(2)/1.D0,0.D0/
      DATA VC(2)/1.D0/
!C
!C       ENTRY POINT. GET CORRECT SUB-INTERVAL NUMBER FOR THIS SET
!C       OF COEFFICIENTS AND THEN GET NORMALIZED CHEBYSHEV TIME
!C       WITHIN THAT SUBINTERVAL.
!C
      DNA=DBLE(NA)
      DT1=DINT(T(1))
      TEMP=DNA*T(1)
      L=IDINT(TEMP-DT1)+1

!C         TC IS THE NORMALIZED CHEBYSHEV TIME (-1 .LE. TC .LE. 1)

      TC=2.D0*(DMOD(TEMP,1.D0)+DT1)-1.D0

!C       CHECK TO SEE WHETHER CHEBYSHEV TIME HAS CHANGED,
!C       AND COMPUTE NEW POLYNOMIAL VALUES IF IT HAS.
!C       (THE ELEMENT PC(2) IS THE VALUE OF T1(TC) AND HENCE
!C       CONTAINS THE VALUE OF TC ON THE PREVIOUS CALL.)

      IF(TC.NE.PC(2)) THEN
        NP=2
        NV=3
        PC(2)=TC
        TWOT=TC+TC
      ENDIF
!C
!C       BE SURE THAT AT LEAST 'NCF' POLYNOMIALS HAVE BEEN EVALUATED
!C       AND ARE STORED IN THE ARRAY 'PC'.
!C
      IF(NP.LT.NCF) THEN
        DO  I=NP+1,NCF
        PC(I)=TWOT*PC(I-1)-PC(I-2)
        END DO 
        NP=NCF
      ENDIF
!C
!C       INTERPOLATE TO GET POSITION FOR EACH COMPONENT
!C
      DO I=1,NCM
      PV(I,1)=0.D0
      DO J=NCF,1,-1
      PV(I,1)=PV(I,1)+PC(J)*BUF(J,I,L)
      END DO
      END DO
      IF(IFL.LE.1) RETURN
!C
!C       IF VELOCITY INTERPOLATION IS WANTED, BE SURE ENOUGH
!C       DERIVATIVE POLYNOMIALS HAVE BEEN GENERATED AND STORED.
!C
      VFAC=(DNA+DNA)/T(2)
      VC(3)=TWOT+TWOT
      IF(NV.LT.NCF) THEN
        DO I=NV+1,NCF
        VC(I)=TWOT*VC(I-1)+PC(I-1)+PC(I-1)-VC(I-2)
        END DO
        NV=NCF
      ENDIF
!C
!C       INTERPOLATE TO GET VELOCITY FOR EACH COMPONENT
!C
      DO I=1,NCM
      PV(I,2)=0.D0
      DO J=NCF,2,-1
      PV(I,2)=PV(I,2)+VC(J)*BUF(J,I,L)
      END DO
      PV(I,2)=PV(I,2)*VFAC
      END DO
!C
      RETURN
    END SUBROUTINE INTERP

!C+++++++++++++++++++++++++

      SUBROUTINE SPLIT(TT,FR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

      DIMENSION FR(2)

      FR(1)=DINT(TT)
      FR(2)=TT-FR(1)
      IF(TT.GE.0.D0 .OR. FR(2).EQ.0.D0) RETURN
      FR(1)=FR(1)-1.D0
      FR(2)=FR(2)+1.D0

      RETURN
    END SUBROUTINE SPLIT

!C++++++++++++++++++++++++++++++++

      SUBROUTINE STATE(ET2,LIST,PV,PNUT, &
           &                 NRFILE,NRECL,KSIZE,NAMFIL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

      SAVE

      INTEGER OLDMAX
      PARAMETER ( OLDMAX = 400)
      INTEGER NMAX
      PARAMETER ( NMAX = 1000)

!      DIMENSION ET2(2),PV(6,11),PNUT(4),T(2),PJD(4),BUF(1500), &
!     & SS(3),CVAL(NMAX),PVSUN(6)
      REAL(8) :: ET2(2),PV(6,11),PNUT(4),T(2),PJD(4),BUF(1500), &
     & SS(3),CVAL(NMAX),PVSUN(6)

      INTEGER LIST(12),IPT(3,13)

      LOGICAL FIRST
      DATA FIRST/.TRUE./

      CHARACTER*6 TTL(14,3),CNAM(NMAX)
      CHARACTER*256, intent(in), optional :: NAMFIL
      INTEGER, intent(in), optional ::  NRFILE, NRECL, KSIZE

      LOGICAL KM,BARY

      INTEGER :: I, J, K, L, NCON, NUMDE
      REAL(8) :: au, emrat
      
      COMMON/EPHHDR/CVAL,SS,AU,EMRAT,NUMDE,NCON,IPT
      COMMON/CHRHDR/CNAM,TTL
      COMMON/STCOMX/KM,BARY,PVSUN

      IF(FIRST) THEN
        FIRST=.FALSE.
        IRECSZ=NRECL*KSIZE
        NCOEFFS=KSIZE/2
        OPEN(NRFILE, &
     &       FILE=NAMFIL, &
     &       ACCESS='DIRECT', &
     &       FORM='UNFORMATTED', &
     &       RECL=IRECSZ, &
     &       STATUS='OLD')
        READ(NRFILE,REC=1)TTL,(CNAM(K),K=1,OLDMAX),SS,NCON,AU,EMRAT, &
     & ((IPT(I,J),I=1,3),J=1,12),NUMDE,(IPT(I,13),I=1,3) &
     & ,(CNAM(L),L=OLDMAX+1,NCON)
        IF(NCON .LE. OLDMAX)THEN
           READ(NRFILE,REC=2)(CVAL(I),I=1,OLDMAX)
        ELSE
           READ(NRFILE,REC=2)(CVAL(I),I=1,NCON)
        ENDIF
        NRL=0
      ENDIF

      IF(ET2(1) .EQ. 0.D0) RETURN

      S=ET2(1)-.5D0
      CALL SPLIT(S,PJD(1))
      CALL SPLIT(ET2(2),PJD(3))
      PJD(1)=PJD(1)+PJD(3)+.5D0
      PJD(2)=PJD(2)+PJD(4)
      CALL SPLIT(PJD(2),PJD(3))
      PJD(1)=PJD(1)+PJD(3)

      IF(PJD(1)+PJD(4).LT.SS(1) .OR. PJD(1)+PJD(4).GT.SS(2)) GO TO 98

      NR=IDINT((PJD(1)-SS(1))/SS(3))+3
      IF(PJD(1).EQ.SS(2)) NR=NR-1

      tmp1 = DBLE(NR-3)*SS(3) + SS(1)
      tmp2 = PJD(1) - tmp1
      T(1) = (tmp2 + PJD(4))/SS(3)

      IF(NR.NE.NRL) THEN
        NRL=NR
        READ(NRFILE,REC=NR,ERR=99)(BUF(K),K=1,NCOEFFS)
      ENDIF

      IF(KM) THEN
         T(2)=SS(3)*86400.D0
         AUFAC=1.D0
      ELSE
         T(2)=SS(3)
         AUFAC=1.D0/AU
      ENDIF

      CALL INTERP(BUF(IPT(1,11)),T,IPT(2,11),3,IPT(3,11),2,PVSUN)
      DO I=1,6
         PVSUN(I)=PVSUN(I)*AUFAC
      ENDDO

      DO I=1,10
         IF(LIST(I).EQ.0) CYCLE
         CALL INTERP(BUF(IPT(1,I)),T,IPT(2,I),3,IPT(3,I), &
     &        LIST(I),PV(1,I))
         DO J=1,6
            IF(I.LE.9 .AND. .NOT.BARY) THEN
               PV(J,I)=PV(J,I)*AUFAC-PVSUN(J)
            ELSE
               PV(J,I)=PV(J,I)*AUFAC
            ENDIF
         ENDDO
      END DO

      IF(LIST(11).GT.0 .AND. IPT(2,12).GT.0) &
     & CALL INTERP(BUF(IPT(1,12)),T,IPT(2,12),2,IPT(3,12), &
     & LIST(11),PNUT)

      IF(LIST(12).GT.0 .AND. IPT(2,13).GT.0) &
     & CALL INTERP(BUF(IPT(1,13)),T,IPT(2,13),3,IPT(3,13), &
     & LIST(12),PV(1,11))

      RETURN

  98  WRITE(*,198)ET2(1)+ET2(2),SS(1),SS(2)
 198  FORMAT(' ***  Requested JED,',f12.2, &
     & ' not within ephemeris limits,',2f12.2,'  ***')

      STOP

   99 WRITE(*,'(2F12.2,A80)')ET2,'ERROR RETURN IN STATE'

      STOP
    END SUBROUTINE STATE

!C+++++++++++++++++++++++++++++

      SUBROUTINE CONST(NAM,VAL,SSS,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

      SAVE

      INTEGER NMAX
      PARAMETER (NMAX = 1000)

      CHARACTER*6 NAM(*),TTL(14,3),CNAM(NMAX)

      DOUBLE PRECISION VAL(*),SSS(3),SS(3),CVAL(NMAX),ZIPS(2)
      DOUBLE PRECISION PVST(6,11),PNUT(4)
      DATA ZIPS/2*0.d0/

      INTEGER IPT(3,13),DENUM,LIST(12)
      logical first
      data first/.true./

      COMMON/EPHHDR/CVAL,SS,AU,EMRAT,DENUM,NCON,IPT
      COMMON/CHRHDR/CNAM,TTL

      IF(FIRST) CALL STATE(ZIPS,LIST,PVST,PNUT)
      first=.false.

      N=NCON

      DO I=1,3
         SSS(I)=SS(I)
      ENDDO

      DO I=1,N
         NAM(I)=CNAM(I)
         VAL(I)=CVAL(I)
      ENDDO

      RETURN
    END SUBROUTINE CONST

end module fshlw


