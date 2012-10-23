module ctnuv
! --- information -----------------------------------------------------
!
!  Estimate the pressure gradient and Coriolis terms of the equation
! of motion.
!
!  HISTORY
!     '03.05.13  H.Hasumi: from COCO3.4
!     '07.04.23  H.Hasumi
!     '07.05.01  H.Hasumi: for McDougall et al. (2003) eq. of state
!     '08.07.10  H.Hasumi: initial/final processing
!     '09.01.15  T.Suzuki: bug fix vertical integration thanks to Hasumi
!     '12.10.23  T.Suzuki: for COCO5.0 
!
! ---------------------------------------------------------------------

  use zocdim, only: &
    &  nxdim, nxydim,  nzdim,  ntdim, &
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
  &    fux,    fuy,   fune,   fuse, &
  &    fvx,    fvy,   fvne,   fvse ) 
  real(8), intent(out)   ::    gxx(nxydim),           gyy(nxydim)
  real(8), intent(inout) ::     gx(nxydim, nzdim),     gy(nxydim, nzdim)
  real(8), intent(inout) ::     xx(nxydim, nzdim),     yy(nxydim, nzdim)
  real(8), intent(in)    ::     uy(nxydim, nzdim),     vy(nxydim, nzdim)
  real(8), intent(in)    ::      r(nxydim, nzdim),     hy(nxydim)
  real(8), intent(in)    ::     ux(nxydim, nzdim),     vx(nxydim, nzdim)
  real(8), intent(in)    ::     fux(nxydim, nzdim),    fvx(nxydim, nzdim)  
  real(8), intent(in)    ::     fuy(nxydim, nzdim),    fvy(nxydim, nzdim)
  real(8), intent(in)    ::    fune(nxydim, nzdim),   fvne(nxydim, nzdim)
  real(8), intent(in)    ::    fuse(nxydim, nzdim),   fvse(nxydim, nzdim)

  real(8), save ::     cof

  real(8) ::     uu(nxydim)       ,     uv(nxydim)
  real(8) ::     vv(nxydim)
  real(8) ::  hvbot(nxydim)
  real(8) ::  dzsig(nxydim, nzdim)

  integer ::     ij,      k
  integer ::  ifpar,  jfpar
  logical, save :: ofirst = .true.

  if (oinit .or. ofinal) then
     return
  end if

  if (ofirst) then
     ofirst = .false.
     cof = -0.25d0 * gravit / rhoo * 1.d-3
  end if

  do ij = ijvstr, ijvend+nxdim+1
     hvbot(ij) = (  (hy(ij) + hy(ij+le)) * dy(ij) &
        &         + (hy(ij+ln) + hy(ij+lne)) * dy(ij+ln)) * &
        &         rym(ij) * 0.25d0 &
        &       + zbot
  end do
  do k = kstr, kstr+kz-1
     do ij = ijvstr, ijvend+nxdim+1
        dzsig(ij, k) = ds(k) * hvbot(ij)
     end do
  end do

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
  do k = kstr+1, kstr+kz-1
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

  do k = kstr+kz, kend
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

! -----------------------------------------------------------------------
  do k = kstr, kstr+kz-1
     do ij = ijvstr, ijvend
         gxx(ij) = gxx(ij) &
     &           + (  (fux(ij+le, k) - fux(ij, k)) * rx &
     &              + (fuy(ij+ln, k) - fuy(ij, k)) * rym(ij)) * &
     &             rxu(ij) * ryu(ij) * dzsig(ij,k) &
     &           + (  fune(ij+lne, k) - fune(ij   , k) &
     &              + fuse(ij+le , k) - fuse(ij+ln, k)) * &
     &             rx * rym(ij) * rxu(ij) * ryu(ij) * dzsig(ij,k)
         gyy(ij) = gyy(ij) &
     &           + (  (fvx(ij+le, k) - fvx(ij, k)) * rx &
     &              + (fvy(ij+ln, k) - fvy(ij, k)) * rym(ij)) * &
     &             rxu(ij) * ryu(ij) * dzsig(ij,k) &
     &           + (  fvne(ij+lne, k) - fvne(ij   , k) &
     &              + fvse(ij+le , k) - fvse(ij+ln, k)) * &
     &             rx * rym(ij) * rxu(ij) * ryu(ij) * dzsig(ij,k)
      enddo
  enddo
  do k = kstr+kz, kend
     do ij = ijvstr, ijvend
         gxx(ij) = gxx(ij) &
     &           + (  (fux(ij+le, k) - fux(ij, k)) * rx &
     &              + (fuy(ij+ln, k) - fuy(ij, k)) * rym(ij)) * &
     &             rxu(ij) * ryu(ij) &
     &           + (  fune(ij+lne, k) - fune(ij   , k) &
     &              + fuse(ij+le , k) - fuse(ij+ln, k)) * &
     &             rx * rym(ij) * rxu(ij) * ryu(ij)
         gyy(ij) = gyy(ij) &
     &           + (  (fvx(ij+le, k) - fvx(ij, k)) * rx &
     &              + (fvy(ij+ln, k) - fvy(ij, k)) * rym(ij)) * &
     &             rxu(ij) * ryu(ij) &
     &           + (  fvne(ij+lne, k) - fvne(ij   , k) &
     &              + fvse(ij+le , k) - fvse(ij+ln, k)) * &
     &             rx * rym(ij) * rxu(ij) * ryu(ij) 
     enddo
  enddo
! -----------------------------------------------------------------------
!      do 300 ij = ijvstr, ijvend+nxdim+1
!         fux (ij, kstr) = fux (ij, kstr) * dzsig(ij, kstr)
!         fuy (ij, kstr) = fuy (ij, kstr) * dzsig(ij, kstr)
!         fune(ij, kstr) = fune(ij, kstr) * dzsig(ij, kstr)
!         fuse(ij, kstr) = fuse(ij, kstr) * dzsig(ij, kstr)
!         fvx (ij, kstr) = fvx (ij, kstr) * dzsig(ij, kstr)
!         fvy (ij, kstr) = fvy (ij, kstr) * dzsig(ij, kstr)
!         fvne(ij, kstr) = fvne(ij, kstr) * dzsig(ij, kstr)
!         fvse(ij, kstr) = fvse(ij, kstr) * dzsig(ij, kstr)
! 300  continue
!      do 330 k = kstr+1, kstr+kz-1
!         do 320 ij = ijvstr, ijvend+nxdim+1
!            fux (ij, kstr) = fux (ij, kstr)
!     &                     + fux (ij, k) * dzsig(ij, k)
!            fuy (ij, kstr) = fuy (ij, kstr)
!     &                     + fuy (ij, k) * dzsig(ij, k)
!            fune(ij, kstr) = fune(ij, kstr)
!     &                     + fune(ij, k) * dzsig(ij, k)
!            fuse(ij, kstr) = fuse(ij, kstr)
!     &                     + fuse(ij, k) * dzsig(ij, k)
!            fvx (ij, kstr) = fvx (ij, kstr)
!     &                     + fvx (ij, k) * dzsig(ij, k)
!            fvy (ij, kstr) = fvy (ij, kstr)
!     &                     + fvy (ij, k) * dzsig(ij, k)
!            fvne(ij, kstr) = fvne(ij, kstr)
!     &                     + fvne(ij, k) * dzsig(ij, k)
!            fvse(ij, kstr) = fvse(ij, kstr)
!     &                     + fvse(ij, k) * dzsig(ij, k)
! 320     continue
! 330  continue
!      do 410 k = kstr+kz, kend
!         do 400 ij = ijvstr, ijvend+nxdim+1
!            fux (ij, kstr) = fux (ij, kstr) + fux (ij, k)
!            fuy (ij, kstr) = fuy (ij, kstr) + fuy (ij, k)
!            fune(ij, kstr) = fune(ij, kstr) + fune(ij, k)
!            fuse(ij, kstr) = fuse(ij, kstr) + fuse(ij, k)
!            fvx (ij, kstr) = fvx (ij, kstr) + fvx (ij, k)
!            fvy (ij, kstr) = fvy (ij, kstr) + fvy (ij, k)
!            fvne(ij, kstr) = fvne(ij, kstr) + fvne(ij, k)
!            fvse(ij, kstr) = fvse(ij, kstr) + fvse(ij, k)
!  400    continue
!  410 continue

  do ij = 1, nxydim
     uu(ij) = 0.d0
     vv(ij) = 0.d0
     uv(ij) = 0.d0
  end do

  do ij = ijvstr, ijvend
         uu(ij) = uy(ij, kstr) * uy(ij, kstr) * &
     &            dzsig(ij, kstr) * amskv(ij, kstr)
         vv(ij) = vy(ij, kstr) * vy(ij, kstr) * &
     &            dzsig(ij, kstr) * amskv(ij, kstr)
         uv(ij) = uy(ij, kstr) * vy(ij, kstr) * &
     &            dzsig(ij, kstr) * amskv(ij, kstr)
  enddo
  do k = kstr+1, kstr+kz-1
     do ij = ijvstr, ijvend
            uu(ij) = uu(ij) &
     &             + uy(ij, k) * uy(ij, k) * &
     &               dzsig(ij, k) * amskv(ij, k)
            vv(ij) = vv(ij) &
     &             + vy(ij, k) * vy(ij, k) * &
     &               dzsig(ij, k) * amskv(ij, k)
            uv(ij) = uv(ij) &
     &             + uy(ij, k) * vy(ij, k) * &
     &               dzsig(ij, k) * amskv(ij, k)
     enddo
  enddo
  do k = kstr+kz, kend
     do  ij = ijvstr, ijvend
            uu(ij) = uu(ij) &
     &             + uy(ij, k) * uy(ij, k) * dzv(ij, k) * amskv(ij, k)
            vv(ij) = vv(ij) &
     &             + vy(ij, k) * vy(ij, k) * dzv(ij, k) * amskv(ij, k)
            uv(ij) = uv(ij) &
     &             + uy(ij, k) * vy(ij, k) * dzv(ij, k) * amskv(ij, k)
     enddo
  enddo
  do ij = ijvstr, ijvend
         gxx(ij) = gxx(ij) + vv(ij) * hyxu(ij) - uv(ij) * hxyu(ij)
         gyy(ij) = gyy(ij) + uu(ij) * hxyu(ij) - uv(ij) * hyxu(ij)
  enddo

!      do 700 ij = ijvstr, ijvend
!         gxx(ij) = gxx(ij)
!     &           + (  (fux(ij+le, kstr) - fux(ij, kstr)) * rx
!     &              + (fuy(ij+ln, kstr) - fuy(ij, kstr)) * rym(ij)) *
!     &             rxu(ij) * ryu(ij)
!     &           + (  fune(ij+lne, kstr) - fune(ij   , kstr)
!     &              + fuse(ij+le , kstr) - fuse(ij+ln, kstr)) *
!     &             rx * rym(ij) * rxu(ij) * ryu(ij)
!     &           + vv(ij) * hyxu(ij) - uv(ij) * hxyu(ij)
!         gyy(ij) = gyy(ij)
!     &           + (  (fvx(ij+le, kstr) - fvx(ij, kstr)) * rx
!     &              + (fvy(ij+ln, kstr) - fvy(ij, kstr)) * rym(ij)) *
!     &             rxu(ij) * ryu(ij)
!     &           + (  fvne(ij+lne, kstr) - fvne(ij   , kstr)
!     &              + fvse(ij+le , kstr) - fvse(ij+ln, kstr)) *
!     &             rx * rym(ij) * rxu(ij) * ryu(ij)
!     &           + uu(ij) * hxyu(ij) - uv(ij) * hyxu(ij)
!  700 continue

  return

end subroutine tnduvd

#ifdef OPT_BBL
! *********************************************************************
subroutine tnduvb( &
  &    gxx,    gyy,     gx,     gy, &
  &     xx,     yy,                 &
  &     uy,     vy,     ty,         &
  &     ux,     vx,                 & 
  &    fux,    fuy,   fune,   fuse, &
  &    fvx,    fvy,   fvne,   fvse )

  use xprst

  implicit none

  real(8), intent(inout) ::     gxx(nxydim)       ,    gyy(nxydim)
  real(8), intent(inout) ::      gx(nxydim, nzdim),     gy(nxydim, nzdim)
  real(8), intent(inout) ::      xx(nxydim, nzdim),     yy(nxydim, nzdim)
  real(8), intent(in)    ::      ux(nxydim, nzdim),     vx(nxydim, nzdim)
  real(8), intent(in)    ::      uy(nxydim, nzdim),     vy(nxydim, nzdim)
  real(8), intent(in)    ::      ty(nxydim, nzdim, ntdim)     


  real(8), intent(in)    ::    fux(nxydim, nzdim),    fvx(nxydim, nzdim)
  real(8), intent(in)    ::    fuy(nxydim, nzdim),    fvy(nxydim, nzdim)
  real(8), intent(in)    ::   fune(nxydim, nzdim),   fvne(nxydim, nzdim)
  real(8), intent(in)    ::   fuse(nxydim, nzdim),   fvse(nxydim, nzdim)

  real(8), save ::  c0b(nxydim), c1b(nxydim), c2b(nxydim), c3b(nxydim)
  real(8), save ::  c4b(nxydim), c5b(nxydim), c6b(nxydim)
  real(8), save ::  d0b(nxydim), d1b(nxydim), d2b(nxydim), d3b(nxydim)
  real(8), save ::  d4b(nxydim), d5b(nxydim), d6b(nxydim), d7b(nxydim)
  real(8), save ::  d8b(nxydim), d9b(nxydim)

  real(8), save ::  dept0(nxydim)
  real(8), save ::  cof
  logical, save ::  ofirst = .true.

  real(8) ::      tl,     sl,   rtmp
  real(8) ::       p,     pe,     pn,    pne
  integer ::      ij,      k,     kk
  integer ::    ijle,   ijln,   ijne

!===== Define statement function 
  real(8) ::    rbbl,     tb,     sb
  real(8) ::  c0, c1, c2, c3, c4, c5, c6
  real(8) ::  d0, d1, d2, d3, d4, d5, d6, d7, d8, d9

  rbbl (tb, sb, c0, c1, c2, c3, c4, c5, c6, &
     &  d0, d1, d2, d3, d4, d5, d6, d7, d8, d9) = &
     &   (c0 + (c1 + (c2 + c3 * tb) * tb) * tb &
     &       + (c4 + c5 * tb + c6 * sb) * sb) &
     & / (d0 + (d1 + (d2 + (d3 + d4 * tb) * tb) * tb) * tb &
     &       + (d5 + (d6 + d7 * tb * tb) * tb &
     &             + (d8 + d9 * tb * tb) * sqrt(sb)) * sb) &
     & - 1.d3
!===== 
  if (oinit .or. ofinal) then
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

  do ij = ijvstr, ijvend

     k = max(nbotv(ij)-1, 1)
     ijle = ij + le
     ijln = ij + ln
     ijne = ij + lne

     tl = ty(ij, k, 1)
     sl = ty(ij, k, 2)
     rtmp = rbbl(tl, sl, &
     &           c0b(ij), c1b(ij), c2b(ij), c3b(ij), &
     &           c4b(ij), c5b(ij), c6b(ij), &
     &           d0b(ij), d1b(ij), d2b(ij), d3b(ij), d4b(ij), &
     &           d5b(ij), d6b(ij), d7b(ij), d8b(ij), d9b(ij))
         p = rtmp * dz0(k) * 0.5d0
     do kk = k+1, nbot(ij)-1
            tl = ty(ij, kk, 1)
            sl = ty(ij, kk, 2)
            rtmp = rbbl(tl, sl, &
     &               c0b(ij), c1b(ij), c2b(ij), c3b(ij), &
     &               c4b(ij), c5b(ij), c6b(ij), &
     &               d0b(ij), d1b(ij), d2b(ij), d3b(ij), d4b(ij), &
     &               d5b(ij), d6b(ij), d7b(ij), d8b(ij), d9b(ij))
            p = p + rtmp * dz0(kk)
     end do
         tl = ty(ij, kend, 1)
         sl = ty(ij, kend, 2)
         rtmp = rbbl(tl, sl, &
     &               c0b(ij), c1b(ij), c2b(ij), c3b(ij), &
     &               c4b(ij), c5b(ij), c6b(ij), &
     &               d0b(ij), d1b(ij), d2b(ij), d3b(ij), d4b(ij), &
     &               d5b(ij), d6b(ij), d7b(ij), d8b(ij), d9b(ij))
         p = p + rtmp * dz0(kend) * 0.5d0

         tl = ty(ijle, k, 1)
         sl = ty(ijle, k, 2)
         rtmp = rbbl(tl, sl, &
     &               c0b(ij), c1b(ij), c2b(ij), c3b(ij), &
     &               c4b(ij), c5b(ij), c6b(ij), &
     &               d0b(ij), d1b(ij), d2b(ij), d3b(ij), d4b(ij), &
     &               d5b(ij), d6b(ij), d7b(ij), d8b(ij), d9b(ij))
         pe = rtmp * dz0(k) * 0.5d0
     do kk = k+1, nbot(ijle)-1
            tl = ty(ijle, kk, 1)
            sl = ty(ijle, kk, 2)
            rtmp = rbbl(tl, sl, &
     &               c0b(ij), c1b(ij), c2b(ij), c3b(ij), &
     &               c4b(ij), c5b(ij), c6b(ij), &
     &               d0b(ij), d1b(ij), d2b(ij), d3b(ij), d4b(ij), &
     &               d5b(ij), d6b(ij), d7b(ij), d8b(ij), d9b(ij))
            pe = pe + rtmp * dz0(kk)
     end do
         tl = ty(ijle, kend, 1)
         sl = ty(ijle, kend, 2)
         rtmp = rbbl(tl, sl, &
     &               c0b(ij), c1b(ij), c2b(ij), c3b(ij), &
     &               c4b(ij), c5b(ij), c6b(ij), &
     &               d0b(ij), d1b(ij), d2b(ij), d3b(ij), d4b(ij), &
     &               d5b(ij), d6b(ij), d7b(ij), d8b(ij), d9b(ij))
         pe = pe + rtmp * dz0(kend) * 0.5d0

         tl = ty(ijln, k, 1)
         sl = ty(ijln, k, 2)
         rtmp = rbbl(tl, sl, &
     &               c0b(ij), c1b(ij), c2b(ij), c3b(ij), &
     &               c4b(ij), c5b(ij), c6b(ij), &
     &               d0b(ij), d1b(ij), d2b(ij), d3b(ij), d4b(ij), &
     &               d5b(ij), d6b(ij), d7b(ij), d8b(ij), d9b(ij))
         pn = rtmp * dz0(k) * 0.5d0
     do kk = k+1, nbot(ijln)-1
            tl = ty(ijln, kk, 1)
            sl = ty(ijln, kk, 2)
            rtmp = rbbl(tl, sl, &
     &               c0b(ij), c1b(ij), c2b(ij), c3b(ij), &
     &               c4b(ij), c5b(ij), c6b(ij), &
     &               d0b(ij), d1b(ij), d2b(ij), d3b(ij), d4b(ij), &
     &               d5b(ij), d6b(ij), d7b(ij), d8b(ij), d9b(ij))
            pn = pn + rtmp * dz0(kk)
     end do
         tl = ty(ijln, kend, 1)
         sl = ty(ijln, kend, 2)
         rtmp = rbbl(tl, sl, &
     &               c0b(ij), c1b(ij), c2b(ij), c3b(ij), &
     &               c4b(ij), c5b(ij), c6b(ij), &
     &               d0b(ij), d1b(ij), d2b(ij), d3b(ij), d4b(ij), &
     &               d5b(ij), d6b(ij), d7b(ij), d8b(ij), d9b(ij))
         pn = pn + rtmp * dz0(kend) * 0.5d0

         tl = ty(ijne, k, 1)
         sl = ty(ijne, k, 2)
         rtmp = rbbl(tl, sl, &
     &               c0b(ij), c1b(ij), c2b(ij), c3b(ij), &
     &               c4b(ij), c5b(ij), c6b(ij), &
     &               d0b(ij), d1b(ij), d2b(ij), d3b(ij), d4b(ij), &
     &               d5b(ij), d6b(ij), d7b(ij), d8b(ij), d9b(ij))
         pne = rtmp * dz0(k) * 0.5d0
     do kk = k+1, nbot(ijne)-1
            tl = ty(ijne, kk, 1)
            sl = ty(ijne, kk, 2)
            rtmp = rbbl(tl, sl, &
     &               c0b(ij), c1b(ij), c2b(ij), c3b(ij), &
     &               c4b(ij), c5b(ij), c6b(ij), &
     &               d0b(ij), d1b(ij), d2b(ij), d3b(ij), d4b(ij), &
     &               d5b(ij), d6b(ij), d7b(ij), d8b(ij), d9b(ij))
            pne = pne + rtmp * dz0(kk)
     end do
         tl = ty(ijne, kend, 1)
         sl = ty(ijne, kend, 2)
         rtmp = rbbl(tl, sl, &
     &               c0b(ij), c1b(ij), c2b(ij), c3b(ij), &
     &               c4b(ij), c5b(ij), c6b(ij), &
     &               d0b(ij), d1b(ij), d2b(ij), d3b(ij), d4b(ij), &
     &               d5b(ij), d6b(ij), d7b(ij), d8b(ij), d9b(ij))
         pne = pne + rtmp * dz0(kend) * 0.5d0

         tl = 0.25d0 * (  ty(ij  , kend, 1) + ty(ijle, kend, 1) &
     &                  + ty(ijln, kend, 1) + ty(ijne, kend, 1))
         sl = 0.25d0 * (  ty(ij  , kend, 2) + ty(ijle, kend, 2) &
     &                  + ty(ijln, kend, 2) + ty(ijne, kend, 2)) 
         rtmp = rbbl(tl, sl, &
     &               c0b(ij), c1b(ij), c2b(ij), c3b(ij), &
     &               c4b(ij), c5b(ij), c6b(ij), &
     &               d0b(ij), d1b(ij), d2b(ij), d3b(ij), d4b(ij), &
     &               d5b(ij), d6b(ij), d7b(ij), d8b(ij), d9b(ij))

         px(ij, kend) = px(ij, k) &
     &                + (  pne + pe - pn - p &
     &                   - rtmp * (  dept0(ijne) + dept0(ijle) &
     &                             - dept0(ijln) - dept0(ij  ))) * &
     &                  cof * rx * rxu(ij) * amskvb(ij)
         gx(ij, kend) = gx(ij, kend) &
     &                + (cor(ij) * vx(ij, kend) + px(ij, kend)) * &
     &                  amskvb(ij)
         xx(ij, kend) = xx(ij, kend) + px(ij, kend) * amskvb(ij) 
         gxx(ij)      = gxx(ij) &
     &                + xx(ij, kend) * dzv(ij, kend) * amskvb(ij)

         py(ij, kend) = py(ij, k) &
     &                + (  pne + pn - pe - p &
     &                   - rtmp * (  dept0(ijne) + dept0(ijln) &
     &                             - dept0(ijle) - dept0(ij  ))) * &
     &                  cof * rym(ij) * ryu(ij) * amskvb(ij)
         gy(ij, kend) = gy(ij, kend) &
     &                - (cor(ij) * ux(ij, kend) - py(ij, kend)) * &
     &                  amskvb(ij)
         yy(ij, kend) = yy(ij, kend) + py(ij, kend) * amskvb(ij) 
         gyy(ij)      = gyy(ij) &
     &                + yy(ij, kend) * dzv(ij, kend) * amskvb(ij)

  end do

  do ij = ijvstr, ijvend
         gxx(ij) = gxx(ij) &
     &           + (  vy(ij, kend) * vy(ij, kend) * hyxu(ij) &
     &              - uy(ij, kend) * vy(ij, kend) * hxyu(ij)) * &
     &             dzv(ij, kend) * amskvb(ij)
         gyy(ij) = gyy(ij) &
     &           + (  uy(ij, kend) * uy(ij, kend) * hxyu(ij) &
     &              - uy(ij, kend) * vy(ij, kend) * hyxu(ij)) * &
     &             dzv(ij, kend) * amskvb(ij)
  end do

  return
end subroutine tnduvb
#endif
end module ctnuv


