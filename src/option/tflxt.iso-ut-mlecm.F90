module tflxt

! --- information -----------------------------------------------------
!
!  HISTORY
!     '03.05.12  H.Hasumi: from COCO3.4
!     '07.04.23  H.Hasumi
!     '07.05.01  H.Hasumi: McDougall et al. (2003) eq. of state
!     '07.09.25  H.Hasumi: arguments of CHEKIN
!     '08.06.11  H.Hasumi: initial/final setup
!     '08.07.10  H.Hasumi: initial/final setup
!     '08.08.??  Y.Komuro: Fox-Kemper et al. (2008) parameterization
!     '09.07.30  Y.Komuro: bug fix (sign of mle fluxes)
!     '10.11.29  Y.Komuro: Coarse-resolution modification applied
!                          following Fox-Kemper et al. (2011, OM).
!     '12.08.01  Y.Komuro: for COCO5.0
!
! ---------------------------------------------------------------------

  use zocdim, only: &
    &  nxdim,  nydim,  nzdim,  ntdim, nxydim, nxyzdm, &
    &   kstr,   kend,     kz, &
    &     nx,     ny,     nz, &
    &  ijstr,  ijend, ijtstr, ijtend, &
    &     le,     lw,     ln,     ls, &
    &    lnw,    lse,    lsw,    lww,    lss, &
    &  oinit, ofinal
  use zocgrd, only: &
    &     dy,    dym,     dz,    dzm,    dzv,     ds,    dsm, &
    &     dx,     rx,     ry,    rym, &
    &     ts,   zbot,    cor, &
    &    hxt,    hxu,    hyt,    hyu,    rxt,    ryt
  use zocmsk, only: &
    &  amskt,  amftx,  amfty,  amftz, &
    &   nbot
  use zocphy, only: &
    & gravit,   rhoo

  implicit none
  private

  real(8), save :: ftx(nxydim, nzdim, ntdim)
  real(8), save :: fty(nxydim, nzdim, ntdim)
  real(8), save :: ftz(nxydim, nzdim, ntdim)
  real(8), save :: tadvx(nxydim, ntdim),  tadvy(nxydim, ntdim)

  public :: flxtrc, chkftx
#ifdef OPT_BBL
  public :: flxtrb
#endif

contains 

subroutine flxtrc( &
  &    adt,  diffz, &
  &     tx,     hx,     ty,     hz, &
  &     uy,     vy,      w,    ahv )

  use ufile

  real(8), intent(out) ::    adt(nxydim, nzdim, ntdim)    
  real(8), intent(out) ::  diffz(nxydim, nzdim)
  real(8), intent(in)  ::     tx(nxydim, nzdim, ntdim)
  real(8), intent(in)  ::     ty(nxydim, nzdim, ntdim)
  real(8), intent(in)  ::     hx(nxydim),     hz(nxydim)
  real(8), intent(in)  ::     uy(nxydim, nzdim),     vy(nxydim, nzdim)
  real(8), intent(in)  ::      w(nxydim, nzdim),    ahv(nxydim, nzdim)

  real(8), save :: rymlls(nxydim), rymrym(nxydim)
  real(8), save ::    dx2,    rx2
  real(8), save ::    eps = 1.d-10
  logical, save :: ofirst = .true.

  real(8) ::    wzc(nxydim, nzdim),    rzm(nxydim, nzdim)
  real(8) ::  hzbot(nxydim)
  real(8) ::  xdzdx(nxydim, nzdim),  ydzdy(nxydim, nzdim)
  real(8) ::  zdzdx(nxydim, nzdim),  zdzdy(nxydim, nzdim)
  real(8) ::  xdtdz(nxydim, nzdim, ntdim),  ydtdz(nxydim, nzdim, ntdim)
  real(8) ::  zdtdx(nxydim, nzdim, ntdim),  zdtdy(nxydim, nzdim, ntdim)
  real(8) ::    tnt(nxydim),    tet(nxydim)
  real(8) ::    tst(nxydim),    ts2(nxydim)
  real(8) ::    twt(nxydim),    tww(nxydim)
  real(8) ::    tnw(nxydim),    tsw(nxydim)
  real(8) ::    tse(nxydim)
  real(8) ::      u,      v
  real(8) ::    cxl,    cyl,    cxd,    cyd,    wut
  real(8) ::     a0,     a1,     a2
  real(8) ::    c00,    c01,    c02,    c10,    c11,    c20
  real(8) ::    d00,    d01,    d02,    d10,    d11,    d20
  integer ::    ij,      k,   ijls,   ijlw,      n
  integer ::    kuu,     ku,     kd
  integer ::  ifpar,  jfpar,  istat

  real(8) ::   tadv
  real(8) ::   vpos,   vneg
  real(8) ::  rzmup,   dzcn,  rzmcn
  real(8) ::  rymup,   dycn,  rymcn
  real(8) ::    tup,    tcn,    tdn
  real(8) ::   tcne,   tcnw,   tcnn,   tcns
  real(8) ::  tdiff,  tcurv
  real(8) ::   tref,  tref1,  tref2
  real(8) ::   tmin,   tmax
  real(8) ::   slp1,   slp2

! for mixed layer eddy parameterization
  real(8) :: psigmx(nxydim, nzdim), psigmy(nxydim, nzdim)
  real(8) ::  xpsiy(nxydim, nzdim),  ypsix(nxydim, nzdim)
  real(8) ::  zpsix(nxydim, nzdim),  zpsiy(nxydim, nzdim)
  real(8) ::   igsy(nxydim, nzdim),   igsx(nxydim, nzdim)

  real(8), save ::    ahh = 0.d0,    ahi = 0.d0,    ahg = 0.d0
  namelist /nmdifh/ ahh
  namelist /nmdifi/ ahi
  namelist /nmdifg/ ahg

  if (oinit .or. ofinal) then
     return
  end if

  if (ofirst) then
     ofirst = .false.
     call rewnml(ifpar, jfpar)
     read(ifpar, nmdifh, iostat=istat)
     call cstnml(jfpar, 'flxtrc', 'nmdifh', istat)
     write(jfpar, nmdifh)
     call rewnml(ifpar, jfpar)
     read(ifpar, nmdifi, iostat=istat)
     call cstnml(jfpar, 'flxtrc', 'nmdifi', istat)
     write(jfpar, nmdifi)
     call rewnml(ifpar, jfpar)
     read(ifpar, nmdifg, iostat=istat)
     call cstnml(jfpar, 'flxtrc', 'nmdifg', istat)
     write(jfpar, nmdifg)

     do ij = ijstr-nxdim, ijend+nxdim
        rymlls(ij) = 1.d0 / dym(ij) / dym(ij-nxdim)
        rymrym(ij) = rymlls(ij) / (dym(ij) + dym(ij-nxdim))
     end do
     dx2 = dx * dx
     rx2 = 1.d0 / dx2
  end if

  call dnsgrd( &
     &  xdzdx,  ydzdy,  zdzdx,  zdzdy, &
     &  xdtdz,  ydtdz,  zdtdx,  zdtdy, &
     &  xpsiy,  ypsix, &
     &  zpsix,  zpsiy, &
     &     ty,     tx,     hz )

  do n = 1, ntdim
     do k = 1, nzdim
        do ij = 1, nxydim
           adt(ij, k, n) = 0.d0
           ftx(ij, k, n) = 0.d0
           fty(ij, k, n) = 0.d0
           ftz(ij, k, n) = 0.d0
        end do
     end do
  end do

  do k = 1, nzdim
     do ij = 1, nxydim
        diffz(ij, k) = 0.d0
     end do
  end do

  do ij = 1, nxydim
     hzbot(ij) = hz(ij) + zbot
  end do
  do k = kstr, kstr+kz-1
     do ij = 1, nxydim
        wzc(ij, k) = w(ij, k) * hzbot(ij)
        rzm(ij, k) = 1.d0 / dsm(k) / hzbot(ij)
     end do
  end do
  do k = kstr+kz, kend
     do ij = 1, nxydim
        wzc(ij, k) = w(ij, k)
        rzm(ij, k) = 1.d0 / dzm(ij, k)
     end do
  end do

  call chekin(wzc, 'WZC', nx, ny, nz, nxyzdm, 'OCN')

  do n = 1, ntdim
     do k = kstr+1, kend
        kuu = k - 2
        ku  = k - 1
        kd  = k + 1
        do ij = ijtstr, ijtend
           diffz(ij, k) = &
              &   (  ahv(ij, k) &
              &    + ahi * (  zdzdx(ij, k) * zdzdx(ij, k) &
              &             + zdzdy(ij, k) * zdzdy(ij, k))) * &
              &   rzm(ij, k) * amftz(ij, k)
           wut = wzc(ij, k) * ts
               
           vpos = 0.5d0 + sign(0.5d0, wut)
           vneg = 0.5d0 - sign(0.5d0, wut)
           rzmup = vpos / dzm(ij, kd) + vneg / dzm(ij, ku)
           dzcn  = vpos * dz(ij, k) + vneg * dz(ij, ku)
           rzmcn = vpos / (dzm(ij, k) + dzm(ij, kd)) &
              &  + vneg / (dzm(ij, k) + dzm(ij, ku))

           tup = vpos * tx(ij, kd, n) + vneg * tx(ij, kuu, n)
           tcn = vpos * tx(ij, k , n) + vneg * tx(ij, ku,  n)
           tdn = vpos * tx(ij, ku, n) + vneg * tx(ij, k ,  n)
           tdiff = tdn - tup
           tcurv = (tdn - tcn) / dzm(ij, k) - (tcn - tup) * rzmup

           a2 = tcurv * rzmcn
           a1 = ( tx(ij, ku, n) - tx(ij, k, n)) / dzm(ij, k) &
              &   + 0.5d0 * (dz(ij, k) - dz(ij, ku)) * a2
           a0 = (  dz(ij, k ) * tx(ij, ku, n) &
              &  + dz(ij, ku) * tx(ij, k , n)) &
              & / dzm(ij, k) * 0.5d0 &
              & - dz(ij, ku) * dz(ij, k) * 0.25d0 * a2
           tadv = (  wut * wut / 3.d0 &
              &    - dz(ij, k) * dz(ij, k) / 12.d0) * a2 &
              & - wut * 0.5d0 * a1 + a0
           tcurv = abs(tcurv) * dzcn

           if (     (tcurv .gt. abs(tdiff)) &
              &      .or. (abs(wut) .lt. eps)) then
              tadv = tcn
           else
              tref = tup + (tcn - tup) / abs(wut) * dzm(ij, k)
              vpos = 0.5d0 + sign(0.5d0, tdiff)
              vneg = 0.5d0 - sign(0.5d0, tdiff)
              tmin = vpos * tcn + vneg * max(tdn, tref)
              tmax = vpos * min(tdn, tref) + vneg * tcn
              tadv = min(max(tadv, tmin), tmax)
           end if

           ftz(ij, k, n) =  &
             & (  diffz(ij, k) * (tx(ij, ku, n) - tx(ij, k, n)) &
             &  - ( ( ahi + ahg ) * &
             &     zdzdx(ij, k) - zpsiy(ij, k) ) * zdtdx(ij, k, n) &
             &  - ( ( ahi + ahg ) * &
             &     zdzdy(ij, k) + zpsix(ij, k) ) * zdtdy(ij, k, n) &
             &  - wzc(ij, k) * tadv ) * amftz(ij, k)
        end do
     end do
  end do

  do n = 1, ntdim
     do k = kstr, kend
        do ij = ijtstr, ijtend+nxdim
           ts2(ij) = tx(ij+lss, k, n) * amskt(ij+lss, k) &
             &         + tx(ij+ls , k, n) * (1.d0 - amskt(ij+lss, k))
           tse(ij) = tx(ij+lse, k, n) * amskt(ij+lse, k) &
             &         + tx(ij+ls , k, n) * (1.d0 - amskt(ij+lse, k))
           tsw(ij) = tx(ij+lsw, k, n) * amskt(ij+lsw, k) &
             &         + tx(ij+ls , k, n) * (1.d0 - amskt(ij+lsw, k))
           tnt(ij) = tx(ij+ln , k, n) * amskt(ij+ln, k) &
             &         + tx(ij    , k, n) * (1.d0 - amskt(ij+ln, k))
           tet(ij) = tx(ij+le , k, n) * amskt(ij+le, k) &
             &         + tx(ij    , k, n) * (1.d0 - amskt(ij+le, k))
           twt(ij) = tx(ij+lw , k, n) * amskt(ij+lw, k) &
             &         + tx(ij    , k, n) * (1.d0 - amskt(ij+lw, k))
        end do
            
        do ij = ijtstr, ijtend+nxdim

           ijls = ij + ls
           u = 0.5d0 * &
             &     (  uy(ij+lsw, k) * dzv(ij+lsw, k) * hyu(ij+lsw) &
             &      + uy(ijls  , k) * dzv(ijls  , k) * hyu(ijls  ))
           v = 0.5d0 * &
             &     (  vy(ij+lsw, k) * dzv(ij+lsw, k) * hxu(ij+lsw) &
             &      + vy(ijls  , k) * dzv(ijls  , k) * hxu(ijls  ))
           cxd = u * ts * 4.d0 &
             &       / (hxt(ij) + hxt(ijls)) / (hyt(ij) + hyt(ijls)) &
             &       / max(eps, amfty(ij, k))
           cyd = v * ts * 4.d0 &
             &       / (hxt(ij) + hxt(ijls)) / (hyt(ij) + hyt(ijls)) &
             &       / max(eps, amfty(ij, k))
               
           vpos = 0.5d0 + sign(0.5d0, cyd)
           vneg = 0.5d0 - sign(0.5d0, cyd)
           rymup = vpos / dym(ij+lss) + vneg / dym(ij)
           rymcn = vpos / (dym(ijls) + dym(ij+lss)) &
              &  + vneg / (dym(ijls) + dym(ij))
           dycn  = vpos * dy(ijls) + vneg * dy(ij)
               
           tup  = vpos * ts2(ij)        + vneg * tnt(ij)
           tcn  = vpos * tx(ijls, k, n) + vneg * tx(ij, k, n)
           tdn  = vpos * tx(ij, k, n)   + vneg * tx(ijls, k, n)
           tcne = vpos * tse(ij)        + vneg * tet(ij)
           tcnw = vpos * tsw(ij)        + vneg * twt(ij)
           tdiff = tdn - tup
           tcurv = (tdn - tcn) * rym(ijls) - (tcn - tup) * rymup

           d02 = tcurv * rymcn
           d00 = (  tx(ij, k, n) * dy(ijls) &
              &       + tx(ijls, k, n) * dy(ij)) * rym(ijls) * 0.5d0 &
              &    - d02 * dy(ij) * dy(ijls) * 0.25d0
           d01 = (tx(ij, k, n) - tx(ijls, k, n)) * rym(ijls) &
              &    - d02 * (dy(ij) - dy(ijls)) * 0.5d0
           d20 = (tcne - 2.d0 * tcn + tcnw) * rx2 * 0.5d0
           d11 = (0.5d0 + sign(0.5d0, u)) * &
              &      (  tsw(ij) - tx(ijls, k, n) &
              &       - twt(ij) + tx(ij, k, n)) * rx * rym(ijls) &
              &     + (0.5d0 - sign(0.5d0, u)) * &
              &       (  tet(ij) - tx(ij, k, n) &
              &        - tse(ij) + tx(ijls, k, n)) * rx * rym(ijls)
           d10 = (0.5d0 + sign(0.5d0, u * v)) * &
              &      (  (tx(ij, k, n) - twt(ij)) * rx &
              &       + d20 * dx &
              &       - d11 * dy(ij) * 0.5d0) &
              &    + (0.5d0 - sign(0.5d0, u * v)) * &
              &      (  (tet(ij) - tx(ij, k, n)) * rx &
              &       - d20 * dx &
              &       - d11 * dy(ij) * 0.5d0)
           tadvy(ij, n) = d00 &
              &     - d01 * 0.5d0 * cyd / 12.d0 &
              &     + d02 * (cyd * cyd / 3.d0 - dy(ij) * dy(ij) / 12.d0) &
              &     - d10 * 0.5d0 * cxd &
              &     + d20 * 0.5d0 * cxd * cxd &
              &     + d11 * 0.25d0 * cxd * cyd
           tcurv = abs(tcurv) * dycn
               
           cxd = cxd * rx
           cyd = cyd * rym(ijls)
               
           d10 = (  (0.5d0 + sign(0.5d0, cxd)) * (tcn - tcnw) &
              &   + (0.5d0 - sign(0.5d0, cxd)) * (tcne - tcn))
           tadvy(ij, n) = tadvy(ij, n) + 0.5d0 * cxd * d10
           if (      (tcurv .gt. abs(tdiff)) &
              & .or. (      (abs(cxd) .lt. eps) &
              &       .and. (abs(cyd) .lt. eps))) then
              tadvy(ij, n) = tcn
           else
              slp1 = (1.d0 + abs(cxd * cyd)) &
                 &       / (abs(cxd) + abs(cyd))
              slp2 = abs(cxd)
              tref1 = tup + (tcn - tup) * slp1
              tref2 = tdn + (tcn - tdn) * slp2
              vpos = 0.5d0 + sign(0.5d0, tdiff)
              vneg = 0.5d0 - sign(0.5d0, tdiff)
              tmin = vpos * tcn + vneg * max(tdn, tref2)
              tmax = vpos * min(tdn, tref1) + vneg * tcn
              tadvy(ij, n) = min(max(tadvy(ij, n), tmin), tmax)
              tadvy(ij, n) = tadvy(ij, n) - 0.5d0 * cxd * d10
           end if

           fty(ij, k, n) = &
             &     (  ( ahh + ahi ) * rym(ijls) &
             &      / ( hyt(ij) + hyt(ijls) ) * &
             &        ( tx(ij, k, n) - tx(ijls, k, n) ) * 2.d0 &
             &      - ( ( ahi - ahg ) &
             &        * ydzdy(ij, k) - ypsix(ij, k) ) * ydtdz(ij, k, n) ) &
             &      * ( hxu(ijls) + hxu(ij+lsw) ) * 0.5d0 * amfty(ij, k) &
             &   - v * tadvy(ij, n)
        end do
     end do
  end do

  do n = 1, ntdim
     do k = kstr, kend
        do ij = ijtstr, ijtend+1
           tww(ij) = tx(ij+lww, k, n) * amskt(ij+lww, k) &
              &        + tx(ij+lw , k, n) * (1.d0 - amskt(ij+lww, k))
           tnw(ij) = tx(ij+lnw, k, n) * amskt(ij+lnw, k) &
              &        + tx(ij+lw , k, n) * (1.d0 - amskt(ij+lnw, k))
           tsw(ij) = tx(ij+lsw, k, n) * amskt(ij+lsw, k) &
              &        + tx(ij+lw , k, n) * (1.d0 - amskt(ij+lsw, k))
           tet(ij) = tx(ij+le , k, n) * amskt(ij+le, k) &
              &        + tx(ij    , k, n) * (1.d0 - amskt(ij+le, k))
           tnt(ij) = tx(ij+ln , k, n) * amskt(ij+ln, k) &
              &        + tx(ij    , k, n) * (1.d0 - amskt(ij+ln, k))
           tst(ij) = tx(ij+ls , k, n) * amskt(ij+ls, k) &
              &        + tx(ij    , k, n) * (1.d0 - amskt(ij+ls, k))
        end do

        do ij = ijtstr, ijtend+1

           ijlw = ij + lw
           ijls = ij + ls
           u   = 0.5d0 * &
              &   (  uy(ijlw  , k) * dzv(ijlw  , k) * hyu(ijlw) &
              &    + uy(ij+lsw, k) * dzv(ij+lsw, k) * hyu(ij+lsw))
           v   = 0.5d0 * &
              &   (  vy(ijlw  , k) * dzv(ijlw  , k) * hxu(ijlw) &
              &    + vy(ij+lsw, k) * dzv(ij+lsw, k) * hxu(ij+lsw))
           cxl = u * ts * 4.d0 &
              &   / (hxt(ij) + hxt(ijlw)) / (hyt(ij) + hyt(ijlw)) &
              &   / max(eps, amftx(ij, k))
           cyl = v * ts * 4.d0 &
              &   / (hxt(ij) + hxt(ijlw)) / (hyt(ij) + hyt(ijlw)) &
              &   / max(eps, amftx(ij, k))

           vpos = 0.5d0 + sign(0.5d0, cxl)
           vneg = 0.5d0 - sign(0.5d0, cxl)

           tup  = vpos * tww(ij)        + vneg * tet(ij)
           tcn  = vpos * tx(ijlw, k, n) + vneg * tx(ij, k, n)
           tdn  = vpos * tx(ij, k, n)   + vneg * tx(ijlw, k, n)
           tcnn = vpos * tnw(ij)        + vneg * tnt(ij)
           tcns = vpos * tsw(ij)        + vneg * tst(ij)
           tdiff = tdn - tup
           tcurv = tdn - 2.d0 * tcn + tup
               
           c10 = (tx(ij, k, n) - tx(ijlw, k, n)) * rx
           c20 = 0.5d0 * tcurv
           c00 = (tx(ij, k, n) + tx(ijlw, k, n)) * 0.5d0 &
              &  - c20 * 0.25d0
           c20 = c20 * rx2
           c11 = (0.5d0 + sign(0.5d0, v)) * &
              &   (  tsw(ij) - tx(ijlw, k, n) &
              &    - tst(ij) + tx(ij, k, n)) * rx * rym(ijls) &
              & + (0.5d0 - sign(0.5d0, v)) * &
              &   (  tnt(ij) - tx(ij, k, n) &
              &    - tnw(ij) + tx(ijlw, k, n)) * rx * rym(ij)
           c02 = (  (tcnn - tcn) * dym(ijls) &
              &   + (tcns - tcn) * dym(ij)) * rymrym(ij)
           c01 = (0.25d0 + sign(0.25d0, u * v)) * &
              &  (  (  (tnw(ij) - tx(ijlw, k, n)) * dym(ijls) &
              &      - (tst(ij) - tx(ij, k, n)) * dym(ij)) * &
              &     rymlls(ij) &
              &   - c02 * (dym(ij) - dym(ijls))) &
              & + (0.25d0 - sign(0.25d0, u * v)) * &
              &   (  (  (tnt(ij) - tx(ij, k, n)) * dym(ijls) &
              &       - (tsw(ij) - tx(ijlw, k, n)) * dym(ij)) * &
              &      rymlls(ij) &
              &    - c02 * (dym(ij) - dym(ijls)))
           tadvx(ij, n) = c00 &
              &  - c10 * 0.5d0 * cxl &
              &  + c20 * (cxl * cxl / 3.d0 - dx2 / 12.d0) &
              &  - c01 * 0.5d0 * cyl &
              &  + c02 * 0.5d0 * cyl * cyl &
              &  + c11 * 0.25d0 * cxl * cyl
               
           cxl = cxl * rx
           cyl = cyl * ry(ij)
               
           c01 =  (0.5d0 + sign(0.5d0, cyl)) * (tcn - tcns) &
              & + (0.5d0 - sign(0.5d0, cyl)) * (tcnn - tcn)
           tadvx(ij, n) = tadvx(ij, n) + 0.5d0 * cyl * c01
           if (      (abs(tcurv) .gt. abs(tdiff)) &
              & .or. (      (abs(cxl) .lt. eps) &
              &       .and. (abs(cyl) .lt. eps))) then
              tadvx(ij, n) = tcn
           else
              slp1 = (1.d0 + abs(cxl * cyl)) &
                 &   / (abs(cxl) + abs(cyl))
              slp2 = abs(cyl)
              tref1 = tup + (tcn - tup) * slp1
              tref2 = tdn + (tcn - tdn) * slp2
              vpos = 0.5d0 + sign(0.5d0, tdiff)
              vneg = 0.5d0 - sign(0.5d0, tdiff)
              tmin = vpos * tcn + vneg * max(tdn, tref2)
              tmax = vpos * min(tdn, tref1) + vneg * tcn
              tadvx(ij, n) = min(max(tadvx(ij, n), tmin), tmax)
              tadvx(ij, n) = tadvx(ij, n) - 0.5d0 * cyl * c01
           end if
               
           ftx(ij, k, n) = &
             &     (  ( ahh + ahi ) * rx &
             &      / ( hxt(ij) + hxt(ijlw) ) &
             &      * ( tx(ij, k, n) - tx(ijlw, k, n) ) * 2.d0 &
             &      - ( ( ahi - ahg ) &
             &      * xdzdx(ij, k) + xpsiy(ij, k) ) * xdtdz(ij, k, n) ) &
             &      * ( hyu(ijlw) + hyu(ij+lsw) ) * 0.5d0 * amftx(ij, k) &
             &   - u * tadvx(ij, n)
        end do
     end do

  end do

  do k = kstr, kend
     do ij=1, nxydim
        igsx(ij, k) = ( ahi - ahg ) * xdzdx(ij, k)
        igsy(ij, k) = ( ahi - ahg ) * ydzdy(ij, k)
     end do
  end do
  call chekin(igsx, 'IGSX', nx, ny, nz, nxyzdm, 'OCN')
  call chekin(igsy, 'IGSY', nx, ny, nz, nxyzdm, 'OCN')

  do n = 1, ntdim
     do k = kstr, kend
        do ij = ijtstr, ijtend
           adt(ij, k, n) = &
              & (  (  (ftx(ij+le, k, n) - ftx(ij, k, n)) * rx &
              &     + (fty(ij+ln, k, n) - fty(ij, k, n)) * ry(ij)) * &
              &    rxt(ij) * ryt(ij) &
              &  + ftz(ij, k, n) - ftz(ij, k+1, n)) / dz(ij, k)
        end do
     end do
  end do

  do k = 1, nzdim
     do ij = 1, nxydim
        psigmx(ij, k) = 0.0d0
        psigmy(ij, k) = 0.0d0
     end do
  end do

  do k = kstr, kend
     do ij = ijtstr, ijtend+nxdim
        psigmx(ij, k) = ahg * ydzdy(ij, k)
        psigmy(ij, k) = - ahg * xdzdx(ij, k)
     end do
  end do

  call chekin(psigmx, 'PSIGMX', nx, ny, nz, nxyzdm, 'OCN')
  call chekin(psigmy, 'PSIGMY', nx, ny, nz, nxyzdm, 'OCN')

  return

end subroutine flxtrc

! *********************************************************************

subroutine dnsgrd( &
  &  xdzdx,  ydzdy,  zdzdx,  zdzdy, &
  &  xdtdz,  ydtdz,  zdtdx,  zdtdy, &
  &  xpsiy,  ypsix, &
  &  zpsix,  zpsiy, &
  &     ty,     tx,     hz )

  use xprst
  use ufile

  real(8), intent(out) ::  xdzdx(nxydim, nzdim),  ydzdy(nxydim, nzdim)
  real(8), intent(out) ::  zdzdx(nxydim, nzdim),  zdzdy(nxydim, nzdim)
  real(8), intent(out) ::  xdtdz(nxydim, nzdim, ntdim)
  real(8), intent(out) ::  ydtdz(nxydim, nzdim, ntdim)
  real(8), intent(out) ::  zdtdx(nxydim, nzdim, ntdim)
  real(8), intent(out) ::  zdtdy(nxydim, nzdim, ntdim)
  real(8), intent(in)  ::     ty(nxydim, nzdim, ntdim)
  real(8), intent(in)  ::     tx(nxydim, nzdim, ntdim)

  real(8), save :: c0(nzdim), c1(nzdim), c2(nzdim)
  real(8), save :: c3(nzdim), c4(nzdim), c5(nzdim), c6(nzdim)
  real(8), save :: d0(nzdim), d1(nzdim), d2(nzdim), d3(nzdim), d4(nzdim)
  real(8), save :: d5(nzdim), d6(nzdim), d7(nzdim), d8(nzdim), d9(nzdim)
  real(8), save :: eps = 1.d-20
  logical, save :: ofirst = .true.

  real(8), save ::  cxpsy(nxydim), cypsx(nxydim) 
  real(8), save ::  czpsx(nxydim), czpsy(nxydim) 
  integer, save ::  kzmin

  real(8) ::      r(nxydim, nzdim)
  real(8) ::   hmld(nxydim), hmld1(nxydim)
  real(8) :: rmavez(nxydim), rmav1(nxydim)
  real(8) ::  dzsig(nxydim, nzdim), dzmsig(nxydim, nzdim)
  real(8) ::     zt(nxydim, nzdim),    ztm(nxydim, nzdim)
  real(8) :: rsigth(nxydim, nzdim)
  real(8) ::    nbv(nxydim),            lf(nxydim)
  real(8) :: xpsiy1(nxydim, nzdim), ypsix1(nxydim, nzdim)
  real(8) :: zpsix1(nxydim, nzdim), zpsiy1(nxydim, nzdim)
  real(8) ::  zmld0(nxydim),         zhmld(nxydim, nzdim)
  real(8) ::  hmldx(nxydim),         hmldy(nxydim)
  integer ::   kmld(nxydim)

  real(8) :: rmavdx(nxydim), rmavdy(nxydim)
!  real(8) ::   muzx(nxydim, nzdim),   muzy(nxydim, nzdim)

  real(8) ::  xpsiy(nxydim, nzdim),  ypsix(nxydim, nzdim)
  real(8) ::  zpsix(nxydim, nzdim),  zpsiy(nxydim, nzdim)
  real(8) ::     hz(nxydim)

  real(8) ::   dtdx(nxydim, nzdim, ntdim),   dtdy(nxydim, nzdim, ntdim)
  real(8) ::  dtfdz(nxydim, nzdim, ntdim)
  real(8) ::     p1,     p2
  real(8) ::     tl,     sl
  real(8) ::     rl,    rlw,    rls,    rlu
  real(8) ::   dzdx,   dzdy
  integer ::     ij,      k,      n
  integer ::  ifpar,  jfpar,  istat

  real(8) ::   muzh,  in2dz, n2min, n2l, hmldt
  real(8) :: rsigdf, rsigbt,  dhmld
  real(8) ::     pi,  omega, cormin

  real(8), save :: slpmax = 1.d-2
  real(8), save ::  cm = 8.0d0,  ce = 0.06d0,  fminlt = 10.0d0
  real(8), save ::  lfmin = 1.0d5,  taumle = 10.0d0,  vscl = 50.0d0
  real(8), save ::  drsig = 0.1d0
  integer, save ::  mz = nz,  mzmin = 1
  integer, save ::  nfltdm = 0,  nfltps = 0,  nfltrm = 0
  logical, save ::  ofltdm = .false.,  ofltps = .false.
  logical, save ::  ofltrm = .false.,  ocoamp = .true.

  namelist /nmslpm/ slpmax
  namelist /nmmlep/ cm, ce, fminlt, mz, ofltdm, nfltdm, &
    &               lfmin, taumle, vscl, ofltps, nfltps, mzmin, &
    &               drsig, ocoamp, ofltrm, nfltrm

  if (ofirst) then
     ofirst = .false.
     call rewnml(ifpar, jfpar)
     read(ifpar, nmslpm, iostat=istat)
     call cstnml(jfpar, 'dnsgrd', 'nmslpm', istat)
     write(jfpar, nmslpm)
     call rewnml(ifpar, jfpar)
     read(ifpar, nmmlep, iostat=istat)
     call cstnml(jfpar, 'dnsgrd', 'nmmlep', istat)
     write(jfpar, nmmlep)

     call secoef( &
        &   c0(kstr), c1(kstr), c2(kstr), c3(kstr), &
        &   c4(kstr), c5(kstr), c6(kstr), &
        &   d0(kstr), d1(kstr), d2(kstr), d3(kstr),  d4(kstr), &
        &   d5(kstr), d6(kstr), d7(kstr), d8(kstr),  d9(kstr))

!    === Mixed layer eddy parameterization
!        from Fox-Kemper and Ferrari(2008) ===
!    *** Coarse-resolution modification is applied.
!        The last row in CXPSY/CYPSX corresponds to ds,
!        and is cancelled out by the row just before it. *** 
     pi = atan( 1.d0 )*4.d0
     omega = 2.d0 * pi / 86400.d0
     cormin = 2.d0 * omega * sin( pi*abs(fminlt)/180.d0 )
     kzmin = mzmin + kstr - 1

     do ij = nxdim+2, nxydim
        if (ocoamp) then
!           dscm = max(dx * 0.5d0 * (hxt(ij) + hxt(ij+lw)), &
!             &        dym(ij+ls) * 0.5d0 * (hyt(ij) + hyt(ij+ls)))
           cxpsy(ij) = ce &
             &       / sqrt( &
             &         max( abs(0.5d0*(cor(ij+lsw) + cor(ij+lw))), &
             &              cormin )**2.0d0 &
             &         + 1.0d0 / (8.64d4 * taumle)**2.0d0 ) &
!             &       * dscm
             &       * dx * 0.5d0 * (hxt(ij) + hxt(ij+lw))
           cypsx(ij) = ce &
             &       / sqrt( &
             &         max( abs(0.5d0*(cor(ij+lsw) + cor(ij+ls))), &
             &              cormin )**2.0d0 &
             &         + 1.0d0 / (8.64d4 * taumle)**2.0d0 ) &
!             &       * dscm
             &       * dym(ij+ls) * 0.5d0 * (hyt(ij) + hyt(ij+ls))        
           czpsy(ij) = ce &
             &       / sqrt( &
             &         max( abs(0.25d0* &
             &              ( cor(ij    ) + cor(ij+lw ) &
             &              + cor(ij+ls ) + cor(ij+lsw))), cormin ) &
             &                                              **2.0d0 &
             &         + 1.0d0 / (8.64d4 * taumle)**2.0d0 ) &
             &         * dx * hxt(ij)
           czpsx(ij) = ce &
             &       / sqrt( &
             &         max( abs(0.25d0* &
             &              ( cor(ij    ) + cor(ij+lw ) &
             &              + cor(ij+ls ) + cor(ij+lsw))), cormin ) &
             &                                              **2.0d0 &
             &         + 1.0d0 / (8.64d4 * taumle)**2.0d0 ) &
             &         * dy(ij) * hyt(ij)
        else
           cxpsy(ij) = ce &
             &       / sqrt( &
             &         max( abs(0.5d0*(cor(ij+lsw) + cor(ij+lw))), &
             &              cormin )**2.0d0 )
           cypsx(ij) = ce &
             &       / sqrt( &
             &         max( abs(0.5d0*(cor(ij+lsw) + cor(ij+ls))), &
             &              cormin )**2.0d0 )
           czpsy(ij) = ce &
             &       / sqrt( &
             &         max( abs(0.25d0* &
             &              ( cor(ij    ) + cor(ij+lw ) &
             &              + cor(ij+ls ) + cor(ij+lsw))), cormin ) &
             &                                              **2.0d0 )
           czpsx(ij) = ce &
             &       / sqrt( &
             &         max( abs(0.25d0* &
             &              ( cor(ij    ) + cor(ij+lw ) &
             &              + cor(ij+ls ) + cor(ij+lsw))), cormin ) &
             &                                              **2.0d0 )
        end if
     end do
  end if

  do n = 1, ntdim
     do k = 1, nzdim
        do ij = 1, nxydim
           dtdx(ij, k, n) = 0.d0
           dtdy(ij, k, n) = 0.d0
           dtfdz(ij, k, n) = 0.d0
        end do
     end do
  end do

  do k = kstr, kend
     do ij = 1, nxydim
        tl = ty(ij, k, 1)
        sl = ty(ij, k, 2)
        p1 = c0(k) &
          &  + (c1(k) + (c2(k) + c3(k) * tl) * tl) * tl &
          &  + (c4(k) + c5(k) * tl + c6(k) * sl) * sl
        p2 = d0(k) &
          &  + (d1(k) + (d2(k) + (d3(k) + d4(k) * tl) * tl) * tl) * tl &
          &  + (d5(k) + (d6(k) + d7(k) * tl * tl) * tl &
          &           + (d8(k) + d9(k) * tl * tl) * sqrt(sl)) * sl
        r(ij, k) = p1 / p2
        p1 = c0(kstr) &
          &  + (c1(kstr) + (c2(kstr) + c3(kstr) * tl) * tl) * tl &
          &  + (c4(kstr) + c5(kstr) * tl + c6(kstr) * sl) * sl
        p2 = d0(kstr) &
          &  + (d1(kstr) + (d2(kstr) + &
          &                (d3(kstr) + d4(kstr) * tl) * tl) * tl) * tl &
          &  + (d5(kstr) + (d6(kstr) + d7(kstr) * tl * tl) * tl &
          &         + (d8(kstr) + d9(kstr) * tl * tl) * sqrt(sl)) * sl
        rsigth(ij, k) = p1 / p2
     end do
  end do

  do k = kstr, kend
     do ij = ijtstr, ijtend+nxdim
        tl = ty(ij, k, 1)
        sl = ty(ij, k, 2)
        p1 = c0(k) &
          &  + (c1(k) + (c2(k) + c3(k) * tl) * tl) * tl &
          &  + (c4(k) + c5(k) * tl + c6(k) * sl) * sl 
        p2 = d0(k) &
          &  + (d1(k) + (d2(k) + (d3(k) + d4(k) * tl) * tl) * tl) * tl &
          &  + (d5(k) + (d6(k) + d7(k) * tl * tl) * tl &
          &           + (d8(k) + d9(k) * tl * tl) * sqrt(sl)) * sl 
        rl = p1 / p2

        tl = ty(ij+lw, k, 1)
        sl = ty(ij+lw, k, 2)
        p1 = c0(k) &
          &  + (c1(k) + (c2(k) + c3(k) * tl) * tl) * tl &
          &  + (c4(k) + c5(k) * tl + c6(k) * sl) * sl
        p2 = d0(k) &
          &  + (d1(k) + (d2(k) + (d3(k) + d4(k) * tl) * tl) * tl) * tl &
          &  + (d5(k) + (d6(k) + d7(k) * tl * tl) * tl &
          &           + (d8(k) + d9(k) * tl * tl) * sqrt(sl)) * sl
        rlw = p1 / p2

        tl = ty(ij+ls, k, 1)
        sl = ty(ij+ls, k, 2)
        p1 = c0(k) &
          &  + (c1(k) + (c2(k) + c3(k) * tl) * tl) * tl &
          &  + (c4(k) + c5(k) * tl + c6(k) * sl) * sl
        p2 = d0(k) &
          &  + (d1(k) + (d2(k) + (d3(k) + d4(k) * tl) * tl) * tl) * tl &
          &  + (d5(k) + (d6(k) + d7(k) * tl * tl) * tl &
          &           + (d8(k) + d9(k) * tl * tl) * sqrt(sl)) * sl
        rls = p1 / p2

        dtdx(ij, k, 1) = (rl - rlw) * rx * 2.d0 &
          &              / (hxt(ij) + hxt(ij+lw)) * &
          &              amskt(ij, k) * amskt(ij+lw, k)
        dtdy(ij, k, 1) = (rl - rls) * rym(ij+ls) * 2.d0 &
          &              / (hyt(ij) + hyt(ij+ls)) * &
          &              amskt(ij, k) * amskt(ij+ls, k)
     end do
  end do

  do k = kstr+1, kend
     do ij = ijtstr-nxdim, ijtend+nxdim
        tl = ty(ij, k, 1)
        sl = ty(ij, k, 2)
        p1 = c0(k) &
           &   + (c1(k) + (c2(k) + c3(k) * tl) * tl) * tl &
           &   + (c4(k) + c5(k) * tl + c6(k) * sl) * sl
        p2 = d0(k) &
           &   + (d1(k) + (d2(k) + (d3(k) + d4(k) * tl) * tl) * tl) * tl &
           &   + (d5(k) + (d6(k) + d7(k) * tl * tl) * tl &
           &            + (d8(k) + d9(k) * tl * tl) * sqrt(sl)) * sl
        rl = p1 / p2

        tl = ty(ij, k-1, 1)
        sl = ty(ij, k-1, 2)
        p1 = c0(k) &
           &   + (c1(k) + (c2(k) + c3(k) * tl) * tl) * tl &
           &   + (c4(k) + c5(k) * tl + c6(k) * sl) * sl
        p2 = d0(k) &
           &   + (d1(k) + (d2(k) + (d3(k) + d4(k) * tl) * tl) * tl) * tl &
           &   + (d5(k) + (d6(k) + d7(k) * tl * tl) * tl &
           &            + (d8(k) + d9(k) * tl * tl) * sqrt(sl)) * sl
        rlu = p1 / p2

        dtfdz(ij, k, 1) = min((rlu - rl) / dzm(ij, k), 0.d0) * &
           &              amftz(ij, k)
     end do
  end do

! === Mixed layer eddy parameterization
!     from Fox-Kemper and Ferrari(2008) ===

  do k = kstr, kstr+kz-1
     do ij = 1, nxydim
!        dzsig (ij, k) = (hz(ij) + zbot) * ds(k)
!        dzmsig(ij, k) = (hz(ij) + zbot) * dsm(k)
        dzsig (ij, k) = zbot * ds(k)
        dzmsig(ij, k) = zbot * dsm(k)
     end do
  end do
  do k = kstr+kz, kend
     do ij = 1, nxydim
        dzsig (ij, k) = dz(ij, k)
        dzmsig(ij, k) = dzm(ij, k)
     end do
  end do
  do k = kend+1, nzdim
     do ij = 1, nxydim
        dzsig (ij, k) = dzsig(ij, kend)
        dzmsig(ij, k) = dzmsig(ij, kend)
     end do
  end do
  
  do ij = 1, nxydim
     zt(ij, kstr) = 0.d0
     ztm(ij, kstr) = 0.5d0 * dzmsig(ij, kstr)
  end do
  do k = kstr, kend
     do ij = 1, nxydim
        zt(ij, k+1) = zt(ij, k) + dzsig(ij, k)
        ztm(ij, k+1) = ztm(ij, k) + dzmsig(ij, k+1)
     end do
  end do

  do k = kstr, kend
     do ij = 1, nxydim
        zhmld(ij, k) = 0.0d0
     end do
  end do
  
  do ij = 1, nxydim
     rmavdx(ij) = 0.0d0
     rmavdy(ij) = 0.0d0
     hmldx(ij) = 0.0d0
     hmldy(ij) = 0.0d0
  end do

! calculating mixed layer depth hmld and mld-averaged buoyancy rmavez
!  do ij = ijtstr-nxdim, ijtend+nxdim
!!  do ij = 1, nxydim
!     in2dz = (-1.0d-3) * gravit / rhoo &
!       &   * dtfdz(ij, kstr+1, 1) * ztm(ij, kstr+1)
!     n2min = (-1.0d-3) * gravit / rhoo &
!       &   * dtfdz(ij, kstr+1, 1)
!     hmld(ij) = ztm(ij, kstr+1)
!     rmavez(ij) = r(ij, kstr) * dzsig(ij, kstr) &
!       &        + 0.5d0 * r(ij, kstr+1) * dzsig(ij, kstr+1)
!     do k = kstr+2, min(nbot(ij), mz)
!        n2l = (-1.0d-3) * gravit / rhoo * dtfdz(ij, k, 1)
!        in2dz = in2dz + n2l * dzmsig(ij, k)
!        n2min = min(n2min, n2l)
!        if (      ( (n2l - n2min) * ztm(ij, k) ) &
!          &  .ge. (cm * in2dz) ) then
!!        if (      ( n2l * ztm(ij, k) ) &
!!          &  .ge. (cm * in2dz) ) then
!           exit
!        end if
!        hmld(ij) = ztm(ij, k)
!        rmavez(ij) = rmavez(ij) + 0.5d0 *  &
!          &        (  r(ij, k-1) * dzsig(ij, k-1) &
!          &         + r(ij, k) * dzsig(ij, k) )
!     end do
!     rmavez(ij) = rmavez(ij) / hmld(ij)
!  end do

!  do ij = ijtstr-nxdim, ijtend+nxdim
  do ij = 1, nxydim
     hmld(ij) = ztm(ij, kstr)
     zhmld(ij, kstr) = ztm(ij, kstr)
     kmld(ij) = kstr
     rmavez(ij) = rsigth(ij, kstr) * ztm(ij, kstr)
     rsigbt = rsigth(ij, kstr) + drsig
     rsigdf = 0.d0
     do k = kstr+1, min(nbot(ij), mz)
        if ( rsigth(ij, k) .ge. rsigbt ) then
           dhmld = dzmsig(ij, k) * (rsigbt - rsigth(ij, k-1)) &
             &   / (rsigth(ij, k) - rsigth(ij, k-1))
           hmld(ij) = hmld(ij) + dhmld
           if (dhmld .le. 0.5d0*dzsig(ij, k-1)) then
              zhmld(ij, k-1) = zhmld(ij, k-1) + dhmld
           else
              zhmld(ij, k-1) = dzsig(ij, k-1)
              zhmld(ij, k) = dhmld - 0.5d0 * dzsig(ij, k-1)
           end if
           rmavez(ij) = rmavez(ij) + &
             &          0.5d0 * (rsigth(ij, k-1)+rsigbt) * dhmld
           rsigdf = drsig
           exit
        end if
        hmld(ij) = hmld(ij) + dzmsig(ij, k)
        zhmld(ij, k-1) = dzsig(ij, k-1)
        zhmld(ij, k) = 0.5d0 * dzsig(ij, k)
        kmld(ij) = k
        rmavez(ij) = rmavez(ij) + 0.5d0 * &
          &        ( rsigth(ij, k-1) * dzsig(ij, k-1) &
          &        + rsigth(ij, k  ) * dzsig(ij, k  ) )
        rsigdf = rsigth(ij, k) - rsigth(ij, kstr)
     end do
     rmavez(ij) = rmavez(ij) / hmld(ij)
!     rmavez(ij) = rsigth(ij, kstr)
     nbv(ij) = sqrt(max(1.0d-3*gravit/rhoo*rsigdf/hmld(ij),0.0d0))
  end do

  do ij = ijtstr-nxdim, ijtend+nxdim
     if (ocoamp) then
        lf(ij) = max(4.0d0*nbv(ij)*hmld(ij)/ &
          &          max( abs( cor(ij   ) + cor(ij+lw ) &
          &                  + cor(ij+ls) + cor(ij+lsw) ), eps), lfmin)
     else
        lf(ij) = 1.0d0
     end if
  end do

  if (ofltdm) then
#ifdef OPT_PARALLEL
     call shift1( &
       &            hmld, &
       &           nxdim,  nydim,      1)
#else
     call stbcsh( &
       &            hmld)
#endif
     do n = 1, nfltdm
!        do ij = ijtstr-nxdim, ijtend+nxdim
        do ij = ijstr, ijend
           if (amskt(ij, kstr) .eq. 1.0d0) then
              hmld1(ij) = ( 4.0d0 * hmld(ij) * amskt(ij, kstr) &
                &         + hmld(ij+ls) * amskt(ij+ls, kstr) &
                &         + hmld(ij+ln) * amskt(ij+ln, kstr) &
                &         + hmld(ij+lw) * amskt(ij+lw, kstr) &
                &         + hmld(ij+le) * amskt(ij+le, kstr) ) &
                &       / ( 4.0d0 * amskt(ij, kstr) &
                &         + amskt(ij+ls, kstr) + amskt(ij+ln, kstr) &
                &         + amskt(ij+lw, kstr) + amskt(ij+le, kstr) )
           else
              hmld1(ij) = hmld(ij)
           end if
        end do
!        do ij = ijtstr-nxdim, ijtend+nxdim
        do ij = ijstr, ijend
           hmld(ij) = hmld1(ij)
        end do
#ifdef OPT_PARALLEL
        call shift1( &
          &            hmld, &
          &           nxdim,  nydim,      1)
#else
        call stbcsh( &
          &            hmld)
#endif
     end do
  end if

  if (ofltrm) then
#ifdef OPT_PARALLEL
     call shift1( &
       &          rmavez, &
       &           nxdim,  nydim,      1)
#else
     call stbcsh( &
       &          rmavez)
#endif
     do n = 1, nfltrm
!        do ij = ijtstr-nxdim, ijtend+nxdim
        do ij = ijstr, ijend
           if (amskt(ij, kstr) .eq. 1.0d0) then
              rmav1(ij) = ( 4.0d0 * rmavez(ij) * amskt(ij, kstr) &
                &         + rmavez(ij+ls) * amskt(ij+ls, kstr) &
                &         + rmavez(ij+ln) * amskt(ij+ln, kstr) &
                &         + rmavez(ij+lw) * amskt(ij+lw, kstr) &
                &         + rmavez(ij+le) * amskt(ij+le, kstr) ) &
                &       / ( 4.0d0 * amskt(ij, kstr) &
                &         + amskt(ij+ls, kstr) + amskt(ij+ln, kstr) &
                &         + amskt(ij+lw, kstr) + amskt(ij+le, kstr) )
           else
              rmav1(ij) = rmavez(ij)
           end if
        end do
!        do ij = ijtstr-nxdim, ijtend+nxdim
        do ij = ijstr, ijend
           rmavez(ij) = rmav1(ij)
        end do
#ifdef OPT_PARALLEL
        call shift1( &
          &          rmavez, &
          &           nxdim,  nydim,      1)
#else
        call stbcsh( &
          &          rmavez)
#endif
     end do
  end if

  do ij = 1, nxydim
     if (kmld(ij) .lt. kzmin) then
        hmld(ij) = 0.0d0
     end if
     zmld0(ij) = kmld(ij) - kstr + 1
  end do
        
  call chekin(   hmld,   'HMLD', nx, ny,  1, nxydim, 'SFC')
  call chekin(  zmld0,   'KMLD', nx, ny,  1, nxydim, 'SFC')
  call chekin( rmavez, 'RMAVEZ', nx, ny,  1, nxydim, 'SFC')
  call chekin(      r,      'R', nx, ny, nz, nxyzdm, 'OCN')
  call chekin(     lf,     'LF', nx, ny,  1, nxydim, 'SFC')

!  do ij = ijtstr, ijtend+nxdim
!     rmavdx(ij) = (-1.0d-3) * gravit / rhoo &
!       &        * rx * 2.d0 / (hxt(ij) + hxt(ij+lw)) &
!       &        * (rmavez(ij) - rmavez(ij+lw)) &
!       &        * amskt(ij, kstr) * amskt(ij+lw, kstr)
!     rmavdy(ij) = (-1.0d-3) * gravit / rhoo &
!       &        * rym(ij+ls) * 2.d0 / (hyt(ij) + hyt(ij+ls)) &
!       &        * (rmavez(ij) - rmavez(ij+ls)) &
!       &        * amskt(ij, kstr) * amskt(ij+ls, kstr)
!  end do

  do k = kstr, mz
     do ij = ijtstr, ijtend+nxdim
!        dhmld = min(zhmld(ij, k), zhmld(ij+lw, k))
        dhmld = max(zhmld(ij, k), zhmld(ij+lw, k))
        rmavdx(ij) = rmavdx(ij) + &
          &          dhmld * (-1.0d-3) * gravit / rhoo &
          &        * rx * 2.d0 / (hxt(ij) + hxt(ij+lw)) &
          &        * (rsigth(ij, k) - rsigth(ij+lw, k))
        hmldx(ij) = hmldx(ij) + dhmld
!        dhmld = min(zhmld(ij, k), zhmld(ij+ls, k))
        dhmld = max(zhmld(ij, k), zhmld(ij+ls, k))
        rmavdy(ij) = rmavdy(ij) + &
          &          dhmld * (-1.0d-3) * gravit / rhoo &
          &        * rym(ij+ls) * 2.d0 / (hyt(ij) + hyt(ij+ls)) &
          &        * (rsigth(ij, k) - rsigth(ij+ls, k))
        hmldy(ij) = hmldy(ij) + dhmld
     end do
  end do

  do ij = ijtstr, ijtend+nxdim
     rmavdx(ij) = rmavdx(ij) / max(hmldx(ij), eps) &
       &        * amskt(ij, kstr) * amskt(ij+lw, kstr)
     rmavdy(ij) = rmavdy(ij) / max(hmldy(ij), eps) &
       &        * amskt(ij, kstr) * amskt(ij+ls, kstr)
  end do

! calculating xpsiy and ypsix
  do k = kstr, kend
     do ij = ijtstr, ijtend+nxdim
!     do ij = nxdim+2, nxydim
        hmldt = max(hmld(ij), hmld(ij+lw))
!        hmldt = min(hmld(ij), hmld(ij+lw))
!        hmldt = 0.5d0 * (hmld(ij) + hmld(ij+lw))
        muzh = ( 1.0d0 - 2.0d0 * &
          &      min(0.5d0*(ztm(ij, k)+ztm(ij+lw, k)) &
          &      / max(hmldt, eps), 1.0d0) &
          &    ) ** 2.0d0
!        muzx(ij, k) = ( 1.0d0 - muzh )
!          &         * ( 1.0d0 + 5.0d0 / 21.0d0 * muzh )
        xpsiy(ij, k) = cxpsy(ij) &
          &          * 2.0d0 / (lf(ij) + lf(ij+lw)) &
          &          * hmldt * hmldt &
          &          * (-1.0d0) * rmavdx(ij) &
          &          * ( 1.0d0 - muzh ) * ( 1.0d0 + 5.0d0 / 21.0d0 * muzh ) &
          &          * amskt(ij, k) * amskt(ij+lw, k)
        xpsiy(ij, k) = sign(1.0d0, xpsiy(ij, k)) * &
          &            min(abs(xpsiy(ij, k)), vscl*dzsig(ij, k))
     end do
  end do
  do k = kstr, kend
     do ij = ijtstr, ijtend+nxdim
!     do ij = nxdim+2, nxydim
        hmldt = max(hmld(ij), hmld(ij+ls))
!        hmldt = min(hmld(ij), hmld(ij+ls))
!        hmldt = 0.5d0 * (hmld(ij) + hmld(ij+ls))
        muzh = ( 1.0d0 - 2.0d0 * &
          &      min(0.5d0*(ztm(ij, k)+ztm(ij+ls, k)) &
          &      / max(hmldt, eps), 1.0d0) &
          &    ) ** 2.0d0
!        muzy(ij, k) = ( 1.0d0 - muzh ) &
!     &              * ( 1.0d0 + 5.0d0 / 21.0d0 * muzh )
        ypsix(ij, k) = cypsx(ij) &
          &          * 2.0d0 / (lf(ij) + lf(ij+ls)) &
          &          * hmldt * hmldt &
          &          * rmavdy(ij) &
          &          * ( 1.0d0 - muzh ) * ( 1.0d0 + 5.0d0 / 21.0d0 * muzh ) &
          &          * amskt(ij, k) * amskt(ij+ls, k)
        ypsix(ij, k) = sign(1.0d0, ypsix(ij, k)) * &
          &            min(abs(ypsix(ij, k)), vscl*dzsig(ij, k))
     end do
  end do

!  do k = kstr+1, kend
!     do ij = ijtstr, ijtend
!        hmldt = hmld(ij)
!        muzh = ( 1.0d0 - 2.0d0 * &
!          &      min(zt(ij, k) / max(hmldt, eps), 1.0d0) &
!          &    ) ** 2.0d0
!!        muzx(ij, k) = ( 1.0d0 - muzh ) &
!!          &         * ( 1.0d0 + 5.0d0 / 21.0d0 * muzh )
!        zpsix(ij, k) = czpsx(ij) &
!          &          / lf(ij) &
!          &          * hmldt * hmldt &
!          &          * 0.5d0 * (rmavdy(ij) + rmavdy(ij+ln)) &
!          &          * ( 1.0d0 - muzh ) * ( 1.0d0 + 5.0d0 / 21.0d0 * muzh ) &
!          &          * amftz(ij, k)
!        zpsiy(ij, k) = czpsy(ij) &
!          &          / lf(ij) &
!          &          * hmldt * hmldt &
!          &          * (-0.5d0) * (rmavdx(ij) + rmavdx(ij+le)) &
!          &          * ( 1.0d0 - muzh ) * ( 1.0d0 + 5.0d0 / 21.0d0 * muzh ) &
!          &          * amftz(ij, k)
!        zpsix(ij, k) = sign(1.0d0, zpsix(ij, k)) * &
!          &            min(abs(zpsix(ij, k)), &
!          &                vscl*0.5d0*(dzsig(ij, k-1)+dzsig(ij, k)))
!        zpsiy(ij, k) = sign(1.0d0, zpsiy(ij, k)) * &
!          &            min(abs(zpsiy(ij, k)), &
!          &                vscl*0.5d0*(dzsig(ij, k-1)+dzsig(ij, k)))
!     end do
!  end do

  if (ofltps) then
#ifdef OPT_PARALLEL
     call shift2( &
       &           xpsiy,  ypsix, &
       &           nxdim,  nydim,  nzdim)
!     call shift2( &
!       &           zpsix,  zpsiy, &
!       &           nxdim,  nydim,  nzdim)
#else
     call stbcav( &
       &           xpsiy,  ypsix)
!     call stbcav(
!       &           zpsix,  zpsiy)
#endif
     do n = 1, nfltps
!        do ij = ijtstr-nxdim, ijtend+nxdim
        do k = kstr, kend
           do ij = ijstr, ijend
              if (amftx(ij, k) .eq. 1.0d0) then
                 xpsiy1(ij, k) = &
                   &           ( 4.0d0 * xpsiy(ij, k) * amftx(ij, k) &
                   &           + xpsiy(ij+ls, k) * amftx(ij+ls, k) &
                   &           + xpsiy(ij+ln, k) * amftx(ij+ln, k) &
                   &           + xpsiy(ij+lw, k) * amftx(ij+lw, k) &
                   &           + xpsiy(ij+le, k) * amftx(ij+le, k) ) &
                   &         / ( 4.0d0 * amftx(ij, k) &
                   &           + amftx(ij+ls, k) + amftx(ij+ln, k) &
                   &           + amftx(ij+lw, k) + amftx(ij+le, k) )
              else
                 xpsiy1(ij, k) = xpsiy(ij, k)
              end if
              if (amfty(ij, k) .eq. 1.0d0) then
                 ypsix1(ij, k) = &
                   &           ( 4.0d0 * ypsix(ij, k) * amfty(ij, k) &
                   &           + ypsix(ij+ls, k) * amfty(ij+ls, k) &
                   &           + ypsix(ij+ln, k) * amfty(ij+ln, k) &
                   &           + ypsix(ij+lw, k) * amfty(ij+lw, k) &
                   &           + ypsix(ij+le, k) * amfty(ij+le, k) ) &
                   &         / ( 4.0d0 * amfty(ij, k) &
                   &           + amfty(ij+ls, k) + amfty(ij+ln, k) &
                   &           + amfty(ij+lw, k) + amfty(ij+le, k) )
              else
                 ypsix1(ij, k) = ypsix(ij, k)
              end if
!              if (amftz(ij, k) .eq. 1.0d0) then
!                 zpsix1(ij, k) = &
!                   &           ( 4.0d0 * zpsix(ij, k) * amftz(ij, k) &
!                   &           + zpsix(ij+ls, k) * amftz(ij+ls, k) &
!                   &           + zpsix(ij+ln, k) * amftz(ij+ln, k) &
!                   &           + zpsix(ij+lw, k) * amftz(ij+lw, k) &
!                   &           + zpsix(ij+le, k) * amftz(ij+le, k) ) &
!                   &         / ( 4.0d0 * amftz(ij, k) &
!                   &           + amftz(ij+ls, k) + amftz(ij+ln, k) &
!                   &           + amftz(ij+lw, k) + amftz(ij+le, k) )
!                 zpsiy1(ij, k) = &
!                   &           ( 4.0d0 * zpsiy(ij, k) * amftz(ij, k) &
!                   &           + zpsiy(ij+ls, k) * amftz(ij+ls, k) &
!                   &           + zpsiy(ij+ln, k) * amftz(ij+ln, k) &
!                   &           + zpsiy(ij+lw, k) * amftz(ij+lw, k) &
!                   &           + zpsiy(ij+le, k) * amftz(ij+le, k) ) &
!                   &         / ( 4.0d0 * amftz(ij, k) &
!                   &           + amftz(ij+ls, k) + amftz(ij+ln, k) &
!                   &           + amftz(ij+lw, k) + amftz(ij+le, k) )
!              else
!                 zpsix1(ij, k) = zpsix(ij, k)
!                 zpsiy1(ij, k) = zpsiy(ij, k)
!              end if
           end do
        end do
        do k = kstr, kend
!           do ij = ijtstr-nxdim, ijtend+nxdim
           do ij = ijstr, ijend
              xpsiy(ij, k) = xpsiy1(ij, k)
              ypsix(ij, k) = ypsix1(ij, k)
!              zpsix(ij, k) = zpsix1(ij, k)
!              zpsiy(ij, k) = zpsiy1(ij, k)
           end do
        end do
#ifdef OPT_PARALLEL
        call shift2( &
          &           xpsiy,  ypsix, &
          &           nxdim,  nydim,  nzdim)
!        call shift2( &
!          &           zpsix,  zpsiy,
!          &           nxdim,  nydim,  nzdim)
#else
        call stbcav( &
          &           xpsiy,  ypsix)
!        call stbcav( &
!          &           zpsix,  zpsiy)
#endif
     end do
  end if
        
  do k = kstr+1, kend
     do ij = ijtstr, ijtend
        zpsix(ij, k) = 0.25d0 * &
          &   ( ( ypsix(ij, k-1) + ypsix(ij+ln, k-1) ) * dz(ij, k-1) &
          &   + ( ypsix(ij, k  ) + ypsix(ij+ln, k  ) ) * dz(ij, k  ) ) &
          &     / dzm(ij, k) * amftz(ij, k)
        zpsiy(ij, k) = 0.25d0 * &
          &   ( ( xpsiy(ij, k-1) + xpsiy(ij+le, k-1) ) * dz(ij, k-1) &
          &   + ( xpsiy(ij, k  ) + xpsiy(ij+le, k  ) ) * dz(ij, k  ) ) &
          &     / dzm(ij, k) * amftz(ij, k)
!        zpsix(ij, k) = sign(1.0d0, zpsix(ij, k)) * &
!          &            min(abs(zpsix(ij, k)), &
!          &                vscl*0.5d0*(dzsig(ij, k-1)+dzsig(ij, k)))
!        zpsiy(ij, k) = sign(1.0d0, zpsiy(ij, k)) * &
!          &            min(abs(zpsiy(ij, k)), &
!          &                vscl*0.5d0*(dzsig(ij, k-1)+dzsig(ij, k)))
     end do
  end do

  call chekin(  cxpsy,  'CXPSY', nx, ny,  1, nxydim, 'SFC')
  call chekin(  cypsx,  'CYPSX', nx, ny,  1, nxydim, 'SFC')
  call chekin(  xpsiy,  'XPSIY', nx, ny, nz, nxyzdm, 'OCN')
  call chekin(  ypsix,  'YPSIX', nx, ny, nz, nxyzdm, 'OCN')
  call chekin(  zpsix,  'ZPSIX', nx, ny, nz, nxyzdm, 'OCN')
  call chekin(  zpsiy,  'ZPSIY', nx, ny, nz, nxyzdm, 'OCN')
  call chekin( rmavdx, 'RMAVDX', nx, ny,  1, nxydim, 'SFC')
  call chekin( rmavdy, 'RMAVDY', nx, ny,  1, nxydim, 'SFC')
!  call chekin(   muzx,   'MUZX', nx, ny, nz, nxyzdm, 'OCN')
!  call chekin(   muzy,   'MUZY', nx, ny, nz, nxyzdm, 'OCN')

  do k = kstr, kend
     do ij = ijtstr, ijtend+nxdim
        dzdx = dtdx(ij, k, 1) * 4.d0 &
           &   / (  dtfdz(ij, k  , 1) + dtfdz(ij+lw, k ,  1) &
           &      + dtfdz(ij, k+1, 1) + dtfdz(ij+lw, k+1, 1) - eps)
        xdzdx(ij, k) = min(slpmax, max(-slpmax, dzdx))
        dzdy = dtdy(ij, k, 1) * 4.d0 * dym(ij+ls) &
           &   / (  (dtfdz(ij, k, 1) + dtfdz(ij, k+1, 1)) * dy(ij) &
           &   + (dtfdz(ij+ls, k, 1) + dtfdz(ij+ls, k+1, 1)) * dy(ij+ls) &
           &      - eps)
        ydzdy(ij, k) = min(slpmax, max(-slpmax, dzdy))
     end do
  end do

  do k = kstr+1, kend
     do ij = ijtstr, ijtend
        dzdx = ((dtdx(ij, k-1, 1)+dtdx(ij+le, k-1, 1))*dz(ij, k-1) &
           &  + (dtdx(ij, k  , 1)+dtdx(ij+le, k  , 1))*dz(ij, k)) * &
           &   0.25d0 / dzm(ij, k) / (dtfdz(ij, k, 1) - eps)
        zdzdx(ij, k) = min(slpmax, max(-slpmax, dzdx)) * &
           &           amftz(ij, k)
        dzdy = &
           &   (  (dtdy(ij, k-1, 1)+dtdy(ij+ln, k-1, 1))*dz(ij, k-1) &
           &    + (dtdy(ij, k  , 1)+dtdy(ij+ln, k  , 1))*dz(ij, k)) * &
           &     0.25d0 / dzm(ij, k) / (dtfdz(ij, k, 1) - eps)
        zdzdy(ij, k) = min(slpmax, max(-slpmax, dzdy)) * &
           &           amftz(ij, k)
     end do
  end do

  do n = 1, ntdim
     do k = 1, nzdim
        do ij = 1, nxydim
           dtdx(ij, k, n) = 0.d0
           dtdy(ij, k, n) = 0.d0
           dtfdz(ij, k, n) = 0.d0
        end do
     end do
  end do

  do n = 1, ntdim
     do k = kstr, kend
        do ij = ijtstr, ijtend+nxdim
           dtdx(ij, k, n) = (tx(ij, k, n) - tx(ij+lw, k, n)) * rx * &
              &           2.d0 / (hxt(ij) + hxt(ij+lw)) * &
              &           amskt(ij, k) * amskt(ij+lw, k) 
           dtdy(ij, k, n) = (tx(ij, k, n) - tx(ij+ls, k, n)) * &
              &           rym(ij+ls) * &
              &           2.d0 / (hyt(ij) + hyt(ij+ls)) * &
              &           amskt(ij, k) * amskt(ij+ls, k)
        end do
     end do
  end do

  do n = 1, ntdim
     do k = kstr+1, kend
        do ij = ijtstr-nxdim, ijtend+nxdim
           dtfdz(ij, k, n) = (tx(ij, k-1, n) - tx(ij, k, n)) &
              &              / dzm(ij, k) * amftz(ij, k)
        end do
     end do
  end do

  do n = 1, ntdim
     do k = kstr, kend
        do ij = ijtstr, ijtend+nxdim
           xdtdz(ij, k, n) = &
              &   (  dtfdz(ij, k, n)   + dtfdz(ij+lw, k, n) &
              &    + dtfdz(ij, k+1, n) + dtfdz(ij+lw, k+1, n)) *  0.25d0
           ydtdz(ij, k, n) = &
              &   (  (dtfdz(ij, k, n)    + dtfdz(ij, k+1, n)) * &
              &      dy(ij) &
              &    + (dtfdz(ij+ls, k, n) + dtfdz(ij+ls, k+1, n)) * &
              &      dy(ij+ls)) * &
              &   rym(ij+ls) * 0.25d0
        end do
     end do
  end do

  do n = 1, ntdim
     do k = kstr+1, kend
        do ij = ijtstr, ijtend
           zdtdx(ij, k, n) = &
              &   (  (dtdx(ij, k-1, n) + dtdx(ij+le, k-1, n)) * &
              &      dz(ij, k-1) &
              &    + (dtdx(ij, k, n) + dtdx(ij+le, k, n)) * &
              &      dz(ij, k)) &
              &   / dzm(ij, k) * 0.25d0 * amftz(ij, k)
           zdtdy(ij, k, n) = &
              &   (  (dtdy(ij, k-1, n) + dtdy(ij+ln, k-1, n)) * &
              &      dz(ij, k-1) &
              &    + (dtdy(ij, k, n) + dtdy(ij+ln, k, n)) * &
              &      dz(ij, k)) &
              &   / dzm(ij, k) * 0.25d0 * amftz(ij, k)
        end do
     end do
  end do

  return

end subroutine dnsgrd
#ifdef OPT_BBL
! *********************************************************************

subroutine flxtrb( &
  &    adt,  diffz, &
  &     tx,     ty,     uy,     vy, &
  &      w,    ahv )

  use ufile

  real(8), intent(out) ::    adt(nxydim, nzdim, ntdim)
  real(8), intent(out) ::  diffz(nxydim, nzdim)
  real(8), intent(in)  ::     tx(nxydim, nzdim, ntdim)
  real(8), intent(in)  ::     ty(nxydim, nzdim, ntdim)
  real(8), intent(in)  ::     uy(nxydim, nzdim),     vy(nxydim, nzdim)
  real(8), intent(in)  ::      w(nxydim, nzdim),    ahv(nxydim, nzdim)

  real(8) ::      u,      v
  integer ::     ij,      k,   ijls,   ijlw,      n
  integer ::  ifpar,  jfpar,  istat

  logical, save :: ofirst = .true.

  real(8), save :: ahhbbl = 0.d0
  namelist /nmbbdh/ ahhbbl

  if (oinit .or. ofinal) then
     return
  end if

  if (ofirst) then
     ofirst = .false.
     call rewnml(ifpar, jfpar)
     read(ifpar, nmbbdh, iostat=istat)
     call cstnml(jfpar, 'flxtrb', 'nmbbdh', istat)
     write(jfpar, nmbbdh)
  end if

  do n = 1, ntdim
     do ij = ijtstr, ijtend
        k = nbot(ij)
        ftz(ij, kend, n) = ftz(ij, k, n)
     end do
  end do

  do n = 1, ntdim
     do ij = ijtstr, ijtend+nxdim
        ijls = ij + ls
        v = 0.5d0 * &
           &   (  vy(ij+lsw, kend) * dzv(ij+lsw, kend) * hxu(ij+lsw) &
           &    + vy(ijls  , kend) * dzv(ijls  , kend) * hxu(ijls  ))
        fty(ij, kend, n) = &
           &     ahhbbl * (tx(ij, kend, n) - tx(ij+ls, kend, n)) * &
           &     rym(ij) * &
           &     (hxu(ij+ls) + hxu(ij+lsw)) / (hyt(ij) + hyt(ij+ls)) * &
           &     amfty(ij, kend) &
           &   - v * tadvy(ij, n)
     end do

     do ij = ijtstr, ijtend+1
        ijlw = ij + lw
        u   = 0.5d0 * &
           &   (  uy(ijlw  , kend) * dzv(ijlw  , kend) * hyu(ijlw) &
           &    + uy(ij+lsw, kend) * dzv(ij+lsw, kend) * hyu(ij+lsw)) 
        ftx(ij, kend, n) = &
           &     ahhbbl * (tx(ij, kend, n) - tx(ij+lw, kend, n)) * rx * &
           &     (hyu(ij+lw) + hyu(ij+lsw)) / (hxt(ij) + hxt(ij+lw)) * &
           &     amftx(ij, kend) &
           &   - u * tadvx(ij, n)
     end do

     do ij = ijtstr, ijtend
        adt(ij, kend, n) = &
           &   (  (  (  ftx(ij+le, kend, n) - ftx(ij, kend, n)) * rx &
           &       + (  fty(ij+ln, kend, n) - fty(ij, kend, n)) * ry(ij)) * &
           &      rxt(ij) * ryt(ij) &
           &    + ftz(ij, kend, n)) / dz(ij, kend)
     end do
  end do

  return
end subroutine flxtrb
#endif
! *********************************************************************

subroutine chkftx

  if (oinit .or. ofinal) then
     return
  end if

  call chekin(   ftx,  'FTX', &
     &            nx,     ny,     nz, nxyzdm, 'OCN')
  call chekin(   fty,  'FTY', &
     &            nx,     ny,     nz, nxyzdm, 'OCN')
  call chekin(   ftz,  'FTZ', &
     &            nx,     ny,     nz, nxyzdm, 'OCN')
  call chekin(   ftx(1, 1, 2),  'FSX', &
     &            nx,     ny,     nz, nxyzdm, 'OCN')
  call chekin(   fty(1, 1, 2),  'FSY', &
     &            nx,     ny,     nz, nxyzdm, 'OCN')
  call chekin(   ftz(1, 1, 2),  'FSZ', &
     &            nx,     ny,     nz, nxyzdm, 'OCN')

  return
end subroutine chkftx

end module tflxt
