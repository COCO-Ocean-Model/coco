module tflxt

! --- information -----------------------------------------------------
!
!  HISTORY
!     '03.05.12  H.Hasumi: from COCO3.4
!     '07.04.23  H.Hasumi
!     '08.06.11  H.Hasumi: initial/final processing
!     '08.07.10  H.Hasumi: initial/final processing
!     '12.06.04  Y.Komuro: for COCO5.0
!
! ---------------------------------------------------------------------

  use zocdim, only: &
    &  nxdim,  nzdim,  ntdim, nxydim, nxyzdm, &
    &   kstr,   kend,     kz, &
    &     nx,     ny,     nz, &
    &  ijstr,  ijend, ijtstr, ijtend, &
    &     le,     lw,     ln,     ls, &
    &    lnw,    lse,    lsw,    lww,    lss, &
    &  oinit, ofinal
  use zocgrd, only: &
    &     dy,    dym,     dz,    dzm,    dzv,    dsm, &
    &     dx,     rx,     ry,    rym, &
    &     ts,   zbot, &
    &    hxt,    hxu,    hyt,    hyu,    rxt,    ryt
  use zocmsk, only: &
    &  amskt,  amftx,  amfty,  amftz, &
    &   nbot

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
  use qckot

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
  real(8) :: fharmx(nxydim), fharmy(nxydim),   harm(nxydim)
  real(8) ::  hzbot(nxydim)
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

  real(8), save ::  ahb = 0.0d0

  namelist /nmdifb/ ahb

  if (oinit .or. ofinal) then
     return
  end if

  if (ofirst) then
     ofirst = .false.
     call rewnml(ifpar, jfpar)
     read(ifpar, nmdifb, iostat=istat)
     call cstnml(jfpar, 'flxtrc', 'nmdifb', istat)
     write(jfpar, nmdifb)

     do ij = ijstr-nxdim, ijend+nxdim
        rymlls(ij) = 1.d0 / dym(ij) / dym(ij-nxdim)
        rymrym(ij) = rymlls(ij) / (dym(ij) + dym(ij-nxdim))
     end do
     dx2 = dx * dx
     rx2 = 1.d0 / dx2
  end if

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

  do n = 1, ntdim
     do k = kstr, kend
        do ij = ijtstr-nxdim, ijtend+nxdim+nxdim
           ijlw = ij + lw
           ijls = ij + ls
           fharmx(ij) = ahb * (hyu(ijlw) + hyu(ij+lsw)) &
             &              / (hxt(ij) + hxt(ijlw)) * &
             &          (tx(ij, k, n) - tx(ijlw, k, n)) * rx * &
             &          amskt(ij, k) * amskt(ijlw, k)
           fharmy(ij) = ahb * (hxu(ijls) + hxu(ij+lsw)) &
             &              / (hyt(ij) + hyt(ijls)) * &
             &          (tx(ij, k, n) - tx(ijls, k, n)) * rym(ij) * &
             &          amskt(ij, k) * amskt(ijls, k)
        end do
        do ij = ijtstr-nxdim, ijtend+nxdim
           harm(ij) = (  (fharmx(ij+le) - fharmx(ij)) * rx &
             &         + (fharmy(ij+ln) - fharmy(ij)) * ry(ij)) * &
             &        rxt(ij) * ryt(ij)
        end do
        do ij = ijtstr, ijtend+nxdim
           ijlw = ij + lw
           ijls = ij + ls
           ftx(ij, k, n) = - (harm(ij) - harm(ijlw)) * rx * &
             &             (hyu(ijlw) + hyu(ij+lsw)) &
!             &             / (hxt(ij) + hxt(ijlw)) * amftx(ij, kstr) &
             &             / (hxt(ij) + hxt(ijlw)) * amftx(ij, k)
           fty(ij, k, n) = - (harm(ij) - harm(ijls)) * rym(ij) * &
             &             (hxu(ijls) + hxu(ij+lsw)) &
!             &             / (hyt(ij) + hyt(ijls)) * amfty(ij, kstr) &
             &             / (hyt(ij) + hyt(ijls)) * amfty(ij, k)
        end do
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

  call chekin(wzc, 'WZC', &
     &     'ocean vertical velocity on sigma coordinate', 'cm/s', &
     & nx, ny, nz, nxyzdm, 'OCLVMT')

  do n = 1, ntdim
     do k = kstr+1, kend
        kuu = k - 2
        ku  = k - 1
        kd  = k + 1
        do ij = ijtstr, ijtend
           diffz(ij, k) = ahv(ij, k) * rzm(ij, k) * amftz(ij, k)
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
             &  - wzc(ij, k) * tadv) * amftz(ij, k)
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

           fty(ij, k, n) = fty(ij, k, n) &
             &           - v * tadvy(ij, n)
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
               
           ftx(ij, k, n) = ftx(ij, k, n) &
             &           - u * tadvx(ij, n)
        end do
     end do

  end do

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

  return

end subroutine flxtrc

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
  use qckot
  if (oinit .or. ofinal) then
     return
  end if

  call chekin(   ftx,  'FTX', &
     &            'ocean zonal heat flux', 'degC cm^3/rad/s', &
     &            nx,     ny,     nz, nxyzdm, 'OCLVTX')
  call chekin(   fty,  'FTY', &
     &            'ocean meridional heat flux', 'degC cm^3/rad/s', &
     &            nx,     ny,     nz, nxyzdm, 'OCLVTY')
  call chekin(   ftz,  'FTZ', &
     &            'ocean vertical heat flux', 'degC cm/s', &
     &            nx,     ny,     nz, nxyzdm, 'OCLVTT')
  call chekin(   ftx(1, 1, 2),  'FSX', &
     &            'ocean zonal salt flux', 'psu cm^3/rad/s', &
     &            nx,     ny,     nz, nxyzdm, 'OCLVTX')
  call chekin(   fty(1, 1, 2),  'FSY', &
     &            'ocean meridional salt flux', 'psu cm^3/rad/s', &
     &            nx,     ny,     nz, nxyzdm, 'OCLVTY')
  call chekin(   ftz(1, 1, 2),  'FSZ', &
     &            'ocean vertical salt flux', 'psu cm/s', &
     &            nx,     ny,     nz, nxyzdm, 'OCLVTT')

  return
end subroutine chkftx

end module tflxt
