module ipmmt

! --- information -----------------------------------------------------
!
!  Elastic-visco-plastic rheology of Hunke and Dukowicz (1997, JPO).
!
!  HISTORY
!     '03.06.05  H.Hasumi: from IcedCOCO3
!     '07.04.23  H.Hasumi
!     '07.05.01  H.Hasumi: sea surface slope term neglected
!     '07.09.25  H.Hasumi: multi-category sea ice
!     '08.06.11  H.Hasumi: initial/final processing
!     '08.07.10  H.Hasumi: initial/final processing
!     '09.01.09  Y.Komuro: bug fix (dimension of SGMXX etc. at FINADD)
!                          masking UX and VX, which are referred to
!                          in calculating ice-ocean stress
!     '09.02.20  Y.Komuro: parallel forward time march
!     '10.04.14  M.Kurogi: (COCO4.4 tripolar code by Dr. Suzuki)
!     '12.07.11  Y.Komuro: for COCO5.0
!
! ---------------------------------------------------------------------

  use zocdim, only: &
    &  nxdim,  nydim,  nzdim, nxydim,   nic,   kstr, &
    & ijtstr, ijtend, ijvstr, ijvend, &
    &     le,     lw,     ln,     ls,    lne,    lsw, &
    &  oinit, ofinal
  use zocgrd, only: &
    &     dx,     dy,     rx,     ry,    rym, &
    &     ts,    its, &
    &    cor, &
    &    hxt,    hyt,   hxyt,   hxyu,   hyxt,   hyxu, &
    &    rxt,    rxu,    ryt,    ryu
  use zocmsk, only: &
    &  amskt,  amskv
  use zocphy, only: &
    &   rhoo,   rhoi, gravit

  implicit none
  private

  public :: pmomnt

contains

subroutine pmomnt( &
  &                   uix,    vix, &
  &                  taux,   tauy, &
  &                    ax,     ay,     az, &
  &                   hix,    hiy,    hiz, &
  &                   hsx,    hsy,    hsz, &
  &                   uiy,    viy,   pice, &
  &                    ux,     vx,     hy,   ptop, &
  &                tauaix, tauaiy, tauaox, tauaoy)

  use dvdif

  real(8), intent(inout) ::    uix(nxydim),           vix(nxydim)
  real(8), intent(out)   ::   taux(nxydim),          tauy(nxydim)
  real(8), intent(in)    ::     ax(nxydim, 0:nic),     ay(nxydim, 0:nic)
  real(8), intent(in)    ::     az(nxydim, 0:nic)
  real(8), intent(in)    ::    hix(nxydim, 0:nic),    hiy(nxydim, 0:nic)
  real(8), intent(in)    ::    hiz(nxydim, 0:nic)
  real(8), intent(in)    ::    hsx(nxydim, 0:nic),    hsy(nxydim, 0:nic)
  real(8), intent(in)    ::    hsz(nxydim, 0:nic)
  real(8), intent(in)    ::    uiy(nxydim),           viy(nxydim)
  real(8), intent(inout) ::   pice(nxydim)
  real(8), intent(in)    ::     ux(nxydim, nzdim),     vx(nxydim, nzdim)
  real(8), intent(in)    ::     hy(nxydim)
  real(8), intent(in)    ::   ptop(nxydim)
  real(8), intent(inout) :: tauaix(nxydim),        tauaiy(nxydim)
  real(8), intent(inout) :: tauaox(nxydim),        tauaoy(nxydim)

  real(8), save ::  sgmxx(nxydim),  sgmyy(nxydim),  sgmxy(nxydim)
  real(8), save ::   rdts

  real(8) ::   zeta(nxydim),    eta(nxydim)
  real(8) ::    emz(nxydim),    epz(nxydim)
  real(8) ::   ecof(nxydim)
  real(8) ::    exx(nxydim),    eyy(nxydim),    exy(nxydim)
  real(8) ::   mice(nxydim),   aice(nxydim)
  real(8) :: avrmsx(nxydim),   avra(nxydim)
  real(8) :: accelu(nxydim), accelv(nxydim)
  real(8) ::   ctau(nxydim),     hh(nxydim)
!      common /work/ zeta, eta, emz, epz, ecof, &
!     &              exx, eyy, exy, &
!     &              mice, aice, avrmsx, avra, &
!     &              accelu, accelv, ctau, hh

  real(8), save ::   caic,   cais,   cioc,   cios
  real(8), save ::   elst(nxydim)
  real(8), save ::   cemz,   cepz,  recc2

  integer ::     ij,      k, isplit, nnsplt
  integer ::   ijle,   ijln,  ijlne
  real(8) ::    rad
  real(8) ::   alph,   beta,   gamm,   delt
  real(8) ::  raabb
  integer ::  ifpar,  jfpar,  istat

  logical, save :: ofirst = .true.,   oeof

  real(8), save ::  thetaa = 0.0d0,  thetao = 25.0d0
  real(8), save ::  cwdrag = 5.0d-3,  elast0 = 0.25d0
  integer, save ::  kglev = 1
  integer, save ::  nsplit = 60
  real(8), save ::  ecc = 2.0d0,  dmin = 2.0d-7,  floss = 17.0d0

  namelist /nmidyn/ ecc, dmin, floss
  namelist /nmiprm/ thetaa, thetao, cwdrag, elast0, kglev
  namelist /nmitsp/ nsplit


  if (oinit) then
#ifdef OPT_TRIPOLE
     call rstadd(sgmxx, oeof, nxdim, nydim, 1, 'SGMXX', 'SFC', &
       &                                        1.d0,  0,  0 )
     call rstadd(sgmyy, oeof, nxdim, nydim, 1, 'SGMYY', 'SFC', &
       &                                        1.d0,  0,  0 )
     call rstadd(sgmxy, oeof, nxdim, nydim, 1, 'SGMXY', 'SFC', &
       &                                        1.d0,  0,  0 )
#else
     call rstadd(sgmxx, oeof, nxdim, nydim, 1, 'SGMXX', 'SFC')
     call rstadd(sgmyy, oeof, nxdim, nydim, 1, 'SGMYY', 'SFC')
     call rstadd(sgmxy, oeof, nxdim, nydim, 1, 'SGMXY', 'SFC')
#endif
     return
  end if

  if (ofinal) then
     call finadd(sgmxx, nxdim, nydim, 1, 'SGMXX', 'SFC')
     call finadd(sgmyy, nxdim, nydim, 1, 'SGMYY', 'SFC')
     call finadd(sgmxy, nxdim, nydim, 1, 'SGMXY', 'SFC')
     return
  end if

  if (ofirst) then
!     ofirst = .false.
     call rewnml(ifpar, jfpar)
     read (ifpar, nmiprm, iostat=istat)
     call cstnml(jfpar, 'pmomnt', 'nmiprm', istat)
     write(jfpar, nmiprm)
     call rewnml(ifpar, jfpar)
     read (ifpar, nmitsp, iostat=istat)
     call cstnml(jfpar, 'pmomnt', 'nmitsp', istat)
     write(jfpar, nmitsp)
     call rewnml(ifpar, jfpar)
     read (ifpar, nmidyn, iostat=istat)
     call cstnml(jfpar, 'pmomnt', 'nmidyn', istat)
     write(jfpar, nmidyn)

     cemz = (1.d0 / ecc / ecc - 1.d0) * 0.5d0
     cepz = (1.d0 / ecc / ecc + 1.d0) * 0.5d0
     recc2 = 0.5d0 / ecc / ecc
!     rdts = dble(nsplit) / dt
     rad = atan(1.d0) * 4.d0 / 180.d0
     kglev = kglev + kstr - 1
     caic = cos(rad * thetaa)
     cais = sin(rad * thetaa)
     cioc = cos(rad * thetao)
     cios = sin(rad * thetao)
     cwdrag = rhoo * cwdrag
     do ij = 1, nxydim
        elst(ij) = 2.d0 * rhoi * elast0 * &
          &        min(dx*dx*hxt(ij)*hxt(ij), &
          &            dy(ij)*dy(ij)*hyt(ij)*hyt(ij)) * &
          &        amskt(ij, kstr)
     end do
  end if
      
  if (its .eq. 2) then
     nnsplt = nsplit
  else
     nnsplt = nsplit / 2
  endif
  rdts = dble(nnsplt) / ts

  do ij = 1, nxydim
     zeta  (ij) = 0.d0
     eta   (ij) = 0.d0
     emz   (ij) = 0.d0
     epz   (ij) = 0.d0
     exx   (ij) = 0.d0
     eyy   (ij) = 0.d0
     exy   (ij) = 0.d0
     avrmsx(ij) = 0.d0
     avra  (ij) = 0.d0
     accelu(ij) = 0.d0
     accelv(ij) = 0.d0
     ctau  (ij) = 0.d0
     ecof  (ij) = 0.d0
     aice  (ij) = 1.d0 - ax(ij, 0)
     mice  (ij) = 0.d0
     tauaix(ij) = tauaix(ij) * aice(ij)
     tauaiy(ij) = tauaiy(ij) * aice(ij)
     tauaox(ij) = tauaox(ij) * ax(ij, 0)
     tauaoy(ij) = tauaoy(ij) * ax(ij, 0)
!     hh    (ij) = hy(ij) &
!       &        + ay(ij) * (rhoi * hiy(ij) + rhos * hsy(ij)) / rhoo
     hh    (ij) = 0.d0
  end do
  do k = 1, nic
     do ij = 1, nxydim
        mice(ij) = mice(ij) + ax(ij, k) * hix(ij, k)
     end do
  end do

!  do ij = 1, nxydim
!     uix(ij) = uiy(ij)
!     vix(ij) = viy(ij)
!  end do
  call strain( &
    &             exx,    eyy,    exy, &
    &             uix,    vix)
  call rheolo( &
    &            zeta,    eta,    emz,    epz,   pice, &
    &            aice,   mice,    exx,    eyy,    exy )

  if (ofirst) then
     ofirst = .false.
     if (oeof) then
        do ij = 1, nxydim
           sgmxx(ij) = 2.d0 * eta(ij) * exx(ij) &
             &       - emz(ij) * (exx(ij) + eyy(ij)) &
             &       - 0.5d0 * pice(ij)
           sgmyy(ij) = 2.d0 * eta(ij) * eyy(ij) &
             &       - emz(ij) * (exx(ij) + eyy(ij)) &
             &       - 0.5d0 * pice(ij)
           sgmxy(ij) = 2.d0 * eta(ij) * exy(ij)
        end do
     end if
  end if

  do ij = ijvstr, ijvend
     avrmsx(ij) = (  mice(ij) + mice(ij+le) &
       &           + mice(ij+ln) + mice(ij+lne)) * 0.25d0 * rhoi * &
       &          amskv(ij, kstr)
     avra(ij)   = (  aice(ij) + aice(ij+le) &
       &           + aice(ij+ln) + aice(ij+lne)) * 0.25d0 
     accelu(ij) = (  caic * tauaix(ij) &
       &           - sign(cais, cor(ij)) * tauaiy(ij) &
       &           - avrmsx(ij) * rx * rxu(ij) * 0.5d0 * &
       &             (  hh(ij+lne) + hh(ij+le) &
       &              - hh(ij+ln) - hh(ij)) * gravit &
       &          ) * amskv(ij, kstr)
     accelv(ij) = (  caic * tauaiy(ij) &
       &           + sign(cais, cor(ij)) * tauaix(ij) &
       &           - avrmsx(ij) * rym(ij) * ryu(ij) * 0.5d0 * &
       &             (  hh(ij+lne) + hh(ij+ln) &
       &              - hh(ij+le) - hh(ij)) * gravit &
       &          ) * amskv(ij, kstr)
  end do

!  do isplit = 1, nsplit
  do isplit = 1, nnsplt

     do ij = ijtstr, ijtend
        if (aice(ij) .eq. 0.d0) then
           sgmxx(ij) = 0.d0
           sgmyy(ij) = 0.d0
           sgmxy(ij) = 0.d0
        else
           ecof(ij) = elst(ij) * mice(ij) * rdts * rdts
           alph = 2.d0 * eta(ij) * rdts + cepz * ecof(ij)
           beta = cemz * ecof(ij)
           gamm = 2.d0 * eta(ij) * &
             &    (sgmxx(ij) * rdts + ecof(ij) * exx(ij)) &
             &  - ecof(ij) * pice(ij) * recc2
           delt = 2.d0 * eta(ij) * &
             &    (sgmyy(ij) * rdts + ecof(ij) * eyy(ij)) &
             &  - ecof(ij) * pice(ij) * recc2
           raabb = 1.d0 / (alph * alph - beta * beta)
           sgmxx(ij) = (alph * gamm - beta * delt) * raabb
           sgmyy(ij) = (alph * delt - beta * gamm) * raabb
           sgmxy(ij) = 2.d0 * eta(ij) * &
             &         (sgmxy(ij) * rdts + ecof(ij) * exy(ij)) &
             &       / (2.d0 * eta(ij) * rdts + ecof(ij))
        end if
     end do

#ifdef OPT_PARALLEL
#ifdef OPT_TRIPOLE
     call shift3( sgmxx,  sgmyy,  sgmxy, &
       &          nxdim,  nydim,      1, &
       &           1.d0,      0,      0 )
#else
     call shift3( &
       &          sgmxx,  sgmyy,  sgmxy, &
       &          nxdim,  nydim,      1)
#endif
#endif

     do ij = ijvstr, ijvend
        if (avrmsx(ij) .eq. 0.d0) then
           uix(ij) = 0.d0
           vix(ij) = 0.d0
        else
           ijle = ij + le
           ijln = ij + ln
           ijlne = ij + lne
           ctau(ij) = avra(ij) * cwdrag * &
             &        sqrt(  (uix(ij) &
             &               - ux(ij, kglev)*amskv(ij, kglev) )**2 &
             &             + (vix(ij) &
             &               - vx(ij, kglev)*amskv(ij, kglev) )**2)
           alph = avrmsx(ij) * rdts &
             &  + ctau(ij) * cioc
           beta = avrmsx(ij) * cor(ij) &
             &  + ctau(ij) * sign(cios, cor(ij))
           gamm = (  (  sgmxx(ijlne) * hyt(ijlne) &
             &        + sgmxx(ijle ) * hyt(ijle ) &
             &        - sgmxx(ijln ) * hyt(ijln ) &
             &        - sgmxx(ij   ) * hyt(ij   )) * rx &
             &     + (  sgmxy(ijlne) * hxt(ijlne) &
             &        + sgmxy(ijln ) * hxt(ijln ) &
             &        - sgmxy(ijle ) * hxt(ijle ) &
             &        - sgmxy(ij   ) * hxt(ij   )) * rym(ij)) * &
             &    0.5d0 * rxu(ij) * ryu(ij) &
             &  + (  sgmxy(ijlne) + sgmxy(ijle) &
             &     + sgmxy(ijln) + sgmxy(ij)) * 0.25d0 * hxyu(ij) &
             &  - (  sgmyy(ijlne) + sgmyy(ijle) &
             &     + sgmyy(ijln) + sgmyy(ij)) * 0.25d0 * hyxu(ij) &
             &  + accelu(ij) &
             &  + ctau(ij) * &
             &    (  ux(ij, kglev)*amskv(ij, kglev) * cioc &
             &     - vx(ij, kglev)*amskv(ij, kglev) &
             &         * sign(cios, cor(ij))) &
             &  + avrmsx(ij) * rdts * uix(ij)
           delt = (  (  sgmxy(ijlne) * hyt(ijlne) &
             &        + sgmxy(ijle ) * hyt(ijle ) &
             &        - sgmxy(ijln ) * hyt(ijln ) &
             &        - sgmxy(ij   ) * hyt(ij   )) * rx &
             &     + (  sgmyy(ijlne) * hxt(ijlne) &
             &        + sgmyy(ijln ) * hxt(ijln ) &
             &        - sgmyy(ijle ) * hxt(ijle ) &
             &        - sgmyy(ij   ) * hxt(ij   )) * rym(ij)) * &
             &    0.5d0 * rxu(ij) * ryu(ij) &
             &  + (  sgmxy(ijlne) + sgmxy(ijle) &
             &     + sgmxy(ijln) + sgmxy(ij)) * 0.25d0 * hyxu(ij) &
             &  - (  sgmxx(ijlne) + sgmxx(ijle) &
             &     + sgmxx(ijln) + sgmxx(ij)) * 0.25d0 * hxyu(ij) &
             &  + accelv(ij) &
             &  + ctau(ij) * &
             &    (  vx(ij, kglev)*amskv(ij, kglev) * cioc &
             &     + ux(ij, kglev)*amskv(ij, kglev) &
             &         * sign(cios, cor(ij))) &
             &  + avrmsx(ij) * rdts * vix(ij)
           raabb = 1.d0 / (alph * alph + beta * beta)
           uix(ij) = (alph * gamm + beta * delt) * raabb
           vix(ij) = (alph * delt - beta * gamm) * raabb
        end if
     end do

#ifdef OPT_PARALLEL
#ifdef OPT_TRIPOLE
     call shift2(   uix,    vix, &
       &          nxdim,  nydim,      1, &
       &          -1.d0,     -1,     -1 )
#else
     call shift2( &
       &            uix,    vix, &
       &          nxdim,  nydim,      1)
#endif
#else
     call stbcvi( &
       &            uix,    vix)
#endif

     call strain( &
       &            exx,    eyy,    exy, &
       &            uix,    vix)
     call rheolo( &
       &            zeta,    eta,    emz,    epz,   pice, &
       &            aice,   mice,    exx,    eyy,    exy )

  end do

  do ij = ijvstr, ijvend
     taux(ij) = (  ctau(ij) * &
       &           (  (uix(ij) &
       &                - ux(ij, kglev)*amskv(ij, kglev)) * cioc &
       &            - (vix(ij) &
       &                - vx(ij, kglev)*amskv(ij, kglev)) * &
       &              sign(cios, cor(ij))) &
       &         + tauaox(ij) * caic &
       &         - tauaoy(ij) * sign(cais, cor(ij))) * &
       &        amskv(ij, kstr)
     tauy(ij) = (  ctau(ij) * &
       &           (  (vix(ij) &
       &                - vx(ij, kglev)*amskv(ij, kglev)) * cioc &
       &            + (uix(ij) &
       &                - ux(ij, kglev)*amskv(ij, kglev)) * &
       &              sign(cios, cor(ij))) &
       &         + tauaoy(ij) * caic &
       &         + tauaox(ij) * sign(cais, cor(ij))) * &
       &        amskv(ij, kstr)
  end do

  call puttao( &
    &            tauaox, tauaoy, &
    &              caic,   cais )

  return
end subroutine pmomnt

! *********************************************************************

subroutine strain( &
  &                   exx,    eyy,    exy, &
  &                    ui,     vi)

  real(8), intent(out) ::    exx(nxydim),    eyy(nxydim),    exy(nxydim)
  real(8), intent(in)  ::     ui(nxydim),     vi(nxydim)

  integer ::     ij,   ijlw,   ijls,  ijlsw

  do ij = ijtstr, ijtend+nxdim+1
     ijlw = ij + lw
     ijls = ij + ls
     ijlsw = ij + lsw
     exx(ij) = (ui(ij) + ui(ijls) - ui(ijlw) - ui(ijlsw)) * &
       &       rx * rxt(ij) * 0.5d0 &
       &     + (vi(ij) + vi(ijls) + vi(ijlw) + vi(ijlsw)) * &
       &       hxyt(ij) * 0.25d0
     eyy(ij) = (vi(ij) + vi(ijlw) - vi(ijls) - vi(ijlsw)) * &
       &       ry(ij) * ryt(ij) * 0.5d0 &
       &     + (ui(ij) + ui(ijlw) + ui(ijls) + ui(ijlsw)) * &
       &       hyxt(ij) * 0.25d0
     exy(ij) = (ui(ij) + ui(ijlw) - ui(ijls) - ui(ijlsw)) * &
       &       ry(ij) * ryt(ij) * 0.25d0 &
       &     - (vi(ij) + vi(ijlw) + vi(ijls) + vi(ijlsw)) * &
       &       hyxt(ij) * 0.125d0 &
       &     + (vi(ij) + vi(ijls) - vi(ijlw) - vi(ijlsw)) * &
       &       rx * rxt(ij) * 0.25d0 &
       &     - (ui(ij) + ui(ijls) + ui(ijlw) + ui(ijlsw)) * &
       &       hxyt(ij) * 0.125d0
  end do

  return
end subroutine strain

! *********************************************************************

subroutine rheolo( &
  &            zeta,    eta,    emz,    epz,   pice, &
  &            aice,   mice,    exx,    eyy,    exy )

  real(8), intent(out) ::   zeta(nxydim),    eta(nxydim)
  real(8), intent(out) ::    emz(nxydim),    epz(nxydim)
  real(8), intent(out) ::   pice(nxydim)
  real(8), intent(in)  ::   aice(nxydim),   mice(nxydim)
  real(8), intent(in)  ::    exx(nxydim),    eyy(nxydim),    exy(nxydim)

  real(8) ::  delta(nxydim)

  real(8), save ::     c1,     c2,     c3,     c4
  logical, save :: ofirst = .true.

  real(8) ::    del
  integer ::     ij
  integer ::  ifpar,  jfpar,  istat

  real(8), save ::  ecc = 2.0d0,  dmin = 2.0d-7,  floss = 17.0d0
  real(8), save ::  p0 = 2.0d5,  cp = 2.0d1

  namelist /nmidyn/ ecc, dmin, floss
  namelist /nmpice/ p0, cp

  if (ofirst) then
     ofirst = .false.
     call rewnml(ifpar, jfpar)
     read (ifpar, nmidyn, iostat=istat)
     call cstnml(jfpar, 'rheolo', 'nmidyn', istat)
     write(jfpar, nmidyn)
     call rewnml(ifpar, jfpar)
     read (ifpar, nmpice, iostat=istat)
     call cstnml(jfpar, 'rheolo', 'nmpice', istat)
     write(jfpar, nmpice)

     c1 = 1.d0 + 1.d0 / ecc / ecc
     c2 = 4.d0 / ecc / ecc
     c3 = 2.d0 * (1.d0 - 1.d0 / ecc / ecc)
     c4 = 1.d0 / ecc / ecc
  end if

  do ij = ijtstr, ijtend+nxdim+1
     pice(ij) = p0 * mice(ij) * exp(- cp * (1.d0 - aice(ij)))
     del = sqrt(  c1 * (exx(ij) * exx(ij) + eyy(ij) * eyy(ij)) &
       &        + c2 * exy(ij) * exy(ij) &
       &        + c3 * exx(ij) * eyy(ij))
     delta(ij) = max(del, dmin)
     zeta(ij) = pice(ij) / delta(ij) * 0.5d0
     eta (ij) = zeta(ij) * c4
     emz (ij) = eta(ij) - zeta(ij)
     epz (ij) = eta(ij) + zeta(ij)
  end do

  return
end subroutine rheolo

end module ipmmt
