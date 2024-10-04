module iprdg

! ---- information ----------------------------------------------------
!
!  Mechanical redistribution of sea ice among thickness categories.
!
!  HISTORY
!     '03.08.06  H.Hasumi
!     '07.09.25  H.Hasumi: for COCO4
!     '07.10.03  H.Hasumi: 1-layer sea ice thermodynamics
!     '08.02.??  Y.Komuro
!     '08.06.11  H.Hasumi: initial/final processing
!     '08.07.10  H.Hasumi: initial/final processing
!     '08.11.22  Y.Komuro: coping with the situation HRIDGE = HIX
!     '12.07.20  Y.Komuro: for COCO5.0
!     '21.05.26  Y.Komuro: snow aging & meltpond parametrization 
!
! ---------------------------------------------------------------------

  use zocdim, only: &
    & nxydim,    nic, ijtstr, ijtend, &
    &     lw,     ls,    lsw, &
    &  oinit, ofinal, myrank
  use zocgrd, only: &
    &     rx,     ry,    hic,     ts, &
    &    hxt,   hxyt,   hyxt,    rxt,    ryt
  use zocphy, only: &
    & gravit,   rhoo,   rhoi,   hfus,    cpi,   dtds

  implicit none

  private

  public :: pridge

contains

subroutine pridge( &
  &                  pice, &
  &                    ax,    hix,    eix,    hsx,    tix, &
  &                   asx,  frlvx,   vmpx,  frmpx,   dsdx,   dsbx, &
  &                    az,    hiz,    eiz,    hsz, &
  &                    ui,     vi )
  use ufile
  use zocite

  integer, parameter :: nrbnd = 3  !! N. of radiation bands; VI, NIR, and IR

  real(8), intent(out)   ::    pice(nxydim)
  real(8), intent(inout) ::      ax(nxydim, 0:nic)
  real(8), intent(inout) ::     hix(nxydim, 0:nic)
  real(8), intent(inout) ::     eix(nxydim, 0:nic)
  real(8), intent(inout) ::     hsx(nxydim, 0:nic)
  real(8), intent(inout) ::     tix(nxydim, 0:nic)
  real(8), intent(inout) ::     asx(nxydim, 0:nic)
  real(8), intent(inout) ::   frlvx(nxydim, 0:nic)
  real(8), intent(inout) ::    vmpx(nxydim, 0:nic)
  real(8), intent(inout) ::   frmpx(nxydim, 0:nic)
  real(8), intent(inout) ::    dsdx(nxydim, 0:nic)
  real(8), intent(inout) ::    dsbx(nxydim, 0:nic)
  real(8), intent(in)    ::      az(nxydim, 0:nic)
  real(8), intent(in)    ::     hiz(nxydim, 0:nic)
  real(8), intent(in)    ::     eiz(nxydim, 0:nic)
  real(8), intent(in)    ::     hsz(nxydim, 0:nic)
  real(8), intent(in)    ::      ui(nxydim),     vi(nxydim)

  real(8) ::   axhix(nxydim, 0:nic),  axhsx(nxydim, 0:nic)
  real(8) ::   axeix(nxydim, 0:nic)
  real(8) ::   axasx(nxydim, 0:nic),  axvmp(nxydim, 0:nic)
  real(8) ::   axflv(nxydim, 0:nic),  axfmp(nxydim, 0:nic)
  real(8) ::   axdsd(nxydim, 0:nic),  axdsb(nxydim, 0:nic)
  real(8) ::     axa(nxydim, 0:nic)
  real(8) ::    divv(nxydim),  edis(nxydim)
  real(8) ::      wa(nxydim, 0:nic),     wn(nxydim, 0:nic)
  real(8) ::      ww(nxydim)
  real(8) ::       g(nxydim, -1:nic)
  real(8) ::       y(nxydim, -1:nic)
  real(8) ::      da(nxydim, nic)
  real(8) ::    dahi(nxydim, nic),   dahs(nxydim, nic)
  real(8) ::    daei(nxydim, nic)
  real(8) ::    daas(nxydim, nic),   davm(nxydim, nic)
  real(8) ::    dafl(nxydim, nic),   dafm(nxydim, nic)
  real(8) ::    dadd(nxydim, nic),   dadb(nxydim, nic)
  real(8) ::  hrdgef(nxydim)
!  common /work/ axhix, axhsx, &
!    &           divv, edis, wa, wn, ww, &
!    &           g, y, da, dahi, dahs, daei

  real(8), save ::      c1,     c2,     c3
  real(8), save ::     gam(nxydim, nic, nic)
  real(8), save ::   pifct
  real(8), save ::     tmi
  logical, save ::  ofirst = .true.

  real(8) ::     exx,    eyy,    exy,    del
  real(8) ::   hrmax,  hrmin
  real(8) ::    hikl,   hskl,   hekl,   vmkl,   ddkl,   dbkl
  real(8) ::    cwan
  integer ::      ij,      k,      l
  integer ::    ijlw,   ijls,  ijlsw
  integer ::   ifpar,  jfpar,  istat

  real(8), save ::  ecc = 2.0d0,  dmin = 2.0d-7,  floss = 17.0d0
  logical, save ::  opt_pice = .false.
  real(8), save ::  cs = 0.5d0,  gridge = 0.15d0,  hridge = 1.0d4
  real(8), save ::  si = 5.0d0
! namelist nmmpnd
  integer, save :: impnd = 0 !! 0: melt pond (MP) parametrization not used
                             !! 1: Holland et al. (2012) MP param.
                             !! 2: Hunke et al. (2013) MP param.
  real(8), save :: hminmp = 10.0d0 !! min. ice thickness for keeping MP [cm]
  real(8), save :: rtdpmp = 80.0d0 !! ratio of melt pond depth to frmp [cm/1]
  real(8), save :: rtmxmp = 0.9d0  !! max. ratio of MP depth to ice thickness
  !! The default dpscl in CICE is 1.0, but we set it to 0.1.
  !! (maybe due to slightly different implementation?)
  real(8), save :: dpscl = 0.1d0  !! permiability scale parameter [ND]
  real(8), save :: rmpcmn(0:2) = & !! minimum water catching rate of MP 
    &                 (/ 0.0d0, 0.15d0, 0.15d0 /) 
  real(8), save :: rmpcmx(0:2) = & !! maximum water catching rate of MP 
    &                 (/ 0.0d0, 0.7d0, 0.85d0 /) 
  real(8), save :: cmpfrz = 3.d-6 !! constant for melt pond freeze-up rate
  real(8), save :: tmpfrz = -2.0d0 !! ref. t for melt pond freeze-up [c]
  real(8), save :: albmpd( nrbnd ) = &  !! deep melt pond shortwave albedo
    &                 (/ 0.4d0, 0.1d0, 0.0d0 /)
  real(8), save :: almpdp(2) = &   !! MP sw albedo, depth dependency [cm]
    &                 (/ 0.5d0, 20.0d0 /)
  real(8), save :: frmpmn = 1.0d-14  !! empirical limiter for frmpx [ND]
  real(8), save :: vmpmin = 1.0d-12  !! empirical limiter for vmpx [cm]

  namelist /nmidyn/ ecc, dmin, floss, opt_pice
  namelist /nmirdg/ cs, gridge, hridge
  namelist /nmislt/ si
  namelist /nmmpnd/   impnd, hminmp, rtdpmp, rtmxmp,  dpscl, &
    &                rmpcmn, rmpcmx, cmpfrz, tmpfrz, albmpd, almpdp, &
    &                frmpmn, vmpmin

  if (oinit .or. ofinal) then
     return
  end if
  
  if (ofirst) then
     ofirst = .false.
     call rewnml(ifpar, jfpar)
     read (ifpar, nmidyn, iostat=istat)
     call cstnml(jfpar, 'pridge', 'nmidyn', istat)
     write(jfpar, nmidyn)
     call rewnml(ifpar, jfpar)
     read (ifpar, nmirdg, iostat=istat)
     call cstnml(jfpar, 'pridge', 'nmirdg', istat)
     write(jfpar, nmirdg)
     call rewnml(ifpar, jfpar)
     read (ifpar, nmislt, iostat=istat)
     call cstnml(jfpar, 'pridge', 'nmislt', istat)
     write(jfpar, nmislt)
     call rewnml(ifpar, jfpar)
     read (ifpar, nmmpnd, iostat=istat)
     call cstnml(jfpar, 'pridge', 'nmmpnd', istat)
     write(jfpar, nmmpnd)

     tmi = dtds * si

     c1 = 1.d0 + 1.d0 / ecc / ecc
     c2 = 4.d0 / ecc / ecc
     c3 = 2.d0 * (1.d0 - 1.d0 / ecc / ecc)

     pifct = 0.5d0 * rhoi / rhoo * gravit * (rhoo - rhoi) * floss

     do l = 1, nic
        do k = 1, nic
           do ij = 1, nxydim
              gam(ij, k, l) = 0.d0
           end do
        end do
     end do
  end if

  do k = 0, nic
     do ij = 1, nxydim
        axhix(ij, k) = ax(ij, k) * hix(ij, k)
        axhsx(ij, k) = ax(ij, k) * hsx(ij, k)
        axeix(ij, k) = ax(ij, k) * eix(ij, k)
        axasx(ij, k) = ax(ij, k) * asx(ij, k)
        axflv(ij, k) = ax(ij, k) * frlvx(ij, k)
        axvmp(ij, k) = ax(ij, k) * vmpx(ij, k)
        axfmp(ij, k) = ax(ij, k) * frmpx(ij, k)
        axdsd(ij, k) = ax(ij, k) * dsdx(ij, k)
        axdsb(ij, k) = ax(ij, k) * dsbx(ij, k)
        g(ij, k) = 0.d0
        wn(ij, k) = 0.d0
     end do
  end do
  do ij = 1, nxydim
     pice(ij) = 0.d0
     g(ij, -1) = 0.d0
     hrdgef(ij) = hridge
  end do

  do ij = ijtstr, ijtend
     ijlw = ij + lw
     ijls = ij + ls
     ijlsw = ij + lsw
     exx = (ui(ij) + ui(ijls) - ui(ijlw) - ui(ijlsw)) * &
       &    rx * rxt(ij) * 0.5d0 &
       & + (vi(ij) + vi(ijls) + vi(ijlw) + vi(ijlsw)) * &
       &    hxyt(ij) * 0.25d0
     eyy = (vi(ij) + vi(ijlw) - vi(ijls) - vi(ijlsw)) * &
       &    ry(ij) * ryt(ij) * 0.5d0 &
       & + (ui(ij) + ui(ijlw) + ui(ijls) + ui(ijlsw)) * &
       &    hyxt(ij) * 0.25d0
     exy = (ui(ij) + ui(ijlw) - ui(ijls) - ui(ijlsw)) * &
       &    ry(ij) * ryt(ij) * 0.25d0 &
       & - (vi(ij) + vi(ijlw) + vi(ijls) + vi(ijlsw)) * &
       &    hyxt(ij) * 0.125d0 &
       & + (vi(ij) + vi(ijls) - vi(ijlw) - vi(ijlsw)) * &
       &    rx * rxt(ij) * 0.25d0 &
       & - (ui(ij) + ui(ijls) + ui(ijlw) + ui(ijlsw)) * &
       &    hxyt(ij) * 0.125d0
     divv(ij) = exx + eyy
     del = sqrt(  c1 * (exx * exx + eyy * eyy) &
       &        + c2 * exy * exy + c3 * exx * eyy)
     edis(ij) = 0.5d0 * cs * (del - abs(divv(ij))) &
       &      - min(divv(ij), 0.d0)
  end do

  do ij = ijtstr, ijtend
     g(ij, 0) = max(ax(ij, 0), 0.0d0)
  end do
  do k = 1, nic
     do ij = ijtstr, ijtend
        g(ij, k) = g(ij, k-1) + ax(ij, k)
     end do
  end do

  do k = -1, nic
     do ij = ijtstr, ijtend
        y(ij, k) = max(0.d0, 1.d0 - g(ij, k) / gridge)**2
     end do
  end do

  do k = 0, nic
     do ij = ijtstr, ijtend
        wa(ij, k) = y(ij, k-1) - y(ij, k)
        if (wa(ij, k) > 0.0d0) then
           hrdgef(ij) = max(hix(ij, k), hridge)
        end if
     end do
  end do

  do l = 1, nic
! === '23.07.12: avoid loop interchange due to a bug
!                 in ES4ve nfort compiler (version 5.0.0 or earlier)
!NEC$ nointerchange
   do k = 1, l
        do ij = ijtstr, ijtend
           hrmax = 2.d0 * sqrt(hrdgef(ij) * hix(ij, k))
           hrmin = 2.d0 * hix(ij, k)
           if (     (hic(l+1) .lt. hrmin) &
             & .or. (hic(l) .gt. hrmax) ) then
              gam(ij, k, l) = 0.d0
           else if (hrmin .ge. hrmax) then
              gam(ij, k, l) = 0.5d0
           else
              gam(ij, k, l) = (  min(hrmax, hic(l+1)) &
                &              - max(hrmin, hic(l))) * 0.5d0 &
                &             / (hrdgef(ij) - hix(ij, k))
           end if
        end do
     end do
  end do

  do l = 1, nic
     do k = 1, l
        do ij = ijtstr, ijtend
           wn(ij, l) = wn(ij, l) + wa(ij, k) * gam(ij, k, l)
        end do
     end do
  end do

  do ij = 1, nxydim
     ww(ij) = 0.d0
  end do
  do k = 0, nic
     do ij = ijtstr, ijtend
        ww(ij) = ww(ij) + wa(ij, k) - wn(ij, k)
     end do
  end do
  do k = 0, nic
     do ij = ijtstr, ijtend
        wa(ij, k) = wa(ij, k) / ww(ij)
        wn(ij, k) = wn(ij, k) / ww(ij)
     end do
  end do

  do k = 1, nic
     do ij = ijtstr, ijtend
         if (edis(ij)*wa(ij, k)*ts > ax(ij, k)) then
            cwan = ax(ij, k) / (edis(ij)*wa(ij, k)*ts)
            do l = k, nic
               wn(ij, l) = wn(ij, l) &
                 &       - (1.0d0 - cwan) * wa(ij, k) * gam(ij, k, l)
            end do
            wa(ij, k) = cwan * wa(ij, k)
         end if
      end do
   end do
 
   do k = 1, nic
     do ij = ijtstr, ijtend
        da(ij, k) = edis(ij) * (wn(ij, k) - wa(ij, k))
     end do
  end do
  do l = 1, nic
     do ij = ijtstr, ijtend
        dahi(ij, l) = - hix(ij, l) * wa(ij, l) * edis(ij)
        dahs(ij, l) = - hsx(ij, l) * wa(ij, l) * edis(ij)
        daei(ij, l) = - eix(ij, l) * wa(ij, l) * edis(ij)
        daas(ij, l) = - asx(ij, l) * wa(ij, l) * edis(ij)
        dafl(ij, l) = - frlvx(ij, l) * wa(ij, l) * edis(ij)
        davm(ij, l) = - vmpx(ij, l) * wa(ij, l) * edis(ij)
        dafm(ij, l) = - frmpx(ij, l) * wa(ij, l) * edis(ij)
        dadd(ij, l) = - dsdx(ij, l) * wa(ij, l) * edis(ij)
        dadb(ij, l) = - dsbx(ij, l) * wa(ij, l) * edis(ij)
        pice(ij) = pice(ij) &
          &      - hix(ij, l) * hix(ij, l) * wa(ij, l) * pifct
        do k = 1, l
           hrmax = 2.d0 * sqrt(hrdgef(ij) * hix(ij, k))
           hrmin = 2.d0 * hix(ij, k)
           hikl = (  min(hrmax, hic(l+1)) &
             &     + max(hrmin, hic(l))) * 0.5d0
           hskl = hsx(ij, k) * (hix(ij, k) + hrmax * 0.5d0) &
             &    / hix(ij, k)
           if (tix(ij, k) .lt. tmi) then
              hekl = hikl * ei(tix(ij, k), si)
           else
              hekl = 0.d0
           end if
           vmkl = vmpx(ij, k) * (hix(ij, k) + hrmax * 0.5d0) &
             &    / hix(ij, k)
           ddkl = dsdx(ij, k) * (hix(ij, k) + hrmax * 0.5d0) &
             &    / hix(ij, k)
           dbkl = dsbx(ij, k) * (hix(ij, k) + hrmax * 0.5d0) &
             &    / hix(ij, k)
           dahi(ij, l) = dahi(ij, l) &
             &         + hikl * wa(ij, k) * gam(ij, k, l) * &
             &           edis(ij)
           dahs(ij, l) = dahs(ij, l) &
             &         + hskl * wa(ij, k) * gam(ij, k, l) * &
             &           edis(ij)
           daei(ij, l) = daei(ij, l) &
             &         + hekl * wa(ij ,k) * gam(ij, k, l) * &
             &           edis(ij)
           daas(ij, l) = daas(ij, l) &
             &         + asx(ij, k) * wa(ij, k) * gam(ij, k, l) * &
             &           edis(ij)
!          frlvx not transported to the destination category;
!           all the ridged ice is classed as deformed ice.
           if (impnd == 1) then
              davm(ij, l) = davm(ij, l) &
                &         + vmkl * wa(ij, k) * gam(ij, k, l) * &
                &           edis(ij)
           end if
!          For Holland et al. (2012) MP parametrization (impnd = 1),
!           the update of frmpx here is just a dummy.
!          For Hunke et al. (2013) MP paramerization (impnd = 2), 
!           vmpx and frmpx not transported to the destination category;
!           the pond water on the ridged ice is regarded to drop to the ocean.
           dadd(ij, l) = dadd(ij, l) &
             &         + ddkl * wa(ij, k) * gam(ij, k, l) * &
             &           edis(ij)
           dadb(ij, l) = dadb(ij, l) &
             &         + dbkl * wa(ij, k) * gam(ij, k, l) * &
             &           edis(ij)
           pice(ij) = pice(ij) &
             &      + hikl * hikl * wa(ij, k) * gam(ij, k, l) * &
             &        pifct
        end do
     end do
  end do

  do k = 1, nic
     do ij = ijtstr, ijtend
        ax(ij, k) = ax(ij, k) + ts * da(ij, k)
!       Store ax before the adjustment
        axa(ij, k) = ax(ij, k)
        ax(ij, k) = min(1.d0, max(0.d0, ax(ij, k)))
     end do
  end do
  do ij = ijtstr, ijtend
     ax(ij, 0) = 1.d0
  end do
  do k = 1, nic
     do ij = ijtstr, ijtend
        ax(ij, 0) = ax(ij, 0) - ax(ij, k)
     end do
  end do
  do ij = ijtstr, ijtend
     if (ax(ij, 0) .lt. 0.d0) then
        do k = 1, nic
           ax(ij, k) = ax(ij, k) / (1.d0 - ax(ij, 0))
        end do
        ax(ij, 0) = 0.d0
     end if
  end do

  do ij = ijtstr, ijtend
     hix(ij, 0) = 0.d0
     hsx(ij, 0) = 0.d0
     eix(ij, 0) = 0.d0
     vmpx(ij, 0) = 0.d0
     dsdx(ij, 0) = 0.d0
     dsbx(ij, 0) = 0.d0
     if (ax(ij, 0) .eq. 1.d0) then
        do k = 1, nic
           ax(ij, k) = 0.d0
           hix(ij, k) = hic(k)
           hsx(ij, k) = 0.d0
           eix(ij, k) = 0.d0
           asx(ij, k) = 0.d0
           frlvx(ij, k) = 1.0d0
           vmpx(ij, k) = 0.d0
           frmpx(ij, k) = 0.0d0
           dsdx(ij, k) = 0.d0
           dsbx(ij, k) = 0.d0
           hix(ij, 0) = hix(ij, 0) + axhix(ij, k)
           hsx(ij, 0) = hsx(ij, 0) + axhsx(ij, k)
           eix(ij, 0) = eix(ij, 0) + axeix(ij, k)
           vmpx(ij, 0) = vmpx(ij, 0) + axvmp(ij, k)
           dsdx(ij, 0) = dsdx(ij, 0) + axdsd(ij, k)
           dsbx(ij, 0) = dsbx(ij, 0) + axdsb(ij, k)
           dahi(ij, k) = 0.d0
           dahs(ij, k) = 0.d0
           daei(ij, k) = 0.d0
           daas(ij, k) = 0.d0
           dafl(ij, k) = 0.d0
           davm(ij, k) = 0.d0
           dafm(ij, k) = 0.d0
           dadd(ij, k) = 0.d0
           dadb(ij, k) = 0.d0
        end do
     end if
  end do

  do k = 1, nic
     do ij = ijtstr, ijtend
        axhix(ij, k) = axhix(ij, k) + ts * dahi(ij, k)
        axhsx(ij, k) = axhsx(ij, k) + ts * dahs(ij, k)
        axeix(ij, k) = axeix(ij, k) + ts * daei(ij, k)
        axasx(ij, k) = axasx(ij, k) + ts * daas(ij, k)
        axflv(ij, k) = axflv(ij, k) + ts * dafl(ij, k)
        axvmp(ij, k) = axvmp(ij, k) + ts * davm(ij, k)
        axfmp(ij, k) = axfmp(ij, k) + ts * dafm(ij, k)
        axdsd(ij, k) = axdsd(ij, k) + ts * dadd(ij, k)
        axdsb(ij, k) = axdsb(ij, k) + ts * dadb(ij, k)
     end do
  end do

! *** ridging process conserves the mass of ice and snow, and the
! *** amount of the thickest category ice is never reduced by ridging.
! *** thus, the situation that axhix(ij, nic) becomes less than zero
! *** needs not be considered in the adjustment below. 
!
! '21.07.19: For variables which are not conserved through the process
!            (asx, frlvx, frmpx, and vmpx for impnd=2), negative
!            ax-values are not transported to the thicker category; 
!            the negative ax-value transport sometimes causes resultant
!            negative value, which is invalid, in the thickest category,
!            since positive counterpart has not been transported there.
  do k = 1, nic-1
     do ij = ijtstr, ijtend
        if (ax(ij, k) .eq. 0.d0) then
           axhix(ij, k+1) = axhix(ij, k+1) + axhix(ij, k)
           axhsx(ij, k+1) = axhsx(ij, k+1) + axhsx(ij, k)
           axeix(ij, k+1) = axeix(ij, k+1) + axeix(ij, k)
           axasx(ij, k+1) = axasx(ij, k+1) + max(axasx(ij, k), 0.d0)
           if (impnd == 1) then
              axvmp(ij, k+1) = axvmp(ij, k+1) + axvmp(ij, k)
           else
              axvmp(ij, k+1) = axvmp(ij, k+1) + max(axvmp(ij, k), 0.d0)
           end if
           axflv(ij, k+1) = axflv(ij, k+1) + max(axflv(ij, k), 0.d0)
           axfmp(ij, k+1) = axfmp(ij, k+1) + max(axfmp(ij, k), 0.d0)
!           axflv(ij, k+1) = axflv(ij, k+1) + axflv(ij, k)
!           axfmp(ij, k+1) = axfmp(ij, k+1) + axfmp(ij, k)
           axdsd(ij, k+1) = axdsd(ij, k+1) + axdsd(ij, k)
           axdsb(ij, k+1) = axdsb(ij, k+1) + axdsb(ij, k)
           axhix(ij, k) = 0.d0
           axhsx(ij, k) = 0.d0
           axeix(ij, k) = 0.d0
           axasx(ij, k) = 0.d0
           axflv(ij, k) = 0.d0
           axvmp(ij, k) = 0.d0
           axfmp(ij, k) = 0.d0
           axdsd(ij, k) = 0.d0
           axdsb(ij, k) = 0.d0
           hix(ij, k) = hic(k)
           hsx(ij, k) = 0.d0
           eix(ij, k) = 0.d0
           asx(ij, k) = 0.d0
           frlvx(ij, k) = 1.d0
           vmpx(ij, k) = 0.d0
           frmpx(ij, k) = 0.d0
           dsdx(ij, k) = 0.d0
           dsbx(ij, k) = 0.d0
        else if (axhix(ij, k) .le. 0.d0) then
           axhix(ij, k+1) = axhix(ij, k+1) + axhix(ij, k)
           axhsx(ij, k+1) = axhsx(ij, k+1) + axhsx(ij, k)
           axeix(ij, k+1) = axeix(ij, k+1) + axeix(ij, k)
           axasx(ij, k+1) = axasx(ij, k+1) + max(axasx(ij, k), 0.d0)
           if (impnd == 1) then
              axvmp(ij, k+1) = axvmp(ij, k+1) + axvmp(ij, k)
           else
              axvmp(ij, k+1) = axvmp(ij, k+1) + max(axvmp(ij, k), 0.d0)
           end if
           axflv(ij, k+1) = axflv(ij, k+1) + max(axflv(ij, k), 0.d0)
           axfmp(ij, k+1) = axfmp(ij, k+1) + max(axfmp(ij, k), 0.d0)
!           axflv(ij, k+1) = axflv(ij, k+1) + axflv(ij, k)
!           axfmp(ij, k+1) = axfmp(ij, k+1) + axfmp(ij, k)
           axdsd(ij, k+1) = axdsd(ij, k+1) + axdsd(ij, k)
           axdsb(ij, k+1) = axdsb(ij, k+1) + axdsb(ij, k)
           ax(ij, k+1) = ax(ij, k+1) + ax(ij, k)
           axa(ij, k+1) = axa(ij, k+1) + ax(ij, k)
           ax(ij, k) = 0.d0
           axhix(ij, k) = 0.d0
           axhsx(ij, k) = 0.d0
           axeix(ij, k) = 0.d0
           axasx(ij, k) = 0.d0
           axflv(ij, k) = 0.d0
           axvmp(ij, k) = 0.d0
           axfmp(ij, k) = 0.d0
           axdsd(ij, k) = 0.d0
           axdsb(ij, k) = 0.d0
           hix(ij, k) = hic(k)
           hsx(ij, k) = 0.d0
           eix(ij, k) = 0.d0
           asx(ij, k) = 0.d0
           frlvx(ij, k) = 1.d0
           vmpx(ij, k) = 0.d0
           frmpx(ij, k) = 0.d0
           dsdx(ij, k) = 0.d0
           dsbx(ij, k) = 0.d0
        else
           if (axhsx(ij, k) .lt. 0.d0) then
              axhsx(ij, k+1) = axhsx(ij, k+1) + axhsx(ij, k)
              axhsx(ij, k) = 0.d0
              hsx(ij, k) = 0.d0
              axasx(ij, k+1) = axasx(ij, k+1) + max(axasx(ij, k), 0.d0)
              axasx(ij, k) = 0.d0
              asx(ij, k) = 0.d0
           end if
           if (axeix(ij, k) .lt. 0.d0) then
              axeix(ij, k+1) = axeix(ij, k+1) + axeix(ij, k)
              axeix(ij, k) = 0.d0
              eix(ij, k) = 0.d0
           end if
           if (axasx(ij, k) .lt. 0.d0) then
!              axasx(ij, k+1) = axasx(ij, k+1) + axasx(ij, k)
              axasx(ij, k) = 0.d0
              asx(ij, k) = 0.d0
           end if
           if (axflv(ij, k) .lt. 0.d0) then
!              axflv(ij, k+1) = axflv(ij, k+1) + axflv(ij, k)
              axflv(ij, k) = 0.d0
              frlvx(ij, k) = 0.d0
           end if
           if (axvmp(ij, k) .lt. 0.d0) then
              if (impnd == 1) then
                 axvmp(ij, k+1) = axvmp(ij, k+1) + axvmp(ij, k)
              end if
              axvmp(ij, k) = 0.d0
              vmpx(ij, k) = 0.d0
              axfmp(ij, k+1) = axfmp(ij, k+1) + max(axfmp(ij, k), 0.d0)
              axfmp(ij, k) = 0.d0
              frmpx(ij, k) = 0.d0
           end if
           if (axfmp(ij, k) .lt. 0.d0) then
!              axfmp(ij, k+1) = axfmp(ij, k+1) + axfmp(ij, k)
              axfmp(ij, k) = 0.d0
              frmpx(ij, k) = 0.d0
!             frmpx update is dummy for Holland MP param. (impnd = 1)
              if (impnd == 2) then  
                 axvmp(ij, k+1) = axvmp(ij, k+1) + max(axvmp(ij, k), 0.d0)
                 axvmp(ij, k) = 0.d0
                 vmpx(ij, k) = 0.d0
              end if
           end if
           if (axdsd(ij, k) .lt. 0.d0) then
              axdsd(ij, k+1) = axdsd(ij, k+1) + axdsd(ij, k)
              axdsd(ij, k) = 0.d0
              dsdx(ij, k) = 0.d0
           end if
           if (axdsb(ij, k) .lt. 0.d0) then
              axdsb(ij, k+1) = axdsb(ij, k+1) + axdsb(ij, k)
              axdsb(ij, k) = 0.d0
              dsbx(ij, k) = 0.d0
           end if
        end if
     end do
  end do

  do k = 1, nic
     do ij = ijtstr, ijtend
        if (ax(ij, k) .gt. 0.d0) then
           hix(ij, k) = axhix(ij, k) / ax(ij, k)
           hsx(ij, k) = axhsx(ij, k) / ax(ij, k)
           eix(ij, k) = axeix(ij, k) / ax(ij, k)
           tix(ij, k) = ti(eix(ij, k)/hix(ij, k), si)
           vmpx(ij, k) = axvmp(ij, k) / ax(ij, k)
           dsdx(ij, k) = axdsd(ij, k) / ax(ij, k)
           dsbx(ij, k) = axdsb(ij, k) / ax(ij, k)
!          Use unadjusted ax to calculate the following variables, 
!           since they are not affected by the adjustment:
!           snow age, level-ice fractaion, and pond fraction.
           if (axa(ij, k) .gt. 0.d0) then
              asx(ij, k) = axasx(ij, k) / axa(ij, k)
              frlvx(ij, k) = axflv(ij, k) / axa(ij, k)
              frmpx(ij, k) = axfmp(ij, k) / axa(ij, k)
           else
              asx(ij, k) = 0.d0
              frlvx(ij, k) = 1.0d0
              frmpx(ij, k) = 0.0d0
!              write(0, *) '### REFRESH ASX/FRLVX/FRMPX (iprdg) ###' !! debug
           end if
        else
           tix(ij, k) = tmi
        end if
!       debug code
!        if ((asx(ij, k) < 0.d0).or.(asx(ij, k) > (1.0d0+1.0d-9))) then
!           write(0,*) '##iprdg; asx##', myrank, ij, k, asx(ij, k)
!        end if        
!        if ((frlvx(ij, k) < 0.d0).or.(frlvx(ij, k) > (1.0d0+1.0d-9))) then
!            write(0,*) '##iprdg; frlvx##', myrank, ij, k, frlvx(ij, k)
!        end if
!        if ((frmpx(ij, k) < 0.d0).or.(frmpx(ij, k) > (1.0d0+1.0d-9))) then
!         write(0,*) '##iprdg; frmpx##', myrank, ij, k, frmpx(ij, k)
!        end if        
     end do
  end do

  return

end subroutine pridge

end module iprdg
