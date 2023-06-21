module ipthm

! ---- information ----------------------------------------------------
!
!  HISTORY
!     '02.10.10  H.Hasumi: from MIROC3.1
!     '03.07.30  H.Hasumi: multi-category sea ice thickness
!     '07.09.25  H.Hasumi: for COCO4
!     '07.10.03  H.Hasumi: 1-layer sea ice thermodynamics
!     '08.04.07  H.Hasumi: condition is modified for new ice formation
!                          on open water
!     '08.06.11  H.Hasumi: initial/final processing
!     '08.07.10  H.Hasumi: initial/final processing
!     '08.12.12  Y.Komuro: linear remapping of Lipscomb(2001)
!     '09.05.25  Y.Komuro: CMIP5 output code included
!                          (basal/lateral melting processes separated)
!     '12.08.01  Y.Komuro: for COCO5.0
!     '13.02.13  Y.Komuro: remove non-parallel code 
!     '13.10.28  Y.Komuro: bug fix (incorrect IMRISF/IGRCON/IMRIBS)
!     '14.04.30  Y.Komuro: bug fix
!     '21.05.26  Y.Komuro: snow aging & meltpond parametrization 
!
! ---------------------------------------------------------------------

  use zocdim, only: &
    & nxydim,  nzdim,   kstr,  ntdim,     nx,     ny,    nic, &
    &  ijstr,  ijend, ijtstr, ijtend, &
    &     lw,     ls,    lsw, &
    &  oinit, ofinal, myrank, nxyidm
  use zocgrd, only: &
    &    hic,     ts,     dt,    cor
  use zocmsk, only: &
    &  amskt
  use zocphy, only: &
    &   rhoo,   rhoi,   rhos,   hfus,    cpo,    cpi,   dtds, &
    &  tmelt, kelvin, gravit

  implicit none

  integer, parameter :: nrbnd = 3

! [namelist parameters] 
! namelist nmsfrc
  logical, save :: oasfrc = .false. !! true if CCSM type snow fraction
  real(8), save :: alsdpt = 1.0d0
                      !! threshold thickness of snow-covered ice, cm
  real(8), save :: alspat = 2.0d0
                      !! snowpatch depth for CCSM type snow fraction, cm
! namelist nmsage
  logical, save :: osage = .false. !! true if snow aging param. is used
  real(8), save :: alssif( nrbnd ) = &  !! albedo of fresh snow on sea ice
    &                 (/ 0.75d0, 0.75d0, 0.0d0 /) 
  real(8), save :: alssio( nrbnd ) = &  !! albedo of old snow on sea ice
    &                 (/ 0.5d0, 0.5d0, 0.0d0 /) 
  real(8), save :: alfmax = 0.999d0 !! maximum value for albfct
  real(8), save :: snrfrs = 1.0d0 !! snowfall for refreshing snow surface [cm]
  real(8), save :: ftage = 5.0d3  !! aging factor, for agefr1
  real(8), save :: tauage = 2.0d6 !! time scale for aging
  real(8), save :: adirt0 = 0.3d0 !! dirt factor for agefr3, normal place
  real(8), save :: adirtc = 0.01d0 !! dirt factor for agefr3, clean place
  real(8), save :: adirts = 0.1d0 !! dirt factor for agefr3, coeff. for dscppm
  real(8), save :: adirtm = 1.0d0 !! dirt factor for agefr3, maximum
  real(8), save :: drsmax = 0.1d0 !! maximum ratio of dust to snow
  logical, save :: oadst = .false.
                      !! using dsdx/dsbx for aging insted of adirt0/adirtc
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
! namelist nmislt
  real(8), save ::     si = 5.0d0  !! sea-ice salinity (psu)
! namelist nmamin
  real(8), save ::  amin = 1.0d-6  !! concentration min.
  real(8), save ::  amax = 1.0d0   !! concentration max.
  integer, save ::  mic = nic

! namelist nmsaab (from matdrv.F in MIROC5.2/MATSIRO parameters)
  real(8), save :: abduvs = 6.36777D2 !! absorp.coef. for dust vis from PARA.bnd29ch111sp
  real(8), save :: abduni = 3.34617D2 !! absorp.coef. for dust NI from PARA.bnd29ch111sp
  real(8), save :: abduir = 8.62054D2 !! absorp.coef. for dust IR from PARA.bnd29ch111sp
  real(8), save :: abbcvs = 7.43200D4 !! absorp.coef. for BC vis from PARA.bnd29ch111sp
  real(8), save :: abbcni = 2.93200D4 !! absorp.coef. for BC NI from PARA.bnd29ch111sp
  real(8), save :: abbcir = 2.47174D3 !! absorp.coef. for BC IR from PARA.bnd29ch111sp
  real(8), save :: wgtvs  = 0.46D0    !! radiation weight of visible band
  real(8), save :: wgtni  = 0.36D0    !! radiation weight of visible band
  real(8), save :: wgtir  = 0.18D0    !! radiation weight of visible band

  namelist /nmsfrc/  oasfrc, alsdpt, alspat
  namelist /nmsage/   osage, alssif, alssio, alfmax, snrfrs,  ftage, &
    &                tauage, adirt0, adirtc, adirts, adirtm, drsmax, &
    &                 oadst
  namelist /nmmpnd/   impnd, hminmp, rtdpmp, rtmxmp,  dpscl, &
    &                rmpcmn, rmpcmx, cmpfrz, tmpfrz, albmpd, almpdp, &
    &                frmpmn, vmpmin
  namelist /nmsaab/ abduvs, abduni, abduir, &
    &               abbcvs, abbcni, abbcir, &
    &                wgtvs,  wgtni,  wgtir
  namelist /nmislt/      si
  namelist /nmamin/    amin,   amax,    mic

  private

  public :: ptherm, ipsage, idfrmp

contains

subroutine ptherm( &
  &                    ax,    hix,    hsx,    eix,    tix, &
  &                   asx,  frlvx,   vmpx,  frmpx,   dsdx,   dsbx, &
  &                  prec,   snow,     ft,     fs,    fdd,    fdb, &
  &                 ftitd, igrfra, igrcon, igrsni, &
  &                igrsfl, inrlat, &
  &                imrsno, imrsmi, imrisf, imribs, &
  &                impinc, impfrz, improf, &
  &                    tx,    tsi, &
  &                   wio,    wao,    was,    wil, &
  &                  evap,   subi,   roff, wiadjs, weadjs, &
  &                  dfdu,   dfbc, &
  &                   qio )

  use ufile
  use qckag
  use qckot
  use zocite

#include "mpif.h"

  real(8), intent(inout) ::     ax(nxydim, 0:nic),    hix(nxydim, 0:nic)
  real(8), intent(inout) ::    hsx(nxydim, 0:nic)
  real(8), intent(inout) ::    eix(nxydim, 0:nic),    tix(nxydim, 0:nic)
  real(8), intent(inout) ::    asx(nxydim, 0:nic)
  real(8), intent(inout) ::  frlvx(nxydim, 0:nic)
  real(8), intent(inout) ::   vmpx(nxydim, 0:nic)
  real(8), intent(inout) ::  frmpx(nxydim, 0:nic)
  real(8), intent(inout) ::   dsdx(nxydim, 0:nic),   dsbx(nxydim, 0:nic)
  real(8), intent(inout) ::   prec(nxydim),   snow(nxydim)
  real(8), intent(inout) ::     ft(nxydim, ntdim),     fs(nxydim)
  real(8), intent(inout) ::    fdd(nxydim),    fdb(nxydim)
  real(8), intent(out)   ::  ftitd(nxydim)
  real(8), intent(out)   :: igrfra(nxydim), igrcon(nxydim), igrsni(nxydim)
  real(8), intent(out)   :: inrlat(nxydim)
  real(8), intent(out)   :: imrsno(nxydim), imrsmi(nxydim)
  real(8), intent(out)   :: imrisf(nxydim), imribs(nxydim)
  real(8), intent(out)   :: igrsfl(nxydim)
  real(8), intent(out)   :: impinc(nxydim, 0:nic), impfrz(nxydim, 0:nic)
  real(8), intent(out)   :: improf(nxydim, 0:nic)
  real(8), intent(inout) ::    wio(nxydim, nic)
  real(8), intent(in)    ::    was(nxydim, nic)
  real(8), intent(in)    ::    wil(nxydim, nic)
  real(8), intent(in)    ::    wao(nxydim)
  real(8), intent(in)    ::   evap(nxydim),   subi(nxydim, nic)
  real(8), intent(in)    ::   roff(nxydim)
  real(8), intent(in)    :: wiadjs(nxydim), weadjs(nxydim)
  real(8), intent(in)    ::     tx(nxydim, nzdim, ntdim)
  real(8), intent(in)    ::    tsi(nxydim, 0:nic)
  real(8), intent(in)    ::   dfdu(nxydim),   dfbc(nxydim)
  real(8), intent(in)    ::    qio(nxydim, nic)

  real(8) ::     az(nxydim, 0:nic)
  real(8) ::  axhix(nxydim, 0:nic)
  real(8) ::  axhsx(nxydim, 0:nic), axhsxn(nxydim, 0:nic)
  real(8) ::  axeix(nxydim, 0:nic), axeixn(nxydim, 0:nic)
  real(8) ::  axvmp(nxydim, 0:nic)
  real(8) ::  axflv(nxydim, 0:nic), axfmp(nxydim, 0:nic)
  real(8) ::  axdsd(nxydim, 0:nic),  axdsb(nxydim, 0:nic)
  real(8) ::    wai(nxydim, nic)
  real(8) ::     wi(nxydim),     ws(nxydim)
  real(8) ::    wen(nxydim),    wsn(nxydim)
  real(8) ::  rmpcc(nxydim)
  real(8) ::   dvmp(nxydim, 0:nic)
  real(8) ::   dfcb(nxydim)
  real(8) :: dsdrhs(nxydim, 0:nic), dsbrhs(nxydim, 0:nic)
  real(8) ::    hiz(nxydim, 0:nic),   vmpz(nxydim, 0:nic)
  real(8) ::  aflrm(nxydim, 0:nic), aflrmc(nxydim)
  real(8) ::   fdtn(nxydim, 0:nic),  fdtcn(nxydim, 0:nic)
  real(8) ::   hicn(nxydim, 0:nic)
  real(8) ::     g0(nxydim, 0:nic),     g1(nxydim, 0:nic)
  real(8) ::    hil(nxydim, 0:nic),    hir(nxydim, 0:nic)
  real(8) ::     da(nxydim, 0:nic),   dahi(nxydim, 0:nic)
  real(8) ::   dahs(nxydim, 0:nic),   daei(nxydim, 0:nic)
  real(8) ::   daas(nxydim, 0:nic),   davm(nxydim, 0:nic)
  real(8) ::   dafl(nxydim, 0:nic),   dafm(nxydim, 0:nic)
  real(8) ::   dadd(nxydim, 0:nic),   dadb(nxydim, 0:nic)
  real(8) :: laxhix(nxydim, 0:nic), laxhsx(nxydim, 0:nic) 
  real(8) :: laxeix(nxydim, 0:nic)
  real(8) :: daxhit(nxydim, 0:nic), daxhib(nxydim, 0:nic)
  real(8) :: laxasx(nxydim, 0:nic), laxvmp(nxydim, 0:nic)
  real(8) :: laxflv(nxydim, 0:nic), laxfmp(nxydim, 0:nic)
  real(8) :: laxdsd(nxydim, 0:nic), laxdsb(nxydim, 0:nic)
!      COMMON /WORK/ AZ, AXHIX, AXHSX, AXHSXN, WAI,
!     &              WI, WS, WEN, WSN, AXEIX, AXEIXN,
!     &              HIZ, AFLRM, AFLRMC, FDTN, FDTCN, HICN,
!     &              G0, G1, HIL, HIR, DA, DAHI,
!     &              DAHS, DAEI, LAXHIX, LAXHSX, LAXEIX,
!     &              DAXHIT, DAXHIB

  real(8), save ::    rri,    rrs, rorirs
  real(8), save ::  rsfus
  real(8), save ::    tmi
  real(8), save :: epsaei
  real(8), save ::   hic0(0:nic+1)
  logical, save :: ofirst = .true.

  real(8) ::   wres
  real(8) ::   hsxo,    dhs,   vmpo,   hfrb
  real(8) :: ax1max
  real(8) ::    phi,   hpnd,   perm,   prhd, dvperm
  real(8) :: delfmp,   cefb
  real(8) ::   etan,   etar, etanrr,    gil,    gir
  real(8) ::     x0,     x1,   gint
  real(8) ::   fahi,   fahs,   faei,   pvol
  real(8) ::   faas,   favm,   fafl,   fafm,   fadd,   fadb
  real(8) ::   eieq,  danew
  integer ::     ij,      k
  integer ::  ifpar,  jfpar,  istat

  real(8) :: hmp, dhmp
  logical :: iscrmp

  real(8), save ::    eps = 1.0d-3,   epsl = 1.0d-6

  if (oinit .or. ofinal) then
     return
  end if

  if (ofirst) then
     ofirst = .false.
     call rewnml(ifpar, jfpar)
     read (ifpar, nmamin, iostat=istat)
     call cstnml(jfpar, 'ptherm', 'nmamin', istat)
     write(jfpar, nmamin)
     call rewnml(ifpar, jfpar)
     read (ifpar, nmislt, iostat=istat)
     call cstnml(jfpar, 'ptherm', 'nmislt', istat)
     write(jfpar, nmislt)
     call rewnml(ifpar, jfpar)
     read (ifpar, nmmpnd, iostat=istat)
     call cstnml(jfpar, 'ptherm', 'nmmpnd', istat)
     write(jfpar, nmmpnd)
     call rewnml(ifpar, jfpar)
     read (ifpar, nmsage, iostat=istat)
     call cstnml(jfpar, 'ptherm', 'nmsage', istat)
     write(jfpar, nmsage)

     tmi = dtds * si
     rri    = rhoo / rhoi
     rrs    = rhoo / rhos
     rorirs = (rhoo - rhoi) / rhos
     rsfus  = rhos * hfus
!     epsaei = amin * abs(ei(tmi+eps, si))
     hic0(0) = 0.0d0
     do k = 1, nic+1
        hic0(k) = hic(k)
     end do
  end if

  do ij = 1, nxydim
     igrfra(ij) = 0.0d0
     igrcon(ij) = 0.0d0
     igrsni(ij) = 0.0d0
     igrsfl(ij) = 0.0d0
     inrlat(ij) = 0.0d0
     imrsno(ij) = 0.0d0
     imrsmi(ij) = 0.0d0
     imrisf(ij) = 0.0d0
     imribs(ij) = 0.0d0
  end do

  do k = 0, nic
     do ij = 1, nxydim
        az(ij, k) = ax(ij, k)
        hiz(ij, k) = hix(ij, k)
        vmpz(ij, k) = vmpx(ij, k)
        aflrm(ij, k) = 0.0d0
        daxhit(ij, k) = 0.0d0
        daxhib(ij, k) = 0.0d0
        dvmp(ij, k) = 0.0d0
        impinc(ij, k) = 0.0d0
        improf(ij, k) = 0.0d0
        impfrz(ij, k) = 0.0d0
     end do
  end do

  do k = 0, nic
     do ij = 1, nxydim
        axhix(ij, k) = ax(ij, k) * hix(ij, k)
        axhsx(ij, k) = ax(ij, k) * hsx(ij, k)
        axeix(ij, k) = ax(ij, k) * eix(ij, k)
        axflv(ij, k) = ax(ij, k) * frlvx(ij, k)
        axvmp(ij, k) = ax(ij, k) * vmpx(ij, k)
!        axfmp(ij, k) = ax(ij, k) * frmpx(ij, k)
        axdsd(ij, k) = ax(ij, k) * dsdx(ij, k)
        axdsb(ij, k) = ax(ij, k) * dsbx(ij, k)
     end do
  end do

  do ij = 1, nxydim
     rmpcc(ij) = rmpcmn(impnd)
  end do

  do k = 1, nic
     do ij = 1, nxydim
        rmpcc(ij) = rmpcc(ij) + &
          &         (rmpcmx(impnd) - rmpcmn(impnd)) * ax(ij, k)
!        impth2(ij) = impth2(ij) - ax(ij, k) * vmpx(ij, k)
     end do
  end do

! *** snowfall ***
  do k = 1, nic
     do ij = ijtstr, ijtend
        if (ax(ij, k) .gt. 0.d0) then
           hsx(ij, k) = axhsx(ij, k) / ax(ij, k) &
             &        + ts * rrs * snow(ij)
           dsdx(ij, k) = axdsd(ij, k) / ax(ij, k) &
             &         + ts * dfdu(ij)
           dsbx(ij, k) = axdsb(ij, k) / ax(ij, k) &
             &         + ts * dfbc(ij)
           axhsx(ij, k) = ax(ij, k) * hsx(ij, k)
           axdsd(ij, k) = ax(ij, k) * dsdx(ij, k)
           axdsb(ij, k) = ax(ij, k) * dsbx(ij, k)
           igrsfl(ij) = igrsfl(ij) + ax(ij, k) * snow(ij) * ts
        else
           hsx(ij, k) = 0.d0
           axhsx(ij, k) = 0.d0
           asx(ij, k) = 0.d0
           axdsd(ij, k) = 0.d0
           dsdx(ij, k) = 0.d0
           axdsb(ij, k) = 0.d0
           dsbx(ij, k) = 0.d0
        end if
        if (hsx(ij, k) .le. 0.d0) then
           dsdrhs(ij, k) = 0.d0
           dsbrhs(ij, k) = 0.d0
        else
           dsdrhs(ij, k) = min(dsdx(ij, k) / hsx(ij, k), drsmax)
           dsbrhs(ij, k) = min(dsbx(ij, k) / hsx(ij, k), drsmax)
        end if
     end do
  end do

! *** catching rainfall ***
  do k=1, nic
     do ij = ijtstr, ijtend
        impinc(ij, k) = impinc(ij, k) + &
          &           rmpcc(ij) * ax(ij, k) * prec(ij) * ts
     end do
  end do

  do ij = ijtstr, ijtend
     snow(ij) = ax(ij, 0) * snow(ij)
     prec(ij) = prec(ij) + snow(ij)
     fdd(ij) = fdd(ij) - ax(ij, 0) * dfdu(ij)
     fdb(ij) = fdb(ij) - ax(ij, 0) * dfbc(ij)
  end do

! *** snow melting ***
  do k = 1, nic
     do ij = ijtstr, ijtend
        wres = axhsx(ij, k) * rsfus / ts + was(ij, k)
        if (ax(ij, k) .gt. 0.d0) then
           if (wres .lt. 0.d0) then
              hsx(ij, k) = 0.d0
              wai(ij, k) = wres
              asx(ij, k) = 0.d0
              dsdx(ij, k) = 0.d0
              dsbx(ij, k) = 0.d0
           else
              hsx(ij, k) = wres * ts / rsfus / ax(ij, k)
              wai(ij, k) = 0.d0
              dsdx(ij, k) = dsdrhs(ij, k) * hsx(ij, k)
              dsbx(ij, k) = dsbrhs(ij, k) * hsx(ij, k)
           end if
        else
           wai(ij, k) = was(ij, k)
        end if
        axhsxn(ij, k) = ax(ij, k) * hsx(ij, k) * amskt(ij, kstr)
        imrsno(ij) = imrsno(ij) - &
          &          ( axhsxn(ij, k) - axhsx(ij, k) ) * rhos
        impinc(ij, k) = impinc(ij, k) - rmpcc(ij) * rhos * &
          &          min((axhsxn(ij, k) - axhsx(ij, k)), 0.0d0)
     end do
  end do

! *** ice top melting ***
  do k = 1, nic
     do ij = ijtstr, ijtend
        wres = rhoi * axeix(ij, k) / ts + wai(ij, k)
        if (ax(ij, k) .gt. 0.d0) then
           if (wres .lt. 0.d0) then
              axeixn(ij, k) = 0.d0
              wio(ij, k) = wio(ij, k) + wres
           else
              axeixn(ij, k) = wres / rhoi * ts
           end if
           daxhit(ij, k) = &
             &   ( axeixn(ij, k) - axeix(ij, k) ) / ei(tix(ij, k), si)
        else
           wio(ij, k) = wio(ij, k) + wai(ij, k)
           axeixn(ij, k) = axeix(ij, k)
           daxhit(ij, k) = 0.d0
        end if
        impinc(ij, k) = impinc(ij, k) &
          &          - min((rmpcc(ij) * rhoi * daxhit(ij, k)), 0.0d0)
     end do
  end do

! *** new ice formation on open water ***
  do ij = ijtstr, ijtend
     axeixn(ij, 0) = wao(ij) * ts / rhoi * amskt(ij, kstr)
     axhsxn(ij, 0) = 0.d0
     tix(ij, 0) = min(tx(ij, kstr, 1)*amskt(ij, kstr), tmi)
     if (      (ax(ij, 0) .gt. 0.d0) &
       & .and. (axeixn(ij, 0) .gt. 0.d0) &
!       & .and. (tix(ij, 0) .lt. tmi)) then
       & .and. (tix(ij, 0) .lt. tmi-eps)) then
        eix(ij, 0) = axeixn(ij, 0) / ax(ij, 0)
        hix(ij, 0) = eix(ij, 0) / ei(tix(ij, 0), si)
        igrfra(ij) = ax(ij, 0) * hix(ij, 0) * rhoi
     else
        eix(ij, 0) = 0.d0
        hix(ij, 0) = 0.d0
     end if
  end do

! *** basal and lateral ice formation/melting processes are divided
! *** in order to apply linear-remapping method
! *** basal ice formation/melting
  do k = 1, nic
     do ij = ijtstr, ijtend
        daxhib(ij, k) = axeixn(ij, k)
        axeixn(ij, k) = axeixn(ij, k) &
          &           + wio(ij, k) * ts / rhoi &
          &           * amskt(ij, kstr)
     end do
  end do
  do k = 1, nic
     do ij = ijtstr, ijtend
        if (axeixn(ij, k) .le. 0.d0 .or. ax(ij, k) .le. 0.d0) then
           if (daxhib(ij, k) .le. 0.d0) then
              daxhib(ij, k) = 0.d0
           else
              daxhib(ij, k) = &
                &   - daxhib(ij, k) / ei(tix(ij, k), si)
           end if
           ax(ij, k) = 0.d0
           hix(ij, k) = hic(k)
           hsx(ij, k) = 0.d0
           tix(ij, k) = tmi
           eix(ij, k) = 0.d0
           asx(ij, k) = 0.d0
           frlvx(ij, k) = 1.d0
           vmpx(ij, k) = 0.d0
           frmpx(ij, k) = 0.d0
           dsdx(ij, k) = 0.d0
           dsbx(ij, k) = 0.d0
        else
           daxhib(ij, k) = &
             &     ( axeixn(ij, k) - daxhib(ij, k) ) &
             &     / ei(tix(ij, k), si)
           eix(ij, k) = axeixn(ij, k) / ax(ij, k)
           hix(ij, k) = eix(ij, k) / ei(tix(ij, k), si)
           vmpx(ij, k) = (axvmp(ij, k) + impinc(ij, k)) / ax(ij, k)
        end if
     end do
  end do

  do k = 1, nic
     do ij = ijtstr, ijtend
        if (az(ij, k) .gt. 0.0d0) then
           imrisf(ij) = imrisf(ij) - rhoi * daxhit(ij, k)
           igrcon(ij) = igrcon(ij) &
             &        + rhoi * max(0.0d0, daxhib(ij, k))
           imribs(ij) = imribs(ij) &
             &        - rhoi * min(0.0d0, daxhib(ij, k))
           imrsmi(ij) = imrsmi(ij) - &
             &          rhos * (ax(ij, k) * hsx(ij, k) - axhsxn(ij, k))
        end if
     end do
  end do

! *** snow-ice formation ***
  do k = 1, nic
     do ij = ijtstr, ijtend
        hsxo = hsx(ij, k)
        hsx(ij, k) = min(hsx(ij, k), rorirs * hix(ij, k))
        dhs = hsxo - hsx(ij, k)
        hix(ij, k) = hix(ij, k) + dhs * rhos / rhoi
        eix(ij, k) = eix(ij, k) + rsfus / rhoi * dhs
        dsdx(ij, k) = dsdrhs(ij, k) * hsx(ij, k)
        dsbx(ij, k) = dsbrhs(ij, k) * hsx(ij, k)
        igrsni(ij) = igrsni(ij) + ax(ij, k) * dhs * rhos
        if (ax(ij, k) .gt. 0.d0) then
           tix(ij, k) = ti(eix(ij, k)/hix(ij, k), si)
        else
           tix(ij, k) = tmi
        end if
     end do
  end do

! *** melt pond freezing ***
  do k = 1, nic
     do ij = ijstr, ijend
        if (ax(ij, k) .gt. 0.0d0) then
           vmpo = vmpx(ij, k)
!          We make the exponent depend on dt, unlike the CICE implementation.
           vmpx(ij, k) = vmpx(ij, k) * &
             &    exp( cmpfrz * dt * &
             &         max(tmpfrz-tsi(ij, k), 0.0d0) / tmpfrz )
           impfrz(ij, k) = - ax(ij, k) * ( vmpx(ij, k) - vmpo )
!          Save the change in vmpx
           dvmp(ij, k) = max( vmpx(ij, k)-vmpz(ij, k), -vmpz(ij, k) )
         end if
     end do
  end do

! *** negative freeboard consideration (virtual) ***
! Runoff here does not change frmpx, following the CICE implementation.
  do k = 1, nic
     do ij = ijtstr, ijtend
        if (frmpx(ij, k) > 0.0d0) then
        vmpo = vmpx(ij, k)
!        CICE uses meltpond depth (= vmpx/frmpx) instead of vmpx as below
!        hfrb = min(vmpx(ij, k) / frmpx(ij, k), &
!          &    ( (rhoo - rhoi) * hix(ij, k) - rhos * hsx(ij,k) ) / rhoo )
!        vmpx(ij, k) = hfrb * frmpx(ij, k)
        vmpx(ij, k) = min(vmpx(ij, k), &
          &    ( (rhoo - rhoi) * hix(ij, k) - rhos * hsx(ij,k) ) / rhoo )
        improf(ij, k) = improf(ij, k) - ax(ij, k) * (vmpx(ij, k) - vmpo)
        end if
     end do
  end do

! *** permiability ***
  if (impnd == 2) then  !! Hunke MP param.
     do k = 1, nic
        do ij = ijstr, ijend
           if (frmpx(ij, k) > 0.0d0) then
              phi = si * (1.0d-3 - 0.054d0 / tix(ij, k))
              if (phi >= 0.05d0) then  !! permiable ice
                 hpnd = max(vmpx(ij, k) / frmpx(ij, k), 0.0d0)
!                In Hunke et al. (2013), perm is proportional to phi**1,
!                 but we follow the CICE Icepack code (ver. 1.2.5).
                 perm = 3.0d-4 * (phi**3)     !! permiability [cm**2]
                 prhd = gravit * rhoo * &
                   &  ( (1.0d0 - 1.0d0/rri) * hix(ij, k) &
                   &   - 1.0d0/rrs * hsx(ij, k) )
                 dvperm = - frmpx(ij, k) * &
                   &      min( hpnd, &
                   &           dpscl * perm * prhd * ts / 1.79d-2 / hix(ij, k) )
                 dvperm = max(-vmpx(ij, k), dvperm)
                 vmpx(ij, k) = vmpx(ij, k) + dvperm
                 dvmp(ij, k) = max( dvmp(ij, k)+dvperm, -vmpz(ij, k) )
                 improf(ij, k) = improf(ij, k) - ax(ij, k) * dvperm
              end if
           end if
        end do
     end do
  end if

! *** update fraction of meltpond ***
  if (impnd == 2) then  !! Hunke MP param.
     do k = 1, nic
        do ij = ijstr, ijend
           if (ax(ij, k) > 0.0d0) then
              if (vmpx(ij, k) <= 0.0d0) then
                 vmpx(ij, k) = 0.0d0
                 frmpx(ij, k) = 0.0d0
              else
!                We solve the following equations for dfmp
!                 under the condition that dvmp >= -vmp:
!                   vmpz + dvmp = (frmpx + delfmp) * (hmp + delhmp)
!                   delhmp = rtdpmp * delfmp,
!                 where hmp is meltpond depth and delhmp is its variation.
                 if (frmpx(ij, k) <= 0.0d0) then  !! MP not exist
                    delfmp = sqrt(max(dvmp(ij, k), 0.0d0)/rtdpmp)
                    frmpx(ij, k) = delfmp
                    iscrmp = .true.
                    hmp = 0.0d0
                    if (delfmp > 0.0d0) then
                       dhmp = dvmp(ij, k)/delfmp
                    else
                       dhmp = 0.0d0
                    end if
                 else
                    cefb = vmpz(ij, k) / frmpx(ij, k) + rtdpmp * frmpx(ij, k)
!                   max function is for avoiding floating invalid due to numerical error
                    delfmp = ( -cefb &
                      &      + sqrt( max( cefb**2 + 4.0d0 * rtdpmp * dvmp(ij, k), 0.0d0 ) ) )&
                      &      / (2.0d0 * rtdpmp)                    
                    hmp = vmpz(ij, k) / frmpx(ij, k)
                    if ((frmpx(ij, k)+delfmp) <= 0.d0) then
                       dhmp = -hmp
                    else
                       dhmp = (vmpz(ij, k)+dvmp(ij, k)) / (frmpx(ij, k)+delfmp) - hmp
                    end if
                    frmpx(ij, k) = frmpx(ij, k) + delfmp
                    iscrmp = .false.
                 end if
                 if (frmpx(ij, k) <= 0.0d0) then
                    improf(ij, k) = improf(ij, k) &
                      &           + ax(ij, k) * max(vmpx(ij, k), 0.0d0)        
                    frmpx(ij, k) = 0.0d0
                    vmpx(ij, k) = 0.0d0
                 else if (frmpx(ij, k) > frlvx(ij, k)) then !! MP water runoff
                    vmpo = vmpx(ij, k)
                    vmpx(ij, k) = vmpx(ij, k) * frlvx(ij, k) / frmpx(ij, k)
                    frmpx(ij, k) = frlvx(ij, k)
                    improf(ij, k) = improf(ij, k) &
                      &           - ax(ij, k) * (vmpx(ij, k) - vmpo)        
                 end if
              end if
           end if
        end do
     end do
  end if

! ****** linear remapping of Lipscomb(2001)
! *** setting flags
  do ij = ijtstr, ijtend
     aflrmc(ij) = amskt(ij, kstr)
     aflrm(ij, 0) = 1.0d0
  end do
 
  do k = 1, nic
     do ij = ijtstr, ijtend
        if ( ( ax(ij, k) .gt. 0.0d0 ) .and. &
          &  ( az(ij, k) .gt. 0.0d0 ) ) then
           aflrm(ij, k) = 1.0d0
        end if
     end do
  end do

! *** growth rate of each categories and category boundaries
  do ij = ijtstr, ijtend
     fdtn(ij, 0) = hix(ij, 0)
  end do
  do k = 1, nic
     do ij = ijtstr, ijtend
        fdtn(ij, k) = ( hix(ij, k) - hiz(ij, k) ) * aflrm(ij, k)
     end do
  end do

  do ij = ijtstr, ijtend
     if ( fdtn(ij, 0) .gt. 0.0d0 ) then
        fdtcn(ij, 1) = fdtn(ij, 0)
     elseif ( aflrm(ij, 1) .eq. 1.0d0 ) then
        fdtcn(ij, 1) = fdtn(ij, 1)
     else
        fdtcn(ij, 1) = 0.0d0
     endif
  end do
  do k = 2, nic
     do ij = ijtstr, ijtend
        if ( ( aflrm(ij, k-1) .eq. 1.0d0 ) &
          &  .and.( aflrm(ij, k) .eq. 1.0d0 ) ) then
           fdtcn(ij, k) = fdtn(ij, k-1) + &
             &          ( fdtn(ij, k) - fdtn(ij, k-1) ) * &
             &          ( hic0(k) - hiz(ij, k-1) ) &
             &          / ( hiz(ij, k) - hiz(ij, k-1) )
        elseif (aflrm(ij, k-1) .eq. 1.0d0) then
           fdtcn(ij, k) = fdtn(ij, k-1)
        elseif (aflrm(ij, k) .eq. 1.0d0) then
           fdtcn(ij, k) = fdtn(ij, k)
        else
           fdtcn(ij, k) = 0.0d0
        endif
     end do
  end do

!  call chekin(  fdtn,  'FDTN', &
!    &             nx,      ny,    nic, nxyidm, 'ICE')
!  call chekin( fdtcn, 'FDTCN', &
!    &             nx,      ny,    nic, nxyidm, 'ICE')
!  call chekin(   hiz,   'HIZ', &
!    &             nx,      ny,    nic, nxyidm, 'ICE')
!  call chekin(   hix,  'HIX1', &
!    &             nx,      ny,    nic, nxyidm, 'ICE')

! *** temporally shift category boundaries
  do k = 1, nic
     do ij = ijtstr, ijtend
        hicn(ij, k) = hic(k) + fdtcn(ij,k)
     end do
  end do

! *** validation check: will not execute remapping when...
! ***  - hicn(ij,k) does not lie between hix(ij,k-1) and hix(ij,k)
! ***  - hicn(ij,k) does not lie between hic0(k-1) and hic0(k+1)
  do k = 1, nic
     do ij = ijtstr, ijtend
        hicn(ij, k) = hic(k) + fdtcn(ij,k)
        if ( (aflrm(ij, k-1) .eq. 1.0d0) .and. &
          &  (hicn(ij, k) .le. hix(ij, k-1)) ) then
           aflrmc(ij) = 0.0d0
        endif
        if ( (aflrm(ij, k) .eq. 1.0d0) .and. &
          &  (hicn(ij, k) .ge. hix(ij, k)) ) then
           aflrmc(ij) = 0.0d0
        endif
        if ( (hicn(ij, k) .le. hic0(k-1)) .or. &
          &  (hicn(ij, k) .ge. hic0(k+1)) ) then
           aflrmc(ij) = 0.0d0
        endif
     end do
  end do

  do k = 0, nic
     do ij = ijtstr, ijtend
        aflrm(ij, k) = aflrm(ij, k) * aflrmc(ij)
     end do
  end do

!  call chekin(   hicn,  'HICN', &
!    &              nx,      ny,    nic, nxyidm, 'ICE')

! *** determine linear distribution function within each category
  do k = 1, nic-1
     do ij = ijtstr, ijtend
        hil(ij, k) = max( hicn(ij, k), &
          &               3.0d0*hix(ij, k) - 2.0d0*hicn(ij, k+1) )
        hir(ij, k) = min( hicn(ij, k+1), &
          &               3.0d0*hix(ij, k) - 2.0d0*hicn(ij, k) )
        etan = hix(ij, k) - hil(ij, k)
!       function max would not applied when ax(ij, k) > 0
        etar = max(hir(ij, k) - hil(ij, k), epsl)
        etanrr = etan / etar
        g1(ij, k) = 12.0d0 * ax(ij, k) / (etar*etar) &
          &         * ( etanrr - 0.5d0) &
          &         * aflrm(ij, k)
        g0(ij, k) = 6.0d0 * ax(ij, k) / etar &
          &         * ( 2.0d0 / 3.0d0 - etanrr) &
          &         * aflrm(ij, k)
     end do
  end do
  do ij = ijtstr, ijtend
     hil(ij, nic) = hicn(ij, nic)
     hir(ij, nic) = 3.0d0*hix(ij, nic) - 2.0d0*hicn(ij, nic)
!     hir(ij, nic) = max(3.0d0*hix(ij, nic) - 2.0d0*hicn(ij, nic), &
!       &                hic(nic))
     etan = hix(ij, nic) - hil(ij, nic)
!    function max would not applied when ax(ij, nic) > 0
     etar = max(hir(ij, nic) - hil(ij, nic), epsl)
     etanrr = etan / etar
     g1(ij, nic) = 12.0d0 * ax(ij, nic) / (etar*etar) &
       &           * ( etanrr - 0.5d0 ) &
       &           * aflrm(ij, nic)
     g0(ij, nic) = 6.0d0 * ax(ij, nic) / etar &
       &           * ( 2.0d0 / 3.0d0 - etanrr) &
       &           * aflrm(ij, nic)
  end do

!  call chekin(    hil,  'HIL', &
!    &              nx,     ny,    nic, nxyidm, 'ICE')
!  call chekin(    hir,  'HIR', &
!    &              nx,     ny,    nic, nxyidm, 'ICE')
!  call chekin(     g0,   'G0', &
!    &              nx,     ny,    nic, nxyidm, 'ICE')
!  call chekin(     g1,   'G1', &
!    &              nx,     ny,    nic, nxyidm, 'ICE')

! *** compute transfer fluxes
! *** area and volume fluxes are calculated from the linear distribution
! *** snow and enthalpy fluxes are propotional to the volume flux
  do k = 1, nic
     do ij = ijtstr, ijtend
        laxhix(ij, k) = ax(ij, k) * hix(ij, k)
        laxhsx(ij, k) = ax(ij, k) * hsx(ij, k)
        laxeix(ij, k) = ax(ij, k) * eix(ij, k)
        laxasx(ij, k) = ax(ij, k) * asx(ij, k)
        laxflv(ij, k) = ax(ij, k) * frlvx(ij, k)
        laxvmp(ij, k) = ax(ij, k) * vmpx(ij, k)
        laxfmp(ij, k) = ax(ij, k) * frmpx(ij, k)
        laxdsd(ij, k) = ax(ij, k) * dsdx(ij, k)
        laxdsb(ij, k) = ax(ij, k) * dsbx(ij, k)
        da(ij, k) = 0.0d0
        dahi(ij, k) = 0.0d0
        dahs(ij, k) = 0.0d0
        daei(ij, k) = 0.0d0
        daas(ij, k) = 0.0d0
        dafl(ij, k) = 0.0d0
        davm(ij, k) = 0.0d0
        dafm(ij, k) = 0.0d0
        dadd(ij, k) = 0.0d0
        dadb(ij, k) = 0.0d0
     end do
  end do

  do ij = ijtstr, ijtend
     da(ij, 0) = 0.0d0
!    for the flux across the lowest boundary,
!    consider only melting situation at the boundary
     if ( ( hicn(ij, 1) .lt. hic(1) ) .and. &
       &  ( aflrm(ij, 1) .eq. 1.0d0 ) ) then
!        if ( hicn(ij, 1) .lt. hic(1) ) then
        if (hir(ij, 1) .le. hic(1)) then
!          all sea ice transfered out of the lowest category
           gint = ax(ij, 1)
        else
!           gil = hil(ij, 1) - hil(ij, 1)
           gil = 0.0d0
           gir = max(hic(1) - hil(ij, 1), gil)
           x0 = gir - gil
           x1 = 0.5d0 * (gir*gir - gil*gil)
           gint = g0(ij, 1)*x0 + g1(ij, 1)*x1
        end if
!       area loss does not cause sea-ice to thicken
        gint = max( &
               min( gint, &
          &         ax(ij, 1)*(1.0d0 - hix(ij, 1)/hiz(ij, 1)) ), &
          &         0.0d0 )
!       only area flux can across the lowest boundary
        da(ij, 1) = da(ij, 1) - gint
        da(ij, 0) = da(ij, 0) + gint
!       snow age, level-ice frac., and melt-pond frac. do not change
        faas = asx(ij, 1) * gint
        fafl = frlvx(ij, 1) * gint
        fafm = frmpx(ij, 1) * gint
        daas(ij, 1) = daas(ij, 1) - faas
        dafl(ij, 1) = dafl(ij, 1) - fafl
        dafm(ij, 1) = dafm(ij, 1) - fafm
     end if
  end do

  do k = 2, nic
     do ij = ijtstr, ijtend
        if ( ( hicn(ij, k) .ge. hic(k) ) .and. &
          &  ( aflrm(ij, k-1) .eq. 1.0d0 ) ) then
!           if ( hicn(ij, k) .ge. hic(k) ) then
           if (hil(ij, k-1) .ge. hic(k)) then
!             all sea ice transfered to the upper category
              gint = ax(ij, k-1)
              fahi = laxhix(ij, k-1)
              pvol = 1.0d0
              fahs = laxhsx(ij, k-1)
              faei = laxeix(ij, k-1)
              faas = laxasx(ij, k-1)
              fafl = laxflv(ij, k-1)
              favm = laxvmp(ij, k-1)
              fafm = laxfmp(ij, k-1)
              fadd = laxdsd(ij, k-1)
              fadb = laxdsb(ij, k-1)
           else
              gil = hic(k) - hil(ij, k-1)
              gir = max(hir(ij, k-1) - hil(ij, k-1), gil)
              x0 = gir - gil
              x1 = 0.5d0 * (gir*gir - gil*gil)
              gint = g0(ij, k-1)*x0 + g1(ij, k-1)*x1
              gint = max(min(gint, ax(ij, k-1)), 0.0d0)
              x0 = 0.5d0 * (gir*gir - gil*gil)
              x1 = (gir*gir*gir - gil*gil*gil) / 3.0d0
              fahi = hil(ij, k-1) * gint + &
                &    g0(ij, k-1)*x0 + g1(ij, k-1)*x1
              fahi = max(min(fahi, laxhix(ij, k-1)), 0.0d0)
              pvol = fahi / laxhix(ij, k-1)
              fahs = laxhsx(ij, k-1) * pvol
              faas = asx(ij, k-1) * gint
              fafl = frlvx(ij, k-1) * gint
              favm = vmpx(ij, k-1) * gint
              fafm = frmpx(ij, k-1) * gint
              fadd = laxdsd(ij, k-1) * pvol
              fadb = laxdsb(ij, k-1) * pvol
              if (tix(ij, k-1) .lt. tmi) then
                 faei = laxeix(ij, k-1) * pvol
!                 faei = fahi * ei(tix(ij, k-1), si)
              else
                 faei = 0.d0
              end if
!              faei = ax(ij, k-1) * eix(ij, k-1) * pvol
           end if
           da(ij, k-1) = da(ij, k-1) - gint
           dahi(ij, k-1) = dahi(ij, k-1) - fahi
           dahs(ij, k-1) = dahs(ij, k-1) - fahs
           daei(ij, k-1) = daei(ij, k-1) - faei
           daas(ij, k-1) = daas(ij, k-1) - faas
           dafl(ij, k-1) = dafl(ij, k-1) - fafl
           davm(ij, k-1) = davm(ij, k-1) - favm
           dafm(ij, k-1) = dafm(ij, k-1) - fafm
           dadd(ij, k-1) = dadd(ij, k-1) - fadd
           dadb(ij, k-1) = dadb(ij, k-1) - fadb
           da(ij, k) = da(ij, k) + gint
           dahi(ij, k) = dahi(ij, k) + fahi
           dahs(ij, k) = dahs(ij, k) + fahs
           daei(ij, k) = daei(ij, k) + faei
           daas(ij, k) = daas(ij, k) + faas
           dafl(ij, k) = dafl(ij, k) + fafl
           davm(ij, k) = davm(ij, k) + favm
           dafm(ij, k) = dafm(ij, k) + fafm
           dadd(ij, k) = dadd(ij, k) + fadd
           dadb(ij, k) = dadb(ij, k) + fadb
        elseif ( ( hicn(ij, k) .lt. hic(k) ) .and. &
          &      ( aflrm(ij, k) .eq. 1.0d0 ) ) then
           if (hir(ij, k) .le. hic(k)) then
!             all sea ice transfered to the lower category
              gint = ax(ij, k)
              fahi = laxhix(ij, k)
              pvol = 1.0d0
              fahs = laxhsx(ij, k)
              faei = laxeix(ij, k)
              faas = laxasx(ij, k)
              fafl = laxflv(ij, k)
              favm = laxvmp(ij, k)
              fafm = laxfmp(ij, k)
              fadd = laxdsd(ij, k)
              fadb = laxdsb(ij, k)
           else
!              gil = hil(ij, k) - hil(ij, k)
              gil = 0.0d0
              gir = max(hic(k) - hil(ij, k), gil)
              x0 = gir - gil
              x1 = 0.5d0 * (gir*gir - gil*gil)
              gint = g0(ij, k)*x0 + g1(ij, k)*x1
              gint = max(min(gint, ax(ij, k)), 0.0d0)
              x0 = 0.5d0 * (gir*gir - gil*gil)
              x1 = (gir*gir*gir - gil*gil*gil) / 3.0d0
              fahi = hil(ij, k) * gint + &
                &    g0(ij, k)*x0 + g1(ij, k)*x1
              fahi = max(min(fahi, laxhix(ij, k)), 0.0d0)
              pvol = fahi / laxhix(ij, k)
              fahs = laxhsx(ij, k) * pvol
              faas = asx(ij, k) * gint
              fafl = frlvx(ij, k) * gint
              favm = vmpx(ij, k) * gint
              fafm = frmpx(ij, k) * gint
              fadd = laxdsd(ij, k) * pvol
              fadb = laxdsb(ij, k) * pvol
              if (tix(ij, k) .lt. tmi) then
                 faei = laxeix(ij, k) * pvol
!                 faei = fahi * ei(tix(ij, k), si)
              else
                 faei = 0.d0
              end if
!              faei = ax(ij, k) * eix(ij, k) * pvol
           endif
           da(ij, k-1) = da(ij, k-1) + gint
           dahi(ij, k-1) = dahi(ij, k-1) + fahi
           dahs(ij, k-1) = dahs(ij, k-1) + fahs
           daei(ij, k-1) = daei(ij, k-1) + faei
           daas(ij, k-1) = daas(ij, k-1) + faas
           dafl(ij, k-1) = dafl(ij, k-1) + fafl
           davm(ij, k-1) = davm(ij, k-1) + favm
           dafm(ij, k-1) = dafm(ij, k-1) + fafm
           dadd(ij, k-1) = dadd(ij, k-1) + fadd
           dadb(ij, k-1) = dadb(ij, k-1) + fadb
           da(ij, k) = da(ij, k) - gint
           dahi(ij, k) = dahi(ij, k) - fahi
           dahs(ij, k) = dahs(ij, k) - fahs
           daei(ij, k) = daei(ij, k) - faei
           daas(ij, k) = daas(ij, k) - faas
           dafl(ij, k) = dafl(ij, k) - fafl
           davm(ij, k) = davm(ij, k) - favm
           dafm(ij, k) = dafm(ij, k) - fafm
           dadd(ij, k) = dadd(ij, k) - fadd
           dadb(ij, k) = dadb(ij, k) - fadb
        end if
     end do
  end do

!  call chekin(     da,   'DA', &
!    &              nx,     ny,    nic, nxyidm, 'ICE')
!  call chekin(   dahi, 'DAHI', &
!    &              nx,     ny,    nic, nxyidm, 'ICE')
!  call chekin(   dahs, 'DAHS', &
!    &              nx,     ny,    nic, nxyidm, 'ICE')
!  call chekin(   daei, 'DAEI', &
!    &              nx,     ny,    nic, nxyidm, 'ICE')

! *** update prediction variables.
  do k = 1, nic
     do ij = ijtstr, ijtend
        ax(ij, k) = ax(ij, k) + da(ij, k) * amskt(ij, kstr)
        laxhix(ij, k) = laxhix(ij, k) + &
          &               dahi(ij, k) * amskt(ij, kstr)
        laxhsx(ij, k) = laxhsx(ij, k) + &
          &               dahs(ij, k) * amskt(ij, kstr)
        laxeix(ij, k) = laxeix(ij, k) + &
          &               daei(ij, k) * amskt(ij, kstr)
        laxasx(ij, k) = laxasx(ij, k) + &
          &               daas(ij, k) * amskt(ij, kstr)
        laxflv(ij, k) = laxflv(ij, k) + &
          &               dafl(ij, k) * amskt(ij, kstr)
        laxvmp(ij, k) = laxvmp(ij, k) + &
          &               davm(ij, k) * amskt(ij, kstr)
        laxfmp(ij, k) = laxfmp(ij, k) + &
          &               dafm(ij, k) * amskt(ij, kstr)
        laxdsd(ij, k) = laxdsd(ij, k) + &
          &               dadd(ij, k) * amskt(ij, kstr)
        laxdsb(ij, k) = laxdsb(ij, k) + &
          &               dadb(ij, k) * amskt(ij, kstr)
        axhsxn(ij, k) = axhsxn(ij, k) + &
          &               dahs(ij, k) * amskt(ij, kstr)
        axeixn(ij, k) = axeixn(ij, k) + &
          &               daei(ij, k) * amskt(ij, kstr)
        if (laxeix(ij, k) .le. 0.d0) then
           ax(ij, k) = 0.d0
           hix(ij, k) = hic(k)
           hsx(ij, k) = 0.d0
           tix(ij, k) = tmi
           eix(ij, k) = 0.d0
           asx(ij, k) = 0.d0
           frlvx(ij, k) = 1.d0
           vmpx(ij, k) = 0.d0
           frmpx(ij, k) = 0.d0
           dsdx(ij, k) = 0.d0
           dsbx(ij, k) = 0.d0
        elseif (ax(ij, k) .le. 0.d0) then
           ax(ij, k) = 0.d0
           hix(ij, k) = hic(k)
           hsx(ij, k) = 0.d0
           tix(ij, k) = tmi
           eix(ij, k) = 0.d0
           asx(ij, k) = 0.d0
           frlvx(ij, k) = 1.d0
           vmpx(ij, k) = 0.d0
           frmpx(ij, k) = 0.d0
           dsdx(ij, k) = 0.d0
           dsbx(ij, k) = 0.d0
        else
           hix(ij, k) = laxhix(ij, k) / ax(ij, k)
           hsx(ij, k) = laxhsx(ij, k) / ax(ij, k)
           eix(ij, k) = laxeix(ij, k) / ax(ij, k)
           asx(ij, k) = laxasx(ij, k) / ax(ij, k)
           frlvx(ij, k) = laxflv(ij, k) / ax(ij, k)
           vmpx(ij, k) = laxvmp(ij, k) / ax(ij, k)
           frmpx(ij, k) = laxfmp(ij, k) / ax(ij, k)
!           tix(ij, k) = ti(eix(ij, k)/hix(ij, k), si)
           dsdx(ij, k) = laxdsd(ij, k) / ax(ij, k)
           dsbx(ij, k) = laxdsb(ij, k) / ax(ij, k)
        end if
!       Some variables can be out of their valid range due to truncation error
        asx(ij, k) = max(0.0d0, min(1.0d0, asx(ij, k)))
        frlvx(ij, k) = max(0.0d0, min(1.0d0, frlvx(ij, k)))
        frmpx(ij, k) = max(0.0d0, min(1.0d0, frmpx(ij, k)))
     end do
  end do

! *** check if variables are in valid range
!  do k = 1, nic
!     do ij = ijtstr, ijtend
!        if (ax(ij, k) .lt. 0.d0) then
!           call rewnml(ifpar, jfpar)
!           write(jfpar, *) '### NEGATIVE AREA (linrmp) ###', &
!             &             ij, k, ax(ij, k)
!           stop
!        end if
!!        if (     (laxhix(ij, k) .lt. 0.d0) &
!!          & .or. (axhsxn(ij, k) .lt. 0.d0)) then
!!           call rewnml(ifpar, jfpar)
!!           write(jfpar, *) '### NEGATIVE ICE/SNOW (linrmp) ###', &
!!             &             ij, k, ax(ij, k), laxhix(ij, k), axhsxn(ij, k)
!!           stop
!!        end if
!!        if (axeixn(ij, k) .lt. 0.d0) then
!!           call rewnml(ifpar, jfpar)
!!           write(jfpar, *) '### NEGATIVE ENTHALPY (linrmp) ###', &
!!             &             ij, k, ax(ij, k), axeixn(ij, k)
!!           stop
!!        end if
!        if (hix(ij, k) .le. 0.d0) then
!           call rewnml(ifpar, jfpar)
!           write(jfpar, *) '### INVALID ICE THICKNESS (linrmp) ###', &
!             &             ij, k, ax(ij, k), hix(ij, k)
!           stop
!        end if
!        if ((frlvx(ij, k) .lt. 0.d0).or.(frlvx(ij, k) .gt. 1.d0)) then
!           call rewnml(ifpar, jfpar)
!           write(jfpar, *) '### INVALID LEVEL-ICE FRACTION (linrmp) ###', &
!             &             ij, k, ax(ij, k), frlvx(ij, k)
!           stop
!        end if
!        if ((frmpx(ij, k) .lt. 0.d0).or.(frmpx(ij, k) .gt. 1.d0)) then
!         call rewnml(ifpar, jfpar)
!         write(jfpar, *) '### INVALID MELT-POND FRACTION (linrmp) ###', &
!           &             ij, k, ax(ij, k), frmpx(ij, k)
!         stop
!      end if
!     end do
!  end do

! ****** end of linear remapping of Lipscomb(2001)

! *** lateral ice formation/melting
  do k = 1, nic
     do ij = ijtstr, ijtend
        laxhix(ij, k) = ax(ij, k) * hix(ij, k)
        laxeix(ij, k) = laxeix(ij, k) &
          &           + wil(ij, k) * ts / rhoi * &
          &             amskt(ij, kstr)
        axeixn(ij, k) = axeixn(ij, k) &
          &           + wil(ij, k) * ts / rhoi * &
          &             amskt(ij, kstr)
        laxflv(ij, k) = ax(ij, k) * frlvx(ij, k)
        laxvmp(ij, k) = ax(ij, k) * vmpx(ij, k)
        laxfmp(ij, k) = ax(ij, k) * frmpx(ij, k)
     end do
  end do
  do k = 1, nic
     do ij = ijtstr, ijtend
        if (ax(ij, k) .gt. 0.d0) then
           danew = ts * wil(ij, k) / eix(ij, k) / rhoi &
             &   * amskt(ij, kstr)
           ax(ij, k) = ax(ij, k) + danew
           if (danew >= 0.0d0) then
!             all the newly-formed ice is classed as level ice
!             laxvmp and laxfmp not changed (vmp and frmp decrease)
              laxflv(ij, k) = laxflv(ij, k) + danew * 1.0d0
           else
!             frlv, vmp, frmp not changed (lax* decrease)
              laxflv(ij, k) = ax(ij, k) * frlvx(ij, k)
              laxvmp(ij, k) = ax(ij, k) * vmpx(ij, k)
              laxfmp(ij, k) = ax(ij, k) * frmpx(ij, k)
           end if
        end if
     end do
  end do
  do k = 1, nic
     do ij = ijtstr, ijtend
        if (laxeix(ij, k) .le. 0.d0) then
           ax(ij, k) = 0.d0
           hix(ij, k) = hic(k)
           hsx(ij, k) = 0.d0
           tix(ij, k) = tmi
           eix(ij, k) = 0.d0
           asx(ij, k) = 0.d0
           frlvx(ij, k) = 1.d0
           vmpx(ij, k) = 0.d0
           frmpx(ij, k) = 0.d0
           dsdx(ij, k) = 0.d0
           dsbx(ij, k) = 0.d0
        else if (ax(ij, k) .le. 0.d0) then
           ax(ij, k) = 0.d0
           hix(ij, k) = hic(k)
           hsx(ij, k) = 0.d0
           tix(ij, k) = tmi
           eix(ij, k) = 0.d0
           asx(ij, k) = 0.d0
           frlvx(ij, k) = 1.d0
           vmpx(ij, k) = 0.d0
           frmpx(ij, k) = 0.d0
           dsdx(ij, k) = 0.d0
           dsbx(ij, k) = 0.d0
        else
           eix(ij, k) = laxeix(ij, k) / ax(ij, k)
           tix(ij, k) = ti(eix(ij, k)/hix(ij, k), si)
           frlvx(ij, k) = laxflv(ij, k) / ax(ij, k)
           vmpx(ij, k) = laxvmp(ij, k) / ax(ij, k)
           frmpx(ij, k) = laxfmp(ij, k) / ax(ij, k)
        end if
        inrlat(ij) = inrlat(ij) + rhoi * &
          &          ( ax(ij, k)*hix(ij, k) - laxhix(ij, k) )
        imrsmi(ij) = imrsmi(ij) - rhos * &
          &          ( ax(ij, k)*hsx(ij, k) - laxhsx(ij, k) )
     end do
  end do

! *** heat and freshwater budget ***
  do ij = 1, nxydim
     ftitd(ij) = 0.d0
  end do
  do ij = ijtstr, ijtend
     wi(ij) = (ax(ij, 0) * hix(ij, 0) - axhix(ij, 0)) / rri / ts
     ws(ij) = (ax(ij, 0) * hsx(ij, 0) - axhsx(ij, 0)) / rrs / ts
     wsn(ij) = (ax(ij, 0) * hsx(ij, 0) - axhsxn(ij, 0)) / rrs / ts
     wen(ij) = (ax(ij, 0) * eix(ij, 0) - axeixn(ij, 0)) / rri / ts
     fs(ij) = 0.d0
  end do
  do k = 1, nic
     do ij = ijtstr, ijtend
        wi(ij) = wi(ij) &
          &    + (ax(ij, k) * hix(ij, k) - axhix(ij, k)) &
          &      / rri / ts
        ws(ij) = ws(ij) &
          &    + (ax(ij, k) * hsx(ij, k) - axhsx(ij, k)) &
          &      / rrs / ts
        wsn(ij) = wsn(ij) &
          &     + (ax(ij, k) * hsx(ij, k) - axhsxn(ij, k)) &
          &       / rrs / ts
        wen(ij) = wen(ij) &
          &     + (ax(ij, k) * eix(ij, k) - axeixn(ij, k)) &
          &       / rri / ts
        fs(ij) = fs(ij) - subi(ij, k) * si * amskt(ij, kstr)
        ftitd(ij) = ftitd(ij) &
          &       - qio(ij, k) * az(ij, k) * amskt(ij, kstr)
        fdd(ij) = fdd(ij) &
          &     + (ax(ij, k) * dsdx(ij, k) - axdsd(ij, k)) &
          &       / ts
        fdb(ij) = fdb(ij) &
          &     + (ax(ij, k) * dsbx(ij, k) - axdsb(ij, k)) &
          &       / ts
     end do
  end do
  do ij = ijtstr, ijtend
     ft(ij, 2) = (  evap(ij) - prec(ij) - roff(ij) &
       &          + ws(ij) + wi(ij) + wiadjs(ij)) * amskt(ij, kstr)
     fs(ij) = fs(ij) + wi(ij) * si * amskt(ij, kstr)
     ft(ij, 1) = - ft(ij, 1) &
       &         + hfus / cpo * (wsn(ij) - snow(ij)) &
       &         + wen(ij) / cpo &
       &         + weadjs(ij) / cpo
     ft(ij, 1) = ft(ij, 1) * amskt(ij, kstr)
     ftitd(ij) = ftitd(ij) + amskt(ij, kstr) * &
       &         ( rhoo * hfus * wi(ij) &
       &         + rhoo * hfus * wsn(ij) )
!       &         + rhoo * wen(ij) )
     fdd(ij) = fdd(ij) * amskt(ij, kstr)
     fdb(ij) = fdb(ij) * amskt(ij, kstr)
  end do

! *** merging newly formed ice into the category 1 ***
  do ij = ijtstr, ijtend
     axhix(ij, 1) = ax(ij, 1) * hix(ij, 1) &
       &          + ax(ij, 0) * hix(ij, 0)
     axeix(ij, 1) = ax(ij, 1) * eix(ij, 1) &
       &          + ax(ij, 0) * eix(ij, 0)
     axhsx(ij, 1) = ax(ij, 1) * hsx(ij, 1)
     axflv(ij, 1) = ax(ij, 1) * frlvx(ij, 1)
     axvmp(ij, 1) = ax(ij, 1) * vmpx(ij, 1)
     axfmp(ij, 1) = ax(ij, 1) * frmpx(ij, 1)
     axdsd(ij, 1) = ax(ij, 1) * dsdx(ij, 1)
     axdsb(ij, 1) = ax(ij, 1) * dsbx(ij, 1)
     hix(ij, 1) = max(hix(ij, 1), hic(1))
     hix(ij, 0) = 0.d0
     eix(ij, 0) = 0.d0
     tix(ij, 0) = tmi
     hsx(ij, 0) = 0.d0
     frlvx(ij, 0) = 1.d0
     vmpx(ij, 0) = 0.d0
     frmpx(ij, 0) = 0.d0
     dsdx(ij, 0) = 0.d0
     dsbx(ij, 0) = 0.d0
  end do
  do ij = ijtstr, ijtend
     danew = ax(ij, 1)
     ax(ij, 1) = axhix(ij, 1) / hix(ij, 1)
!     ax1max = az(ij, 0) + az(ij, 1)
     ax1max = az(ij, 0) + az(ij, 1) + da(ij, 0) + da(ij, 1) 
     if (ax(ij, 1) .gt. ax1max) then
        danew = ax1max - danew
        ax(ij, 1) = ax1max
        hix(ij, 1) = axhix(ij, 1) / ax1max
        eix(ij, 1) = axeix(ij, 1) / ax1max
        tix(ij, 1) = ti(eix(ij, 1) / hix(ij, 1), si)
        hsx(ij, 1) = axhsx(ij, 1) / ax1max
        vmpx(ij, 1) = axvmp(ij, 1) / ax1max
        dsdx(ij, 1) = axdsd(ij, 1) / ax1max
        dsbx(ij, 1) = axdsb(ij, 1) / ax1max
        if (danew > 0.0d0) then  !! new ice formed
           frlvx(ij, 1) = (axflv(ij, 1) + danew * 1.0d0) / ax1max
           frmpx(ij, 1) = axfmp(ij, 1) / ax1max
        end if
!       frlvx and frmpx do not change when danew < 0 (adjustment)
!       asx does not change
     else if (ax(ij, 1) .gt. 0.d0) then
        danew = ax(ij, 1) - danew
        hsx(ij, 1) = axhsx(ij, 1) / ax(ij, 1)
        eix(ij, 1) = axeix(ij, 1) / ax(ij, 1)
        tix(ij, 1) = ti(eix(ij, 1) / hix(ij, 1), si)
        vmpx(ij, 1) = axvmp(ij, 1) / ax(ij, 1)
        dsdx(ij, 1) = axdsd(ij, 1) / ax(ij, 1)
        dsbx(ij, 1) = axdsb(ij, 1) / ax(ij, 1)
        if (danew > 0.0d0) then  !! new ice formed
           frlvx(ij, 1) = (axflv(ij, 1) + danew * 1.0d0) / ax(ij, 1)
           frmpx(ij, 1) = axfmp(ij, 1) / ax(ij, 1)
        end if
!       frlvx and frmpx do not change when danew < 0 (adjustment)
!       asx does not change
     end if
  end do

! for cmip5 output: unit conversion
  do ij = ijtstr, ijtend
     igrfra(ij) = igrfra(ij) / ts
     igrcon(ij) = igrcon(ij) / ts
     igrsni(ij) = igrsni(ij) / ts
     igrsfl(ij) = igrsfl(ij) / ts
     inrlat(ij) = inrlat(ij) / ts
     imrsno(ij) = imrsno(ij) / ts
     imrsmi(ij) = imrsmi(ij) / ts
     imrisf(ij) = imrisf(ij) / ts
     imribs(ij) = imribs(ij) / ts
  end do

!  do k = 1, nic
!     do ij = 1, nxydim
!        impth2(ij) = impth2(ij) + ax(ij, k) * vmpx(ij, k)
!     end do
!  end do

  do k = 1, nic
     do ij = ijtstr, ijtend
        impinc(ij, k) = impinc(ij, k) / ts
        impfrz(ij, k) = impfrz(ij, k) / ts
        improf(ij, k) = improf(ij, k) / ts
     end do
  end do

  call cofpwi( &
    &               ws,     wi)

! debug code
!  do k = 1, nic
!     do ij = ijtstr, ijtend
!        if ((frlvx(ij, k).lt.0.d0).or.(frlvx(ij, k) > (1.d0+1.0d-9))) then
!           write(0,*) '##ipthm; frlvx##', myrank, ij, k, frlvx(ij, k)
!           stop
!        end if        
!        if ((frmpx(ij, k).lt.0.d0).or.(frmpx(ij, k) > (1.d0+1.0d-9))) then
!           write(0,*) '##ipthm; frmpx##', myrank, ij, k, frmpx(ij, k)
!        end if
!     end do
!  end do

  return
end subroutine ptherm

!***********************************************************************
subroutine ipsage( &
  &                  asx, &
  &                  hsx,    tsi,    snow,   dsdx,   dsbx)

  use ufile
  use qckot

  real(8), intent(inout) ::    asx(nxydim, 0:nic)
  real(8), intent(in)    ::    hsx(nxydim, 0:nic),    tsi(nxydim, 0:nic)
  real(8), intent(in)    ::   dsdx(nxydim, 0:nic),   dsbx(nxydim, 0:nic)
  real(8), intent(in)    ::   snow(nxydim)

  real(8), save ::  ildir(nxydim)
  real(8), save :: wabdst, wabblc

  real(8) ::   rafr3(nxydim, 0:nic)
  real(8) ::  dscppm(nxydim, 0:nic), dsmppm(nxydim, 0:nic)
  real(8) ::    dscb(nxydim, 0:nic),   dstm(nxydim, 0:nic)

  logical, save :: ofirst = .true.

  real(8) ::  pi, omega, cor50s, cort, abdst, abblc
  real(8) ::  agesn, agefr1, agefr2, agefr3, agefct, rfrfct, snwfal
  integer :: ij, k
  integer :: ifpar,  jfpar,  istat

  if (oinit .or. ofinal) then
     return
  end if

  if (ofirst) then
     ofirst = .false.
     call rewnml(ifpar, jfpar)
     read (ifpar, nmsage, iostat=istat)
     call cstnml(jfpar, 'ptherm', 'nmsage', istat)
     write(jfpar, nmsage)
     call rewnml(ifpar, jfpar)
     read (ifpar, nmsaab, iostat=istat)
     call cstnml(jfpar, 'ptherm', 'nmsaab', istat)
     write(jfpar, nmsaab)

     do ij = 1, nxydim
        ildir(ij) = 0.0d0
     end do

     pi = atan( 1.d0 )*4.d0
     omega = 2.d0 * pi / 86400.d0
     cor50s = - 2.d0 * omega * sin( pi*50.d0/180.d0 )
     do ij = ijstr, ijend
        cort=(cor(ij)+cor(ij+lw)+cor(ij+lsw)+cor(ij+ls))*0.25d0
        if (cort .le. cor50s) then
           ildir(ij) = 1.0d0
        end if
     end do
     abdst = abduvs*wgtvs + abduni*wgtni + abduir*wgtir
     abblc = abbcvs*wgtvs + abbcni*wgtni + abbcir*wgtir
     wabdst = abdst / (abdst + abblc)                  
     wabblc = abblc / (abdst + abblc) 
  end if

  do k = 1, nic
     do ij = ijstr, ijend
        dscb(ij, k) = wabdst*dsdx(ij, k) + wabblc*dsbx(ij, k)
        dstm(ij, k) = dsdx(ij, k) + dsbx(ij, k)
        if (hsx(ij, k) .gt. 0.d0) then
!!         Multiplier 1.D6 is for converting into ppmw
           dscppm(ij, k) = 1.d6 * dscb(ij, k) / (hsx(ij, k) * rhos)
           dsmppm(ij, k) = 1.d6 * dstm(ij, k) / (hsx(ij, k) * rhos)
        else
           dscppm(ij, k) = 0.d0
           dsmppm(ij, k) = 0.d0
        end if
     end do
  end do

  do k = 1, nic
     do ij = ijstr, ijend
!           Yang et al. (1997) snow aging
        if (hsx(ij, k) .le. 0.d0) then !! no snow
           asx(ij, k) = 0.0d0  !! 0.d0 = fresh snow
           rafr3(ij, k) = 0.d0
        else
           agesn  = asx(ij, k) / ( 1.d0 - asx(ij, k) )
           agefr1 = exp( ftage * ( 1.d0 / tmelt &
             &                   - 1.d0 / (tsi(ij, k)+kelvin) ) )
           agefr2 = min( agefr1**10, 1.0d0 )
           if (oadst) then
              agefr3 = min( adirtc + adirts * dscppm(ij, k), &
                &           adirtm )
           else
              agefr3 = adirt0 * (1.0d0 - ildir(ij)) &
                &    + adirtc * ildir(ij)
           end if
           rafr3(ij, k) = agefr3
           agefct = agefr1 + agefr2 + agefr3
           snwfal = snow(ij) * dt 
           rfrfct = max(0.0d0, 1.0d0 - max(snwfal, 0.0d0) / snrfrs)

           agesn = ( agesn + agefct * dt / tauage) * rfrfct
           asx(ij, k) = max(0.0d0, &
             &              min( alfmax, agesn / (1.0d0 + agesn) ) )
        end if
     end do
  end do

  call chekin( rafr3(1,1), 'RAFR3', &
   &          'r3 in Yang parameterization', '', &
   &              nx,     ny,    nic, nxydim * nic, 'OCICET')
  call chekin( dscb(1,1), 'DSCB', &
   &          'Concentration of dust, combined', 'g/cm^2', &
   &              nx,     ny,    nic, nxydim * nic, 'OCICET')
  call chekin( dstm(1,1), 'DSTM', &
   &          'Concentration of dust, total', 'g/cm^2', &
   &              nx,     ny,    nic, nxydim * nic, 'OCICET')
  call chekin( dscppm(1,1), 'OCDST', &
   &          'Conc. of conbined dust / snow', 'ppmw', &
   &              nx,     ny,    nic, nxydim * nic, 'OCICET')
  call chekin( dsmppm(1,1), 'OCDSTM', &
   &          'Conc. of dust(mass) / snow', 'ppmw', &
   &              nx,     ny,    nic, nxydim * nic, 'OCICET')

  return
end subroutine ipsage

!***********************************************************************
subroutine idfrmp( &
  &                frlvx,   vmpx,  frmpx, &
  &               improf, &
  &                   ax,    hix,    hsx)

  use ufile

  real(8), intent(inout) ::  frlvx(nxydim, 0:nic)  !! ratio of level-ice, 0 - 1. 
  real(8), intent(inout) ::   vmpx(nxydim, 0:nic)  !! volume per unit area of ice [cm]
  real(8), intent(inout) ::  frmpx(nxydim, 0:nic)  !! ratio covered by MP, max. 1 
  real(8), intent(inout) :: improf(nxydim, 0:nic)
  real(8), intent(in)    ::     ax(nxydim, 0:nic)
  real(8), intent(in)    ::    hix(nxydim, 0:nic),    hsx(nxydim, 0:nic)

  real(8), save ::  csfrc
  logical, save :: ofirst = .true.

  real(8) ::  axvmp(nxydim,   nic)
  real(8) :: fbarei, hmp
  logical :: oromp
  integer :: ij, k
  integer :: ifpar,  jfpar,  istat

  if (ofinal) then
     return
  end if

  if (ofirst) then
     ofirst = .false.
     call rewnml(ifpar, jfpar)
     read (ifpar, nmsfrc, iostat=istat)
     call cstnml(jfpar, 'ptherm', 'nmsfrc', istat)
     write(jfpar, nmsfrc)
     call rewnml(ifpar, jfpar)
     read (ifpar, nmsage, iostat=istat)
     call cstnml(jfpar, 'ptherm', 'nmsage', istat)
     write(jfpar, nmsage)
     call rewnml(ifpar, jfpar)
     read (ifpar, nmmpnd, iostat=istat)
     call cstnml(jfpar, 'ptherm', 'nmmpnd', istat)
     write(jfpar, nmmpnd)
      
     if (oasfrc) then
        csfrc = 1.0d0
     else
        csfrc = 0.0d0
     end if
  end if

  if (impnd == 0) then
     do k = 1, nic
        do ij = ijstr, ijend
           frmpx(ij, k) = 0.0d0
           vmpx(ij, k) = 0.0d0
        end do
     end do
     return
  end if

  do k = 1, nic
     do ij = 1, nxydim
        axvmp(ij, k) = ax(ij, k) * vmpx(ij, k)
     end do
  end do

  if (impnd == 1) then  !! Holland MP param.
     do k = 1, nic
        do ij = ijstr, ijend
           oromp = .false.
           if (hix(ij, k) .lt. hminmp) then
              vmpx(ij, k) = 0.0d0
           end if
           fbarei = csfrc * (1.0d0 - hsx(ij, k)/(alspat + hsx(ij, k))) &
             &    + (1.0d0 - csfrc) * &
             &      ( 0.5d0 + sign(0.5d0, alspat - hsx(ij, k)))
           frmpx(ij, k) = sqrt(vmpx(ij,k)/rtdpmp)
           if (frmpx(ij, k) .gt. fbarei) then
              frmpx(ij, k) = fbarei
              oromp = .true.
           end if
           hmp = rtdpmp * frmpx(ij, k)
           if (hmp .gt. rtmxmp*hix(ij, k)) then
              hmp = rtmxmp*hix(ij, k)
              frmpx(ij, k) = hmp/rtdpmp
              oromp = .true.
           end if
           if (oromp) then
              vmpx(ij, k) = frmpx(ij, k) * hmp
           end if
        end do
     end do
  else if (impnd == 2) then   !! Hunke MP param.
     do k = 1, nic
        do ij = ijstr, ijend
           oromp = .false.
           if (frmpx(ij, k) <= 0.0d0) then
              frmpx(ij, k) = 0.0d0
              vmpx(ij, k) = 0.0d0
           else if (frmpx(ij, k) > frlvx(ij, k)) then
              vmpx(ij, k) = vmpx(ij, k) * frlvx(ij, k) / frmpx(ij, k)
              frmpx(ij, k) = frlvx(ij, k)
           end if
        end do
     end do
  end if

! Applying a limiter to vmpx/frmpx
  do k = 1, nic
     do ij = ijstr, ijend
        if ((frmpx(ij, k) < frmpmn).or.(vmpx(ij, k) < vmpmin)) then
           frmpx(ij, k) = 0.0d0
           vmpx(ij, k) = 0.0d0
        end if
     end do
  end do

  if (oinit) then
     return
  end if
  
  do k = 1, nic
     do ij = 1, nxydim
        improf(ij, k) = improf(ij, k) - &
          &           ( ax(ij, k) * vmpx(ij, k) - axvmp(ij, k) ) &
          &           / ts * amskt(ij, kstr)
     end do
  end do

  return
end subroutine idfrmp

end module ipthm
