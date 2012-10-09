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
!
! ---------------------------------------------------------------------

  use zocdim, only: &
    & nxydim,  nzdim,   kstr,  ntdim,    nic, ijtstr, ijtend, &
    &  oinit, ofinal
  use zocgrd, only: &
    &    hic,     ts
  use zocmsk, only: &
    &  amskt
  use zocphy, only: &
    &   rhoo,   rhoi,   rhos,   hfus,    cpo,    cpi,   dtds

  implicit none

  private

  public :: ptherm

contains

subroutine ptherm( &
  &                    ax,    hix,    hsx, &
  &                   eix,    tix, &
  &                  prec,   snow,     ft,     fs, &
  &                 ftitd, igrfra, igrcon, igrsni, &
  &                inrlat, imrsno, imrisf, imribs, &
  &                    tx, &
  &                   wio,    wao,    was,    wil, &
  &                  evap,   subi,   roff, adjlat, &
  &                   qio )

  use ufile
  use qckot

#ifdef OPT_PARALLEL
#include "mpif.h"
#endif

  real(8), intent(inout) ::     ax(nxydim, 0:nic),    hix(nxydim, 0:nic)
  real(8), intent(inout) ::    hsx(nxydim, 0:nic)
  real(8), intent(inout) ::    eix(nxydim, 0:nic),    tix(nxydim, 0:nic)
  real(8), intent(inout) ::   prec(nxydim),   snow(nxydim)
  real(8), intent(inout) ::     ft(nxydim, ntdim),     fs(nxydim)
  real(8), intent(out)   ::  ftitd(nxydim)
  real(8), intent(out)   :: igrfra(nxydim), igrcon(nxydim), igrsni(nxydim)
  real(8), intent(out)   :: inrlat(nxydim)
  real(8), intent(out)   :: imrsno(nxydim), imrisf(nxydim), imribs(nxydim)
  real(8), intent(inout) ::    wio(nxydim, nic)
  real(8), intent(in)    ::    was(nxydim, nic)
  real(8), intent(in)    ::    wil(nxydim, nic)
  real(8), intent(in)    ::    wao(nxydim)
  real(8), intent(in)    ::   evap(nxydim),   subi(nxydim, nic)
  real(8), intent(in)    ::   roff(nxydim), adjlat(nxydim)
  real(8), intent(in)    ::     tx(nxydim, nzdim, ntdim)
  real(8), intent(in)    ::    qio(nxydim, nic)

  real(8) ::     az(nxydim, 0:nic)
  real(8) ::  axhix(nxydim, 0:nic)
  real(8) ::  axhsx(nxydim, 0:nic), axhsxn(nxydim, 0:nic)
  real(8) ::  axeix(nxydim, 0:nic), axeixn(nxydim, 0:nic)
  real(8) ::    wai(nxydim, nic)
  real(8) ::     wi(nxydim),     ws(nxydim)
  real(8) ::    wen(nxydim),    wsn(nxydim)
  real(8) ::    hiz(nxydim, 0:nic)
  real(8) ::  aflrm(nxydim, 0:nic), aflrmc(nxydim)
  real(8) ::   fdtn(nxydim, 0:nic),  fdtcn(nxydim, 0:nic)
  real(8) ::   hicn(nxydim, 0:nic)
  real(8) ::     g0(nxydim, 0:nic),     g1(nxydim, 0:nic)
  real(8) ::    hil(nxydim, 0:nic),    hir(nxydim, 0:nic)
  real(8) ::     da(nxydim, 0:nic),   dahi(nxydim, 0:nic)
  real(8) ::   dahs(nxydim, 0:nic),   daei(nxydim, 0:nic)
  real(8) :: laxhix(nxydim, 0:nic), laxhsx(nxydim, 0:nic) 
  real(8) :: laxeix(nxydim, 0:nic)
  real(8) :: daxeit(nxydim, 0:nic), daxeib(nxydim, 0:nic)
!      COMMON /WORK/ AZ, AXHIX, AXHSX, AXHSXN, WAI,
!     &              WI, WS, WEN, WSN, AXEIX, AXEIXN,
!     &              HIZ, AFLRM, AFLRMC, FDTN, FDTCN, HICN,
!     &              G0, G1, HIL, HIR, DA, DAHI,
!     &              DAHS, DAEI, LAXHIX, LAXHSX, LAXEIX,
!     &              DAXEIT, DAXEIB

  real(8), save ::    rri,    rrs, rorirs
  real(8), save ::  rsfus
  real(8), save ::    tmi
  real(8), save :: epsaei
  real(8), save ::   hic0(0:nic+1)
  logical, save :: ofirst = .true.

  real(8) ::   wres
  real(8) ::   hsxo,    dhs
  real(8) :: ax1max
  real(8) ::   etan,   etar, etanrr,   gil,    gir
  real(8) ::     x0,     x1,   gint
  real(8) ::   fahi,   fahs,   faei,   pvol
  real(8) :: daxhix, daxeix, daxhit, daxhib,   eieq
  integer ::     ij,      k
  integer ::  ifpar,  jfpar,  istat

  real(8), save ::    eps = 1.0d-3,   epsl = 1.0d-6

  real(8), save ::  amin = 1.0d-6,  amax = 1.0d0,  si = 5.0d0
  integer, save ::  mic = nic

  namelist /nmamin/ amin, amax, mic
  namelist /nmislt/ si

!===== define statement function 
#include "zocite.F90"
!===== 

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

     tmi = dtds * si
     rri    = rhoo / rhoi
     rrs    = rhoo / rhos
     rorirs = (rhoo - rhoi) / rhos
     rsfus  = rhos * hfus
     epsaei = amin * abs(ei(tmi+eps, si))
     hic0(0) = 0.0d0
     do k = 1, nic+1
        hic0(k) = hic(k)
     end do
  end if

  do ij = 1, nxydim
     igrfra(ij) = 0.0d0
     igrcon(ij) = 0.0d0
     igrsni(ij) = 0.0d0
     inrlat(ij) = 0.0d0
     imrsno(ij) = 0.0d0
     imrisf(ij) = 0.0d0
     imribs(ij) = 0.0d0
  end do

  do k = 0, nic
     do ij = 1, nxydim
        az(ij, k) = ax(ij, k)
        hiz(ij, k) = hix(ij, k)
        aflrm(ij, k) = 0.0d0
        daxeit(ij, k) = 0.0d0
        daxeib(ij, k) = 0.0d0
     end do
  end do

  do k = 0, nic
     do ij = 1, nxydim
        axhix(ij, k) = ax(ij, k) * hix(ij, k)
        axhsx(ij, k) = ax(ij, k) * hsx(ij, k)
        axeix(ij, k) = ax(ij, k) * eix(ij, k)
     end do
  end do

! *** snowfall ***
  do k = 1, nic
     do ij = ijtstr, ijtend
        if (ax(ij, k) .gt. 0.d0) then
           hsx(ij, k) = axhsx(ij, k) / ax(ij, k) &
             &        + ts * rrs * snow(ij)
           axhsx(ij, k) = ax(ij, k) * hsx(ij, k)
        else
           hsx(ij, k) = 0.d0
           axhsx(ij, k) = 0.d0
        end if
     end do
  end do
  do ij = ijtstr, ijtend
     snow(ij) = ax(ij, 0) * snow(ij)
     prec(ij) = prec(ij) + snow(ij)
  end do

! *** snow melting ***
  do k = 1, nic
     do ij = ijtstr, ijtend
        wres = axhsx(ij, k) * rsfus / ts + was(ij, k)
        if (ax(ij, k) .gt. 0.d0) then
           if (wres .lt. 0.d0) then
              hsx(ij, k) = 0.d0
              wai(ij, k) = wres
           else
              hsx(ij, k) = wres * ts / rsfus / ax(ij, k)
              wai(ij, k) = 0.d0
           end if
        else
           wai(ij, k) = was(ij, k)
        end if
        axhsxn(ij, k) = ax(ij, k) * hsx(ij, k) * amskt(ij, kstr)
        imrsno(ij) = imrsno(ij) - &
          &          ( axhsxn(ij, k) - axhsx(ij, k) )
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
        else
           wio(ij, k) = wio(ij, k) + wai(ij, k)
           axeixn(ij, k) = axeix(ij, k)
        end if
        daxeit(ij, k) = axeixn(ij, k) - axeix(ij, k)
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
        daxeib(ij, k) = axeixn(ij, k)
        axeixn(ij, k) = axeixn(ij, k) &
          &           + wio(ij, k) * ts / rhoi &
          &           * amskt(ij, kstr)
     end do
  end do
  do k = 1, nic
     do ij = ijtstr, ijtend
        if (axeixn(ij, k) .le. 0.d0) then
           ax(ij, k) = 0.d0
           hix(ij, k) = hic(k)
           hsx(ij, k) = 0.d0
           tix(ij, k) = tmi
           eix(ij, k) = 0.d0
        else
           eix(ij, k) = axeixn(ij, k) / ax(ij, k)
           hix(ij, k) = eix(ij, k) / ei(tix(ij, k), si)
        end if
        daxeib(ij, k) = ax(ij, k) * eix(ij, k) - daxeib(ij, k)
     end do
  end do

  do k = 1, nic
     do ij = ijtstr, ijtend
        if (az(ij, k) .gt. 0.0d0) then
           daxhix = ax(ij, k) * hix(ij, k) - axhix(ij, k)
           daxeix = daxeit(ij, k) + daxeib(ij, k)
           if ( abs(daxeix) .lt. epsaei ) then
              daxhit = daxeit(ij, k) / ei(tix(ij, k), si)
              daxhib = daxeib(ij, k) / ei(tix(ij, k), si)
           else
              daxhit = daxhix * daxeit(ij, k) / daxeix
              daxhib = daxhix * daxeib(ij, k) / daxeix
           endif
           imrisf(ij) = imrisf(ij) - rhoi * daxhit
           igrcon(ij) = igrcon(ij) + rhoi * max(0.0d0, daxhib)
           imribs(ij) = imribs(ij) - rhoi * min(0.0d0, daxhib)
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
        igrsni(ij) = igrsni(ij) + ax(ij, k) * dhs * rhos
        if (ax(ij, k) .gt. 0.d0) then
           tix(ij, k) = ti(eix(ij, k)/hix(ij, k), si)
        else
           tix(ij, k) = tmi
        end if
     end do
  end do

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
        da(ij, k) = 0.0d0
        dahi(ij, k) = 0.0d0
        dahs(ij, k) = 0.0d0
        daei(ij, k) = 0.0d0
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
        gint = min( gint, &
          &         ax(ij, 1)*(1.0d0 - hix(ij, 1)/hiz(ij, 1)) )
!       only area flux can across the lowest boundary
        da(ij, 1) = da(ij, 1) - gint
        da(ij, 0) = da(ij, 0) + gint
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
           else
              gil = hic(k) - hil(ij, k-1)
              gir = max(hir(ij, k-1) - hil(ij, k-1), gil)
              x0 = gir - gil
              x1 = 0.5d0 * (gir*gir - gil*gil)
              gint = g0(ij, k-1)*x0 + g1(ij, k-1)*x1
              x0 = 0.5d0 * (gir*gir - gil*gil)
              x1 = (gir*gir*gir - gil*gil*gil) / 3.0d0
              fahi = hil(ij, k-1) * gint + &
                &    g0(ij, k-1)*x0 + g1(ij, k-1)*x1
              pvol = fahi / laxhix(ij, k-1)
              fahs = laxhsx(ij, k-1) * pvol
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
           da(ij, k) = da(ij, k) + gint
           dahi(ij, k) = dahi(ij, k) + fahi
           dahs(ij, k) = dahs(ij, k) + fahs
           daei(ij, k) = daei(ij, k) + faei
        elseif ( ( hicn(ij, k) .lt. hic(k) ) .and. &
          &      ( aflrm(ij, k) .eq. 1.0d0 ) ) then
           if (hir(ij, k) .le. hic(k)) then
!             all sea ice transfered to the lower category
              gint = ax(ij, k)
              fahi = laxhix(ij, k)
              pvol = 1.0d0
              fahs = laxhsx(ij, k)
              faei = laxeix(ij, k)
           else
!              gil = hil(ij, k) - hil(ij, k)
              gil = 0.0d0
              gir = max(hic(k) - hil(ij, k), gil)
              x0 = gir - gil
              x1 = 0.5d0 * (gir*gir - gil*gil)
              gint = g0(ij, k)*x0 + g1(ij, k)*x1
              gint = min(gint, ax(ij, k))
              x0 = 0.5d0 * (gir*gir - gil*gil)
              x1 = (gir*gir*gir - gil*gil*gil) / 3.0d0
              fahi = hil(ij, k) * gint + &
                &    g0(ij, k)*x0 + g1(ij, k)*x1
              fahi = min(fahi, laxhix(ij, k))
              pvol = fahi / laxhix(ij, k)
              fahs = laxhsx(ij, k) * pvol
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
           da(ij, k) = da(ij, k) - gint
           dahi(ij, k) = dahi(ij, k) - fahi
           dahs(ij, k) = dahs(ij, k) - fahs
           daei(ij, k) = daei(ij, k) - faei
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
        elseif (ax(ij, k) .le. 0.d0) then
           ax(ij, k) = 0.d0
           hix(ij, k) = hic(k)
           hsx(ij, k) = 0.d0
           tix(ij, k) = tmi
           eix(ij, k) = 0.d0
        else
           hix(ij, k) = laxhix(ij, k) / ax(ij, k)
           hsx(ij, k) = laxhsx(ij, k) / ax(ij, k)
           eix(ij, k) = laxeix(ij, k) / ax(ij, k)
!           tix(ij, k) = ti(eix(ij, k)/hix(ij, k), si)
        endif
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
     end do
  end do
  do k = 1, nic
     do ij = ijtstr, ijtend
        if (ax(ij, k) .gt. 0.d0) then
           ax(ij, k) = ax(ij, k) &
             &       + ts * wil(ij, k) / eix(ij, k) / rhoi &
             &       * amskt(ij, kstr)
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
        else if (ax(ij, k) .le. 0.d0) then
           ax(ij, k) = 0.d0
           hix(ij, k) = hic(k)
           hsx(ij, k) = 0.d0
           tix(ij, k) = tmi
           eix(ij, k) = 0.d0
        else
           eix(ij, k) = laxeix(ij, k) / ax(ij, k)
           tix(ij, k) = ti(eix(ij, k)/hix(ij, k), si)
        end if
        inrlat(ij) = inrlat(ij) + rhoi *  &
          &          ( ax(ij, k)*hix(ij, k) - laxhix(ij, k) )
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
     end do
  end do
  do ij = ijtstr, ijtend
     ft(ij, 2) = (  evap(ij) - prec(ij) - roff(ij) &
       &          + ws(ij) + wi(ij)) * amskt(ij, kstr)
     fs(ij) = fs(ij) + wi(ij) * si * amskt(ij, kstr)
     ft(ij, 1) = - ft(ij, 1) &
       &         + hfus / cpo * (wsn(ij) - snow(ij)) &
       &         + wen(ij) / cpo &
       &         + adjlat(ij) * hfus / cpo 
     ft(ij, 1) = ft(ij, 1) * amskt(ij, kstr)
     ftitd(ij) = ftitd(ij) + amskt(ij, kstr) * &
       &         ( rhoo * hfus * wi(ij) &
       &         + rhoo * hfus * wsn(ij) )
!       &         + rhoo * wen(ij) )
  end do

! *** merging newly formed ice into the category 1 ***
  do ij = ijtstr, ijtend
     axhix(ij, 1) = ax(ij, 1) * hix(ij, 1) &
       &          + ax(ij, 0) * hix(ij, 0)
     axeix(ij, 1) = ax(ij, 1) * eix(ij, 1) &
       &          + ax(ij, 0) * eix(ij, 0)
     axhsx(ij, 1) = ax(ij, 1) * hsx(ij, 1)
     hix(ij, 1) = max(hix(ij, 1), hic(1))
     hix(ij, 0) = 0.d0
     eix(ij, 0) = 0.d0
     tix(ij, 0) = tmi
     hsx(ij, 0) = 0.d0
  end do
  do ij = ijtstr, ijtend
     ax(ij, 1) = axhix(ij, 1) / hix(ij, 1)
!     ax1max = az(ij, 0) + az(ij, 1)
     ax1max = az(ij, 0) + az(ij, 1) + da(ij, 0) + da(ij, 1) 
     if (ax(ij, 1) .gt. ax1max) then
        ax(ij, 1) = ax1max
        hix(ij, 1) = axhix(ij, 1) / ax1max
        eix(ij, 1) = axeix(ij, 1) / ax1max
        tix(ij, 1) = ti(eix(ij, 1) / hix(ij, 1), si)
        hsx(ij, 1) = axhsx(ij, 1) / ax1max
     else if (ax(ij, 1) .gt. 0.d0) then
        hsx(ij, 1) = axhsx(ij, 1) / ax(ij, 1)
        eix(ij, 1) = axeix(ij, 1) / ax(ij, 1)
        tix(ij, 1) = ti(eix(ij, 1) / hix(ij, 1), si)
     end if
  end do

! for cmip5 output: unit conversion
  do ij = ijtstr, ijtend
     igrfra(ij) = igrfra(ij) / ts
     igrcon(ij) = igrcon(ij) / ts
     igrsni(ij) = igrsni(ij) / ts
     inrlat(ij) = inrlat(ij) / ts
     imrsno(ij) = imrsno(ij) / ts
     imrisf(ij) = imrisf(ij) / ts
     imribs(ij) = imribs(ij) / ts
  end do
  
  return
end subroutine ptherm

end module ipthm
