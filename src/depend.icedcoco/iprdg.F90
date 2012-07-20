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
!
! ---------------------------------------------------------------------

  use zocdim, only: &
    & nxydim,    nic, ijtstr, ijtend, &
    &     lw,     ls,    lsw, &
    &  oinit, ofinal
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
  &                    az,    hiz,    eiz,    hsz, &
  &                    ui,     vi )

  real(8), intent(out)   ::    pice(nxydim)
  real(8), intent(inout) ::      ax(nxydim, 0:nic)
  real(8), intent(inout) ::     hix(nxydim, 0:nic)
  real(8), intent(inout) ::     eix(nxydim, 0:nic)
  real(8), intent(inout) ::     hsx(nxydim, 0:nic)
  real(8), intent(inout) ::     tix(nxydim, 0:nic)
  real(8), intent(in)    ::      az(nxydim, 0:nic)
  real(8), intent(in)    ::     hiz(nxydim, 0:nic)
  real(8), intent(in)    ::     eiz(nxydim, 0:nic)
  real(8), intent(in)    ::     hsz(nxydim, 0:nic)
  real(8), intent(in)    ::      ui(nxydim),     vi(nxydim)

  real(8) ::   axhix(nxydim, 0:nic),  axhsx(nxydim, 0:nic)
  real(8) ::   axeix(nxydim, 0:nic)
  real(8) ::    divv(nxydim),  edis(nxydim)
  real(8) ::      wa(nxydim, 0:nic),     wn(nxydim, 0:nic)
  real(8) ::      ww(nxydim)
  real(8) ::       g(nxydim, -1:nic)
  real(8) ::       y(nxydim, -1:nic)
  real(8) ::      da(nxydim, nic)
  real(8) ::    dahi(nxydim, nic),   dahs(nxydim, nic)
  real(8) ::    daei(nxydim, nic)
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
  real(8) ::    hikl,   hskl,   hekl
  integer ::      ij,      k,      l
  integer ::    ijlw,   ijls,  ijlsw
  integer ::   ifpar,  jfpar,  istat

  real(8), save ::  ecc = 2.0d0,  dmin = 2.0d-7,  floss = 17.0d0
  real(8), save ::  cs = 0.5d0,  gridge = 0.15d0,  hridge = 1.0d4
  real(8), save ::  si = 5.0d0

  namelist /nmidyn/ ecc, dmin, floss
  namelist /nmirdg/ cs, gridge, hridge
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
        g(ij, k) = 0.d0
        wn(ij, k) = 0.d0
     end do
  end do
  do ij = 1, nxydim
     pice(ij) = 0.d0
     g(ij, -1) = 0.d0
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
     end do
  end do

  do l = 1, nic
     do k = 1, l
        do ij = ijtstr, ijtend
           hrmax = 2.d0 * sqrt(hridge * hix(ij, k))
           hrmin = 2.d0 * hix(ij, k)
           if (     (hic(l+1) .lt. hrmin) &
             & .or. (hic(l) .gt. hrmax) &
             & .or. (hrmin .ge. hrmax) ) then
              gam(ij, k, l) = 0.d0
           else
              gam(ij, k, l) = (  min(hrmax, hic(l+1)) &
                &              - max(hrmin, hic(l))) * 0.5d0 &
                &             / (hridge - hix(ij, k))
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
        da(ij, k) = edis(ij) * (wn(ij, k) - wa(ij, k))
     end do
  end do
  do l = 1, nic
     do ij = ijtstr, ijtend
        dahi(ij, l) = - hix(ij, l) * wa(ij, l) * edis(ij)
        dahs(ij, l) = - hsx(ij, l) * wa(ij, l) * edis(ij)
        daei(ij, l) = - eix(ij, l) * wa(ij, l) * edis(ij)
        pice(ij) = pice(ij) &
          &      - hix(ij, l) * hix(ij, l) * wa(ij, l) * pifct
        do k = 1, l
           hrmax = 2.d0 * sqrt(hridge * hix(ij, k))
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
           dahi(ij, l) = dahi(ij, l) &
             &         + hikl * wa(ij, k) * gam(ij, k, l) * &
             &           edis(ij)
           dahs(ij, l) = dahs(ij, l) &
             &         + hskl * wa(ij, k) * gam(ij, k, l) * &
             &           edis(ij)
           daei(ij, l) = daei(ij, l) &
             &         + hekl * wa(ij ,k) * gam(ij, k, l) * &
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
     if (ax(ij, 0) .eq. 1.d0) then
        do k = 1, nic
           ax(ij, k) = 0.d0
           hix(ij, k) = hic(k)
           hsx(ij, k) = 0.d0
           eix(ij, k) = 0.d0
           hix(ij, 0) = hix(ij, 0) + axhix(ij, k)
           hsx(ij, 0) = hsx(ij, 0) + axhsx(ij, k)
           eix(ij, 0) = eix(ij, 0) + axeix(ij, k)
           dahi(ij, k) = 0.d0
           dahs(ij, k) = 0.d0
           daei(ij, k) = 0.d0
        end do
     end if
  end do

  do k = 1, nic
     do ij = ijtstr, ijtend
        axhix(ij, k) = axhix(ij, k) + ts * dahi(ij, k)
        axhsx(ij, k) = axhsx(ij, k) + ts * dahs(ij, k)
        axeix(ij, k) = axeix(ij, k) + ts * daei(ij, k)
     end do
  end do

! *** ridging process conserves the mass of ice and snow, and the
! *** amount of the thickest category ice is never reduced by ridging.
! *** thus, the situation that axhix(ij, nic) becomes less than zero
! *** needs not be considered in the adjustment below. 
  do k = 1, nic-1
     do ij = ijtstr, ijtend
        if (ax(ij, k) .eq. 0.d0) then
           axhix(ij, k+1) = axhix(ij, k+1) + axhix(ij, k)
           axhsx(ij, k+1) = axhsx(ij, k+1) + axhsx(ij, k)
           axeix(ij, k+1) = axeix(ij, k+1) + axeix(ij, k)
           axhix(ij, k) = 0.d0
           axhsx(ij, k) = 0.d0
           axeix(ij, k) = 0.d0
           hix(ij, k) = hic(k)
           hsx(ij, k) = 0.d0
           eix(ij, k) = 0.d0
        else if (axhix(ij, k) .lt. 0.d0) then
           axhix(ij, k+1) = axhix(ij, k+1) + axhix(ij, k)
           axhsx(ij, k+1) = axhsx(ij, k+1) + axhsx(ij, k)
           axeix(ij, k+1) = axeix(ij, k+1) + axeix(ij, k)
           ax(ij, k+1) = ax(ij, k+1) + ax(ij, k)
           ax(ij, k) = 0.d0
           axhix(ij, k) = 0.d0
           axhsx(ij, k) = 0.d0
           axeix(ij, k) = 0.d0
           hix(ij, k) = hic(k)
           hsx(ij, k) = 0.d0
           eix(ij, k) = 0.d0
        else
           if (axhsx(ij, k) .lt. 0.d0) then
              axhsx(ij, k+1) = axhsx(ij, k+1) + axhsx(ij, k)
              axhsx(ij, k) = 0.d0
              hsx(ij, k) = 0.d0
           end if
           if (axeix(ij, k) .lt. 0.d0) then
              axeix(ij, k+1) = axeix(ij, k+1) + axeix(ij, k)
              axeix(ij, k) = 0.d0
              eix(ij, k) = 0.d0
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
        else
           tix(ij, k) = tmi
        end if
     end do
  end do

  return

end subroutine pridge

end module iprdg
