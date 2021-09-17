module ipadv

! ---- information ----------------------------------------------------
!
!  Horizontal advection of ice category fraction, ice volume, and
! snow volume. The current choice of the advection algorithm is the
! simplest first-order upstream. No explict diffusion is taken into
! account. Advection algorithm must be positivity preserving.
!
!  HISTORY
!     '03.08.06  H.Hasumi
!     '07.09.25  H.Hasumi: for COCO4
!     '07.10.03  H.Hasumi: 1-layer sea ice thermodynamics
!     '08.??.??  Y.Komuro
!     '08.06.11  H.Hasumi: initial/final processing
!     '08.07.10  H.Hasumi: initial/final processing
!     '09.05.25  Y.Komuro: CMIP5 output code included
!     '09.09.04  Y.Komuro: extra output code (FEX/FEY)
!     '12.07.20  Y.Komuro: for COCO5.0
!     '21.05.26  Y.Komuro: snow aging & meltpond parametrization 
!
! ---------------------------------------------------------------------

  use zocdim, only: &
    &  nxdim, nxydim,   kstr,    nic, ijtstr, ijtend, &
    &     lw,     le,     ln,     ls,    lsw, &
    &  oinit, ofinal
  use zocgrd, only: &
    &     rx,     ry,    hic,     ts, &
    &    hxu,    hyu,    rxt,    ryt
  use zocmsk, only: &
    &  amskt
  use zocphy, only: &
    &    cpi,   hfus,   dtds

  implicit none

  private

  public :: padvct

contains

subroutine padvct( &
  &                    ax,    hix,    eix,    hsx,    tix, &
  &                   asx,  frlvx,   vmpx,  frmpx,   dsdx,   dsbx, &
  &                    az,    hiz,    eiz,    hsz, &
  &                   fix,    fiy,    fsx,    fsy, &
  &                   fex,    fey, &
  &                   uiy,    viy)
  use ufile
  use zocite

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
  real(8), intent(out)   ::      az(nxydim, 0:nic)
  real(8), intent(out)   ::     hiz(nxydim, 0:nic)
  real(8), intent(out)   ::     eiz(nxydim, 0:nic)
  real(8), intent(out)   ::     hsz(nxydim, 0:nic)
  real(8), intent(out)   ::     fix(nxydim, 0:nic),    fiy(nxydim, 0:nic)
  real(8), intent(out)   ::     fsx(nxydim, 0:nic),    fsy(nxydim, 0:nic)
  real(8), intent(out)   ::     fex(nxydim, 0:nic),    fey(nxydim, 0:nic)
  real(8), intent(in)    ::     uiy(nxydim),    viy(nxydim)

  real(8) ::     fax(nxydim,   nic),    fay(nxydim,   nic)
  real(8) ::    fasx(nxydim,   nic),   fasy(nxydim,   nic)
  real(8) ::    fflx(nxydim,   nic),   ffly(nxydim,   nic)
  real(8) ::    fvmx(nxydim,   nic),   fvmy(nxydim,   nic)
  real(8) ::    ffmx(nxydim,   nic),   ffmy(nxydim,   nic)
  real(8) ::    fddx(nxydim,   nic),   fddy(nxydim,   nic)
  real(8) ::    fdbx(nxydim,   nic),   fdby(nxydim,   nic)
  real(8) ::   axhix(nxydim, 0:nic),  axhsx(nxydim, 0:nic)
  real(8) ::   axeix(nxydim, 0:nic)
  real(8) ::   axasx(nxydim, 0:nic),  axvmp(nxydim, 0:nic)
  real(8) ::   axflv(nxydim, 0:nic),  axfmp(nxydim, 0:nic)
  real(8) ::   axdsd(nxydim, 0:nic),  axdsb(nxydim, 0:nic)
!  common /work/ fax, fay, axhix, axhsx, axeix

  real(8), save ::    tmi
  logical, save :: ofirst = .true.

  real(8) ::       u,     up,     um
  real(8) ::       v,     vp,     vm
  integer ::      ij,      k
  integer ::    ijls,   ijlw,   ijln,   ijle
  integer ::   ifpar,  jfpar,  istat

  real(8), save ::  si = 5.0d0

  namelist /nmislt/ si


  if (oinit .or. ofinal) then
     return
  end if

  if (ofirst) then
     ofirst = .false.
     call rewnml(ifpar, jfpar)
     read (ifpar, nmislt, iostat=istat)
     call cstnml(jfpar, 'padvct', 'nmislt', istat)
     write(jfpar, nmislt)
     tmi = dtds * si
  end if

  do k = 0, nic
     do ij = 1, nxydim
        az (ij, k) = ax (ij, k)
        hiz(ij, k) = hix(ij, k)
        eiz(ij, k) = eix(ij, k)
        hsz(ij, k) = hsx(ij, k)
        axhix(ij, k) = ax(ij, k) * hix(ij, k)
        axeix(ij, k) = ax(ij, k) * eix(ij, k)
        axhsx(ij, k) = ax(ij, k) * hsx(ij, k)
        axasx(ij, k) = ax(ij, k) * asx(ij, k)
        axflv(ij, k) = ax(ij, k) * frlvx(ij, k)
        axvmp(ij, k) = ax(ij, k) * vmpx(ij, k)
        axfmp(ij, k) = ax(ij, k) * frmpx(ij, k)
        axdsd(ij, k) = ax(ij, k) * dsdx(ij, k)
        axdsb(ij, k) = ax(ij, k) * dsbx(ij, k)
        fix(ij, k) = 0.d0
        fiy(ij, k) = 0.d0
        fsx(ij, k) = 0.d0
        fsy(ij, k) = 0.d0
        fex(ij, k) = 0.d0
        fey(ij, k) = 0.d0
     end do
  end do

  do k = 1, nic
     do ij = 1, nxydim
        fax(ij, k) = 0.d0
        fay(ij, k) = 0.d0
        fasx(ij, k) = 0.d0
        fasy(ij, k) = 0.d0
        fflx(ij, k) = 0.d0
        ffly(ij, k) = 0.d0
        fvmx(ij, k) = 0.d0
        fvmy(ij, k) = 0.d0
        ffmx(ij, k) = 0.d0
        ffmy(ij, k) = 0.d0
        fddx(ij, k) = 0.d0
        fddy(ij, k) = 0.d0
        fdbx(ij, k) = 0.d0
        fdby(ij, k) = 0.d0
     end do
  end do

  do k = 1, nic

     do ij = ijtstr, ijtend+1
        ijlw = ij + lw
        u  = (  uiy(ijlw  ) * hyu(ijlw  ) &
          &   + uiy(ij+lsw) * hyu(ij+lsw)) * 0.25d0
        up = u + abs(u)
        um = u - abs(u)
        fax(ij, k) = - (  up * ax(ijlw, k) &
          &          + um * ax(ij, k)) * &
          &         amskt(ij, kstr) * amskt(ijlw, kstr)
        fix(ij, k) = - (  up * axhix(ijlw, k) &
          &          + um * axhix(ij, k)) * &
          &         amskt(ij, kstr) * amskt(ijlw, kstr)
        fex(ij, k) = - (  up * axeix(ijlw, k) &
          &          + um * axeix(ij, k)) * &
          &         amskt(ij, kstr) * amskt(ijlw, kstr)
        fsx(ij, k) = - (  up * axhsx(ijlw, k) &
          &          + um * axhsx(ij, k)) * &
          &         amskt(ij, kstr) * amskt(ijlw, kstr)
        fasx(ij, k) = - (  up * axasx(ijlw, k) &
          &          + um * axasx(ij, k)) * &
          &         amskt(ij, kstr) * amskt(ijlw, kstr)
        fflx(ij, k) = - (  up * axflv(ijlw, k) &
          &          + um * axflv(ij, k)) * &
          &         amskt(ij, kstr) * amskt(ijlw, kstr)
        fvmx(ij, k) = - (  up * axvmp(ijlw, k) &
          &          + um * axvmp(ij, k)) * &
          &         amskt(ij, kstr) * amskt(ijlw, kstr)
        ffmx(ij, k) = - (  up * axfmp(ijlw, k) &
          &          + um * axfmp(ij, k)) * &
          &         amskt(ij, kstr) * amskt(ijlw, kstr)
        fddx(ij, k) = - (  up * axdsd(ijlw, k) &
          &          + um * axdsd(ij, k)) * &
          &         amskt(ij, kstr) * amskt(ijlw, kstr)
        fdbx(ij, k) = - (  up * axdsb(ijlw, k) &
          &          + um * axdsb(ij, k)) * &
          &         amskt(ij, kstr) * amskt(ijlw, kstr)
     end do

     do ij = ijtstr, ijtend+nxdim
        ijls = ij + ls
        v  = (  viy(ijls  ) * hxu(ijls  ) &
          &   + viy(ij+lsw) * hxu(ij+lsw)) * 0.25d0
        vp = v + abs(v)
        vm = v - abs(v)
        fay(ij, k) = - (  vp * ax(ijls, k) &
          &          + vm * ax(ij, k)) * &
          &         amskt(ij, kstr) * amskt(ijls, kstr)
        fiy(ij, k) = - (  vp * axhix(ijls, k) &
          &          + vm * axhix(ij, k)) * &
          &         amskt(ij, kstr) * amskt(ijls, kstr)
        fey(ij, k) = - (  vp * axeix(ijls, k) &
          &          + vm * axeix(ij, k)) * &
          &         amskt(ij, kstr) * amskt(ijls, kstr)
        fsy(ij, k) = - (  vp * axhsx(ijls, k) &
          &          + vm * axhsx(ij, k)) * &
          &         amskt(ij, kstr) * amskt(ijls, kstr)
        fasy(ij, k) = - (  vp * axasx(ijls, k) &
          &          + vm * axasx(ij, k)) * &
          &         amskt(ij, kstr) * amskt(ijls, kstr)
        ffly(ij, k) = - (  vp * axflv(ijls, k) &
          &          + vm * axflv(ij, k)) * &
          &         amskt(ij, kstr) * amskt(ijls, kstr)
        fvmy(ij, k) = - (  vp * axvmp(ijls, k) &
          &          + vm * axvmp(ij, k)) * &
          &         amskt(ij, kstr) * amskt(ijls, kstr)
        ffmy(ij, k) = - (  vp * axfmp(ijls, k) &
          &          + vm * axfmp(ij, k)) * &
          &         amskt(ij, kstr) * amskt(ijls, kstr)
        fddy(ij, k) = - (  vp * axdsd(ijls, k) &
          &          + vm * axdsd(ij, k)) * &
          &         amskt(ij, kstr) * amskt(ijls, kstr)
        fdby(ij, k) = - (  vp * axdsb(ijls, k) &
          &          + vm * axdsb(ij, k)) * &
          &         amskt(ij, kstr) * amskt(ijls, kstr)
     end do

     do ij = ijtstr, ijtend
        ijle = ij + le
        ijln = ij + ln
        ax(ij, k) = ax(ij, k) &
          &       + ts * (  (fax(ijle, k) - fax(ij, k)) * rx &
          &               + (fay(ijln, k) - fay(ij, k)) * ry(ij)) * &
          &         rxt(ij) * ryt(ij) * amskt(ij, kstr)
        axhix(ij, k) = axhix(ij, k) &
          &          + ts * (  (fix(ijle, k) - fix(ij, k)) * rx &
          &                  + (fiy(ijln, k) - fiy(ij, k)) * ry(ij)) * &
          &            rxt(ij) * ryt(ij) * amskt(ij, kstr)
        axeix(ij, k) = axeix(ij, k) &
          &          + ts * (  (fex(ijle, k) - fex(ij, k)) * rx &
          &                  + (fey(ijln, k) - fey(ij, k)) * ry(ij)) * &
          &            rxt(ij) * ryt(ij) * amskt(ij, kstr)
        axhsx(ij, k) = axhsx(ij, k) &
          &          + ts * (  (fsx(ijle, k) - fsx(ij, k)) * rx &
          &                  + (fsy(ijln, k) - fsy(ij, k)) * ry(ij)) * &
          &            rxt(ij) * ryt(ij) * amskt(ij, kstr)
        axasx(ij, k) = axasx(ij, k) &
          &          + ts * (  (fasx(ijle, k) - fasx(ij, k)) * rx &
          &              + (fasy(ijln, k) - fasy(ij, k)) * ry(ij)) * &
          &            rxt(ij) * ryt(ij) * amskt(ij, kstr)
        axflv(ij, k) = axflv(ij, k) &
          &          + ts * (  (fflx(ijle, k) - fflx(ij, k)) * rx &
          &              + (ffly(ijln, k) - ffly(ij, k)) * ry(ij)) * &
          &            rxt(ij) * ryt(ij) * amskt(ij, kstr)
        axvmp(ij, k) = axvmp(ij, k) &
          &          + ts * (  (fvmx(ijle, k) - fvmx(ij, k)) * rx &
          &              + (fvmy(ijln, k) - fvmy(ij, k)) * ry(ij)) * &
          &            rxt(ij) * ryt(ij) * amskt(ij, kstr)
        axfmp(ij, k) = axfmp(ij, k) &
          &          + ts * (  (ffmx(ijle, k) - ffmx(ij, k)) * rx &
          &              + (ffmy(ijln, k) - ffmy(ij, k)) * ry(ij)) * &
          &            rxt(ij) * ryt(ij) * amskt(ij, kstr)
        axdsd(ij, k) = axdsd(ij, k) &
          &          + ts * (  (fddx(ijle, k) - fddx(ij, k)) * rx &
          &              + (fddy(ijln, k) - fddy(ij, k)) * ry(ij)) * &
          &            rxt(ij) * ryt(ij) * amskt(ij, kstr)
        axdsb(ij, k) = axdsb(ij, k) &
          &          + ts * (  (fdbx(ijle, k) - fdbx(ij, k)) * rx &
          &              + (fdby(ijln, k) - fdby(ij, k)) * ry(ij)) * &
          &            rxt(ij) * ryt(ij) * amskt(ij, kstr)
     end do
   
  end do

  do k = 1, nic
     do ij = ijtstr, ijtend
        if (ax(ij, k) .gt. 0.d0) then
           hix(ij, k) = axhix(ij, k) / ax(ij, k)
           eix(ij, k) = axeix(ij, k) / ax(ij, k)
           hsx(ij, k) = axhsx(ij, k) / ax(ij, k)
           tix(ij, k) = ti(eix(ij, k)/hix(ij, k), si)
           asx(ij, k) = axasx(ij, k) / ax(ij, k)
           frlvx(ij, k) = axflv(ij, k) / ax(ij, k)
           vmpx(ij, k) = axvmp(ij, k) / ax(ij, k)
           frmpx(ij, k) = axfmp(ij, k) / ax(ij, k)
           dsdx(ij, k) = axdsd(ij, k) / ax(ij, k)
           dsbx(ij, k) = axdsb(ij, k) / ax(ij, k)
        else
           hix(ij, k) = hic(k)
           eix(ij, k) = 0.d0
           hsx(ij, k) = 0.d0
           tix(ij, k) = tmi
           asx(ij, k) = 0.d0
           frlvx(ij, k) = 1.d0
           vmpx(ij, k) = 0.d0
           frmpx(ij, k) = 0.d0
           dsdx(ij, k) = 0.d0
           dsbx(ij, k) = 0.d0
        end if
     end do
  end do

! ax(ij, 0) can be negative
  do ij = ijtstr, ijtend
     ax(ij, 0) = 1.d0
  end do
  do k = 1, nic
     do ij = ijtstr, ijtend
        ax(ij, 0) = ax(ij, 0) - ax(ij, k)
     end do
  end do

!  call rewnml(ifpar, jfpar)
!  do k = 1, nic
!     do ij = ijtstr, ijtend
!        if (ax(ij, k) .lt. 0.d0) then
!           write(jfpar, *) '### NEGATIVE AREA (ipadv) ###', &
!             &             ij, k, ax(ij, k)
!           stop
!        end if
!        if (     (axhix(ij, k) .lt. 0.d0) &
!          & .or. (axhsx(ij, k) .lt. 0.d0)) then
!           write(jfpar, *) '### NEGATIVE ICE/SNOW (ipadv) ###', &
!             &             ij, k, ax(ij, k), axhix(ij, k), axhsx(ij, k)
!           stop
!        end if
!        if (hix(ij, k) .le. 0.d0) then
!           write(jfpar, *) '### INVALID ICE THICKNESS (ipadv) ###', &
!             &             ij, k, ax(ij, k), hix(ij, k)
!           stop
!        end if
!     end do
!  end do

  return

end subroutine padvct

end module ipadv
