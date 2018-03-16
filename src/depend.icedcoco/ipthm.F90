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
!     '09.05.25  Y.Komuro: CMIP5 output code included
!                          (basal/lateral melting processes separated)
!     '12.07.30  Y.Komuro: for COCO5.0
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
  use zocite

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
  real(8) :: laxhix(nxydim, 0:nic), laxeix(nxydim, 0:nic)
  real(8) :: daxeit(nxydim, 0:nic), daxeib(nxydim, 0:nic)
!      common /work/ az, axhix, axhsx, axhsxn, wai,
!     &              wi, ws, wen, wsn, axeix, axeixn,
!     &              laxhix, laxeix, daxeit, daxeib

  real(8), save ::    rri,    rrs, rorirs
  real(8), save ::  rsfus
  real(8), save ::    tmi
  real(8), save :: epsaei
  logical, save :: ofirst = .true.

  real(8) ::   wres
  real(8) ::   hsxo,    dhs
  real(8) :: ax1max
  real(8) :: daxhix, daxeix, daxhit, daxhib,   eieq
  integer ::     ij,      k
  integer ::  ifpar,  jfpar,  istat

  real(8), save ::    eps = 1.0d-3,   epsl = 1.0d-6

  real(8), save ::  amin = 1.0d-6,  amax = 1.0d0,  si = 5.0d0
  integer, save ::  mic = nic

  namelist /nmamin/ amin, amax, mic
  namelist /nmislt/ si


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

! *** lateral ice formation/melting
  do k = 1, nic
     do ij = ijtstr, ijtend
        laxhix(ij, k) = ax(ij, k) * hix(ij, k)
        laxeix(ij, k) = ax(ij, k) * eix(ij, k) &
          &           + wil(ij, k) * ts / rhoi &
          &           * amskt(ij, kstr)
        axeixn(ij, k) = axeixn(ij, k) &
          &           + wil(ij, k) * ts / rhoi &
          &           * amskt(ij, kstr)
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
     ax1max = az(ij, 0) + az(ij, 1)
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
