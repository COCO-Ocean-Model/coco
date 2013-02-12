module iprdc

! --- information -----------------------------------------------------
!
!  HISTORY
!     '02.10.10  H.Hasumi: from MIROC3.1
!     '07.04.23  H.Hasumi
!     '07.09.25  H.Hasumi: multi-category sea ice
!     '07.09.28  H.Hasumi: 1-layer sea ice thermodynamics
!     '08.07.10  H.Hasumi: initial/final processing
!     '09.02.20  Y.Komuro: parallel forward time march
!     '09.05.25  Y.Komuro: CMIP5 output code included
!     '09.09.04  Y.Komuro: extra output code (FEX/FEY)
!     '09.09.26  Y.Komuro: bug fix
!     '10.04.14  M.Kurogi: (COCO4.4 tripolar code by Dr. Suzuki)
!     '12.07.30  Y.Komuro: for COCO5.0
!     '13.02.12  Y.Komuro: remove non-parallel code 
!
! ---------------------------------------------------------------------

  use zocdim, only: &
    & nxydim,  nxdim,  nydim,  nzdim,  ntdim,    nic, nxyidm, &
    &     nx,     ny, &
    &  oinit, ofinal, &
    & myrank, ijnode
  use zocphy, only: &
    & gravit,   rhoi,   rhos

  implicit none

  private

  public :: predci

contains

subroutine predci( &
  &                    ax,    hix,    uix,    vix,    tix,    hsx, &
  &                    ft,     fs,   taux,   tauy,   ptop, &
  &                    ay,    hiy,    uiy,    viy,    tiy,    hsy, &
  &                    tx,     ux,     vx,     hx,     hy, &
  &                   qao,    qai,    qii,    qio,  swabs, &
  &                   wev,    wsb, &
  &                  prec,   snow,   roff,   soff, &
  &                tauaix, tauaiy, tauaox, tauaoy )

  use ictrn
  use ifhea
  use ifwtr
  use ipadv
  use ipmmt
  use iprdg
  use ipthm
  use iptmp
  use ucloc
  use qckot
  use bshft

  real(8), intent(inout) ::     ax(nxydim, 0:nic)
  real(8), intent(inout) ::    hix(nxydim, 0:nic)
  real(8), intent(inout) ::    uix(nxydim)
  real(8), intent(inout) ::    vix(nxydim)
  real(8), intent(inout) ::    tix(nxydim, 0:nic)
  real(8), intent(inout) ::    hsx(nxydim, 0:nic)
  
  real(8), intent(out)   ::     ft(nxydim, ntdim),     fs(nxydim)
  real(8), intent(out)   ::   taux(nxydim),   tauy(nxydim)
  real(8), intent(out)   ::   ptop(nxydim)

  real(8), intent(in)    ::     ay(nxydim, 0:nic)
  real(8), intent(in)    ::    hiy(nxydim, 0:nic)
  real(8), intent(in)    ::    uiy(nxydim)
  real(8), intent(in)    ::    viy(nxydim)
  real(8), intent(in)    ::    tiy(nxydim, 0:nic)
  real(8), intent(in)    ::    hsy(nxydim, 0:nic)
  real(8), intent(in)    ::     tx(nxydim, nzdim, ntdim)
  real(8), intent(in)    ::     ux(nxydim, nzdim),     vx(nxydim, nzdim)
  real(8), intent(in)    ::     hx(nxydim),     hy(nxydim)
  real(8), intent(inout) ::    qio(nxydim, nic),    qai(nxydim, nic)
  real(8), intent(inout) ::    qii(nxydim, nic)
  real(8), intent(inout) ::    qao(nxydim)
  real(8), intent(in)    ::  swabs(nxydim)
  real(8), intent(in)    ::    wev(nxydim)
  real(8), intent(inout) ::    wsb(nxydim, nic)
  real(8), intent(inout) ::   prec(nxydim),   snow(nxydim)
  real(8), intent(in)    ::   roff(nxydim),   soff(nxydim)
  real(8), intent(inout) :: tauaix(nxydim), tauaiy(nxydim)
  real(8), intent(inout) :: tauaox(nxydim), tauaoy(nxydim)

  real(8), save ::    wao(nxydim)
  real(8), save ::    wio(nxydim, nic),    was(nxydim, nic)
  real(8), save ::    wil(nxydim, nic)
  real(8), save ::   subi(nxydim, nic)
  real(8), save ::   evap(nxydim), adjlat(nxydim)
  real(8), save ::     az(nxydim, 0:nic),    hiz(nxydim, 0:nic)
  real(8), save ::    hsz(nxydim, 0:nic)
  real(8), save ::   pice(nxydim)

  real(8), save ::    eix(nxydim, 0:nic),    eiz(nxydim, 0:nic)

  real(8), save ::    fix(nxydim, 0:nic),    fiy(nxydim, 0:nic)
  real(8), save ::    fsx(nxydim, 0:nic),    fsy(nxydim, 0:nic)
  real(8), save ::    fex(nxydim, 0:nic),    fey(nxydim, 0:nic)
  real(8), save ::  ftitd(nxydim)
  real(8), save :: igrfra(nxydim), igrcon(nxydim), igrsni(nxydim)
  real(8), save :: inrlat(nxydim)
  real(8), save :: imrsno(nxydim), imrisf(nxydim), imribs(nxydim)

!! for check
!  real(8) :: imrtot(nxydim), igrtot(nxydim), itrtot(nxydim)

  integer ::     ij,      l

  call clcstr('ICE')

  if (      (myrank .ge. ijnode) &
    & .and. (.not. oinit) .and. (.not. ofinal)) then 
     return
  end if

  do ij = 1, nxydim
     evap(ij) = wev(ij)
     fs  (ij) = 0.d0
     taux(ij) = 0.d0
     tauy(ij) = 0.d0
     az  (ij, 0) = 1.d0
     hiz (ij, 0) = 0.d0
     hsz (ij, 0) = 0.d0
  end do
  do l = 1, nic
     do ij = 1, nxydim
        subi(ij, l) = 0.d0
        az  (ij, l) = 0.d0
        hiz (ij, l) = 0.d0
        hsz (ij, l) = 0.d0
     end do
  end do
  do l = 1, 2
     do ij = 1, nxydim
        ft(ij, l) = 0.d0
     end do
  end do

  call icetmp( &
    &            eix, &
    &            tix, &
    &            qao,    qai,    qio,    qii, &
    &             ax,    hix )

  call fiheat( &
    &             ft, &
    &            wao,    wio,    was,    wil, &
    &             ax,     tx,     hx, &
    &            qao,    qai,    qio,    qii,  swabs )
  call fwater( &
    &             ax,    hix,    hsx, &
    &           prec,   snow, &
    &           evap,   subi, adjlat, &
    &            wev,    wsb,   soff )
  call ptherm( &
    &             ax,    hix,    hsx, &
    &            eix,    tix, &
    &           prec,   snow,     ft,     fs, &
    &          ftitd, igrfra, igrcon, igrsni, &
    &         inrlat, imrsno, imrisf, imribs, &
    &             tx, &
    &            wio,    wao,    was,    wil, &
    &           evap,   subi,   roff, adjlat, &
    &            qio )
  call ictrns( &
    &             ax,    hix,    hsx,    eix,    tix, &
    &             ft,     fs )

#ifdef OPT_TRIPOLE
  call shift3( &
    &             ax,    hix,    hsx, &
    &          nxdim,  nydim,  nic+1, &
    &          1.0d0,      0,      0 )
  call shift2( &
    &            eix,    tix, &
    &          nxdim,  nydim,  nic+1, &
    &          1.0d0,      0,      0 )
#else
  call shift3( &
    &             ax,    hix,    hsx, &
    &          nxdim,  nydim,  nic+1 )
  call shift2( &
    &            eix,    tix, &
    &          nxdim,  nydim,  nic+1 )
#endif

  call padvct( &
    &             ax,    hix,    eix,    hsx,    tix, &
    &             az,    hiz,    eiz,    hsz, &
    &            fix,    fiy,    fsx,    fsy, &
    &            fex,    fey, &
    &            uix,    vix )
  call pridge( &
    &           pice, &
    &             ax,    hix,    eix,    hsx,    tix, &
    &             az,    hiz,    eiz,    hsz, &
    &            uix,    vix )
  call icadjs( &
    &             ax,    hix,    hsx,    eix )
  call ictrns( &
    &             ax,    hix,    hsx,    eix,    tix, &
    &             ft,     fs )
#ifdef OPT_TRIPOLE
  call shift3( &
    &             ax,    hix,    hsx, &
    &          nxdim,  nydim,  nic+1, &
    &          1.0d0,      0,      0 )
  call shift2( &
    &            eix,    tix, &
    &          nxdim,  nydim,  nic+1, &
    &          1.0d0,      0,      0 )
#else
  call shift3( &
    &             ax,    hix,    hsx, &
    &          nxdim,  nydim,  nic+1 )
  call shift2( &
    &            eix,    tix, &
    &          nxdim,  nydim,  nic+1 )
#endif

  call clcstr('ICEDYN')
  call pmomnt( &
    &            uix,    vix, &
    &           taux,   tauy, &
    &             ax,     ay,     az, &
    &            hix,    hiy,    hiz, &
    &            hsx,    hsy,    hsz, &
    &            uiy,    viy,   pice, &
    &             ux,     vx,     hy,   ptop, &
    &         tauaix, tauaiy, tauaox, tauaoy )
#ifdef OPT_TRIPOLE
  call shift2( &
    &           taux,    tauy, &
    &          nxdim,   nydim,     1, &
    &          -1.d0,      -1,    -1 )
#else
  call shift2( &
    &           taux,   tauy, &
    &          nxdim,  nydim,      1 )
#endif
  call clcend('ICEDYN')

  if (myrank .ge. ijnode) then
     return
  end if

  do l = 1, nic
     do ij = 1, nxydim
        ptop(ij) = ptop(ij) &
          &      + gravit * ay(ij, l) * &
          &        (rhoi * hiy(ij, l) + rhos * hsy(ij, l))
     end do
  end do

  call clcend('ICE')

! output section for CMIP5
! FIX, FIY: eastward/northward sea ice transport
!           unit [cm^3(ice)/s]
  call chekin(   fix,  'FIX', &
    &             nx,     ny,    nic, nxyidm, 'ICE')
  call chekin(   fiy,  'FIY', &
    &             nx,     ny,    nic, nxyidm, 'ICE')
! FISX, FISY: eastward/northward snow transport
!             unit [cm^3(snow)/s]
  call chekin(   fsx, 'FISX', &
    &             nx,     ny,    nic, nxyidm, 'ICE')
  call chekin(   fsy, 'FISY', &
    &             nx,     ny,    nic, nxyidm, 'ICE')
! FTITD: heat flux into sea water due to sea ice thermodynamics
!        unit [erg/cm^2/s]
  call chekin(  ftitd,  'FTITD', &
    &              nx,     ny,      1, nxydim, 'SFC')
! IGRFRA: frazil sea ice growth rate, unit [g/cm^2/s]
  call chekin( igrfra, 'IGRFRA', &
    &              nx,     ny,      1, nxydim, 'SFC')
! IGRCON: congelation sea ice growth rate, unit [g/cm^2/s]
  call chekin( igrcon, 'IGRCON', &
    &              nx,     ny,      1, nxydim, 'SFC')
! IGRSNI: snow-ice formation rate, unit [g/cm^2/s]
  call chekin( igrsni, 'IGRSNI', &
    &              nx,     ny,      1, nxydim, 'SFC')
! INRLAT: lateral sea ice net growth rate, unit [g/cm^2/s]
!  (COCO4.4 represents lateral melting process but not freezing,
!   thus this value will be zero or negative.)
  call chekin( inrlat, 'INRLAT', &
    &              nx,     ny,      1, nxydim, 'SFC')
! IMRSNO: snow melt rate, unit [g/cm^2/s]
  call chekin( imrsno, 'IMRSNO', &
    &              nx,     ny,      1, nxydim, 'SFC')
! IMRISF: rate of melt at upper surface of sea ice, unit [g/cm^2/s]
  call chekin( imrisf, 'IMRISF', &
    &              nx,     ny,      1, nxydim, 'SFC')
! IMRIBS: rate of melt at sea ice base, unit [g/cm^2/s]
  call chekin( imribs, 'IMRIBS', &
    &              nx,     ny,      1, nxydim, 'SFC')

! extra output section
! FEX, FEY: eastward/northward thermal transport by sea ice
!           unit [(erg/g)*(cm^3(ice)/s)]
  call chekin(   fex,  'FIEX', &
    &             nx,     ny,    nic, nxyidm, 'ICE')
  call chekin(   fey,  'FIEY', &
    &             nx,     ny,    nic, nxyidm, 'ICE')

!!     for check: not necessary for CMIP5 output
!  do ij = 1, nxydim
!     igrtot(ij) = igrfra(ij)+igrcon(ij)+igrsni(ij)
!     imrtot(ij) = imrsno(ij)+imrisf(ij)+imribs(ij)-inrlat(ij)
!     itrtot(ij) = igrtot(ij)-imrtot(ij)
!  end do
!  call chekin( igrtot, 'IGRTOT', &
!    &              nx,     ny,      1, nxydim, 'SFC')
!  call chekin( imrtot, 'IMRTOT', &
!    &              nx,     ny,      1, nxydim, 'SFC')
!  call chekin( itrtot, 'ITRTOT', &
!    &              nx,     ny,      1, nxydim, 'SFC')

  return

end subroutine predci

end module iprdc
