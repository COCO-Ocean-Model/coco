module qckag

! --- information -----------------------------------------------------
!
!  Aggregate variables from various modules and pass them to CHEKIN
!
!  HISTORY
!     '21.03.26  Y.Komuro: for COCO5.0
!
! ---------------------------------------------------------------------

  use zocdim, only: &
    &     nx,     ny, nxydim,  nzdim,  ntdim,   kstr,    nic, oinit, ofinal
  use zocmsk, only: &
    &  amskt,  amskv

  implicit none

! Variables for chkstk and related subroutines
  real(8), save ::  ofprec(nxydim) = 0.0d0, ofsnow(nxydim) = 0.0d0
  real(8), save ::  ofevap(nxydim) = 0.0d0, ofroff(nxydim) = 0.0d0
  real(8), save ::  ofsoff(nxydim) = 0.0d0, ofwiws(nxydim) = 0.0d0
  real(8), save ::  ofwnml(nxydim) = 0.0d0, ofsrst(nxydim) = 0.0d0
  real(8), save ::     sst(nxydim) = 0.0d0
  real(8), save ::  oftocn(nxydim) = 0.0d0, ofwocn(nxydim) = 0.0d0
  real(8), save ::  ofsocn(nxydim) = 0.0d0
  real(8), save ::  oftaux(nxydim) = 0.0d0, oftauy(nxydim) = 0.0d0
  real(8), save ::  swcnv1(nxydim) = 0.0d0, oswabs(nxydim) = 0.0d0
  logical, save ::  of_acc_enter_data=.true.
  private

  public :: chksfx, cofpfw, cofpwi, cofpnw, cofpsr, cofpsf, &
    &       cofptu, copswc

contains

subroutine acc_enter_data
  implicit none
  if(.not. of_acc_enter_data) return  
  of_acc_enter_data=.false.
  !$acc enter data create(ofprec, ofsnow, ofevap, ofroff)
  !$acc enter data create(ofsoff, ofwiws, ofwnml, ofsrst, sst)
  !$acc enter data create(oftocn, ofwocn, ofsocn)
  !$acc enter data create(oftaux, oftauy, swcnv1, oswabs)
end subroutine acc_enter_data
    
subroutine chksfx
  use zocphy, only:    cpo,  rhoo
  use qckot
  
  real(8) :: ofithm(nxydim), ofhfds(nxydim)
  real(8) :: ofvsfc(nxydim), ofsfdi(nxydim)
  integer ::     ij
  logical, save :: of=.true.

  call acc_enter_data
  if (of) then
     of=.false.
     !$acc enter data create(ofithm, ofhfds, ofvsfc, ofsfdi)
  end if

  !$acc kernels default(present)
  do ij = 1, nxydim
!    OFITHM: calculated by residual, not from OFWIWS,
!            because ice-related FW is also added in ICTRNS
     ofithm(ij) = ( ofwocn(ij) - &
       &       ( ofprec(ij) + ofsnow(ij) + ofevap(ij) &
       &       + ofroff(ij) + ofsoff(ij) ) ) &
       &       * amskt(ij, kstr)
!    OFHFDS: net downwelling heat flux
     ofhfds(ij) = ( rhoo * cpo * &
       &            ( oftocn(ij) + sst(ij) * ofwocn(ij)) &
       &          + oswabs(ij) * (1.0d0 - swcnv1(ij)) ) &
       &        * amskt(ij, kstr)
!    OFVSFC,OFSFDI: 1 [psu g(water)/cm2 s] = 1.0D-3 [g(salt)/cm2 s]
     ofvsfc(ij) = 1.0d-3 * rhoo * ofsrst(ij) * amskt(ij, kstr)
     ofsfdi(ij) = 1.0d-3 * rhoo * ofsocn(ij) * amskt(ij, kstr)
  end do
  !$acc end kernels
  
  call chekin( ofprec, 'OFPREC', &
    &          'Precipitation entering into the ocean', 'cm/s', &
    &              nx,     ny,      1, nxydim, 'OCSFCT')      
  call chekin( ofsnow, 'OFSNOW', &
    &          'Snow directly entering into the ocean', 'cm/s', &
    &              nx,     ny,      1, nxydim, 'OCSFCT')      
  call chekin( ofroff, 'OFROFF', &
    &          'Runoff entering into the ocean', 'cm/s', &
    &              nx,     ny,      1, nxydim, 'OCSFCT')      
  call chekin( ofsoff, 'OFSOFF', &
    &          'Snow runoff entering into the ocean', 'cm/s', &
    &              nx,     ny,      1, nxydim, 'OCSFCT')      
  call chekin( ofevap, 'OFEVAP', &
    &          'Evaporation from the ocean (positive downward)', &
    &          'cm/s', &
    &              nx,     ny,      1, nxydim, 'OCSFCT')      
  call chekin( ofithm, 'OFITHM', &
    &  'Fw induced by sea-ice thermodynamics (positive downward)', &
    &          'cm/s', &
    &              nx,     ny,      1, nxydim, 'OCSFCT')      
  call chekin( ofwnml, 'OFWNML', &
    &  'Fw for normalizing surface height (positive downward)', &
    &          'cm/s', &
    &              nx,     ny,      1, nxydim, 'OCSFCT')      
  call chekin( ofhfds, 'OFHFDS', &
    &  'Total heat flux into the ocean (positive downward)', &
    &          'erg/cm2', &
    &              nx,     ny,      1, nxydim, 'OCSFCT')      
  call chekin( ofvsfc, 'OFSVSF', &
    &          'Fs by SSS restoring, entering into the ocean', &
    &          'g(salt)/cm2/s', &
    &              nx,     ny,      1, nxydim, 'OCSFCT')      
  call chekin( ofsfdi, 'OFSFDI', &
    &          'Fs by ice processes, entering into the ocean', &
    &          'g(salt)/cm2/s', &
    &              nx,     ny,      1, nxydim, 'OCSFCT')      
  call chekin( oftaux, 'OFTAUX', &
    &          'X-stress entering into the ocean', 'dyn/cm2', &
    &              nx,     ny,      1, nxydim, 'OCSFCV')      
  call chekin( oftauy, 'OFTAUY', &
    &          'Y-stress entering into the ocean', 'dyn/cm2', &
    &              nx,     ny,      1, nxydim, 'OCSFCV')      
  
  return
end subroutine chksfx
! ======================================================================
subroutine cofpfw( &
  &                 prec,   snow,   roff,   soff,   evap, &
  &                   ax)
  
  real(8) ::   prec(nxydim),   snow(nxydim)
  real(8) ::   evap(nxydim),   roff(nxydim)
  real(8) ::   soff(nxydim)
  real(8) ::     ax(nxydim, 0:nic)
  integer ::     ij
  
  if (oinit .or. ofinal) return
  call acc_enter_data 
  !$acc kernels default(present)
  do ij = 1, nxydim
!    fluxes are positive when entering the ocean (i.e., downward)
!    SNOW = SNOW_original + SOFF_original
!    SNOW, SOFF: portion falling to ice-top must be removed
!                i.e., SNOW(SOFF) = (1-A) * SNOW(SOFF)_original
     ofprec(ij) = prec(ij) * amskt(ij, kstr)
     ofsnow(ij) = ax(ij, 0) * &
       &          (snow(ij) - soff(ij)) * amskt(ij, kstr)
     ofroff(ij) = roff(ij) * amskt(ij, kstr)
     ofsoff(ij) = ax(ij, 0) * soff(ij) * amskt(ij, kstr)
     ofevap(ij) = - evap(ij) * amskt(ij, kstr)

!    FWNMD has been added to PREC/WEV in NMLPER.
!    That portion must be removed from the output value.
     ofprec(ij) = ofprec(ij) - max(ofwnml(ij), 0.0d0) &
       &        * amskt(ij, kstr)
     ofevap(ij) = ofevap(ij) - min(ofwnml(ij), 0.0d0) &
       &        * amskt(ij, kstr)
  end do
  !$acc end kernels

  return
end subroutine cofpfw
! ======================================================================
subroutine cofpwi( &
  &                   ws,     wi)

  real(8) ::     wi(nxydim),     ws(nxydim)
  integer ::     ij

  call acc_enter_data 
  !$acc kernels default(present)
  do ij = 1, nxydim
!    fluxes are positive when entering the ocean (i.e., downward)
     ofwiws(ij) = - (wi(ij) + ws(ij)) * amskt(ij, kstr)
  end do
  !$acc end kernels
  
  return
end subroutine cofpwi
! ======================================================================
subroutine cofpnw( &
  &                fwnmd)

  real(8) ::  fwnmd(nxydim)
  integer ::     ij

  call acc_enter_data
  !$acc kernels default(present)
  do ij = 1, nxydim
!    fluxes are positive when entering the ocean (i.e., downward)
     ofwnml(ij) = - fwnmd(ij) * amskt(ij, kstr)
  end do
  !$acc end kernels
  
  return
end subroutine cofpnw
! ======================================================================
subroutine cofpsr( &
  &                fsrst)

  real(8) ::  fsrst(nxydim)
  integer ::     ij

  call acc_enter_data
  !$acc kernels default(present)
  do ij = 1, nxydim
!    fluxes are positive when entering the ocean (i.e., downward)
     ofsrst(ij) = fsrst(ij) * amskt(ij, kstr)
  end do
  !$acc end kernels
  return
end subroutine cofpsr
! ======================================================================
subroutine cofpsf( &
  &                   tx,     ft,     fs,  swabs)

  real(8) ::     tx(nxydim, nzdim, ntdim)
  real(8) ::     ft(nxydim, ntdim)
  real(8) ::     fs(nxydim)
  real(8) ::  swabs(nxydim)
  integer ::     ij

  call acc_enter_data
  !$acc kernels default(present)
  do ij = 1, nxydim
     sst(ij) = tx(ij, kstr, 1) * amskt(ij, kstr)
!    FT(IJ, 1): positive downward
!    FT(IJ, 2), FS: positive upward
!    OF?OCN: positive downward
!    SWABS: total SW flux into the ocn
!           OFTOCN includes SW only into the uppermost level
     oftocn(ij) = ft(ij, 1) * amskt(ij, kstr)
     ofwocn(ij) = -ft(ij, 2) * amskt(ij, kstr)
     ofsocn(ij) = -fs(ij) * amskt(ij, kstr)
     oswabs(ij) = swabs(ij) * amskt(ij, kstr)
  end do
  !$acc end kernels
  return
end subroutine cofpsf
! ======================================================================
subroutine cofptu( &
  &                 taux,   tauy)

  real(8) ::   taux(nxydim),   tauy(nxydim)
  integer ::     ij

  call acc_enter_data
  !$acc kernels default(present)
  do ij = 1, nxydim
!    fluxes are positive when entering the ocean (i.e., downward)
     oftaux(ij) = taux(ij) * amskv(ij, kstr)
     oftauy(ij) = tauy(ij) * amskv(ij, kstr)
  end do
  !$acc end kernels
  return
end subroutine cofptu
! =====================================================================
subroutine copswc( &
  &               swconv)

  real(8) :: swconv(nxydim)
  integer ::     ij

  call acc_enter_data
  !$acc kernels default(present)
  do ij = 1, nxydim
     swcnv1(ij) = swconv(ij)
  end do
  !$acc end kernels
  return
end subroutine copswc
  
end module qckag
