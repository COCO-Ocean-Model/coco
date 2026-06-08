module dvdif

! --- information -----------------------------------------------------
!
!  Setting vertical diffusion and viscosity coefficients
!
!  HISTORY
!     '97.03.23  H.Hasumi
!     '98.05.08  H.Hasumi: revision of the consistency with the
!                          enhanced diffusivity version for the
!                          parameterization of the vertical convection
!     '98.05.11  H.Hasumi: namelist
!     '02.05.29  H.Nakano: tracer dimension
!     '07.04.23  H.Hasumi
!     '08.06.11  H.Hasumi: initial/final processing
!     '08.07.10  H.Hasumi: initial/final processing
!     '12.08.02  Y.Komuro: for COCO5.0
!
! ---------------------------------------------------------------------

  use zocdim, only: &
    & nxydim,  nzdim,  ntdim, &
    &   kstr,   kend,     nz, &
    &  ijstr,  ijend, &
    &  oinit, ofinal
  use zocmsk, only: &
#ifdef OPT_BBL
    & amsktb, amskvb, nbotv,  &
#endif
    &   nbot

  implicit none
#include "coco.h"
  private

  public :: vdiff, puttao
#ifdef OPT_BBL
  public :: vdiffb
#endif

contains

subroutine vdiff( &
  &                  amv,    ahv, &
  &                   uy,     vy,      r,   taux,   tauy, &
  &                   ty,     hy )
  use ufile
  real(8), intent(out) ::     amv(nxydim, nzdim),    ahv(nxydim, nzdim)
  real(8), intent(in)  ::      uy(nxydim, nzdim),     vy(nxydim, nzdim)
  real(8), intent(in)  ::       r(nxydim, nzdim)
  real(8), intent(in)  ::    taux(nxydim)       ,   tauy(nxydim)
  real(8), intent(in)  ::      ty(nxydim, nzdim, ntdim)
  real(8), intent(in)  ::      hy(nxydim)

  integer ::     ij,      k
  integer ::  ifpar,  jfpar,  istat

  logical, save :: ofirst = .true.

  real(8), save ::  amv0(nz) = 0.d0,  ahv0(nz) = 0.d0

  namelist /nmvisv/ amv0
  namelist /nmdifv/ ahv0

  if (oinit .or. ofinal) then
     return
  end if

  if (ofirst) then
     ofirst = .false.
     READ_NAMELIST( nmvisv )
     READ_NAMELIST( nmdifv )
  end if

  do k = kstr, kend
     do ij = 1, nxydim
        amv(ij, k) = amv0(k-kstr+1)
        ahv(ij, k) = ahv0(k-kstr+1)
     end do
  end do

  return

end subroutine vdiff
! =====================================================================
subroutine puttao( &
  &                 taox,   taoy, &
  &                 caic,   cais)

  use bshft

  real(8), intent(in) ::   taox(nxydim),   taoy(nxydim)
  real(8), intent(in) ::   caic,   cais

! dummy routine

  return

end subroutine puttao


#ifdef OPT_BBL
! *********************************************************************

subroutine vdiffb( &
  &                   amv,    ahv )
  
! --- information -----------------------------------------------------
!
!  Vertical viscosity and diffusion coefficients for the bottom
! boundary layer.
!
!  HISTORY
!     '01.02.08  H.Hasumi
!     '12.08.02  Y.Komuro: for COCO5.0
!
! ---------------------------------------------------------------------
  use ufile
  real(8), intent(inout) ::    amv(nxydim, nzdim),    ahv(nxydim, nzdim)

  integer ::     ij
  integer ::     kt,     kv
  integer ::  ifpar,  jfpar,  istat

  logical, save :: ofirst = .true.

  real(8), save :: amvbbl = 0.d0, ahvbbl = 0.d0
  namelist /nmbbdv/ amvbbl, ahvbbl

  if (oinit .or. ofinal) then
     return
  end if

  if (ofirst) then
     ofirst = .false.
     READ_NAMELIST( nmbbdv )
  end if

  do ij = ijstr, ijend
     kt = nbot(ij)
     kv = nbotv(ij)
     amv(ij, kv) = max(amvbbl, amv(ij, kv)) * amskvb(ij) &
       &         + amv(ij, kv) * (1.d0 - amskvb(ij))
     ahv(ij, kt) = max(ahvbbl, ahv(ij, kt)) * amsktb(ij) &
       &         + ahv(ij, kt) * (1.d0 - amsktb(ij))
  end do

  do ij = ijstr, ijend
     kt = nbot(ij)
     kv = nbotv(ij)
     amv(ij, kend) = amv(ij, kv)
     ahv(ij, kend) = ahv(ij, kt)
  end do

  return

end subroutine vdiffb
#endif

end module dvdif
