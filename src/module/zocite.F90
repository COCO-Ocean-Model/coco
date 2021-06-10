module zocite
! sea ice thermal energy

  use zocphy, only: &
    &    cpi,   dtds,   hfus

  implicit none

  private
  public :: ei, ti

contains

function ei(tice, sice)
  real(8)             ::     ei
  real(8), intent(in) ::   tice,   sice
  real(8), parameter  ::    eps = 1.0d-2

  ei = cpi * (dtds * sice - tice) &
     & + hfus * (1.d0 - dtds * sice / tice)
  ei = sign(max(abs(ei), eps), ei)

  return

end function ei

!#######################################################################

function ti(eice, sice)
  real(8)             ::     ti
  real(8), intent(in) ::   eice,   sice

  ti = (  hfus + cpi * dtds * sice - eice &
     & - sqrt(  (hfus + cpi * dtds * sice - eice)**2 &
     &          - 4.d0 * cpi * hfus * dtds * sice)) &
     &   / cpi * 0.5d0

  return

end function ti

end module zocite
