module zqsat   ! Clausis-Clapeyron

  use zocphy, only: &
    &   epsv,    es0,     el,  emelt,  tqice,   rvap,  tmelt

  implicit none

  private
  public :: fqsat, fdqsat

contains

function fqsat(t, p)   ! saturation water vapour mixing ratio
  !$acc routine seq
  real(8) ::  fqsat  
  real(8), intent(in), value ::      t,      p
  
  fqsat = epsv * es0 / p &
     &    * exp( (el+emelt/2.d0*(1.d0-sign(1.d0,t-tqice))) &
     &           / rvap *( 1.d0/tmelt - 1.d0/t ) )
  return

end function fqsat

!#######################################################################

function fdqsat(t, qs)  ! d(qsat)/d(t)
  !$acc routine seq
  real(8) :: fdqsat
  real(8), intent(in), value :: t, qs
  
  fdqsat = (el+emelt/2.d0*(1.d0-sign(1.d0,t-tmelt))) &
     &     * qs / ( rvap * t*t )
  return

end function fdqsat

end module zqsat
