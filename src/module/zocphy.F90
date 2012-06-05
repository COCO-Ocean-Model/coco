module zocphy

  implicit none
  public
  save

! Physical constants in SI units
  real(8), parameter :: ER = 6.37D6, GRAV = 9.8D0
  real(8), parameter :: CP = 1004.6, RAIR = 287.04
  real(8), parameter :: EL = 2.5D6, EMELT = 3.4D5
  real(8), parameter :: CPVAP = 1810., RVAP = 461.
  real(8), parameter :: DWATR = 1.D3
  real(8), parameter :: ES0 = 611.
  real(8), parameter :: STB = 5.67D-8
  real(8), parameter :: FKARM = 0.4
  real(8), parameter :: TMELT = 273.15D0
  real(8), parameter :: CPWATR = 4200., CPICE = 2000.
  real(8), parameter :: TFRZS = 271.35, TQICE = 273.15
  real(8), parameter :: EPSV = RAIR/RVAP, EPSVT = 1.D0/EPSV-1.D0

! Physical constants in cgs units
  real(8), parameter :: gravit = 9.8D+2
  real(8), parameter :: rhoo = 1.D0, rhoi = 0.9D0, rhos = 0.33D0, rhow = 1.D0
  real(8), parameter :: ckarm = 0.4D0
  real(8), parameter :: kelvin = 273.15D0
  real(8), parameter :: hfus = EMELT * 1.D4
  real(8), parameter :: cpo = 3.990D7, cpi= 2.093D7
  real(8), parameter :: cdi = 2.04D5, cds = 0.31D5
  real(8), parameter :: dtds = -0.0543D0, dtdz = -7.59D-6

end module zocphy
