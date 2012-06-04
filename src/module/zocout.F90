module zocout

use zocdim, only: &
 & nxyzdm,   nxyz

implicit none
public
save

integer, parameter :: nwork = 25*nxyzdm

real(8) ::  wrkout(nwork)
real(8) ::  dbleou(nxyz)
real(4) ::  snglou(nxyz)

#ifdef OPT_PARALLEL
real(8), allocatable ::  dbloug(:, :, :)
real(4), allocatable ::  sngoug(:, :, :)
#endif

end module zocout
