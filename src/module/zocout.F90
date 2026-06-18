module zocout

  use zocdim, only: &
    & nxyzdm,   nxyz

  implicit none
  public
  save

  integer, parameter :: cnwork = 60, cnwrks = 15
  integer, parameter :: nncmax = 10
  integer, parameter :: nchmax = 4
  integer, parameter :: nwork = cnwork*nxyzdm

  real(8) ::  wrkout(nwork)
  logical, allocatable :: owrksg(:)
  real(8), allocatable :: dbleou(:)
  real(4), allocatable :: snglou(:)
  integer :: nworks
#ifndef OPT_IO_COCOMPI
  real(8), allocatable ::  dbloug(:, :, :)
  real(4), allocatable ::  sngoug(:, :, :)
#endif

  integer(4) :: loglev = 0

end module zocout
