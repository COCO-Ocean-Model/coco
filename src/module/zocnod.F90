module zocnod

  implicit none
  public
  save

  integer :: mpi_comm_ogcm
  integer :: nprocs, myrank, igrank
  integer :: ijnode,  iroot
  integer ::   ierr
  integer ::  irank,    iup,  idown
  integer ::  jrank,    jup,  jdown
#ifdef OPT_TRIPOLE
  integer ::  jupe,    jupw
#endif
  integer :: ijtstr, ijtend
  integer :: ijvstr, ijvend

end module zocnod
