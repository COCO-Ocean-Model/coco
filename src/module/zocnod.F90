module zocnod

implicit none
public
save

integer :: nprocs, myrank
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
