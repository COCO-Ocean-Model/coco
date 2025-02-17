module zocdim
 
  use zocnod
 
  implicit none
 
  integer, parameter :: nxg    =    360, nyg    =    256, nz     =     63
  integer, parameter :: kz     =     10
  integer, parameter :: nic    =      5
  integer, parameter :: inodes =     36, jnodes =     32
 
  integer, parameter :: igstr  =      3, jgstr  =      3, kstr   =      2
  integer, parameter :: nxgdim = nxg + 2*(igstr-1)
  integer, parameter :: nygdim = nyg + 2*(jgstr-1)
  integer, parameter ::  nzdim =  nz + 2*(kstr-1)
  integer, parameter :: igend = igstr + nxg-1
  integer, parameter :: jgend = jgstr + nyg-1
  integer, parameter ::  kend =  kstr + nz-1
 
  integer, parameter :: nxyg = nxg*nyg
  integer, parameter :: nxygdm = nxgdim*nygdim
  integer(8), parameter :: iz = nz, izdim = nzdim
  integer(8), parameter :: ixyg = nxyg, ixygdm = nxygdm
  integer(8), parameter :: nxyzg = ixyg*iz, nxyzgd = ixygdm*izdim
 
  integer, parameter :: nx = nxg/inodes
  integer, parameter :: ny = nyg/jnodes
  integer, parameter :: istr   =      3, jstr   =      3
  integer, parameter :: nxdim = nx + 2*(istr-1)
  integer, parameter :: nydim = ny + 2*(jstr-1)
  integer, parameter :: iend = istr + nx -1
  integer, parameter :: jend = jstr + ny -1
  integer, parameter :: ijstr = istr + (jstr-1)*nxdim
  integer, parameter :: ijend = iend + (jend-1)*nxdim
 
! *** Write here the number of tracer.
! *** The minimum number is 2, for temperature and salinity.
  integer, parameter :: ntdim =      2
 
! ***
  integer, parameter :: nxy = nx*ny, nxyz = nxy*nz
  integer, parameter :: nxydim = nxdim*nydim, nxyzdm = nxydim*nzdim
  integer, parameter :: nztdim = nzdim*ntdim
  integer, parameter :: nxyidm = nxydim*(nic+1)
 
  integer, parameter :: nwrk = nxyzdm*(ntdim*3+2)
 
  integer, parameter :: lng = nxgdim, lsg = -nxgdim
  integer, parameter :: lneg = nxgdim+1, lnwg = nxgdim-1
  integer, parameter :: lseg = -nxgdim+1, lswg = -nxgdim-1
  integer, parameter :: lssg = -2 * nxgdim
 
  integer, parameter :: le = 1, lw = -1, ln = nxdim, ls = -nxdim
  integer, parameter :: lne = nxdim+1, lnw = nxdim-1
  integer, parameter :: lse = -nxdim+1, lsw = -nxdim-1
  integer, parameter :: lww = -2, lss = -2 * nxdim
 
  integer, parameter :: icomm = istr-1, jcomm = jstr-1
 
 
  logical, save :: oinit, ofinal
 
end module zocdim
