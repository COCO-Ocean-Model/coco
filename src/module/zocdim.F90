module zocdim
 
  use zocnod
 
  implicit none
 
  integer, parameter :: nxg    =    360, nyg    =    256, nz     =     63
  integer, parameter :: kz     =     10
  integer, parameter :: nic    =      5
  integer, parameter :: nxgdim =    364, nygdim =    260, nzdim  =     65
  integer, parameter :: igstr  =      3, jgstr  =      3, kstr   =      2
  integer, parameter :: igend  =    362, jgend  =    258, kend   =     64
 
  integer, parameter :: nxyg = nxg*nyg, nxyzg = nxyg*nz
  integer, parameter :: nxygdm = nxgdim*nygdim, nxyzgd = nxygdm*nzdim
 
  integer, parameter :: nx     =     10, ny     =      8
  integer, parameter :: nxdim  =     14, nydim  =     12
  integer, parameter :: istr   =      3, jstr   =      3
  integer, parameter :: iend   =     12, jend   =     10
  integer, parameter :: ijstr  =     31, ijend  =    138
 
  integer, parameter :: ntdim =      2
 
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
 
  integer, parameter :: inodes =     36, jnodes =     32
  integer, parameter :: icomm = istr-1, jcomm = jstr-1
 
 
  logical, save :: oinit, ofinal
 
end module zocdim
