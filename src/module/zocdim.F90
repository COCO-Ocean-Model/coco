module zocdim
 
  use zocnod
 
  implicit none
 
  integer, parameter :: nxg    =     32, nyg    =     32, nz     =     28
  integer, parameter :: kz     =      5
  integer, parameter :: nic    =      5
  integer, parameter :: nxgdim =     36, nygdim =     36, nzdim  =     30
  integer, parameter :: igstr  =      3, jgstr  =      3, kstr   =      2
  integer, parameter :: igend  =     34, jgend  =     34, kend   =     29
  integer, parameter :: ijgstr =     75, ijgend =   1222
 
  integer, parameter :: nxyg = nxg*nyg, nxyzg = nxyg*nz
  integer, parameter :: nxygdm = nxgdim*nygdim, nxyzgd = nxygdm*nzdim
 
  integer, parameter :: nx     =     16, ny     =     16
  integer, parameter :: nxdim  =     20, nydim  =     20
  integer, parameter :: istr   =      3, jstr   =      3
  integer, parameter :: iend   =     18, jend   =     18
  integer, parameter :: ijstr  =     43, ijend  =    358
 
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
 
  integer, parameter :: inodes =      2, jnodes =      2
  integer, parameter :: icomm = istr-1, jcomm = jstr-1
 
  integer, parameter :: ijtgst =    112, ijtged =   1185
  integer, parameter :: ijvgst =    112, ijvged =   1148
 
  logical, save :: oinit, ofinal
 
end module zocdim
