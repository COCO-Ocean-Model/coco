module zocgrd

use zocdim, only: &
 & nxydim,  nzdim,    nic

implicit none
public
save

real(8) ::     dy(nxydim),    dym(nxydim)
real(8) ::     dz(nxydim, nzdim) ,    dzm(nxydim, nzdim)
real(8) ::    dzv(nxydim, nzdim)
real(8) ::    dz0(nzdim),      ds(nzdim) ,    dsm(nzdim)
real(8) ::     dt,     dx
real(8) ::     rx
real(8) ::     ry(nxydim),    rym(nxydim)
real(8) ::     rs(nzdim),    rsm(nzdim)
real(8) ::    hic(nic+1)

real(8) ::     tt,     ts,    tss
integer ::     nt,    its,   itst,   ntss
integer :: ieuler

real(8) ::    rea,   zbot
real(8) ::    cor(nxydim),   dept(nxydim),  rdepv(nxydim)
real(8) ::    hxt(nxydim),    hxu(nxydim)
real(8) ::    hyt(nxydim),    hyu(nxydim)
real(8) ::   hxyt(nxydim),   hxyu(nxydim)
real(8) ::   hyxt(nxydim),   hyxu(nxydim)
real(8) ::    rxt(nxydim),    rxu(nxydim)
real(8) ::    ryt(nxydim),    ryu(nxydim)

end module zocgrd
