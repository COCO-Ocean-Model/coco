module zocmsk

use zocdim, only: &
 & nxydim,  nzdim

implicit none
public
save

real(8) ::  amskt(nxydim, nzdim)
real(8) ::  amftx(nxydim, nzdim),  amfty(nxydim, nzdim)
real(8) ::  amftz(nxydim, nzdim)
real(8) ::  amskv(nxydim, nzdim)
real(8) ::  amfvx(nxydim, nzdim),  amfvy(nxydim, nzdim)
real(8) ::  amfvz(nxydim, nzdim)
real(8) ::  amskb(nxydim, nzdim)
integer ::  nbot(nxydim)

#ifdef OPT_BBL
real(8) ::  amsktb(nxydim), amskt0(nxydim), amskt1(nxydim)
real(8) ::  amskvb(nxydim), amskv0(nxydim), amskv1(nxydim)
integer ::  nbotv(nxydim)
#endif

end module zocmsk
