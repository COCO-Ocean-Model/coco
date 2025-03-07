module binst

! --- information -----------------------------------------------------
!
!  HISTORY
!     '02.10.10  H.Hasumi: from COCO3.4
!     '03.05.27  H.Hasumi: bug fix for non-par
!     '07.04.23  H.Hasumi
!     '09.02.23  Y.Komuro: call DDENST(diagnose R) instead of OVTURN
!     '10.04.14  M.kurogi
!     '10.04.14  M.Kurogi: (COCO4.4 tripolar code by Dr. Suzuki)
!     '12.10.22  T.Suzuki: for COCO5.0
!
! ---------------------------------------------------------------------
  implicit none
  private
  public :: iniset

contains

  subroutine iniset(                                &
   &                 uadv,   vadv,   wadv,      r,  &
   &                   ub,     vb,     tb,          &
   &                   hb,   ubtb,   vbtb,          &
   &                    w,    amv,    ahv)

  use zocdim, only : nxdim, nydim, nzdim, ntdim, nxydim, nxyzdm, nztdim
!      use zocgrd
!      use zocphy
  use zocmsk, only : amskv
  use tovtr
  use dvlva
  use bstbc
  use bchmk
  use bshft

  implicit none

  real(8), intent(out)   :: uadv(nxyzdm), vadv(nxyzdm), wadv(nxyzdm, 9), r(nxyzdm)
  real(8), intent(inout) ::   ub(nxyzdm),   vb(nxyzdm),   tb(nxyzdm, ntdim)
  real(8), intent(inout) ::   hb(nxydim), ubtb(nxydim),  vbtb(nxydim)
  real(8), intent(inout) ::    w(nxyzdm),  amv(nxyzdm),   ahv(nxyzdm)

  real(8) ::      uz(nxydim),     vz(nxydim)
! common /work/ uz, vz

  integer ::   ijk, ij, k

  call shift_pack_begin
  call shift1(  hb,       nxdim, nydim,     1,  1.d0,  0,  0)
  call shift2(ubtb, vbtb, nxdim, nydim,     1, -1.d0, -1, -1)
  call shift2(  ub,   vb, nxdim, nydim, nzdim, -1.d0, -1, -1)
  call shift1(   w,       nxdim, nydim, nzdim,  1.d0,  0,  0)
  call shift_pack_end
  call shift_unpack(  hb, 1)
  call shift_unpack(ubtb, 2)
  call shift_unpack(vbtb, 3)
  call shift_unpack(  ub, 4)
  call shift_unpack(  vb, 5)
  call shift_unpack(   w, 6)
    
  !$acc kernels default(present)
  do ijk = 1, nxyzdm
     r(ijk) = 0.0d+0
  end do
  !$acc end kernels
  
  call ddenst(  r,   tb)
  
  call shift_pack_begin
  call shift1( r, nxdim, nydim,  nzdim, 1.d0, 0, 0)
  call shift1(tb, nxdim, nydim, nztdim, 1.d0, 0, 0)
  call shift_pack_end
  call shift_unpack( r, 1)
  call shift_unpack(tb, 2)

  call stbctr( tb,      r)
#ifdef OPT_BBL
  call rmmskv
  call admkvb
  do k=1, nzdim
     do ij=1,nxydim
        ijk= ij + (k-1)*nxydim
        ub(ijk)=ub(ijk)*amskv(ij,k)
        vb(ijk)=vb(ijk)*amskv(ij,k)
     end do
  end do
  call admskv

  call rmmskt
  call admktb
#endif
  call velvad( &
   &              uadv,   vadv,   wadv, &
   &                ub,     vb,      w)
#ifdef OPT_BBL
  call velvab( wadv,   w )
  call stbbvt(   ub,  vb )
  call shift2(ub, vb, nxdim, nydim, nzdim, -1.d0, -1,  -1 )
#endif

  call shift_pack_begin
  call shift2(uadv, vadv, nxdim, nydim, nzdim, -1.d0,  0,  0)
  call shift1(wadv(:,1),  nxdim, nydim, nzdim,  1.d0, -1, -1)
  call shift1(wadv(:,2),  nxdim, nydim, nzdim,  1.d0, -1, -1)
  call shift1(wadv(:,3),  nxdim, nydim, nzdim,  1.d0, -1, -1)
  call shift1(wadv(:,4),  nxdim, nydim, nzdim,  1.d0, -1, -1)
  call shift1(wadv(:,5),  nxdim, nydim, nzdim,  1.d0, -1, -1)
  call shift1(wadv(:,6),  nxdim, nydim, nzdim,  1.d0, -1, -1)
  call shift1(wadv(:,7),  nxdim, nydim, nzdim,  1.d0, -1, -1)
  call shift1(wadv(:,8),  nxdim, nydim, nzdim,  1.d0, -1, -1)
  call shift1(wadv(:,9),  nxdim, nydim, nzdim,  1.d0, -1, -1)
  call shift1(amv,        nxdim, nydim, nzdim,  1.d0, -1, -1)
  call shift1(ahv,        nxdim, nydim, nzdim,  1.d0,  0,  0)
  call shift_pack_end
  call shift_unpack(     uadv,  1)
  call shift_unpack(     vadv,  2)
  call shift_unpack(wadv(:,1),  3)
  call shift_unpack(wadv(:,2),  4)
  call shift_unpack(wadv(:,3),  5)
  call shift_unpack(wadv(:,4),  6)
  call shift_unpack(wadv(:,5),  7)
  call shift_unpack(wadv(:,6),  8)
  call shift_unpack(wadv(:,7),  9)
  call shift_unpack(wadv(:,8), 10)
  call shift_unpack(wadv(:,9), 11)
  call shift_unpack(      amv, 12)
  call shift_unpack(      ahv, 13)

#ifdef OPT_TRIPOLE
  call excngw(wadv)
#endif
  return
end subroutine iniset
end module binst
