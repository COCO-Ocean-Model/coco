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

#ifdef OPT_TRIPOLE
  call shift1(    hb,                       &
   &               nxdim,  nydim,      1,   &
   &                1.d0,      0,      0 )
  call shift2(  ubtb,   vbtb,               &
   &               nxdim,  nydim,      1,   &
   &               -1.d0,     -1,     -1 )
  call shift2(    ub,     vb,               &
   &               nxdim,  nydim,  nzdim,   &
   &               -1.d0,     -1,     -1 )
  call shift1(     w,                       &
   &               nxdim,  nydim,  nzdim,   &
   &                1.d0,      0,      0 )
#else
  call shift3(                              &
   &                  hb,   ubtb,   vbtb,   &
   &               nxdim,  nydim,      1)
  call shift3(                              &
   &                  ub,     vb,      w,   &
   &               nxdim,  nydim,  nzdim)
#endif
  !$acc kernels default(present)
  do ijk = 1, nxyzdm
     r(ijk) = 0.0d+0
  end do
  !$acc end kernels
  
  call ddenst(  r,   tb)
  
#ifdef OPT_TRIPOLE
  call shift1(     r, &
   &             nxdim,  nydim,  nzdim, &
   &              1.d0,      0,      0 )
  call shift1(    tb, &
   &             nxdim,  nydim, nztdim, &
   &              1.d0,      0,      0 )
#else
  call shift1( &
   &                 r, &
   &             nxdim,  nydim,  nzdim)
  call shift1( &
   &                tb, &
   &             nxdim,  nydim, nztdim)
#endif

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
#ifdef OPT_TRIPOLE
  call shift2(    ub,     vb, &
   &                nxdim,  nydim,  nzdim, &
   &                -1.d0,     -1,     -1 )
#else
  call shift2( &
   &                   ub,     vb, &
   &                nxdim,  nydim,  nzdim)
#endif
#endif
#ifdef OPT_TRIPOLE
  call shift2( &
   &                 uadv,  vadv, &
   &                nxdim,  nydim,  nzdim, &
   &                -1.d0,      0,      0)
  call shift3( &
   &                 wadv(1, 1),   wadv(1, 2),   wadv(1, 3), &
   &                nxdim,  nydim,  nzdim, &
   &                 1.d0,     -1,     -1)
  call shift3( &
   &                 wadv(1, 4),   wadv(1, 5),   wadv(1, 6), &
   &                nxdim,  nydim,  nzdim, &
   &                 1.d0,     -1,     -1)
  call shift3( &
   &                 wadv(1, 7),   wadv(1, 8),   wadv(1, 9), &
   &                nxdim,  nydim,  nzdim, &
   &                 1.d0,     -1,     -1)
  call excngw(wadv)

  call shift1(   amv, &
   &             nxdim,  nydim,  nzdim, &
   &              1.d0,     -1,     -1 )
  call shift1(   ahv, &
   &             nxdim,  nydim,  nzdim, &
   &              1.d0,      0,      0 )
#else
  call shift2( &
   &              uadv,   vadv, &
   &             nxdim,  nydim,  nzdim)
  call shift3( &
   &              wadv(1, 1),   wadv(1, 2),   wadv(1, 3), &
   &             nxdim,  nydim,  nzdim)
  call shift3( &
   &              wadv(1, 4),   wadv(1, 5),   wadv(1, 6), &
   &             nxdim,  nydim,  nzdim)
  call shift3( &
   &              wadv(1, 7),   wadv(1, 8),   wadv(1, 9), &
   &             nxdim,  nydim,  nzdim)
  call shift2( &
   &               amv,    ahv, &
   &             nxdim,  nydim,  nzdim)
#endif
  return
end subroutine iniset
end module binst
