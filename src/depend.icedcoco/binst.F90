module binst

! --- information -----------------------------------------------------
!
!  HISTORY
!     '02.10.10  H.Hasumi: from COCO3.4
!     '07.04.23  H.Hasumi
!     '09.02.23  Y.Komuro: call DDENST(diagnose R) instead of OVTURN
!                          bug fix (SHIFT UX/VX after STBBVT)
!     '10.04.14  M.kurogi
!     '10.04.14  M.Kurogi: (COCO4.4 tripolar code by Dr. Suzuki)
!     '12.07.06  M.kurogi: for COCO5.0
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
    use zocmsk, only : amskv
    use tovtr
    use dvlva
    use bstbc
    use bchmk
    use bshft

    implicit none

    real(8), intent(out)   :: uadv(nxyzdm), vadv(nxyzdm), wadv(nxyzdm), r(nxyzdm)
    real(8), intent(inout) ::   ub(nxyzdm),   vb(nxyzdm),   tb(nxyzdm, ntdim)
    real(8), intent(inout) ::   hb(nxydim), ubtb(nxydim),  vbtb(nxydim)
    real(8), intent(inout) ::    w(nxyzdm),  amv(nxyzdm),   ahv(nxyzdm)

    real(8) ::      uz(nxydim),     vz(nxydim)
    integer ::   ijk, ij, k

    call shift1(hb,         nxdim, nydim,      1,  1.d0,  0,  0 )
    call shift2(ubtb, vbtb, nxdim, nydim,      1, -1.d0, -1, -1 )
    call shift2(  ub,   vb, nxdim, nydim,  nzdim, -1.d0, -1, -1 )
    call shift1(   w,       nxdim, nydim,  nzdim,  1.d0,  0,  0 )

    do ijk = 1, nxyzdm
       r(ijk) = 0.0d+0
    end do

    call ddenst(                           &
         &                      r,         &
         &                      tb)

#ifdef OPT_BBL
    call stbbtr(                           &
         &                  tb)
#endif

    call shift1( r, nxdim, nydim, nzdim,  1.d0,  0,  0 )
    call shift1(tb, nxdim, nydim, nztdim, 1.d0,  0,  0 )

    call stbctr(                           &
         &                  tb,      r)

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

    call rmmskt
    call admktb
    call admskv
#endif
    call velvad(                            &
         &           uadv,   vadv,   wadv,  &
         &             ub,     vb,      w)

#ifdef OPT_BBL
    call stbbvt(ub, vb)

    call shift2(ub, vb, nxdim, nydim, nzdim, -1.d0, -1, -1)
    call rmmskt
    call admkt1
    call rmmskv
    call admkv1
#endif


    call shift1(uadv, nxdim, nydim, nzdim, -1.d0,  0,  -1 )
    call shift1(vadv, nxdim, nydim, nzdim, -1.d0, -1,   0 )
    call shift1(wadv, nxdim, nydim, nzdim,  1.d0, -1,  -1 )
    call shift1( amv, nxdim, nydim, nzdim,  1.d0, -1,  -1 )
    call shift1( ahv, nxdim, nydim, nzdim,  1.d0,  0,   0 )
    return
  end subroutine iniset
end module binst

