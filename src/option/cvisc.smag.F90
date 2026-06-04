
module cvisc

! --- information -----------------------------------------------------
!
!  Viscosity term of the equation of motion.
!
!  HISTORY
!     '03.04.21  H.Hasumi: from COCO3.4
!     '07.04.23  H.Hasumi: Formulation changed
!     '07.04.16  H.Hasumi
!     '08.06.11  H.Hasumi: initial/final processing
!     '08.07.10  H.Hasumi: initial/final processing
!     '09.01.06  T.Suzuki: bug fix (FUZ/FVZ added to GX/GY in VSCVLB)
!     '12.07.04  H.Tatebe: for COCO5.0 in F90
! ---------------------------------------------------------------------

  use zocdim,  only :  nxydim,  nzdim

  implicit none
#include "coco.h"
  private

  real(8),     save  ::    sxx(nxydim),          syy(nxydim)
  real(8),     save  ::    sxy(nxydim),          syx(nxydim)
  real(8),     save  ::    bxx(nxydim),          byy(nxydim)
  real(8),     save  ::    bxy(nxydim),          byx(nxydim)
  real(8),     save  ::  hvbot(nxydim)
  logical,     save  ::  ofirst, ofirst_bbl       
  character(len=64)  ::  chead(1:16)
  data ofirst, ofirst_bbl / .true., .true. /

!---- temporary arrays
  real(8),     allocatable,    dimension(:,:)     ::  buf2,    g2d

  public  ::  vscvel
#ifdef OPT_BBL
  public  ::  vscvlb
#endif

contains

  subroutine vscvel(                                                  &
         &      gx,     gy,     xx,     yy,                           &
         &      uy,     ux,     vy,     vx,                           &
         &      hx,     hy,                                           &
         &     amv,   taux,   tauy        )
    use zocdim,  only :                                               &
         &  nxydim,  nxdim,  nzdim,                                   &
         &    kstr,   kend,     kz,                                   &
         &   ijstr,  ijend, ijvstr, ijvend,                           &
         &      le,     lw,     ln,     ls,                           &
         &     lne,    lsw,    lnw,    lse,                           &
         &   oinit,  ofinal
    use zocgrd,  only :                                               &
         &    zbot,     dx,                                           &  
         &      dy,    dzv,    dzm,    rea,                           &  
         &      rx,     ry,    rym,     rs,    rsm,                   & 
         &     hxt,    hyt,   hxyt,   hyxt,                           &
         &     rxu,    ryu,    rxt,    ryt   
    use zocmsk,  only :  amskv,  amfvx,  amfvy
    use ufile
    
    implicit none

    real(8),   intent(inout)  ::     gx(nxydim,nzdim),    gy(nxydim,nzdim)
    real(8),   intent(inout)  ::     xx(nxydim,nzdim),    yy(nxydim,nzdim)
    real(8),   intent(in)     ::     uy(nxydim,nzdim),    ux(nxydim,nzdim)
    real(8),   intent(in)     ::     vy(nxydim,nzdim),    vx(nxydim,nzdim)
    real(8),   intent(in)     ::     hx(nxydim),          hy(nxydim)
    real(8),   intent(in)     ::    amv(nxydim,nzdim)
    real(8),   intent(in)     ::   taux(nxydim),        tauy(nxydim)

!---- local variables
    real(8)            ::    fux(nxydim,nzdim),    fvx(nxydim,nzdim)
    real(8)            ::    fuy(nxydim,nzdim),    fvy(nxydim,nzdim)
    real(8)            ::    fuz(nxydim,nzdim),    fvz(nxydim,nzdim)
    real(8),     save  ::     rz(nxydim,nzdim),    rzm(nxydim,nzdim)
    real(8),     save  ::    szx(nxydim,nzdim),    szy(nxydim,nzdim)
    real(8)            ::     db(nxydim)
    real(8)            ::     hu(nxydim),          hv(nxydim)
    real(8)            ::    hux(nxydim),         hvx(nxydim)
    real(8)            ::    huy(nxydim),         hvy(nxydim)
    
    real(8)            ::    exx,    eyy,    exy,     ez,    dd
    integer(4)         ::     ij,      k,      i,      j
    integer(4)         ::   ijln,   ijls,   ijle,   ijlw
    integer(4)         ::   ijnw,   ijse,   ijsw,   ijne
    integer(4)         ::  ifpar,  jfpar,  istat

    real(8),     save  ::    asm(nxydim), asmb(nxydim)
    real(8),     save  ::    csm=0.d0, csmb=0.d0
    real(8)            ::     pi
    namelist /nmsmag/ csm, csmb


    real(8)            ::   wte,   wtw,   wtn,   wts
    real(8)            ::  wtne,  wtnw,  wtse,  wtsw
    real(8)            :: sxxne(nxydim), sxxnw(nxydim)
    real(8)            :: sxxse(nxydim), sxxsw(nxydim)
    real(8)            :: sxyne(nxydim), sxynw(nxydim)
    real(8)            :: sxyse(nxydim), sxysw(nxydim)
    real(8)            :: syyne(nxydim), syynw(nxydim)
    real(8)            :: syyse(nxydim), syysw(nxydim)
    real(8)            ::  dbsw(nxydim),  dbnw(nxydim)
    real(8)            ::  dbse(nxydim),  dbne(nxydim)
    real(8)            ::  szx0,   szy0
    real(8)            :: bxxsw,  bxxnw,  bxxse,  bxxne
    real(8)            :: byysw,  byynw,  byyse,  byyne
    real(8)            :: bxysw,  bxynw,  bxyse,  bxyne
    real(8)            ::  huxw,   huxe,   huys,   huyn
    real(8)            ::  hvxw,   hvxe,   hvys,   hvyn
    logical,     save  ::    opslvis = .false.
    namelist /nmpslv/ opslvis

    if ( oinit .or. ofinal ) then
       return
    end if

    if ( ofirst ) then
       ofirst = .false.
       
       !$acc enter data create(asm, asmb)
       !$acc enter data create(sxx,syy)
       !$acc enter data create(sxy,syx)
       !$acc enter data create(bxx,byy)
       !$acc enter data create(bxy,byx)
       !$acc enter data create(hvbot)

       !$acc enter data create(fux,fvx)
       !$acc enter data create(fuy,fvy)
       !$acc enter data create(fuz,fvz)
       !$acc enter data create(rz, rzm)
       !$acc enter data create(szx,szy)
       !$acc enter data create(db)
       !$acc enter data create(hu,hv)
       !$acc enter data create(hux,hvx)
       !$acc enter data create(huy,hvy)

       !$acc enter data create( sxxne, sxxnw )
       !$acc enter data create( sxxse, sxxsw )
       !$acc enter data create( sxyne, sxynw )
       !$acc enter data create( sxyse, sxysw )
       !$acc enter data create( syyne, syynw )
       !$acc enter data create( syyse, syysw )
       !$acc enter data create(  dbsw,  dbnw )
       !$acc enter data create(  dbse,  dbne )
    
       call rewnml( ifpar, jfpar )
       READ_NAMELIST( nmsmag )
       READ_NAMELIST( nmpslv )

       pi = 4.d0 * atan(1.d0)
       !$acc kernels default(present)
       do ij = 1, nxydim
          asm(ij) = ( csm * min( dx*hxt(ij), dy(ij)*hyt(ij)) / pi ) ** 2
          asmb(ij) = csmb * min( dx*hxt(ij), dy(ij)*hyt(ij))**2        &
   &              / pi / sqrt(8.d0)
       end do
       !$acc end kernels
    end if

    !$acc kernels default(present)
    do k = 1, nzdim
       do ij = 1, nxydim
          fux(ij, k) = 0.d0
          fvx(ij, k) = 0.d0
          fuy(ij, k) = 0.d0
          fvy(ij, k) = 0.d0
          fuz(ij, k) = 0.d0
          fvz(ij, k) = 0.d0
       end do
    end do
    
    do ij = ijvstr, ijvend
       fuz(ij, kstr) = taux(ij)
       fvz(ij, kstr) = tauy(ij)
    end do

    do ij = ijvstr, ijvend
       hvbot(ij) = (  (hx(ij)    + hx(ij+le) ) * dy(ij)               &
    &               + (hx(ij+ln) + hx(ij+lne)) * dy(ij+ln) ) *        &
    &               rym(ij) * 0.25d0                                  &
    &             + zbot
       hvbot(ij) = 1.d0 / hvbot(ij)
    end do
    do k = kstr, kstr+kz-1
       do ij = ijvstr, ijvend
          rz (ij, k) = 1.d0 * rs (k) * hvbot(ij)
          rzm(ij, k) = 1.d0 * rsm(k) * hvbot(ij)
       end do
    end do
    do k = kstr+kz, kend
       do ij = ijvstr, ijvend
          rz (ij, k) = 1.d0 / dzv(ij, k)
          rzm(ij, k) = 1.d0 / dzm(ij, k)
       end do
    end do
    
    do k = kstr+1, kend
       do ij = ijvstr, ijvend
          fuz(ij, k) = amv(ij, k) * rzm(ij, k) * ( ux(ij, k-1) - ux(ij, k) )
          fvz(ij, k) = amv(ij, k) * rzm(ij, k) * ( vx(ij, k-1) - vx(ij, k) )
       end do
    end do
    
    do k = kstr, kend
       do ij = ijvstr, ijvend
          gx(ij, k) = gx(ij, k) + (fuz(ij, k) - fuz(ij, k+1)) * rz(ij, k)
          gy(ij, k) = gy(ij, k) + (fvz(ij, k) - fvz(ij, k+1)) * rz(ij, k)
       end do
    end do
    !$acc end kernels
    
    if (opslvis) then
       !$acc kernels default(present)
       do k = kstr, kend
          sxxsw(:) = 0.d0
          sxxnw(:) = 0.d0
          sxxse(:) = 0.d0
          sxxne(:) = 0.d0
          syysw(:) = 0.d0
          syynw(:) = 0.d0
          syyse(:) = 0.d0
          syyne(:) = 0.d0
          sxysw(:) = 0.d0
          sxynw(:) = 0.d0
          sxyse(:) = 0.d0
          sxyne(:) = 0.d0

          do ij = ijvstr-nxdim-1, ijvend+nxdim+1
             ijls = ij + ls
             ijlw = ij + lw
             ijln = ij + ln
             ijle = ij + le
             ijnw = ij + lnw
             ijse = ij + lse
             ijsw = ij + lsw
             ijne = ij + lne

             if (dzv(ij, k) .ge. dzv(ijlw, k)) then
                wtw = dzv(ijlw, k) / dzv(ij, k)
             else
                wtw = 2.d0 - dzv(ij, k) / dzv(ijlw, k)
             end if
             if (dzv(ij, k) .ge. dzv(ijls, k)) then
                wts = dzv(ijls, k) / dzv(ij, k)
             else
                wts = 2.d0 - dzv(ij, k) / dzv(ijls, k)
             end if
             if (dzv(ij, k) .ge. dzv(ijln, k)) then
                wtn = dzv(ijln, k) / dzv(ij, k)
             else
                wtn = 2.d0 - dzv(ij, k) / dzv(ijln, k)
             end if
             if (dzv(ij, k) .ge. dzv(ijle, k)) then
                wte = dzv(ijle, k) / dzv(ij, k)
             else
                wte = 2.d0 - dzv(ij, k) / dzv(ijle, k)
             end if
             if (dzv(ij, k) .ge. dzv(ijsw, k)) then
                wtsw = dzv(ijsw, k) / dzv(ij, k)
             else
                wtsw = 2.d0 - dzv(ij, k) / dzv(ijsw, k)
             end if
             if (dzv(ij, k) .ge. dzv(ijnw, k)) then
                wtnw = dzv(ijnw, k) / dzv(ij, k)
             else
                wtnw = 2.d0 - dzv(ij, k) / dzv(ijnw, k)
             end if
             if (dzv(ij, k) .ge. dzv(ijse, k)) then
                wtse = dzv(ijse, k) / dzv(ij, k)
             else
                wtse = 2.d0 - dzv(ij, k) / dzv(ijse, k)
             end if
             if (dzv(ij, k) .ge. dzv(ijne, k)) then
                wtne = dzv(ijne, k) / dzv(ij, k)
             else
                wtne = 2.d0 - dzv(ij, k) / dzv(ijne, k)
             end if

             exx =    (ux(ij,   k)       + ux(ijls, k) * wts                               &
                  & -  ux(ijlw, k) * wtw - ux(ijsw, k) * wtsw) *  0.5d0 * rx * rxt(ij)     &
                  & + (vx(ij,   k)       + vx(ijls, k) * wts                               &
                  & +  vx(ijlw, k) * wtw + vx(ijsw, k) * wtsw) * 0.25d0 * hxyt(ij)
             eyy =    (vx(ij,   k)       + vx(ijlw, k) * wtw                               &
                  & -  vx(ijls, k) * wts - vx(ijsw, k) * wtsw) *  0.5d0 * ry(ij) * ryt(ij) &
                  & + (ux(ij,   k)       + ux(ijlw, k) * wtw                               &
                  & +  ux(ijls, k) * wts + ux(ijsw, k) * wtsw) * 0.25d0 * hyxt(ij)
             exy =    (ux(ij,   k)       + ux(ijlw, k) * wtw                               &
                  & -  ux(ijls, k) * wts - ux(ijsw, k) * wtsw) *  0.5d0 * ry(ij) * ryt(ij) &
                  & + (vx(ij,   k)       + vx(ijls, k) * wts                               &
                  & -  vx(ijlw, k) * wtw - vx(ijsw, k) * wtsw) *  0.5d0 * rx * rxt(ij)     &
                  & - (ux(ij,   k)       + ux(ijlw, k) * wtw                               &
                  & +  ux(ijls, k) * wts + ux(ijsw, k) * wtsw) * 0.25d0 * hxyt(ij)         &
                  & - (vx(ij,   k)       + vx(ijls, k) * wts                               &
                  & +  vx(ijlw, k) * wtw + vx(ijsw, k) * wtsw) * 0.25d0 * hyxt(ij)
             dd = sqrt((exx - eyy) * (exx - eyy) + exy * exy)
             sxxsw(ij) = asm(ij) * dd * (exx - eyy)
             syysw(ij) = asm(ij) * dd * (eyy - exx)
             sxysw(ij) = asm(ij) * dd * exy

             dbsw(ij) = sqrt(sqrt((exx - eyy) * (exx - eyy) + exy * exy))
             bxxsw = asmb(ij) * dbsw(ij) * (exx - eyy)
             byysw = asmb(ij) * dbsw(ij) * (eyy - exx)
             bxysw = asmb(ij) * dbsw(ij) * exy

             exx =    (ux(ij,   k)       + ux(ijln, k) * wtn                                   &
                  & -  ux(ijlw, k) * wtw - ux(ijnw, k) * wtnw) *  0.5d0 * rx * rxt(ijln)       &
                  & + (vx(ij,   k)       + vx(ijln, k) * wtn                                   &
                  & +  vx(ijlw, k) * wtw + vx(ijnw, k) * wtnw) * 0.25d0 * hxyt(ijln)
             eyy =    (vx(ijnw, k) * wtnw+ vx(ijln, k) * wtn                                   &
                  & -  vx(ijlw, k) * wtw - vx(ij,   k)       ) *  0.5d0 * ry(ijln) * ryt(ijln) &
                  & + (ux(ijnw, k) * wtnw+ ux(ijln, k) * wtn                                   &
                  & +  ux(ijlw, k) * wtw + ux(ij,   k)       ) * 0.25d0 * hyxt(ijln)
             exy = (   ux(ijnw, k) * wtnw+ ux(ijln, k) * wtn                                   &
                  & -  ux(ijlw, k) * wtw - ux(ij  , k)       ) *  0.5d0 * ry(ijln) * ryt(ijln) &
                  & + (vx(ij,   k)       + vx(ijln, k) * wtn                                   &
                  & -  vx(ijlw, k) * wtw - vx(ijnw, k) * wtnw) *  0.5d0 * rx * rxt(ijln)       &
                  & - (ux(ijnw, k) * wtnw+ ux(ijln, k) * wtn                                   &
                  & +  ux(ijlw, k) * wtw + ux(ij  , k)       ) * 0.25d0 * hxyt(ijln)           &
                  & - (vx(ij,   k)       + vx(ijln, k) * wtn                                   &
                  & +  vx(ijlw, k) * wtw + vx(ijnw, k) * wtnw) * 0.25d0 * hyxt(ijln)
             dd = sqrt((exx - eyy) * (exx - eyy) + exy * exy)
             sxxnw(ij) = asm(ijln) * dd * (exx - eyy)
             syynw(ij) = asm(ijln) * dd * (eyy - exx)
             sxynw(ij) = asm(ijln) * dd * exy

             dbnw(ij) = sqrt(sqrt((exx - eyy) * (exx - eyy) + exy * exy))
             bxxnw = asmb(ijln) * dbnw(ij) * (exx - eyy)
             byynw = asmb(ijln) * dbnw(ij) * (eyy - exx)
             bxynw = asmb(ijln) * dbnw(ij) * exy
             
             exx =    (ux(ijse, k) * wtse+ ux(ijle, k) * wte                                   &
                  & -  ux(ijls, k) * wts - ux(ij,   k)       ) *  0.5d0 * rx * rxt(ijle)       &
                  & + (vx(ijse, k) * wtse+ vx(ijle, k) * wte                                   &
                  & +  vx(ijls, k) * wts + vx(ij,   k)       ) * 0.25d0 * hxyt(ijle)
             eyy =    (vx(ij,   k)       + vx(ijle, k) * wte                                   &
                  & -  vx(ijls, k) * wts - vx(ijse, k) * wtse) *  0.5d0 * ry(ijle) * ryt(ijle) &
                  & + (ux(ij,   k)       + ux(ijle, k) * wte                                   &
                  & +  ux(ijls, k) * wts + ux(ijse, k) * wtse) * 0.25d0 * hyxt(ijle)
             exy =    (ux(ij,   k)       + ux(ijle, k) * wte                                   &
                  & -  ux(ijls, k) * wts - ux(ijse, k) * wtse) *  0.5d0 * ry(ijle) * ryt(ijle) &
                  & + (vx(ijse, k) * wtse+ vx(ijle, k) * wte                                   &
                  & -  vx(ijls, k) * wts - vx(ij,   k)       ) *  0.5d0 * rx * rxt(ijle)       &
                  & - (ux(ij,   k)       + ux(ijle, k) * wte                                   &
                  & +  ux(ijls, k) * wts + ux(ijse, k) * wtse) * 0.25d0 * hxyt(ijle)           &
                  & - (vx(ijse, k) * wtse+ vx(ijle, k) * wte                                   &
                  & +  vx(ijls, k) * wts + vx(ij,   k)       ) * 0.25d0 * hyxt(ijle)
             dd = sqrt((exx - eyy) * (exx - eyy) + exy * exy)
             sxxse(ij) = asm(ijle) * dd * (exx - eyy)
             syyse(ij) = asm(ijle) * dd * (eyy - exx)
             sxyse(ij) = asm(ijle) * dd * exy

             dbse(ij) = sqrt(sqrt((exx - eyy) * (exx - eyy) + exy * exy))
             bxxse = asmb(ijle) * dbse(ij) * (exx - eyy)
             byyse = asmb(ijle) * dbse(ij) * (eyy - exx)
             bxyse = asmb(ijle) * dbse(ij) * exy
             
             exx =    (ux(ijne, k) * wtne+ ux(ijle, k) * wte                                   &
                  & -  ux(ijln, k) * wtn - ux(ij,   k)       ) *  0.5d0 * rx * rxt(ijne)       &
                  & + (vx(ijne, k) * wtne+ vx(ijle, k) * wte                                   &
                  & +  vx(ijln, k) * wtn + vx(ij,   k)       ) * 0.25d0 * hxyt(ijne)
             eyy =    (vx(ijne, k) * wtne+ vx(ijln, k) * wtn                                   &
                  & -  vx(ijle, k) * wte - vx(ij,   k)       ) *  0.5d0 * ry(ijne) * ryt(ijne) &
                  & + (ux(ijne, k) * wtne+ ux(ijln, k) * wtn                                   &
                  & +  ux(ijle, k) * wte + ux(ij,   k)       ) * 0.25d0 * hyxt(ijne)
             exy =    (ux(ijne, k) * wtne+ ux(ijln, k) * wtn                                   &
                  & -  ux(ijle, k) * wte - ux(ij,   k)       ) *  0.5d0 * ry(ijne) * ryt(ijne) &
                  & + (vx(ijne, k) * wtne+ vx(ijle, k) * wte                                   &
                  & -  vx(ijln, k) * wtn - vx(ij,   k)       ) *  0.5d0 * rx * rxt(ijne)       &
                  & - (ux(ijne, k) * wtne+ ux(ijln, k) * wtn                                   &
                  & +  ux(ijle, k) * wte + ux(ij,   k)       ) * 0.25d0 * hxyt(ijne)           &
                  & - (vx(ijne, k)       + vx(ijle, k) * wte                                   &
                  & +  vx(ijln, k) * wtn + vx(ij,   k)       ) * 0.25d0 * hyxt(ijne)
             dd = sqrt((exx - eyy) * (exx - eyy) + exy * exy)
             sxxne(ij) = asm(ijne) * dd * (exx - eyy)
             syyne(ij) = asm(ijne) * dd * (eyy - exx)
             sxyne(ij) = asm(ijne) * dd * exy

             dbne(ij) = sqrt(sqrt((exx - eyy) * (exx - eyy) + exy * exy))
             bxxne = asmb(ijne) * dbne(ij) * (exx - eyy)
             byyne = asmb(ijne) * dbne(ij) * (eyy - exx)
             bxyne = asmb(ijne) * dbne(ij) * exy

             huxw = ( bxxsw * hyt(ij)   * hyt(ij)   &
                  & + bxxnw * hyt(ijln) * hyt(ijln) ) * 0.5d0
             huxe = ( bxxse * hyt(ijle) * hyt(ijle) &
                  & + bxxne * hyt(ijne) * hyt(ijne) ) * 0.5d0
             hvxw = ( bxysw * hyt(ij)   * hyt(ij)   &
                  & + bxynw * hyt(ijln) * hyt(ijln) ) * 0.5d0
             hvxe = ( bxyse * hyt(ijle) * hyt(ijle) &
                  & + bxyne * hyt(ijne) * hyt(ijne) ) * 0.5d0
             huys = ( bxysw * hxt(ij)   * hxt(ij)   &
                  & + bxyse * hxt(ijle) * hxt(ijle) ) * 0.5d0
             huyn = ( bxynw * hxt(ijln) * hxt(ijln) &
                  & + bxyne * hxt(ijne) * hxt(ijne) ) * 0.5d0
             hvys = ( byysw * hxt(ij)   * hxt(ij)   &
                  & + byyse * hxt(ijle) * hxt(ijle) ) * 0.5d0
             hvyn = ( byynw * hxt(ijln) * hxt(ijln) &
                  & + byyne * hxt(ijne) * hxt(ijne) ) * 0.5d0

             hu(ij) = - ((huxe - huxw) * rx      * ryu(ij)  &
                  &   +  (huyn - huys) * rym(ij) * rxu(ij)) &
                  &                    * rxu(ij) * ryu(ij) * amskv(ij, k)            
             hv(ij) = - ((hvxe - hvxw) * rx      * ryu(ij)  &
                  &   +  (hvyn - hvys) * rym(ij) * rxu(ij)) &
                  &                    * rxu(ij) * ryu(ij) * amskv(ij, k)

          end do

          if (csmb /= 0.d0) then
             do ij = ijvstr, ijvend

                ijls = ij + ls
                ijlw = ij + lw
                ijln = ij + ln
                ijle = ij + le
                ijnw = ij + lnw
                ijse = ij + lse
                ijsw = ij + lsw
                ijne = ij + lne

                if (dzv(ij, k) .ge. dzv(ijlw, k)) then
                   wtw = dzv(ijlw, k) / dzv(ij, k)
                else
                   wtw = 2.d0 - dzv(ij, k) / dzv(ijlw, k)
                end if
                if (dzv(ij, k) .ge. dzv(ijls, k)) then
                   wts = dzv(ijls, k) / dzv(ij, k)
                else
                   wts = 2.d0 - dzv(ij, k) / dzv(ijls, k)
                end if
                if (dzv(ij, k) .ge. dzv(ijln, k)) then
                   wtn = dzv(ijln, k) / dzv(ij, k)
                else
                   wtn = 2.d0 - dzv(ij, k) / dzv(ijln, k)
                end if
                if (dzv(ij, k) .ge. dzv(ijle, k)) then
                   wte = dzv(ijle, k) / dzv(ij, k)
                else
                   wte = 2.d0 - dzv(ij, k) / dzv(ijle, k)
                end if
                if (dzv(ij, k) .ge. dzv(ijsw, k)) then
                   wtsw = dzv(ijsw, k) / dzv(ij, k)
                else
                   wtsw = 2.d0 - dzv(ij, k) / dzv(ijsw, k)
                end if
                if (dzv(ij, k) .ge. dzv(ijnw, k)) then
                   wtnw = dzv(ijnw, k) / dzv(ij, k)
                else
                   wtnw = 2.d0 - dzv(ij, k) / dzv(ijnw, k)
                end if
                if (dzv(ij, k) .ge. dzv(ijse, k)) then
                   wtse = dzv(ijse, k) / dzv(ij, k)
                else
                   wtse = 2.d0 - dzv(ij, k) / dzv(ijse, k)
                end if
                if (dzv(ij, k) .ge. dzv(ijne, k)) then
                   wtne = dzv(ijne, k) / dzv(ij, k)
                else
                   wtne = 2.d0 - dzv(ij, k) / dzv(ijne, k)
                end if

                exx =    (hu(ij)         + hu(ijls) * wts                               &
                     & -  hu(ijlw) * wtw - hu(ijsw) * wtsw) *  0.5d0 * rx * rxt(ij)     &
                     & + (hv(ij)         + hv(ijls) * wts                               &
                     & +  hv(ijlw) * wtw + hv(ijsw) * wtsw) * 0.25d0 * hxyt(ij)
                eyy =    (hv(ij)         + hv(ijlw) * wtw                               &
                     & -  hv(ijls) * wts - hv(ijsw) * wtsw) *  0.5d0 * ry(ij) * ryt(ij) &
                     & + (hu(ij)         + hu(ijlw) * wtw                               &
                     & +  hu(ijls) * wts + hu(ijsw) * wtsw) * 0.25d0 * hyxt(ij)
                exy =    (hu(ij)         + hu(ijlw) * wtw                               &
                     & -  hu(ijls) * wts - hu(ijsw) * wtsw) *  0.5d0 * ry(ij) * ryt(ij) &
                     & + (hv(ij)         + hv(ijls) * wts                               &
                     & -  hv(ijlw) * wtw - hv(ijsw) * wtsw) *  0.5d0 * rx * rxt(ij)     &
                     & - (hu(ij)         + hu(ijlw) * wtw                               &
                     & +  hu(ijls) * wts + hu(ijsw) * wtsw) * 0.25d0 * hxyt(ij)         &
                     & - (hv(ij)         + hv(ijls) * wts                               &
                     & +  hv(ijlw) * wtw + hv(ijsw) * wtsw) * 0.25d0 * hyxt(ij)

                sxxsw(ij) = sxxsw(ij) + asmb(ij) * dbsw(ij) * (exx - eyy)
                syysw(ij) = syysw(ij) + asmb(ij) * dbsw(ij) * (eyy - exx)
                sxysw(ij) = sxysw(ij) + asmb(ij) * dbsw(ij) * exy

                exx =    (hu(ij)         + hu(ijln) * wtn                                   &
                     & -  hu(ijlw) * wtw - hu(ijnw) * wtnw) *  0.5d0 * rx * rxt(ijln)       &
                     & + (hv(ij)         + hv(ijln) * wtn                                   &
                     & +  hv(ijlw) * wtw + hv(ijnw) * wtnw) * 0.25d0 * hxyt(ijln)
                eyy =    (hv(ijnw) * wtnw+ hv(ijln) * wtn                                   &
                     & -  hv(ijlw) * wtw - hv(ij)         ) *  0.5d0 * ry(ijln) * ryt(ijln) &
                     & + (hu(ijnw) * wtnw+ hu(ijln) * wtn                                   &
                     & +  hu(ijlw) * wtw + hu(ij)         ) * 0.25d0 * hyxt(ijln)
                exy = (   hu(ijnw) * wtnw+ hu(ijln) * wtn                                   &
                     & -  hu(ijlw) * wtw - hu(ij)         ) *  0.5d0 * ry(ijln) * ryt(ijln) &
                     & + (hv(ij)         + hv(ijln) * wtn                                   &
                     & -  hv(ijlw) * wtw - hv(ijnw) * wtnw) *  0.5d0 * rx * rxt(ijln)       &
                     & - (hu(ijnw) * wtnw+ hu(ijln) * wtn                                   &
                     & +  hu(ijlw) * wtw + hu(ij)         ) * 0.25d0 * hxyt(ijln)           &
                     & - (hv(ij)         + hv(ijln) * wtn                                   &
                     & +  hv(ijlw) * wtw + hv(ijnw) * wtnw) * 0.25d0 * hyxt(ijln)

                sxxnw(ij) = sxxnw(ij) + asmb(ijln) * dbnw(ij) * (exx - eyy)
                syynw(ij) = syynw(ij) + asmb(ijln) * dbnw(ij) * (eyy - exx)
                sxynw(ij) = sxynw(ij) + asmb(ijln) * dbnw(ij) * exy

                exx =    (hu(ijse) * wtse+ hu(ijle) * wte                                   &
                     & -  hu(ijls) * wts - hu(ij)         ) *  0.5d0 * rx * rxt(ijle)       &
                     & + (hv(ijse) * wtse+ hv(ijle) * wte                                   &
                     & +  hv(ijls) * wts + hv(ij)         ) * 0.25d0 * hxyt(ijle)
                eyy =    (hv(ij)         + hv(ijle) * wte                                   &
                     & -  hv(ijls) * wts - hv(ijse) * wtse) *  0.5d0 * ry(ijle) * ryt(ijle) &
                     & + (hu(ij)         + hu(ijle) * wte                                   &
                     & +  hu(ijls) * wts + hu(ijse) * wtse) * 0.25d0 * hyxt(ijle)
                exy =    (hu(ij)         + hu(ijle) * wte                                   &
                     & -  hu(ijls) * wts - hu(ijse) * wtse) *  0.5d0 * ry(ijle) * ryt(ijle) &
                     & + (hv(ijse) * wtse+ hv(ijle) * wte                                   &
                     & -  hv(ijls) * wts - hv(ij)         ) *  0.5d0 * rx * rxt(ijle)       &
                     & - (hu(ij)         + hu(ijle) * wte                                   &
                     & +  hu(ijls) * wts + hu(ijse) * wtse) * 0.25d0 * hxyt(ijle)           &
                     & - (hv(ijse) * wtse+ hv(ijle) * wte                                   &
                     & +  hv(ijls) * wts + hv(ij)         ) * 0.25d0 * hyxt(ijle)

                sxxse(ij) = sxxse(ij) + asmb(ijle) * dbse(ij) * (exx - eyy)
                syyse(ij) = syyse(ij) + asmb(ijle) * dbse(ij) * (eyy - exx)
                sxyse(ij) = sxyse(ij) + asmb(ijle) * dbse(ij) * exy

                exx =    (hu(ijne) * wtne+ hu(ijle) * wte                                   &
                     & -  hu(ijln) * wtn - hu(ij)         ) *  0.5d0 * rx * rxt(ijne)       &
                     & + (hv(ijne) * wtne+ hv(ijle) * wte                                   &
                     & +  hv(ijln) * wtn + hv(ij)         ) * 0.25d0 * hxyt(ijne)
                eyy =    (hv(ijne) * wtne+ hv(ijln) * wtn                                   &
                     & -  hv(ijle) * wte - hv(ij)         ) *  0.5d0 * ry(ijne) * ryt(ijne) &
                     & + (hu(ijne) * wtne+ hu(ijln) * wtn                                   &
                     & +  hu(ijle) * wte + hu(ij)         ) * 0.25d0 * hyxt(ijne)
                exy =    (hu(ijne) * wtne+ hu(ijln) * wtn                                   &
                     & -  hu(ijle) * wte - hu(ij)         ) *  0.5d0 * ry(ijne) * ryt(ijne) &
                     & + (hv(ijne) * wtne+ hv(ijle) * wte                                   &
                     & -  hv(ijln) * wtn - hv(ij)         ) *  0.5d0 * rx * rxt(ijne)       &
                     & - (hu(ijne) * wtne+ hu(ijln) * wtn                                   &
                     & +  hu(ijle) * wte + hu(ij)         ) * 0.25d0 * hxyt(ijne)           &
                     & - (hv(ijne)       + hv(ijle) * wte                                   &
                     & +  hv(ijln) * wtn + hv(ij)         ) * 0.25d0 * hyxt(ijne)

                sxxne(ij) = sxxne(ij) + asmb(ijne) * dbne(ij) * (exx - eyy)
                syyne(ij) = syyne(ij) + asmb(ijne) * dbne(ij) * (eyy - exx)
                sxyne(ij) = sxyne(ij) + asmb(ijne) * dbne(ij) * exy

             end do
          end if
          
          do ij = ijvstr, ijvend
             
             ijln = ij + ln
             ijle = ij + le
             ijne = ij + lne
             ez =   ( (amv(ij, k) + amv(ij, k+1)) * 0.5d0 / rea &
                  & + (amv(ij, k) - amv(ij, k+1)) * rz(ij, k) )
             szx0 = - ux(ij, k) * ez
             szy0 = - vx(ij, k) * ez
             
             gx(ij, k) = (gx(ij, k)                                                    &
                  & + ((sxxne(ij) * hyt(ijne) * hyt(ijne)                              &
                  & +   sxxse(ij) * hyt(ijle) * hyt(ijle)                              &
                  & -   sxxnw(ij) * hyt(ijln) * hyt(ijln)                              &
                  & -   sxxsw(ij) * hyt(ij)   * hyt(ij) ) * 0.5d0 * rx * ryu(ij)       &
                  & +  (sxyne(ij) * hxt(ijne) * hxt(ijne)                              &
                  & +   sxynw(ij) * hxt(ijln) * hxt(ijln)                              &
                  & -   sxyse(ij) * hxt(ijle) * hxt(ijle)                              &
                  & -   sxysw(ij) * hxt(ij)   * hxt(ij) ) * 0.5d0 * rym(ij) * rxu(ij)) &
                  & * rxu(ij) * ryu(ij) + szx0 / rea) * amskv(ij, k)

             gy(ij, k) = (gy(ij, k)                                                    &
                  & + ((sxyne(ij) * hyt(ijne) * hyt(ijne)                              &
                  & +   sxyse(ij) * hyt(ijle) * hyt(ijle)                              &
                  & -   sxynw(ij) * hyt(ijln) * hyt(ijln)                              &
                  & -   sxysw(ij) * hyt(ij)   * hyt(ij) ) * 0.5d0 * rx * ryu(ij)       &
                  & +  (syyne(ij) * hxt(ijne) * hxt(ijne)                              &
                  & +   syynw(ij) * hxt(ijln) * hxt(ijln)                              &
                  & -   syyse(ij) * hxt(ijle) * hxt(ijle)                              &
                  & -   syysw(ij) * hxt(ij)   * hxt(ij) ) * 0.5d0 * rym(ij) * rxu(ij)) &
                  & * rxu(ij) * ryu(ij) + szy0 / rea) * amskv(ij, k)
          end do
       end do
       !$acc end kernels
    else
       !$acc kernels default(present)
       do k = kstr, kend
          do ij = 1, nxydim
             sxx(ij) = 0.d0
             syy(ij) = 0.d0
             sxy(ij) = 0.d0
          end do

          do ij = ijstr-nxdim-1, ijend+nxdim+nxdim+2

             ijls  = ij + ls
             ijlw  = ij + lw
             ijsw = ij + lsw
             exx =    (ux(ij,   k) + ux(ijls, k)                              &
                  & -  ux(ijlw, k) - ux(ijsw, k)) *  0.5d0 * rx * rxt(ij)     &
                  & + (vx(ij,   k) + vx(ijls, k)                              &
                  & +  vx(ijlw, k) + vx(ijsw, k)) * 0.25d0 * hxyt(ij)
             eyy =    (vx(ij,   k) + vx(ijlw, k)                              & 
                  & -  vx(ijls, k) - vx(ijsw, k)) *  0.5d0 * ry(ij) * ryt(ij) &
                  & + (ux(ij,   k) + ux(ijlw, k)                              &
                  & +  ux(ijls, k) + ux(ijsw, k)) * 0.25d0 * hyxt(ij)
             exy =    (ux(ij,   k) + ux(ijlw, k)                              &
                  & -  ux(ijls, k) - ux(ijsw, k)) *  0.5d0 * ry(ij) * ryt(ij) &
                  & + (vx(ij,   k) + vx(ijls, k)                              &
                  & -  vx(ijlw, k) - vx(ijsw, k)) *  0.5d0 * rx * rxt(ij)     &
                  & - (ux(ij,   k) + ux(ijlw, k)                              &
                  & +  ux(ijls, k) + ux(ijsw, k)) * 0.25d0 * hxyt(ij)         &
                  & - (vx(ij,   k) + vx(ijlw, k)                              &
                  & +  vx(ijls, k) + vx(ijsw, k)) * 0.25d0 * hyxt(ij)

             dd = sqrt((exx - eyy) * (exx - eyy) + exy * exy)
             sxx(ij) = asm(ij) * dd * (exx - eyy)
             syy(ij) = asm(ij) * dd * (eyy - exx)
             sxy(ij) = asm(ij) * dd * exy

             db(ij) = sqrt(sqrt((exx - eyy) * (exx - eyy) + exy * exy))

             bxx(ij) = asmb(ij) * db(ij) * (exx - eyy)
             byy(ij) = asmb(ij) * db(ij) * (eyy - exx)
             bxy(ij) = asmb(ij) * db(ij) * exy

          end do

          do ij = ijstr-nxdim-1, ijend+nxdim+2
             ijln = ij + ln
             hux(ij) = (bxx(ij)   * hyt(ij)   * hyt(ij) &
                  &  +  bxx(ijln) * hyt(ijln) * hyt(ijln)) * 0.5d0
             hvx(ij) = (bxy(ij)   * hyt(ij)   * hyt(ij) &
                  &  +  bxy(ijln) * hyt(ijln) * hyt(ijln)) * 0.5d0
          end do

          do ij = ijstr-nxdim-1, ijend+nxdim+nxdim+1
             ijle = ij + le
             huy(ij) = (bxy(ij)   * hxt(ij)   * hxt(ij) &
                  &  +  bxy(ijle) * hxt(ijle) * hxt(ijle)) * 0.5d0
             hvy(ij) = (byy(ij)   * hxt(ij)   * hxt(ij) &
                  &  +  byy(ijle) * hxt(ijle) * hxt(ijle)) * 0.5d0
          end do

          do ij = ijstr-nxdim-1, ijend+nxdim+1
             hu(ij) = - ((hux(ij+le) - hux(ij)) * rx      * ryu(ij)  &
                  &   +  (huy(ij+ln) - huy(ij)) * rym(ij) * rxu(ij)) &
                  &                             * rxu(ij) * ryu(ij) * amskv(ij, k)            
             hv(ij) = - ((hvx(ij+le) - hvx(ij)) * rx      * ryu(ij)  &
                  &   +  (hvy(ij+ln) - hvy(ij)) * rym(ij) * rxu(ij)) &
                  &                             * rxu(ij) * ryu(ij) * amskv(ij, k)
          end do

          do ij = ijstr, ijend+nxdim+1
             ijls = ij + ls
             ijlw = ij + lw
             ijsw = ij + lsw
             exx =    (hu(ij)   + hu(ijls)                              &
                  & -  hu(ijlw) - hu(ijsw)) *  0.5d0 * rx * rxt(ij)     &
                  & + (hv(ij)   + hv(ijls)                              &
                  & +  hv(ijlw) + hv(ijsw)) * 0.25d0 * hxyt(ij)                                            
             eyy =    (hv(ij)   + hv(ijlw)                              & 
                  & -  hv(ijls) - hv(ijsw)) *  0.5d0 * ry(ij) * ryt(ij) &
                  & + (hu(ij)   + hu(ijlw)                              &
                  & +  hu(ijls) + hu(ijsw)) * 0.25d0 * hyxt(ij)                                            
             exy =    (hu(ij)   + hu(ijlw )                             &
                  & -  hu(ijls) - hu(ijsw)) *  0.5d0 * ry(ij) * ryt(ij) &
                  & + (hv(ij)   + hv(ijls)                              &
                  & -  hv(ijlw) - hv(ijsw)) *  0.5d0 * rx * rxt(ij)     &
                  & - (hu(ij)   + hu(ijlw)                              &
                  & +  hu(ijls) + hu(ijsw)) * 0.25d0 * hxyt(ij)         &
                  & - (hv(ij)   + hv(ijlw)                              &
                  & +  hv(ijls) + hv(ijsw)) * 0.25d0 * hyxt(ij)

             sxx(ij) = sxx(ij) + asmb(ij) * db(ij) * (exx - eyy)
             syy(ij) = syy(ij) + asmb(ij) * db(ij) * (eyy - exx)
             sxy(ij) = sxy(ij) + asmb(ij) * db(ij) * exy
          end do

          do ij = ijvstr, ijvend+1
             ijln = ij + ln
             fux(ij, k) = (sxx(ij)   * hyt(ij)   * hyt(ij) &
                  &     +  sxx(ijln) * hyt(ijln) * hyt(ijln)) * 0.5d0 * amfvx(ij, k)
             fvx(ij, k) = (sxy(ij)   * hyt(ij)   * hyt(ij) &
                  &     +  sxy(ijln) * hyt(ijln) * hyt(ijln)) * 0.5d0 * amfvx(ij, k)
          end do

          do ij = ijvstr, ijvend+nxdim
             ijle = ij + le
             fuy(ij, k) = (sxy(ij)   * hxt(ij)   * hxt(ij) &
                  &     +  sxy(ijle) * hxt(ijle) * hxt(ijle)) * 0.5d0 * amfvy(ij, k)
             fvy(ij, k) = (syy(ij)   * hxt(ij)   * hxt(ij) &
                  &     +  syy(ijle) * hxt(ijle) * hxt(ijle)) * 0.5d0 * amfvy(ij, k)
          end do

          do ij = ijvstr, ijvend
             ez =    ((amv(ij, k) + amv(ij, k+1)) * 0.5d0 / rea &
                  & + (amv(ij, k) - amv(ij, k+1)) * rz(ij, k))
             szx(ij, k) = - ux(ij, k) * ez
             szy(ij, k) = - vx(ij, k) * ez
          end do
       end do

       do k = kstr, kstr+kz-1
          do ij = ijvstr, ijvend
             gx(ij, k) = (gx(ij, k)                                       &
                  &  + ((fux(ij+le, k) - fux(ij, k)) * rx * ryu(ij)       & 
                  &  +  (fuy(ij+ln, k) - fuy(ij, k)) * rym(ij) * rxu(ij)) &
                  &                                  * rxu(ij) * ryu(ij)  &
                  &  +   szx(ij, k) / rea) * amskv(ij, k)
             gy(ij, k) = (gy(ij, k)                                       &
                  &  + ((fvx(ij+le, k) - fvx(ij, k)) * rx * ryu(ij)       &
                  &  +  (fvy(ij+ln, k) - fvy(ij, k)) * rym(ij) * rxu(ij)) &
                  &                                  * rxu(ij) * ryu(ij)  &
                  &  +   szy(ij, k) / rea ) * amskv(ij, k)
          end do
       end do

       do k = kstr+kz, kend
          do ij = ijvstr, ijvend
             gx(ij, k) = (gx(ij, k)                                                  &
                  &  + ((fux(ij+le, k) - fux(ij, k)) * rx * ryu(ij)                  &
                  &  +  (fuy(ij+ln, k) - fuy(ij, k)) * rym(ij) * rxu(ij))            &
                  &                                  * rxu(ij) * ryu(ij) * rz(ij, k) &
                  &  +   szx(ij, k) / rea ) * amskv(ij, k)
             gy(ij, k) = (gy(ij, k)                                                  &
                  &  + ((fvx(ij+le, k) - fvx(ij, k)) * rx * ryu(ij)                  &
                  &  +  (fvy(ij+ln, k) - fvy(ij, k)) * rym(ij) * rxu(ij))            &
                  &                                  * rxu(ij) * ryu(ij) * rz(ij, k) &
                  &  +   szy(ij, k) / rea ) * amskv(ij, k)
          end do
       end do
       !$acc end kernels
    end if

    !$acc kernels default(present)
    do k = kstr, kend
       do ij = ijvstr, ijvend
          xx(ij, k) = gx(ij, k)
          yy(ij, k) = gy(ij, k)
       end do
    end do
    !$acc end kernels
  end subroutine vscvel


#ifdef OPT_BBL
! *********************************************************************
! --- information -----------------------------------------------------
!
!  Visocity term of the BBL momentum eqs.
!
! ---------------------------------------------------------------------
  subroutine vscvlb(                                                  &
         &      gx,     gy,     xx,     yy,                           &
         &      uy,     ux,     vy,     vx,                           &
         &     amv  )
    use zocdim,  only :                                               &
         &  nxydim,  nxdim,  nzdim,                                   &
         &    kstr,   kend,                                           &
         &   ijstr,  ijend, ijvstr, ijvend,                           &
         &      le,     lw,     ln,     ls,   lsw,                    &
         &   oinit,  ofinal
    use zocgrd,  only :                                               &
         &    zbot,     dx,                                           &  
         &      dy,    dzv,    dzm,    rea,                           &  
         &      rx,     ry,    rym,     rs,   rsm,                    & 
         &     hxt,    hyt,   hxyt,   hyxt,                           &
         &     rxu,    ryu,    rxt,    ryt
    use zocmsk,  only :  amskvb,  amfvx, amfvy, nbotv
    use ufile
  
    implicit none

    real(8),   intent(inout)  ::     gx(nxydim,nzdim),    gy(nxydim,nzdim)
    real(8),   intent(inout)  ::     xx(nxydim,nzdim),    yy(nxydim,nzdim)
    real(8),   intent(in)     ::     uy(nxydim,nzdim),    ux(nxydim,nzdim)
    real(8),   intent(in)     ::     vy(nxydim,nzdim),    vx(nxydim,nzdim)
    real(8),   intent(in)     ::    amv(nxydim,nzdim)

!---- local variables
    real(8)            ::    fux(nxydim),    fvx(nxydim)
    real(8)            ::    fuy(nxydim),    fvy(nxydim)
    real(8)            ::    fuz(nxydim),    fvz(nxydim)
    real(8),     save  ::     rz(nxydim),    rzm(nxydim)
    real(8)            ::    szx(nxydim),    szy(nxydim)
    real(8)            ::     db(nxydim)
    real(8)            ::     hu(nxydim),    hv(nxydim)
    real(8)            ::    hux(nxydim),   hvx(nxydim)
    real(8)            ::    huy(nxydim),   hvy(nxydim)

    real(8)            ::    exx,    eyy,    exy,     ez,    dd
    integer(4)         ::     ij,      k,      i,      j
    integer(4)         ::   ijln,   ijls,   ijle,   ijlw
    integer(4)         ::   ijnw,   ijse,   ijsw,   ijne
    integer(4)         ::  ifpar,  jfpar,  istat
    integer(4)         ::    kup

    real(8),     save  ::    asm(nxydim),  asmb(nxydim)
    real(8),     save  :: csmbbl=0.d0, csmbbb=0.d0
    real(8)            ::     pi
    namelist /nmsmgb/ csmbbl, csmbbb

    if ( oinit .or. ofinal ) then
       return
    end if
    
    if ( ofirst_bbl ) then
       ofirst_bbl = .false.
       READ_NAMELIST( nmsmgb )

       pi = 4.d0 * atan(1.d0)
       do ij = 1, nxydim
          asm(ij) = (csmbbl * min(dx*hxt(ij), dy(ij)*hyt(ij)) / pi)**2
          asmb(ij) = csmbbb * min(dx*hxt(ij), dy(ij)*hyt(ij))**2       &
     &            / pi / sqrt(8.d0)

          rz (ij) = 1.d0 / dzv(ij, kend)
          rzm(ij) = 1.d0 / dzm(ij, kend)
       end do

    end if
       
    k = kend

    do ij = 1, nxydim
       fux(ij) = 0.d0
       fvx(ij) = 0.d0
       fuy(ij) = 0.d0
       fvy(ij) = 0.d0
       sxx(ij) = 0.d0
       syy(ij) = 0.d0
       sxy(ij) = 0.d0
    end do

    do ij = ijvstr, ijvend
       kup = max( nbotv(ij)-1, 1 )
       fuz(ij) = amv(ij, kend) * rzm(ij) * (ux(ij, kup) - ux(ij, kend))
       fvz(ij) = amv(ij, kend) * rzm(ij) * (vx(ij, kup) - vx(ij, kend))
    end do

    do ij = ijstr-nxdim-1, ijend+nxdim+nxdim+2
       ijls = ij + ls
       ijlw = ij + lw
       ijsw = ij + lsw
       exx = (  ux(ij  , k) + ux(ijls, k)                            &
    &         - ux(ijlw, k) - ux(ijsw, k)) * 0.5d0 *                 &
    &        rx * rxt(ij)                                            &
    &      + (  vx(ij  , k) + vx(ijls, k)                            &
    &         + vx(ijlw, k) + vx(ijsw, k)) * 0.25d0 *                &
    &        hxyt(ij)     
       eyy = (  vx(ij  , k) + vx(ijlw, k)                            &
    &         - vx(ijls, k) - vx(ijsw, k)) * 0.5d0 *                 &
    &        ry(ij) * ryt(ij)                                        &
    &      + (  ux(ij  , k) + ux(ijlw, k)                            &
    &         + ux(ijls, k) + ux(ijsw, k)) * 0.25d0 *                &
    &        hyxt(ij)
       exy = (  ux(ij  , k) + ux(ijlw, k)                            &
    &         - ux(ijls, k) - ux(ijsw, k)) * 0.5d0 *                 &
    &        ry(ij) * ryt(ij)                                        &
    &      + (  vx(ij  , k) + vx(ijls, k)                            &
    &         - vx(ijlw, k) - vx(ijsw, k)) * 0.5d0 *                 &
    &        rx * rxt(ij)                                            &
    &      - (  ux(ij  , k) + ux(ijlw, k)                            &
    &         + ux(ijls, k) + ux(ijsw, k)) * 0.25d0 *                &    
    &        hxyt(ij)                                                &
    &      - (  vx(ij  , k) + vx(ijlw, k)                            &
    &         + vx(ijls, k) + vx(ijsw, k)) * 0.25d0 *                &
    &        hyxt(ij)

       dd = sqrt((exx - eyy) * (exx - eyy) + exy * exy)
       sxx(ij) = asm(ij) * dd * (exx - eyy)
       syy(ij) = asm(ij) * dd * (eyy - exx)
       sxy(ij) = asm(ij) * dd * exy
       
       db(ij) = sqrt(sqrt((exx - eyy) * (exx - eyy) + exy * exy))
       bxx(ij) = asmb(ij) * db(ij) * (exx - eyy)
       byy(ij) = asmb(ij) * db(ij) * (eyy - exx)
       bxy(ij) = asmb(ij) * db(ij) * exy
    end do
     
    do ij = ijstr-nxdim-1, ijend+nxdim+2
       ijln = ij + ln
       hux(ij) = (  bxx(ij) * hyt(ij) * hyt(ij)                       &
    &             + bxx(ijln) * hyt(ijln) * hyt(ijln)) * 0.5d0
       hvx(ij) = (  bxy(ij) * hyt(ij) * hyt(ij)                       &
    &             + bxy(ijln) * hyt(ijln) * hyt(ijln)) * 0.5d0
    end do

    do ij = ijstr-nxdim-1, ijend+nxdim+nxdim+1
       ijle = ij + le
       huy(ij) = (  bxy(ij) * hxt(ij) * hxt(ij)                       &
    &             + bxy(ijle) * hxt(ijle) * hxt(ijle)) * 0.5d0
       hvy(ij) = (  byy(ij) * hxt(ij) * hxt(ij)                       &
    &             + byy(ijle) * hxt(ijle) * hxt(ijle)) * 0.5d0
    end do

    do ij = ijstr-nxdim-1, ijend+nxdim+1
       hu(ij) = - (  (hux(ij+le) - hux(ij)) *                         &
    &                 rx * ryu(ij)                                    &
    &              + (huy(ij+ln) - huy(ij)) *                         &
    &                 rym(ij) * rxu(ij)) *                            &
    &              rxu(ij) * ryu(ij) * amskvb(ij)
       hv(ij) = - (  (hvx(ij+le) - hvx(ij)) *                         &
    &                 rx * ryu(ij)                                    &
    &              + (hvy(ij+ln) - hvy(ij)) *                         &
    &                 rym(ij) * rxu(ij)) *                            &
    &              rxu(ij) * ryu(ij) * amskvb(ij)
    end do

    do ij = ijstr, ijend+nxdim+1
       ijls = ij + ls
       ijlw = ij + lw
       ijsw = ij + lsw
       exx = (  hu(ij  ) + hu(ijls)                                  &
    &         - hu(ijlw) - hu(ijsw)) * 0.5d0 *                       &
    &         rx * rxt(ij)                                           &
    &      + (  hv(ij  ) + hv(ijls)                                  &
    &         + hv(ijlw) + hv(ijsw)) * 0.25d0 *                      &
    &         hxyt(ij)
        eyy = (  hv(ij  ) + hv(ijlw)                                 &
    &          - hv(ijls) - hv(ijsw)) * 0.5d0 *                      &
    &         ry(ij) * ryt(ij)                                       &
    &       + (  hu(ij  ) + hu(ijlw)                                 &
    &          + hu(ijls) + hu(ijsw)) * 0.25d0 *                     &
    &         hyxt(ij)
        exy = (  hu(ij  ) + hu(ijlw)                                 &
    &          - hu(ijls) - hu(ijsw)) * 0.5d0 *                      &
    &         ry(ij) * ryt(ij)                                       &
    &       + (  hv(ij  ) + hv(ijls)                                 &
    &          - hv(ijlw) - hv(ijsw)) * 0.5d0 *                      &
    &         rx * rxt(ij)                                           &
    &       - (  hu(ij  ) + hu(ijlw)                                 &
    &          + hu(ijls) + hu(ijsw)) * 0.25d0 *                     &
    &         hxyt(ij)                                               &
    &       - (  hv(ij  ) + hv(ijlw)                                 & 
    &          + hv(ijls) + hv(ijsw)) * 0.25d0 *                     &
    &         hyxt(ij)

        sxx(ij) = sxx(ij) + asmb(ij) * db(ij) * (exx - eyy)
        syy(ij) = syy(ij) + asmb(ij) * db(ij) * (eyy - exx)
        sxy(ij) = sxy(ij) + asmb(ij) * db(ij) * exy
        
     end do
     
     do ij = ijvstr, ijvend+1
        ijln = ij + ln
        fux(ij) = (  sxx(ij)   * hyt(ij) * hyt(ij)                    &
    &              + sxx(ijln) * hyt(ijln) * hyt(ijln))               &
    &           * 0.5d0 * amfvx(ij, k)
        fvx(ij) = (  sxy(ij) * hyt(ij) * hyt(ij)                      &
    &              + sxy(ijln) * hyt(ijln) * hyt(ijln))               &
    &           * 0.5d0 * amfvx(ij, k)
     end do
     do ij = ijvstr, ijvend+nxdim
        ijle = ij + le
        fuy(ij) = (  sxy(ij) * hxt(ij) * hxt(ij)                      &
    &              + sxy(ijle) * hxt(ijle) * hxt(ijle))               &
    &           * 0.5d0 * amfvy(ij, k)
        fvy(ij) = (  syy(ij) * hxt(ij) * hxt(ij)                      &
    &              + syy(ijle) * hxt(ijle) * hxt(ijle))               &
    &           * 0.5d0 * amfvy(ij, k)
     end do

     do ij = ijvstr, ijvend
        szx(ij) = - ux(ij, k) * amv(ij, k) / rea
        szy(ij) = - vx(ij, k) * amv(ij, k) / rea
     end do

     do ij = ijvstr, ijvend
        gx(ij, k) = (  gx(ij, k)                                      &
    &                + fuz(ij) * rz(ij)                               &
    &                + (  (fux(ij+le) - fux(ij)) *                    &
    &                     rx * ryu(ij)                                &
    &                   + (fuy(ij+ln) - fuy(ij)) *                    &
    &                     rym(ij) * rxu(ij)                           &
    &                  ) * rxu(ij) * ryu(ij) * rz(ij)                 &  
    &                + szx(ij) / rea ) * amskvb(ij)
        gy(ij, k) = (  gy(ij, k)                                      &
    &                + fvz(ij) * rz(ij)                               &
    &                + (  (fvx(ij+le) - fvx(ij)) *                    &
    &                     rx * ryu(ij)                                &
    &                   + (fvy(ij+ln) - fvy(ij)) *                    &
    &                     rym(ij) * rxu(ij)                           &
    &                  ) * rxu(ij) * ryu(ij) * rz(ij)                 &
    &                + szy(ij) / rea ) * amskvb(ij)
     end do

     do ij = ijvstr, ijvend
        xx(ij, kend) = gx(ij, kend)
        yy(ij, kend) = gy(ij, kend)
     end do
     
  end subroutine vscvlb

#endif
  
end module cvisc

