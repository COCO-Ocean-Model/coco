
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

  private
  public  ::  vscvel, vscvlb  !  used in aprdc

  real(8),     save  ::    sxx(nxydim),          syy(nxydim)
  real(8),     save  ::    sxy(nxydim),          syx(nxydim)
  real(8),     save  ::  hvbot(nxydim)
  logical,     save  ::  ofirst, ofirst_bbl       
  character(len=64)  ::  chead(1:16)
  data ofirst, ofirst_bbl / .true., .true. /

!---- temporary arrays
  real(8),     allocatable,    dimension(:,:)     ::  buf2,    g2d

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
         &      le,     lw,     ln,     ls,    lne,    lsw,           &
         &   oinit,  ofinal
    use zocgrd,  only :                                               &
         &    zbot,     dx,                                           &  
         &      dy,    dzv,    dzm,    rea,                           &  
         &      rx,     ry,    rym,     rs,    rsm,                   & 
         &     hxt,    hyt,   hxyt,   hyxt,                           &
         &     rxu,    ryu,    rxt,    ryt   
    use zocmsk,  only :  amskv,  amfvx,  amfvy
    
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
    real(8)            ::     dd(nxydim,nzdim)
    real(8)            ::     hu(nxydim,nzdim),    hv(nxydim,nzdim)

    real(8)            ::    exx,    eyy,    exy,     ez
    integer(4)         ::     ij,      k,      i,      j
    integer(4)         ::   ijln,   ijls,   ijle,   ijlw
    integer(4)         ::   ijnw,   ijse,   ijsw,   ijlsw
    integer(4)         ::  ifpar,  jfpar,  istat

    real(8),     save  ::    asm(nxydim)
    real(8),     save  ::    csmb
    real(8)            ::     pi
    namelist /nmsmab/ csmb
    data csmb / 0.d0 /

    if ( oinit .or. ofinal ) then
       return
    end if

    if ( ofirst ) then

       ofirst = .false.
       call rewnml( ifpar, jfpar )
       read( ifpar, nmsmab, iostat = istat )
       call cstnml( jfpar, 'vscvel', 'nmsmab', istat )
       write( jfpar, nmsmab )

       pi = 4.d0 * atan(1.d0)
       do ij = 1, nxydim
          asm(ij) = csmb * min( dx*hxt(ij), dy(ij)*hyt(ij))**2        &
   &              / pi / sqrt(8.d0)
       end do
      
    end if

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

    do k = kstr, kend
       do ij = ijstr-nxdim-1, ijend+nxdim+nxdim+2

          ijls  = ij + ls
          ijlw  = ij + lw
          ijlsw = ij + lsw
          exx = (  ux(ij  , k) + ux(ijls , k)                         &
    &            - ux(ijlw, k) - ux(ijlsw, k)) * 0.5d0 *              &
    &            rx * rxt(ij)                                         &
    &         + (  vx(ij  , k) + vx(ijls , k)                         &
    &            + vx(ijlw, k) + vx(ijlsw, k)) * 0.25d0 *             &
    &            hxyt(ij)                         
          eyy = (  vx(ij  , k) + vx(ijlw , k)                         & 
    &            - vx(ijls, k) - vx(ijlsw, k)) * 0.5d0 *              &
    &           ry(ij) * ryt(ij)                                      &
    &         + (  ux(ij  , k) + ux(ijlw , k)                         &
    &            + ux(ijls, k) + ux(ijlsw, k)) * 0.25d0 *             &
    &            hyxt(ij)
          exy = (  ux(ij  , k) + ux(ijlw , k)                         &
    &            - ux(ijls, k) - ux(ijlsw, k)) * 0.5d0 *              &
    &           ry(ij) * ryt(ij)                                      &
    &         + (  vx(ij  , k) + vx(ijls , k)                         &
    &            - vx(ijlw, k) - vx(ijlsw, k)) * 0.5d0 *              &
    &           rx * rxt(ij)                                          &
    &         - (  ux(ij  , k) + ux(ijlw , k)                         &
    &            + ux(ijls, k) + ux(ijlsw, k)) * 0.25d0 *             &
    &           hxyt(ij)                                              &
    &         - (  vx(ij  , k) + vx(ijlw , k)                         &
    &            + vx(ijls, k) + vx(ijlsw, k)) * 0.25d0 *             &
    &           hyxt(ij)

          dd(ij, k) = sqrt(sqrt(                                      &
    &                  (exx - eyy) * (exx - eyy) + exy * exy))

          sxx(ij) = asm(ij) * dd(ij, k) * (exx - eyy)
          syy(ij) = asm(ij) * dd(ij, k) * (eyy - exx)
          sxy(ij) = asm(ij) * dd(ij, k) * exy

       end do

       do ij = ijstr-nxdim-1, ijend+nxdim+2
          ijln = ij + ln
          fux(ij, k) = (  sxx(ij)   * hyt(ij)   * hyt(ij)             &
    &                   + sxx(ijln) * hyt(ijln) * hyt(ijln))          & 
    &                  * 0.5d0
          fvx(ij, k) = (  sxy(ij)   * hyt(ij)   * hyt(ij)             &
    &                   + sxy(ijln) * hyt(ijln) * hyt(ijln))          &
    &                  * 0.5d0
       end do

       do ij = ijstr-nxdim-1, ijend+nxdim+nxdim+1
          ijle = ij + le
          fuy(ij, k) = (  sxy(ij)   * hxt(ij)   * hxt(ij)             &
    &                   + sxy(ijle) * hxt(ijle) * hxt(ijle))          &
    &                * 0.5d0
          fvy(ij, k) = (  syy(ij)   * hxt(ij)   * hxt(ij)             &
    &                   + syy(ijle) * hxt(ijle) * hxt(ijle))          &
    &                * 0.5d0
       end do

    end do

    do k = kstr, kend
       do ij = ijstr-nxdim-1, ijend+nxdim+1
          hu(ij, k) = - (  (fux(ij+le, k) - fux(ij, k)) *             &
    &                       rx * ryu(ij)                              &
    &                    + (fuy(ij+ln, k) - fuy(ij, k)) *             &   
    &                       rym(ij) * rxu(ij)) *                      &
    &                     rxu(ij) * ryu(ij) * amskv(ij, k)            
          hv(ij, k) = - (  (fvx(ij+le, k) - fvx(ij, k)) *             &
    &                       rx * ryu(ij)                              &
    &                    + (fvy(ij+ln, k) - fvy(ij, k)) *             &
    &                       rym(ij) * rxu(ij)) *                      &
    &                  rxu(ij) * ryu(ij) * amskv(ij, k)
       end do
    end do

    do k = kstr, kend
       do ij = ijstr, ijend+nxdim+1
          ijls = ij + ls
          ijlw = ij + lw
          ijlsw = ij + lsw
          exx = (  hu(ij  , k) + hu(ijls , k)                          &
    &            - hu(ijlw, k) - hu(ijlsw, k)) * 0.5d0 *               &
    &           rx * rxt(ij)                                           &
    &         + (  hv(ij  , k) + hv(ijls , k)                          &
    &            + hv(ijlw, k) + hv(ijlsw, k)) * 0.25d0 *              &
    &           hxyt(ij)                                            
          eyy = (  hv(ij  , k) + hv(ijlw , k)                          & 
    &            - hv(ijls, k) - hv(ijlsw, k)) * 0.5d0 *               &
    &           ry(ij) * ryt(ij)                                       &
    &         + (  hu(ij  , k) + hu(ijlw , k)                          &
    &            + hu(ijls, k) + hu(ijlsw, k)) * 0.25d0 *              &
    &           hyxt(ij)                                            
          exy = (  hu(ij  , k) + hu(ijlw , k)                          &
    &            - hu(ijls, k) - hu(ijlsw, k)) * 0.5d0 *               &
    &           ry(ij) * ryt(ij)                                       &
    &         + (  hv(ij  , k) + hv(ijls , k)                          &
    &            - hv(ijlw, k) - hv(ijlsw, k)) * 0.5d0 *               &
    &           rx * rxt(ij)                                           &
    &         - (  hu(ij  , k) + hu(ijlw , k)                          &
    &            + hu(ijls, k) + hu(ijlsw, k)) * 0.25d0 *              &
    &           hxyt(ij)                                               &
    &         - (  hv(ij  , k) + hv(ijlw , k)                          &
    &            + hv(ijls, k) + hv(ijlsw, k)) * 0.25d0 *              &
    &           hyxt(ij)

          sxx(ij) = asm(ij) * dd(ij, k) * (exx - eyy)
          syy(ij) = asm(ij) * dd(ij, k) * (eyy - exx)
          sxy(ij) = asm(ij) * dd(ij, k) * exy
       end do

       do ij = ijvstr, ijvend+1
          ijln = ij + ln
          fux(ij, k) = (  sxx(ij)   * hyt(ij)   * hyt(ij)             &
    &                   + sxx(ijln) * hyt(ijln) * hyt(ijln))          & 
    &                * 0.5d0 * amfvx(ij, k)
          fvx(ij, k) = (  sxy(ij)   * hyt(ij)   * hyt(ij)             &
    &                   + sxy(ijln) * hyt(ijln) * hyt(ijln))          &
    &                * 0.5d0 * amfvx(ij, k)
       end do

       do ij = ijvstr, ijvend+nxdim
          ijle = ij + le
          fuy(ij, k) = (  sxy(ij)   * hxt(ij)   * hxt(ij)             &
    &                   + sxy(ijle) * hxt(ijle) * hxt(ijle))          &
    &                * 0.5d0 * amfvy(ij, k)                           
          fvy(ij, k) = (  syy(ij)   * hxt(ij)   * hxt(ij)             &
    &                   + syy(ijle) * hxt(ijle) * hxt(ijle))          &
    &                * 0.5d0 * amfvy(ij, k)
       end do

    end do

    do k = kstr, kend
       do ij = ijvstr, ijvend
          ez = (  (amv(ij, k) + amv(ij, k+1)) * 0.5d0 / rea           &
     &          + (amv(ij, k) - amv(ij, k+1)) * rz(ij, k))
          szx(ij, k) = - ux(ij, k) * ez
          szy(ij, k) = - vx(ij, k) * ez
       end do
    end do

    do k = kstr, kstr+kz-1
       do ij = ijvstr, ijvend
          gx(ij, k) = (  gx(ij, k)                                    &
    &                + (  (fux(ij+le, k) - fux(ij, k)) *              &
    &                      rx * ryu(ij)                               & 
    &                   + (fuy(ij+ln, k) - fuy(ij, k)) *              &
    &                      rym(ij) * rxu(ij)) * rxu(ij) * ryu(ij)     &
    &                + szx(ij, k) / rea ) * amskv(ij, k)
          gy(ij, k) = (  gy(ij, k)                                    &
    &                + (  (fvx(ij+le, k) - fvx(ij, k)) *              &
    &                      rx * ryu(ij)                               &
    &                   + (fvy(ij+ln, k) - fvy(ij, k)) *              & 
    &                      rym(ij) * rxu(ij)) * rxu(ij) * ryu(ij)     &
    &                + szy(ij, k) / rea ) * amskv(ij, k)
       end do
    end do
      
    do k = kstr+kz, kend
       do ij = ijvstr, ijvend
          gx(ij, k) = (  gx(ij, k)                                   &
    &                + (  (fux(ij+le, k) - fux(ij, k)) *             &
    &                      rx * ryu(ij)                              &
    &                   + (fuy(ij+ln, k) - fuy(ij, k)) *             &
    &                      rym(ij) * rxu(ij)                         &
    &                   ) * rxu(ij) * ryu(ij) * rz(ij, k)            &
    &                + szx(ij, k) / rea ) * amskv(ij, k)
          gy(ij, k) = (  gy(ij, k)                                   &
    &                + (  (fvx(ij+le, k) - fvx(ij, k)) *             &
    &                     rx * ryu(ij)                               &
    &                   + (fvy(ij+ln, k) - fvy(ij, k)) *             & 
    &                     rym(ij) * rxu(ij)                          &
    &                   ) * rxu(ij) * ryu(ij) * rz(ij, k)            &
    &                + szy(ij, k) / rea ) * amskv(ij, k)
       end do
    end do

    do k = kstr, kend
       do ij = ijvstr, ijvend
          xx(ij, k) = gx(ij, k)
          yy(ij, k) = gy(ij, k)
       end do
    end do

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
    use zocmsk,  only :  amskv,   amskvb,  amfvx, amfvy, nbotv
    
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
    real(8)            ::     dd(nxydim)
    real(8)            ::     hu(nxydim),    hv(nxydim)
    real(8)            ::    exx,    eyy,    exy,     ez
    integer(4)         ::     ij,      k,      i,      j
    integer(4)         ::   ijln,   ijls,   ijle,   ijlw
    integer(4)         ::   ijnw,   ijse,   ijsw,   ijlsw
    integer(4)         ::  ifpar,  jfpar,  istat
    integer(4)         ::    kup

    real(8),     save  ::    asm(nxydim)
    real(8),     save  :: csmbbb
    real(8)            ::     pi
    namelist /nmsmmb/ csmbbb
    data csmbbb / 0.d0 /

    if ( oinit .or. ofinal ) then
       return
    end if
    
    if ( ofirst_bbl ) then

       ofirst_bbl = .false.
       call rewnml( ifpar, jfpar )
       read( ifpar, nmsmmb, iostat = istat )
       call cstnml( jfpar, 'vscvlb', 'nmsmgb', istat )
       write( jfpar, nmsmmb )

       pi = 4.d0 * atan(1.d0)
       do ij = 1, nxydim
          asm(ij) = csmbbb * min(dx*hxt(ij), dy(ij)*hyt(ij))**2       &
     &            / pi / sqrt(8.d0)

          rz (ij) = 1.d0 / dzv(ij, kend)
          rzm(ij) = 1.d0 / dzm(ij, kend)
       end do

    end if
       
    k = kend

    do ij = ijvstr, ijvend
       kup = max( nbotv(ij)-1, 1 )
       fuz(ij) = amv(ij, kend) * rzm(ij) * (ux(ij, kup) - ux(ij, kend))
       fvz(ij) = amv(ij, kend) * rzm(ij) * (vx(ij, kup) - vx(ij, kend))
    end do

    do ij = ijstr-nxdim-1, ijend+nxdim+nxdim+2
       ijls = ij + ls
       ijlw = ij + lw
       ijlsw = ij + lsw
       exx = (  ux(ij  , k) + ux(ijls , k)                            &
    &         - ux(ijlw, k) - ux(ijlsw, k)) * 0.5d0 *                 &
    &        rx * rxt(ij)                                             &
    &      + (  vx(ij  , k) + vx(ijls , k)                            &
    &         + vx(ijlw, k) + vx(ijlsw, k)) * 0.25d0 *                &
    &        hxyt(ij)     
       eyy = (  vx(ij  , k) + vx(ijlw , k)                            &
    &         - vx(ijls, k) - vx(ijlsw, k)) * 0.5d0 *                 &
    &        ry(ij) * ryt(ij)                                         &
    &      + (  ux(ij  , k) + ux(ijlw , k)                            &
    &         + ux(ijls, k) + ux(ijlsw, k)) * 0.25d0 *                &
    &        hyxt(ij)
       exy = (  ux(ij  , k) + ux(ijlw , k)                            &
    &         - ux(ijls, k) - ux(ijlsw, k)) * 0.5d0 *                 &
    &        ry(ij) * ryt(ij)                                         &
    &      + (  vx(ij  , k) + vx(ijls , k)                            &
    &         - vx(ijlw, k) - vx(ijlsw, k)) * 0.5d0 *                 &
    &        rx * rxt(ij)                                             &
    &      - (  ux(ij  , k) + ux(ijlw , k)                            &
    &         + ux(ijls, k) + ux(ijlsw, k)) * 0.25d0 *                &    
    &        hxyt(ij)                                                 &
    &      - (  vx(ij  , k) + vx(ijlw , k)                            &
    &         + vx(ijls, k) + vx(ijlsw, k)) * 0.25d0 *                &
    &        hyxt(ij)

       dd(ij) = sqrt(sqrt((exx - eyy) * (exx - eyy) + exy * exy))

       sxx(ij) = asm(ij) * dd(ij) * (exx - eyy)
       syy(ij) = asm(ij) * dd(ij) * (eyy - exx)
       sxy(ij) = asm(ij) * dd(ij) * exy

    end do
     
    do ij = ijstr-nxdim-1, ijend+nxdim+2
       ijln = ij + ln
       fux(ij) = (  sxx(ij) * hyt(ij) * hyt(ij)                       &
    &             + sxx(ijln) * hyt(ijln) * hyt(ijln)) * 0.5d0
       fvx(ij) = (  sxy(ij) * hyt(ij) * hyt(ij)                       &
    &             + sxy(ijln) * hyt(ijln) * hyt(ijln)) * 0.5d0
    end do

    do ij = ijstr-nxdim-1, ijend+nxdim+nxdim+1
       ijle = ij + le
       fuy(ij) = (  sxy(ij) * hxt(ij) * hxt(ij)                       &
    &             + sxy(ijle) * hxt(ijle) * hxt(ijle)) * 0.5d0
       fvy(ij) = (  syy(ij) * hxt(ij) * hxt(ij)                       &
    &             + syy(ijle) * hxt(ijle) * hxt(ijle)) * 0.5d0
    end do

    do ij = ijstr-nxdim-1, ijend+nxdim+1
       hu(ij) = - (  (fux(ij+le) - fux(ij)) *                         &
    &                 rx * ryu(ij)                                    &
    &              + (fuy(ij+ln) - fuy(ij)) *                         &
    &                 rym(ij) * rxu(ij)) *                            &
    &              rxu(ij) * ryu(ij) * amskv(ij, k)
       hv(ij) = - (  (fvx(ij+le) - fvx(ij)) *                         &
    &                 rx * ryu(ij)                                    &
    &              + (fvy(ij+ln) - fvy(ij)) *                         &
    &                 rym(ij) * rxu(ij)) *                            &
    &              rxu(ij) * ryu(ij) * amskv(ij, k)
    end do

    do ij = ijstr, ijend+nxdim+1
       ijls = ij + ls
       ijlw = ij + lw
       ijlsw = ij + lsw
       exx = (  hu(ij  ) + hu(ijls )                                  &
    &         - hu(ijlw) - hu(ijlsw)) * 0.5d0 *                       &
    &         rx * rxt(ij)                                            &
    &      + (  hv(ij  ) + hv(ijls )                                  &
    &         + hv(ijlw) + hv(ijlsw)) * 0.25d0 *                      &
    &         hxyt(ij)
        eyy = (  hv(ij  ) + hv(ijlw )                                 &
    &          - hv(ijls) - hv(ijlsw)) * 0.5d0 *                      &
    &         ry(ij) * ryt(ij)                                        &
    &       + (  hu(ij  ) + hu(ijlw )                                 &
    &          + hu(ijls) + hu(ijlsw)) * 0.25d0 *                     &
    &         hyxt(ij)
        exy = (  hu(ij  ) + hu(ijlw )                                 &
    &          - hu(ijls) - hu(ijlsw)) * 0.5d0 *                      &
    &         ry(ij) * ryt(ij)                                        &
    &       + (  hv(ij  ) + hv(ijls )                                 &
    &          - hv(ijlw) - hv(ijlsw)) * 0.5d0 *                      &
    &         rx * rxt(ij)                                            &
    &       - (  hu(ij  ) + hu(ijlw )                                 &
    &          + hu(ijls) + hu(ijlsw)) * 0.25d0 *                     &
    &         hxyt(ij)                                                &
    &       - (  hv(ij  ) + hv(ijlw )                                 & 
    &          + hv(ijls) + hv(ijlsw)) * 0.25d0 *                     &
    &         hyxt(ij)

        sxx(ij) = asm(ij) * dd(ij) * (exx - eyy)
        syy(ij) = asm(ij) * dd(ij) * (eyy - exx)
        sxy(ij) = asm(ij) * dd(ij) * exy
        
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
    &                + szx(ij) / rea ) * amskv(ij,k)
        gy(ij, k) = (  gy(ij, k)                                      &
    &                + fvz(ij) * rz(ij)                               &
    &                + (  (fvx(ij+le) - fvx(ij)) *                    &
    &                     rx * ryu(ij)                                &
    &                   + (fvy(ij+ln) - fvy(ij)) *                    &
    &                     rym(ij) * rxu(ij)                           &
    &                  ) * rxu(ij) * ryu(ij) * rz(ij)                 &
    &                + szy(ij) / rea ) * amskv(ij,k)
     end do

     do ij = ijvstr, ijvend
        xx(ij, kend) = gx(ij, kend)
        yy(ij, kend) = gy(ij, kend)
     end do
     
  end subroutine vscvlb

#endif
  
end module cvisc

