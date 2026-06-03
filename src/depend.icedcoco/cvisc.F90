
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
!     '12.06.28  H.Tatebe: for COCO5.0 in F90
! ---------------------------------------------------------------------

  use zocdim,  only :  nxydim,  nzdim

  implicit none
#include "coco.h"
  private

  real(8),     save  ::    sxx(nxydim),          syy(nxydim)
  real(8),     save  ::    sxy(nxydim),          syx(nxydim)
  real(8),     save  ::    szx(nxydim),          szy(nxydim)
  real(8),     save  ::  hvbot(nxydim)
  logical,     save  ::  ofirst, ofirst_bbl       
  character(len=16)  ::  chead(1:64)
  data ofirst, ofirst_bbl / .true., .true. /
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
         &     nxg,    nyg, nxgdim, nygdim,                           &
         &  nxydim,  nxdim,  nydim,  nzdim,                           &
         &    kstr,   kend,     kz,                                   &
         &   ijstr,  ijend, ijvstr, ijvend,                           &
         &   igstr,  jgstr,  igend,  jgend,                           &
         &      le,     lw,     ln,     ls,                           &
         &     lnw,    lne,    lse,    lsw,                           &
         &   oinit,  ofinal
    use zocgrd,  only :                                               &
         &    zbot,                                                   &  
         &      dy,    dzv,    dzm,    hxt,   rea,                    &  
         &      rx,     ry,    rym,     rs,   rsm,                    & 
         &     hxu,    hyu,   hxyu,   hyxu,                           &
         &     rxu,    ryu   
    use zocnod,  only :  iroot,  myrank
    use zocmsk,  only :  amskv,  amfvx,  amfvy
    use zocfil,  only :  ncf
    use ufile
#ifdef OPT_IO_COCOMPI
    use mpiio
#else
    use bgs2d
#endif
    use bshft
 
    implicit none
#include "mpif.h"

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
    real(8)            ::    smx(nxydim,nzdim),    smy(nxydim,nzdim)
    real(8),     save  ::     rz(nxydim,nzdim),    rzm(nxydim,nzdim)
    integer(4)         ::     ij,      k,   i,      j
    integer(4)         ::   ijln,   ijls,   ijle,   ijlw
    integer(4)         ::   ijnw,   ijse,   ijsw,   ijne
    integer(4)         ::  ifpar,  jfpar,  istat

    real(8),     save  ::  amh
    real(8),     save  ::  amhmod(nxydim)
    integer(4)         ::  iam,    nfamh
    character(len=ncf) ::  cfamh
    namelist /nmvish/ amh
    namelist /nmcvis/ iam, cfamh
    data amh, iam, cfamh / 0.d0, -1, 'not-specified'/
    logical,     save  ::  opslvis = .false.
    namelist /nmpslv/ opslvis

#ifdef OPT_IO_COCOMPI
    integer :: mpi_fh
    integer :: icread
    integer (kind = mpi_offset_kind) :: disp
#else
    real(8), allocatable :: buf2(:, :),  g2d(:, :)
#endif

    real(8)            ::    sxxe(nxydim),        sxxw(nxydim)
    real(8)            ::    sxyn(nxydim),        sxys(nxydim)
    real(8)            ::    syxe(nxydim),        syxw(nxydim)
    real(8)            ::    syyn(nxydim),        syys(nxydim)
    real(8)            ::   smx2d(nxydim),       smy2d(nxydim)
    real(8)            ::      wt

    if ( oinit .or. ofinal ) then
       return
    end if

    if ( ofirst ) then
       ofirst = .false.
       READ_NAMELIST( nmvish )
       READ_NAMELIST( nmcvis )
       READ_NAMELIST( nmpslv )

       if ( iam < 0 ) then
!---- spatially constant
          write(jfpar,*) 'spatially constant AMH'
          amhmod(1:nxydim) = amh
       end if
       if ( iam == 0 ) then
!---- zonal resolution dependent (cvisc.clat setting)
          write(jfpar,*) 'dx dependent AMH'
          do ij = 1, nxydim
             amhmod(ij) = amh * hxt(ij) / rea
          end do
       end if
       if ( iam == 1 ) then
!---- given by input file
          write(jfpar,*) 'AMH is given by a file'          
#ifdef OPT_IO_COCOMPI
          call mpi_filopn(mpi_fh, cfamh, 'READ')
          disp=0
          call mpi_read_chead(chead, mpi_fh, disp, icread)
          call mpi_read_2d(amhmod, mpi_fh  , disp)
          call mpi_filcls(mpi_fh)
#else
          allocate ( buf2(1:nxg,1:nyg), g2d(1:nxgdim,1:nygdim) )
          buf2(1:nxg,1:nyg) = 0.d0
          g2d (1:nxgdim,1:nygdim) = 0.d0
          if ( myrank == iroot ) then
             call filopn( nfamh, cfamh, 'READ' )
!---- for MIROC
!               CALL IFLOPN(
!     O                     NFAMH,    IERR,
!     I                     CFAMH,  'READ', 'UNFORMATTED') 
             rewind( nfamh )
             read( nfamh ) chead
             read( nfamh ) buf2
             do j = 1, nyg
                do i = 1, nxg
                   g2d(igstr+i-1, jgstr+j-1) = buf2(i, j)
                end do
             end do
             call filcls( nfamh )
!---- for MIROC
!               CLOSE(NFAMH)
          end if
          call scatter_2d( amhmod, g2d )
          deallocate ( buf2, g2d )
#endif

#ifdef OPT_TRIPOLE
          call shift1(amhmod,                                      &
    &                  nxdim,  nydim,      1,                      &
    &                   1.d0,      0,      0 )
#else
          call shift1(                                             &
    &                 amhmod,                                      &
    &                  nxdim,  nydim,      1)
#endif
          
       end if
       
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

    if (opslvis) then ! partial step lateral viscosity
       do k = kstr, kend
          do ij = ijvstr, ijvend
             ijls = ij + ls
             ijlw = ij + lw
             ijln = ij + ln
             ijle = ij + le
             ijnw = ij + lnw
             ijse = ij + lse
             ijsw = ij + lsw
             ijne = ij + lne

             if (k .ge. kstr+kz) then
                if (dzv(ij, k) .ge. dzv(ijlw, k)) then
                   wt = dzv(ijlw, k) / dzv(ij, k)
                else
                   wt = 2.d0 - dzv(ij, k) / dzv(ijlw, k)
                end if
             else
                wt = 1.d0
             end if
             sxxw(ij) = (ux(ij, k) - ux(ijlw, k) * wt) * amhmod(ij) *       &
                  &                 rx / (hxu(ij) + hxu(ijlw)) * 2.d0                  &
                  &              + (vx(ij, k) + vx(ijlw, k) * wt) * amhmod(ij) *       &
                  &              (hxyu(ij) + hxyu(ijlw)) * 0.25d0

             if (k .ge. kstr+kz) then
                if (dzv(ijle, k) .ge. dzv(ij, k)) then
                   wt = 2.d0 - dzv(ij, k) / dzv(ijle, k)
                else
                   wt = dzv(ijle, k) / dzv(ij, k)
                endif
             else
                wt = 1.d0
             end if
             sxxe(ij) = (ux(ijle, k) * wt - ux(ij, k)) * amhmod(ij) *       &
                  &                    rx / (hxu(ijle) + hxu(ij)) * 2.d0               &
                  &                 + (vx(ijle, k) * wt + vx(ij, k)) * amhmod(ij) *    &
                  &                 (hxyu(ijle) + hxyu(ij)) * 0.25d0

             if (k .ge. kstr+kz) then
                if (dzv(ij, k) .ge. dzv(ijls, k)) then
                   wt = dzv(ijls, k) / dzv(ij, k)
                else
                   wt = 2.d0 - dzv(ij, k) / dzv(ijls, k)
                endif
             else
                wt = 1.d0
             end if
             sxys(ij) = (ux(ij, k) - ux(ijls, k) * wt) * amhmod(ij) *       &
                  &                   ry(ij) / (hyu(ij) + hyu(ijls)) * 2.d0            &
                  &                 - (vx(ij, k) + vx(ijls, k) * wt) * amhmod(ij) *    &
                  &                 (hyxu(ij) + hyxu(ijls)) * 0.25d0

             if (k .ge. kstr+kz) then
                if (dzv(ijln, k) .ge. dzv(ij, k)) then
                   wt = 2.d0 - dzv(ij, k) / dzv(ijln, k)
                else
                   wt = dzv(ijln, k) / dzv(ij, k)
                endif
             else
                wt = 1.d0
             end if
             sxyn(ij) = (ux(ijln, k) * wt - ux(ij, k)) * amhmod(ij) *       &
                  &                    ry(ijln) / (hyu(ijln) + hyu(ij)) * 2.d0         &
                  &                 - (vx(ijln, k) * wt + vx(ij, k)) * amhmod(ij) *    &
                  &                 (hyxu(ijln) + hyxu(ij)) * 0.25d0

             if (k .ge. kstr+kz) then
                if (dzv(ij, k) .ge. dzv(ijls, k)) then
                   wt = dzv(ijls, k) / dzv(ij, k)
                else
                   wt = 2.d0 - dzv(ij, k) / dzv(ijls, k)
                endif
             else
                wt = 1.d0
             end if
             syys(ij) = (vx(ij, k) - vx(ijls, k) * wt) * amhmod(ij) *       &
                  &                    ry(ij) / (hyu(ij) + hyu(ijls)) * 2.d0           &
                  &                 + (ux(ij, k) + ux(ijls, k) * wt) * amhmod(ij) *    &
                  &                 (hyxu(ij) + hyxu(ijls)) * 0.25d0

             if (k .ge. kstr+kz) then
                if (dzv(ijln, k) .ge. dzv(ij, k)) then
                   wt = 2.d0 - dzv(ij, k) / dzv(ijln, k)
                else
                   wt = dzv(ijln, k) / dzv(ij, k)
                endif
             else
                wt = 1.d0
             end if
             syyn(ij) = (vx(ijln, k) * wt - vx(ij, k)) * amhmod(ij) *       &
                  &                    ry(ijln) / (hyu(ijln) + hyu(ij)) * 2.d0         &
                  &                 + (ux(ijln, k) * wt + ux(ij, k)) * amhmod(ij) *    &
                  &                 (hyxu(ijln) + hyxu(ij)) * 0.25d0

             if (k .ge. kstr+kz) then
                if (dzv(ij, k) .ge. dzv(ijlw, k)) then
                   wt = dzv(ijlw, k) / dzv(ij, k)
                else
                   wt = 2.d0 - dzv(ij, k) / dzv(ijlw, k)
                endif
             else
                wt = 1.d0
             end if
             syxw(ij) = (vx(ij, k) - vx(ijlw, k) * wt) * amhmod(ij) *       &
                  &                 rx / (hxu(ij) + hxu(ijlw)) * 2.d0                  &
                  &              - (ux(ij, k) + ux(ijlw, k) * wt) * amhmod(ij) *       &
                  &              (hxyu(ij) + hxyu(ijlw)) * 0.25d0

             if (k .ge. kstr+kz) then
                if (dzv(ijle, k) .ge. dzv(ij, k)) then
                   wt = 2.d0 - dzv(ij, k) / dzv(ijle, k)
                else
                   wt = dzv(ijle, k) / dzv(ij, k)
                endif
             else
                wt = 1.d0
             end if
             syxe(ij) = (vx(ijle, k) * wt - vx(ij, k)) * amhmod(ij) *       &
                  &                 rx / (hxu(ijle) + hxu(ij)) * 2.d0                  &
                  &              - (ux(ijle, k) * wt + ux(ij, k)) * amhmod(ij) *       &
                  &              (hxyu(ijle) + hxyu(ij)) * 0.25d0

             smx2d(ij) = - ux(ij, k) *                                      &
                  &              (  (amv(ij, k) + amv(ij, k+1)) * 0.5d0 / rea          &
                  &               + (amv(ij, k) - amv(ij, k+1)) * rz(ij, k)) / rea
             smy2d(ij) = - vx(ij, k) *                                      &
                  &              (  (amv(ij, k) + amv(ij, k+1)) * 0.5d0 / rea          &
                  &               + (amv(ij, k) - amv(ij, k+1)) * rz(ij, k)) / rea

             gx(ij, k) = (  gx(ij, k)                                       &
                  &          + (  (  (hyu(ij+le) + hyu(ij)) * sxxe(ij)                 &
                  &                - (hyu(ij) + hyu(ij+lw)) * sxxw(ij)) * rx           &
                  &             + (  (hxu(ij+ln) + hxu(ij)) * sxyn(ij)                 &
                  &                - (hxu(ij) + hxu(ij+ls)) * sxys(ij)) * rym(ij)      &
                  &               ) * rxu(ij) * ryu(ij) * 0.5d0                        &
                  &                + hxyu(ij) * 0.5d0 * (syxe(ij) + syxw(ij))          &
                  &                - hyxu(ij) * 0.5d0 * (syyn(ij) + syys(ij))          &
                  &               + smx2d(ij)) * amskv(ij, k)
             gy(ij, k) = (  gy(ij, k)                                       &
                  &          + (  (  (hyu(ij+le) + hyu(ij)) * syxe(ij)                 &
                  &                - (hyu(ij) + hyu(ij+lw)) * syxw(ij)) * rx           &
                  &             + (  (hxu(ij+ln) + hxu(ij)) * syyn(ij)                 &
                  &                - (hxu(ij) + hxu(ij+ls)) * syys(ij)) * rym(ij)      &
                  &               ) * rxu(ij) * ryu(ij) * 0.5d0                        &
                  &                + hyxu(ij) * 0.5d0 * (sxyn(ij) + sxys(ij))          &
                  &                - hxyu(ij) * 0.5d0 * (sxxe(ij) + sxxw(ij))          &
                  &               + smy2d(ij)) * amskv(ij, k)

          end do
       end do
    else ! previous version
       do k = kstr, kend
          do ij = ijstr, ijend+nxdim+1
             ijls = ij + ls
             ijlw = ij + lw
             ijln = ij + ln
             ijle = ij + le
             ijnw = ij + lnw
             ijse = ij + lse
             ijsw = ij + lsw
             sxx(ij) = (ux(ij, k) - ux(ijlw, k)) * amhmod(ij) *          &
                  &               rx / (hxu(ij) + hxu(ijlw)) * 2.d0                 &
                  &             + (vx(ij, k) + vx(ijlw, k)) * amhmod(ij) *          &
                  &               (hxyu(ij) + hxyu(ijlw)) * 0.25d0                  &
                  &             - (  vx(ijln, k) + vx(ijnw, k)                      &
                  &                - vx(ijls, k) - vx(ijsw, k)) * amhmod(ij)        &
                  &               / (hyu(ij) + hyu(ijlw)) / (dy(ij) + dy(ijln))     &
                  &             - (ux(ij, k) + ux(ijlw, k)) * amhmod(ij) *          &
                  &               (hyxu(ij) + hyxu(ijlw)) * 0.25d0                  
             sxy(ij) = (ux(ij, k) - ux(ijls, k)) * amhmod(ij) *          &
                  &               ry(ij) / (hyu(ij) + hyu(ijls)) * 2.d0             &
                  &             - (ux(ij, k) + ux(ijls, k)) * amhmod(ij) *          &
                  &               (hxyu(ij) + hxyu(ijls)) * 0.25d0                  &
                  &             + (  vx(ijle, k) + vx(ijse, k)                      &
                  &                - vx(ijlw, k) - vx(ijsw, k)) * amhmod(ij)        &
                  &               / (hxu(ij) + hxu(ijls)) * rx * 0.5d0              &
                  &             - (vx(ij, k) + vx(ijls, k)) * amhmod(ij) *          &
                  &               (hyxu(ij) + hyxu(ijls)) * 0.25d0                  
             syy(ij) = (vx(ij, k) - vx(ijls, k)) * amhmod(ij) *          &
                  &               ry(ij) / (hyu(ij) + hyu(ijls)) * 2.d0             &
                  &             + (ux(ij, k) + ux(ijls, k)) * amhmod(ij) *          &
                  &               (hyxu(ij) + hyxu(ijls)) * 0.25d0                  &
                  &             - (  ux(ijle, k) + ux(ijse, k)                      &
                  &                - ux(ijlw, k) - ux(ijsw, k)) * amhmod(ij)        &
                  &               / (hxu(ij) + hxu(ijls)) * rx * 0.5d0              &
                  &             - (vx(ij, k) + vx(ijls, k)) * amhmod(ij) *          &
                  &               (hxyu(ij) + hxyu(ijls)) * 0.25d0                  
             syx(ij) = (vx(ij, k) - vx(ijlw, k)) * amhmod(ij) *          &
                  &               rx / (hxu(ij) + hxu(ijlw)) * 2.d0                 &
                  &             - (vx(ij, k) + vx(ijlw, k)) * amhmod(ij) *          &
                  &               (hyxu(ij) + hyxu(ijlw)) * 0.25d0                  &
                  &             + (  ux(ijln, k) + ux(ijnw, k)                      &
                  &                - ux(ijls, k) - ux(ijsw, k)) * amhmod(ij)        &
                  &               / (hyu(ij) + hyu(ijlw)) / (dy(ij) + dy(ijln))     &
                  &             - (ux(ij, k) + ux(ijlw, k)) * amhmod(ij) *          &
                  &               (hxyu(ij) + hxyu(ijlw)) * 0.25d0                  
             szx(ij) = - ux(ij, k) *                                     &
                  &                 (  (amv(ij, k) + amv(ij, k+1)) * 0.5d0 / rea    &
                  &                  + (amv(ij, k) - amv(ij, k+1)) * rz(ij, k))
             szy(ij) = - vx(ij, k) *                                     &
                  &                 (  (amv(ij, k) + amv(ij, k+1)) * 0.5d0 / rea    &
                  &                  + (amv(ij, k) - amv(ij, k+1)) * rz(ij, k))
          end do

          do ij = ijvstr, ijvend+1
             ijlw = ij + lw
             fux(ij, k) = sxx(ij) * (hyu(ij) + hyu(ijlw)) *              &
                  &                            (hyu(ij) + hyu(ijlw)) * 0.25d0 *     &
                  &                  amfvx(ij, k)
             fvx(ij, k) = syx(ij) * (hyu(ij) + hyu(ijlw)) *              &
                  &                            (hyu(ij) + hyu(ijlw)) * 0.25d0 *     &
                  &                   amfvx(ij, k)
          end do

          do ij = ijvstr, ijvend+nxdim
             ijls = ij + ls
             fuy(ij, k) = sxy(ij) * (hxu(ij) + hxu(ijls)) *              &
                  &                            (hxu(ij) + hxu(ijls)) * 0.25d0 *     &
                  &                  amfvy(ij, k)
             fvy(ij, k) = syy(ij) * (hxu(ij) + hxu(ijls)) *              &
                  &                            (hxu(ij) + hxu(ijls)) * 0.25d0 *     &
                  &                  amfvy(ij, k)
          end do

          do ij = ijvstr, ijvend
             smx(ij, k) = szx(ij) / rea
             smy(ij, k) = szy(ij) / rea
          end do

       end do

       do k = kstr, kstr+kz-1
          do ij = ijvstr, ijvend
             gx(ij, k) = (  gx(ij, k)                                    &
                  &                + (  (fux(ij+le, k) - fux(ij, k)) *              &
                  &                      rx * ryu(ij)                               & 
                  &                   + (fuy(ij+ln, k) - fuy(ij, k)) *              &
                  &                      rym(ij) * rxu(ij)) * rxu(ij) * ryu(ij)     &
                  &                + smx(ij, k) ) * amskv(ij, k)
             gy(ij, k) = (  gy(ij, k)                                    &
                  &                + (  (fvx(ij+le, k) - fvx(ij, k)) *              &
                  &                      rx * ryu(ij)                               &
                  &                   + (fvy(ij+ln, k) - fvy(ij, k)) *              & 
                  &                      rym(ij) * rxu(ij)) * rxu(ij) * ryu(ij)     &
                  &                + smy(ij, k) ) * amskv(ij, k)
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
                  &                + smx(ij, k) ) * amskv(ij, k)
             gy(ij, k) = (  gy(ij, k)                                   &
                  &                + (  (fvx(ij+le, k) - fvx(ij, k)) *             &
                  &                     rx * ryu(ij)                               &
                  &                   + (fvy(ij+ln, k) - fvy(ij, k)) *             & 
                  &                     rym(ij) * rxu(ij)                          &
                  &                   ) * rxu(ij) * ryu(ij) * rz(ij, k)            &
                  &                + smy(ij, k) ) * amskv(ij, k)
          end do
       end do
    end if
 
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
         &     nxg,    nyg, nxgdim, nygdim,                           &
         &  nxydim,  nxdim,  nydim,  nzdim,                           &
         &    kstr,   kend,     kz,                                   &
         &   ijstr,  ijend, ijvstr, ijvend,                           &
         &   igstr,  jgstr,  igend,  jgend,                           &
         &      le,     lw,     ln,     ls,                           &
         &     lnw,    lse,    lsw,                                   &
         &   oinit,  ofinal
    use zocgrd,  only :                                               &
         &    zbot,                                                   &  
         &      dy,    dzv,    dzm,    hxt,   rea,                    &  
         &      rx,     ry,    rym,     rs,   rsm,                    & 
         &     hxu,    hyu,   hxyu,   hyxu,                           &
         &     rxu,    ryu   
    use zocnod,  only :   iroot,  myrank
    use zocmsk,  only :  amskvb,  amfvx, amfvy, nbotv
    use zocfil,  only :  ncf
    use ufile
#ifdef OPT_IO_COCOMPI
    use mpiio
#else
    use bgs2d
#endif
    use bshft

    implicit none
#include "mpif.h"

    real(8),   intent(inout)  ::     gx(nxydim,nzdim),    gy(nxydim,nzdim)
    real(8),   intent(inout)  ::     xx(nxydim,nzdim),    yy(nxydim,nzdim)
    real(8),   intent(in)     ::     uy(nxydim,nzdim),    ux(nxydim,nzdim)
    real(8),   intent(in)     ::     vy(nxydim,nzdim),    vx(nxydim,nzdim)
    real(8),   intent(in)     ::    amv(nxydim,nzdim)

!---- local variables
    real(8)            ::    fux(nxydim),    fvx(nxydim)
    real(8)            ::    fuy(nxydim),    fvy(nxydim)
    real(8)            ::    fuz(nxydim),    fvz(nxydim)
    real(8)            ::    smx(nxydim),    smy(nxydim)
    real(8),     save  ::     rz(nxydim),    rzm(nxydim)
    integer(4)         ::     ij,      k,   i,      j
    integer(4)         ::   ijln,   ijls,   ijle,   ijlw
    integer(4)         ::   ijnw,   ijse,   ijsw
    integer(4)         ::    kup
    integer(4)         ::  ifpar,  jfpar,  istat

    real(8),     save  ::  amhmod(nxydim)
    real(8),     save  ::  amhbbl
    integer(4)         ::  iam,    nfamh
    character(len=ncf) ::  cfamh
    namelist /nmbbvh/ amhbbl, iam, cfamh
    data amhbbl, iam, cfamh / 0.d0, 0, 'not-specified' /

#ifdef OPT_IO_COCOMPI
    integer :: mpi_fh
    integer :: icread
    integer (kind = mpi_offset_kind) :: disp
#else
    real(8), allocatable :: buf2(:, :),  g2d(:, :)
#endif


    if ( oinit .or. ofinal ) then
       return
    end if
    
    if ( ofirst_bbl ) then

       ofirst_bbl = .false.
       READ_NAMELIST( nmbbvh )

       do ij = 1, nxydim
          rz   (ij) = 1.d0 / dzv(ij, kend)
          rzm  (ij) = 1.d0 / dzm(ij, kend)
       end do
       
       if ( iam == 0 ) then

          write(jfpar, *) ' amhbbl = ', amhbbl
          do ij = 1, nxydim
             amhmod(ij) = amhbbl * hxt(ij) / rea
          end do
          
       else

#ifdef OPT_IO_COCOMPI
          call mpi_filopn(mpi_fh, cfamh, 'read')
          disp=0
          call mpi_read_chead(chead, mpi_fh, disp, icread)
          call mpi_read_2d(amhmod, mpi_fh  , disp)
          call mpi_filcls(mpi_fh)
#else
!---- reading file of viscosity coefficient (i,j)
          allocate ( buf2(1:nxg,1:nyg) )
          allocate ( g2d(1:nxgdim,1:nygdim) )
          if ( myrank == iroot ) then
             call filopn( nfamh, cfamh, 'read' )
!---- for MIROC
!               CALL IFLOPN(
!     O                     NFAMH,    IERR,
!     I                     CFAMH,  'READ', 'UNFORMATTED') 
             rewind( nfamh )
             read( nfamh ) chead
             read( nfamh ) buf2
             do j = 1, nyg
                do i = 1, nxg
                   g2d(igstr+i-1, jgstr+j-1) = buf2(i, j)
                end do
             end do
             call filcls( nfamh )
!---- for MIROC
!               CLOSE(NFAMH)
          end if
          call scatter_2d( amhmod, g2d )
          deallocate ( buf2, g2d )
#endif
#ifdef OPT_TRIPOLE
          call shift1(amhmod,                                      &
    &                  nxdim,  nydim,      1,                      &
    &                   1.d0,      0,      0 )
#else
          call shift1(                                             &
    &                 amhmod,                                      &
    &                  nxdim,  nydim,      1)
#endif
       end if
       
    end if

    k = kend

    do ij = 1, nxydim
       fux(ij) = 0.d0
       fvx(ij) = 0.d0
       fuy(ij) = 0.d0
       fvy(ij) = 0.d0
    end do

    do ij = ijvstr, ijvend
       kup = max( nbotv(ij)-1, 1 )
       fuz(ij) = amv(ij, kend) * rzm(ij) * (ux(ij, kup) - ux(ij, kend))
       fvz(ij) = amv(ij, kend) * rzm(ij) * (vx(ij, kup) - vx(ij, kend))
    end do

    do ij = ijstr, ijend+nxdim+1
       ijls = ij + ls
       ijlw = ij + lw
       ijln = ij + ln
       ijle = ij + le
       ijnw = ij + lnw
       ijse = ij + lse
       ijsw = ij + lsw
       sxx(ij) = (ux(ij, k) - ux(ijlw, k)) * amhmod(ij) *             &
    &            rx / (hxu(ij) + hxu(ijlw)) * 2.d0                    &
    &          + (vx(ij, k) + vx(ijlw, k)) * amhmod(ij) *             &
    &            (hxyu(ij) + hxyu(ijlw)) * 0.25d0                     &
    &          - (  vx(ijln, k) + vx(ijnw, k)                         &
    &             - vx(ijls, k) - vx(ijsw, k)) * amhmod(ij)           &
    &            / (hyu(ij) + hyu(ijlw)) / (dy(ij) + dy(ijln))        &
    &          - (ux(ij, k) + ux(ijlw, k)) * amhmod(ij) *             & 
    &            (hyxu(ij) + hyxu(ijlw)) * 0.25d0
       sxy(ij) = (ux(ij, k) - ux(ijls, k)) * amhmod(ij) *             &
    &            ry(ij) / (hyu(ij) + hyu(ijls)) * 2.d0                &
    &          - (ux(ij, k) + ux(ijls, k)) * amhmod(ij) *             &
    &           (hxyu(ij) + hxyu(ijls)) * 0.25d0                      &
    &          + (  vx(ijle, k) + vx(ijse, k)                         &
    &             - vx(ijlw, k) - vx(ijsw, k)) * amhmod(ij)           &
    &            / (hxu(ij) + hxu(ijls)) * rx * 0.5d0                 &
    &          - (vx(ij, k) + vx(ijls, k)) * amhmod(ij) *             & 
    &            (hyxu(ij) + hyxu(ijls)) * 0.25d0
       syy(ij) = (vx(ij, k) - vx(ijls, k)) * amhmod(ij) *             &
    &            ry(ij) / (hyu(ij) + hyu(ijls)) * 2.d0                &
    &          + (ux(ij, k) + ux(ijls, k)) * amhmod(ij) *             &
    &            (hyxu(ij) + hyxu(ijls)) * 0.25d0                     &
    &          - (  ux(ijle, k) + ux(ijse, k)                         &
    &             - ux(ijlw, k) - ux(ijsw, k)) * amhmod(ij)           &
    &            / (hxu(ij) + hxu(ijls)) * rx * 0.5d0                 &
    &          - (vx(ij, k) + vx(ijls, k)) * amhmod(ij) *             &
    &            (hxyu(ij) + hxyu(ijls)) * 0.25d0                     
       syx(ij) = (vx(ij, k) - vx(ijlw, k)) * amhmod(ij) *             &
    &            rx / (hxu(ij) + hxu(ijlw)) * 2.d0                    &
    &          - (vx(ij, k) + vx(ijlw, k)) * amhmod(ij) *             &
    &            (hyxu(ij) + hyxu(ijlw)) * 0.25d0                     &
    &          + (  ux(ijln, k) + ux(ijnw, k)                         &
    &             - ux(ijls, k) - ux(ijsw, k)) * amhmod(ij)           &
    &            / (hyu(ij) + hyu(ijlw)) / (dy(ij) + dy(ijln))        &
    &          - (ux(ij, k) + ux(ijlw, k)) * amhmod(ij) *             &
    &            (hxyu(ij) + hxyu(ijlw)) * 0.25d0
       szx(ij) = - ux(ij, k) * amv(ij, k) / rea
       szy(ij) = - vx(ij, k) * amv(ij, k) / rea
    end do
     
    do ij = ijvstr, ijvend+1
       ijlw = ij + lw
       fux(ij) = sxx(ij) * (hyu(ij) + hyu(ijlw)) *                    &
    &                      (hyu(ij) + hyu(ijlw)) * 0.25d0 *           &
    &            amfvx(ij, k)
       fvx(ij) = syx(ij) * (hyu(ij) + hyu(ijlw)) *                    &
    &                      (hyu(ij) + hyu(ijlw)) * 0.25d0 *           &
    &            amfvx(ij, k)
    end do
    do ij = ijvstr, ijvend+nxdim
       ijls = ij + ls
       fuy(ij) = sxy(ij) * (hxu(ij) + hxu(ijls)) *                    &
    &                      (hxu(ij) + hxu(ijls)) * 0.25d0 *           &
    &            amfvy(ij, k) 
       fvy(ij) = syy(ij) * (hxu(ij) + hxu(ijls)) *                    &
    &                      (hxu(ij) + hxu(ijls)) * 0.25d0 *           &
    &            amfvy(ij, k)
    end do
    
    do ij = ijvstr, ijvend
       smx(ij) = szx(ij) / rea
       smy(ij) = szy(ij) / rea
    end do

    do ij = ijvstr, ijvend
       gx(ij, k) = (  gx(ij, k)                                      &
    &               + fuz(ij) * rz(ij)                               &
    &               + (  (fux(ij+le) - fux(ij)) *                    &
    &                    rx * ryu(ij)                                &
    &                  + (fuy(ij+ln) - fuy(ij)) *                    &
    &                    rym(ij) * rxu(ij)                           &
    &                 ) * rxu(ij) * ryu(ij) * rz(ij)                 &  
    &               + smx(ij)) * amskvb(ij)
       gy(ij, k) = (  gy(ij, k)                                      &
    &               + fvz(ij) * rz(ij)                               &
    &               + (  (fvx(ij+le) - fvx(ij)) *                    &
    &                    rx * ryu(ij)                                &
    &                  + (fvy(ij+ln) - fvy(ij)) *                    &
    &                    rym(ij) * rxu(ij)                           &
    &                 ) * rxu(ij) * ryu(ij) * rz(ij)                 &
    &               + smy(ij)) * amskvb(ij)
    end do

    do ij = ijvstr, ijvend
       xx(ij, kend) = gx(ij, kend)
       yy(ij, kend) = gy(ij, kend)
    end do
    
  end subroutine vscvlb

#endif
  
end module cvisc

