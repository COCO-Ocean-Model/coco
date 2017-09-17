module tflxt

! --- information -----------------------------------------------------
!  Second order moment scheme for tracer advection of Prather (1988, JGR)
!
!  HISTORY
!     '07.07.10  H.Tatebe
!     '07.10.24  H.Hasumi: for COCO4.3
!     '08.06.11  H.Hasumi: initial/final processing
!     '08.07.10  H.Hasumi: initial/final processing
!     '09.04.02  H.Tatebe: bug fix
!     '09.04.08  H.Tatebe: bug fix 2
!     '10.04.14  M.Kurogi: (COCO4.4 tripolar code by Dr. Suzuki)
!     '12.08.02  Y.Komuro: for COCO5.0
!     '13.02.13  Y.Komuro: remove non-parallel code 
!     '13.09.24  s.urakawa: bug fix (overshoot limiter)
!     '15.03.25  M.Kurogi: scalar tuning
!
! ---------------------------------------------------------------------

  use zocdim, only: &
    &  nxdim,  nydim,  nzdim,  ntdim, nxydim, nxyzdm, nztdim, &
    & nxgdim, nygdim, &
    &   kstr,   kend,     kz, &
    &     nx,     ny,     nz,    nxg,    nyg, &
    & ijtstr, ijtend, &
    &  igstr,  jgstr, &
    &     le,     lw,     ln,     ls, &
    &    lsw, &
    &  oinit, ofinal, &
    & myrank,  iroot
  use zocgrd, only: &
    &     dy,    dym,     dz,    dz0,    dzm,    dzv,     ds,    dsm, &
    &     dx,     rx,     ry,    rym, &
    &     ts,   zbot,    cor, &
    &    hxt,    hxu,    hyt,    hyu,    rxt,    ryt
  use zocmsk, only: &
#ifdef OPT_BBL
    & amsktb, &
#endif
    &  amskt,  amftx,  amfty,  amftz, &
    &   nbot
  use zocfil, only: &
    &    ncf
  use zocphy, only: &
    & gravit,   rhoo

  use brstt

  implicit none
  private

  real(8), save :: ftx(nxydim, nzdim, ntdim)
  real(8), save :: fty(nxydim, nzdim, ntdim)
  real(8), save :: ftz(nxydim, nzdim, ntdim)

  public :: flxtrc, chkftx
#ifdef OPT_BBL
  public :: flxtrb
#endif

contains 

subroutine flxtrc( &
  &    adt,  diffz, &
  &     tx,     hx,     ty,     hz, &
  &     uy,     vy,      w,    ahv )

  use bstbc
  use ufile
#ifdef OPT_IO_COCOMPI
  use mpiio
#else
  use bgs3d
#endif
  use qckot
  use bshft
  implicit none
#include "mpif.h"

  real(8), intent(out)    ::    adt(nxydim, nzdim, ntdim)    
  real(8), intent(out)    ::  diffz(nxydim, nzdim)
  real(8), intent(inout)  ::     tx(nxydim, nzdim, ntdim)
  real(8), intent(in)     ::     ty(nxydim, nzdim, ntdim)
  real(8), intent(in)     ::     hx(nxydim),     hz(nxydim)
  real(8), intent(in)     ::     uy(nxydim, nzdim),     vy(nxydim, nzdim)
  real(8), intent(in)     ::      w(nxydim, nzdim),    ahv(nxydim, nzdim)

  real(8) ::    wzc(nxydim, nzdim),    rzm(nxydim, nzdim)
  real(8) ::  hzbot(nxydim)

  real(8) ::  xdzdx(nxydim, nzdim),  ydzdy(nxydim, nzdim)
  real(8) ::  zdzdx(nxydim, nzdim),  zdzdy(nxydim, nzdim)
  real(8) ::  xdtdz(nxydim, nzdim, ntdim),  ydtdz(nxydim, nzdim, ntdim)
  real(8) ::  zdtdx(nxydim, nzdim, ntdim),  zdtdy(nxydim, nzdim, ntdim)

! ---- spatially varying isopycnal diffusion coefficient
  real(8), save ::   ahh3d(nxydim, nzdim),  ahi3d(nxydim, nzdim)
  real(8), save ::   ahg3d(nxydim, nzdim) 

!---- for second order moment
!---- bug fix (save these variables)
  real(8), save ::  s0 (nxydim, nzdim, ntdim)=0.d0
  real(8), save ::  sm (nxydim, nzdim, ntdim)=0.d0
  real(8), save ::  sx (nxydim, nzdim, ntdim), sxx(nxydim, nzdim, ntdim)
  real(8), save ::  sy (nxydim, nzdim, ntdim), syy(nxydim, nzdim, ntdim)
  real(8), save ::  sz (nxydim, nzdim, ntdim), szz(nxydim, nzdim, ntdim)
  real(8), save ::  sxy(nxydim, nzdim, ntdim), sxz(nxydim, nzdim, ntdim)
  real(8), save ::  syz(nxydim, nzdim, ntdim)   

  real(8), save :: alf(nxydim)=0.d0
  real(8), save ::  f0(nxydim)=0.d0
  real(8), save ::  fm(nxydim)=0.d0
  real(8), save ::  fx(nxydim)=0.d0, fxx(nxydim)=0.d0
  real(8), save ::  fy(nxydim)=0.d0, fyy(nxydim)=0.d0
  real(8), save ::  fz(nxydim)=0.d0, fzz(nxydim)=0.d0
  real(8), save :: fxy(nxydim)=0.d0, fxz(nxydim)=0.d0
  real(8), save :: fyz(nxydim)=0.d0

  integer :: ij1, ij2
!                             L2(kb), dble,  , 3D var num
  integer, parameter :: iblock=512000/ 8 /nz / 30

  real(8), save :: zalf(iblock,nzdim)=0.d0
  real(8), save ::  zf0(iblock,nzdim)=0.d0
  real(8), save ::  zfm(iblock,nzdim)=0.d0
  real(8), save ::  zfx(iblock,nzdim)=0.d0, zfxx(iblock,nzdim)=0.d0
  real(8), save ::  zfy(iblock,nzdim)=0.d0, zfyy(iblock,nzdim)=0.d0
  real(8), save ::  zfz(iblock,nzdim)=0.d0, zfzz(iblock,nzdim)=0.d0
  real(8), save :: zfxy(iblock,nzdim)=0.d0, zfxz(iblock,nzdim)=0.d0
  real(8), save :: zfyz(iblock,nzdim)=0.d0

  real(8), save ::  vlmx(nxydim, nzdim)=0.d0, vlmy(nxydim, nzdim)=0.d0
  real(8), save ::  vlmz(nxydim)=0.d0
  real(8), save ::  r(nxyzdm)
  real(8), save ::  uv(nxydim, nzdim)

  real(8) ::  s0m,    s1m,    s0p,    sxp
  real(8) ::  alfq,   alf1,   alf1q
  real(8) ::  u,      v,      tmp

  real(8), save ::  eps,   sq3,   ci3,   tsiv

  integer ::     ij,     k,      n,      i,     j
  integer ::   ijlw,   ijlsw,  ijle
  integer ::   ijln,   ijls
  integer ::    kuu,     ku,     kd
  integer ::  ifpar,  jfpar,  istat

  logical, save :: ofirst = .true.,   oeof

! for mixed layer eddy parameterization
  real(8) :: psigmx(nxydim, nzdim), psigmy(nxydim, nzdim)
  real(8) ::  xpsiy(nxydim, nzdim),  ypsix(nxydim, nzdim)
  real(8) ::  zpsix(nxydim, nzdim),  zpsiy(nxydim, nzdim)
  real(8) ::   igsy(nxydim, nzdim),   igsx(nxydim, nzdim)

  real(8), save ::  ahh = 0.0d0,  ahi = 0.0d0,  ahg = 0.0d0

  namelist /nmdifh/ ahh
  namelist /nmdifi/ ahi
  namelist /nmdifg/ ahg

!---- 
#ifdef OPT_IO_COCOMPI
 integer :: mpi_fh
 integer :: icread
 integer (kind = mpi_offset_kind) :: disp
#else
  real(8) ::  buf3(nxg, nyg, nz)
  real(8) ::  g3d(nxgdim, nygdim, nzdim)
#endif

!---- file name of isotropic diffusion and thickness diffusion
  character(len=ncf) ::  cfahi = 'not-specified'
  character(len=ncf) ::  cfahg = 'not-specified'
  character(len= 16) ::  chead(64) 
  integer :: iah = 0
  integer :: nfahi, nfahg

  namelist /nmcah/ cfahi, cfahg, iah

#ifdef OPT_BBL
  real(8), save :: ahhbbl = 0.0d0

  namelist /nmbbdh/ ahhbbl
#endif
  real(8) :: ss, d1, d2
  integer :: l, ll(nxydim, nzdim)

  if (oinit) then
     do n = 1, ntdim
#ifdef OPT_TRIPOLE
        call rstadd(sx(1, 1, n), oeof, &
          &         nxdim, nydim, nzdim, 'SX', 'OCN', &
          &                            1.d0,  0,  0 )
        call rstadd(sy(1, 1, n), oeof, &
          &         nxdim, nydim, nzdim, 'SY', 'OCN', &
          &                           -1.d0,  0,  0 )
        call rstadd(sz(1, 1, n), oeof, &
          &         nxdim, nydim, nzdim, 'SZ', 'OCN', &
          &                            1.d0,  0,  0 )
        call rstadd(sxx(1, 1, n), oeof, &
          &         nxdim, nydim, nzdim, 'SXX', 'OCN', &
          &                             1.d0,  0,  0 )
        call rstadd(syy(1, 1, n), oeof, &
          &         nxdim, nydim, nzdim, 'SYY', 'OCN', &
          &                             1.d0,  0,  0 )
        call rstadd(szz(1, 1, n), oeof, &
          &         nxdim, nydim, nzdim, 'SZZ', 'OCN', &
          &                             1.d0,  0,  0 )
        call rstadd(sxy(1, 1, n), oeof, &
          &         nxdim, nydim, nzdim, 'SXY', 'OCN', &
          &                            -1.d0,  0,  0 )
        call rstadd(sxz(1, 1, n), oeof, &
          &         nxdim, nydim, nzdim, 'SXZ', 'OCN', &
          &                             1.d0,  0,  0 )
        call rstadd(syz(1, 1, n), oeof, &
          &         nxdim, nydim, nzdim, 'SYZ', 'OCN', &
          &                            -1.d0,  0,  0 )
#else
        call rstadd(sx(1, 1, n), oeof, &
          &         nxdim, nydim, nzdim, 'SX', 'OCN')
        call rstadd(sy(1, 1, n), oeof, &
          &         nxdim, nydim, nzdim, 'SY', 'OCN')
        call rstadd(sz(1, 1, n), oeof, &
          &         nxdim, nydim, nzdim, 'SZ', 'OCN')
        call rstadd(sxx(1, 1, n), oeof, &
          &         nxdim, nydim, nzdim, 'SXX', 'OCN')
        call rstadd(syy(1, 1, n), oeof, &
          &         nxdim, nydim, nzdim, 'SYY', 'OCN')
        call rstadd(szz(1, 1, n), oeof, &
          &         nxdim, nydim, nzdim, 'SZZ', 'OCN')
        call rstadd(sxy(1, 1, n), oeof, &
          &         nxdim, nydim, nzdim, 'SXY', 'OCN')
        call rstadd(sxz(1, 1, n), oeof, &
          &         nxdim, nydim, nzdim, 'SXZ', 'OCN')
        call rstadd(syz(1, 1, n), oeof, &
          &         nxdim, nydim, nzdim, 'SYZ', 'OCN')
#endif
     end do
     return
  end if

  if (ofinal) then
     do n = 1, ntdim
        call finadd( sx(1, 1, n), nxdim, nydim, nzdim, 'SX',  'OCN')
        call finadd( sy(1, 1, n), nxdim, nydim, nzdim, 'SY',  'OCN')
        call finadd( sz(1, 1, n), nxdim, nydim, nzdim, 'SZ',  'OCN')
        call finadd(sxx(1, 1, n), nxdim, nydim, nzdim, 'SXX', 'OCN')
        call finadd(syy(1, 1, n), nxdim, nydim, nzdim, 'SYY', 'OCN')
        call finadd(szz(1, 1, n), nxdim, nydim, nzdim, 'SZZ', 'OCN')
        call finadd(sxy(1, 1, n), nxdim, nydim, nzdim, 'SXY', 'OCN')
        call finadd(sxz(1, 1, n), nxdim, nydim, nzdim, 'SXZ', 'OCN')
        call finadd(syz(1, 1, n), nxdim, nydim, nzdim, 'SYZ', 'OCN')
     end do
     return
  end if

  if (ofirst) then
     ofirst = .false.
     call rewnml(ifpar, jfpar)
     read(ifpar, nmdifh, iostat=istat)
     call cstnml(jfpar, 'flxtrc', 'nmdifh', istat)
     write(jfpar, nmdifh)
     call rewnml(ifpar, jfpar)
     read(ifpar, nmdifi, iostat=istat)
     call cstnml(jfpar, 'flxtrc', 'nmdifi', istat)
     write(jfpar, nmdifi)
     call rewnml(ifpar, jfpar)
     read(ifpar, nmdifg, iostat=istat)
     call cstnml(jfpar, 'flxtrc', 'nmdifg', istat)
     write(jfpar, nmdifg)

     eps = 1.d-20
     sq3 = sqrt( 3.d0 )
     ci3 = 1.d0 / 3.d0

!----- initialization
     if (oeof) then
        do n = 1, ntdim
           do k = 1, nzdim
              do ij = 1, nxydim
                 sx (ij, k, n) = 0.d0
                 sy (ij, k, n) = 0.d0
                 sz (ij, k, n) = 0.d0
                 sxx(ij, k, n) = 0.d0
                 syy(ij, k, n) = 0.d0
                 szz(ij, k, n) = 0.d0
                 sxy(ij, k, n) = 0.d0
                 sxz(ij, k, n) = 0.d0
                 syz(ij, k, n) = 0.d0
              end do
           end do
        end do
     end if

!    ---- area normal to u defined on UV-grid
     do k = kstr, kstr+kz-1
        do ij = 1, nxydim

           vlmx(ij, k) =  hxu(ij) * zbot * 0.5d0 * ds(k) * dx
           vlmy(ij, k) =  hyu(ij) * zbot * 0.5d0 * ds(k)
           
        end do
     end do

     do k = kstr+kz, kend
        do ij = 1, nxydim
            
           vlmx(ij, k) =  hxu(ij) * dzv(ij, k) * 0.5d0 * dx
           vlmy(ij, k) =  hyu(ij) * dzv(ij, k) * 0.5d0

        end do
     end do

     do ij = 1, nxydim
        vlmx(ij, kstr-1) = vlmx(ij, kstr)
        vlmx(ij, kend+1) = vlmx(ij, kend)
        vlmy(ij, kstr-1) = vlmy(ij, kstr)
        vlmy(ij, kend+1) = vlmy(ij, kend)
     end do

     do ij = 1, nxydim
         
        vlmz(ij) = hxt(ij) * dx * hyt(ij) * dy(ij)
        
     end do

!---- 
     call rewnml(ifpar, jfpar)
     read(ifpar, nmcah, iostat=istat)
     call cstnml(jfpar, 'flxtrc', 'nmcah', istat)
     write(jfpar, nmcah)

     if ( iah .eq. 0 ) then

        write(jfpar, *) 'Background horizontal diffusion :', ahh
        write(jfpar, *) 'Isopycnal diffusion             :', ahi
        write(jfpar, *) 'G-M thickness diffusion         :', ahg

        do k = 1, nzdim
           do ij = 1, nxydim
              ahi3d(ij, k) = ahi
              ahg3d(ij, k) = ahg
              ahh3d(ij, k) = ahh
           end do
        end do

     else

        do k = 1, nzdim
           do ij = 1, nxydim
              ahh3d(ij, k) = ahh
           end do
        end do

        write(jfpar, *) '  file name of ahi: ', cfahi
        write(jfpar, *) '  file name of ahg: ', cfahg

!       ---- reading diffusion coefficient
#ifdef OPT_IO_COCOMPI
!---- ahi                                                                                                                                        
        call mpi_filopn(mpi_fh, cfahi, 'READ')
        disp=0
        call mpi_read_chead(chead, mpi_fh, disp, icread)
        call mpi_read_3d(ahi3d, mpi_fh  , disp)
        call mpi_filcls(mpi_fh)

!---- ahg                                                                                                                                        
        call mpi_filopn(mpi_fh, cfahg, 'READ')
        disp=0
        call mpi_read_chead(chead, mpi_fh, disp, icread)
        call mpi_read_3d(ahg3d, mpi_fh  , disp)
        call mpi_filcls(mpi_fh)
#else
!       ---- ahi
        if ( myrank .eq. iroot ) then

           call filopn( nfahi, cfahi, 'READ' )
           rewind( nfahi )
           read( nfahi ) chead
           read( nfahi ) buf3

           do k = 1, nz
              do j = 1, nyg
                 do i = 1, nxg

                    g3d(igstr+i-1, jgstr+j-1, kstr+k-1) = buf3(i, j, k )

                 end do
              end do
           end do
           call filcls( nfahi )

        end if
        call scatter_3d( ahi3d, g3d )

!       ---- ahg
        if ( myrank .eq. iroot ) then

           call filopn( nfahg, cfahg, 'READ' )
           rewind( nfahg )
           read( nfahg ) chead
           read( nfahg ) buf3

           do k = 1, nz
              do j = 1, nyg
                 do i = 1, nxg

                    g3d(igstr+i-1, jgstr+j-1, kstr+k-1) = buf3(i, j, k )

                 end do
              end do
           end do
           call filcls( nfahg )

        end if
        call scatter_3d( ahg3d, g3d )
#endif
     end if

#ifdef OPT_BBL
     call rewnml(ifpar, jfpar)
     read(ifpar, nmbbdh, iostat=istat)
     call cstnml(jfpar, 'flxtrc', 'nmbbdh', istat)
     write(jfpar, nmbbdh)
#endif

  end if

  tsiv = 1.d0 / ts

  call dnsgrd( &
     &  xdzdx,  ydzdy,  zdzdx,  zdzdy, &
     &  xdtdz,  ydtdz,  zdtdx,  zdtdy, &
     &  xpsiy,  ypsix, &
     &  zpsix,  zpsiy, &
     &     ty,     tx,     hz )

  do n = 1, ntdim
     do k = 1, nzdim
        do ij = 1, nxydim
           adt (ij, k, n) = 0.d0
           ftx (ij, k, n) = 0.d0
           fty (ij, k, n) = 0.d0
           ftz (ij, k, n) = 0.d0
        end do
     end do
  end do

  do k = 1, nzdim
     do ij = 1, nxydim
        diffz(ij, k) = 0.d0
     end do
  end do

  do ij = 1, nxydim
     hzbot(ij) = hz(ij) + zbot
  end do

! ---- vertical velocity on sigma coordinate
  do k = kstr, kstr+kz-1
     do ij = 1, nxydim
        wzc(ij, k) = w(ij, k) * hzbot(ij)
        rzm(ij, k) = 1.d0 / dsm(k) / hzbot(ij)
     end do
  end do
  do k = kstr+kz, kend
     do ij = 1, nxydim
        wzc(ij, k) = w(ij, k)
        rzm(ij, k) = 1.d0 / dzm(ij, k)
     end do
  end do

  call chekin(wzc, 'WZC', nx, ny, nz, nxyzdm, 'OCN')

! ======  GM  isopycnal and diapycnal diffusion  ======
! ---- z diffusion flux of GM
  do n = 1, ntdim
     do k = kstr+1, kend
        kuu = k - 2
        ku  = k - 1
        kd  = k + 1
        do ij = ijtstr, ijtend

           diffz(ij, k) = &
             & (  ahv(ij, k) &
             &  + ahi3d(ij, k) * (  zdzdx(ij, k) * zdzdx(ij, k) &
             &                    + zdzdy(ij, k) * zdzdy(ij, k) ) ) * &
             & rzm(ij, k) * amftz(ij, k)

           ftz(ij, k, n) =  &
             & (  diffz(ij, k) * (tx(ij, ku, n) - tx(ij, k, n)) &
             &  - ( ( ahi3d(ij, k) + ahg3d(ij, k) ) * &
             &     zdzdx(ij, k) - zpsiy(ij, k) ) * zdtdx(ij, k, n) &
             &  - ( ( ahi3d(ij, k) + ahg3d(ij, k) ) * &
             &     zdzdy(ij, k) + zpsix(ij, k) ) * zdtdy(ij, k, n) &
             &  ) * amftz(ij, k)

        end do
     end do
  end do

! ---- y diffusion flux of GM
  do k = kstr, kend
  do n = 1, ntdim
        do ij = ijtstr, ijtend+nxdim

           ijls = ij + ls

           fty(ij, k, n) = &
             &     (  ( ahh3d(ij, k) + ahi3d(ij, k) ) * rym(ijls) &
             &      / ( hyt(ij) + hyt(ijls) ) * &
             &        ( tx(ij, k, n) - tx(ijls, k, n) ) * 2.d0 &
             &      - ( ( ahi3d(ij, k) - ahg3d(ij, k) ) &
             &        * ydzdy(ij, k) - ypsix(ij, k) ) * ydtdz(ij, k, n) ) &
             &      * ( hxu(ijls) + hxu(ij+lsw) ) * 0.5d0 * amfty(ij, k)

        end do

! ---- x diffusion flux of GM

        do ij = ijtstr, ijtend+1

           ijlw = ij + lw
           
           ftx(ij, k, n) = &
             &     (  ( ahh3d(ij, k) + ahi3d(ij, k) ) * rx &
             &      / ( hxt(ij) + hxt(ijlw) ) &
             &      * ( tx(ij, k, n) - tx(ijlw, k, n) ) * 2.d0 &
             &      - ( ( ahi3d(ij, k) - ahg3d(ij, k) ) &
             &      * xdzdx(ij, k) + xpsiy(ij, k) ) * xdtdz(ij, k, n) ) &
             &      * ( hyu(ijlw) + hyu(ij+lsw) ) * 0.5d0 * amftx(ij, k)

        end do

! ---- divergence of diffusion fluxes

        do ij = ijtstr, ijtend

           adt(ij, k, n) = &
             & (  (  (ftx(ij+le, k, n) - ftx(ij, k, n)) * rx &
             &     + (fty(ij+ln, k, n) - fty(ij, k, n)) * ry(ij)) * &
             &    rxt(ij) * ryt(ij) &
             &  + ftz(ij, k, n) - ftz(ij, k+1, n)) / dz(ij, k)

        end do
  end do
  end do

!---- diffusion in BBL
#ifdef OPT_BBL

  do n = 1, ntdim
     do ij = ijtstr, ijtend
        k = nbot(ij)
        ftz(ij, kend, n) = ftz(ij, k, n)
     end do
  end do

  do n = 1, ntdim

     do ij = ijtstr, ijtend+nxdim
        ijls = ij + ls
        fty(ij, kend, n) = &
          &   ahhbbl * (tx(ij, kend, n) - tx(ijls, kend, n)) * &
          &   rym(ij) * &
          &   (hxu(ijls) + hxu(ij+lsw)) / (hyt(ij) + hyt(ijls)) * &
          &   amfty(ij, kend)
     end do
     
     do ij = ijtstr, ijtend+1
        ijlw = ij + lw
        ftx(ij, kend, n) = &
          &   ahhbbl * (tx(ij, kend, n) - tx(ijlw, kend, n)) * rx * &
          &   (hyu(ijlw) + hyu(ij+lsw)) / (hxt(ij) + hxt(ijlw)) * &
          &   amftx(ij, kend) 
     end do

  end do

#endif

  do k = kstr, kend
     do ij=1, nxydim
        igsx(ij, k) = ( ahi3d(ij, k) - ahg3d(ij, k) ) * xdzdx(ij, k)
        igsy(ij, k) = ( ahi3d(ij, k) - ahg3d(ij, k) ) * ydzdy(ij, k)
     end do
  end do
  call chekin(igsx, 'IGSX', nx, ny, nz, nxyzdm, 'OCN')
  call chekin(igsy, 'IGSY', nx, ny, nz, nxyzdm, 'OCN')

  do k = 1, nzdim
     do ij = 1, nxydim
        psigmx(ij, k) = 0.0d0
        psigmy(ij, k) = 0.0d0
     end do
  end do

  do k = kstr, kend
     do ij = ijtstr, ijtend+nxdim
        psigmx(ij, k) = ahg3d(ij, k) * ydzdy(ij, k)
        psigmy(ij, k) = - ahg3d(ij, k) * xdzdx(ij, k)
     end do
  end do

  call chekin(psigmx, 'PSIGMX', nx, ny, nz, nxyzdm, 'OCN')
  call chekin(psigmy, 'PSIGMY', nx, ny, nz, nxyzdm, 'OCN')

!---- SOM 
!---- mass contained in a tracer grid

  do k = kstr, kend
!---- in X-direction

     do ij = ijtstr-nxdim-1, ijtend+nxdim+1

        ijlw  = ij + lw
        ijlsw = ij + lsw
               
        uv(ij, k) = ( uy(ijlw , k) * vlmy(ijlw, k) &
          &         + uy(ijlsw, k) * vlmy(ijlsw, k) ) &
          &       * amskt(ij, k) * amskt(ijlw, k) * dy(ijlw)

              !if u>0 LL=-1, else LL=0
               ll(ij,k)= - nint(0.5d0 +dsign(0.5D0, uv(ij,k))) 
     end do

     do n = 1, ntdim
        if(k .le. kstr+kz-1) then
           do ij = 1, nxydim
              sm(ij, k, n) = vlmz(ij) * hzbot(ij) * ds(k)
           end do
        else
           do ij = 1, nxydim
              sm(ij, k, n) = vlmz(ij) * dz(ij, k)
           end do
        end if

        do ij = 1, nxydim
           s0(ij, k, n) = sm(ij, k, n) * tx(ij, k, n)               
        end do

        if(k == kend) then
           do ij = 1, nxydim
              sm(ij, kstr-1, n) = sm(ij, kstr, n)
              sm(ij, kend+1, n) = sm(ij, kend, n)
              s0(ij, kstr-1, n) = s0(ij, kstr, n)
              s0(ij, kend+1, n) = s0(ij, kend, n)
           end do
        end if

! ---- undershoot limiter (Method B) of Morales Maqueda and Holloway (2006)

        do ij = ijtstr-nxdim-1, ijtend+nxdim+1

           ijlw = ij + lw * nint( amskt(ij+lw, k) )
           ijle = ij + le * nint( amskt(ij+le, k) )


           s0m  = s0(ij, k, n) &
             &  - min( s0(ijlw, k, n) / sm(ijlw, k, n), &
             &         s0(ij  , k, n) / sm(ij  , k, n), &
             &         s0(ijle, k, n) / sm(ijle, k, n)  ) &
             &  * sm(ij, k, n)


           s1m  = sq3 * s0m
           sx(ij, k, n) = min( s1m, max( - s1m, sx(ij, k, n) ) )

           s0p = s0m * s0m
           sxp = sqrt( max( s0p - sx(ij, k, n) * sx(ij, k, n) * ci3, &
             &              eps ) )

           d1=max( abs( sx(ij, k, n) ) - s0m, sxx(ij, k, n) ) 
           d2=max( s0m - sxp, sxx(ij, k, n) ) 
           ss=0.5D0 +dsign(0.5d0, abs( sx(ij, k, n) ) - 1.5d0 * s0m)
           sxx(ij, k, n) = min( s0m + sxp, d1*(1.d0-ss) + d2*ss )
           sxy(ij, k, n) = min( s0m, max( - s0m, sxy(IJ, K, N) ) )
           sxz(ij, k, n) = min( s0m, max( - s0m, sxz(IJ, K, N) ) )

!---- overshoot limiter (Method B) of Morales Maqueda and Holloway (2006)

           s0m = - s0(ij, k, n) &
             &   + max( s0(ijlw, k, n) / sm(ijlw, k, n), &
             &          s0(ij  , k, n) / sm(ij  , k, n), &
             &          s0(ijle, k, n) / sm(ijle, k, n)  ) &
             &   * sm(ij, k, n)

           s1m  = sq3 * s0m
           sx(ij, k, n) = min( s1m, max( - s1m, sx(ij, k, n) ) )

           s0p = s0m * s0m
           sxp = sqrt( max( s0p - sx(ij, k, n) * sx(ij, k, n) * ci3, &
             &              eps ) )

!               d1=max( abs( sx(ij, k, n) ) - s0m, sxx(ij, k, n) ) 
!               d2=max( s0m - sxp, sxx(ij, k, n) )
!               ss=0.5d0 +dsign(0.5d0, abs( sx(ij, k, n) ) - 1.5d0 * s0m)
!               sxx(ij, k, n) = min( s0m + sxp, d1*(1.d0-ss) + d2*ss )
               d1 = max( -s0m-sxp, sxx(ij, k, n) )
!               d2 = max( -s0m-sxp, sxx(ij, k, n) )
               ss=0.5d0 +dsign(0.5d0, abs( sx(ij, k, n) ) - 1.5d0 * s0m)
               sxx(ij, k, n) = min( (s0m-abs(sx(ij, k, n)))*(1.d0-ss) &
     &                                + (-s0m+sxp)*ss,                &
     &                              d1 )

               sxy(ij, k, n) = min( s0m, max( - s0m, sxy(ij, k, n) ) )
               sxz(ij, k, n) = min( s0m, max( - s0m, sxz(ij, k, n) ) )
          
        end do

!---- calculating ALF and MASS between box (I-1,J,K) <---> (I,J,K)

        do ij = ijtstr-nxdim-1, ijtend+nxdim+1
           ijlw  = ij + lw

           l=ll(ij,k)
           ss=dble(-2*l -1)
               
           fm(ijlw) = ss*uv(ij, k) * ts
           alf(ij)  = fm(ijlw) / sm(ij+l, k, n)              
        end do


!---- calculating flux and moments between box (I-1,J,K) <---> (I,J,K)                  

        do ij = ijtstr-nxdim-1, ijtend+nxdim+1
           ijlw  = ij + lw

           l=ll(ij,k)
           ss=dble(-2*l -1)

           alfq  = alf(ij) * alf(ij)
           alf1  = 1.d0 - alf(ij)
           alf1q = alf1 * alf1

           f0(ijlw) =                    &
     &       alf(ij) * ( s0(ij+l, k, n)  &
     &             +ss* alf1 * (         &
     &                  sx(ij+l, k, n)   &                        
     &       +ss* ( alf1 - alf(ij) ) * sxx(ij+l, k, n) ) )


           fx(ijlw) =                &
     &       alfq * (  sx(ij+l, k, n) &
     &       +ss* 3.d0 * alf1 * sxx(ij+l, k, n) )

           fxx(ijlw) = alf(ij) * alfq * sxx(ij+l, k, n)

           fy(ijlw) = alf(ij)  * ( sy(ij+l, k, n) &
     &               +ss* alf1 * sxy(ij+l, k, n) )

           fz(ijlw) = alf(ij) * ( sz(ij+l, k, n)  &
     &              + ss* alf1 * sxz(ij+l, k, n) )


           fxy(ijlw) = alfq    * sxy(ij+l, k, n)
           fxz(ijlw) = alfq    * sxz(ij+l, k, n)
           fyy(ijlw) = alf(ij) * syy(ij+l, k, n)
           fzz(ijlw) = alf(ij) * szz(ij+l, k, n)
           fyz(ijlw) = alf(ij) * syz(ij+l, k, n)

           ftx (ij, k, n) = ftx(ij, k, n) &
             &            -ss*f0(ijlw) * tsiv * ry(ijlw)          
        end do

!---- calculating flux and moments between box (I-1,J,K) <---> (I,J,K)                  

        do ij = ijtstr-nxdim-1, ijtend+nxdim+1

           ijlw  = ij + lw
           l=ll(ij,k)
           ss=dble(-2*l -1)


           alfq  = alf(ij) * alf(ij)
           alf1  = 1.d0 - alf(ij)
           alf1q = alf1 * alf1

           sm (ij+l, k, n) = sm(ij+l, k, n) - fm(ijlw)
           s0 (ij+l, k, n) = s0(ij+l, k, n) - f0(ijlw)

           sx (ij+l, k, n) =                 &
     &         alf1q * ( sx(ij+l, k, n)      &
     &               - ss* 3.d0 * alf(ij) * sxx(ij+l, k, n) )

           sxx(ij+l, k, n) = alf1 * alf1q * sxx(ij+l, k, n)
           sy (ij+l, k, n) = sy (ij+l, k, n) - fy(ijlw)
           syy(ij+l, k, n) = syy(ij+l, k, n) - fyy(ijlw)
           sz (ij+l, k, n) = sz (ij+l, k, n) - fz(ijlw)
           szz(ij+l, k, n) = szz(ij+l, k, n) - fzz(ijlw)
           sxy(ij+l, k, n) = alf1q * sxy(ij+l, k, n)
           sxz(ij+l, k, n) = alf1q * sxz(ij+l, k, n)
           syz(ij+l, k, n) = syz(ij+l, k, n) - fyz(ijlw)
        end do

!---- put the temporary moments (fi) into appropriate neighboring boxes
        do ij = ijtstr-nxdim-1, ijtend+nxdim+1

           ijlw  = ij + lw
           l=-(ll(ij,k)+1)

           sm(ij+l, k, n) = sm(ij+l, k, n) + fm(ijlw)
           alf(ij)     = fm(ijlw) / sm(ij+l, k, n)

        end do

        do ij = ijtstr-nxdim-1, ijtend+nxdim+1

           ijlw  = ij + lw
           l=-(ll(ij,k)+1)
           ss=dble(-2*ll(ij,k) -1)


           alf1 = 1.d0 - alf(ij)

           tmp = ss* alf(ij) * s0(ij+l, k, n) &
     &         - ss* alf1 * f0(ijlw)

           s0(ij+l, k, n) = s0(ij+l, k, n) + f0(ijlw)


           sxx(ij+l, k, n) =                              &
     &          alf(ij) * alf(ij) * fxx(ijlw)             &
     &        + alf1 * alf1 * sxx(ij+l, k, n)             &
     &        + 5.d0 * (  alf(ij) * alf1 * ss*(           &
     &                   sx(ij+l, k, n)  - fx(ijlw)    )  &
     &                  -ss* ( alf1 - alf(ij) ) * tmp  )


           sx (ij+l, k, n) =                                &
     &          alf(ij) * fx(ijlw) + alf1 * sx(ij+l, k, n)  &
     &        + 3.d0 * tmp


           sxy(ij+l, k, n) =         &
     &          alf(ij) * fxy(ijlw)  &
     &     + alf1 * sxy(ij+l, k, n)  &
     &    + 3.d0 * ss* ( alf(ij) * sy(ij+l, k, n) - alf1 * fy(ijlw) )


           sxz(ij+l, k, n) =         &
     &        alf(ij)  * fxz(ijlw)   &
     &      + alf1 * sxz(ij+l, k, n) &
     &      + 3.d0 * ss* ( alf(ij) * sz(ij+l, k, n) - alf1 * fz(ijlw) )


           sy (ij+l, k, n) = sy (ij+l, k, n) + fy(ijlw)
           sz (ij+l, k, n) = sz (ij+l, k, n) + fz(ijlw)

           syy(ij+l, k, n) = syy(ij+l, k, n) + fyy(ijlw)
           szz(ij+l, k, n) = szz(ij+l, k, n) + fzz(ijlw)
           syz(ij+l, k, n) = syz(ij+l, k, n) + fyz(ijlw)
        end do
     end do ! n
!---- y-direction
     do ij = ijtstr-nxdim-1, ijtend+nxdim+1

        ijls  = ij + ls
        ijlsw = ij + lsw

        uv(ij, k) = ( vy(ijls , k) * vlmx(ijls, k)    &
     &              + vy(ijlsw, k) * vlmx(ijlsw, k) ) &
     &            * amskt(ij, k) * amskt(ijls, k) 

        ll(ij,k)= - nint(0.5d0 +dsign(0.5d0, uv(ij,k))) 
     end do

     do n = 1, ntdim
!---- undershoot limiter (Method B) of Morales Maqueda and Holloway (2006)

        do ij = ijtstr-nxdim-1, ijtend+nxdim+1

           ijls  = ij + ls * nint( amskt(ij+ls, k) )
           ijln  = ij + ln * nint( amskt(ij+ln, k) )
                
           s0m = s0(ij, k, n)                            &
     &         - min( s0(ijls, k, n) / sm(ijls, k, n),   &
     &                s0(ij  , k, n) / sm(ij  , k, n),   &
     &                s0(ijln, k, n) / sm(ijln, k, n)  ) &
     &         * sm(ij, k, n)

           s1m = sq3 * s0m
           sy(ij, k, n) = min( s1m, max( - s1m, sy(ij, k, n) ) )

           s0p = s0m * s0m
           sxp = sqrt( max( s0p - sy(ij, k, n) * sy(ij, k, n) * ci3,  &
     &                          eps ) )

           d1=max( abs( sy(ij, k, n) ) - s0m, syy(ij, k, n) ) 
           d2=max( s0m - sxp, syy(ij, k, n) )
           ss=0.5d0 +dsign(0.5d0, abs( sy(ij, k, n) ) - 1.5d0 * s0m)
           syy(ij, k, n) = min( s0m + sxp, d1*(1.d0-ss) + d2*ss )

             
           sxy(ij, k, n) = min( s0m, max( - s0m, sxy(ij, k, n) ) )
           syz(ij, k, n) = min( s0m, max( - s0m, syz(ij, k, n) ) )

!---- overshoot limiter (Method B) of Morales Maqueda and Holloway (2006)

           s0m = - s0(ij, k, n)                             &
     &           + max( s0(ijls, k, n) / sm(ijls, k, n),    &
     &                  s0(ij  , k, n) / sm(ij  , k, n),    &
     &                  s0(ijln, k, n) / sm(ijln, k, n)  )  &
     &           * sm(ij, k, n)

           s1m = sq3 * s0m
           sy(ij, k, n) = min( s1m, max( - s1m, sy(ij, k, n) ) )

           s0p = s0m * s0m
           sxp = sqrt( max( s0p - sy(ij, k, n) * sy(ij, k, n) * ci3,  &
     &                      eps ) )
               

!          d1=max( abs( sy(ij, k, n) ) - s0m, syy(ij, k, n) ) 
!          d2=max( s0m - sxp, syy(ij, k, n) )
!          ss=0.5d0 +dsign(0.5d0, abs( sy(ij, k, n) ) - 1.5d0 * s0m)
!          syy(ij, k, n) = min( s0m + sxp, d1*(1.d0-ss) + d2*ss )

           d1 = max( -s0m-sxp, syy(ij, k, n) )
!          d2 = max( -s0m-sxp, syy(ij, k, n) )
           ss=0.5d0 +dsign(0.5d0, abs( sy(ij, k, n) ) - 1.5d0 * s0m)
           syy(ij, k, n) = min( (s0m-abs(sy(ij, k, n)))*(1.d0-ss)    &
     &                            + (-s0m+sxp)*ss,                   &
     &                          d1 )
             
           sxy(ij, k, n) = min( s0m, max( - s0m, sxy(ij, k, n) ) )
           syz(ij, k, n) = min( s0m, max( - s0m, syz(ij, k, n) ) )

        end do
     end do ! n
  end do ! k

!---- bug fix 2
#ifdef OPT_TRIPOLE
     call shift1( sy, &
       &           nxdim, nydim, nztdim, &
       &           -1.d0,     0,      0 )
     call shift1( syy, &
       &           nxdim, nydim, nztdim, &
       &            1.d0,     0,      0 )
     call shift1( sxy, &
       &           nxdim, nydim, nztdim, &
       &           -1.d0,     0,      0 )
     call shift1( syz, &
       &           nxdim, nydim, nztdim, &
       &           -1.d0,     0,      0 )
#else
     call shift1( sy , nxdim, nydim, nztdim )
     call shift1( syy, nxdim, nydim, nztdim )
     call shift1( sxy, nxdim, nydim, nztdim )
     call shift1( syz, nxdim, nydim, nztdim )
#endif

  do k = kstr, kend
     do n=1, ntdim
!---- calculating ALF  and MASS between box (I,J-1,K) <---> (I,J,K)
        do ij = ijtstr-nxdim-1, ijtend+nxdim+1
           ijls  = ij + ls
           l=ll(ij,k) *(-ls)
           ss=dble(-2*ll(ij,k) -1)

           fm(ijls) = ss* uv(ij, k) * ts
           alf(ij)  = fm(ijls) / sm(ij+l, k, n)
        end do

!---- calculating flux between box (I,J-1,K) <---> (I,J,K)
        do ij = ijtstr-nxdim-1, ijtend+nxdim+1
           ijls  = ij + ls
           l=ll(ij,k) *(-ls)
           ss=dble(-2*ll(ij,k) -1)

           alfq  = alf(ij) * alf(ij)
           alf1  = 1.d0 - alf(ij)
           alf1q = alf1 * alf1


           f0(ijls) =                        &
     &         alf(ij) * ( s0(ij+l, k, n)    &
     &            + alf1 * ss* (             &
     &               sy(ij+l, k, n)          &
     &        + ss* ( alf1 - alf(ij) ) * syy(ij+l, k, n) ) )


           fy(ijls) =  alfq * (              &
     &         sy(ij+l, k, n)   + ss* 3.d0 * alf1 * syy(ij+l, k, n) ) 

           fyy(ijls) = alf(ij) * alfq * syy(ij+l, k, n)

           fx(ijls) = alf(ij) * ( sx(ij+l, k, n)  &
     &               + ss* alf1 * sxy(ij+l, k, n) )

           fz(ijls) = alf(ij) * ( sz(ij+l, k, n)  &
     &             + ss * alf1 * syz(ij+l, k, n) )


           fxy(ijls) = alfq       * sxy(ij+l, k, n)
           fyz(ijls) = alfq       * syz(ij+l, k, n)
           fxx(ijls) = alf(ij)    * sxx(ij+l, k, n)
           fzz(ijls) = alf(ij)    * szz(ij+l, k, n)
           fxz(ijls) = alf(ij)    * sxz(ij+l, k, n)

           fty (ij, k, n) = fty(ij, k, n) &
     &             - ss * f0(ijls) * tsiv * rx
        end do

!---- calculating moments
        do ij = ijtstr-nxdim-1, ijtend+nxdim+1

           ijls  = ij + ls
           l=ll(ij,k) *(-ls)
           ss=dble(-2*ll(ij,k) -1)

           alfq  = alf(ij) * alf(ij)
           alf1  = 1.d0 - alf(ij)
           alf1q = alf1 * alf1

           sm (ij+l, k, n) = sm(ij+l, k, n) - fm(ijls)
           s0 (ij+l, k, n) = s0(ij+l, k, n) - f0(ijls)

           sy (ij+l, k, n) = alf1q * ( &
     &              sy(ij+l, k, n)     &
     &             - ss* 3.d0 * alf(ij) * syy(ij+l, k, n) )

           syy(ij+l, k, n) = alf1 * alf1q * syy(ij+l, k, n)
           sx (ij+l, k, n) = sx (ij+l, k, n) -  fx(ijls)
           sxx(ij+l, k, n) = sxx(ij+l, k, n) - fxx(ijls)
           sz (ij+l, k, n) = sz (ij+l, k, n) -  fz(ijls)
           szz(ij+l, k, n) = szz(ij+l, k, n) - fzz(ijls)
           sxy(ij+l, k, n) = alf1q * sxy(ij+l, k, n)
           syz(ij+l, k, n) = alf1q * syz(ij+l, k, n)
           sxz(ij+l, k, n) = sxz(ij+l, k, n) - fxz(ijls)

        end do

!---- put the temporary moments (fi) into appropriate neighboring boxes

        do ij = ijtstr-nxdim-1, ijtend+nxdim+1

           ijls  = ij + ls
           l=(ll(ij,k)+1)* ls
           ss=dble(-2*ll(ij,k) -1)

           sm(ij+l, k, n) = sm(ij+l, k, n) + fm(ijls)
           alf(ij)     = fm(ijls) / sm(ij+l, k, n)

        end do

        do ij = ijtstr-nxdim-1, ijtend+nxdim+1

           ijls  = ij + ls
           l=(ll(ij,k)+1)* ls
           ss=dble(-2*ll(ij,k) -1)
           alf1 = 1.d0 - alf(ij)

           tmp =   ss*alf(ij) * s0(ij+l, k, n)   - ss*alf1 * f0(ijls)
           s0(ij+l, k, n) = s0(ij+l, k, n) + f0(ijls)

           syy(ij+l, k, n) =                 &
     &        alf(ij) * alf(ij) * fyy(ijls)  &
     &       + alf1 * alf1 * syy(ij+l, k, n) &
     &     + 5.d0 * (  alf(ij) * alf1 * ss*( & 
     &           sy(ij+l, k, n) - fy(ijls) ) &
              - ss*( alf1 - alf(ij) ) * tmp )


           sy (ij+l, k, n) = &
     &      alf(ij) * fy(ijls) + alf1 * sy(ij+l, k, n) + 3.d0 * tmp

           sxy(ij+l, k, n) =                         &
     &        alf(ij) * fxy(ijls)                    &
     &      + alf1 * sxy(ij+l, k, n)                 &
     &      + 3.d0 * ss* ( alf(ij) * sx(ij+l, k, n)  &
     &                      - alf1 * fx(ijls) )

           syz(ij+l, k, n) =                        &
     &        alf(ij) * fyz(ijls)                   & 
     &      + alf1 * syz(ij+l, k, n)                &
     &      + 3.d0 * ss*( alf(ij) * sz(ij+l, k, n)  &
     &                    - alf1 * fz(ijls) )

           sx (ij+l, k, n) = sx (ij+l, k, n) + fx(ijls)
           sz (ij+l, k, n) = sz (ij+l, k, n) + fz(ijls)
           sxx(ij+l, k, n) = sxx(ij+l, k, n) + fxx(ijls)
           szz(ij+l, k, n) = szz(ij+l, k, n) + fzz(ijls)
           sxz(ij+l, k, n) = sxz(ij+l, k, n) + fxz(ijls)
        end do
     end do ! n
  end do ! k

#ifdef OPT_BBL
!---- keeping BBL variables consistent at two levels
  call stbbt2( s0 )
  call stbbt2( sm )
  call stbbt2( sx )
  call stbbt2( sy )
  call stbbt2( sz )
  call stbbt2( sxx )
  call stbbt2( syy )
  call stbbt2( szz )
  call stbbt2( sxy )
  call stbbt2( sxz )
  call stbbt2( syz )
#endif

!---- Z-direction

  do k = kstr, kend
     do ij = ijtstr, ijtend

        uv(ij, k) = - wzc(ij, k) * vlmz(ij) &
     &       * amskt(ij, k) * amskt(ij, k-1)
               
     end do
  end do

#ifdef OPT_BBL
  do ij = ijtstr, ijtend
     k = nbot(ij)
     uv(ij, k) = - wzc(ij, k) * vlmz(ij) * amsktb(ij)
  end do
#endif

  do k = kstr, kend
     do ij = ijtstr, ijtend
          !if w>0 ll=-1, else l=0
           ll(ij,k)= - nint(0.5d0 +dsign(0.5d0, uv(ij,k))) 
     end do
  end do


  do n = 1, ntdim
     do ij1 = ijtstr, ijtend, iblock       !cashe blocking
        ij2 =min(ij1 + iblock-1, ijtend) 
!---- undershoot limiter (Method B) of Morales Maqueda and Holloway (2006)
        do k = kstr, kend
           ku = max( k - 1, kstr )
           kd = min( k + 1, kend )
            
           do ij = ij1, ij2
              s0m = s0(ij, k, n)                          &
     &            - min( s0(ij, ku, n) / sm(ij, ku, n),   &
     &                   s0(ij, k , n) / sm(ij, k , n),   &
     &                   s0(ij, kd, n) / sm(ij, kd, n)  ) &
     &            * sm(ij, k, n)

              s1m = sq3 * s0m
              sz(ij, k, n) = min( s1m, max( - s1m, sz(ij, k, n) ) )
              s0p = s0m * s0m
              sxp = sqrt( max( s0p - sz(ij, k, n) * sz(ij, k, n) * ci3,  &
     &                         eps ) )


              d1=max( abs( sz(ij, k, n) ) - s0m, szz(ij, k, n) )
              d2=max( s0m - sxp, szz(ij, k, n) )
              ss=0.5d0 +dsign(0.5d0, abs( sz(ij, k, n) ) - 1.5d0 * s0m)
              szz(ij, k, n) = min( s0m + sxp, d1*(1.d0-ss) + d2*ss )

             
              sxz(ij, k, n) = min( s0m, max( - s0m, sxz(ij, k, n) ) )
              syz(ij, k, n) = min( s0m, max( - s0m, syz(ij, k, n) ) )


!---- overshoot limiter (Method B) of Morales Maqueda and Holloway (2006)

              s0m = - s0(ij, k, n)                          &
     &              + max( s0(ij, ku, n) / sm(ij, ku, n),   &
     &                     s0(ij, k,  n) / sm(ij, k , n),   &
     &                     s0(ij, kd, n) / sm(ij, kd, n)  ) &
     &              * sm(ij, k, n)

              s1m = sq3 * s0m
              sz(ij, k, n) = min( s1m, max( - s1m, sz(ij, k, n) ) )
              s0p = s0m * s0m
              sxp = sqrt( max( s0p - sz(ij, k, n) * sz(ij, k, n) * ci3,  &
     &                         eps ) )
               

!             d1=abs( sz(ij, k, n) ) - s0m
!             d2=s0m - sxp
!             ss=0.5d0 +dsign(0.5d0, abs( sz(ij, k, n) ) - 1.5d0 * s0m)
!             szz(ij, k, n) = min( s0m + sxp, 
!     &               max( d1*(1.d0-ss) + d2*ss,  szz(ij, k, n) ) )

              d1 = max( -s0m-sxp, szz(ij, k, n) )
!             d2 = max( -s0m-sxp, szz(ij, k, n) )
              ss=0.5d0 +dsign(0.5d0, abs( sz(ij, k, n) ) - 1.5d0 * s0m)
              szz(ij, k, n) = min( (s0m-abs(sz(ij, k, n)))*(1.d0-ss) &
     &                               + (-s0m+sxp)*ss,  &
     &                             d1 )

              sxz(ij, k, n) = min( s0m, max( - s0m, sxz(ij, k, n) ) )
              syz(ij, k, n) = min( s0m, max( - s0m, syz(ij, k, n) ) )
           end do
        end do

!---- calculating ALF
        do k = kstr, kend
           ku = k - 1
           do ij = ij1, ij2

              l=ll(ij,k)
              ss=dble(-2*l -1)

              zfm(ij-ij1+1, ku) = ss*uv(ij, k) * ts
              zalf(ij-ij1+1, k) = zfm(ij-ij1+1, ku) / sm(ij, k+l, n)
           end do
        end do

!---- calculating flux between box (I,J,K-1) <---> (I,J,K)
        do k = kstr, kend
           ku = k - 1
           do ij = ij1, ij2
              alfq  = zalf(ij-ij1+1, k) * zalf(ij-ij1+1, k)
              alf1  = 1.d0 - zalf(ij-ij1+1, k)
              alf1q = alf1 * alf1

              l=ll(ij,k)
              ss=dble(-2*l -1)

              zf0(ij-ij1+1, ku) =                           &
     &          zalf(ij-ij1+1, k) * ( s0(ij, k+l, n)        &
     &               + ss*alf1 * (                          &
     &                     sz(ij, k+l, n)                   &
     &                   + ss*( alf1 - zalf(ij-ij1+1, k) )  &
     &                   * szz(ij, k+l, n) ) )

              zfz(ij-ij1+1, ku) =                           &
     &          alfq * ( sz(ij, k+l, n)                     &
     &                 + ss*3.d0 * alf1 * szz(ij, k+l, n) ) 

              zfzz(ij-ij1+1, ku) = zalf(ij-ij1+1, k)  &
     &                     * alfq * szz(ij, k+l, n)

              zfx(ij-ij1+1, ku) = zalf(ij-ij1+1, k)*( sx(ij, k+l, n)  &
     &                              + ss*alf1 * sxz(ij, k+l, n) )

              zfy(ij-ij1+1, ku) = zalf(ij-ij1+1, k)*( sy(ij, k+l, n) &
     &                                + ss*alf1 * syz(ij, k+l, n) )


              zfxz(ij-ij1+1, ku) = alfq            * sxz(ij, k+l, n)
              zfyz(ij-ij1+1, ku) = alfq            * syz(ij, k+l, n)
              zfxx(ij-ij1+1, ku) = zalf(ij-ij1+1, k)*sxx(ij, k+l, n)
              zfyy(ij-ij1+1, ku) = zalf(ij-ij1+1, k)*syy(ij, k+l, n)
              zfxy(ij-ij1+1, ku) = zalf(ij-ij1+1, k)*sxy(ij, k+l, n)
                  
              ftz (ij, k, n) = ftz(ij, k, n)  &
     &                   +ss* zf0(ij-ij1+1, ku) * tsiv / vlmz(ij)
           end do
        end do
!---- calculating flux between box (i,j,k-1) <---> (i,j,k)
        do k = kstr, kend
           ku = k - 1
           do ij = ij1, ij2
              l=ll(ij,k)
              ss=dble(-2*l -1)

              alfq  = zalf(ij-ij1+1, k) * zalf(ij-ij1+1, k)
              alf1  = 1.d0 - zalf(ij-ij1+1, k)
              alf1q = alf1 * alf1

              sm (ij, k+l, n) = sm(ij, k+l, n) - zfm(ij-ij1+1, ku)
              s0 (ij, k+l, n) = s0(ij, k+l, n) - zf0(ij-ij1+1, ku)

              sz (ij, k+l, n) =                                  &
     &          alf1q * ( sz(ij, k+l, n)                         &
     &         -ss* 3.d0 * zalf(ij-ij1+1, k) * szz(ij, k+l, n) )


              szz(ij, k+l, n) = alf1 * alf1q * szz(ij, k+l, n)
              sx (ij, k+l, n) = sx (ij, k+l, n) - zfx(ij-ij1+1, ku)
              sxx(ij, k+l, n) = sxx(ij, k+l, n) - zfxx(ij-ij1+1, ku)
              sy (ij, k+l, n) = sy (ij, k+l, n) - zfy(ij-ij1+1, ku)
              syy(ij, k+l, n) = syy(ij, k+l, n) - zfyy(ij-ij1+1, ku)

              sxz(ij, k+l, n) = alf1q * sxz(ij, k+l, n)
              syz(ij, k+l, n) = alf1q * syz(ij, k+l, n)
              sxy(ij, k+l, n) = sxy(ij, k+l, n) - zfxy(ij-ij1+1, ku)
          
           end do
        end do

!---- put the temporary moments (fi) into appropriate neighboring boxes
        do k = kstr, kend
           ku = k - 1
           do ij = ij1, ij2
              l=-(ll(ij,k)+1)
              sm(ij, k+l, n) = sm(ij, k+l, n) + zfm(ij-ij1+1, ku)
              zalf(ij-ij1+1, k)   = zfm(ij-ij1+1, ku) / sm(ij, k+l, n)
           end do
        end do

        do k = kstr, kend
           ku = k - 1
           do ij = ij1, ij2
              l=-(ll(ij,k)+1)
              ss=dble(-2*ll(ij,k) -1)
              alf1 = 1.d0 - zalf(ij-ij1+1, k)

              tmp = ss* zalf(ij-ij1+1, k) * s0(ij, k+l, n)  &
     &              -ss* alf1 * zf0(ij-ij1+1, ku)

              s0(ij, k+l, n) = s0(ij, k+l, n) + zf0(ij-ij1+1, ku)

              szz(ij, k+l, n) =                                           &
     &          zalf(ij-ij1+1, k) * zalf(ij-ij1+1, k)*zfzz(ij-ij1+1, ku)  &
     &             + alf1 * alf1 * szz(ij, k+l, n)                        &
     &             + 5.d0 * (  zalf(ij-ij1+1, k)                          &
     &             * alf1 * ss*(   sz(ij, k+l, n) - zfz(ij-ij1+1, ku) )   &
     &                     -ss* ( alf1 - zalf(ij-ij1+1, k) ) * tmp )


              sz (ij, k+l, n) =                            &
     &              zalf(ij-ij1+1, k) * zfz(ij-ij1+1, ku)  &
     &            + alf1 * sz(ij, k+l, n)                  &
     &            + 3.d0 * tmp

              sxz(ij, k+l, n) =                                     &
     &              zalf(ij-ij1+1, k) * zfxz(ij-ij1+1, ku)          &
     &            + alf1 * sxz(ij, k+l, n)                          &
     &            + 3.d0 * ss*( zalf(ij-ij1+1, k) * sx(ij, k+l, n)  &
     &                      - alf1 * zfx(ij-ij1+1, ku) )


              syz(ij, k+l, n) =                                     &
     &              zalf(ij-ij1+1, k) * zfyz(ij-ij1+1, ku)          &
     &            + alf1 * syz(ij, k+l, n)                          &
     &            + 3.d0 * ss*( zalf(ij-ij1+1, k) * sy(ij, k+l, n)  &
     &                      - alf1 * zfy(ij-ij1+1, ku) )

              sx (ij, k+l, n) = sx (ij, k+l, n) +  zfx(ij-ij1+1, ku)
              sy (ij, k+l, n) = sy (ij, k+l, n) +  zfy(ij-ij1+1, ku)
              sxx(ij, k+l, n) = sxx(ij, k+l, n) + zfxx(ij-ij1+1, ku)
              syy(ij, k+l, n) = syy(ij, k+l, n) + zfyy(ij-ij1+1, ku)
              sxy(ij, k+l, n) = sxy(ij, k+l, n) + zfxy(ij-ij1+1, ku)
           end do
        end do
     end do !blocking
  end do

!---- tx
  do n = 1, ntdim
     do k = kstr, kend
        do ij = ijtstr, ijtend

!           tx(ij, k, n) = s0(ij, k, n) / sm(ij, k, n)

               adt(ij, k, n) =                                       &
     &         (  (  (ftx(ij+le, k, n) - ftx(ij, k, n)) * rx         &
     &             + (fty(ij+ln, k, n) - fty(ij, k, n)) * ry(ij)) *  &
     &            rxt(ij) * ryt(ij)                                  &
     &          + ftz(ij, k, n) - ftz(ij, k+1, n)) / dz(ij, k)

        end do
     end do
  end do

#ifdef OPT_BBL
  do n = 1, ntdim
     do ij = ijtstr, ijtend
        k = nbot(ij)
        ftz(ij, kend, n) = ftz(ij, k, n)
     end do
  end do

  do n = 1, ntdim
     do ij = ijtstr, ijtend
        adt(ij, kend, n) = &
             & ( ( (ftx(ij+le, kend, n) - ftx(ij, kend, n)) * rx &
             &   + (fty(ij+ln, kend, n) - fty(ij, kend, n)) * ry(ij)) * &
             &  rxt(ij) * ryt(ij) &
             & + ftz(ij, kend, n) ) / dz(ij, kend)
     end do
  end do
#endif

#ifdef OPT_BBL
!---- keeping BBL variables consistent at two levels
  call stbbtr( s0 )
  call stbbtr( sm )
  call stbbtr( sx )
  call stbbtr( sy )
  call stbbtr( sz )
  call stbbtr( sxx )
  call stbbtr( syy )
  call stbbtr( szz )
  call stbbtr( sxy )
  call stbbtr( sxz )
  call stbbtr( syz )
#endif


#ifdef OPT_TRIPOLE
  call shift1(    sx, &
    &          nxdim, nydim, nztdim, &
    &           1.d0,     0,      0 )
  call shift1(    sy, &
    &          nxdim, nydim, nztdim, &
    &          -1.d0,     0,      0 )
  call shift1(    sz, &
    &          nxdim, nydim, nztdim, &
    &           1.d0,     0,      0 )
  call shift1(   sxx, &
    &          nxdim, nydim, nztdim, &
    &           1.d0,     0,      0 )
  call shift1(   syy, &
    &          nxdim, nydim, nztdim, &
    &           1.d0,     0,      0 )
  call shift1(   szz, &
    &          nxdim, nydim, nztdim, &
    &           1.d0,     0,      0 )
  call shift1(   sxy, &
    &          nxdim, nydim, nztdim, &
    &          -1.d0,     0,      0 )
  call shift1(   sxz, &
    &          nxdim, nydim, nztdim, &
    &           1.d0,     0,      0 )
  call shift1(   syz, &
    &          nxdim, nydim, nztdim, &
    &          -1.d0,     0,      0 )
#else
!  call shift1( s0 , nxdim, nydim, nztdim )
  call shift1( sx , nxdim, nydim, nztdim )
  call shift1( sy , nxdim, nydim, nztdim )
  call shift1( sz , nxdim, nydim, nztdim )
  call shift1( sxx, nxdim, nydim, nztdim )
  call shift1( syy, nxdim, nydim, nztdim )
  call shift1( szz, nxdim, nydim, nztdim )
  call shift1( sxy, nxdim, nydim, nztdim )
  call shift1( sxz, nxdim, nydim, nztdim )
  call shift1( syz, nxdim, nydim, nztdim )
#endif

  return

end subroutine flxtrc

! *********************************************************************

subroutine dnsgrd( &
  &  xdzdx,  ydzdy,  zdzdx,  zdzdy, &
  &  xdtdz,  ydtdz,  zdtdx,  zdtdy, &
  &  xpsiy,  ypsix, &
  &  zpsix,  zpsiy, &
  &     ty,     tx,     hz )

  use bshft
  use qckot
  use ufile
  use xprst

  real(8), intent(out) ::  xdzdx(nxydim, nzdim),  ydzdy(nxydim, nzdim)
  real(8), intent(out) ::  zdzdx(nxydim, nzdim),  zdzdy(nxydim, nzdim)
  real(8), intent(out) ::  xdtdz(nxydim, nzdim, ntdim)
  real(8), intent(out) ::  ydtdz(nxydim, nzdim, ntdim)
  real(8), intent(out) ::  zdtdx(nxydim, nzdim, ntdim)
  real(8), intent(out) ::  zdtdy(nxydim, nzdim, ntdim)
  real(8), intent(in)  ::     ty(nxydim, nzdim, ntdim)
  real(8), intent(in)  ::     tx(nxydim, nzdim, ntdim)

  real(8), save :: c0(nzdim), c1(nzdim), c2(nzdim)
  real(8), save :: c3(nzdim), c4(nzdim), c5(nzdim), c6(nzdim)
  real(8), save :: d0(nzdim), d1(nzdim), d2(nzdim), d3(nzdim), d4(nzdim)
  real(8), save :: d5(nzdim), d6(nzdim), d7(nzdim), d8(nzdim), d9(nzdim)
  real(8), save :: eps = 1.d-20
  logical, save :: ofirst = .true.

  real(8), save ::  cxpsy(nxydim), cypsx(nxydim) 
  real(8), save ::  czpsx(nxydim), czpsy(nxydim) 
  integer, save ::  kzmin

  real(8) ::      r(nxydim, nzdim)
  real(8) ::   hmld(nxydim), hmld1(nxydim)
  real(8) :: rmavez(nxydim), rmav1(nxydim)
  real(8) ::  dzsig(nxydim, nzdim), dzmsig(nxydim, nzdim)
  real(8) ::     zt(nxydim, nzdim),    ztm(nxydim, nzdim)
  real(8) :: rsigth(nxydim, nzdim)
  real(8) ::    nbv(nxydim),            lf(nxydim)
  real(8) :: xpsiy1(nxydim, nzdim), ypsix1(nxydim, nzdim)
  real(8) :: zpsix1(nxydim, nzdim), zpsiy1(nxydim, nzdim)
  real(8) ::  zmld0(nxydim),         zhmld(nxydim, nzdim)
  real(8) ::  hmldx(nxydim),         hmldy(nxydim)
  real(8) :: xpsiyz(nxydim, nzdim), ypsixz(nxydim, nzdim)
  integer ::   kmld(nxydim)

  real(8) :: rmavdx(nxydim), rmavdy(nxydim)
  real(8) ::   muzx(nxydim, nzdim),   muzy(nxydim, nzdim)

  real(8) ::  xpsiy(nxydim, nzdim),  ypsix(nxydim, nzdim)
  real(8) ::  zpsix(nxydim, nzdim),  zpsiy(nxydim, nzdim)
  real(8) ::     hz(nxydim)

  real(8) ::   dtdx(nxydim, nzdim, ntdim),   dtdy(nxydim, nzdim, ntdim)
  real(8) ::  dtfdz(nxydim, nzdim, ntdim)
  real(8) ::     p1,     p2
  real(8) ::     tl,     sl
  real(8) ::     rl,    rlw,    rls,    rlu
  real(8) ::   dzdx,   dzdy
  integer ::     ij,      k,      n
  integer ::  ifpar,  jfpar,  istat

  real(8) ::   muzh,  in2dz, n2min, n2l, hmldt
  real(8) :: rsigdf, rsigbt,  dhmld, cpsi
  real(8) ::     pi,  omega, cormin

  real(8), save :: slpmax = 1.d-2
  real(8), save ::  cm = 8.0d0,  ce = 0.06d0,  fminlt = 10.0d0
  real(8), save ::  lfmin = 1.0d5,  taumle = 10.0d0,  vscl = 50.0d0
  real(8), save ::  drsig = 0.1d0
  integer, save ::  mz = nz,  mzmin = 1
  integer, save ::  nfltdm = 0,  nfltps = 0,  nfltrm = 0
  logical, save ::  ofltdm = .false.,  ofltps = .false.
  logical, save ::  ofltrm = .false.,  ocoamp = .true.

  namelist /nmslpm/ slpmax
  namelist /nmmlep/ cm, ce, fminlt, mz, ofltdm, nfltdm, &
    &               lfmin, taumle, vscl, ofltps, nfltps, mzmin, &
    &               drsig, ocoamp, ofltrm, nfltrm

  if (ofirst) then
     ofirst = .false.
     call rewnml(ifpar, jfpar)
     read(ifpar, nmslpm, iostat=istat)
     call cstnml(jfpar, 'dnsgrd', 'nmslpm', istat)
     write(jfpar, nmslpm)
     call rewnml(ifpar, jfpar)
     read(ifpar, nmmlep, iostat=istat)
     call cstnml(jfpar, 'dnsgrd', 'nmmlep', istat)
     write(jfpar, nmmlep)

     call secoef( &
        &   c0(kstr), c1(kstr), c2(kstr), c3(kstr), &
        &   c4(kstr), c5(kstr), c6(kstr), &
        &   d0(kstr), d1(kstr), d2(kstr), d3(kstr),  d4(kstr), &
        &   d5(kstr), d6(kstr), d7(kstr), d8(kstr),  d9(kstr))

!    === Mixed layer eddy parameterization
!        from Fox-Kemper and Ferrari(2008) ===
!    *** Coarse-resolution modification is applied.
!        The last row in CXPSY/CYPSX corresponds to ds,
!        and is cancelled out by the row just before it. *** 
     pi = atan( 1.d0 )*4.d0
     omega = 2.d0 * pi / 86400.d0
     cormin = 2.d0 * omega * sin( pi*abs(fminlt)/180.d0 )
     kzmin = mzmin + kstr - 1

     do ij = nxdim+2, nxydim
        if (ocoamp) then
           cxpsy(ij) = ce &
             &       / sqrt( &
             &         max( abs(0.5d0*(cor(ij+lsw) + cor(ij+lw))), &
             &              cormin )**2.0d0 &
             &         + 1.0d0 / (8.64d4 * taumle)**2.0d0 ) &
             &       * dx * 0.5d0 * (hxt(ij) + hxt(ij+lw))
           cypsx(ij) = ce &
             &       / sqrt( &
             &         max( abs(0.5d0*(cor(ij+lsw) + cor(ij+ls))), &
             &              cormin )**2.0d0 &
             &         + 1.0d0 / (8.64d4 * taumle)**2.0d0 ) &
             &       * dym(ij+ls) * 0.5d0 * (hyt(ij) + hyt(ij+ls))        
           czpsy(ij) = ce &
             &       / sqrt( &
             &         max( abs(0.25d0* &
             &              ( cor(ij    ) + cor(ij+lw ) &
             &              + cor(ij+ls ) + cor(ij+lsw))), cormin ) &
             &                                              **2.0d0 &
             &         + 1.0d0 / (8.64d4 * taumle)**2.0d0 ) &
             &         * dx * hxt(ij)
           czpsx(ij) = ce &
             &       / sqrt( &
             &         max( abs(0.25d0* &
             &              ( cor(ij    ) + cor(ij+lw ) &
             &              + cor(ij+ls ) + cor(ij+lsw))), cormin ) &
             &                                              **2.0d0 &
             &         + 1.0d0 / (8.64d4 * taumle)**2.0d0 ) &
             &         * dy(ij) * hyt(ij)
        else
           cxpsy(ij) = ce &
             &       / sqrt( &
             &         max( abs(0.5d0*(cor(ij+lsw) + cor(ij+lw))), &
             &              cormin )**2.0d0 )
           cypsx(ij) = ce &
             &       / sqrt( &
             &         max( abs(0.5d0*(cor(ij+lsw) + cor(ij+ls))), &
             &              cormin )**2.0d0 )
           czpsy(ij) = ce &
             &       / sqrt( &
             &         max( abs(0.25d0* &
             &              ( cor(ij    ) + cor(ij+lw ) &
             &              + cor(ij+ls ) + cor(ij+lsw))), cormin ) &
             &                                              **2.0d0 )
           czpsx(ij) = ce &
             &       / sqrt( &
             &         max( abs(0.25d0* &
             &              ( cor(ij    ) + cor(ij+lw ) &
             &              + cor(ij+ls ) + cor(ij+lsw))), cormin ) &
             &                                              **2.0d0 )
        end if
     end do
  end if

  do n = 1, ntdim
     do k = 1, nzdim
        do ij = 1, nxydim
           dtdx(ij, k, n) = 0.d0
           dtdy(ij, k, n) = 0.d0
           dtfdz(ij, k, n) = 0.d0
        end do
     end do
  end do

  do k = kstr, kend
     do ij = 1, nxydim
        tl = ty(ij, k, 1)
        sl = ty(ij, k, 2)
        p1 = c0(k) &
          &  + (c1(k) + (c2(k) + c3(k) * tl) * tl) * tl &
          &  + (c4(k) + c5(k) * tl + c6(k) * sl) * sl
        p2 = d0(k) &
          &  + (d1(k) + (d2(k) + (d3(k) + d4(k) * tl) * tl) * tl) * tl &
          &  + (d5(k) + (d6(k) + d7(k) * tl * tl) * tl &
          &           + (d8(k) + d9(k) * tl * tl) * sqrt(sl)) * sl
        r(ij, k) = p1 / p2
        p1 = c0(kstr) &
          &  + (c1(kstr) + (c2(kstr) + c3(kstr) * tl) * tl) * tl &
          &  + (c4(kstr) + c5(kstr) * tl + c6(kstr) * sl) * sl
        p2 = d0(kstr) &
          &  + (d1(kstr) + (d2(kstr) + &
          &                (d3(kstr) + d4(kstr) * tl) * tl) * tl) * tl &
          &  + (d5(kstr) + (d6(kstr) + d7(kstr) * tl * tl) * tl &
          &         + (d8(kstr) + d9(kstr) * tl * tl) * sqrt(sl)) * sl
        rsigth(ij, k) = p1 / p2
     end do
  end do

  do k = kstr, kend
     do ij = ijtstr, ijtend+nxdim
        tl = ty(ij, k, 1)
        sl = ty(ij, k, 2)
        p1 = c0(k) &
          &  + (c1(k) + (c2(k) + c3(k) * tl) * tl) * tl &
          &  + (c4(k) + c5(k) * tl + c6(k) * sl) * sl 
        p2 = d0(k) &
          &  + (d1(k) + (d2(k) + (d3(k) + d4(k) * tl) * tl) * tl) * tl &
          &  + (d5(k) + (d6(k) + d7(k) * tl * tl) * tl &
          &           + (d8(k) + d9(k) * tl * tl) * sqrt(sl)) * sl 
        rl = p1 / p2

        tl = ty(ij+lw, k, 1)
        sl = ty(ij+lw, k, 2)
        p1 = c0(k) &
          &  + (c1(k) + (c2(k) + c3(k) * tl) * tl) * tl &
          &  + (c4(k) + c5(k) * tl + c6(k) * sl) * sl
        p2 = d0(k) &
          &  + (d1(k) + (d2(k) + (d3(k) + d4(k) * tl) * tl) * tl) * tl &
          &  + (d5(k) + (d6(k) + d7(k) * tl * tl) * tl &
          &           + (d8(k) + d9(k) * tl * tl) * sqrt(sl)) * sl
        rlw = p1 / p2

        tl = ty(ij+ls, k, 1)
        sl = ty(ij+ls, k, 2)
        p1 = c0(k) &
          &  + (c1(k) + (c2(k) + c3(k) * tl) * tl) * tl &
          &  + (c4(k) + c5(k) * tl + c6(k) * sl) * sl
        p2 = d0(k) &
          &  + (d1(k) + (d2(k) + (d3(k) + d4(k) * tl) * tl) * tl) * tl &
          &  + (d5(k) + (d6(k) + d7(k) * tl * tl) * tl &
          &           + (d8(k) + d9(k) * tl * tl) * sqrt(sl)) * sl
        rls = p1 / p2

        dtdx(ij, k, 1) = (rl - rlw) * rx * 2.d0 &
          &              / (hxt(ij) + hxt(ij+lw)) * &
          &              amskt(ij, k) * amskt(ij+lw, k)
        dtdy(ij, k, 1) = (rl - rls) * rym(ij+ls) * 2.d0 &
          &              / (hyt(ij) + hyt(ij+ls)) * &
          &              amskt(ij, k) * amskt(ij+ls, k)
     end do
  end do

  do k = kstr+1, kend
     do ij = ijtstr-nxdim, ijtend+nxdim
        tl = ty(ij, k, 1)
        sl = ty(ij, k, 2)
        p1 = c0(k) &
           &   + (c1(k) + (c2(k) + c3(k) * tl) * tl) * tl &
           &   + (c4(k) + c5(k) * tl + c6(k) * sl) * sl
        p2 = d0(k) &
           &   + (d1(k) + (d2(k) + (d3(k) + d4(k) * tl) * tl) * tl) * tl &
           &   + (d5(k) + (d6(k) + d7(k) * tl * tl) * tl &
           &            + (d8(k) + d9(k) * tl * tl) * sqrt(sl)) * sl
        rl = p1 / p2

        tl = ty(ij, k-1, 1)
        sl = ty(ij, k-1, 2)
        p1 = c0(k) &
           &   + (c1(k) + (c2(k) + c3(k) * tl) * tl) * tl &
           &   + (c4(k) + c5(k) * tl + c6(k) * sl) * sl
        p2 = d0(k) &
           &   + (d1(k) + (d2(k) + (d3(k) + d4(k) * tl) * tl) * tl) * tl &
           &   + (d5(k) + (d6(k) + d7(k) * tl * tl) * tl &
           &            + (d8(k) + d9(k) * tl * tl) * sqrt(sl)) * sl
        rlu = p1 / p2

        dtfdz(ij, k, 1) = min((rlu - rl) / dzm(ij, k), 0.d0) * &
           &              amftz(ij, k)
     end do
  end do

! === Mixed layer eddy parameterization
!     from Fox-Kemper and Ferrari(2008) ===

  do k = kstr, kstr+kz-1
     do ij = 1, nxydim
        dzsig (ij, k) = zbot * ds(k)
        dzmsig(ij, k) = zbot * dsm(k)
     end do
  end do
  do k = kstr+kz, kend
     do ij = 1, nxydim
        dzsig (ij, k) = dz(ij, k)
        dzmsig(ij, k) = dzm(ij, k)
     end do
  end do
  do k = kend+1, nzdim
     do ij = 1, nxydim
        dzsig (ij, k) = dzsig(ij, kend)
        dzmsig(ij, k) = dzmsig(ij, kend)
     end do
  end do
  
  do ij = 1, nxydim
     zt(ij, kstr) = 0.d0
     ztm(ij, kstr) = 0.5d0 * dzsig(ij, kstr)
  end do
  do k = kstr, kend
     do ij = 1, nxydim
        zt(ij, k+1) = zt(ij, k) + dzsig(ij, k)
        ztm(ij, k+1) = ztm(ij, k) + dzmsig(ij, k+1)
     end do
  end do

  do k = kstr, kend
     do ij = 1, nxydim
        zhmld(ij, k) = 0.0d0
     end do
  end do
  
  do ij = 1, nxydim
     rmavdx(ij) = 0.0d0
     rmavdy(ij) = 0.0d0
     hmldx(ij) = 0.0d0
     hmldy(ij) = 0.0d0
  end do

! calculating mixed layer depth hmld and mld-averaged buoyancy rmavez

!  determine HMLD by N^2
!  do ij = 1, nxydim
!     in2dz = (-1.0d-3) * gravit / rhoo &
!       &   * dtfdz(ij, kstr+1, 1) * ztm(ij, kstr+1)
!     n2min = (-1.0d-3) * gravit / rhoo &
!       &   * dtfdz(ij, kstr+1, 1)
!     hmld(ij) = ztm(ij, kstr+1)
!     rmavez(ij) = r(ij, kstr) * dzsig(ij, kstr) &
!       &        + 0.5d0 * r(ij, kstr+1) * dzsig(ij, kstr+1)
!     do k = kstr+2, min(nbot(ij), mz)
!        n2l = (-1.0d-3) * gravit / rhoo * dtfdz(ij, k, 1)
!        in2dz = in2dz + n2l * dzmsig(ij, k)
!        n2min = min(n2min, n2l)
!        if (      ( (n2l - n2min) * ztm(ij, k) ) &
!          &  .ge. (cm * in2dz) ) then
!           exit
!        end if
!        hmld(ij) = ztm(ij, k)
!        rmavez(ij) = rmavez(ij) + 0.5d0 *  &
!          &        (  r(ij, k-1) * dzsig(ij, k-1) &
!          &         + r(ij, k) * dzsig(ij, k) )
!     end do
!     rmavez(ij) = rmavez(ij) / hmld(ij)
!  end do

! determine HMLD by sigma_theta
  do ij = 1, nxydim
     hmld(ij) = ztm(ij, kstr)
     zhmld(ij, kstr) = ztm(ij, kstr)
     kmld(ij) = kstr
     rmavez(ij) = rsigth(ij, kstr) * ztm(ij, kstr)
     rsigbt = rsigth(ij, kstr) + drsig
     rsigdf = 0.d0
     do k = kstr+1, min(nbot(ij), mz+kstr-1)
        if ( rsigth(ij, k) .ge. rsigbt ) then
           dhmld = dzmsig(ij, k) * (rsigbt - rsigth(ij, k-1)) &
             &   / (rsigth(ij, k) - rsigth(ij, k-1))
           hmld(ij) = hmld(ij) + dhmld
           if (dhmld .le. 0.5d0*dzsig(ij, k-1)) then
              zhmld(ij, k-1) = zhmld(ij, k-1) + dhmld
           else
              zhmld(ij, k-1) = dzsig(ij, k-1)
              zhmld(ij, k) = dhmld - 0.5d0 * dzsig(ij, k-1)
           end if
           rmavez(ij) = rmavez(ij) + &
             &          0.5d0 * (rsigth(ij, k-1)+rsigbt) * dhmld
           rsigdf = drsig
           exit
        end if
        hmld(ij) = hmld(ij) + dzmsig(ij, k)
        zhmld(ij, k-1) = dzsig(ij, k-1)
        zhmld(ij, k) = 0.5d0 * dzsig(ij, k)
        kmld(ij) = k
        rmavez(ij) = rmavez(ij) + 0.5d0 * &
          &        ( rsigth(ij, k-1) * dzsig(ij, k-1) &
          &        + rsigth(ij, k  ) * dzsig(ij, k  ) )
        rsigdf = rsigth(ij, k) - rsigth(ij, kstr)
     end do
     rmavez(ij) = rmavez(ij) / hmld(ij)
     nbv(ij) = sqrt(max(1.0d-3*gravit/rhoo*rsigdf/hmld(ij),0.0d0))
  end do

  do ij = ijtstr-nxdim, ijtend+nxdim
     if (ocoamp) then
        lf(ij) = max(4.0d0*nbv(ij)*hmld(ij)/ &
          &          max( abs( cor(ij   ) + cor(ij+lw ) &
          &                  + cor(ij+ls) + cor(ij+lsw) ), eps), lfmin)
     else
        lf(ij) = 1.0d0
     end if
  end do

!  if (ofltdm) then
!#ifdef OPT_TRIPOLE
!     call shift1( &
!       &            hmld, &
!       &           nxdim,  nydim,      1, &
!       &           1.0d0,      0,      0 )
!#else
!     call shift1( &
!       &            hmld, &
!       &           nxdim,  nydim,      1)
!#endif
!     do n = 1, nfltdm
!        do ij = ijstr, ijend
!           if (amskt(ij, kstr) .eq. 1.0d0) then
!              hmld1(ij) = ( 4.0d0 * hmld(ij) * amskt(ij, kstr) &
!                &         + hmld(ij+ls) * amskt(ij+ls, kstr) &
!                &         + hmld(ij+ln) * amskt(ij+ln, kstr) &
!                &         + hmld(ij+lw) * amskt(ij+lw, kstr) &
!                &         + hmld(ij+le) * amskt(ij+le, kstr) ) &
!                &       / ( 4.0d0 * amskt(ij, kstr) &
!                &         + amskt(ij+ls, kstr) + amskt(ij+ln, kstr) &
!                &         + amskt(ij+lw, kstr) + amskt(ij+le, kstr) )
!           else
!              hmld1(ij) = hmld(ij)
!           end if
!        end do
!        do ij = ijstr, ijend
!           hmld(ij) = hmld1(ij)
!        end do
!#ifdef OPT_TRIPOLE
!     call shift1( &
!       &            hmld, &
!       &           nxdim,  nydim,      1, &
!       &           1.0d0,      0,      0 )
!#else
!     call shift1( &
!       &            hmld, &
!       &           nxdim,  nydim,      1)
!#endif
!     end do
!  end if

!  if (ofltrm) then
!#ifdef OPT_TRIPOLE
!     call shift1( &
!       &          rmavez, &
!       &           nxdim,  nydim,      1, &
!       &           1.0d0,      0,      0 )
!#else
!     call shift1( &
!       &          rmavez, &
!       &           nxdim,  nydim,      1)
!#endif
!     do n = 1, nfltrm
!        do ij = ijstr, ijend
!           if (amskt(ij, kstr) .eq. 1.0d0) then
!              rmav1(ij) = ( 4.0d0 * rmavez(ij) * amskt(ij, kstr) &
!                &         + rmavez(ij+ls) * amskt(ij+ls, kstr) &
!                &         + rmavez(ij+ln) * amskt(ij+ln, kstr) &
!                &         + rmavez(ij+lw) * amskt(ij+lw, kstr) &
!                &         + rmavez(ij+le) * amskt(ij+le, kstr) ) &
!                &       / ( 4.0d0 * amskt(ij, kstr) &
!                &         + amskt(ij+ls, kstr) + amskt(ij+ln, kstr) &
!                &         + amskt(ij+lw, kstr) + amskt(ij+le, kstr) )
!           else
!              rmav1(ij) = rmavez(ij)
!           end if
!        end do
!        do ij = ijstr, ijend
!           rmavez(ij) = rmav1(ij)
!        end do
!#ifdef OPT_TRIPOLE
!        call shift1( &
!          &          rmavez, &
!          &           nxdim,  nydim,      1, &
!          &           1.0d0,      0,      0 )
!#else
!        call shift1( &
!          &          rmavez, &
!          &           nxdim,  nydim,      1)
!#endif
!     end do
!  end if

  do ij = 1, nxydim
     if (kmld(ij) .lt. kzmin) then
        hmld(ij) = 0.0d0
     end if
     zmld0(ij) = kmld(ij) - kstr + 1
  end do

! for diagnosis
!  call chekin(   hmld,   'HMLD', nx, ny,  1, nxydim, 'SFC')
!  call chekin(  zmld0,   'KMLD', nx, ny,  1, nxydim, 'SFC')
!  call chekin( rmavez, 'RMAVEZ', nx, ny,  1, nxydim, 'SFC')
!  call chekin(      r,      'R', nx, ny, nz, nxyzdm, 'OCN')
!  call chekin(     lf,     'LF', nx, ny,  1, nxydim, 'SFC')

  do k = kstr, mz+kstr-1
     do ij = ijtstr, ijtend+nxdim
        dhmld = min(zhmld(ij, k), zhmld(ij+lw, k))
        rmavdx(ij) = rmavdx(ij) + &
          &          dhmld * (-1.0d-3) * gravit / rhoo &
          &        * rx * 2.d0 / (hxt(ij) + hxt(ij+lw)) &
          &        * (rsigth(ij, k) - rsigth(ij+lw, k))
        hmldx(ij) = hmldx(ij) + dhmld
        dhmld = min(zhmld(ij, k), zhmld(ij+ls, k))
        rmavdy(ij) = rmavdy(ij) + &
          &          dhmld * (-1.0d-3) * gravit / rhoo &
          &        * rym(ij+ls) * 2.d0 / (hyt(ij) + hyt(ij+ls)) &
          &        * (rsigth(ij, k) - rsigth(ij+ls, k))
        hmldy(ij) = hmldy(ij) + dhmld
     end do
  end do

  do ij = ijtstr, ijtend+nxdim
     rmavdx(ij) = rmavdx(ij) / max(hmldx(ij), eps) &
       &        * amskt(ij, kstr) * amskt(ij+lw, kstr)
     rmavdy(ij) = rmavdy(ij) / max(hmldy(ij), eps) &
       &        * amskt(ij, kstr) * amskt(ij+ls, kstr)
  end do

! calculating xpsiy and ypsix
  do k = kstr, kend
     do ij = ijtstr, ijtend+nxdim
        hmldt = min(hmld(ij), hmld(ij+lw))
        cpsi = cxpsy(ij) &
          &  * 2.0d0 / (lf(ij) + lf(ij+lw)) &
          &  * hmldt * hmldt &
          &  * (-1.0d0) * rmavdx(ij)
        muzh = ( 1.0d0 - 2.0d0 * &
          &      min(min(ztm(ij, k), ztm(ij+lw, k)) &
          &      / max(hmldt, eps), 1.0d0) &
          &    ) ** 2.0d0
        muzx(ij, k) = ( 1.0d0 - muzh ) &
          &         * ( 1.0d0 + 5.0d0 / 21.0d0 * muzh )
        xpsiy(ij, k) = cpsi &
          &          * ( 1.0d0 - muzh ) * ( 1.0d0 + 5.0d0 / 21.0d0 * muzh ) &
          &          * amskt(ij, k) * amskt(ij+lw, k)
        xpsiy(ij, k) = sign(1.0d0, xpsiy(ij, k)) * &
          &            min(abs(xpsiy(ij, k)), vscl*dzsig(ij, k))
        muzh = ( 1.0d0 - 2.0d0 * &
          &      min(min(zt(ij, k+1), zt(ij+lw, k+1)) &
          &      / max(hmldt, eps), 1.0d0) &
          &    ) ** 2.0d0
        xpsiyz(ij, k+1) = cpsi &
          &          * ( 1.0d0 - muzh ) * ( 1.0d0 + 5.0d0 / 21.0d0 * muzh ) &
          &          * amskt(ij, k+1) * amskt(ij+lw, k+1)
        xpsiyz(ij, k+1) = sign(1.0d0, xpsiyz(ij, k+1)) * &
          &            min(abs(xpsiyz(ij, k+1)), vscl*dzmsig(ij, k+1))
     end do
  end do
  do k = kstr, kend
     do ij = ijtstr, ijtend+nxdim
        hmldt = min(hmld(ij), hmld(ij+ls))
        cpsi = cypsx(ij) &
          &  * 2.0d0 / (lf(ij) + lf(ij+ls)) &
          &  * hmldt * hmldt &
          &  * rmavdy(ij)
        muzh = ( 1.0d0 - 2.0d0 * &
          &      min(min(ztm(ij, k), ztm(ij+ls, k)) &
          &      / max(hmldt, eps), 1.0d0) &
          &    ) ** 2.0d0
        muzy(ij, k) = ( 1.0d0 - muzh ) &
          &         * ( 1.0d0 + 5.0d0 / 21.0d0 * muzh )
        ypsix(ij, k) = cpsi &
          &          * ( 1.0d0 - muzh ) * ( 1.0d0 + 5.0d0 / 21.0d0 * muzh ) &
          &          * amskt(ij, k) * amskt(ij+ls, k)
        ypsix(ij, k) = sign(1.0d0, ypsix(ij, k)) * &
          &            min(abs(ypsix(ij, k)), vscl*dzsig(ij, k))
        muzh = ( 1.0d0 - 2.0d0 * &
          &      min(min(zt(ij, k+1), zt(ij+ls, k+1)) &
          &      / max(hmldt, eps), 1.0d0) &
          &    ) ** 2.0d0
        ypsixz(ij, k+1) = cpsi &
          &          * ( 1.0d0 - muzh ) * ( 1.0d0 + 5.0d0 / 21.0d0 * muzh ) &
          &          * amskt(ij, k+1) * amskt(ij+ls, k+1)
        ypsixz(ij, k+1) = sign(1.0d0, ypsixz(ij, k+1)) * &
          &            min(abs(ypsixz(ij, k+1)), vscl*dzmsig(ij, k+1))
     end do
  end do

!  if (ofltps) then
!#ifdef OPT_TRIPOLE
!     call shift1( &
!       &           xpsiy, &
!       &           nxdim,  nydim,  nzdim, &
!       &          -1.0d0,      1,      0 )
!     call shift1( &
!       &           ypsix, &
!       &           nxdim,  nydim,  nzdim, &
!       &          -1.0d0,      0,      1 )
!#else
!     call shift2( &
!       &           xpsiy,  ypsix, &
!       &           nxdim,  nydim,  nzdim)
!#endif
!     do n = 1, nfltps
!        do k = kstr, kend
!           do ij = ijstr, ijend
!              if (amftx(ij, k) .eq. 1.0d0) then
!                 xpsiy1(ij, k) = &
!                   &           ( 4.0d0 * xpsiy(ij, k) * amftx(ij, k) &
!                   &           + xpsiy(ij+ls, k) * amftx(ij+ls, k) &
!                   &           + xpsiy(ij+ln, k) * amftx(ij+ln, k) &
!                   &           + xpsiy(ij+lw, k) * amftx(ij+lw, k) &
!                   &           + xpsiy(ij+le, k) * amftx(ij+le, k) ) &
!                   &         / ( 4.0d0 * amftx(ij, k) &
!                   &           + amftx(ij+ls, k) + amftx(ij+ln, k) &
!                   &           + amftx(ij+lw, k) + amftx(ij+le, k) )
!              else
!                 xpsiy1(ij, k) = xpsiy(ij, k)
!              end if
!              if (amfty(ij, k) .eq. 1.0d0) then
!                 ypsix1(ij, k) = &
!                   &           ( 4.0d0 * ypsix(ij, k) * amfty(ij, k) &
!                   &           + ypsix(ij+ls, k) * amfty(ij+ls, k) &
!                   &           + ypsix(ij+ln, k) * amfty(ij+ln, k) &
!                   &           + ypsix(ij+lw, k) * amfty(ij+lw, k) &
!                   &           + ypsix(ij+le, k) * amfty(ij+le, k) ) &
!                   &         / ( 4.0d0 * amfty(ij, k) &
!                   &           + amfty(ij+ls, k) + amfty(ij+ln, k) &
!                   &           + amfty(ij+lw, k) + amfty(ij+le, k) )
!              else
!                 ypsix1(ij, k) = ypsix(ij, k)
!              end if
!           end do
!        end do
!        do k = kstr, kend
!           do ij = ijstr, ijend
!              xpsiy(ij, k) = xpsiy1(ij, k)
!              ypsix(ij, k) = ypsix1(ij, k)
!           end do
!        end do
!#ifdef OPT_TRIPOLE
!     call shift1( &
!       &           xpsiy, &
!       &           nxdim,  nydim,  nzdim, &
!       &          -1.0d0,      1,      0 )
!     call shift1( &
!       &           ypsix, &
!       &           nxdim,  nydim,  nzdim, &
!       &          -1.0d0,      0,      1 )
!#else
!     call shift2( &
!       &           xpsiy,  ypsix, &
!       &           nxdim,  nydim,  nzdim)
!#endif
!     end do
!  end if
        
  do k = kstr+1, kend
     do ij = ijtstr, ijtend
        zpsix(ij, k) = 0.5d0 * &
          &   ( ypsixz(ij, k) + ypsixz(ij+ln, k) ) * amftz(ij, k)
        zpsiy(ij, k) = 0.5d0 * &
          &   ( xpsiyz(ij, k) + xpsiyz(ij+le, k) ) * amftz(ij, k)
     end do
  end do

! for diagnosis
!  call chekin(  cxpsy,  'CXPSY', nx, ny,  1, nxydim, 'SFC')
!  call chekin(  cypsx,  'CYPSX', nx, ny,  1, nxydim, 'SFC')
!  call chekin(  xpsiy,  'XPSIY', nx, ny, nz, nxyzdm, 'OCN')
!  call chekin(  ypsix,  'YPSIX', nx, ny, nz, nxyzdm, 'OCN')
!  call chekin(  zpsix,  'ZPSIX', nx, ny, nz, nxyzdm, 'OCN')
!  call chekin(  zpsiy,  'ZPSIY', nx, ny, nz, nxyzdm, 'OCN')
!  call chekin( rmavdx, 'RMAVDX', nx, ny,  1, nxydim, 'SFC')
!  call chekin( rmavdy, 'RMAVDY', nx, ny,  1, nxydim, 'SFC')
!  call chekin(   muzx,   'MUZX', nx, ny, nz, nxyzdm, 'OCN')
!  call chekin(   muzy,   'MUZY', nx, ny, nz, nxyzdm, 'OCN')

  do k = kstr, kend
     do ij = ijtstr, ijtend+nxdim
        dzdx = dtdx(ij, k, 1) * 4.d0 &
           &   / (  dtfdz(ij, k  , 1) + dtfdz(ij+lw, k ,  1) &
           &      + dtfdz(ij, k+1, 1) + dtfdz(ij+lw, k+1, 1) - eps)
        xdzdx(ij, k) = min(slpmax, max(-slpmax, dzdx))
        dzdy = dtdy(ij, k, 1) * 4.d0 * dym(ij+ls) &
           &   / (  (dtfdz(ij, k, 1) + dtfdz(ij, k+1, 1)) * dy(ij) &
           &   + (dtfdz(ij+ls, k, 1) + dtfdz(ij+ls, k+1, 1)) * dy(ij+ls) &
           &      - eps)
        ydzdy(ij, k) = min(slpmax, max(-slpmax, dzdy))
     end do
  end do

  do k = kstr+1, kend
     do ij = ijtstr, ijtend
        dzdx = ((dtdx(ij, k-1, 1)+dtdx(ij+le, k-1, 1))*dz(ij, k-1) &
           &  + (dtdx(ij, k  , 1)+dtdx(ij+le, k  , 1))*dz(ij, k)) * &
           &   0.25d0 / dzm(ij, k) / (dtfdz(ij, k, 1) - eps)
        zdzdx(ij, k) = min(slpmax, max(-slpmax, dzdx)) * &
           &           amftz(ij, k)
        dzdy = &
           &   (  (dtdy(ij, k-1, 1)+dtdy(ij+ln, k-1, 1))*dz(ij, k-1) &
           &    + (dtdy(ij, k  , 1)+dtdy(ij+ln, k  , 1))*dz(ij, k)) * &
           &     0.25d0 / dzm(ij, k) / (dtfdz(ij, k, 1) - eps)
        zdzdy(ij, k) = min(slpmax, max(-slpmax, dzdy)) * &
           &           amftz(ij, k)
     end do
  end do

  do n = 1, ntdim
     do k = 1, nzdim
        do ij = 1, nxydim
           dtdx(ij, k, n) = 0.d0
           dtdy(ij, k, n) = 0.d0
           dtfdz(ij, k, n) = 0.d0
        end do
     end do
  end do

  do n = 1, ntdim
     do k = kstr, kend
        do ij = ijtstr, ijtend+nxdim
           dtdx(ij, k, n) = (tx(ij, k, n) - tx(ij+lw, k, n)) * rx * &
              &           2.d0 / (hxt(ij) + hxt(ij+lw)) * &
              &           amskt(ij, k) * amskt(ij+lw, k) 
           dtdy(ij, k, n) = (tx(ij, k, n) - tx(ij+ls, k, n)) * &
              &           rym(ij+ls) * &
              &           2.d0 / (hyt(ij) + hyt(ij+ls)) * &
              &           amskt(ij, k) * amskt(ij+ls, k)
        end do
     end do
  end do

  do n = 1, ntdim
     do k = kstr+1, kend
        do ij = ijtstr-nxdim, ijtend+nxdim
           dtfdz(ij, k, n) = (tx(ij, k-1, n) - tx(ij, k, n)) &
              &              / dzm(ij, k) * amftz(ij, k)
        end do
     end do
  end do

  do n = 1, ntdim
     do k = kstr, kend
        do ij = ijtstr, ijtend+nxdim
           xdtdz(ij, k, n) = &
              &   (  dtfdz(ij, k, n)   + dtfdz(ij+lw, k, n) &
              &    + dtfdz(ij, k+1, n) + dtfdz(ij+lw, k+1, n)) *  0.25d0
           ydtdz(ij, k, n) = &
              &   (  (dtfdz(ij, k, n)    + dtfdz(ij, k+1, n)) * &
              &      dy(ij) &
              &    + (dtfdz(ij+ls, k, n) + dtfdz(ij+ls, k+1, n)) * &
              &      dy(ij+ls)) * &
              &   rym(ij+ls) * 0.25d0
        end do
     end do
  end do

  do n = 1, ntdim
     do k = kstr+1, kend
        do ij = ijtstr, ijtend
           zdtdx(ij, k, n) = &
              &   (  (dtdx(ij, k-1, n) + dtdx(ij+le, k-1, n)) * &
              &      dz(ij, k-1) &
              &    + (dtdx(ij, k, n) + dtdx(ij+le, k, n)) * &
              &      dz(ij, k)) &
              &   / dzm(ij, k) * 0.25d0 * amftz(ij, k)
           zdtdy(ij, k, n) = &
              &   (  (dtdy(ij, k-1, n) + dtdy(ij+ln, k-1, n)) * &
              &      dz(ij, k-1) &
              &    + (dtdy(ij, k, n) + dtdy(ij+ln, k, n)) * &
              &      dz(ij, k)) &
              &   / dzm(ij, k) * 0.25d0 * amftz(ij, k)
        end do
     end do
  end do

  return

end subroutine dnsgrd
#ifdef OPT_BBL
! *********************************************************************

subroutine flxtrb( &
  &    adt,  diffz, &
  &     tx,     ty,     uy,     vy, &
  &      w,    ahv )
      
  real(8), intent(out) ::    adt(nxydim, nzdim, ntdim)
  real(8), intent(out) ::  diffz(nxydim, nzdim)
  real(8), intent(in)  ::     tx(nxydim, nzdim, ntdim)
  real(8), intent(in)  ::     ty(nxydim, nzdim, ntdim)
  real(8), intent(in)  ::     uy(nxydim, nzdim),     vy(nxydim, nzdim)
  real(8), intent(in)  ::      w(nxydim, nzdim),    ahv(nxydim, nzdim)

  integer ::  idummy

! ---- only a dummy routine to preserve a consistency
  idummy = 1

  return

end subroutine flxtrb
#endif
! *********************************************************************

subroutine chkftx
  use qckot

  if (oinit .or. ofinal) then
     return
  end if

  call chekin(   ftx,  'FTX', &
     &            nx,     ny,     nz, nxyzdm, 'OCN')
  call chekin(   fty,  'FTY', &
     &            nx,     ny,     nz, nxyzdm, 'OCN')
  call chekin(   ftz,  'FTZ', &
     &            nx,     ny,     nz, nxyzdm, 'OCN')
  call chekin(   ftx(1, 1, 2),  'FSX', &
     &            nx,     ny,     nz, nxyzdm, 'OCN')
  call chekin(   fty(1, 1, 2),  'FSY', &
     &            nx,     ny,     nz, nxyzdm, 'OCN')
  call chekin(   ftz(1, 1, 2),  'FSZ', &
     &            nx,     ny,     nz, nxyzdm, 'OCN')

  return
end subroutine chkftx

end module tflxt
