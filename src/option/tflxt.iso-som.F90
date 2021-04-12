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
!     '15.04.07  M.Kurogi: for MPI-IO
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
    &     ts,   zbot, &
    &    hxt,    hxu,    hyt,    hyu,    rxt,    ryt
  use zocmsk, only: &
#ifdef OPT_BBL
    & amsktb, &
#endif
    &  amskt,  amftx,  amfty,  amftz, &
    &   nbot
  use zocfil, only: &
    &    ncf
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
  real(8), save ::  s0 (nxydim, nzdim, ntdim)
  real(8), save ::  sm (nxydim, nzdim, ntdim)
  real(8), save ::  sx (nxydim, nzdim, ntdim), sxx(nxydim, nzdim, ntdim)
  real(8), save ::  sy (nxydim, nzdim, ntdim), syy(nxydim, nzdim, ntdim)
  real(8), save ::  sz (nxydim, nzdim, ntdim), szz(nxydim, nzdim, ntdim)
  real(8), save ::  sxy(nxydim, nzdim, ntdim), sxz(nxydim, nzdim, ntdim)
  real(8), save ::  syz(nxydim, nzdim, ntdim)   

  real(8), save ::  f0 (nxydim, nzdim)
  real(8), save ::  fm (nxydim, nzdim)
  real(8), save ::  fx (nxydim, nzdim), fxx(nxydim, nzdim)
  real(8), save ::  fy (nxydim, nzdim), fyy(nxydim, nzdim) 
  real(8), save ::  fz (nxydim, nzdim), fzz(nxydim, nzdim) 
  real(8), save ::  fxy(nxydim, nzdim), fxz(nxydim, nzdim)
  real(8), save ::  fyz(nxydim, nzdim) 

  real(8), save ::  vlmx(nxydim, nzdim), vlmy(nxydim, nzdim)
  real(8), save ::  vlmz(nxydim)
  real(8), save ::  r(nxyzdm)
  real(8), save ::  alf(nxydim, nzdim), uv(nxydim, nzdim)

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
     do n = 1, ntdim
        do k = 1, nzdim
           do ij = 1, nxydim
              f0 (ij, k)    = 0.d0
              fm (ij, k)    = 0.d0
              fx (ij, k)    = 0.d0
              fy (ij, k)    = 0.d0
              fz (ij, k)    = 0.d0
              fxx(ij, k)    = 0.d0
              fyy(ij, k)    = 0.d0
              fzz(ij, k)    = 0.d0
              fxy(ij, k)    = 0.d0
              fxz(ij, k)    = 0.d0
              fyz(ij, k)    = 0.d0
              s0 (ij, k, n) = 0.d0
              sm (ij, k, n) = 0.d0
              vlmx(ij, k)   = 0.d0
              vlmy(ij, k)   = 0.d0
              vlmz(ij)      = 0.d0
              alf(ij, k)    = 0.d0
           end do
        end do
     end do

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

#ifdef OPT_TRIPOLE
     call shift2( ahi3d,  ahg3d, &
          &       nxdim,  nydim,  nzdim, &
          &        1.d0,      0,      0 )
#else
     call shift2( &
          &       ahi3d,  ahg3d, &
          &       nxdim,  nydim,  nzdim )
#endif

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
     &     ty,     tx )

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

  call chekin(wzc, 'WZC', &
     &     'ocean vertical velocity on sigma coordinate', 'cm/s', &
     & nx, ny, nz, nxyzdm, 'OCLVMT')

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

           ftz(ij, k, n) = &
             & (  diffz(ij, k) * (tx(ij, ku, n) - tx(ij, k, n)) &
             &  - ( ahi3d(ij, k) + ahg3d(ij, k) ) * &
             &    (  zdzdx(ij, k) * zdtdx(ij, k, n) &
             &     + zdzdy(ij, k) * zdtdy(ij, k, n) ) &
             & ) * amftz(ij, k)

        end do
     end do
  end do

! ---- y diffusion flux of GM
  do n = 1, ntdim

     do k = kstr, kend
        do ij = ijtstr, ijtend+nxdim

           ijls = ij + ls

           fty(ij, k, n) = &
             &   (  ( ahh3d(ij, k) + ahi3d(ij, k) ) * rym(ijls) &
             &    / ( hyt(ij) + hyt(ijls) ) * &
             &      ( tx(ij, k, n) - tx(ijls, k, n) ) * 2.d0 &
             &    - ( ahi3d(ij, k) - ahg3d(ij, k) ) &
             &    * ydzdy(ij, k) * ydtdz(ij, k, n) ) &
             &    * ( hxu(ijls) + hxu(ij+lsw) ) * 0.5d0 * amfty(ij, k)

        end do
     end do

  end do

! ---- x diffusion flux of GM
  do n = 1, ntdim

     do k = kstr, kend
        do ij = ijtstr, ijtend+1

           ijlw = ij + lw
           
           ftx(ij, k, n) = &
             &   (  ( ahh3d(ij, k) + ahi3d(ij, k) ) * rx &
             &    / ( hxt(ij) + hxt(ijlw) ) &
             &    * ( tx(ij, k, n) - tx(ijlw, k, n) ) * 2.d0 &
             &    - ( ahi3d(ij, k) - ahg3d(ij, k) ) &
             &    * xdzdx(ij, k) * xdtdz(ij, k, n) ) &
             &    * ( hyu(ijlw) + hyu(ij+lsw) ) * 0.5d0 * amftx(ij, k)

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

! ---- SOM 

! ---- mass contained in a tracer grid
  do n = 1, ntdim

     do k = kstr, kstr+kz-1
        do ij = 1, nxydim

           sm(ij, k, n) = vlmz(ij) * hzbot(ij) * ds(k)

        end do
     end do

     do k = kstr+kz, kend
        do ij = 1, nxydim

           sm(ij, k, n) = vlmz(ij) * dz(ij, k)

        end do
     end do

     do k = kstr, kend
        do ij = 1, nxydim

           s0(ij, k, n) = sm(ij, k, n) * tx(ij, k, n)
               
        end do
     end do

     do ij = 1, nxydim
        sm(ij, kstr-1, n) = sm(ij, kstr, n)
        sm(ij, kend+1, n) = sm(ij, kend, n)
        s0(ij, kstr-1, n) = s0(ij, kstr, n)
        s0(ij, kend+1, n) = s0(ij, kend, n)
     end do

  end do

! ---- in X-direction

  do k = kstr, kend
     do ij = ijtstr-nxdim-1, ijtend+nxdim+1

        ijlw  = ij + lw
        ijlsw = ij + lsw
               
        uv(ij, k) = ( uy(ijlw , k) * vlmy(ijlw, k) &
          &         + uy(ijlsw, k) * vlmy(ijlsw, k) ) &
          &       * amskt(ij, k) * amskt(ijlw, k) * dy(ijlw)

     end do
  end do

  do n = 1, ntdim

! ---- undershoot limiter (Method B) of Morales Maqueda and Holloway (2006)
     do k = kstr, kend
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
          
           if ( abs( sx(ij, k, n) ) < 1.5d0 * s0m ) then
             
              sxx(ij, k, n) = min( s0m + sxp, &
                &                  max( abs( sx(ij, k, n) ) - s0m, &
                &                       sxx(ij, k, n) ) ) 

           else
                  
              sxx(ij, k, n) = min( s0m + sxp, &
                &                  max( s0m - sxp, sxx(ij, k, n) ) )
             
           end if

           sxy(ij, k, n) = min( s0m, max( - s0m, sxy(ij, k, n) ) )
           sxz(ij, k, n) = min( s0m, max( - s0m, sxz(ij, k, n) ) )
          
        end do
     end do

!    ---- overshoot limiter (Method B) of Morales Maqueda and Holloway (2006)
     do k = kstr, kend
        do ij = ijtstr-nxdim-1, ijtend+nxdim+1

           ijlw = ij + lw * nint( amskt(ij+lw, k) )
           ijle = ij + le * nint( amskt(ij+le, k) )

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
          
           if ( abs( sx(ij, k, n) ) < 1.5d0 * s0m ) then
             
!             sxx(ij, k, n) = min( s0m + sxp, &
!               &                  max( abs( sx(ij, k, n) ) - s0m, &
!               &                       sxx(ij, k, n) ) )
              sxx(ij, k, n) = min( s0m - abs( sx(ij, k, n) ), &
                &                  max( - s0m - sxp,          &
                &                       sxx(ij, k, n) ) )

           else
                  
!             sxx(ij, k, n) = min( s0m + sxp, &
!               &                  max( s0m - sxp, sxx(ij, k, n) ) )
              sxx(ij, k, n) = min( - s0m + sxp,       &
                &                  max( - s0m - sxp,  &
                &                        sxx(ij, k, n) ) )            
           end if

           sxy(ij, k, n) = min( s0m, max( - s0m, sxy(ij, k, n) ) )
           sxz(ij, k, n) = min( s0m, max( - s0m, sxz(ij, k, n) ) )
          
        end do
     end do

!---- bug fix 2
!     call shift1( sx (:,:,n), nxdim, nydim, nzdim )
!     call shift1( sxx(:,:,n), nxdim, nydim, nzdim )
!     call shift1( sxy(:,:,n), nxdim, nydim, nzdim )
!     call shift1( sxz(:,:,n), nxdim, nydim, nzdim )

!    ---- calculating ALF and MASS between box (i-1,j,k) <---> (i,j,k)
     do k = kstr, kend
        do ij = ijtstr-nxdim-1, ijtend+nxdim+1

           ijlw  = ij + lw
               
           if ( uv(ij, k) .gt. 0.d0 ) then
              fm(ijlw, k) = uv(ij, k) * ts 
              alf(ij, k)  = fm(ijlw, k) / sm(ijlw, k, n)
           else
              fm(ijlw, k) = - uv(ij, k) * ts
              alf(ij, k)  = fm(ijlw, k) / sm(ij  , k, n)
           end if

        end do
     end do

!    ---- calculating flux and moments between box (i-1,j,k) <---> (i,j,k)
     do k = kstr, kend
        do ij = ijtstr-nxdim-1, ijtend+nxdim+1

           ijlw  = ij + lw

!          ---- flux from (i-1) to (i),  when u > 0
           if ( uv(ij, k) .gt. 0.d0 ) then

              alfq  = alf(ij, k) * alf(ij, k)
              alf1  = 1.d0 - alf(ij, k)
              alf1q = alf1 * alf1

!             ---- moments be transported
              f0(ijlw, k) = &
                &      alf(ij, k) * ( s0(ijlw, k, n) &
                &        + alf1 * ( &
                &          sx(ijlw, k, n) &
                &          + ( alf1 - alf(ij, k) ) * sxx(ijlw, k, n) ) )

              fx(ijlw, k) = &
                &           alfq * ( sx(ijlw, k, n) &
                &         + 3.d0 * alf1 * sxx(ijlw, k, n) ) 

              fxx(ijlw, k) = alf(ij, k) * alfq * sxx(ijlw, k, n)

              fy (ijlw, k) = alf(ij, k) * ( sy(ijlw, k, n) &
                &                         + alf1 * sxy(ijlw, k, n) )
              fz (ijlw, k) = alf(ij, k) * ( sz(ijlw, k, n) & 
                &                         + alf1 * sxz(ijlw, k, n) )

              fxy(ijlw, k) = alfq       * sxy(ijlw, k, n)
              fxz(ijlw, k) = alfq       * sxz(ijlw, k, n)
              fyy(ijlw, k) = alf(ij, k) * syy(ijlw, k, n)
              fzz(ijlw, k) = alf(ij, k) * szz(ijlw, k, n)
              fyz(ijlw, k) = alf(ij, k) * syz(ijlw, k, n)

              ftx (ij, k, n) = ftx(ij, k, n) &
                &            - f0(ijlw, k) * tsiv * ry(ijlw)

           else
             
!             ---- flux from (i) to (i-1),  when u < 0

              alfq  = alf(ij, k) * alf(ij, k)
              alf1  = 1.d0 - alf(ij, k)
              alf1q = alf1 * alf1
             
              f0 (ijlw, k) = &
                &       alf(ij, k) * ( s0(ij, k, n) &
                &         - alf1 * (  sx(ij, k, n) &
                &           - ( alf1 - alf(ij, k) ) &
                &           * sxx(ij, k, n) ) )
              fx (ijlw, k) = alfq * ( &
                &            sx(ij, k, n) - 3.d0 * alf1 * sxx(ij, k, n) )
              fxx(ijlw, k) = alf(ij, k) * alfq * sxx(ij, k, n)
                  
              fy (ijlw, k) = alf(ij, k)  * ( sy(ij, k, n) &
                &                          - alf1 * sxy(ij, k, n) )
              fz (ijlw, k) = alf(ij, k)  * ( sz(ij, k, n) &
                &                          - alf1 * sxz(ij, k, n) )

              fxy(ijlw, k) = alfq       * sxy(ij, k, n)
              fxz(ijlw, k) = alfq       * sxz(ij, k, n)
              fyy(ijlw, k) = alf(ij, k) * syy(ij, k, n)
              fzz(ijlw, k) = alf(ij, k) * szz(ij, k, n)
              fyz(ijlw, k) = alf(ij, k) * syz(ij, k, n)
             
              ftx (ij, k, n) = ftx(ij, k, n) &
                &            + f0(ijlw, k) * tsiv * ry(ijlw)

           end if
          
        end do
     end do

!    ---- calculating flux and moments between box (i-1,j,k) <---> (i,j,k)
     do k = kstr, kend
        do ij = ijtstr-nxdim-1, ijtend+nxdim+1

           ijlw  = ij + lw

!          ---- flux from (i-1) to (i),  when u > 0
           if ( uv(ij, k) .gt. 0.d0 ) then

              alfq  = alf(ij, k) * alf(ij, k)
              alf1  = 1.d0 - alf(ij, k)
              alf1q = alf1 * alf1

!             ---- moments remaining 
              sm (ijlw, k, n) = sm(ijlw, k, n) - fm(ijlw, k)
              s0 (ijlw, k, n) = s0(ijlw, k, n) - f0(ijlw, k)

              sx (ijlw, k, n) = &
                &     alf1q * ( sx(ijlw, k, n) &
                &             - 3.d0 * alf(ij, k) * sxx(ijlw, k, n) )
              sxx(ijlw, k, n) = alf1 * alf1q * sxx(ijlw, k, n)

              sy (ijlw, k, n) = sy (ijlw, k, n) - fy (ijlw, k)
              syy(ijlw, k, n) = syy(ijlw, k, n) - fyy(ijlw, k)

              sz (ijlw, k, n) = sz (ijlw, k, n) - fz (ijlw, k)
              szz(ijlw, k, n) = szz(ijlw, k, n) - fzz(ijlw, k)

              sxy(ijlw, k, n) = alf1q * sxy(ijlw, k, n)
              sxz(ijlw, k, n) = alf1q * sxz(ijlw, k, n)
              syz(ijlw, k, n) = syz(ijlw, k, n) - fyz(ijlw, k)
             
           else
             
!             ---- flux from (i) to (i-1),  when u < 0

              alfq  = alf(ij, k) * alf(ij, k)
              alf1  = 1.d0 - alf(ij, k)
              alf1q = alf1 * alf1
             
              sm (ij, k, n) = sm(ij, k, n) - fm(ijlw, k)
              s0 (ij, k, n) = s0(ij, k, n) - f0(ijlw, k)

              sx (ij, k, n) = alf1q * ( &
                &    sx(ij, k, n) + 3.d0 * alf(ij, k) * sxx(ij, k, n)  )
              sxx(ij, k, n) = alf1 * alf1q * sxx(ij, k, n)

              sy (ij, k, n) = sy (ij, k, n) - fy (ijlw, k)
              syy(ij, k, n) = syy(ij, k, n) - fyy(ijlw, k)

              sz (ij, k, n) = sz (ij, k, n) - fz (ijlw, k)
              szz(ij, k, n) = szz(ij, k, n) - fzz(ijlw, k)

              sxy(ij, k, n) = alf1q * sxy(ij, k, n)
              sxz(ij, k, n) = alf1q * sxz(ij, k, n)
              syz(ij, k, n) = syz(ij, k, n) - fyz(ijlw, k)
             
           end if
          
        end do
     end do

!    ---- put the temporary moments (fi) into appropriate neighboring boxes
     do k = kstr, kend
        do ij = ijtstr-nxdim-1, ijtend+nxdim+1

           ijlw  = ij + lw

           if ( uv(ij, k) .gt. 0.d0 ) then

              sm(ij  , k, n) = sm(ij, k, n) + fm(ijlw, k)
              alf(ij, k)     = fm(ijlw, k) / sm(ij, k, n)

           else

              sm(ijlw, k, n) = sm(ijlw, k, n) + fm(ijlw, k)
              alf(ij, k)     = fm(ijlw, k) / sm(ijlw, k, n)

           end if

        end do
     end do

     do k = kstr, kend
        do ij = ijtstr-nxdim-1, ijtend+nxdim+1

           ijlw  = ij + lw

!          ---- flux from (i-1) to (i),  when u > 0
           if ( uv(ij, k) .gt. 0.d0 ) then

              alf1 = 1.d0 - alf(ij, k)

              tmp = alf(ij, k) * s0(ij, k, n) &
                & - alf1 * f0(ijlw, k)

              s0(ij, k, n) = s0(ij, k, n) + f0(ijlw, k)
              sxx(ij, k, n) = &
                &    alf(ij, k)  * alf(ij, k) * fxx(ijlw, k) &
                &  + alf1 * alf1 * sxx(ij, k, n) &
                &  + 5.d0 * (  alf(ij, k) * alf1 * ( & 
                &                sx(ij, k, n) - fx(ijlw, k) ) &
                &            - ( alf1 - alf(ij, k) ) * tmp )
              sx (ij, k, n) = &
                &    alf(ij, k) * fx(ijlw, k) &
                &  + alf1 * sx(ij, k, n) + 3.d0 * tmp 

              sxy(ij, k, n) = &
                &    alf(ij, k)  * fxy(ijlw, k) &
                &  + alf1 * sxy(ij, k, n) &
                &  + 3.d0 * ( alf(ij, k) * sy(ij, k, n) &
                &           - alf1 * fy(ijlw, k) )
              sxz(ij, k, n) = &
                &    alf(ij, k)  * fxz(ijlw, k) &
                &  + alf1 * sxz(ij, k, n) &
                &  + 3.d0 * ( alf(ij, k) * sz(ij, k, n) &
                &           - alf1 * fz(ijlw, k) )

              sy (ij, k, n) = sy (ij, k, n) + fy (ijlw, k)
              syy(ij, k, n) = syy(ij, k, n) + fyy(ijlw, k)

              sz (ij, k, n) = sz (ij, k, n) + fz (ijlw, k)
              szz(ij, k, n) = szz(ij, k, n) + fzz(ijlw, k)

              syz(ij, k, n) = syz(ij, k, n) + fyz(ijlw, k)

!          ---- flux from (i) to (i-1),  when u < 0
           else

              alf1 = 1.d0 - alf(ij, k)

              tmp = - alf(ij, k) * s0(ijlw, k, n) &
                &   + alf1 * f0(ijlw, k)

              s0(ijlw, k, n) = s0(ijlw, k, n) + f0(ijlw, k)
              sxx(ijlw, k, n) = &
                &      alf(ij, k) * alf(ij, k) * fxx(ijlw, k) &
                &    + alf1 * alf1 * sxx(ijlw, k, n) &
                &    + 5.d0 * (  alf(ij, k) * alf1 * ( - sx(ijlw, k, n) &
                &                               + fx(ijlw, k)    ) &
                &              + ( alf1 - alf(ij, k) ) * tmp  )
              sx (ijlw, k, n) = &
                &      alf(ij, k) * fx(ijlw, k) + alf1 * sx(ijlw, k, n) &
                &    + 3.d0 * tmp

              sxy(ijlw, k, n) = &
                &      alf(ij, k) * fxy(ijlw, k) &
                &    + alf1 * sxy(ijlw, k, n) &
                &    + 3.d0 * ( alf1 * fy(ijlw, k) &
                &             - alf(ij, k) * sy(ijlw, k, n) )
              sxz(ijlw, k, n) = &
                &      alf(ij, k) * fxz(ijlw, k) &
                &    + alf1 * sxz(ijlw, k, n) &
                &    + 3.d0 * ( alf1 * fz(ijlw, k) &
                &             - alf(ij, k) * sz(ijlw, k, n) )

              sy (ijlw, k, n) = sy (ijlw, k, n) + fy (ijlw, k)
              sz (ijlw, k, n) = sz (ijlw, k, n) + fz (ijlw, k)

              syy(ijlw, k, n) = syy(ijlw, k, n) + fyy(ijlw, k)
              szz(ijlw, k, n) = szz(ijlw, k, n) + fzz(ijlw, k)
              syz(ijlw, k, n) = syz(ijlw, k, n) + fyz(ijlw, k)

           end if
           
        end do
     end do
     
  end do

! ---- Y-direction

  do k = kstr, kend
     do ij = ijtstr-nxdim-1, ijtend+nxdim+1

        ijls  = ij + ls
        ijlsw = ij + lsw

        uv(ij, k) = ( vy(ijls , k) * vlmx(ijls , k) &
          &         + vy(ijlsw, k) * vlmx(ijlsw, k) ) &
          &       * amskt(ij, k) * amskt(ijls, k) 

     end do
  end do

  do n = 1, ntdim

!    ---- undershoot limiter (Method B) of Morales Maqueda and Holloway (2006)
     do k = kstr, kend
        do ij = ijtstr-nxdim-1, ijtend+nxdim+1

           ijls  = ij + ls * nint( amskt(ij+ls, k) )
           ijln  = ij + ln * nint( amskt(ij+ln, k) )
                
           s0m = s0(ij, k, n) &
             & - min( s0(ijls, k, n) / sm(ijls, k, n), &
             &        s0(ij  , k, n) / sm(ij  , k, n), &
             &        s0(ijln, k, n) / sm(ijln, k, n)  ) &
             & * sm(ij, k, n)

           s1m = sq3 * s0m
           sy(ij, k, n) = min( s1m, max( - s1m, sy(ij, k, n) ) )

           s0p = s0m * s0m
           sxp = sqrt( max( s0p - sy(ij, k, n) * sy(ij, k, n) * ci3, &
             &              eps ) )
               
           if ( abs( sy(ij, k, n) ) < 1.5d0 * s0m ) then

              syy(ij, k, n) = min( s0m + sxp, &
                &                  max( abs( sy(ij, k, n) ) - s0m, &
                &                       syy(ij, k, n) ) )

           else
                  
              syy(ij, k, n) = min( s0m + sxp, &
                &                  max( s0m - sxp, syy(ij, k, n) ) )

           end if
             
           sxy(ij, k, n) = min( s0m, max( - s0m, sxy(ij, k, n) ) )
           syz(ij, k, n) = min( s0m, max( - s0m, syz(ij, k, n) ) )

        end do
     end do

!    ---- overshoot limiter (Method B) of Morales Maqueda and Holloway (2006)
     do k = kstr, kend
        do ij = ijtstr-nxdim-1, ijtend+nxdim+1

           ijls  = ij + ls * nint( amskt(ij+ls, k) )
           ijln  = ij + ln * nint( amskt(ij+ln, k) )

           s0m = - s0(ij, k, n) &
             &   + max( s0(ijls, k, n) / sm(ijls, k, n), &
             &          s0(ij  , k, n) / sm(ij  , k, n), &
             &          s0(ijln, k, n) / sm(ijln, k, n)  ) &
             &   * sm(ij, k, n)

           s1m = sq3 * s0m
           sy(ij, k, n) = min( s1m, max( - s1m, sy(ij, k, n) ) )

           s0p = s0m * s0m
           sxp = sqrt( max( s0p - sy(ij, k, n) * sy(ij, k, n) * ci3, &
             &              eps ) )
               
           if ( abs( sy(ij, k, n) ) < 1.5d0 * s0m ) then

!             syy(ij, k, n) = min( s0m + sxp, &
!               &                  max( abs( sy(ij, k, n) ) - s0m, &
!               &                       syy(ij, k, n) ) )
              syy(ij, k, n) = min( s0m - abs( sy(ij, k, n) ),  &
                &                  max( - s0m - sxp,           &
                &                       syy(ij, k, n) ) )

           else
                  
!             syy(ij, k, n) = min( s0m + sxp, &
!               &                  max( s0m - sxp, syy(ij, k, n) ) )
              syy(ij, k, n) = min( - s0m + sxp,           &
                &                  max( - s0m - sxp,      &
                &                       syy(ij, k, n) ) )

           end if
             
           sxy(ij, k, n) = min( s0m, max( - s0m, sxy(ij, k, n) ) )
           syz(ij, k, n) = min( s0m, max( - s0m, syz(ij, k, n) ) )

        end do
     end do

!---- bug fix 2
#ifdef OPT_TRIPOLE
     call shift1( sy(:,:,n), &
       &           nxdim,  nydim,  nzdim, &
       &           -1.d0,      0,      0 )
     call shift1( syy(:,:,n), &
       &           nxdim,  nydim,  nzdim, &
       &            1.d0,      0,      0 )
     call shift1( sxy(:,:,n), &
       &           nxdim,  nydim,  nzdim, &
       &           -1.d0,      0,      0 )
     call shift1( syz(:,:,n), &
       &           nxdim,  nydim,  nzdim, &
       &           -1.d0,      0,      0 )
#else
     call shift1( sy (:,:,n), nxdim, nydim, nzdim )
     call shift1( syy(:,:,n), nxdim, nydim, nzdim )
     call shift1( sxy(:,:,n), nxdim, nydim, nzdim )
     call shift1( syz(:,:,n), nxdim, nydim, nzdim )
#endif

!    ---- calculating ALF  and MASS between box (i,j-1,k) <---> (i,j,k)
     do k = kstr, kend
        do ij = ijtstr-nxdim-1, ijtend+nxdim+1

           ijls  = ij + ls

           if ( uv(ij, k) .gt. 0.d0 ) then
              fm(ijls, k) = uv(ij, k) * ts
              alf(ij, k)  = fm(ijls, k) / sm(ijls, k, n)
           else
              fm(ijls, k) = - uv(ij, k) * ts
              alf(ij, k)  = fm(ijls, k) / sm(ij, k, n)
           end if

        end do
     end do

!    ---- calculating flux between box (i,j-1,k) <---> (i,j,k)
     do k = kstr, kend
        do ij = ijtstr-nxdim-1, ijtend+nxdim+1

           ijls  = ij + ls

!          ---- flux from (j-1) to (j),  when v > 0
           if ( uv(ij, k) .gt. 0.d0 ) then

              alfq  = alf(ij, k) * alf(ij, k)
              alf1  = 1.d0 - alf(ij, k)
              alf1q = alf1 * alf1

!             ---- moments be transported
              f0(ijls, k) =  &
                &   alf(ij, k) * ( s0(ijls, k, n) &
                &       + alf1 * ( &
                &          sy(ijls, k, n) &
                &        + ( alf1 - alf(ij, k) ) * syy(ijls, k, n) ) )

              fy(ijls, k) = &
                &   alfq * ( sy(ijls, k, n) &
                &          + 3.d0 * alf1 * syy(ijls, k, n) ) 

              fyy(ijls, k) = alf(ij, k) * alfq * syy(ijls, k, n)

              fx (ijls, k) = alf(ij, k) * ( sx(ijls, k, n) &
                &                         + alf1 * sxy(ijls, k, n) )
              fz (ijls, k) = alf(ij, k) * ( sz(ijls, k, n) &
                &                         + alf1 * syz(ijls, k, n) )

              fxy(ijls, k) = alfq       * sxy(ijls, k, n)
              fyz(ijls, k) = alfq       * syz(ijls, k, n)
              fxx(ijls, k) = alf(ij, k) * sxx(ijls, k, n)
              fzz(ijls, k) = alf(ij, k) * szz(ijls, k, n)
              fxz(ijls, k) = alf(ij, k) * sxz(ijls, k, n)
             
              fty (ij, k, n) = fty(ij, k, n) &
                &            - f0(ijls, k) * tsiv * rx

           else
             
!             ---- flux from (j) to (j-1),  when v < 0

              alfq  = alf(ij, k) * alf(ij, k)
              alf1  = 1.d0 - alf(ij, k)
              alf1q = alf1 * alf1
             
              f0 (ijls, k) = &
                &    alf(ij, k) * ( s0(ij, k, n) &
                &          - alf1 * (  sy(ij, k, n) &
                &                    - ( alf1 - alf(ij, k) ) &
                &                     * syy(ij, k, n) ) )
              fy (ijls, k) = alfq * ( &
                &    sy(ij, k, n) - 3.d0 * alf1 * syy(ij, k, n) )
              fyy(ijls, k) = alf(ij, k) * alfq * syy(ij, k, n)
                  
              fx (ijls, k) = alf(ij, k) * ( sx(ij, k, n) &
                &                         - alf1 * sxy(ij, k, n) )
              fz (ijls, k) = alf(ij, k) * ( sz(ij, k, n) &
                &                         - alf1 * syz(ij, k, n) )

              fxy(ijls, k) = alfq       * sxy(ij, k, n)
              fyz(ijls, k) = alfq       * syz(ij, k, n)
              fxx(ijls, k) = alf(ij, k) * sxx(ij, k, n)
              fzz(ijls, k) = alf(ij, k) * szz(ij, k, n)
              fxz(ijls, k) = alf(ij, k) * sxz(ij, k, n)
             
              fty (ij, k, n) = fty(ij, k, n) &
                &            + f0(ijls, k) * tsiv * rx

           end if
          
        end do
     end do

!    ---- calculating moments
     do k = kstr, kend
        do ij = ijtstr-nxdim-1, ijtend+nxdim+1

           ijls  = ij + ls

!          ---- flux from (j-1) to (j),  when v > 0
           if ( uv(ij, k) .gt. 0.d0 ) then

              alfq  = alf(ij, k) * alf(ij, k)
              alf1  = 1.d0 - alf(ij, k)
              alf1q = alf1 * alf1

!             ---- moments remaining 
              sm (ijls, k, n) = sm(ijls, k, n) - fm(ijls, k)
              s0 (ijls, k, n) = s0(ijls, k, n) - f0(ijls, k)

              sy (ijls, k, n) = &
                &     alf1q * ( sy(ijls, k, n) &
                &             - 3.d0 * alf(ij, k) * syy(ijls, k, n) )
              syy(ijls, k, n) = alf1 * alf1q * syy(ijls, k, n)

              sx (ijls, k, n) = sx (ijls, k, n) - fx (ijls, k)
              sxx(ijls, k, n) = sxx(ijls, k, n) - fxx(ijls, k)

              sz (ijls, k, n) = sz (ijls, k, n) - fz (ijls, k)
              szz(ijls, k, n) = szz(ijls, k, n) - fzz(ijls, k)

              sxy(ijls, k, n) = alf1q * sxy(ijls, k, n)
              syz(ijls, k, n) = alf1q * syz(ijls, k, n)
              sxz(ijls, k, n) = sxz(ijls, k, n) - fxz(ijls, k)
             
           else
             
!             ---- flux from (j) to (j-1),  when v < 0

              alfq  = alf(ij, k) * alf(ij, k)
              alf1  = 1.d0 - alf(ij, k)
              alf1q = alf1 * alf1
             
              sm (ij, k, n) = sm(ij, k, n) - fm(ijls, k)
              s0 (ij, k, n) = s0(ij, k, n) - f0(ijls, k)

              sy (ij, k, n) = alf1q * ( &
                &         sy(ij, k, n) &
                &       + 3.d0 * alf(ij, k) * syy(ij, k, n)  )
              syy(ij, k, n) = alf1 * alf1q * syy(ij, k, n)

              sx (ij, k, n) = sx (ij, k, n) - fx (ijls, k)
              sxx(ij, k, n) = sxx(ij, k, n) - fxx(ijls, k)

              sz (ij, k, n) = sz (ij, k, n) - fz (ijls, k)
              szz(ij, k, n) = szz(ij, k, n) - fzz(ijls, k)

              sxy(ij, k, n) = alf1q * sxy(ij, k, n)
              syz(ij, k, n) = alf1q * syz(ij, k, n)
              sxz(ij, k, n) = sxz(ij, k, n) - fxz(ijls, k)
             
           end if

        end do
     end do

!    ---- put the temporary moments (fi) into appropriate neighboring boxes
     do k = kstr, kend
        do ij = ijtstr-nxdim-1, ijtend+nxdim+1

           ijls  = ij + ls

           if ( uv(ij, k) .gt. 0.d0 ) then
              sm(ij, k, n) = sm(ij, k, n) + fm(ijls, k)
              alf(ij, k)   = fm(ijls, k) / sm(ij, k, n)
           else
              sm(ijls, k, n) = sm(ijls, k, n) + fm(ijls, k)
              alf(ij, k)     = fm(ijls, k) / sm(ijls, k, n)
           end if

        end do
     end do

     do k = kstr, kend
        do ij = ijtstr-nxdim-1, ijtend+nxdim+1

           ijls  = ij + ls

!          ---- flux from (j-1) to (j),  when v > 0
           if ( uv(ij, k) .gt. 0.d0 ) then

              alf1 = 1.d0 - alf(ij, k)

              tmp = alf(ij, k) * s0(ij, k, n) - alf1 * f0(ijls, k)

              s0(ij, k, n) = s0(ij, k, n) + f0(ijls, k)
              syy(ij, k, n) = &
                &    alf(ij, k) * alf(ij, k) * fyy(ijls, k) &
                &  + alf1 * alf1 * syy(ij, k, n) &
                &  + 5.d0 * (  alf(ij, k) * alf1 * ( &
                &                sy(ij, k, n) - fy(ijls, k) ) &
                &            - ( alf1 - alf(ij, k) ) * tmp )
              sy (ij, k, n) = &
                &    alf(ij, k) * fy(ijls, k) &
                &  + alf1 * sy(ij, k, n) + 3.d0 * tmp

              sxy(ij, k, n) = &
                &    alf(ij, k) * fxy(ijls, k) &
                &  + alf1 * sxy(ij, k, n) &
                &  + 3.d0 * ( alf(ij, k) * sx(ij, k, n) &
                &           - alf1 * fx(ijls, k) )
              syz(ij, k, n) = &
                &    alf(ij, k) * fyz(ijls, k) &
                &  + alf1 * syz(ij, k, n) &
                &  + 3.d0 * ( alf(ij, k) * sz(ij, k, n) &
                &           - alf1 * fz(ijls, k) )

              sx (ij, k, n) = sx (ij, k, n) + fx (ijls, k)
              sxx(ij, k, n) = sxx(ij, k, n) + fxx(ijls, k)

              sz (ij, k, n) = sz (ij, k, n) + fz (ijls, k)
              szz(ij, k, n) = szz(ij, k, n) + fzz(ijls, k)

              sxz(ij, k, n) = sxz(ij, k, n) + fxz(ijls, k)

!          ---- flux from (j) to (j-1),  when v < 0
           else

              alf1 = 1.d0 - alf(ij, k)

              tmp = - alf(ij, k) * s0(ijls, k, n) &
                &   + alf1 * f0(ijls, k)

              s0(ijls, k, n) = s0(ijls, k, n) + f0(ijls, k)
              syy(ijls, k, n) = &
                &      alf(ij, k) * alf(ij, k) * fyy(ijls, k) &
                &    + alf1 * alf1 * syy(ijls, k, n) &
                &    + 5.d0 * (  alf(ij, k) * alf1 * ( &
                &              - sy(ijls, k, n) + fy(ijls, k) ) &
                &              + ( alf1 - alf(ij, k) ) * tmp  )
              sy (ijls, k, n) = &
                &      alf(ij, k) * fy(ijls, k) + alf1 * sy(ijls, k, n) &
                &    + 3.d0 * tmp

              sxy(ijls, k, n) = &
                &      alf(ij, k) * fxy(ijls, k) &
                &    + alf1 * sxy(ijls, k, n) &
                &    + 3.d0 * ( alf1       * fx(ijls, k) &
                &             - alf(ij, k) * sx(ijls, k, n) )
              syz(ijls, k, n) = &
                &      alf(ij, k) * fyz(ijls, k) &
                &    + alf1 * syz(ijls, k, n) &
                &    + 3.d0 * ( alf1       * fz(ijls, k) &
                &             - alf(ij, k) * sz(ijls, k, n) )

              sx (ijls, k, n) = sx (ijls, k, n) + fx (ijls, k)
              sz (ijls, k, n) = sz (ijls, k, n) + fz (ijls, k)
              sxx(ijls, k, n) = sxx(ijls, k, n) + fxx(ijls, k)
              szz(ijls, k, n) = szz(ijls, k, n) + fzz(ijls, k)
              sxz(ijls, k, n) = sxz(ijls, k, n) + fxz(ijls, k)

           end if

        end do
     end do

  end do

#ifdef OPT_BBL
! ---- keeping BBL variables consistent at two levels
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

! ---- Z-direction

  do k = kstr, kend
     do ij = ijtstr, ijtend

!        uv(ij, k) = - wzc(ij, k) * vlmz(ij)
        uv(ij, k) = - wzc(ij, k) * vlmz(ij) &
          &         * amskt(ij, k) * amskt(ij, k-1)
               
     end do
  end do

#ifdef OPT_BBL
  do ij = ijtstr, ijtend
     k = nbot(ij)
     uv(ij, k) = - wzc(ij, k) * vlmz(ij) * amsktb(ij)
  end do
#endif

  do n = 1, ntdim

!    ---- undershoot limiter (Method B) of Morales Maqueda and Holloway (2006)
     do k = kstr, kend

        ku = max( k - 1, kstr )
        kd = min( k + 1, kend )
            
        do ij = ijtstr, ijtend

           s0m = s0(ij, k, n) &
             & - min( s0(ij, ku, n) / sm(ij, ku, n), &
             &        s0(ij, k , n) / sm(ij, k , n), &
             &        s0(ij, kd, n) / sm(ij, kd, n)  ) & 
             & * sm(ij, k, n)

           s1m = sq3 * s0m
           sz(ij, k, n) = min( s1m, max( - s1m, sz(ij, k, n) ) )
           s0p = s0m * s0m
           sxp = sqrt( max( s0p - sz(ij, k, n) * sz(ij, k, n) * ci3, &
             &              eps ) )
               
           if ( abs( sz(ij, k, n) ) < 1.5d0 * s0m ) then

              szz(ij, k, n) = min( s0m + sxp, &
                &                  max( abs( sz(ij, k, n) ) - s0m, &
                &                       szz(ij, k, n) ) )

           else
                  
              szz(ij, k, n) = min( s0m + sxp, &
                &                  max( s0m - sxp, szz(ij, k, n) ) )

           end if
             
           sxz(ij, k, n) = min( s0m, max( - s0m, sxz(ij, k, n) ) )
           syz(ij, k, n) = min( s0m, max( - s0m, syz(ij, k, n) ) )

        end do
     end do

!    ---- overshoot limiter (Method B) of Morales Maqueda and Holloway (2006)
     do k = kstr, kend

        ku = max( k - 1, kstr )
        kd = min( k + 1, kend )

        do ij = ijtstr, ijtend

           s0m = - s0(ij, k, n) &
             &   + max( s0(ij, ku, n) / sm(ij, ku, n), &
             &          s0(ij, k,  n) / sm(ij, k , n), &
             &          s0(ij, kd, n) / sm(ij, kd, n)  ) &
             &   * sm(ij, k, n)

           s1m = sq3 * s0m
           sz(ij, k, n) = min( s1m, max( - s1m, sz(ij, k, n) ) )
           s0p = s0m * s0m
           sxp = sqrt( max( s0p - sz(ij, k, n) * sz(ij, k, n) * ci3, &
             &              eps ) )
               
           if ( abs( sz(ij, k, n) ) < 1.5d0 * s0m ) then

!             szz(ij, k, n) = min( s0m + sxp, &
!               &                  max( abs( sz(ij, k, n) ) - s0m, &
!               &                       szz(ij, k, n) ) )
              szz(ij, k, n) = min( s0m - abs( sz(ij, k, n) ),  &
                &                  max( - s0m - sxp,           &
                &                       szz(ij, k, n) ) )

           else
                  
!             szz(ij, k, n) = min( s0m + sxp, &
!               &                  max( s0m - sxp, szz(ij, k, n) ) )
              szz(ij, k, n) = min( - s0m + sxp,          &
     &                             max( - s0m - sxp,     &
     &                                  szz(ij, k, n) ) )

           end if
             
           sxz(ij, k, n) = min( s0m, max( - s0m, sxz(ij, k, n) ) )
           syz(ij, k, n) = min( s0m, max( - s0m, syz(ij, k, n) ) )
           
        end do
     end do

!    ---- calculating ALF
     do k = kstr, kend

        ku = k - 1

        do ij = ijtstr, ijtend

           if ( uv(ij, k) .gt. 0.d0 ) then
              fm(ij, ku) = uv(ij, k) * ts
              alf(ij, k) = fm(ij, ku) / sm(ij, ku, n)
           else
              fm(ij, ku) = - uv(ij, k) * ts
              alf(ij, k) = fm(ij, ku) / sm(ij, k, n)
           end if

        end do
     end do

!    ---- calculating flux between box (i,j,k-1) <---> (i,j,k)
     do k = kstr, kend

        ku = k - 1

        do ij = ijtstr, ijtend

           if ( uv(ij, k) .gt. 0.d0 ) then

!             ---- flux from (k-1) to (k),  when w > 0
              alfq  = alf(ij, k) * alf(ij, k)
              alf1  = 1.d0 - alf(ij, k)
              alf1q = alf1 * alf1

!             ---- moments be transported
              f0(ij, ku) = &
                &   alf(ij, k) * ( s0(ij, ku, n) &
                &         + alf1 * ( &
                &              sz(ij, ku, n) &
                &            + ( alf1 - alf(ij, k) ) &
                &            * szz(ij, ku, n) ) )

              fz(ij, ku) = &
                &   alfq * ( sz(ij, ku, n) &
                &          + 3.d0 * alf1 * szz(ij, ku, n) ) 

              fzz(ij, ku) = alf(ij, k) * alfq * szz(ij, ku, n)

              fx (ij, ku) = alf(ij, k) * ( sx(ij, ku, n) &
                &                        + alf1 * sxz(ij, ku, n) )
              fy (ij, ku) = alf(ij, k) * ( sy(ij, ku, n) &
                &                        + alf1 * syz(ij, ku, n) )

              fxz(ij, ku) = alfq       * sxz(ij, ku, n)
              fyz(ij, ku) = alfq       * syz(ij, ku, n)
              fxx(ij, ku) = alf(ij, k) * sxx(ij, ku, n)
              fyy(ij, ku) = alf(ij, k) * syy(ij, ku, n)
              fxy(ij, ku) = alf(ij, k) * sxy(ij, ku, n)
                  
              ftz (ij, k, n) = ftz(ij, k, n) &
                &            + f0(ij, ku) * tsiv / vlmz(ij)

           else
             
!             ---- flux from (k) to (k-1),  when w < 0

              alfq  = alf(ij, k) * alf(ij, k)
              alf1  = 1.d0 - alf(ij, k)
              alf1q = alf1 * alf1
             
              f0 (ij, ku) = &
                &    alf(ij, k) * ( s0(ij, k, n) &
                &          - alf1 * (  sz(ij, k, n) &
                &                     - ( alf1 - alf(ij, k) ) &
                &                     * szz(ij, k, n) ) )
              fz (ij, ku) = alfq * ( &
                &    sz(ij, k, n) - 3.d0 * alf1 * szz(ij, k, n) )
              fzz(ij, ku) = alf(ij, k) * alfq * szz(ij, k, n)
                  
              fx (ij, ku) = alf(ij, k)  * ( sx(ij, k, n) &
                &                         - alf1 * sxz(ij, k, n) )
              fy (ij, ku) = alf(ij, k)  * ( sy(ij, k, n) &
                &                         - alf1 * syz(ij, k, n) )

              fxz(ij, ku) = alfq       * sxz(ij, k, n)
              fyz(ij, ku) = alfq       * syz(ij, k, n)
              fxx(ij, ku) = alf(ij, k) * sxx(ij, k, n)
              fyy(ij, ku) = alf(ij, k) * syy(ij, k, n)
              fxy(ij, ku) = alf(ij, k) * sxy(ij, k, n)
             
              ftz (ij, k, n) = ftz(ij, k, n) &
                &            - f0(ij, ku) * tsiv / vlmz(ij)

           end if
          
        end do
     end do

!    ---- calculating flux between box (i,j,k-1) <---> (i,j,k)
     do k = kstr, kend

        ku = k - 1

        do ij = ijtstr, ijtend

           if ( uv(ij, k) .gt. 0.d0 ) then

!             ---- flux from (k-1) to (k),  when w > 0
              alfq  = alf(ij, k) * alf(ij, k)
              alf1  = 1.d0 - alf(ij, k)
              alf1q = alf1 * alf1

!             ---- moments remaining 
              sm (ij, ku, n) = sm(ij, ku, n) - fm(ij, ku)
              s0 (ij, ku, n) = s0(ij, ku, n) - f0(ij, ku)

              sz (ij, ku, n) = &
                &     alf1q * ( sz(ij, ku, n) &
                &             - 3.d0 * alf(ij, k) * szz(ij, ku, n) )
              szz(ij, ku, n) = alf1 * alf1q * szz(ij, ku, n)

              sx (ij, ku, n) = sx (ij, ku, n) - fx (ij, ku)
              sxx(ij, ku, n) = sxx(ij, ku, n) - fxx(ij, ku)

              sy (ij, ku, n) = sy (ij, ku, n) - fy (ij, ku)
              syy(ij, ku, n) = syy(ij, ku, n) - fyy(ij, ku)

              sxz(ij, ku, n) = alf1q * sxz(ij, ku, n)
              syz(ij, ku, n) = alf1q * syz(ij, ku, n)
              sxy(ij, ku, n) = sxy(ij, ku, n) - fxy(ij, ku)
             
           else
             
!             ---- flux from (k) to (k-1),  when w < 0

              alfq  = alf(ij, k) * alf(ij, k)
              alf1  = 1.d0 - alf(ij, k)
              alf1q = alf1 * alf1
             
              sm (ij, k, n) = sm(ij, k, n) - fm(ij, ku)
              s0 (ij, k, n) = s0(ij, k, n) - f0(ij, ku)

              sz (ij, k, n) = alf1q * ( &
                &     sz(ij, k, n) &
                &   + 3.d0 * alf(ij, k) * szz(ij, k, n)  )
              szz(ij, k, n) = alf1 * alf1q * szz(ij, k, n)

              sx (ij, k, n) = sx (ij, k, n) - fx (ij, ku)
              sxx(ij, k, n) = sxx(ij, k, n) - fxx(ij, ku)

              sy (ij, k, n) = sy (ij, k, n) - fy (ij, ku)
              syy(ij, k, n) = syy(ij, k, n) - fyy(ij, ku)

              sxz(ij, k, n) = alf1q * sxz(ij, k, n)
              syz(ij, k, n) = alf1q * syz(ij, k, n)
              sxy(ij, k, n) = sxy(ij, k, n) - fxy(ij, ku)
             
           end if
          
        end do
     end do

!    ---- put the temporary moments (fi) into appropriate neighboring boxes
     do k = kstr, kend

        ku = k - 1

        do ij = ijtstr, ijtend

           if ( uv(ij, k) .gt. 0.d0 ) then
              sm(ij, k, n) = sm(ij, k, n) + fm(ij, ku)
              alf(ij, k)   = fm(ij, ku) / sm(ij, k, n)
           else
              sm(ij, ku, n) = sm(ij, ku, n) + fm(ij, ku)
              alf(ij, k)    = fm(ij, ku) / sm(ij, ku, n)
           end if

        end do
     end do

     do k = kstr, kend

        ku = k - 1

        do ij = ijtstr, ijtend

!          ---- flux from (k-1) to (k),  when w > 0
           if ( uv(ij, k) .gt. 0.d0 ) then

              alf1 = 1.d0 - alf(ij, k)

              tmp = alf(ij, k) * s0(ij, k, n) - alf1 * f0(ij, ku)

              s0(ij, k, n) = s0(ij, k, n) + f0(ij, ku)

              szz(ij, k, n) = &
                &    alf(ij, k) * alf(ij, k) * fzz(ij, ku) &
                &  + alf1 * alf1 * szz(ij, k, n) &
                &  + 5.d0 * (  alf(ij, k) * alf1 * ( &
                &                sz(ij, k, n) - fz(ij, ku) ) &
                &            - ( alf1 - alf(ij, k) ) * tmp )
              sz (ij, k, n) = &
                &    alf(ij, k) * fz(ij, ku) &
                &  + alf1       * sz(ij, k, n) + 3.d0 * tmp

              sxz(ij, k, n) = &
                &    alf(ij, k) * fxz(ij, ku) &
                &  + alf1       * sxz(ij, k, n) &
                &  + 3.d0 * ( alf(ij, k) * sx(ij, k, n) &
                &           - alf1 * fx(ij, ku) )
              syz(ij, k, n) = &
                &    alf(ij, k) * fyz(ij, ku) &
                &  + alf1       * syz(ij, k, n) &
                &  + 3.d0 * ( alf(ij, k) * sy(ij, k, n) &
                &           - alf1 * fy(ij, ku) )

              sx (ij, k, n) = sx (ij, k, n) + fx (ij, ku)
              sxx(ij, k, n) = sxx(ij, k, n) + fxx(ij, ku)

              sy (ij, k, n) = sy (ij, k, n) + fy (ij, ku)
              syy(ij, k, n) = syy(ij, k, n) + fyy(ij, ku)

              sxy(ij, k, n) = sxy(ij, k, n) + fxy(ij, ku)

!          ---- flux from (k) to (k-1),  when w < 0
           else

              alf1 = 1.d0 - alf(ij, k)

              tmp = - alf(ij, k) * s0(ij, ku, n) + alf1 * f0(ij, ku)

              s0(ij, ku, n) = s0(ij, ku, n) + f0(ij, ku)

              szz(ij, ku, n) = &
                &      alf(ij, k) * alf(ij, k) * fzz(ij, ku) &
                &    + alf1 * alf1 * szz(ij, ku, n) &
                &    + 5.d0 * (  alf(ij, k) * alf1 * ( - sz(ij, ku, n) &
                &                                      + fz(ij, ku)    ) &
                &              + ( alf1 - alf(ij, k) ) * tmp  )
              sz (ij, ku, n) = &
                &      alf(ij, k) * fz(ij, ku) + alf1 * sz(ij, ku, n) &
                &    + 3.d0 * tmp

              sxz(ij, ku, n) = &
                &      alf(ij, k) * fxz(ij, ku) + alf1 * sxz(ij, ku, n) &
                &    + 3.d0 * ( alf1       * fx(ij, ku) &
                &             - alf(ij, k) * sx(ij, ku, n) )
              syz(ij, ku, n) = &
                &      alf(ij, k) * fyz(ij, ku) + alf1 * syz(ij, ku, n) &
                &    + 3.d0 * ( alf1       * fy(ij, ku) &
                &             - alf(ij, k) * sy(ij, ku, n) )

              sx (ij, ku, n) = sx (ij, ku, n) + fx (ij, ku)
              sy (ij, ku, n) = sy (ij, ku, n) + fy (ij, ku)
              sxx(ij, ku, n) = sxx(ij, ku, n) + fxx(ij, ku)
              syy(ij, ku, n) = syy(ij, ku, n) + fyy(ij, ku)
              sxy(ij, ku, n) = sxy(ij, ku, n) + fxy(ij, ku)

           end if

        end do
     end do

  end do

! ---- tx
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
!  call stbbtr( s0 )
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
  &     ty,     tx )

  use xprst
  use ufile

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

  real(8) ::   dtdx(nxydim, nzdim, ntdim),   dtdy(nxydim, nzdim, ntdim)
  real(8) ::   dtdz(nxydim, nzdim, ntdim)
  real(8) ::     p1,     p2
  real(8) ::     tl,     sl
  real(8) ::     rl,    rlw,    rls,    rlu
  real(8) ::   dzdx,   dzdy
  integer ::     ij,      k,      n
  integer ::  ifpar,  jfpar,  istat

  real(8), save ::   slpz(nzdim)
  real(8) ::   dzm0(nzdim)

  real(8), save :: slpmax = 1.d-2
  integer :: kslp = -999
  
  namelist /nmslpm/ slpmax, kslp

  if (ofirst) then
     ofirst = .false.
     call rewnml(ifpar, jfpar)
     read(ifpar, nmslpm, iostat=istat)
     call cstnml(jfpar, 'dnsgrd', 'nmslpm', istat)
     write(jfpar, nmslpm)

     call secoef( &
        &   c0(kstr), c1(kstr), c2(kstr), c3(kstr), &
        &   c4(kstr), c5(kstr), c6(kstr), &
        &   d0(kstr), d1(kstr), d2(kstr), d3(kstr),  d4(kstr), &
        &   d5(kstr), d6(kstr), d7(kstr), d8(kstr),  d9(kstr))

     if ( kslp < 0 ) then

        do k = 1, nzdim
           slpz(k) = slpmax
        end do

     else

!       ---- vertically varying slpmax (like OPA)
        do k = 1, kstr-1
           dzm0(k) = 0.5d0 * dz0(k)
        end do
        do k = kstr, nzdim
           dzm0(k) = 0.5d0 * (dz0(k-1) + dz0(k))
        end do

        p1 = 0.d0
        do k = kstr, kstr+kslp-1
           p1 = p1 + dzm0(k)
        end do

        p2 = 0.d0
        do k = kstr, kstr+kslp-1
           p2 = p2 + dzm0(k)
!           slpz(k) = slpmax * exp( - 10.d0 * ( p1 - p2 ) / p1 )
           slpz(k) = ( slpmax - 0.d0 ) / ( p1 - 0.d0 ) * p2
        end do

        do k = 1, kstr-1
           slpz(k) = slpz(kstr)
        end do

        do k = kstr+kslp, nzdim
           slpz(k) = slpz(kstr+kslp-1)
        end do

        write(jfpar, *) 'depth dependent Maximum slope:'
        do k = 1, nzdim
           write(jfpar, '(i6,f10.3)') k, slpz(k)
        end do

     end if

  end if

  do n = 1, ntdim
     do k = 1, nzdim
        do ij = 1, nxydim
           dtdx(ij, k, n) = 0.d0
           dtdy(ij, k, n) = 0.d0
           dtdz(ij, k, n) = 0.d0
        end do
     end do
  end do

  do k = kstr, kend
     do ij = ijtstr, ijtend+nxdim
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

        tl = ty(ij+lw, k, 1)
        sl = ty(ij+lw, k, 2)
        p1 = c0(k) &
           &   + (c1(k) + (c2(k) + c3(k) * tl) * tl) * tl &
           &   + (c4(k) + c5(k) * tl + c6(k) * sl) * sl
        p2 = d0(k) &
           &   + (d1(k) + (d2(k) + (d3(k) + d4(k) * tl) * tl) * tl) * tl &
           &   + (d5(k) + (d6(k) + d7(k) * tl * tl) * tl &
           &            + (d8(k) + d9(k) * tl * tl) * sqrt(sl)) * sl
        rlw = p1 / p2

        tl = ty(ij+ls, k, 1)
        sl = ty(ij+ls, k, 2)
        p1 = c0(k) &
           &   + (c1(k) + (c2(k) + c3(k) * tl) * tl) * tl &
           &   + (c4(k) + c5(k) * tl + c6(k) * sl) * sl
        p2 = d0(k) &
           &   + (d1(k) + (d2(k) + (d3(k) + d4(k) * tl) * tl) * tl) * tl &
           &   + (d5(k) + (d6(k) + d7(k) * tl * tl) * tl &
           &            + (d8(k) + d9(k) * tl * tl) * sqrt(sl)) * sl
        rls = p1 / p2

        dtdx(ij, k, 1) = (rl - rlw) * rx * 2.d0 &
           &             / (hxt(ij) + hxt(ij+lw)) * &
           &             amskt(ij, k) * amskt(ij+lw, k)
        dtdy(ij, k, 1) = (rl - rls) * rym(ij+ls) * 2.d0 &
           &             / (hyt(ij) + hyt(ij+ls)) * &
           &             amskt(ij, k) * amskt(ij+ls, k)
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

        dtdz(ij, k, 1) = min((rlu - rl) / dzm(ij, k), 0.d0) * &
           &             amftz(ij, k)
     end do
  end do

  do k = kstr, kend
     do ij = ijtstr, ijtend+nxdim
        dzdx = dtdx(ij, k, 1) * 4.d0 &
           &   / (  dtdz(ij, k  , 1) + dtdz(ij+lw, k ,  1) &
           &      + dtdz(ij, k+1, 1) + dtdz(ij+lw, k+1, 1) - eps)
        xdzdx(ij, k) = min(slpz(k), max(-slpz(k), dzdx))
        dzdy = dtdy(ij, k, 1) * 4.d0 * dym(ij+ls) &
           &   / (  (dtdz(ij, k, 1) + dtdz(ij, k+1, 1)) * dy(ij) &
           &   + (dtdz(ij+ls, k, 1) + dtdz(ij+ls, k+1, 1)) * dy(ij+ls) &
           &      - eps)
        ydzdy(ij, k) = min(slpz(k), max(-slpz(k), dzdy))
     end do
  end do

  do k = kstr+1, kend
     do ij = ijtstr, ijtend
        dzdx = ((dtdx(ij, k-1, 1)+dtdx(ij+le, k-1, 1))*dz(ij, k-1) &
           &  + (dtdx(ij, k  , 1)+dtdx(ij+le, k  , 1))*dz(ij, k)) * &
           &   0.25d0 / dzm(ij, k) / (dtdz(ij, k, 1) - eps)
        zdzdx(ij, k) = min(slpz(k), max(-slpz(k), dzdx)) * &
           &           amftz(ij, k)
        dzdy = &
           &   (  (dtdy(ij, k-1, 1)+dtdy(ij+ln, k-1, 1))*dz(ij, k-1) &
           &    + (dtdy(ij, k  , 1)+dtdy(ij+ln, k  , 1))*dz(ij, k)) * &
           &     0.25d0 / dzm(ij, k) / (dtdz(ij, k, 1) - eps)
        zdzdy(ij, k) = min(slpz(k), max(-slpz(k), dzdy)) * &
           &           amftz(ij, k)
     end do
  end do

  do n = 1, ntdim
     do k = 1, nzdim
        do ij = 1, nxydim
           dtdx(ij, k, n) = 0.d0
           dtdy(ij, k, n) = 0.d0
           dtdz(ij, k, n) = 0.d0
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
           dtdz(ij, k, n) = (tx(ij, k-1, n) - tx(ij, k, n)) &
              &             / dzm(ij, k) * amftz(ij, k)
        end do
     end do
  end do

  do n = 1, ntdim
     do k = kstr, kend
        do ij = ijtstr, ijtend+nxdim
           xdtdz(ij, k, n) = &
              &   (  dtdz(ij, k, n)   + dtdz(ij+lw, k, n) &
              &    + dtdz(ij, k+1, n) + dtdz(ij+lw, k+1, n)) *  0.25d0
           ydtdz(ij, k, n) = &
              &   (  (dtdz(ij, k, n)    + dtdz(ij, k+1, n)) * &
              &      dy(ij) &
              &    + (dtdz(ij+ls, k, n) + dtdz(ij+ls, k+1, n)) * &
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
     &            'ocean zonal heat flux', 'degC cm^3/rad/s', &
     &            nx,     ny,     nz, nxyzdm, 'OCLVTT')
  call chekin(   fty,  'FTY', &
     &            'ocean meridional heat flux', 'degC cm^3/rad/s', &
     &            nx,     ny,     nz, nxyzdm, 'OCLVTT')
  call chekin(   ftz,  'FTZ', &
     &            'ocean vertical heat flux', 'degC cm/s', &
     &            nx,     ny,     nz, nxyzdm, 'OCLVTT')
  call chekin(   ftx(1, 1, 2),  'FSX', &
     &            'ocean zonal salt flux', 'psu cm^3/rad/s', &
     &            nx,     ny,     nz, nxyzdm, 'OCLVTT')
  call chekin(   fty(1, 1, 2),  'FSY', &
     &            'ocean meridional salt flux', 'psu cm^3/rad/s', &
     &            nx,     ny,     nz, nxyzdm, 'OCLVTT')
  call chekin(   ftz(1, 1, 2),  'FSZ', &
     &            'ocean vertical salt flux', 'psu cm/s', &
     &            nx,     ny,     nz, nxyzdm, 'OCLVTT')

  return
end subroutine chkftx

end module tflxt
