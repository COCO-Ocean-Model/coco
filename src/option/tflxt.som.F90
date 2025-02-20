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
    &  ijstr,  ijend, &
    &  igstr,  jgstr, &
    &     le,     lw,     ln,     ls, &
    &    lsw, &
    &  oinit, ofinal, &
    & myrank,  iroot
  use zocgrd, only: &
    &     dy,    dym,     dz,    dz0,    dzm,    dzv,     ds,    dsm, &
    &     dx,     rx,     ry,    rym, &
    &     ts,   zbot, &
    &    cor, &
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
  use zocphy, only: &
    & gravit,   rhoo

  implicit none
  private

  real(8), save :: ftx(nxydim, nzdim, ntdim)
  real(8), save :: fty(nxydim, nzdim, ntdim)
  real(8), save :: ftz(nxydim, nzdim, ntdim)
! --- for diffusive flux
  real(8), save :: ftxd(nxydim, nzdim, ntdim)
  real(8), save :: ftyd(nxydim, nzdim, ntdim)
  real(8), save :: ftzd(nxydim, nzdim, ntdim)

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
  real(8) :: fharmx(nxydim), fharmy(nxydim),   harm(nxydim)
  real(8) ::  hzbot(nxydim)
  real(8) ::     dh(nxydim)

  real(8) ::  xdzdx(nxydim, nzdim),  ydzdy(nxydim, nzdim)
  real(8) ::  zdzdx(nxydim, nzdim),  zdzdy(nxydim, nzdim)
  real(8) ::  xdtdz(nxydim, nzdim, ntdim),  ydtdz(nxydim, nzdim, ntdim)
  real(8) ::  zdtdx(nxydim, nzdim, ntdim),  zdtdy(nxydim, nzdim, ntdim)

! --- for flux output
  real(8) ::   adt2(nxydim, nzdim, ntdim)    
  real(8) ::   adtd(nxydim, nzdim, ntdim)    
  real(8) ::  adtah(nxydim, nzdim, ntdim)    
  real(8) ::  adtgm(nxydim, nzdim, ntdim)    
  real(8) ::  adtis(nxydim, nzdim, ntdim)    
  real(8) ::   ftx2(nxydim, nzdim, ntdim)
  real(8) ::   fty2(nxydim, nzdim, ntdim)
  real(8) ::   ftz2(nxydim, nzdim, ntdim)
  real(8) ::  ftxah(nxydim, nzdim, ntdim)
  real(8) ::  ftyah(nxydim, nzdim, ntdim)
  real(8) ::  ftxgm(nxydim, nzdim, ntdim)
  real(8) ::  ftygm(nxydim, nzdim, ntdim)
  real(8) ::  ftzgm(nxydim, nzdim, ntdim)
  real(8) ::  ftxis(nxydim, nzdim, ntdim)
  real(8) ::  ftyis(nxydim, nzdim, ntdim)
  real(8) ::  ftzis(nxydim, nzdim, ntdim)

! ---- spatially varying isopycnal diffusion coefficient
  real(8), save ::   ahh3d(nxydim, nzdim),  ahi3d(nxydim, nzdim)
  real(8), save ::   ahg3d(nxydim, nzdim) 

!---- bolus velocity output (for CMIP6)
  real(8) :: ublsx(nxydim, nzdim), vblsy(nxydim, nzdim)
  real(8) :: ublsw(nxydim, nzdim), vblsw(nxydim, nzdim)

!---- for second order moment
!---- bug fix (save these variables)
  real(8), save ::  s0 (nxydim, nzdim, ntdim)=0.d0
  real(8), save ::  sm (nxydim, nzdim, ntdim)=0.d0
  real(8), save ::  sx (nxydim, nzdim, ntdim), sxx(nxydim, nzdim, ntdim)
  real(8), save ::  sy (nxydim, nzdim, ntdim), syy(nxydim, nzdim, ntdim)
  real(8), save ::  sz (nxydim, nzdim, ntdim), szz(nxydim, nzdim, ntdim)
  real(8), save ::  sxy(nxydim, nzdim, ntdim), sxz(nxydim, nzdim, ntdim)
  real(8), save ::  syz(nxydim, nzdim, ntdim)   

  real(8), save ::  f0 (nxydim, nzdim)=0.d0
  real(8), save ::  fm (nxydim, nzdim)=0.d0
  real(8), save ::  fx (nxydim, nzdim)=0.d0, fxx(nxydim, nzdim)=0.d0
  real(8), save ::  fy (nxydim, nzdim)=0.d0, fyy(nxydim, nzdim)=0.d0
  real(8), save ::  fz (nxydim, nzdim)=0.d0, fzz(nxydim, nzdim)=0.d0
  real(8), save ::  fxy(nxydim, nzdim)=0.d0, fxz(nxydim, nzdim)=0.d0
  real(8), save ::  fyz(nxydim, nzdim)=0.d0

  real(8), save ::  vlmx(nxydim, nzdim)=0.d0, vlmy(nxydim, nzdim)=0.d0
  real(8), save ::  vlmz(nxydim)=0.d0
  real(8), save ::  r(nxyzdm)
  real(8), save ::  alf(nxydim, nzdim)=0.d0, uv(nxydim, nzdim)

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

  real(8), save ::    ahb = 0.0d0
  real(8), save ::    ahh = 0.0d0,    ahi = 0.0d0,    ahg = 0.0d0

  namelist /nmdifb/ ahb
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

!---- latitudinally varying GM diffusivity
  integer, save ::  isvgm = -1
  real(8), save ::  ahgno = 1.d7, nlats =  40.d0, nlatn =  50.d0
  real(8), save ::  ahgso = 1.d7, slatn = -40.d0, slats = -50.d0
  real(8) :: pi, lat
  real(8) :: cort, omega

  namelist /nmsvgm/ isvgm
  namelist /nmdifn/ ahgno, nlats, nlatn
  namelist /nmdifs/ ahgso, slatn, slats

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
     !$acc update self(sx,sy,sz, sxx,syy,szz, sxy,sxz,syz)
     return
  end if

  if (ofirst) then
     ofirst = .false.
     call rewnml(ifpar, jfpar)
     read(ifpar, nmdifb, iostat=istat)
     call cstnml(jfpar, 'flxtrc', 'nmdifb', istat)
     write(jfpar, nmdifb)
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
     call rewnml(ifpar, jfpar)
     read(ifpar, nmsvgm, iostat=istat)
     call cstnml(jfpar, 'flxtrc', 'nmsvgm', istat)
     write(jfpar, nmsvgm)
     call rewnml(ifpar, jfpar)
     read(ifpar, nmdifn, iostat=istat)
     call cstnml(jfpar, 'flxtrc', 'nmdifn', istat)
     write(jfpar, nmdifn)
     call rewnml(ifpar, jfpar)
     read(ifpar, nmdifs, iostat=istat)
     call cstnml(jfpar, 'flxtrc', 'nmdifs', istat)
     write(jfpar, nmdifs)
     !$acc enter data create(ahi3d, ahg3d, ahh3d)
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
        
        if ( isvgm > 0 ) then
           write(jfpar, *) 'latitudinally varying GM diffusivity is used.'
           pi = atan( 1.d0 )*4.d0
           omega = 2.d0 * pi / 86400.d0
           do ij = ijstr, ijend
              cort = (cor(ij) + cor(ij+lw) + cor(ij+lsw) + cor(ij+ls)) * 0.25d0
              lat = asin( cort * 0.5d0 / omega ) * 180.d0 / pi
              if (lat .ge. slatn .and. lat .le. nlats) then
                 ahg3d(ij, :) = ahg
              else if (lat .ge. slats .and. lat .lt. slatn) then
                 ahg3d(ij, :) = (ahgso * (slatn - lat) + &
                      &            ahg * (lat - slats)) / (slatn - slats)
              else if (lat .gt. nlats .and. lat .le. nlatn) then
                 ahg3d(ij, :) = (ahgno * (lat - nlats) + &
                      &            ahg * (nlatn - lat)) / (nlatn - nlats)
              else if (lat .le. slats) then
                 ahg3d(ij, :) = ahgso
              else
                 ahg3d(ij, :) = ahgno
              endif
           end do
        end if

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
     !$acc update device(ahh3d, ahi3d, ahg3d)
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
     
     !$acc enter data create(ftx,fty,ftz, ftxd,ftyd,ftzd)
     !$acc enter data create(wzc,    rzm)
     !$acc enter data create(fharmx, fharmy,   harm)
     !$acc enter data create(hzbot)
     !$acc enter data create(dh)

     !$acc enter data create(xdzdx,  ydzdy)
     !$acc enter data create(zdzdx,  zdzdy)
     !$acc enter data create(xdtdz,  ydtdz)
     !$acc enter data create(zdtdx,  zdtdy)

     !$acc enter data create(adt2)
     !$acc enter data create(adtd)
     !$acc enter data create(adtah)
     !$acc enter data create(adtgm)
     !$acc enter data create(adtis)
     !$acc enter data create(ftx2)
     !$acc enter data create(fty2)
     !$acc enter data create(ftz2)
     !$acc enter data create(ftxah)
     !$acc enter data create(ftyah)
     !$acc enter data create(ftxgm)
     !$acc enter data create(ftygm)
     !$acc enter data create(ftzgm)
     !$acc enter data create(ftxis)
     !$acc enter data create(ftyis)
     !$acc enter data create(ftzis)

     !$acc enter data create(ublsx, vblsy)
     !$acc enter data create(ublsw, vblsw)

     !$acc enter data copyin(s0)
     !$acc enter data copyin(sm)

     !$acc enter data copyin(sx,sy,sz, sxx,syy,szz, sxy,sxz,syz)

     !$acc enter data copyin(f0)
     !$acc enter data copyin(fm)
     !$acc enter data copyin(fx , fxx)
     !$acc enter data copyin(fy , fyy)
     !$acc enter data copyin(fz , fzz)
     !$acc enter data copyin(fxy, fxz)
     !$acc enter data copyin(fyz)
     
     !$acc enter data copyin(vlmx, vlmy)
     !$acc enter data copyin(vlmz)     

     !$acc enter data create(r)
     !$acc enter data copyin(alf)
     !$acc enter data create(uv)

     !$acc enter data create(psigmx, psigmy)
     !$acc enter data create(xpsiy,  ypsix)
     !$acc enter data create(zpsix,  zpsiy)
     !$acc enter data create( igsy,   igsx)  
  end if

  tsiv = 1.d0 / ts

  call dnsgrd( &
     &  xdzdx,  ydzdy,  zdzdx,  zdzdy, &
     &  xdtdz,  ydtdz,  zdtdx,  zdtdy, &
     &  xpsiy,  ypsix, &
     &  zpsix,  zpsiy, &
     &     ty,     tx,     hz )

!$acc kernels default(present)
!$omp parallel
!$omp do
  do n = 1, ntdim
     do k = 1, nzdim
        do ij = 1, nxydim
           adt (ij, k, n) = 0.d0
           ftx (ij, k, n) = 0.d0
           fty (ij, k, n) = 0.d0
           ftz (ij, k, n) = 0.d0
           adt2 (ij, k, n) = 0.d0
           ftx2 (ij, k, n) = 0.d0
           fty2 (ij, k, n) = 0.d0
           ftz2 (ij, k, n) = 0.d0
           adtd (ij, k, n) = 0.d0
           ftxd (ij, k, n) = 0.d0
           ftyd (ij, k, n) = 0.d0
           ftzd (ij, k, n) = 0.d0
           ftxgm(ij, k, n) = 0.d0
           ftygm(ij, k, n) = 0.d0
           ftzgm(ij, k, n) = 0.d0
           ftxis(ij, k, n) = 0.d0
           ftyis(ij, k, n) = 0.d0
           ftzis(ij, k, n) = 0.d0
           ftxah(ij, k, n) = 0.d0
           ftyah(ij, k, n) = 0.d0
           adtgm(ij, k, n) = 0.d0
           adtis(ij, k, n) = 0.d0
           adtah(ij, k, n) = 0.d0
        end do
     end do
  end do
!$omp end do

!$omp do
  do k = 1, nzdim
     do ij = 1, nxydim
        diffz(ij, k) = 0.d0
     end do
  end do
!$omp end do
!$omp end parallel
!$acc end kernels

  do n = 1, ntdim
!$acc kernels default(present)
!$omp parallel private( &
!$omp k, ij, ijls, ijlw, &
!$omp fharmx, fharmy, harm &
!$omp )
!$omp do
     do k = kstr, kend
        do ij = ijtstr-nxdim, ijtend+nxdim+nxdim
           ijlw = ij + lw
           ijls = ij + ls
           fharmx(ij) = ahb * (hyu(ijlw) + hyu(ij+lsw)) &
             &              / (hxt(ij) + hxt(ijlw)) * &
             &          (tx(ij, k, n) - tx(ijlw, k, n)) * rx * &
             &          amskt(ij, k) * amskt(ijlw, k)
           fharmy(ij) = ahb * (hxu(ijls) + hxu(ij+lsw)) &
             &              / (hyt(ij) + hyt(ijls)) * &
             &          (tx(ij, k, n) - tx(ijls, k, n)) * rym(ij) * &
             &          amskt(ij, k) * amskt(ijls, k)
        end do
        do ij = ijtstr-nxdim, ijtend+nxdim
           harm(ij) = (  (fharmx(ij+le) - fharmx(ij)) * rx &
             &         + (fharmy(ij+ln) - fharmy(ij)) * ry(ij)) * &
             &        rxt(ij) * ryt(ij)
        end do
        do ij = ijtstr, ijtend+nxdim
           ijlw = ij + lw
           ijls = ij + ls
           ftx(ij, k, n) = - (harm(ij) - harm(ijlw)) * rx * &
             &             (hyu(ijlw) + hyu(ij+lsw)) &
!             &             / (hxt(ij) + hxt(ijlw)) * amftx(ij, kstr) &
             &             / (hxt(ij) + hxt(ijlw)) * amftx(ij, k)
           fty(ij, k, n) = - (harm(ij) - harm(ijls)) * rym(ij) * &
             &             (hxu(ijls) + hxu(ij+lsw)) &
!             &             / (hyt(ij) + hyt(ijls)) * amfty(ij, kstr) &
             &             / (hyt(ij) + hyt(ijls)) * amfty(ij, k)
        end do
     end do
!$omp end do
!$omp end parallel
!$acc end kernels
  end do

!$acc kernels default(present)
!$omp parallel
!$omp do
  do ij = 1, nxydim
     hzbot(ij) = hz(ij) + zbot
  end do
!$omp end do

! ---- vertical velocity on sigma coordinate
!$omp do
  do k = kstr, kstr+kz-1
     do ij = 1, nxydim
        wzc(ij, k) = w(ij, k) * hzbot(ij)
        rzm(ij, k) = 1.d0 / dsm(k) / hzbot(ij)
     end do
  end do
!$omp end do
!$omp do
  do k = kstr+kz, kend
     do ij = 1, nxydim
        wzc(ij, k) = w(ij, k)
        rzm(ij, k) = 1.d0 / dzm(ij, k)
     end do
  end do
!$omp end do
!$omp end parallel
!$acc end kernels

  call chekin(wzc, 'WZC', &
     &     'ocean vertical velocity on sigma coordinate', 'cm/s', &
     & nx, ny, nz, nxyzdm, 'OCLVMT')

! ======  GM  isopycnal and diapycnal diffusion  ======
! ---- z diffusion flux of GM
!$acc kernels default(present)
!$omp parallel private(k, n, ij, kuu, ku, kd, ijls, ijlw, ijlsw)
!$omp do
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
           ftzd(ij, k, n) = ftz(ij, k, n)
           ftzgm(ij, k, n) =  &
             &   - ahg3d(ij, k) * &
             &   ( zdzdx(ij, k) * zdtdx(ij, k, n) &
             &   + zdzdy(ij, k) * zdtdy(ij, k, n) ) &
             &   * amftz(ij, k)
           ftzis(ij, k, n) =  &
             &  (  ahi3d(ij, k) * (  zdzdx(ij, k) * zdzdx(ij, k) &
             &                     + zdzdy(ij, k) * zdzdy(ij, k) ) * &
             &       rzm(ij, k) * (tx(ij, ku, n) - tx(ij, k, n)) &
             &   - ahi3d(ij, k) * &
             &   ( zdzdx(ij, k) * zdtdx(ij, k, n) &
             &   + zdzdy(ij, k) * zdtdy(ij, k, n) ) ) &
             &   * amftz(ij, k)

        end do
     end do
  end do
!$omp end do

! ---- y diffusion flux of GM
!$omp do
  do n = 1, ntdim

     do k = kstr, kend
        do ij = ijtstr, ijtend+nxdim

           ijls = ij + ls

           fty(ij, k, n) = fty(ij, k, n) + &
             &     (  ( ahh3d(ij, k) + ahi3d(ij, k) ) * rym(ijls) &
             &      / ( hyt(ij) + hyt(ijls) ) * &
             &        ( tx(ij, k, n) - tx(ijls, k, n) ) * 2.d0 &
             &      - ( ( ahi3d(ij, k) - ahg3d(ij, k) ) &
             &        * ydzdy(ij, k) - ypsix(ij, k) ) * ydtdz(ij, k, n) ) &
             &      * ( hxu(ijls) + hxu(ij+lsw) ) * 0.5d0 * amfty(ij, k)
           ftyd(ij, k, n) = fty(ij, k, n)
           ftyah(ij, k, n) = &
             &     (  ahh3d(ij, k) * rym(ijls) &
             &      / ( hyt(ij) + hyt(ijls) ) * &
             &        ( tx(ij, k, n) - tx(ijls, k, n) ) * 2.d0 ) &
             &      * ( hxu(ijls) + hxu(ij+lsw) ) * 0.5d0 * amfty(ij, k)
           ftygm(ij, k, n) = &
             &       (  ahg3d(ij, k) &
             &        * ydzdy(ij, k) * ydtdz(ij, k, n) ) &
             &        * ( hxu(ijls) + hxu(ij+lsw) ) * 0.5d0 * amfty(ij, k)
           ftyis(ij, k, n) = &
             &     (  ahi3d(ij, k) * rym(ijls) &
             &      / ( hyt(ij) + hyt(ijls) ) * &
             &        ( tx(ij, k, n) - tx(ijls, k, n) ) * 2.d0 &
             &      - ahi3d(ij, k) &
             &      * ydzdy(ij, k) * ydtdz(ij, k, n) ) &
             &      * ( hxu(ijls) + hxu(ij+lsw) ) * 0.5d0 * amfty(ij, k)
        end do
     end do

  end do
!$omp end do

! ---- x diffusion flux of GM
!$omp do
  do n = 1, ntdim

     do k = kstr, kend
        do ij = ijtstr, ijtend+1

           ijlw = ij + lw
           
           ftx(ij, k, n) = ftx(ij, k, n) + &
             &     (  ( ahh3d(ij, k) + ahi3d(ij, k) ) * rx &
             &      / ( hxt(ij) + hxt(ijlw) ) &
             &      * ( tx(ij, k, n) - tx(ijlw, k, n) ) * 2.d0 &
             &      - ( ( ahi3d(ij, k) - ahg3d(ij, k) ) &
             &      * xdzdx(ij, k) + xpsiy(ij, k) ) * xdtdz(ij, k, n) ) &
             &      * ( hyu(ijlw) + hyu(ij+lsw) ) * 0.5d0 * amftx(ij, k)
           ftxd(ij, k, n) = ftx(ij, k, n)
           ftxah(ij, k, n) = &
             &     (  ahh3d(ij, k) * rx &
             &      / ( hxt(ij) + hxt(ijlw) ) &
             &      * ( tx(ij, k, n) - tx(ijlw, k, n) ) * 2.d0 ) &
             &      * ( hyu(ijlw) + hyu(ij+lsw) ) * 0.5d0 * amftx(ij, k)
           ftxgm(ij, k, n) = &
             &     (  ahg3d(ij, k) &
             &      * xdzdx(ij, k) * xdtdz(ij, k, n) ) &
             &      * ( hyu(ijlw) + hyu(ij+lsw) ) * 0.5d0 * amftx(ij, k)
           ftxis(ij, k, n) = &
             &     (  ahi3d(ij, k) * rx &
             &      / ( hxt(ij) + hxt(ijlw) ) &
             &      * ( tx(ij, k, n) - tx(ijlw, k, n) ) * 2.d0 &
             &      - ahi3d(ij, k) &
             &      * xdzdx(ij, k) * xdtdz(ij, k, n) ) &
             &      * ( hyu(ijlw) + hyu(ij+lsw) ) * 0.5d0 * amftx(ij, k)

        end do

! ---- divergence of diffusion fluxes

!        do ij = ijtstr, ijtend
!
!           adt(ij, k, n) = &
!             & (  (  (ftx(ij+le, k, n) - ftx(ij, k, n)) * rx &
!             &     + (fty(ij+ln, k, n) - fty(ij, k, n)) * ry(ij)) * &
!             &    rxt(ij) * ryt(ij) &
!             &  + ftz(ij, k, n) - ftz(ij, k+1, n)) / dz(ij, k)
!
!        end do
     end do

  end do
!$omp end do

!---- bolus velocity (for CMIP6 output)
!$omp do
  do k = 1, nzdim
     do ij = 1, nxydim
        ublsx(ij, k) = 0.0d0
        ublsw(ij, k) = 0.0d0
        vblsy(ij, k) = 0.0d0
        vblsw(ij, k) = 0.0d0
     end do
  end do
!$omp end do

!$omp do
  do k = kstr+1, kend
     do ij = ijtstr, ijtend+nxdim
        ublsw(ij, k) = ( ahg3d(ij, k-1) * xdzdx(ij, k-1) &
             &         - ahg3d(ij, k  ) * xdzdx(ij, k  ) ) &
             &       / dzm(ij, k) &
             &       * (1.0d0- &
             &         (1.0d0-amftz(ij,k))*(1.0d0-amftz(ij+lw,k)))
        vblsw(ij, k) = ( ahg3d(ij, k-1) * ydzdy(ij, k-1) &
             &         - ahg3d(ij, k  ) * ydzdy(ij, k  ) ) &
             &       / dzm(ij, k) &
             &       * (1.0d0- &
             &         (1.0d0-amftz(ij,k))*(1.0d0-amftz(ij+ls,k)))
     end do
  end do
!$omp end do
!$omp do
  do k = kstr, kend
     do ij = ijtstr, ijtend+nxdim
        ublsx(ij, k) = 0.5d0 * (ublsw(ij, k) + ublsw(ij, k+1))
        vblsy(ij, k) = 0.5d0 * (vblsw(ij, k) + vblsw(ij, k+1))
     end do
  end do
!$omp end do

!---- diffusion in BBL
#ifdef OPT_BBL

!$omp do
  do n = 1, ntdim
     do ij = ijtstr, ijtend
        k = nbot(ij)
        ftz(ij, kend, n) = ftz(ij, k, n)
     end do
  end do
!$omp end do

!$omp do
  do n = 1, ntdim

     do ij = ijtstr, ijtend+nxdim
        ijls = ij + ls
        fty(ij, kend, n) = &
          &   ahhbbl * (tx(ij, kend, n) - tx(ijls, kend, n)) * &
          &   rym(ij) * &
          &   (hxu(ijls) + hxu(ij+lsw)) / (hyt(ij) + hyt(ijls)) * &
          &   amfty(ij, kend)
        ftyah(ij, kend, n) = fty(ij, kend, n)
     end do
     
     do ij = ijtstr, ijtend+1
        ijlw = ij + lw
        ftx(ij, kend, n) = &
          &   ahhbbl * (tx(ij, kend, n) - tx(ijlw, kend, n)) * rx * &
          &   (hyu(ijlw) + hyu(ij+lsw)) / (hxt(ij) + hxt(ijlw)) * &
          &   amftx(ij, kend) 
        ftxah(ij, kend, n) = ftx(ij, kend, n)
     end do

  end do
!$omp end do

#endif

!$omp do
  do k = kstr, kend
     do ij = ijtstr, ijtend+nxdim
        igsx(ij, k) = ( ahi3d(ij, k) - ahg3d(ij, k) ) * xdzdx(ij, k)
        igsy(ij, k) = ( ahi3d(ij, k) - ahg3d(ij, k) ) * ydzdy(ij, k)
     end do
  end do
!$omp end do
!  call chekin(igsx, 'IGSX', nx, ny, nz, nxyzdm, 'OCN')
!  call chekin(igsy, 'IGSY', nx, ny, nz, nxyzdm, 'OCN')

!$omp do
  do k = 1, nzdim
     do ij = 1, nxydim
        psigmx(ij, k) = 0.0d0
        psigmy(ij, k) = 0.0d0
     end do
  end do
!$omp end do

!$omp do
  do k = kstr, kend
     do ij = ijtstr, ijtend+nxdim
        psigmx(ij, k) = ahg3d(ij, k) * ydzdy(ij, k)
        psigmy(ij, k) = - ahg3d(ij, k) * xdzdx(ij, k)
     end do
  end do
!$omp end do
!$omp end parallel
!$acc end kernels
  
!  call chekin(psigmx, 'PSIGMX', nx, ny, nz, nxyzdm, 'OCN')
!  call chekin(psigmy, 'PSIGMY', nx, ny, nz, nxyzdm, 'OCN')


! ---- SOM 

! ---- mass contained in a tracer grid

!$acc kernels default(present)
!$omp do
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
!$omp end do

! ---- in X-direction

!$omp do
  do k = kstr, kend
     do ij = ijtstr-nxdim-1, ijtend+nxdim+1

        ijlw  = ij + lw
        ijlsw = ij + lsw
               
        uv(ij, k) = ( uy(ijlw , k) * vlmy(ijlw, k) &
          &         + uy(ijlsw, k) * vlmy(ijlsw, k) ) &
          &       * amskt(ij, k) * amskt(ijlw, k) * dy(ijlw)

     end do
  end do
!$omp end do
!$omp end parallel
!$acc end kernels
  
  do n = 1, ntdim

! ---- undershoot limiter (Method B) of Morales Maqueda and Holloway (2006)
!$acc kernels default(present)
!!*POPTION PARALLEL
!$omp parallel do private( &
!$omp ij, ijlw, ijle, k, s0m, s1m, s0p, sxp, &
!$omp alfq, alf1, alf1q, tmp &
!$omp )
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

!    ---- overshoot limiter (Method B) of Morales Maqueda and Holloway (2006)
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

!---- bug fix 2
!     call shift1( sx (:,:,n), nxdim, nydim, nzdim )
!     call shift1( sxx(:,:,n), nxdim, nydim, nzdim )
!     call shift1( sxy(:,:,n), nxdim, nydim, nzdim )
!     call shift1( sxz(:,:,n), nxdim, nydim, nzdim )

!    ---- calculating ALF and MASS between box (i-1,j,k) <---> (i,j,k)
        do ij = ijtstr-nxdim-1, ijtend+nxdim+1

           ijlw  = ij + lw
               
           if ( uv(ij, k) .gt. 0.d0 ) then
              fm(ijlw, k) = uv(ij, k) * ts 
              alf(ij, k)  = fm(ijlw, k) / sm(ijlw, k, n)
           else
              fm(ijlw, k) = - uv(ij, k) * ts
              alf(ij, k)  = fm(ijlw, k) / sm(ij  , k, n)
           end if

!    ---- calculating flux and moments between box (i-1,j,k) <---> (i,j,k)

           alfq  = alf(ij, k) * alf(ij, k)
           alf1  = 1.d0 - alf(ij, k)
           alf1q = alf1 * alf1

           if ( uv(ij, k) .gt. 0.d0 ) then
!          ---- flux from (i-1) to (i),  when u > 0

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

              ftx2(ij, k, n) =-f0(ijlw, k) * tsiv * ry(ijlw)
              ftx (ij, k, n) = ftx(ij, k, n) &
                &            - f0(ijlw, k) * tsiv * ry(ijlw)

           else
!          ---- flux from (i) to (i-1),  when u < 0
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
             
              ftx2(ij, k, n) = f0(ijlw, k) * tsiv * ry(ijlw)
              ftx (ij, k, n) = ftx(ij, k, n) &
                &            + f0(ijlw, k) * tsiv * ry(ijlw)

           end if
          
        end do

!    ---- calculating flux and moments between box (i-1,j,k) <---> (i,j,k)
        do ij = ijtstr-nxdim-1, ijtend+nxdim+1

           ijlw  = ij + lw

           alfq  = alf(ij, k) * alf(ij, k)
           alf1  = 1.d0 - alf(ij, k)
           alf1q = alf1 * alf1
           if ( uv(ij, k) .le. 0.d0 ) then
!          ---- flux from (i) to (i-1),  when u < 0
             
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

        do ij = ijtstr-nxdim-1, ijtend+nxdim+1

           ijlw  = ij + lw

           alfq  = alf(ij, k) * alf(ij, k)
           alf1  = 1.d0 - alf(ij, k)
           alf1q = alf1 * alf1
           if ( uv(ij, k) .gt. 0.d0 ) then
!          ---- flux from (i-1) to (i),  when u > 0

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
           end if
        end do

!    ---- put the temporary moments (fi) into appropriate neighboring boxes
        do ij = ijtstr-nxdim-1, ijtend+nxdim+1
           ijlw  = ij + lw
           if ( uv(ij, k) .gt. 0.d0 ) then
              sm(ij  , k, n) = sm(ij, k, n) + fm(ijlw, k)
              alf(ij, k)     = fm(ijlw, k) / sm(ij, k, n)
           end if
        end do
        do ij = ijtstr-nxdim-1, ijtend+nxdim+1
           ijlw  = ij + lw
           if ( uv(ij, k) .le. 0.d0 ) then
              sm(ijlw, k, n) = sm(ijlw, k, n) + fm(ijlw, k)
              alf(ij, k)     = fm(ijlw, k) / sm(ijlw, k, n)
           end if
        end do

!       ---- flux from (i-1) to (i),  when u > 0
        do ij = ijtstr-nxdim-1, ijtend+nxdim+1
           ijlw  = ij + lw
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

           end if
        end do
!       ---- flux from (i) to (i-1),  when u < 0
        do ij = ijtstr-nxdim-1, ijtend+nxdim+1
           ijlw  = ij + lw
           if ( uv(ij, k) .le. 0.d0 ) then

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
!$omp end parallel do
!$acc end kernels
  end do

! ---- Y-direction
!$acc kernels default(present)
!$omp parallel do &
!$omp private( ij, k, ijls, ijlsw )
  do k = kstr, kend
     do ij = ijtstr-nxdim-1, ijtend+nxdim+1

        ijls  = ij + ls
        ijlsw = ij + lsw
        uv(ij, k) = ( vy(ijls , k) * vlmx(ijls , k) &
          &         + vy(ijlsw, k) * vlmx(ijlsw, k) ) &
          &       * amskt(ij, k) * amskt(ijls, k) 

     end do
  end do
!$omp end parallel do
!$acc end kernels
  do n = 1, ntdim

!    ---- undershoot limiter (Method B) of Morales Maqueda and Holloway (2006)
!$acc kernels default(present)
!!*POPTION PARALLEL
!$omp parallel private( &
!$omp ij, ijls, ijln, k, s0m, s1m, s0p, sxp, &
!$omp alfq, alf1, alf1q, tmp &
!$omp )
!$omp do
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

!    ---- overshoot limiter (Method B) of Morales Maqueda and Holloway (2006)
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
!$omp end do
!$omp end parallel
!$acc end kernels

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
     call shift1( s0(:,:,n), &
       &           nxdim,  nydim,  nzdim, &
       &            1.d0,      0,      0 )
#else
     call shift1( sy (:,:,n), nxdim, nydim, nzdim )
     call shift1( syy(:,:,n), nxdim, nydim, nzdim )
     call shift1( sxy(:,:,n), nxdim, nydim, nzdim )
     call shift1( syz(:,:,n), nxdim, nydim, nzdim )
     call shift1( s0 (:,:,n), nxdim, nydim, nzdim )
#endif

!    ---- calculating ALF  and MASS between box (i,j-1,k) <---> (i,j,k)
!$acc kernels default(present)
!!*POPTION PARALLEL
!$omp parallel private( &
!$omp ij, ijls, ijln, k, s0m, s1m, s0p, sxp, &
!$omp alfq, alf1, alf1q, tmp &
!$omp )
!$omp do
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

!    ---- calculating flux between box (i,j-1,k) <---> (i,j,k)
        do ij = ijtstr-nxdim-1, ijtend+nxdim+1

           ijls  = ij + ls

           alfq  = alf(ij, k) * alf(ij, k)
           alf1  = 1.d0 - alf(ij, k)
           alf1q = alf1 * alf1
           
           if ( uv(ij, k) .gt. 0.d0 ) then
!          ---- flux from (j-1) to (j),  when v > 0

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
             
              fty2(ij, k, n) =-f0(ijls, k) * tsiv * rx
              fty (ij, k, n) = fty(ij, k, n) &
                &            - f0(ijls, k) * tsiv * rx

           else
!          ---- flux from (j) to (j-1),  when v < 0

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
             
              fty2(ij, k, n) = f0(ijls, k) * tsiv * rx
              fty (ij, k, n) = fty(ij, k, n) &
                &            + f0(ijls, k) * tsiv * rx

           end if
          
        end do

!    ---- calculating moments
        do ij = ijtstr-nxdim-1, ijtend+nxdim+1

           ijls  = ij + ls

           if ( uv(ij, k) .le. 0.d0 ) then
!          ---- flux from (j) to (j-1),  when v < 0
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
        
        do ij = ijtstr-nxdim-1, ijtend+nxdim+1

           ijls  = ij + ls
           
           if ( uv(ij, k) .gt. 0.d0 ) then
!          ---- flux from (j-1) to (j),  when v > 0

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
             
           end if
             

        end do

!    ---- put the temporary moments (fi) into appropriate neighboring boxes
        do ij = ijtstr-nxdim-1, ijtend+nxdim+1
           ijls  = ij + ls
           if ( uv(ij, k) .gt. 0.d0 ) then
              sm(ij, k, n) = sm(ij, k, n) + fm(ijls, k)
              alf(ij, k)   = fm(ijls, k) / sm(ij, k, n)
           end if
        end do
        do ij = ijtstr-nxdim-1, ijtend+nxdim+1
           ijls  = ij + ls
           if ( uv(ij, k) .le. 0.d0 ) then
              sm(ijls, k, n) = sm(ijls, k, n) + fm(ijls, k)
              alf(ij, k)     = fm(ijls, k) / sm(ijls, k, n)
           end if
        end do

!       ---- flux from (j-1) to (j),  when v > 0
        do ij = ijtstr-nxdim-1, ijtend+nxdim+1

           ijls  = ij + ls

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
           end if
        end do
!       ---- flux from (j) to (j-1),  when v < 0
        do ij = ijtstr-nxdim-1, ijtend+nxdim+1

           ijls  = ij + ls

           if ( uv(ij, k) .le. 0.d0 ) then

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
!$omp end do
!$omp end parallel
!$acc end kernels
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
!$acc kernels default(present)
!$omp parallel do
  do k = kstr, kend
     do ij = ijtstr, ijtend

!        uv(ij, k) = - wzc(ij, k) * vlmz(ij)
        uv(ij, k) = - wzc(ij, k) * vlmz(ij) &
          &         * amskt(ij, k) * amskt(ij, k-1)
               
     end do
  end do
!$omp end parallel do
!$acc end kernels
  
#ifdef OPT_BBL
!$omp parallel do private(ij, k)
  do ij = ijtstr, ijtend
     k = nbot(ij)
     uv(ij, k) = - wzc(ij, k) * vlmz(ij) * amsktb(ij) &
          &      +  uv(ij, k) * (1.d0 - amsktb(ij))
  end do
!$omp end parallel do
#endif

  do n = 1, ntdim

!    ---- undershoot limiter (Method B) of Morales Maqueda and Holloway (2006)
!$acc kernels default(present)
!!*POPTION PARALLEL
!$omp parallel private( &
!$omp ij, k, ku, kd, s0m, s1m, s0p, sxp, &
!$omp alfq, alf1, alf1q, tmp &
!$omp )
!$omp do
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

!    ---- overshoot limiter (Method B) of Morales Maqueda and Holloway (2006)
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
!$omp end do

!    ---- calculating ALF
!$omp do
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

!    ---- calculating flux between box (i,j,k-1) <---> (i,j,k)

           alfq  = alf(ij, k) * alf(ij, k)
           alf1  = 1.d0 - alf(ij, k)
           alf1q = alf1 * alf1

!          ---- flux from (k-1) to (k),  when w > 0
           if ( uv(ij, k) .gt. 0.d0 ) then

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
                  
              ftz2(ij, k, n) = f0(ij, ku) * tsiv / vlmz(ij)
              ftz (ij, k, n) = ftz(ij, k, n) &
                &            + f0(ij, ku) * tsiv / vlmz(ij)

!          ---- flux from (k) to (k-1),  when w < 0
           else
              
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
             
              ftz2(ij, k, n) =-f0(ij, ku) * tsiv / vlmz(ij)
              ftz (ij, k, n) = ftz(ij, k, n) &
                &            - f0(ij, ku) * tsiv / vlmz(ij)

           end if
          
        end do
     end do
!$omp end do

!    ---- calculating flux between box (i,j,k-1) <---> (i,j,k)
!$omp do
     do k = kstr, kend

        ku = k - 1

        do ij = ijtstr, ijtend

           alfq  = alf(ij, k) * alf(ij, k)
           alf1  = 1.d0 - alf(ij, k)
           alf1q = alf1 * alf1
           if ( uv(ij, k) .gt. 0.d0 ) then
!          ---- flux from (k-1) to (k),  when w > 0

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
!          ---- flux from (k) to (k-1),  when w < 0
             
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
!$omp end do

!    ---- put the temporary moments (fi) into appropriate neighboring boxes
!$omp do
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

           alf1 = 1.d0 - alf(ij, k)
           if ( uv(ij, k) .gt. 0.d0 ) then
!          ---- flux from (k-1) to (k),  when w > 0

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

           else
!          ---- flux from (k) to (k-1),  when w < 0

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
!$omp end do
!$omp end parallel
!$acc end kernels
  end do

  do n = 1, ntdim
!$acc kernels default(present)
!$omp parallel
!$omp do
     do k = kstr, kend
        do ij = ijtstr, ijtend

! ---- tx
!           tx(ij, k, n) = s0(ij, k, n) / sm(ij, k, n)

           adt(ij, k, n) =                                       &
     &         (  (  (ftx(ij+le, k, n) - ftx(ij, k, n)) * rx         &
     &             + (fty(ij+ln, k, n) - fty(ij, k, n)) * ry(ij)) *  &
     &            rxt(ij) * ryt(ij)                                  &
     &          + ftz(ij, k, n) - ftz(ij, k+1, n)) / dz(ij, k)

           adt2(ij, k, n) = &
             & ( ( (ftx2(ij+le, k, n) - ftx2(ij, k,   n)) * rx &
             &   + (fty2(ij+ln, k, n) - fty2(ij, k,   n)) * ry(ij)) * &
             &      rxt(ij) * ryt(ij) &
             &    + ftz2(ij,    k, n) - ftz2(ij, k+1, n)) / dz(ij, k)
           adtd(ij, k, n) = &
             & ( ( (ftxd(ij+le, k, n) - ftxd(ij, k,   n)) * rx &
             &   + (ftyd(ij+ln, k, n) - ftyd(ij, k,   n)) * ry(ij)) * &
             &      rxt(ij) * ryt(ij) &
             &    + ftzd(ij,    k, n) - ftzd(ij, k+1, n)) / dz(ij, k)
           adtah(ij, k, n) = &
             & ( ( (ftxah(ij+le, k, n) - ftxah(ij, k,   n)) * rx &
             &   + (ftyah(ij+ln, k, n) - ftyah(ij, k,   n)) * ry(ij)) * &
             &      rxt(ij) * ryt(ij) ) / dz(ij, k)
           adtgm(ij, k, n) = &
             & ( ( (ftxgm(ij+le, k, n) - ftxgm(ij, k,   n)) * rx &
             &   + (ftygm(ij+ln, k, n) - ftygm(ij, k,   n)) * ry(ij)) * &
             &      rxt(ij) * ryt(ij) &
             &    + ftzgm(ij,    k, n) - ftzgm(ij, k+1, n)) / dz(ij, k)
           adtis(ij, k, n) = &
             & ( ( (ftxis(ij+le, k, n) - ftxis(ij, k,   n)) * rx &
             &   + (ftyis(ij+ln, k, n) - ftyis(ij, k,   n)) * ry(ij)) * &
             &      rxt(ij) * ryt(ij) &
             &    + ftzis(ij,    k, n) - ftzis(ij, k+1, n)) / dz(ij, k)
        end do
     end do
!$omp end do
!$omp end parallel
!$acc end kernels     
  end do

#ifdef OPT_BBL
  do n = 1, ntdim
     do ij = ijtstr, ijtend
        k = nbot(ij)
        ftz(ij, kend, n) = ftz(ij, k, n)
        ftz2 (ij, kend, n) = ftz2 (ij, k, n)
        ftzd (ij, kend, n) = ftzd (ij, k, n)
        ftzgm(ij, kend, n) = ftzgm(ij, k, n)
        ftzis(ij, kend, n) = ftzis(ij, k, n)
     end do
  end do

  do n = 1, ntdim
     do ij = ijtstr, ijtend
        adt(ij, kend, n) = &
             & ( ( (ftx(ij+le, kend, n) - ftx(ij, kend, n)) * rx &
             &   + (fty(ij+ln, kend, n) - fty(ij, kend, n)) * ry(ij)) * &
             &  rxt(ij) * ryt(ij) &
             & + ftz(ij, kend, n) ) / dz(ij, kend)
        adt2(ij, kend, n) = &
             & ( ( (ftx2(ij+le, kend, n) - ftx2(ij, kend, n)) * rx &
             &   + (fty2(ij+ln, kend, n) - fty2(ij, kend, n)) * ry(ij)) * &
             &  rxt(ij) * ryt(ij) &
             & + ftz2(ij, kend, n) ) / dz(ij, kend)
        adtd(ij, kend, n) = &
             & ( ( (ftxd(ij+le, kend, n) - ftxd(ij, kend, n)) * rx &
             &   + (ftyd(ij+ln, kend, n) - ftyd(ij, kend, n)) * ry(ij)) * &
             &  rxt(ij) * ryt(ij) &
             & + ftzd(ij, kend, n) ) / dz(ij, kend)
        adtgm(ij, kend, n) = &
             & ( ( (ftxgm(ij+le, kend, n) - ftxgm(ij, kend, n)) * rx &
             &   + (ftygm(ij+ln, kend, n) - ftygm(ij, kend, n)) * ry(ij)) * &
             &  rxt(ij) * ryt(ij) &
             & + ftzgm(ij, kend, n) ) / dz(ij, kend)
        adtis(ij, kend, n) = &
             & ( ( (ftxis(ij+le, kend, n) - ftxis(ij, kend, n)) * rx &
             &   + (ftyis(ij+ln, kend, n) - ftyis(ij, kend, n)) * ry(ij)) * &
             &  rxt(ij) * ryt(ij) &
             & + ftzis(ij, kend, n) ) / dz(ij, kend)
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

!---- for CMIP6 output
  call chekin(  ublsx, 'UBOLUS', &
  &            'G-M bolus velocity, x-dir.', 'cm/s', &
  &                nx,     ny,     nz, nxyzdm, 'OCLVTT')
  call chekin(  vblsy, 'VBOLUS', &
  &            'G-M bolus velocity, y-dir.', 'cm/s', &
  &                nx,     ny,     nz, nxyzdm, 'OCLVTT')

  call chekin(ftxgm(1, 1, 1), 'FTXGM', &
  &            'GM zonal heat flux', 'degC cm3/rad/s', &
  &            nx, ny, nz, nxyzdm, 'OCLVTX')
  call chekin(ftygm(1, 1, 1), 'FTYGM', &
  &            'GM meridional heat flux', 'degC cm3/rad/s', &
  &            nx, ny, nz, nxyzdm, 'OCLVTY')
  call chekin(ftxis(1, 1, 1), 'FTXIS', &
  &            'isopycnal zonal heat flux', 'degC cm3/rad/s', &
  &            nx, ny, nz, nxyzdm, 'OCLVTX')
  call chekin(ftyis(1, 1, 1), 'FTYIS', &
  &            'isopycnal meridional heat flux','degC cm3/rad/s', &
  &            nx, ny, nz, nxyzdm, 'OCLVTY')
  call chekin(ftxgm(1, 1, 2), 'FSXGM', &
  &            'GM zonal salt flux', 'psu cm3/rad/s', &
  &            nx, ny, nz, nxyzdm, 'OCLVTX')
  call chekin(ftygm(1, 1, 2), 'FSYGM', &
  &            'GM meridional salt flux', 'psu cm3/rad/s', &
  &            nx, ny, nz, nxyzdm, 'OCLVTY')
  call chekin(ftxis(1, 1, 2), 'FSXIS', &
  &            'isopycnal zonal salt flux', 'psu cm3/rad/s', &
  &            nx, ny, nz, nxyzdm, 'OCLVTX')
  call chekin(ftyis(1, 1, 2), 'FSYIS', &
  &            'isopycnal meridional salt flux', 'psu cm3/rad/s', &
  &            nx, ny, nz, nxyzdm, 'OCLVTY')
  call chekin(ftx2(1, 1, 1), 'FTX2', &
  &            'ocean zonal heat flux', 'degC cm3/rad/s', &
  &            nx, ny, nz, nxyzdm, 'OCLVTX')
  call chekin(fty2(1, 1, 1), 'FTY2', &
  &            'ocean meridional heat flux', 'degC cm3/rad/s', &
  &            nx, ny, nz, nxyzdm, 'OCLVTY')
  call chekin(ftz2(1, 1, 1), 'FTZ2', &
  &            'ocean vertical heat flux', 'degC cm3/rad/s', &
  &            nx, ny, nz, nxyzdm, 'OCLVMT')
  call chekin(ftx2(1, 1, 2), 'FSX2', &
  &            'ocean zonal salt flux', 'psu cm3/rad/s', &
  &            nx, ny, nz, nxyzdm, 'OCLVTX')
  call chekin(fty2(1, 1, 2), 'FSY2', &
  &            'ocean meridional salt flux', 'psu cm3/rad/s', &
  &            nx, ny, nz, nxyzdm, 'OCLVTY')
  call chekin(ftz2(1, 1, 2), 'FSZ2', &
  &            'ocean vertical salt flux', 'psu cm3/rad/s', &
  &            nx, ny, nz, nxyzdm, 'OCLVMT')

  call chekin(adtah(1, 1, 1), 'DTDTAHH', &
  &            'tendency of temp due to AHH', 'K/sec', &
  &            nx, ny, nz, nxyzdm, 'OCLVTT')
  call chekin(adtah(1, 1, 2), 'DSDTAHH', &
  &            'tendency of salt due to AHH', 'psu/sec', &
  &            nx, ny, nz, nxyzdm, 'OCLVTT')
  
  call chekin(adtgm(1, 1, 1), 'DTDTAHG', &
  &            'tendency of temp due to GM', 'K/sec', &
  &            nx, ny, nz, nxyzdm, 'OCLVTT')
  call chekin(adtgm(1, 1, 2), 'DSDTAHG', &
  &            'tendency of salt due to GM', 'psu/sec', &
  &            nx, ny, nz, nxyzdm, 'OCLVTT')
  call chekin(adtis(1, 1, 1), 'DTDTAHI', &
  &           'tendency of temp due to isopycnal diff', 'K/sec', &
  &            nx, ny, nz, nxyzdm, 'OCLVTT')
  call chekin(adtis(1, 1, 2), 'DSDTAHI', &
  &          'tendency of salt due to isopycnal diff', 'psu/sec', &
  &            nx, ny, nz, nxyzdm, 'OCLVTT')

  call chekin(adtd(1, 1, 1), 'DTDTD', &
  &            'tendency of temp by diffusion', 'K/sec', &
  &            nx, ny, nz, nxyzdm, 'OCLVTT')
  call chekin(adtd(1, 1, 2), 'DSDTD', &
  &            'tendency of salt by diffusion', 'psu/sec', &
  &            nx, ny, nz, nxyzdm, 'OCLVTT')

  !$acc kernels default(present)
  do ij = ijtstr, ijtend
     dh(ij) = hx(ij) - hz(ij)
  end do
  do n = 1, ntdim
     do k = kstr, kstr+kz-1
        do ij = ijtstr, ijtend
           adt2(ij, k, n) = adt2(ij, k, n) &
                &           - dh(ij) / zbot * tsiv * tx(ij, k, n)
        end do
     end do
  end do
  !$acc end kernels
  call chekin(adt2(1, 1, 1), 'DTDTV', &
       &            'tendency of temp by advection', 'K/sec', &
       &            nx, ny, nz, nxyzdm, 'OCLVTT')
  call chekin(adt2(1, 1, 2), 'DSDTV', &
       &            'tendency of salt by advection', 'psu/sec', &
       &            nx, ny, nz, nxyzdm, 'OCLVTT')

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
  logical, save :: ofirst = .true., ofirst2 = .true.

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
  logical, save ::  omlep  = .false.

  namelist /nmslpm/ slpmax
  namelist /nmmlep/ cm, ce, fminlt, mz, ofltdm, nfltdm, &
    &               lfmin, taumle, vscl, ofltps, nfltps, mzmin, &
    &               drsig, ocoamp, ofltrm, nfltrm, omlep

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
     !$acc enter data copyin(c0,c1,c2,c3,c4,c5,c6)
     !$acc enter data copyin(d0,d1,d2,d3,d4,d5,d6,d7,d8,d9)

     !$acc enter data create(cxpsy, cypsx)
     !$acc enter data create(czpsx, czpsy)

     !$acc enter data create(r)
     !$acc enter data create(   hmld, hmld1)
     !$acc enter data create( rmavez, rmav1)
     !$acc enter data create(  dzsig, dzmsig)
     !$acc enter data create(     zt,    ztm)
     !$acc enter data create( rsigth)
     !$acc enter data create(    nbv,     lf)
     !$acc enter data create( xpsiy1, ypsix1)
     !$acc enter data create( zpsix1, zpsiy1)
     !$acc enter data create(  zmld0,  zhmld)
     !$acc enter data create(  hmldx,  hmldy)
     !$acc enter data create( xpsiyz, ypsixz)
     !$acc enter data create(  kmld)

     !$acc enter data create( rmavdx, rmavdy)
     !$acc enter data create(   muzx,   muzy)

     !$acc enter data create(  xpsiy,  ypsix)
     !$acc enter data create(  zpsix,  zpsiy)
     !$acc enter data create(     hz)

     !$acc enter data create(   dtdx,   dtdy)
     !$acc enter data create(  dtfdz)       
  end if

!$acc kernels default(present)
!$omp parallel private(n, k, ij, &
!$omp tl, sl, p1, p2, rl, rlw, rls, rlu, dzdx, dzdy)
!$omp do
  do k = 1, nzdim
  do n = 1, ntdim
     do ij = 1, nxydim
        dtdx(ij, k, n) = 0.d0
        dtdy(ij, k, n) = 0.d0
        dtfdz(ij, k, n) = 0.d0
     end do
  end do
  end do

!$omp do
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
!$omp end do nowait

!$omp do
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

!$omp do
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

!$omp do
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

!$omp do
  do k = 1, nzdim
  do n = 1, ntdim
     do ij = 1, nxydim
        dtdx(ij, k, n) = 0.d0
        dtdy(ij, k, n) = 0.d0
        dtfdz(ij, k, n) = 0.d0
     end do
  end do
  end do

!$omp do
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
!$omp end do nowait

!$omp do
  do n = 1, ntdim
     do k = kstr+1, kend
        do ij = ijtstr-nxdim, ijtend+nxdim
           dtfdz(ij, k, n) = (tx(ij, k-1, n) - tx(ij, k, n)) &
              &              / dzm(ij, k) * amftz(ij, k)
        end do
     end do
  end do

!$omp do
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
!$omp end do nowait

!$omp do
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
!$omp end parallel

  xpsiy(:,:) = 0.d0
  ypsix(:,:) = 0.d0
  zpsix(:,:) = 0.d0
  zpsiy(:,:) = 0.d0
!$acc end kernels
  if ( .not. omlep ) return

  if (ofirst2) then
     ofirst2 = .false.
!    === Mixed layer eddy parameterization
!        from Fox-Kemper and Ferrari(2008) ===
!    *** Coarse-resolution modification is applied.
!        The last row in CXPSY/CYPSX corresponds to ds,
!        and is cancelled out by the row just before it. *** 
     pi = atan( 1.d0 )*4.d0
     omega = 2.d0 * pi / 86400.d0
     cormin = 2.d0 * omega * sin( pi*abs(fminlt)/180.d0 )
     kzmin = mzmin + kstr - 1

     !$acc kernels default(present)
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
     !$acc end kernels
  end if

!$acc kernels default(present)
!$omp do
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
!$acc end kernels
  
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
  !$acc kernels default(present)
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
  !$acc end kernels

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

  !$acc kernels default(present)
  do ij = 1, nxydim
     if (kmld(ij) .lt. kzmin) then
        hmld(ij) = 0.0d0
     end if
     zmld0(ij) = kmld(ij) - kstr + 1
  end do
  !$acc end kernels
  
! for diagnosis
!  call chekin(   hmld,   'HMLD', nx, ny,  1, nxydim, 'SFC')
!  call chekin(  zmld0,   'KMLD', nx, ny,  1, nxydim, 'SFC')
!  call chekin( rmavez, 'RMAVEZ', nx, ny,  1, nxydim, 'SFC')
!  call chekin(      r,      'R', nx, ny, nz, nxyzdm, 'OCN')
!  call chekin(     lf,     'LF', nx, ny,  1, nxydim, 'SFC')

  !$acc kernels default(present)
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
  !$acc end kernels

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

  !$acc kernels default(present)
  do k = kstr+1, kend
     do ij = ijtstr, ijtend
        zpsix(ij, k) = 0.5d0 * &
          &   ( ypsixz(ij, k) + ypsixz(ij+ln, k) ) * amftz(ij, k)
        zpsiy(ij, k) = 0.5d0 * &
          &   ( xpsiyz(ij, k) + xpsiyz(ij+le, k) ) * amftz(ij, k)
     end do
  end do
  !$acc end kernels
  
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
     &            nx,     ny,     nz, nxyzdm, 'OCLVTX')
  call chekin(   fty,  'FTY', &
     &            'ocean meridional heat flux', 'degC cm^3/rad/s', &
     &            nx,     ny,     nz, nxyzdm, 'OCLVTY')
  call chekin(   ftz,  'FTZ', &
     &            'ocean vertical heat flux', 'degC cm/s', &
     &            nx,     ny,     nz, nxyzdm, 'OCLVMT')
  call chekin(   ftx(1, 1, 2),  'FSX', &
     &            'ocean zonal salt flux', 'psu cm^3/rad/s', &
     &            nx,     ny,     nz, nxyzdm, 'OCLVTX')
  call chekin(   fty(1, 1, 2),  'FSY', &
     &            'ocean meridional salt flux', 'psu cm^3/rad/s', &
     &            nx,     ny,     nz, nxyzdm, 'OCLVTY')
  call chekin(   ftz(1, 1, 2),  'FSZ', &
     &            'ocean vertical salt flux', 'psu cm/s', &
     &            nx,     ny,     nz, nxyzdm, 'OCLVMT')

  call chekin(   ftxd,  'FTXD', &
     &            'ocean zonal diffusive heat flux', 'degC cm^3/rad/s', &
     &            nx,     ny,     nz, nxyzdm, 'OCLVTX')
  call chekin(   ftyd,  'FTYD', &
     &            'ocean meridional diffusive heat flux', 'degC cm^3/rad/s', &
     &            nx,     ny,     nz, nxyzdm, 'OCLVTY')
  call chekin(   ftzd,  'FTZD', &
     &            'ocean vertical diffusive heat flux', 'degC cm/s', &
     &            nx,     ny,     nz, nxyzdm, 'OCLVMT')
  call chekin(   ftxd(1, 1, 2),  'FSXD', &
     &            'ocean zonal diffusive salt flux', 'psu cm^3/rad/s', &
     &            nx,     ny,     nz, nxyzdm, 'OCLVTX')
  call chekin(   ftyd(1, 1, 2),  'FSYD', &
     &            'ocean meridional diffusive salt flux', 'psu cm^3/rad/s', &
     &            nx,     ny,     nz, nxyzdm, 'OCLVTY')
  call chekin(   ftzd(1, 1, 2),  'FSZD', &
     &            'ocean vertical diffusive salt flux', 'psu cm/s', &
     &            nx,     ny,     nz, nxyzdm, 'OCLVMT')

  return
end subroutine chkftx

end module tflxt
