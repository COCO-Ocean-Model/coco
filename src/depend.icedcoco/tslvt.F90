module tslvt

! --- information -----------------------------------------------------
!
!  Time integration of the tracer equations.
!
!  HISTORY
!     '02.10.10  H.Hasumi: from MIROC3.1-OMIP
!     '07.04.23  H.Hasumi
!     '07.09.25  H.Hasumi: multi-category sea ice
!     '08.06.11  H.Hasumi: initial/final processing
!     '08.07.10  H.Hasumi: initial/final processing
!     '12.01.30  Y.Komuro: (change surface water flux by T. Suzuki)
!     '12.06.15  H.Tatebe: for COCO5.0 in F90
!     '13.04.22  T.Suzuki: bug fix (change surface water)
! ---------------------------------------------------------------------

  use zocdim, only :                                   &
       &   nxydim ,  nzdim,                            &
       &   ijtstr,  ijtend,                            &
       &     kstr,    kend,     nz
  use zocmsk,  only :                                  &
       &     nbot
  use zocfil, only :                                   &
       &      ncf

  implicit none
  private

  real(8),    save  ::  rgamma(nzdim)
  real(8),    save  ::  swconv(nxydim, nzdim)

  real(8),    save  ::      aa(nxydim, nzdim)
  real(8),    save  ::      ab(nxydim, nzdim)
  real(8),    save  ::      ac(nxydim, nzdim)

  real(8),    save  ::      dh(nxydim)  
  real(8),    save  ::   hzbot(nxydim),   hxbot(nxydim)
  real(8),    save  ::       gamma(nz),  swcnv1(nzdim)

  real(8),        save  ::   garea(nxydim),  sdmp2d(nxydim)
  real(8),        save  ::    sdmp,  dsmax
  logical,        save  ::  osrstr, osrsti, osrnml, ofsdmp
  character(len=ncf)    ::  cfsdmp

  real(8)           ::     dzb,    dzt
  real(8)           ::   radup,  raddn,  depth
  real(8)           ::  tswcnv
  real(8)           ::     rrr,  zeta1,  zeta2
  integer(4)        ::   ifpar,  jfpar,  istat

!---- for geothermal heating
  real(8),        save  ::   gthm(nxydim)
  logical,        save  ::  ogthm = .false.
!----

  real(8),    save  ::    smin = 5.d0
  integer(4), save  ::  mixsss = 0

  public  ::  slvtrc     !   aprdc.F
  public  ::  svtset     !   aocea.F

  namelist /nmacct/ gamma
  namelist /nmswab/   rrr,  zeta1,  zeta2
  namelist /nmsrst/  sdmp, osrstr, osrsti, osrnml, dsmax, &
    &              ofsdmp, cfsdmp

  namelist /nmmixsss/smin, mixsss

!  data gamma / nz*1.d0 /
  data rrr, zeta1, zeta2 / 5.8d-1, 3.5d+1, 2.3d+3 /
  data sdmp, dsmax / 0.d0, 999.d0 /
  data osrstr, osrsti, osrnml / .false., .false., .false. /
  data ofsdmp / .false. /
  data cfsdmp / 'not-specified' /

contains

! **********************************************************************

  subroutine svtset

    use zocdim,  only :                                &
         &  nxydim,  nxdim,  nydim,  nzdim,  ntdim,    &
         &    istr,   iend,   jstr,   jend,            &
         &    kstr,   kend,     kz,                    &
         &   ijstr,  ijend, ijtstr,  ijtend,           &
#ifndef OPT_IO_COCOMPI
         &     nxg,    nyg, nxgdim,  nygdim,           &
         &   igstr,  jgstr,                            &
#endif
         &   oinit,  ofinal
    use zocgrd,  only :                                &
         &      dx,     dy,    dz0,     hxt,    hyt
    use zocmsk,  only :                                &
         &   amskt,   nbot
    use zocphy,  only :                                &
         &     cpo,    rhoo

    use ifhea
    use bshft
    use qckag
    use ufile
    use zocfil,  only :                                &
         &     ncf
#ifdef OPT_IO_COCOMPI
    use mpiio
#else
    use bgs2d
    use zocnod,  only :  iroot,  myrank
#endif

    implicit none
#include "mpif.h"

    integer(4)        ::      ij,      k,      n,     l
    integer(4)        ::       i,      j
    logical,   save   ::  ofirst = .true.

!---- for geothermal heating
    character(len = ncf)  :: cfgthm
    character(len = 16)   :: chead(64) 
#ifdef OPT_IO_COCOMPI
    integer :: mpi_fh
    integer :: icread
    integer (kind = mpi_offset_kind) :: disp
#else
    integer ::  nfgthm, nfsdmp
    real(8) ::  buf2(nxg, nyg)
    real(8) ::  g2d(nxgdim, nygdim)
#endif
    data cfgthm /'not-specified'/
    namelist /nmgthm/ ogthm, cfgthm
!----

    if ( ofirst ) then

       ofirst = .false.
       gamma(1:nz) = 1.d0
       call rewnml( ifpar, jfpar )
       write(jfpar, *) '*** svtset ***'
       read(ifpar, nmacct, iostat = istat )
       call cstnml( jfpar, 'svtset', 'nmacct', istat )
       write(jfpar, nmacct)
       write(jfpar, *) ' level     gamma'
       do k = 1, nz
          write(jfpar, '(i6,f10.4)') k, gamma(k)
       end do
       call rewnml(ifpar, jfpar)
       read(ifpar, nmswab, iostat = istat )
       call cstnml(jfpar, 'svtset', 'nmswab', istat )
       write(jfpar, nmswab)
       call rewnml(ifpar, jfpar)
       read(ifpar, nmsrst, iostat = istat )
       call cstnml(jfpar, 'svtset', 'nmsrst', istat )
       write(jfpar, nmsrst)

       if ( osrstr ) then
!----- extended SSS restoring
         if ( .not.ofsdmp ) then
            do ij = 1, nxydim
               sdmp2d(ij) = 1.d0 / 8.64d4 / sdmp
            end do
         else
!----- reading SDMP coefficient from file
            write(jfpar, *) '  file name of gthm: ', cfsdmp
#ifdef OPT_IO_COCOMPI
            call mpi_filopn(mpi_fh, cfsdmp, 'READ')
            disp=0
            call mpi_read_chead(chead, mpi_fh, disp, icread)
            call mpi_read_2d(sdmp2d, mpi_fh  , disp)
            call mpi_filcls(mpi_fh)
#else   
            if ( myrank .eq. iroot ) then   
               call filopn( nfsdmp, cfsdmp, 'READ' )
               rewind( nfsdmp )
               read( nfsdmp ) chead
               read( nfsdmp ) buf2
               do j = 1, nyg
                  do i = 1, nxg
                     g2d(igstr+i-1, jgstr+j-1) = buf2(i, j)
                  end do
               end do
               call filcls( nfsdmp )               
            end if
            call scatter_2d( sdmp2d, g2d )
#endif
#ifdef OPT_TRIPOLE
            call shift1(sdmp2d, &
              &          nxdim,  nydim,      1, &
              &           1.d0,      0,      0 )
#else
            call shift1( &
              &         sdmp2d, &
              &          nxdim,  nydim,      1)
#endif
         end if

         do ij = 1, nxydim
            garea(ij) = 0.0d0
         end do
         if (osrnml) then
            write(jfpar, *) &
         &  '*** Normalization of SSS-restoring flux is applied. ***'
            do j=jstr, jend
               do i=istr, iend
                  ij = nxdim*(j-1) + i
                  garea(ij) = dx * dy(ij) * hxt(ij) * hyt(ij) &
                    &         * amskt(ij, kstr)
               end do
            end do
         endif         
       end if

       do k = kstr, kend
          rgamma(k) = 1.d0 / gamma(k-kstr+1)
       end do
  
       depth = 0.d0
       radup = 1.d0
       do k = kstr, kend
          depth = depth + dz0(k)
          raddn = rrr * exp(- depth / zeta1)                               &
       &        + (1.d0 - rrr) * exp(- depth / zeta2)
          swcnv1(k) = radup - raddn
          radup = raddn
       end do

       do k = 1, nzdim
          do ij = 1, nxydim
             swconv(ij, k) = 0.d0
          end do
       end do
       
       do ij = ijtstr, ijtend
          tswcnv = 0.d0
          do k = kstr, nbot(ij)-1
             swconv(ij, k) = swcnv1(k)
             tswcnv = tswcnv + swcnv1(k)
          end do
          if (nbot(ij) .ge. kstr) then
             swconv(ij, nbot(ij)) = 1.d0 - tswcnv
          end if
       end do
       
       call putswc( swconv(1, kstr) )
       call copswc( swconv(1, kstr) )
       
       do ij = 1, nxydim
          swconv(ij, kstr) = 0.d0
       end do
       
       do k = kstr+1, kend
          do ij = ijtstr, ijtend
             swconv(ij, k) = swconv(ij, k) / rhoo / cpo
          end do
       end do

!---- for geothermal heating
       call rewnml(ifpar, jfpar)
       read(ifpar, nmgthm, iostat=istat)
       call cstnml(jfpar, 'svtset', 'nmgthm', istat)
       write(jfpar, nmgthm)
       write(jfpar, *) '  file name of gthm: ', cfgthm

       if ( ogthm ) then
#ifdef OPT_IO_COCOMPI
          call mpi_filopn(mpi_fh, cfgthm, 'READ')
          disp=0
          call mpi_read_chead(chead, mpi_fh, disp, icread)
          call mpi_read_2d(gthm, mpi_fh  , disp)
          call mpi_filcls(mpi_fh)
#else

          if ( myrank .eq. iroot ) then

             call filopn( nfgthm, cfgthm, 'READ' )
             rewind( nfgthm )
             read( nfgthm ) chead
             read( nfgthm ) buf2
             do j = 1, nyg
                do i = 1, nxg
                   g2d(igstr+i-1, jgstr+j-1) = buf2(i, j)
                end do
             end do
             call filcls( nfgthm )
             
          end if
          call scatter_2d( gthm, g2d )
#endif
       end if
!----
       call rewnml(ifpar, jfpar)
       read(ifpar, nmmixsss, iostat = istat )
       call cstnml(jfpar, 'svtset', 'nmmixsss', istat )
       write(jfpar, nmmixsss)
       
    end if

  end subroutine svtset

! **********************************************************************

  subroutine slvtrc(                                   &
         &      tx,     hx,                            &
         &     adt,  diffz,                            &
         &      ft,  swabs,     fs,     hz,   ssfc,    &
         &      ax  )

    use zocdim,  only :                                &
         &      nx,     ny, nxydim,  nzdim,  ntdim,    &
         &    kstr,   kend,     kz,                    &
         &   ijstr,  ijend, ijtstr,  ijtend,           &
         &     nic,                                    &
         &  inodes, jnodes,                            &
         &   oinit,  ofinal
    use zocgrd,  only :                                &
         &      dz,    dz0,     ds,   zbot,            &
         &      ts
    use zocmsk,  only :                                &
#ifdef OPT_BBL
         &   amsktb,                                   &
#endif
         &   amskt,  nbot
    use zocnod,  only :                                &
         &   iroot,  ierr,  myrank
    use zocphy,  only :                                &
         &     cpo,    rhoo
    use utrdg
    use qckag
    use qckot

    implicit none

#include "mpif.h"

    real(8),   intent(inout)  ::     tx(nxydim, nzdim, ntdim)
    real(8),   intent(inout)  ::     hx(nxydim)
    real(8),   intent(inout)  ::    adt(nxydim, nzdim, ntdim)
    real(8),   intent(in)     ::  diffz(nxydim, nzdim)
    real(8),   intent(in)     ::     ft(nxydim, ntdim)
    real(8),   intent(in)     ::  swabs(nxydim)
    real(8),   intent(in)     ::     fs(nxydim)
    real(8),   intent(in)     ::     hz(nxydim)
    real(8),   intent(in)     ::   ssfc(nxydim)
    real(8),   intent(in)     ::     ax(nxydim, 0:nic)

    real(8)                   ::  fsrst(nxydim)
    real(8)                   :: vwteqg(inodes*jnodes)
    real(8)                   :: vareag(inodes*jnodes)
    real(8)                   :: vwteqt, vareat, tarea, dsss, fsnml

    integer(4)                ::     ij,    k,     n,      i,     j
    
    if ( oinit .or. ofinal ) then
       return
    end if

   
    do ij = ijtstr, ijtend
       hxbot(ij) = hx(ij) + zbot
    end do

    do k = 1, nzdim
       do ij = 1, nxydim
          aa(ij, k) = 0.d0
          ab(ij, k) = 0.d0
          ac(ij, k) = 0.d0
       end do
    end do
    
    do ij = 1, nxydim
       dh(ij) = hx(ij) - hz(ij)
    end do
    
    do k = kstr, kstr+kz-1
       do ij = ijtstr, ijtend
          aa(ij, k) = - ts * diffz(ij, k) / dz(ij, k)
          ac(ij, k) = - ts * diffz(ij, k+1) / dz(ij, k)
          ab(ij, k) = (hx(ij) + zbot) / zbot - aa(ij, k) - ac(ij, k)
       end do
    end do
    
    do n = 1, ntdim
       do k = kstr, kstr+kz-1
          do ij = ijtstr, ijtend
             adt(ij, k, n) = adt(ij, k, n)                            &
    &                      - dh(ij) / zbot / ts * tx(ij, k, n)
          end do
       end do
    end do
    do k = kstr+kz, kend
       do ij = ijtstr, ijtend
          aa(ij, k) = - ts * diffz(ij, k) / dz(ij, k)
          ac(ij, k) = - ts * diffz(ij, k+1) / dz(ij, k)
          ab(ij, k) = 1.d0 - aa(ij, k) - ac(ij, k)
       end do
    end do
    
    do ij = ijstr, ijend
!       adt(ij, kstr, 1) = adt(ij, kstr, 1)                            &
!    &                   + tx(ij, kstr, 1) * ft(ij, 2) / zbot          &
!    &                    * amskt(ij, kstr)
       adt(ij, kstr, 1) = adt(ij, kstr, 1)                            &
    &                   - tx(ij, kstr, 1) *ft(ij, 2)                  &
    &                     /dz(ij, kstr) * amskt(ij, kstr)         
       adt(ij, kstr, 2) = adt(ij, kstr, 2)                            &
    &                    - fs(ij)                                     &
    &                     /dz(ij, kstr) * amskt(ij, kstr)        
    end do

    call thomas( adt, ac, aa, ab )

    do n = 1, ntdim
       do k = kstr, kend
          do ij = ijtstr, ijtend
             adt(ij, k, n) = adt(ij, k, n) * rgamma(k) * amskt(ij, k)
          end do
       end do
    end do
#ifdef OPT_BBL
    do n = 1, ntdim
       do ij = ijtstr, ijtend
          k = max(nbot(ij), kstr)
          adt(ij, k, n) = adt(ij, k, n) *                                  &
    &                (  1.d0                                               &
    &                 + amsktb(ij) * (  gamma(k-kstr+1) * rgamma(kend)     &
    &                                 - 1.d0) )
       end do
    end do
#endif

    do n = 1, ntdim
       do k = kstr, kend
          do ij = ijtstr, ijtend
             tx(ij, k, n) = tx(ij, k, n) + ts * adt(ij, k, n)
          end do
       end do
    end do

!     '12.01.30: removed 
!      DO IJ = IJTSTR, IJTEND
!         HZBOT(IJ) = HX(IJ) + ZBOT
!         DH(IJ) = - TS * FT(IJ, 2) * AMSKT(IJ, KSTR)
!         HX(IJ) = HX(IJ) + DH(IJ)
!         HXBOT(IJ) = HX(IJ) + ZBOT
!         IF (DH(IJ) .GE. 0.D0) THEN
!            DZB = 0.D0
!            DZT = DH(IJ) * DS(KSTR+KZ-1)
!            DO K = KSTR+KZ-1, KSTR+1, -1
!               DO N = 1, NTDIM
!                  TX(IJ, K, N) = (  TX(IJ, K  , N) * DS(K) * HZBOT(IJ)
!     &                            - TX(IJ, K  , N) * DZB
!     &                            + TX(IJ, K-1, N) * DZT)
!     &                           / HXBOT(IJ) / DS(K)
!               END DO
!               DZB = DZT
!               DZT = DZB + DH(IJ) * DS(K-1)
!            END DO
!            TX(IJ, KSTR, 2) = (  TX(IJ, KSTR, 2) * DS(KSTR) * HZBOT(IJ)
!     &                         - TX(IJ, KSTR, 2) * DZB)
!     &                        / HXBOT(IJ) / DS(KSTR)
!            TX(IJ, KSTR, 1) = (  TX(IJ, KSTR, 1) * DS(KSTR) * HZBOT(IJ)
!     &                         - TX(IJ, KSTR, 1) * DZB
!     &                         + TX(IJ, KSTR, 1) * DZT
!     &                         + TS * FT(IJ, 1))
!     &                        / HXBOT(IJ) / DS(KSTR)
!            DO N = 3, NTDIM
!               TX(IJ, KSTR, N) = (  TX(IJ, KSTR, N)*DS(KSTR)*HZBOT(IJ)
!     &                            - TX(IJ, KSTR, N) * DZB
!     &                            + TX(IJ, KSTR, N) * DZT
!     &                            + TS * FT(IJ, N))
!     &                           / HXBOT(IJ) / DS(KSTR)
!            END DO
!         ELSE
!            DZT = - DH(IJ)
!            DZB = DZT + DH(IJ) * DS(KSTR)
!            TX(IJ, KSTR, 2) = (  TX(IJ, KSTR, 2) * DS(KSTR) * HZBOT(IJ)
!     &                         + TX(IJ, KSTR+1, 2) * DZB)
!     &                        / HXBOT(IJ) / DS(KSTR)
!            TX(IJ, KSTR, 1) = (  TX(IJ, KSTR, 1) * DS(KSTR) * HZBOT(IJ)
!     &                         + TX(IJ, KSTR+1, 1) * DZB
!     &                         - TX(IJ, KSTR, 1) * DZT
!     &                         + TS * FT(IJ, 1))
!     &                        / HXBOT(IJ) / DS(KSTR)
!            DO N = 3, NTDIM
!               TX(IJ, KSTR, N) = (  TX(IJ, KSTR, N)*DS(KSTR)*HZBOT(IJ)
!     &                            + TX(IJ, KSTR+1, N) * DZB
!     &                            - TX(IJ, KSTR  , N) * DZT
!     &                            + TS * FT(IJ, N))
!     &                           / HXBOT(IJ) / DS(KSTR)  
!            END DO            
!            DO K = KSTR+1, KSTR+KZ-2
!               DZT = DZB
!               DZB = DZT + DH(IJ) * DS(K)
!               DO N = 1, NTDIM
!                  TX(IJ, K, N) = (  TX(IJ, K, N) * DS(K) * HZBOT(IJ)
!     &                            + TX(IJ, K+1, N) * DZB
!     &                            - TX(IJ, K  , N) * DZT)
!     &                           / HXBOT(IJ) / DS(K)
!               END DO
!            END DO
!         END IF
!      END DO

    call cofpsf( &
      &             tx,     ft,     fs,  swabs)

    do ij = ijtstr, ijtend
       tx(ij, kstr, 1) = tx(ij, kstr, 1)                              &
    &                  + ts * ft(ij, 1) / hxbot(ij) / ds(kstr) 
    end do

    do k = kstr, kstr+kz-1
       do ij = ijtstr, ijtend
          tx(ij, k, 1) = tx(ij, k, 1)                                 &
    &                  + ts * swconv(ij, k) * swabs(ij)               &
    &                  / hxbot(ij) / ds(k)
       end do
    end do
    do k = kstr+kz, kend
       do ij = ijtstr, ijtend
          tx(ij, k, 1) = tx(ij, k, 1)                                &
    &                  + ts * swconv(ij, k) * swabs(ij) / dz(ij, k)
       end do
    end do

!    do ij = ijtstr, ijtend
!       tx(ij, kstr, 2) = tx(ij, kstr, 2)                             &
!    &                  - ts * fs(ij) / hxbot(ij) / ds(kstr)
!    end do

    do n = 3, ntdim
       do ij = ijtstr, ijtend
          tx(ij, kstr, n) = tx(ij, kstr, n)                          &
    &                     + ts * ft(ij, n) / hxbot(ij) / ds(kstr)
       end do
    end do

#ifdef OPT_SRST
    do ij = 1, nxydim
       fsrst(ij) = 0.0d0
    end do

    if (osrstr) then
!       call tmintp(ssfc, 10)
       vwteqt = 0.0d0
       vareat = 0.0d0
       if (osrsti) then
          do ij = ijtstr, ijtend
             dsss = (ssfc(ij) - tx(ij, kstr, 2)) * amskt(ij, kstr)
             fsrst(ij) = sdmp2d(ij) &
               &       * sign(min(abs(dsss), dsmax), dsss)
             vwteqt = vwteqt + fsrst(ij) * garea(ij)
             vareat = vareat + garea(ij)
          end do
       else
          do ij = ijtstr, ijtend
             if (ax(ij, 0) .eq. 1.d0) then
                dsss = (ssfc(ij) - tx(ij, kstr, 2)) * amskt(ij, kstr)
                fsrst(ij) = sdmp2d(ij) &
                  &       * sign(min(abs(dsss), dsmax), dsss)
                vwteqt = vwteqt + fsrst(ij) * garea(ij)
                vareat = vareat + garea(ij)
             end if
          end do
       end if

       if (osrnml) then
          fsnml = 0.0d0
          tarea = 0.0d0
          do i = 1, inodes*jnodes
             vwteqg(i) = 0.0d0
             vareag(i) = 0.0d0
          end do

          call mpi_gather( &
            &    vwteqt, 1, mpi_real8, vwteqg(1), 1, mpi_real8, &
            &    iroot, mpi_comm_world, ierr)
          call mpi_gather( &
            &    vareat, 1, mpi_real8, vareag(1), 1, mpi_real8, &
            &    iroot, mpi_comm_world, ierr)
          if (myrank .eq. iroot) then
             do i = 1, inodes*jnodes
                fsnml = fsnml + vwteqg(i)
                tarea = tarea + vareag(i)
             end do
          end if
          call mpi_bcast( &
            &    fsnml, 1, mpi_real8, &
            &    iroot, mpi_comm_world, ierr)
          call mpi_bcast( &
            &    tarea, 1, mpi_real8, &
            &    iroot, mpi_comm_world, ierr)
          fsnml = fsnml / tarea
          if (osrsti) then
             do ij = ijtstr, ijtend
                fsrst(ij) = ( fsrst(ij) - fsnml ) * amskt(ij, kstr)
             end do
          else
             do ij = ijtstr, ijtend
                if (ax(ij, 0) .eq. 1.d0) then
                   fsrst(ij) = ( fsrst(ij) - fsnml ) * amskt(ij, kstr)
                end if
             end do
          end if
       end if

       do ij = ijtstr, ijtend
          tx(ij, kstr, 2) = tx(ij, kstr, 2) + ts * fsrst(ij)
       end do

       do ij = ijtstr, ijtend
          fsrst(ij) = fsrst(ij) * hxbot(ij) * ds(kstr)
       end do

       call cofpsr( &
         &           fsrst)

       call chekin( fsrst, 'FSRST', &
         &          'SSS resotring flux', 'psu cm/s', &
         &          nx,     ny,      1, nxydim, 'OCSFCT') 
    end if
#endif

    if (ogthm) then
       do ij = ijtstr, ijtend
          k = nbot(ij)
          tx(ij, k, 1) = tx(ij, k, 1) &
               & + ts * gthm(ij) / dz(ij, k) / rhoo / cpo
       end do
    end if

!---- mixing salinity in sigma-layers to avoid extremely low SSS
    if ( mixsss > 0 ) then
       call tmixss( tx(1, 1, 2) )
    end if

  end subroutine slvtrc

  subroutine tmixss(                                   & !! mix sea surface
         &      tx  )

    use zocdim,  only :                                &
         &  nxydim,  nzdim,                            &
         &    kstr,     kz
    use zocgrd,  only :                                &
         &     dz0
    use zocmsk,  only :                                &
         &   amskt

    implicit none

    real(8),    intent(inout) :: tx(nxydim, nzdim)

    real(8)                   ::  smean(nxydim), ssum(nxydim)
    integer(4)                ::   kmix(nxydim)
    integer(4)                ::     ij,    k

    depth = dz0(kstr)
    do ij = 1, nxydim
       ssum(ij) = tx(ij, kstr) * dz0(kstr)
       smean(ij) = tx(ij, kstr)
       kmix(ij) = kstr
    end do

    do k = kstr+1, kstr+kz-1
       depth = depth + dz0(k)
       do ij = 1, nxydim
          if ( smean(ij) <= smin ) then
             ssum(ij) = ssum(ij) + dz0(k) * tx(ij, k)
             smean(ij) = ssum(ij) / depth
             kmix(ij) = k
          end if
       end do
    end do

    do k = kstr, kstr+kz-1
       do ij = 1, nxydim
          if ( k <= kmix(ij) ) then
             tx(ij, k) = smean(ij) * amskt(ij, k) &
                  & + tx(ij, k) * (1.d0 - amskt(ij, k))
          end if
       end do
    end do
   
  end subroutine tmixss
  
end module tslvt



