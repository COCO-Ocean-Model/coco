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
!     '14.10.15  M.Kurogi: previous form for nesting code
!     '15.04.02  M.Kurogi: for vector machine
! ---------------------------------------------------------------------

  use zocdim, only :                                   &
       &   nxydim ,  nzdim,                            &
       &   ijtstr,  ijtend,                            &
       &     kstr,    kend,     nz
  use zocmsk,  only :                                  &
       &     nbot

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

  real(8),    save  ::    sdmp
  logical,    save  ::  osrstr, osrsti

  real(8)           ::     dzb(nxydim),    dzt(nxydim)
  real(8)           ::   radup,  raddn,  depth
  real(8)           ::  tswcnv
  real(8)           ::     rrr,  zeta1,  zeta2
  integer(4)        ::   ifpar,  jfpar,  istat

  public  ::  slvtrc     !   aprdc.F
  public  ::  svtset     !   aocea.F

  namelist /nmacct/ gamma
  namelist /nmswab/   rrr,  zeta1,  zeta2
  namelist /nmsrst/  sdmp, osrstr, osrsti

!  data gamma / nz*1.d0 /
  data rrr, zeta1, zeta2 / 5.8d-1, 3.5d+1, 2.3d+3 /
  data sdmp / 0.d0 /
  data osrstr, osrsti / .false., .false. /

contains

! **********************************************************************

  subroutine svtset

    use zocdim,  only :                                &
         &  nxydim,  nzdim,  ntdim,                    &
         &    kstr,   kend,     kz,                    &
         &   ijstr,  ijend, ijtstr,  ijtend,           &
         &     nic,                                    &
         &   oinit,  ofinal
    use zocgrd,  only :                                &
         &      dz,    dz0,     ds,   zbot,            &
         &      ts
    use zocmsk,  only :                                &
#ifdef OPT_BBL
         & amsktb
#endif
         &   amskt,  nbot
    use zocphy,  only :                                &
         &     cpo,    rhoo

    use ifhea
    use ufile

    implicit none

    integer(4)        ::      ij,      k,      n
    logical,   save   ::  ofirst = .true.

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
       write(jfpar, nmswab)
       call rewnml(ifpar, jfpar)
       read(ifpar, nmsrst, iostat = istat )
       write(jfpar, nmsrst)

       if ( osrstr ) then
          sdmp = 1.d0 / 86400.d0 / sdmp
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
       
       do ij = 1, nxydim
          swconv(ij, kstr) = 0.d0
       end do
       
       do k = kstr+1, kend
          do ij = ijtstr, ijtend
             swconv(ij, k) = swconv(ij, k) / rhoo / cpo
          end do
       end do
       
    end if

  end subroutine svtset

! **********************************************************************

  subroutine slvtrc(                                   &
         &      tx,     hx,                            &
         &     adt,  diffz,                            &
         &      ft,  swabs,     fs,     hz,   ssfc,    &
         &      ax  )

    use zocdim,  only :                                &
         &  nxydim,  nzdim,  ntdim,                    &
         &    kstr,   kend,     kz,                    &
         &   ijstr,  ijend, ijtstr,  ijtend,           &
         &     nic,                                    &
         &   oinit,  ofinal
    use zocgrd,  only :                                &
         &      dz,    dz0,     ds,   zbot,            &
         &      ts
    use zocmsk,  only :                                &
#ifdef OPT_BBL
         & amsktb
#endif
         &   amskt,  nbot
    use utrdg

    implicit none

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

    integer(4)                ::     ij,    k,     n
    
    if ( oinit .or. ofinal ) then
       return
    end if
  
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
          k = max(nbot(ij), k)
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

    do ij = ijtstr, ijtend
       hzbot(ij) = hx(ij) + zbot
       dh(ij) = - ts * ft(ij, 2) * amskt(ij, kstr)
       hx(ij) = hx(ij) + dh(ij)
       hxbot(ij) = hx(ij) + zbot
       if (dh(ij) .ge. 0.d0) then
          dzb(ij) = 0.d0
          dzt(ij) = dh(ij) * ds(kstr+kz-1)
       else
          dzt(ij) = -dh(ij)
          dzb(ij) = dzt(ij) + dh(ij) * ds(kstr)
       end if
    end do

    do k = kstr+kz-1, kstr+1, -1
       do ij = ijtstr, ijtend
          if (dh(ij) .ge. 0.d0) then
             do n = 1, ntdim
                tx(ij, k, n) = (  tx(ij, k  , n) * ds(k) * hzbot(ij) &
    &                           - tx(ij, k  , n) * dzb(ij)           &
    &                           + tx(ij, k-1, n) * dzt(ij))          &
    &                          / hxbot(ij) / ds(k)
             end do
             dzb(ij) = dzt(ij)
             dzt(ij) = dzb(ij) + dh(ij) * ds(k-1)
          end if
       end do
    end do
    do ij = ijtstr, ijtend
       if (dh(ij) .ge. 0.d0) then
          tx(ij, kstr, 2) = (  tx(ij, kstr, 2) * ds(kstr) * hzbot(ij) &
    &                        - tx(ij, kstr, 2) * dzb(ij))             &
    &                       / hxbot(ij) / ds(kstr)
          tx(ij, kstr, 1) = (  tx(ij, kstr, 1) * ds(kstr) * hzbot(ij) &
    &                        - tx(ij, kstr, 1) * dzb(ij)              &
    &                        + tx(ij, kstr, 1) * dzt(ij)              &
    &                        + ts * ft(ij, 1))                        &
    &                       / hxbot(ij) / ds(kstr)
          do n = 3, ntdim
             tx(ij, kstr, n) = (  tx(ij, kstr, n)*ds(kstr)*hzbot(ij) &
    &                           - tx(ij, kstr, n) * dzb(ij)          &
    &                           + tx(ij, kstr, n) * dzt(ij)          &
    &                           + ts * ft(ij, n))                    &
    &                          / hxbot(ij) / ds(kstr)
          end do
       else
          tx(ij, kstr, 2) = (  tx(ij, kstr, 2) * ds(kstr) * hzbot(ij) &
    &                        + tx(ij, kstr+1, 2) * dzb(ij))           &
    &                       / hxbot(ij) / ds(kstr)
          tx(ij, kstr, 1) = (  tx(ij, kstr, 1) * ds(kstr) * hzbot(ij) &
    &                        + tx(ij, kstr+1, 1) * dzb(ij)            &
    &                        - tx(ij, kstr, 1) * dzt(ij)              &
    &                        + ts * ft(ij, 1))                        &
    &                       / hxbot(ij) / ds(kstr)
          do n = 3, ntdim
             tx(ij, kstr, n) = (  tx(ij, kstr, n)*ds(kstr)*hzbot(ij) &
    &                           + tx(ij, kstr+1, n) * dzb(ij)        &
    &                           - tx(ij, kstr  , n) * dzt(ij)        &
    &                           + ts * ft(ij, n))                    &
    &                          / hxbot(ij) / ds(kstr)  
          end do
       end if
   end do

   do k = kstr+1, kstr+kz-2
      do ij = ijtstr, ijtend
         if (dh(ij) .lt. 0.d0) then
             dzt(ij) = dzb(ij)
             dzb(ij) = dzt(ij) + dh(ij) * ds(k)
             do n = 1, ntdim
                tx(ij, k, n) = (  tx(ij, k, n) * ds(k) * hzbot(ij) &
    &                           + tx(ij, k+1, n) * dzb(ij)         &
    &                           - tx(ij, k  , n) * dzt(ij))        &
    &                          / hxbot(ij) / ds(k)
             end do
          end if
       end do
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

    do ij = ijtstr, ijtend
       tx(ij, kstr, 2) = tx(ij, kstr, 2)                             &
    &                  - ts * fs(ij) / hxbot(ij) / ds(kstr)
    end do

#ifdef OPT_SRST
    if ( osrstr ) then
       if ( osrsti ) then
          do ij = ijtstr, ijtend
             tx(ij, kstr, 2)                                          &
           &       = tx(ij, kstr, 2)                                  &
           &       + ts * sdmp * (ssfc(ij) - tx(ij, kstr, 2))         &
           &       * amskt(ij, kstr)
          end do
       else
          do ij = ijtstr, ijtend
             if ( ax(ij, 0) == 1.d0 ) then
                tx(ij, kstr, 2)                                       &
          &           = tx(ij, kstr, 2)                               &
          &           + ts * sdmp * (ssfc(ij) - tx(ij, kstr, 2))      &
          &           * amskt(ij, kstr)
             end if
          end do
       end if
    end if
#endif

  end subroutine slvtrc
  
end module tslvt



