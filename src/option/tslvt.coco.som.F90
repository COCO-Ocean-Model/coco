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
!     '12.09.20  H.Tatebe: for COCO5.0 in F90
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

  real(8),    save  ::      aa(nxydim, nzdim)
  real(8),    save  ::      ab(nxydim, nzdim)
  real(8),    save  ::      ac(nxydim, nzdim)

  real(8),    save  ::      dh(nxydim)  
  real(8),    save  ::   hzbot(nxydim),   hxbot(nxydim)
  real(8),    save  ::       gamma(nz)

  real(8)           ::     dzb,    dzt
  integer(4)        ::   ifpar,  jfpar,  istat

  public  ::  slvtrc     !   aprdc.F
  public  ::  svtset     !   aocea.F

  namelist /nmacct/ gamma

contains

! **********************************************************************

  subroutine svtset

    use zocdim,  only  :    kstr,    kend
    use ufile

    implicit none

    integer(4)        ::       k
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

       do k = kstr, kend
          rgamma(k) = 1.d0 / gamma(k-kstr+1)
       end do
  
    end if

  end subroutine svtset

! **********************************************************************

  subroutine slvtrc(                                   &
         &      tx,     hx,                            &
         &     adt,  diffz,                            &
         &      ft,     hz  )

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
         &   amskt,  amsktb,  nbot
    use utrdg

    implicit none

    real(8),   intent(inout)  ::     tx(nxydim, nzdim, ntdim)
    real(8),   intent(inout)  ::     hx(nxydim)
    real(8),   intent(inout)  ::    adt(nxydim, nzdim, ntdim)
    real(8),   intent(in)     ::  diffz(nxydim, nzdim)
    real(8),   intent(in)     ::     ft(nxydim, ntdim)
    real(8),   intent(in)     ::     hz(nxydim)

    integer(4)                ::     ij,    k,     n
    
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
    
!---- when SOM is used, the lines below should be commented.
!    do ij = 1, nxydim
!       dh(ij) = hx(ij) - hz(ij)
!    end do
!----
    
    do k = kstr, kstr+kz-1
       do ij = ijtstr, ijtend
          aa(ij, k) = - ts * diffz(ij, k) / dz(ij, k)
          ac(ij, k) = - ts * diffz(ij, k+1) / dz(ij, k)
          ab(ij, k) = (hx(ij) + zbot) / zbot - aa(ij, k) - ac(ij, k)
       end do
    end do

!---- when SOM is used, the lines below should be commented.    
!    do n = 1, ntdim
!       do k = kstr, kstr+kz-1
!          do ij = ijtstr, ijtend
!             adt(ij, k, n) = adt(ij, k, n)                            &
!    &                      - dh(ij) / zbot / ts * tx(ij, k, n)
!          end do
!       end do
!    end do
!----

    do k = kstr+kz, kend
       do ij = ijtstr, ijtend
          aa(ij, k) = - ts * diffz(ij, k) / dz(ij, k)
          ac(ij, k) = - ts * diffz(ij, k+1) / dz(ij, k)
          ab(ij, k) = 1.d0 - aa(ij, k) - ac(ij, k)
       end do
    end do
    
    do ij = ijstr, ijend
       adt(ij, kstr, 1) = adt(ij, kstr, 1)                            &
    &                   + tx(ij, kstr, 1) * ft(ij, 2) / zbot          &
    &                    * amskt(ij, kstr)
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

    do ij = ijtstr, ijtend
       tx(ij, kstr, 1) = tx(ij, kstr, 1)                              &
    &                  + ts * ft(ij, 1) / hxbot(ij) / ds(kstr) 
    end do

    do n = 3, ntdim
       do ij = ijtstr, ijtend
          tx(ij, kstr, n) = tx(ij, kstr, n)                          &
    &                     + ts * ft(ij, n) / hxbot(ij) / ds(kstr)
       end do
    end do


  end subroutine slvtrc
  
end module tslvt



