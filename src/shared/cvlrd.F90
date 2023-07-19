module cvlrd

! --- information -----------------------------------------------------
!
!  Time integration of the equation of motion with semi-implicit
! treatment of the Coriolis term
!
!  HISTORY
!     '99.08.13  H.Hasumi: from CCSR2-MASK
!     '01.01.30  H.Hasumi: for partial step bottom topography
!     '01.05.02  H.Hasumi
!     '07.04.23  H.Hasumi
!     '08.06.11  H.Hasumi: initial/final processing
!     '08.07.10  H.Hasumi: initial/final processing
!     '12.08.05  H.Tatebe: for COCO5.0 in F90
!
! ---------------------------------------------------------------------

  use zocdim,  only :  nxydim,  nzdim

  implicit none

  private
  public  ::  velrds

  real(8),     save  ::     aa(nxydim,nzdim),      ac(nxydim,nzdim)
  real(8),     save  ::    abr(nxydim,nzdim),     abi(nxydim,nzdim)
  real(8),     save  ::    adr(nxydim,nzdim),     adi(nxydim,nzdim)
  real(8),     save  :: tdiffz(nxydim,nzdim)
  real(8),     save  ::  hvbot(nxydim)
  real(8),     save  ::    cf1,     cf2
  logical,     save  ::  ofirst
  data ofirst / .true. /

contains

  subroutine velrds(                                                  &
         &       ux,     vx,                                          &
         &       gx,     gy,    amv,                                  &
         &       hx   )

    use zocdim,  only :                                               &
         &  nxydim,  nzdim,                                           &
         &    kstr,   kend,     kz,                                   &
         &  ijvstr, ijvend,                                           &
         &      ln,     le,    lne,                                   &
         &   oinit,  ofinal
    use zocgrd,  only :                                               &
         &    zbot,     ts,    cor,                                   &  
         &      dy,    dym,    dzv,    dzm,    rs,     rsm
    use zocmsk,  only :  amskv,  amfvz
    use utrdg
    use ufile

    implicit none

    real(8),   intent(inout)  ::     ux(nxydim,nzdim),    vx(nxydim,nzdim)
    real(8),   intent(in)     ::     gx(nxydim,nzdim),    gy(nxydim,nzdim)
    real(8),   intent(in)     ::    amv(nxydim,nzdim),    hx(nxydim)

!---- local variables
    integer(4)         ::     ij,      k
    integer(4)         ::  ifpar,  jfpar,  istat

    real(8),    save   ::   acc,    aimp
    namelist /nmaccv/ acc
    namelist /nmimpl/ aimp
    data acc  / 1.0d0 /
    data aimp / 0.5d0 /

    if ( oinit .or. ofinal ) then
       return
    end if

    if ( ofirst ) then

       ofirst = .false.
       call rewnml(ifpar, jfpar)
       read(ifpar, nmaccv, iostat = istat )
       call cstnml( jfpar, 'velrds', 'nmaccv', istat )
       write( jfpar, nmaccv )

       call rewnml(ifpar, jfpar)
       read(ifpar, nmimpl, iostat = istat )
       call cstnml( jfpar, 'velrds', 'nmimpl', istat )
       write( jfpar, nmimpl )

       cf1 = aimp   / acc
       cf2 = 1.0d+0 / acc

    end if

!$omp parallel do
    do k = 1, nzdim
       do ij = 1, nxydim
          aa    (ij, k) = 0.d0
          abr   (ij, k) = 0.d0
          abi   (ij, k) = 0.d0
          ac    (ij, k) = 0.d0
          adr   (ij, k) = 0.d0
          adi   (ij, k) = 0.d0
          tdiffz(ij, k) = 0.d0
       end do
    end do
!$omp end parallel do
    
!$omp parallel do
    do ij = ijvstr, ijvend
       hvbot(ij) = (   (hx(ij)    + hx(ij+le) ) * dy(ij)              &
    &                + (hx(ij+ln) + hx(ij+lne)) * dy(ij+ln) )         &
    &              / dym(ij) * 0.25d0                                 &
    &             + zbot
       hvbot(ij) = 1.d0 / hvbot(ij)
    end do
!$omp end parallel do

!$omp parallel do
    do k = kstr+1, kstr+kz-1
       do ij = ijvstr, ijvend
          tdiffz(ij, k) = ts * cf2 * amv(ij, k) * hvbot(ij) * rsm(k) * amfvz(ij, k)
       end do
    end do
!$omp end parallel do

!$omp parallel do
    do k = kstr+kz, kend
       do ij = ijvstr, ijvend
          tdiffz(ij, k) = ts * cf2 * amv(ij, k) / dzm(ij, k) * amfvz(ij, k)
       end do
    end do
!$omp end parallel do

!$omp parallel do
    do k = kstr, kstr+kz-1
       do ij = ijvstr, ijvend
          aa (ij, k) = - tdiffz(ij, k  ) * hvbot(ij) * rs(k)
          ac (ij, k) = - tdiffz(ij, k+1) * hvbot(ij) * rs(k)
          abr(ij, k) = 1.d0 - aa(ij, k) - ac(ij, k)
          abi(ij, k) = ts * cf1 * cor(ij)
          adr(ij, k) = gx(ij, k)
          adi(ij, k) = gy(ij, k)
       end do
    end do
!$omp end parallel do

!$omp parallel do
    do k = kstr+kz, kend
       do ij = ijvstr, ijvend
          aa (ij, k) = - tdiffz(ij, k  ) / dzv(ij, k)
          ac (ij, k) = - tdiffz(ij, k+1) / dzv(ij, k)
          abr(ij, k) = 1.d0 - aa(ij, k) - ac(ij, k)
          abi(ij, k) = ts * cf1 * cor(ij)
          adr(ij, k) = gx(ij, k)
          adi(ij, k) = gy(ij, k)
       end do
    end do
!$omp end parallel do
    
    call thmasc( adr, adi, aa, abr, abi, ac )

!$omp parallel do
    do k = kstr, kend
       do ij = ijvstr, ijvend
          ux(ij, k) = (  ux(ij, k) + ts * cf2 * adr(ij, k)) * amskv(ij, k)
          vx(ij, k) = (  vx(ij, k) + ts * cf2 * adi(ij, k)) * amskv(ij, k)
       end do
    end do
!$omp end parallel do

  end subroutine velrds

end module cvlrd


