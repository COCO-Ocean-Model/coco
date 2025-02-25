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

#ifdef ADF_
!$acc data copyin(dy, dym) 
!$acc data copyin(rsm) 
!$acc data copyin(amfvz) 
!$acc data copyin(dzm) 
!$acc data copyin(rs) 
!$acc data copyin(cor) 
!$acc data copyin(dzv) 
!$acc data copyin(amskv) 
!$acc data copy(ux, vx) 
!$acc data copyin(gx, gy, amv, hx) 
!$acc data copy(tdiffz) 
!$acc data copy(hvbot) 
#ifdef ADF_
!$acc data copy(aa, ac) 
!$acc data copy(abr, abi) 
!$acc data copy(adr, adi) 
#endif ! ADF_
#endif

#ifdef ACC_
!$acc data copy(tdiffz, adr, adi, aa, abi, abr, ac) 
!$acc kernels 
#elif  OMP_
!$omp parallel do private(k, ij)
#endif
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
#ifdef ACC_
!$acc end kernels
!$acc end data ! copy(tdiffz, adr, adi, aa, abi, abr, ac)
#elif  OMP_
!$omp end parallel do
#endif
    
#ifdef ACC_
!$acc data copyin(hx) 
!$acc data copyin(dy, dym) 
!$acc data copy(hvbot) 
!$acc kernels 
#elif  OMP_
!$omp parallel do private(ij)
#endif
    do ij = ijvstr, ijvend
       hvbot(ij) = (   (hx(ij)    + hx(ij+le) ) * dy(ij)              &
    &                + (hx(ij+ln) + hx(ij+lne)) * dy(ij+ln) )         &
    &              / dym(ij) * 0.25d0                                 &
    &             + zbot
       hvbot(ij) = 1.d0 / hvbot(ij)
    end do
#ifdef ACC_
!$acc end kernels
!$acc end data ! copyin(hx)
!$acc end data ! copyin(dy, dym)
!$acc end data ! copy(hvbot)
#elif  OMP_
!$omp end parallel do
#endif

#ifdef ACC_
!$acc data copyin(hvbot) 
!$acc data copy(tdiffz) 
!$acc data copyin(amv) 
!$acc data copyin(rsm, amfvz) 
!$acc kernels 
#elif  OMP_
!$omp parallel do private(k, ij)
#endif
    do k = kstr+1, kstr+kz-1
       do ij = ijvstr, ijvend
          tdiffz(ij, k) = ts * cf2 * amv(ij, k) * hvbot(ij) * rsm(k) * amfvz(ij, k)
       end do
    end do
#ifdef ACC_
!$acc end kernels
!$acc end data ! copy(tdiffz)
!$acc end data ! copyin(hvbot)
!$acc end data ! copyin(amv)
!$acc end data ! copyin(rsm, amfvz)
#elif  OMP_
!$omp end parallel do
#endif

#ifdef ACC_
!$acc data copy(tdiffz) 
!$acc data copyin(amv) 
!$acc data copyin(amfvz, dzm) 
!$acc kernels 
#elif  OMP_
!$omp parallel do private(k, ij)
#endif
    do k = kstr+kz, kend
       do ij = ijvstr, ijvend
          tdiffz(ij, k) = ts * cf2 * amv(ij, k) / dzm(ij, k) * amfvz(ij, k)
       end do
    end do
#ifdef ACC_
!$acc end kernels
!$acc end data ! copy(tdiffz)
!$acc end data ! copyin(amv)
!$acc end data ! copyin(amfvz, dzm)
#elif  OMP_
!$omp end parallel do
#endif

#ifdef ACC_
!$acc data copyin(tdiffz, hvbot) 
!$acc data copy(adi, abi, adr, abr, aa, ac) 
!$acc data copyin(cor, rs) 
!$acc data copyin(gy, gx) 
!$acc kernels 
#elif  OMP_
!$omp parallel do private(k, ij)
#endif
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
#ifdef ACC_
!$acc end kernels
!$acc end data ! copyin(gy, gx)
!$acc end data ! copyin(cor, rs)
!$acc end data ! copyin(tdiffz, hvbot)
!$acc end data ! copy(adi, abi, adr, abr, aa, ac)
#elif  OMP_
!$omp end parallel do
#endif

#ifdef ACC_
!$acc data copy(adi, abi, adr, abr, aa, ac) 
!$acc data copyin(tdiffz) 
!$acc data copyin(gy, gx) 
!$acc data copyin(dzv, cor) 
!$acc kernels 
#elif  OMP_
!$omp parallel do private(k, ij)
#endif
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
#ifdef ACC_
!$acc end kernels
!$acc end data ! copy(adi, abi, adr, abr, aa, ac)
!$acc end data ! copyin(tdiffz)
!$acc end data ! copyin(gy, gx)
!$acc end data ! copyin(dzv, cor)
#elif  OMP_
!$omp end parallel do
#endif
    
#ifdef ADF_
!$acc data copyin(aa, ac) 
!$acc data copy(adr, adi, adr, abi) 
#endif
    call thmasc( adr, adi, aa, abr, abi, ac )
#ifdef ADF_
!$acc end data ! copyin(aa, ac)
!$acc end data ! copy(adr, adi, adr, abi)  
#endif

#ifdef ACC_
!$acc data copy(vx, ux) 
!$acc data copyin(adi, adr) 
!$acc data copyin(amskv) 
!$acc kernels 
#elif  OMP_
!$omp parallel do private(k, ij)
#endif
    do k = kstr, kend
       do ij = ijvstr, ijvend
          ux(ij, k) = (  ux(ij, k) + ts * cf2 * adr(ij, k)) * amskv(ij, k)
          vx(ij, k) = (  vx(ij, k) + ts * cf2 * adi(ij, k)) * amskv(ij, k)
       end do
    end do
#ifdef ACC_
!$acc end kernels
!$acc end data ! copy(vx, ux)
!$acc end data ! copyin(adi, adr)
!$acc end data ! copyin(amskv)
#elif  OMP_
!$omp end parallel do
#endif

#ifdef ADF_
!$acc end data ! copyin(dy, dym)
!$acc end data ! copyin(rsm)
!$acc end data ! copyin(amfvz)
!$acc end data ! copyin(dzm)
!$acc end data ! copyin(rs)
!$acc end data ! copyin(cor)
!$acc end data ! copyin(dzv)
!$acc end data ! copyin(amskv)
!$acc end data ! copy(ux, vx)
!$acc end data ! copyin(gx, gy, amv, hx)
!$acc end data ! copy(tdiffz)
!$acc end data ! copy(hvbot)
#ifdef ADF_
!$acc end data ! copy(aa, ac)
!$acc end data ! copy(abr, abi)
!$acc end data ! copy(adr, adi)
#endif
#endif

  end subroutine velrds

end module cvlrd


