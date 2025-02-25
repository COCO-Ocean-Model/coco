module csrvl

! --- information -----------------------------------------------------
!
!  Estimate the bottom friction term of the equation of motion
!
!  HISTORY
!     '99.08.13  H.Hasumi: from CCSR2-MASK
!     '01.01.30  H.Hasumi: for partial step bottom topography
!     '07.04.23  H.Hasumi
!     '07.11.28  Y.Komuro: give different MZ for each hemisphere
!     '08.06.11  H.Hasumi: initial/final processing
!     '08.07.10  H.Hasumi: initial/final processing
!     '12.08.02  Y.Komuro: for COCO5.0
!
! ---------------------------------------------------------------------

  use zocdim, only: &
    & nxydim,  nzdim,   kstr,   kend,     nz, &
    & ijvstr, ijvend, &
    &  oinit, ofinal
  use zocgrd, only: &
    &    dzv,    cor
  use zocmsk, only: &
#ifdef OPT_BBL
    & amskvb,  nbotv, &
#endif
    &  amskb

  implicit none

  private

  public :: srcvel
#ifdef OPT_BBL
  public :: srcvlb
#endif

contains
  
subroutine srcvel( &
  &                    gx,     gy,     xx,     yy, &
  &                    ux,     vx)
  use ufile

  real(8), intent(out) ::     gx(nxydim, nzdim),     gy(nxydim, nzdim)
  real(8), intent(out) ::     xx(nxydim, nzdim),     yy(nxydim, nzdim)
  real(8), intent(in)  ::     ux(nxydim, nzdim),     vx(nxydim, nzdim)

  real(8) ::    abv
  integer ::     ij,      k
  integer ::  ifpar,  jfpar,  istat

  logical, save :: ofirst = .true.

  real(8), save :: btmfrc = 0.0d0

  namelist /nmbtmf/ btmfrc

  if (oinit .or. ofinal) then
     return
  end if

  if (ofirst) then
     ofirst = .false.
     call rewnml(ifpar, jfpar)
     write(jfpar, *) '*** svlset  ***'
     read (ifpar, nmbtmf, iostat=istat)
     call cstnml(jfpar, 'srcvel', 'nmbtmf', istat)
     write(jfpar, nmbtmf)
  end if

#ifdef ADF_
!$acc data copyin(amskb, dzv) 
!$acc data copy(gx, gy, xx, yy) 
!$acc data copyin(ux, vx) 
#endif

#ifdef ACC_
!$acc data copy(gx, gy, xx, yy)
!$acc kernels 
#elif  OMP_
!$omp parallel do private(k, ij)
#endif
  do k = 1, nzdim
     do ij = 1, nxydim
        gx(ij, k) = 0.d0
        gy(ij, k) = 0.d0
        xx(ij, k) = 0.d0
        yy(ij, k) = 0.d0
     end do
  end do
#ifdef ACC_
!$acc end kernels
!$acc end data ! copy(gx, gy, xx, yy)
#elif  OMP_
!$omp end parallel do
#endif

#ifdef ACC_
!$acc data copyin(vx, ux) 
!$acc data copy(gy, gx, yy, xx) 
!$acc data copyin(amskb, dzv) 
!$acc kernels 
#elif  OMP_
!$omp parallel do private(k, ij, abv)
#endif
  do k = kstr, kend
     do ij = ijvstr, ijvend
        abv = - btmfrc / dzv(ij, k) * &
          &     sqrt(ux(ij, k) * ux(ij, k) + vx(ij, k) * vx(ij, k))
        gx(ij, k) = abv * ux(ij, k) * amskb(ij, k)
        gy(ij, k) = abv * vx(ij, k) * amskb(ij, k)
        xx(ij, k) = gx(ij, k)
        yy(ij, k) = gy(ij, k)
     end do
  end do
#ifdef ACC_
!$acc end kernels
!$acc end data ! copyin(amskb, dzv)
!$acc end data ! copy(gy, gx, yy, xx)
!$acc end data ! copyin(vx, ux)
#elif  OMP_
!$omp end parallel do
#endif

#ifdef ADF_
!$acc end data ! copyin(ux, vx)
!$acc end data ! copy(gx, gy, xx, yy)
!$acc end data ! copyin(amskb, dzv)
#endif

  return

end subroutine srcvel
#ifdef OPT_BBL
! *********************************************************************

subroutine srcvlb( &
  &                   gx,     gy,     xx,     yy, &
  &                   ux,     vx)

! --- information -----------------------------------------------------
!
!  Bottom friction term for the BBL momemtum equations. See Nakano and
! Suginohara (2002, JPO).
!
!  HISTORY
!     '01.02.08  H.Hasumi
!     '12.08.02  Y.Komuro: for COCO5.0
!
! ---------------------------------------------------------------------
  use ufile

  real(8), intent(out) ::      gx(nxydim, nzdim),     gy(nxydim, nzdim)
  real(8), intent(out) ::      xx(nxydim, nzdim),     yy(nxydim, nzdim)
  real(8), intent(in)  ::      ux(nxydim, nzdim),     vx(nxydim, nzdim)

#ifndef ADC_SRCVLB
  real(8), save ::    rfrc(nxydim)
#endif

  real(8) ::    abv
  integer ::     ij,      k
  integer ::  ifpar,  jfpar,  istat

  logical, save :: ofirst = .true.

  real(8), save :: btmfrc = 0.0d0,  rayfrc = 1.0d0
  integer, save :: mzn = nz, mzs = nz

  namelist /nmbtmf/ btmfrc
  namelist /nmbbrf/ rayfrc, mzn, mzs

  if (oinit .or. ofinal) then
     return
  end if

  if (ofirst) then
     ofirst = .false.
     call rewnml(ifpar, jfpar)
     write(jfpar, *) '*** svbset  ***'
     read (ifpar, nmbtmf, iostat=istat)
     call cstnml(jfpar, 'srcvlb', 'nmbtmf', istat)
     write(jfpar, nmbtmf)
     call rewnml(ifpar, jfpar)
     read (ifpar, nmbbrf, iostat=istat)
     call cstnml(jfpar, 'srcvlb', 'nmbbrf', istat)
     write(jfpar, nmbbrf)

     do ij = 1, nxydim
        if (((cor(ij).ge.0.0d0).and.(nbotv(ij).le.mzn+kstr-1)) .or. &
          & ((cor(ij).lt.0.0d0).and.(nbotv(ij).le.mzs+kstr-1))) then
           rfrc(ij) = rayfrc * abs(cor(ij))
        else
           rfrc(ij) = 0.d0
        end if
     end do
#ifdef ADC_SRCVLB
!$acc update device(rfrc)
#endif
  end if

#ifdef ADF_
!$acc data copy(xx, yy, gx, gy) 
!$acc data copyin(ux, vx) 
!$acc data copyin(amskvb, dzv) 
!$acc data copyin(rfrc) 
#endif

#ifdef ACC_
!$acc kernels 
#elif  OMP_
!$omp parallel do private(ij, abv)
#endif
  do ij = ijvstr, ijvend
     abv = - btmfrc / dzv(ij, kend) * &
       &     sqrt(  ux(ij, kend) * ux(ij, kend) &
       &          + vx(ij, kend) * vx(ij, kend)) &
       &   - rfrc(ij)
     gx(ij, kend) = abv * ux(ij, kend) * amskvb(ij)
     gy(ij, kend) = abv * vx(ij, kend) * amskvb(ij)
     xx(ij, kend) = gx(ij, kend)
     yy(ij, kend) = gy(ij, kend)
  end do
#ifdef ACC_
!$acc end kernels
#elif  OMP_
!$omp end parallel do
#endif

#ifdef ADF_
!$acc end data ! copy(xx, yy, gx, gy) 
!$acc end data ! copyin(ux, vx)
!$acc end data ! copyin(amskvb, dzv)
!$acc end data ! copyin(rfrc)
#endif

  return

end subroutine srcvlb
#endif

end module csrvl
