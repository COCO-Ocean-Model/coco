module sfcng

! --- information -----------------------------------------------------
!
!  Estimate the surface forcing fluxes from the data given in the
! boundary condition files
!
!  HISTORY
!     '99.08.18  H.Hasumi
!     '01.01.30  H.Hasumi: for partial step bottom topography
!     '02.05.29  H.Nakano: tracer dimension
!     '02.06.02  H.Hasumi
!     '07.04.23  H.Hasumi
!     '12.09.20  H.Tatebe: Rewrite in F95 format
!
! ---------------------------------------------------------------------

  implicit none
  private

  public  ::  sfcflx
#ifdef OPT_BODY
  public  ::  bdyflx
#endif

contains

  subroutine sfcflx(                                                  &
    &          ft,   taux,   tauy,                                    &
    &           t  )

    use zocdim,  only  :  nxydim, nzdim, ntdim, ijstr, ijend,         &
                            kstr
    use zocgrd,  only  :      dz
    use zocmsk,  only  :   amskt
    use utint

    implicit none

    real(8),    intent(inout)  ::     ft(nxydim,ntdim)
    real(8),    intent(inout)  ::   taux(nxydim)
    real(8),    intent(inout)  ::   tauy(nxydim)
    real(8),    intent(in)     ::      t(nxydim,nzdim,ntdim)

!---- local variables
    real(8)        ::    tsfc(nxydim, ntdim), tdmp(nxydim, ntdim)
    integer(4)     ::      ij,      l,      n

    call tmintp(  taux, 1 )
    call tmintp(  tauy, 2 )

    do l = 1, ntdim
       n = l + 2
       call tmintp( tsfc(1, l), n )
    end do
    do l = 1, ntdim
       n = l + 2 + ntdim
       call tmintp( tdmp(1, l), n )
    end do

    do ij = ijstr, ijend
       ft(ij, 1) = tdmp(ij, 1) * (tsfc(ij, 1) - t(ij, kstr, 1)) *     &
    &              dz(ij, kstr) * amskt(ij, kstr)
       ft(ij, 2) = tdmp(ij, 2) * (tsfc(ij, 2) - t(ij, kstr, 2)) *     &
    &              dz(ij, kstr) / t(ij, kstr, 2) * amskt(ij, kstr)
    end do

    do l = 3, ntdim
       do ij = ijstr, ijend
          ft(ij, l) = tdmp(ij, l) * (tsfc(ij, l) - t(ij, kstr, l)) *  &
   &                  dz(ij, kstr) * amskt(ij, kstr)
       end do
    end do

  end subroutine sfcflx

! --- information -----------------------------------------------------
!
!  Estimate the surface forcing fluxes from the data given in the
! boundary condition files
!
!  HISTORY
!     '99.08.25  H.Hasumi
!     '02.06.02  H.Hasumi: tracer dimension
!     '12.09.20  H.Tatebe: Rewrite in F95 format
! ---------------------------------------------------------------------

#ifdef OPT_BODY

  subroutine bdyflx(  tq,  t )
    
    use zocdim,   only  :                                             &
    &      nxydim,   nzdim,   ntdim,                                  &
    &      ijtstr,  ijtend,    kstr,   kend
    use zocmsk,   only  :                                             &
    &       amskt
    use utint

    implicit none

    real(8),    intent(inout)  ::  tq(nxydim, nzdim, ntdim)
    real(8),    intent(inout)  ::  t (nxydim, nzdim, ntdim)


    real(8)      ::    tbdy(nxydim, nzdim, ntdim)
    real(8)      ::    tdmb(nxydim, nzdim, ntdim)
    integer(4)   ::      ij,      k,      l,      n

    do l = 1, 2
       n = l
       call tmintb(  tbdy(1, 1, l),  n  )
    end do

    do l = 1, 2
       n = l + 2
       call tmintb(  tdmb(1, 1, l),  n  )
    end do

    do l = 1, 2
       do k = kstr, kend
          do ij = ijtstr, ijtend
             tq(ij, k, l) = tdmb(ij, k, l) *                          &
    &                      (tbdy(ij, k, l) - t(ij, k, l)) * amskt(ij, k)
          end do
       end do
    end do

  end subroutine bdyflx

#endif

end module sfcng
