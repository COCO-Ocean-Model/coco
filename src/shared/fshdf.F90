module fshdf

! --- information -----------------------------------------------------
!
!  HISTORY
!     '03.04.24  H.Hasumi: from COCO3.4
!     '07.04.23  H.Hasumi
!     '08.06.11  H.Hasumi: initial/final processing
!     '08.07.10  H.Hasumi: initial/final processing
!     '12.06.29  H.Tatebe: for COCO5.0 in F90
! ---------------------------------------------------------------------

  use zocdim,  only  :  nxydim,  nzdim,  ntdim

  implicit none

  private
  public  ::  shdiff  !  used in aprdc.F

  real(8),     save  ::    tsh(nxydim,nzdim,ntdim)
  real(8),     save  ::    fhx(nxydim),    fhy(nxydim)
  real(8),     save  ::    ftx(nxydim),    fty(nxydim)
  real(8),     save  :: rhxbot(nxydim)

  real(8),     save  ::    ash
  data ash / 0.d0 /

  logical,     save  ::  ofirst
  data ofirst / .true. /

contains
  
  subroutine shdiff( tx, hx )

    use zocdim,  only :                                               &
         &  nxydim,  nxdim,  nzdim,                                   &
         &    kstr,   kend,     kz,                                   &
         &  ijtstr, ijtend,                                           &
         &      le,     lw,     ln,     ls,                           &
         &     lnw,    lne,    lse,    lsw,                           &
         &   oinit,  ofinal
    use zocgrd,  only :                                               &
         &    zbot,     ts,                                           &  
         &      rx,     ry,    rym,                                   & 
         &     hxu,    hyu,    hxt,    hyt,                           &
         &     rxt,    ryt
    use zocmsk,  only  :  amskt
    use ufile
    
    implicit none

    real(8),   intent(inout)  ::     hx(nxydim),    tx(nxydim,nzdim,ntdim)

!---- local
    integer(4)         ::     ij,      k,      n
    integer(4)         ::  ifpar,  jfpar,  istat

    namelist /nmdfsh/ ash

    
    if ( oinit .or. ofinal ) then
       return
    end if

    if (ofirst) then
       ofirst = .false.
       call rewnml( ifpar, jfpar )
       read(ifpar, nmdfsh, iostat = istat )
       call cstnml( jfpar, 'shdiff', 'nmdfsh', istat )
       write( jfpar, nmdfsh )
    end if

    do n = 1, ntdim
       do k = kstr, kstr+kz-1
          do ij = 1, nxydim
             tsh(ij, k, n) = tx(ij, k, n) * ( hx(ij) + zbot )
          end do
       end do
    end do

    do ij = ijtstr, ijtend+nxdim
       fhx(ij) = ash * (hx(ij) - hx(ij+lw)) * rx *                    &
    &             ( hyu(ij+lw) + hyu(ij+lsw) )                        &
    &           / ( hxt(ij)    + hxt(ij+lw)  ) *                      &
    &           amskt(ij, kstr) * amskt(ij+lw, kstr)
       fhy(ij) = ash * (hx(ij) - hx(ij+ls)) * rym(ij+ls) *            &
    &             ( hxu(ij+ls) + hxu(ij+lsw) )                        &
    &           / ( hyt(ij)    + hyt(ij+ls)  ) *                      &
    &           amskt(ij, kstr) * amskt(ij+ls, kstr)
    end do

    do ij = ijtstr, ijtend
       hx(ij) = hx(ij)                                                &
    &         + ts * (  (fhx(ij+le) - fhx(ij)) * rx                   & 
    &                 + (fhy(ij+ln) - fhy(ij)) * ry(ij)) *            &
    &           rxt(ij) * ryt(ij) * amskt(ij, kstr)
    end do

    do ij = ijtstr, ijtend
       rhxbot(ij) = 1.d0 / ( hx(ij) + zbot )
    end do

    do n = 1, ntdim
       do k = kstr, kstr+kz-1
          do ij = ijtstr, ijtend+nxdim
             ftx(ij) =                                                &
    &         (  (fhx(ij) + abs(fhx(ij))) * tx(ij, k, n)              &
    &          + (fhx(ij) - abs(fhx(ij))) * tx(ij+lw, k, n)) *        &
    &         0.5d0 * amskt(ij, k) * amskt(ij+lw, k)               
             fty(ij) =                                                &
    &         (  (fhy(ij) + abs(fhy(ij))) * tx(ij, k, n)              &
    &          + (fhy(ij) - abs(fhy(ij))) * tx(ij+ls, k, n)) *        &
    &         0.5d0 * amskt(ij, k) * amskt(ij+ls, k)
          end do

          do ij = ijtstr, ijtend
             tsh(ij, k, n) =                                          &
    &        tsh(ij, k, n)                                            &
    &        + ts * (  (ftx(ij+le) - ftx(ij)) * rx                    &
    &                + (fty(ij+ln) - fty(ij)) * ry(ij)) *             &
    &          rxt(ij) * ryt(ij) * amskt(ij, k)
          end do

          do ij = ijtstr, ijtend
             tx(ij, k, n) = tsh(ij, k, n) * rhxbot(ij)
          end do

       end do
    end do

  end subroutine shdiff
 
end module fshdf


