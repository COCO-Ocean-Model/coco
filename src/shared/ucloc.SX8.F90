module ucloc

! --- information -----------------------------------------------------
!
!  Monitoring CPU time on Sun OS (Sun FORTRAN, using the intrinsic
! function "ETIME")
!
!  HISTORY
!     '97.03.18  H.Hasumi: From AGCM5.4 developed by A.Numaguti
!     '07.04.23  H.Hasumi
!     '12.12.06  H.Tatebe: for COCO5.0 in F90
!
! ---------------------------------------------------------------------

  implicit none
  private

  integer(4),    parameter          ::  nclmax = 100
  real(8),                    save  ::  cputim(nclmax), vputim(nclmax)
  real(8),                    save  ::  cpuold(nclmax), vpuold(nclmax)
  real(8),                    save  ::  cput,   vput        
  integer(4),                 save  ::  ifpar,  jfpar
  integer(4),                 save  ::  nclock
  integer(4),                 save  ::  ofirst
  character(16),              save  ::  htitle(nclmax)
                              
!---- used in yclock
  real(8),                    save  ::  cput2,  vput2
  real(8),                    save  ::  cput0

  public  ::  clcout,  clcstr,  clcend  !  ued in icedcoco.F, aprdc.F, iprdc,F
  
  data nclock / 0 /
  data ofirst / .true. /
  data htitle / nclmax*' ' /
  data cput0 / 0.d0 /

contains

  subroutine clcout

    use ufile

    implicit none

    integer(4)  ::  ic

    call rewnml(ifpar, jfpar)

    do ic = 1, nclock
       if ( htitle(ic) /= ' ' ) then
          write(jfpar, '(1x,a16,2f15.6)') htitle(ic), cputim(ic), vputim(ic)
       end if
    end do
    
    call yclock(cput, vput)
    write(jfpar, '(1x,a16,2f15.6)') ' total time = ', cput, vput
    
  end subroutine clcout
  
! **********************************************************************

  subroutine clcstr( httl )

    implicit none

    character(*),           intent(in)     ::  httl

    integer(4)  ::  ic

    if (ofirst) then
       ofirst = .false.
       call yclocl
    end if

    call yclock( cput, vput )
    do ic = 1, nclock
       if ( htitle(ic) == httl ) then
          cpuold(ic) = cput
          vpuold(ic) = vput
          return
       end if
    end do
    if ( nclock <= nclmax ) then
       nclock = nclock + 1
       htitle(nclock) = httl
       cpuold(nclock) = cput
       vpuold(nclock) = vput
    end if
    
  end subroutine clcstr

! **********************************************************************

  subroutine clcend( httl )

    implicit none

    character(*),           intent(in)     ::  httl

    integer(4)  ::  ic

    call yclock( cput, vput )
    do ic = 1, nclock
       if (htitle(ic) == httl) then
          cputim(ic) = cputim(ic) + cput - cpuold(ic)
          vputim(ic) = vputim(ic) + vput - vpuold(ic)
          return
       end if
    end do

  end subroutine clcend

! **********************************************************************

  subroutine yclock( cput2, vput2 )

    implicit none

    real(8),                intent(inout)  ::  cput2,   vput2

    real(8) ::  cput1

    call clock( cput1 )

    cput2 = cput1 - cput0
    vput2 = 0.d0

  end subroutine yclock

! **********************************************************************

  subroutine yclocl

    implicit none

    call clock( cput0 )

  end subroutine yclocl

end module ucloc


