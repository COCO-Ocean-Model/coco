module ucloc

! --- information -----------------------------------------------------
!
!  Monitoring CPU time on Sun OS (Sun FORTRAN, using the intrinsic
! function "ETIME")
!
!  HISTORY
!     '97.03.18  H.Hasumi: From AGCM5.4 developed by A.Numaguti
!     '07.04.23  H.Hasumi
!     '12.06.15  H.Tatebe: for COCO5.0 in F90
!
! ---------------------------------------------------------------------

  implicit none
  private
#include "mpif.h"
  integer(4),    parameter          ::  nclmax = 100
  real(8),                    save  ::  cputim(nclmax), vputim(nclmax)
  real(8),                    save  ::  cpuold(nclmax), vpuold(nclmax)
  real(8),                    save  ::  cput,   vput        
  integer(4),                 save  ::  ifpar,  jfpar
  integer(4),                 save  ::  nclock
  logical,                    save  ::  ofirst
  character(16),              save  ::  htitle(nclmax)
                              
!---- used in yclock
  real(8),                    save  ::  cput2,  vput2
  real(8),                    save  ::  ticks,  tick0,  tusr0
  real(4),                    save  ::  tarray(1:2)

  public  ::  clcout,  clcstr,  clcend  !  ued in icedcoco.F, aprdc.F, iprdc,F
  
  data tick0, tusr0 / 0.d0, 0.d0 /
  data nclock / 0 /
  data ofirst / .true. /
  data htitle / nclmax*' ' /

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
!51  format(' ', a16, 2f15.6)
    
    call yclock(cput, vput)
    write(jfpar, '(1x,a16,2f15.6)') ' total time = ', cput, vput
    
  end subroutine clcout
  
! **********************************************************************

  subroutine clcstr( httl )
    use zocnod, only: mpi_comm_ogcm
    implicit none
    character(*),           intent(in)     ::  httl

    integer(4)  ::  ic
    integer :: ierr
    call mpi_barrier(mpi_comm_ogcm, ierr)

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

    real(4)                                ::  etime

    ticks = etime( tarray )
    cput2 = mpi_wtime()-tick0
    vput2 = mpi_wtime()-tusr0

  end subroutine yclock

! **********************************************************************

  subroutine yclocl

    implicit none

    real(4)                                ::  etime

    tick0 = mpi_wtime()
    tusr0 = mpi_wtime()

  end subroutine yclocl

end module ucloc


