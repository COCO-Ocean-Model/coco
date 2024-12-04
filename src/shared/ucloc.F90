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

  integer(4),    parameter          ::  nclmax = 100
  real(8),                    save  ::  cputim(nclmax)=0.d0, wcltim(nclmax)=0.d0
  real(8),                    save  ::  cpuold(nclmax), wclold(nclmax)
  real(8),                    save  ::  cput, wclt
  integer(4),                 save  ::  ifpar,  jfpar
  integer(4),                 save  ::  nclock = 0
  logical,                    save  ::  ofirst = .true.
  character(16),              save  ::  htitle(nclmax)
                              
!---- used in yclock
  real(8),                    save  ::  cput2,  cput0
  real(8),                    save  ::  wclt2,  wclt0
  integer(8),                 save  ::    crt,    cmx

  public  ::  clcout,  clcstr,  clcend  !  ued in icedcoco.F, aprdc.F, iprdc,F
  
  data htitle / nclmax*' ' /

contains

  subroutine clcout
    use ufile
    implicit none
    integer(4)  ::  ic

    call rewnml(ifpar, jfpar)
    call yclock(cput, wclt)
    
    write(jfpar, '(17x,2a25)') 'CPU time (cpu_time)', 'real time (system_clock)'
    do ic = 1, nclock
       if ( htitle(ic) /= ' ' ) then
          write(jfpar, '(1x,a16,2(f15.6,"s (",f5.2,"%)"))') htitle(ic), &
               & cputim(ic), cputim(ic)/cput*100.d0, &
               & wcltim(ic), wcltim(ic)/wclt*100.d0
       end if
    end do
    
    write(jfpar, '(1x,a16,2(f15.6,"s",9x))') ' total time = ', cput, wclt
    
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

    call yclock( cput, wclt )
    do ic = 1, nclock
       if ( htitle(ic) == httl ) then
          cpuold(ic) = cput
          wclold(ic) = wclt
          return
       end if
    end do
    if ( nclock <= nclmax ) then
       nclock = nclock + 1
       htitle(nclock) = httl
       cpuold(nclock) = cput
       wclold(nclock) = wclt
    end if
    
  end subroutine clcstr

! **********************************************************************

  subroutine clcend( httl )
    implicit none
    character(*),           intent(in)     ::  httl
    integer(4)  ::  ic

    call yclock( cput, wclt )
    do ic = 1, nclock
       if (htitle(ic) == httl) then
          cputim(ic) = cputim(ic) + cput - cpuold(ic)
          wcltim(ic) = wcltim(ic) + wclt - wclold(ic)
          return
       end if
    end do

  end subroutine clcend

! **********************************************************************

  subroutine yclock( cput2, wclt2 )
    implicit none
    real(8),                intent(inout)  ::  cput2,  wclt2
    integer(8)                             ::    cnt

    call cpu_time(cput2)
    cput2 = cput2 - cput0
    
    call system_clock(cnt)
    wclt2 = dble(cnt) / dble(crt) - wclt0
    if (wclt2 < 0) then
       wclt2 = wclt2 + dble(cmx) / dble(crt)
    end if

  end subroutine yclock

! **********************************************************************

  subroutine yclocl
    implicit none
    integer(8)                             ::    cnt
    call cpu_time(cput0)
    call system_clock(cnt, crt, cmx)
    wclt0 = dble(cnt) / dble(crt)
    
  end subroutine yclocl

end module ucloc


