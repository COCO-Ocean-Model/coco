module atmct
! --- information -----------------------------------------------------
!
!  Controling the scheme of time integration and the timing of output
!
!  HISTORY
!     '00.05.30  H.Hasumi: parallelized COCO3
!     '00.12.07  H.Hasumi: combine par and non-par routines
!     '01.12.07  H.Hasumi
!     '03.06.04  H.Hasumi: allow to specify the time step by the unit
!                          other than 'hour'
!     '07.04.23  H.Hasumi
!     '10.04.14  M.Kurogi: staggered time stepping
!     '12.07.18  M.kurogi: for COCO5.0
! ---------------------------------------------------------------------

 use zocfil, only : nfomax, ncf
 implicit none

 private
 public :: tmstup, tmstpc
 real(8), save :: dt
 integer, save :: ionext(6, nfomax)
 integer, save :: irnext(6)
 integer, save :: iflout(nfomax)=0
 real(8), save :: tostrt(nfomax),  toend(nfomax), tonext(nfomax)
 real(8), save :: trnext
 integer, save :: irintv=1, iurint=1, irsrwd=1
 integer, save :: iointv(nfomax), iuintv(nfomax)
 integer, save :: ntsplt=0
 integer, save :: nohitm


contains
  subroutine tmstup(                                                           &
   &             tstrt,   tend,    dtt)
  use zocnod, only : ierr
  use ufile
  use ucaln

  implicit none
#include "mpif.h"

  real(8), intent(out) ::   tstrt,   tend,    dtt
  integer ::  i, iitem
  character (6) :: cunit(6) =                                                  &
   & (/ 'year', 'month', 'day', 'hour', 'minute', 'second' /)
  integer   :: iohitm
  character (16) :: citem(nfomax)

  integer :: itstrt(6)=(/ 0, 0, 0, 0, 0, 0 /)
  integer ::  itend(6)=(/ 0, 0, 0, 0, 0, 0 /)
  integer ::  icaln=0
  integer, save :: iostrt(6, nfomax), ioend(6, nfomax)
  integer :: iutstp=0
  integer :: idate(6) =(/ 0, 1, 1, 0, 0, 0 /)
  real(8) ::  tmstp=1.d10

  integer :: iodstr(6)=(/ 0, 0, 0, 0, 0, 0 /)
  integer :: iodend(6)=(/ 0, 0, 0, 0, 0, 0 /)
  integer :: iodint=1, iudint=1, iodavr=1, iodsng=1

  integer :: iohstr(6), iohend(6), iohint, iuhint, iohavr, iohsng
  integer :: ioxstr, ioxend, ioystr, ioyend, iozstr, iozend
  character (ncf) :: cohfil
  character (16) :: cohitm='   not-specified'
  integer :: ifpar,  jfpar
  integer :: istat 

  namelist /nmtime/ itstrt, itend, tmstp, iutstp, ntsplt
  namelist /nmcaln/ icaln
  namelist /nmrstr/ irintv, iurint, irsrwd
  namelist /nmdout/ iodstr, iodend, iodint, iudint, iodavr, iodsng
  namelist /nmhist/ cohitm, cohfil,                                            &
   &                  iohstr, iohend, iohint, iuhint, iohavr, iohsng,          &
   &                  ioxstr, ioxend, ioystr, ioyend, iozstr, iozend


  call rewnml(ifpar, jfpar)
  read(ifpar, nmtime, iostat=istat)
  call cstnml(jfpar, 'tmstup', 'nmtime', istat)

  call rewnml(ifpar, jfpar)
  read(ifpar, nmcaln, iostat=istat)
  call cstnml(jfpar, 'tmstup', 'nmcaln', istat)

  call rewnml(ifpar, jfpar)
  read(ifpar, nmrstr, iostat=istat)
  call cstnml(jfpar, 'tmstup', 'nmrstr', istat)

  call rewnml(ifpar, jfpar)
  read(ifpar, nmdout, iostat=istat)
  call cstnml(jfpar, 'tmstup', 'nmdout', istat)

  do iitem = 1, nfomax
     do i = 1, 6
        iostrt(i, iitem) = iodstr(i)
        ioend (i, iitem) = iodend(i)
     end do
        iointv(iitem) = iodint
        iuintv(iitem) = iudint
  end do

  call rewnml(ifpar, jfpar)
  nohitm = 0

  do 
      cohitm = 'not-specified'
      iohstr(1) = -1
      iohend(1) = -1
      iohint = -1
      iuhint = -1
      read(ifpar, nmhist, iostat=istat)
      if(istat < 0) exit

      if (cohitm(1:13) /= 'not-specified') then
         nohitm = nohitm + 1
         iohitm = nohitm
         iflout(iohitm) = 1
         citem(iohitm) = cohitm
         if (iohstr(1) >= 0) then
            do i = 1, 6
               iostrt(i, iohitm) = iohstr(i)
            end do
         end if
         if (iohend(1) >= 0) then
            do i = 1, 6
               ioend(i, iohitm) = iohend(i)
            end do
         end if
         if (iohint > 0) then
            iointv(iohitm) = iohint
         end if
         if (iuhint > 0) then
            iuintv(iohitm) = iuhint
         end if
      end if
  end do




!     ***** reading time control parameters *****

  call calndr(                                                                 &
   &             icaln)
  call cyh2ss(                                                                 &
   &             tstrt,                                                        &
   &            itstrt)
  call cyh2ss(                                                                 &
   &              tend,                                                        &
   &             itend)

  do i = 1, 6
     irnext(i) = itstrt(i)
  end do
  do iitem = 1, nohitm
     call cyh2ss(                                                              &
   &            tostrt(iitem),                                                 &
   &            iostrt(1, iitem))
     call cyh2ss(                                                              &
   &             toend(iitem),                                                 &
   &             ioend(1, iitem))
     if (tostrt(iitem) < tstrt) then
         tostrt(iitem) = tstrt
        call css2yh(                                                           &
   &               iostrt(1, iitem),                                           &
   &               tostrt(iitem))
     end if
     do i = 1, 6
        ionext(i, iitem) = iostrt(i, iitem)
     end do
  end do

  if (iutstp == 0) then
     dt = anint(tmstp * 3.6d+3)
     dtt = dt
  else if ((iutstp >= 1) .and. (iutstp <= 6)) then
     idate(iutstp) = idate(iutstp) + nint(tmstp)
     call cyh2ss(                                                              &
   &                   dt,                                                     &
   &                idate)
     dtt = dt
  else
     write(jfpar, *) '*** tmstup: no such unit of time ***'

     call mpi_abort(mpi_comm_world, 1, ierr)
  end if

  do iitem = 1, nohitm
     if ((iuintv(iitem) >= 1) .and. (iuintv(iitem) <= 6)) then
        ionext(iuintv(iitem), iitem)                                           &
   &  = ionext(iuintv(iitem), iitem) + iointv(iitem)
        call cyh2ss(                                                           &
   &             tonext(iitem),                                                &
   &             ionext(1, iitem))
        call css2yh(                                                           &
   &             ionext(1, iitem),                                             &
   &             tonext(iitem))
     else if (iflout(iitem) == 1) then
        write(jfpar, *) '*** tmstup: no such unit of time ***'
        call mpi_abort(mpi_comm_world, 1, ierr)
     end if
  end do

  if ((iurint >= 1) .and. (iurint <= 6)) then
     irnext(iurint) = irnext(iurint) + irintv
     call cyh2ss(                                                              &
   &             trnext,                                                       &
   &             irnext)
     call css2yh(                                                              &
   &             irnext,                                                       &
   &             trnext)
  else
     write(jfpar, *) '*** tmstup: no such unit of time ***'
     call mpi_abort(mpi_comm_world, 1, ierr)
  end if

  write(jfpar, *) '***** time control parameters *****'
  write(jfpar, *) '***** time staggering  version ****'
  write(jfpar, *) 
  write(jfpar, *) '     start time  :', itstrt
  write(jfpar, *) '       end time  :', itend
  if (iutstp == 0) iutstp = 4
  write(jfpar, *) '      time step  :', tmstp, ' ', cunit(iutstp)
  write(jfpar, *) '     time split  :', ntsplt
  write(jfpar, *)
  write(jfpar, *) 'restart interval :', irintv, ' ', cunit(iurint)
  write(jfpar, *)

  do iitem = 1, nohitm
     if (iflout(iitem) == 1) then
        write(jfpar, *)
        write(jfpar, *) 'output item      :', citem(iitem)
        write(jfpar, *) 'output starts at :',                                  &
   &                     (iostrt(i, iitem), i = 1, 6)
        write(jfpar, *) 'output  ends  at :',                                  &
   &                     (ioend(i, iitem), i = 1, 6)        
        write(jfpar, *) 'output interval  :', iointv(iitem),                   &
   &                     ' ', cunit(iuintv(iitem))
     end if
  end do

  if( mod(ntsplt,2) /= 0) then
     write(jfpar, *) "error: ntsplt mut be even number !"
     call mpi_abort(mpi_comm_world, 1, ierr)
  end if

  return
  end subroutine tmstup
! =====================================================================

  subroutine tmstpc (                                                          &
   &                itst,     ts,    its,                                      &
   &                ntss,    tss,                                              &
   &              oflout, oflstk, orsout, orsrwd,                              &
   &                  nt,     tt)

  use ucaln

  implicit none

  integer, intent(out) ::    its,   itst,   ntss
  real(8), intent(out) ::    ts,    tss
  logical, intent(out) :: oflout(nfomax), oflstk(nfomax)
  logical, intent(out) :: orsout, orsrwd
  integer, intent (in) :: nt
  real(8), intent (in) :: tt
  integer :: iitem

  itst = 2
  ts  = dt
  its = 1
  ntss = ntsplt 
  tss  = ts / dble(ntss)

  do iitem = 1, nohitm
     if (iflout(iitem) == 1) then
        if (      (tt > tostrt(iitem))                                         &
   &        .and. (tt < toend(iitem)+dt)) then
           oflstk(iitem) = .true.
           if (tt >= tonext(iitem)) then
              oflout(iitem) = .true.
              ionext(iuintv(iitem), iitem)                                     &
   &        = ionext(iuintv(iitem), iitem) + iointv(iitem)
              call cyh2ss(                                                     &
   &                tonext(iitem),                                             &
   &                ionext(1, iitem))
              call css2yh(                                                     &
   &                ionext(1, iitem),                                          &
   &                tonext(iitem))
           else
              oflout(iitem) = .false.
           end if
        else
           oflstk(iitem) = .false.
           oflout(iitem) = .false.
        end if
     else
        oflstk(iitem) = .false.
        oflout(iitem) = .false.
     end if
  end do

  if (tt >= trnext) then
     orsout = .true.
     if (irsrwd == 1) then
        orsrwd = .true.
     else
        orsrwd = .false.
     end if
     irnext(iurint) = irnext(iurint) + irintv
     call cyh2ss(                                                              &
   &             trnext,                                                       &
   &             irnext)
     call css2yh(                                                              &
   &             irnext,                                                       &
   &             trnext)
  else
     orsout = .false.
  end if

  return
  end subroutine tmstpc
end module atmct
