module ufile
! --- information -----------------------------------------------------
!
!  Open and close files.
!
!  HISTORY
!     '01.12.06  H.Hasumi
!     '07.04.23  H.Hasumi
!     '12.10.06  M.kurogi: for COCO5.0
!     '15.04.08  M.Kurogi: for MPI-IO
!
! ---------------------------------------------------------------------

 use zocfil, only : nfmax
 use zocnod, only : mpi_comm_ogcm
 implicit none
 private
 public ::  &
#ifdef OPT_IO_COCOMPI
 & mpi_filopn, mpi_filcls, mpi_filcls_all, &
#endif
 & filopn, filcls, rewnml, cstnml, stop_msg


 logical, save   ::  opn(nfmax)
 data opn / nfmax*.false. /

 integer, parameter :: mpfmax=1000
 logical, save :: mpiopn(mpfmax)=.false.
 integer, save :: mpifh(mpfmax) 
 integer, save :: mpf=0
contains
! =====================================================================
 subroutine filopn(                                                            &
  &                      nf,                                                   &
  &                      cf,   cact)
 use zocfil
 implicit none
#include "mpif.h"
 character,  intent(in ) ::    cf*(ncf),   cact*(*)
 integer,    intent(out) ::    nf
 logical   ::   oex
 integer   ::     i,   ierr

 if (cact .eq. 'READ') then
    inquire(file=cf, exist=oex)
    if (.not. oex) then
       i = index(cf, ' ') - 1
       write(nfstdo, *) '### FILE "', cf(1:i), '" DOES NOT EXIST ###'
       call mpi_abort(mpi_comm_ogcm, 1, ierr)
    end if
 else if (cf(1:13) .eq. 'not-specified') then
    write(nfstdo, *) '### NAME NOT SPECIFIED FOR OUTPUT FILE ###'
    call mpi_abort(mpi_comm_ogcm, 1, ierr)
 end if

 do i = nfmin, nfmax
    if (.not. opn(i)) then
       nf = i
       opn(nf) = .true.
       go to 123
    end if
 end do
 write(nfstdo, *) '### MAXIMUM FILE NUMBER EXCEEDED ###'
 call mpi_abort(mpi_comm_ogcm, 1, ierr)

 123  continue
 open(unit=nf, file=cf, form='unformatted', access='sequential')
 rewind(unit=nf)
 i = index(cf, ' ') - 1
 write(nfstdo, *) '*** FILE "', cf(1:i), '" OPENED FOR UNIT', nf, '***'

 return
 end subroutine filopn

#ifdef OPT_IO_COCOMPI
 subroutine mpi_filopn(                                                            &
  &                      nf,                                                   &
  &                      cf,   cact)
 use zocfil
 implicit none
#include "mpif.h"
 character,  intent(in ) ::    cf*(ncf),   cact*(*)
 integer,    intent(out) ::    nf
 logical   ::   oex
 integer   ::     i,   ierr

 if (cact .eq. 'READ') then
    inquire(file=cf, exist=oex)
    if (.not. oex) then
       i = index(cf, ' ') - 1
       write(nfstdo, *) '### FILE "', cf(1:i), '" DOES NOT EXIST ###'
       call mpi_abort(mpi_comm_ogcm, 1, ierr)
    end if
 else if (cf(1:13) .eq. 'not-specified') then
    write(nfstdo, *) '### NAME NOT SPECIFIED FOR OUTPUT FILE ###'
    call mpi_abort(mpi_comm_ogcm, 1, ierr)
 end if

 if (cact == 'READ') then
   call mpi_file_open(mpi_comm_ogcm, cf, mpi_mode_rdonly, &
  &   mpi_info_null, nf, ierr)
 else
!    inquire(file=cf, exist=oex)
!    if (oex) call mpi_file_delete(cf,mpi_info_null)
   call mpi_file_open(mpi_comm_ogcm, cf, mpi_mode_create + mpi_mode_wronly, &
  &   mpi_info_null, nf, ierr)
 end if

 i = index(cf, ' ') - 1
 write(nfstdo, *) '*** FILE "', cf(1:i), '" OPENED BY MPI FOR UNIT', nf, '***'

  mpf=mpf+1
  mpifh(mpf)=nf
  mpiopn(mpf)=.true.

 return
 end subroutine mpi_filopn
#endif
! =====================================================================

 subroutine filcls(                                                            &
  &                    nf)
 use zocfil, only : nfstdo
 implicit none
 integer, intent(in) :: nf       
 close(unit=nf)
 opn(nf) = .false.
 write(nfstdo, *) '*** FILE UNIT', nf, 'CLOSED ***'

 end subroutine filcls

#ifdef OPT_IO_COCOMPI
 subroutine mpi_filcls(                                                        &
  &                    nf)
 use zocfil, only : nfstdo
 implicit none
 integer, intent(in) :: nf       
 integer :: ierr, n

 write(nfstdo, *) '*** FILE UNIT', nf, 'CLOSED BY MPI***'
 do n=1,mpf
   if(mpifh(n) == nf ) mpiopn(n)=.false.
 end do

 call mpi_file_close(nf, ierr)


 end subroutine mpi_filcls


 subroutine mpi_filcls_all
 use zocfil, only : nfstdo
 implicit none
 integer :: ierr, n

 do n=1,mpf
    if(mpiopn(n)) then
       write(nfstdo, *) '*** FILE UNIT', mpifh(n), 'CLOSED BY MPI***'
       call mpi_file_close(mpifh(n), ierr)
    else
       write(nfstdo, *) '*** FILE UNIT', mpifh(n), 'already closed'
    end if
 end do

 end subroutine mpi_filcls_all
#endif
! *********************************************************************

 subroutine rewnml(                                                            &
 &                       ifile,  jfile)
 use zocfil, only : nfparm, nfstdo
 implicit none
 integer, intent(out) :: ifile,  jfile

 ifile = nfparm
 jfile = nfstdo
 rewind(ifile, err=999)
 return

 999  write(jfile, *) '*** ERROR IN REWINDING THE NAMELIST FILE ***'
 return
 end subroutine rewnml

! *********************************************************************

 subroutine cstnml(                                                            &
  &                        jfile,  cnsbr,  cnnml,  istat)
 implicit none
 integer,      intent(in) ::  jfile, istat 
 character(*), intent(in) ::  cnsbr,  cnnml

 if (istat > 0) then ! an error occured while reading
    write(jfile, *)                                                            &
   &          '*** ERROR OCCURS while reading namelist:', cnnml,               &
   &          ' at subroutine:', cnsbr, ' ***'
    stop
 end if

 return
 end subroutine cstnml

 subroutine stop_msg(cmsg, cfil, line)
   implicit none
   character(len=*), intent(in) :: cmsg, cfil
   integer,          intent(in) :: line
   integer :: ifpar, jfpar

   call rewnml(ifpar, jfpar)
   write(jfpar,'(a)')'===================================================='
   write(jfpar,'(a,i4,a)')'program stopped at line ', line, ' of '//cfil
   write(jfpar,'(a)')cmsg
   write(jfpar,'(a)')'===================================================='
   stop

 end subroutine stop_msg
end module ufile
