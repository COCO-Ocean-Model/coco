module ufile
! --- information -----------------------------------------------------
!
!  Open and close files.
!
!  HISTORY
!     '01.12.06  H.Hasumi
!     '07.04.23  H.Hasumi
!     '12.10.06  M.kurogi: for COCO5.0
!
! ---------------------------------------------------------------------

 use zocfil, only : nfmax
 implicit none
 private
 public :: filopn, filcls, rewnml, cstnml
 logical, save   ::  opn(nfmax)
 data opn / nfmax*.false. /

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

 if (cact .eq. 'read') then
    inquire(file=cf, exist=oex)
    if (.not. oex) then
       i = index(cf, ' ') - 1
       write(nfstdo, *) '### file "', cf(1:i), '" does not exist ###'
       call mpi_abort(mpi_comm_world, 1, ierr)
    end if
 else if (cf(1:13) .eq. 'not-specified') then
    write(nfstdo, *) '### name not specified for output file ###'
    call mpi_abort(mpi_comm_world, 1, ierr)
 end if

 do i = nfmin, nfmax
    if (.not. opn(i)) then
       nf = i
       opn(nf) = .true.
       go to 123
    end if
 end do
 write(nfstdo, *) '### maximum file number exceeded ###'
 call mpi_abort(mpi_comm_world, 1, ierr)

 123  continue
 open(unit=nf, file=cf, form='unformatted', access='sequential')
 rewind(unit=nf)
 i = index(cf, ' ') - 1
 write(nfstdo, *) '*** file "', cf(1:i), '" opened for unit', nf, '***'

 return
 end subroutine filopn

! =====================================================================

 subroutine filcls(                                                            &
  &                    nf)
 use zocfil, only : nfstdo
 implicit none
 integer, intent(in) :: nf       
 close(unit=nf)
 opn(nf) = .false.
 write(nfstdo, *) '*** file unit', nf, 'closed ***'

 end subroutine filcls

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

 999  write(jfile, *) '*** error in rewinding the namelist file ***'
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
   &          '*** error occurs while reading namelist:', cnnml,               &
   &          ' at subroutine:', cnsbr, ' ***'
    stop
 end if

 return
 end subroutine cstnml
end module ufile
