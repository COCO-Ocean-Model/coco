program offline

! --- information -----------------------------------------------------
!
!  Off-line tracer model based on COCO4
!
!  HISTORY
!     '04.10.14  H.Hasumi: from COCO4.0
!     '06.01.30  H.Hasumi: change the timing of calling MPI_INIT
!     '07.04.24  H.Hasumi
!     '14.11.06  M.Watanabe: Tripolar + ecosys_o_fe
!
! ---------------------------------------------------------------------

   use zocdim, only: &
     &  nxdim,  nydim,  nzdim,  ntdim, &
     &   nwrk, &
     & myrank,   ierr
   use zocfil, only: &
     & nfstdo, nfomax,    ncf
   use zocout, only: &
     & loglev

   use aocea
   use atmct
   use brstt
   use ucloc
   use ufile
   use ucaln
   
   implicit none

#include "mpif.h"

   real(8) ::      t(nxdim, nydim, nzdim, ntdim)
   real(8) ::     ft(nxdim, nydim, ntdim)
   real(8) ::      u(nxdim, nydim, nzdim),      v(nxdim, nydim, nzdim)
   real(8) ::     ha(nxdim, nydim),            hb(nxdim, nydim)
   real(8) ::    ahv(nxdim, nydim, nzdim)

   real(8) ::  wrk(nwrk)

   real(8) ::      dt,     tt,     ts,    tss
   integer ::     nt,   itst,    its,   ntss
   real(8) ::   tstrt,   tend
   logical :: oflout(nfomax), oflstk(nfomax), orsout, orsrwd
   integer ::    ijk
   integer ::  ifpar,  jfpar,  istat
   integer :: lenstd
   integer :: idate(6)

   character(len=ncf) :: crun = '(RUN NAME WAS NOT SET)'
   character(len=ncf) :: cstdo = 'STDOUT'

   namelist /nmrun/ crun
   namelist /nmstdo/ cstdo
   namelist /nmlog/ loglev

! *** initial setup ***
   call mpi_init(ierr)
   call clcstr('SETUP')
   call rewnml(ifpar, jfpar)
   open(unit=ifpar, file='PARAMET', status='old', &
     &  access='sequential', form='formatted')
   call parset
   call rewnml(ifpar, jfpar)
   read (ifpar, nmstdo, iostat=istat)
   call cstnml(jfpar, __FILE__, __LINE__ -1, istat)

   lenstd = index(cstdo, ' ')
   write(cstdo(lenstd:lenstd+5), '(a1,i5.5)') '.', myrank
   call rewnml(ifpar, jfpar)
   open(unit=jfpar, file=cstdo, &
     &  access='sequential', form='formatted')
   write(jfpar, *) 'MESSAGE OUTPUT FOR RANK', myrank
   call rewnml(ifpar, jfpar)
   read (ifpar, nmrun, iostat=istat)
   call cstnml(jfpar, __FILE__, __LINE__ -1, istat)
   write(jfpar, *) 'Run name :'//crun
   call rewnml(ifpar, jfpar)
   read (ifpar, nmlog, iostat=istat)
   call cstnml(jfpar, __FILE__, __LINE__ -1, istat)
   write(nfstdo, '(a,i2)') ' log level : ', loglev
   do ijk = 1, nwrk
     wrk(ijk) = 0.d0
   end do

   call tmstup( tstrt,  tend,    dt )
   call restrt( tstrt,     t,    ft ) ! passive tracers' initial condition
   tt = tstrt
   nt = 0
   call ocstup(     t,     u,     v,    hb,   ahv,  &
   &               tt,    dt )
   call clcend('SETUP')

! *** main loop ***

   call clcstr('MAIN')
   do ijk = 1, nwrk
      wrk(ijk) = 0.d0
   end do
   do
      tt = tt + dt
      nt = nt + 1
      if (loglev > 1) then
         call css2yh(idate, tt)
         call rewnml(ifpar,jfpar)
         write(jfpar,'(A,I6,5(1X,I2.2))') ' *** time ***:', idate
         call flush(jfpar)
      end if
      call tmstpc( &
         &          itst,     ts,    its,         &
         &          ntss,    tss,                 &
         &        oflout, oflstk, orsout, orsrwd, &
         &            nt,     tt)
      call ocean ( &
         &             t,     ft,                         & !! (out)
         &             u,      v,     ha,     hb,    ahv, & !! (out)
         &            nt,     tt,   itst,     ts,    its, & !! (in )
         &          ntss,    tss,                         & !! (in )
         &        oflout, oflstk)                           !! (in )
      if (tt >= tend) exit
      call finout( &
         &            tt,     nt,         &
         &             t,     ft,     hb, &
         &        orsout, orsrwd)
   end do
   call clcend('MAIN')

! *** handling for termination ***

   call clcstr('FINOUT')
   orsout = .true.
   call finout( &
         &            tt,     nt,         &
         &             t,     ft,     hb, &
         &        orsout, orsrwd)
   call clcend('FINOUT')
   call clcout
   call parfin

end program offline
! *********************************************************************

subroutine parset

! --- information -----------------------------------------------------
!
!  Initialization for inter-node communication using MPI
!
!  HISTORY
!     '00.05.30  H.Hasumi: from COCO2
!
! ---------------------------------------------------------------------
  use zocdim, only: &
    & nxg, nyg, inodes, jnodes, &
#ifdef OPT_TRIPOLE
    &   jupe,   jupw,  jupfy, jdownfy, &
#endif
    &  nprocs, myrank, ijnode,  iroot,   ierr, &
    &  irank,    iup,  idown, &
    &  jrank,    jup,  jdown, &
    &  igrank, mpi_comm_ogcm

  use bgs2d
  use bgs3d
  use ufile

  implicit none

#include "mpif.h"

  integer ::  ifpar,  jfpar,  istat
#ifdef OPT_TRIPOLE
  integer ::  ii
#endif

  integer :: ndroot = 0
  namelist /nmroot/ ndroot
  integer :: key

!!call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world, nprocs, ierr)
  call mpi_comm_rank(mpi_comm_world, myrank, ierr)
  call mpi_comm_split(mpi_comm_world, 1, igrank, mpi_comm_ogcm, ierr)
  call mpi_comm_rank(mpi_comm_ogcm, myrank, ierr)
  
  call rewnml(ifpar, jfpar)
  read (ifpar, nmroot, iostat=istat)
  call cstnml(jfpar, 'parset', 'nmroot', istat)
!  write(jfpar, nmroot)

  ijnode = inodes * jnodes
  iroot  = ndroot
  
  if (ijnode > nprocs) then
     write(jfpar, *) &
     &   '### NODE NUMBER ERROR                   ###'
     write(jfpar, *) &
     &   '### INODES*JNODES must be = or < NPROCS ###'
     write(jfpar, *) 'INODES:', inodes
     write(jfpar, *) 'JNODES:', jnodes
     write(jfpar, *) 'NPROCS:', nprocs
     call mpi_finalize(ierr)
     stop
  end if

  if ( (mod(nxg, inodes) .ne. 0) .or. (mod(nyg, jnodes) .ne. 0) ) then
     write(*, *) ' #### NODE NUMBER ERROR #### '
     write(*, *) ' The size of sub-regions must be identical. '
     call mpi_finalize(ierr)
     stop
  end if

#ifdef OPT_TRIPOLE
  if ( (mod(inodes, 2) .ne. 0) .and. (inodes .ne. 1) ) then
     write(*, *) ' #### NODE NUMBER ERROR #### '
     write(*,*)'#### For tripolar grid, inodes shold be 1 or even.'
     call mpi_finalize(ierr)
     stop
  end if
#endif

  if ((ndroot < 0) .or. (ndroot > nprocs-1)) then
     write(jfpar, *) &
     &   '### ERROR in the choice of NDROOT          ###'
     write(jfpar, *) 'NDROOT:', ndroot
     write(jfpar, *) 'NPROCS:', nprocs
     write(jfpar, *) &
     &   '### NDROOT must lie between 0 and NPROCS-1 ###'
     call mpi_finalize(ierr)
     stop
  end if

  irank = mod(myrank, inodes)
  jrank = myrank / inodes

  if (irank == 0) then
     idown = inodes * (jrank + 1) - 1
  else
     idown = inodes * jrank + irank - 1
  end if
  if (irank == inodes-1) then
     iup = inodes * jrank
  else
     iup = inodes * jrank + irank + 1
  end if

  if (jrank == 0) then
     jdown = mpi_proc_null
  else
     jdown = inodes * (jrank - 1) + irank
  end if

#ifdef OPT_TRIPOLE
  if (jrank == jnodes-1) then
     if (inodes >= 2) then
        ii = mod(inodes,2)
        if (ii /= 0) then
           write(jfpar, *) 'INODES:', inodes
           write(jfpar, *) &
             &     '### INODES must be 1 or multiples of 2 ###'
           call mpi_finalize(ierr)
           stop
        end if

     ii = inodes / 2
     if (irank < ii) then
        jupe = inodes * jnodes - 1 - irank
        jupw = mpi_proc_null
     else
        jupe = mpi_proc_null
        jupw = inodes * jnodes  -1 - irank
     end if
  else
     jupe = myrank
     jupw = myrank
  endif
  else
     jupe = mpi_proc_null
     jupw = mpi_proc_null
  end if

  if (jrank .eq. 0) then
     jdownfy = inodes * (jnodes - 1) + irank
  else
     jdownfy = mpi_proc_null
  end if
  if (jrank .eq. jnodes-1) then
     jupfy = irank
  else
     jupfy = mpi_proc_null
  end if

#endif
  if (jrank == jnodes-1) then
     jup = mpi_proc_null
  else
     jup = inodes * (jrank + 1) + irank
  end if

  if (myrank >= inodes*jnodes) then
     iup = mpi_proc_null
     idown = mpi_proc_null
     jup = mpi_proc_null
     jdown = mpi_proc_null
  end if

  call gs2dst
#ifndef OPT_IO_COCOMPI
  call gs3dst
#endif
!!call gsidst
  return

end subroutine parset

! *********************************************************************

subroutine parfin

! --- information -----------------------------------------------------
!
!  Finalization for parallelization
!
!  HISTORY
!     '00.05.30  H.Hasumi: from COCO3
!
! ---------------------------------------------------------------------
  use zocnod, only: &
    &   ierr

  implicit none

#include "mpif.h"

  call mpi_finalize(ierr)

  return
end subroutine parfin

! *********************************************************************

subroutine cstnml( &
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
