program coco

! --- information -----------------------------------------------------
!
!  COCO4: free surface, z-sigma hybrid vertical coordinate, general
!  curvilinear horizonal coordinate, primitive equation OGCM
!
!  HISTORY
!     '03.04.21  H.Hasumi: from COCO3.4
!     '06.01.30  H.Hasumi: change the timing of calling MPI_INIT
!     '06.04.23  H.Hasumi
!     '08.07.10  H.Hasumi: initial/final processing
!     '09.09.02  Y.Komuro: STDOUT filename format changed
!     '10.04.14  M.Kurogi: (COCO4.4 tripolar code by Dr. Suzuki)
!     '12.12.06  Y.Komuro: for COCO5.0
!
! ---------------------------------------------------------------------

  use zocdim, only: &
    &  nxdim,  nydim,  nzdim,  ntdim, &
#ifdef OPT_PARALLEL
    & myrank,   ierr, &
#endif
    & ofinal,  oinit
  use zocfil, only: &
    & nfstdo, nfomax,    ncf

  use aocea
  use atmct
  use bfrch
  use brstt
  use ucloc
  use ufile

  implicit none

#ifdef OPT_PARALLEL
#include "mpif.h"
#endif

  real(8) ::     ta(nxdim, nydim, nzdim, ntdim)
  real(8) ::     tb(nxdim, nydim, nzdim, ntdim)
  real(8) ::     ua(nxdim, nydim, nzdim),     ub(nxdim, nydim, nzdim)
  real(8) ::     va(nxdim, nydim, nzdim),     vb(nxdim, nydim, nzdim)
  real(8) ::     ha(nxdim, nydim),            hb(nxdim, nydim)
  real(8) ::   ubta(nxdim, nydim),          ubtb(nxdim, nydim)
  real(8) ::   vbta(nxdim, nydim),          vbtb(nxdim, nydim)
  real(8) ::      w(nxdim, nydim, nzdim),      r(nxdim, nydim, nzdim)
  real(8) ::    amv(nxdim, nydim, nzdim),    ahv(nxdim, nydim, nzdim)
  real(8) ::     ft(nxdim, nydim, ntdim)            

  real(8) ::      dt,     tt,     ts,    tss
  integer ::      nt,   itst,    its,   ntss
  real(8) ::   tstrt,   tend
  logical ::  oflout(nfomax), oflstk(nfomax), orsout, orsrwd
  integer ::     ijk
  integer ::   ifpar,  jfpar,  istat
#ifdef OPT_PARALLEL
  integer ::  lenstd
#endif

  character(len=ncf) :: crun = '(RUN NAME WAS NOT SET)'
  character(len=ncf) :: cstdo = 'STDOUT'

  namelist /nmrun/ crun
  namelist /nmstdo/ cstdo

! *** Initial setup ***

#ifdef OPT_PARALLEL
  call mpi_init(ierr)
#endif
  call clcstr('SETUP')
  call rewnml(ifpar, jfpar)
  open(unit=ifpar, file='PARAMET', status='old', &
    &  access='sequential', form='formatted')
#ifdef OPT_PARALLEL
  call parset
#endif
  call rewnml(ifpar, jfpar)
  read (ifpar, nmstdo, iostat=istat)
  call cstnml(jfpar, 'coco', 'nmstdo', istat)
!  write(jfpar, nmstdo)
#ifdef OPT_PARALLEL
  lenstd = index(cstdo, ' ')
  write(cstdo(lenstd:lenstd+3), '(a1,i3.3)') '.', myrank
#endif
  call rewnml(ifpar, jfpar)
  open(unit=jfpar, file=cstdo, &
    &  access='sequential', form='formatted')
#ifdef OPT_PARALLEL
  write(jfpar, *) 'MESSAGE OUTPUT FOR RANK', myrank
#endif
  call rewnml(ifpar, jfpar)
  read (ifpar, nmrun, iostat=istat)
  call cstnml(jfpar, 'coco', 'nmrun', istat)
!  write(jfpar, nmrun)
  write(nfstdo, *) 'Run name :'//crun

  oinit = .false.
  ofinal = .false.
  call tmstup( &
    &           tstrt,   tend,     dt )
  call restrt( &
    &           tstrt, &
    &              ub,     vb,     tb, &
    &              hb,   ubtb,   vbtb, &
    &               w,    amv,    ahv, &
    &              ft )
  call ocstup( &
    &              ub,     vb,     tb, &
    &              hb,   ubtb,   vbtb,      w,      r, &
    &             amv,    ahv, &
    &              ft, &
    &              dt )
  oinit = .true.
  call ocean ( &
    &              ua,     va,     ta, &
    &              ha,   ubta,   vbta, &
    &              ub,     vb,     tb, &
    &              hb,   ubtb,   vbtb, &
    &               w,      r,    amv,    ahv, &
    &              ft, &
    &              nt,     tt,   itst,     ts,    its, &
    &            ntss,    tss, &
    &          oflout, oflstk )
  oinit = .false.
  call forsto( &
    &              ua,     va,     ta, &
    &              ha,   ubta,   vbta, &
    &              ub,     vb,     tb, &
    &              hb,   ubtb,   vbtb )
  call clcend('SETUP')

! *** Main loop ***

  call clcstr('MAIN')
  tt = tstrt
  nt = 0
  do
     tt = tt + dt
     nt = nt + 1
     call tmstpc( &
       &            itst,     ts,    its, &
       &            ntss,    tss, &
       &          oflout, oflstk, orsout, orsrwd, &
       &              nt,     tt )
     call ocean ( &
       &              ua,     va,     ta, &
       &              ha,   ubta,   vbta, &
       &              ub,     vb,     tb, &
       &              hb,   ubtb,   vbtb, &
       &               w,      r,    amv,    ahv, &
       &              ft, &
       &              nt,     tt,   itst,     ts,    its, &
       &            ntss,    tss, &
       &          oflout, oflstk )
     if (tt >= tend) then
        exit
     end if
     if (itst == 3) then
        call finout( &
          &              tt,     nt, &
          &              ua,     va,     ta, &
          &              ha,   ubta,   vbta, &
          &               w,    amv,    ahv, &
          &              ft, &
          &          orsout, orsrwd )
     else
        call finout( &
          &              tt,     nt, &
          &              ub,     vb,     tb, &
          &              hb,   ubtb,   vbtb, &
          &               w,    amv,    ahv, &
          &              ft, &
          &          orsout, orsrwd )
     end if
     if (orsout) then
        ofinal = .true.
        call ocean ( &
          &           ua,     va,     ta, &
          &           ha,   ubta,   vbta, &
          &           ub,     vb,     tb, &
          &           hb,   ubtb,   vbtb, &
          &            w,      r,    amv,    ahv, &
          &           ft, &
          &           nt,     tt,   itst,     ts,    its, &
          &         ntss,    tss, &
          &       oflout, oflstk )
        ofinal = .false.
     end if
  end do
  call clcend('MAIN')

! *** Handling for termination ***

  call clcstr('FINOUT')
  ofinal = .true.
  orsout = .true.
  if (itst == 3) then
     call finout( &
       &              tt,     nt, &
       &              ua,     va,     ta, &
       &              ha,   ubta,   vbta, &
       &               w,    amv,    ahv, &
       &              ft, &
       &          orsout, orsrwd )
  else
     call finout( &
       &              tt,     nt, &
       &              ub,     vb,     tb, &
       &              hb,   ubtb,   vbtb, &
       &               w,    amv,    ahv, &
       &              ft, &
       &          orsout, orsrwd )
  end if
  call ocean ( &
    &                 ua,     va,     ta, &
    &                 ha,   ubta,   vbta, &
    &                 ub,     vb,     tb, &
    &                 hb,   ubtb,   vbtb, &
    &                  w,      r,    amv,    ahv, &
    &                 ft, &
    &                 nt,     tt,   itst,     ts,    its, &
    &               ntss,    tss, &
    &             oflout, oflstk )
  call clcend('FINOUT')
  call clcout
#ifdef OPT_PARALLEL
  call parfin
#endif

end program coco

#ifdef OPT_PARALLEL
! *********************************************************************

subroutine parset

! --- information -----------------------------------------------------
!
!  Initialization for inter-node communication using MPI
!
!  HISTORY
!     '00.05.30  H.Hasumi: from COCO2
!     '12.12.06  Y.Komuro: for COCO5.0
!
! ---------------------------------------------------------------------

  use zocdim, only: &
    & inodes, jnodes, &
#ifdef OPT_TRIPOLE
    &   jupe,   jupw, &
#endif
    & nprocs, myrank, ijnode,  iroot,   ierr, &
    &  irank,    iup,  idown, &
    &  jrank,    jup,  jdown

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

!!!  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world, nprocs, ierr)
  call mpi_comm_rank(mpi_comm_world, myrank, ierr)

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
             &      '### INODES must be 1 or multiples of 2 ###'
           call mpi_finalize(ierr)
           stop
        endif
         
        ii = inodes / 2
        if (irank < ii) then
           jupe = inodes * jnodes - 1 - irank
           jupw = mpi_proc_null
        else
           jupe = mpi_proc_null
           jupw = inodes * jnodes  -1 - irank  
        endif
     else
        jupe = myrank
        jupw = myrank
     endif
  else
     jupe = mpi_proc_null
     jupw = mpi_proc_null
  endif
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
  call gs3dst

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
!     '12.12.06  Y.Komuro: for COCO5.0
!
! ---------------------------------------------------------------------

  use zocnod, only: &
    &   ierr   

  implicit none

#include "mpif.h"

  call mpi_finalize(ierr)

  return
end subroutine parfin
#endif

