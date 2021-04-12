program icedcoco

! --- information -----------------------------------------------------
!
!  COCO3.4 coupled with sea ice
!
!  HISTORY
!     '02.10.10  H.Hasumi: from COCO3.4
!     '06.01.30  H.Hasumi: change the timing of calling MPI_INIT
!     '07.04.23  H.Hasumi
!     '07.09.25  H.Hasumi: multi-category sea ice
!     '07.09.27  H.Hasumi: 1-layer sea ice thermodynamics
!     '08.06.11  H.Hasumi: initial/final processing
!     '08.07.10  H.Hasumi: initial/final processing
!     '09.09.02  Y.Komuro: STDOUT filename format changed
!     '10.04.14  M.kurogi: (COCO4.4 tripolar code by Dr. Suzuki)
!     '12.10.16  Y.Komuro: for COCO5.0
!     '13.02.13  Y.Komuro: remove non-parallel code 
!
! ---------------------------------------------------------------------

  use zocdim, only: &
    &  nxdim,  nydim,  nzdim,  ntdim,  nic, &
    & myrank,   ierr, &
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

#include "mpif.h"

  real(8) ::     ta(nxdim, nydim, nzdim, ntdim)
  real(8) ::     tb(nxdim, nydim, nzdim, ntdim)
  real(8) ::     ua(nxdim, nydim, nzdim),     ub(nxdim, nydim, nzdim)
  real(8) ::     va(nxdim, nydim, nzdim),     vb(nxdim, nydim, nzdim)
  real(8) ::     ha(nxdim, nydim),            hb(nxdim, nydim)
  real(8) ::   ubta(nxdim, nydim),          ubtb(nxdim, nydim)
  real(8) ::   vbta(nxdim, nydim),          vbtb(nxdim, nydim)
  real(8) ::      w(nxdim, nydim, nzdim),      r(nxdim, nydim, nzdim)
  real(8) ::    amv(nxdim, nydim, nzdim),    ahv(nxdim, nydim, nzdim)

  real(8) ::     aa(nxdim, nydim, 0:nic),     ab(nxdim, nydim, 0:nic)
  real(8) ::    hia(nxdim, nydim, 0:nic),    hib(nxdim, nydim, 0:nic)
  real(8) ::    uia(nxdim, nydim),    uib(nxdim, nydim)
  real(8) ::    via(nxdim, nydim),    vib(nxdim, nydim)
  real(8) ::    tia(nxdim, nydim, 0:nic),    tib(nxdim, nydim, 0:nic)
  real(8) ::    hsa(nxdim, nydim, 0:nic),    hsb(nxdim, nydim, 0:nic)
  real(8) ::    tsi(nxdim, nydim, 0:nic)

  real(8) ::     ft(nxdim, nydim, ntdim)
  real(8) ::  swabs(nxdim, nydim),     fs(nxdim, nydim)
  real(8) ::   taux(nxdim, nydim),   tauy(nxdim, nydim)
  real(8) ::   ptop(nxdim, nydim)

  real(8) ::     dt,     tt,     ts,    tss
  integer ::     nt,   itst,    its,   ntss
  real(8) ::  tstrt,   tend
  logical :: oflout(nfomax), oflstk(nfomax), orsout, orsrwd
  integer ::    ijk
  integer ::  ifpar,  jfpar,  istat
  integer :: lenstd

  character(len=ncf) :: crun = '(RUN NAME WAS NOT SET)'
  character(len=ncf) :: cstdo = 'STDOUT'

  namelist /nmrun/ crun
  namelist /nmstdo/ cstdo

! *** Initial setup ***

  call mpi_init(ierr)
  call clcstr('SETUP')
  call rewnml(ifpar, jfpar)
  open(unit=ifpar, file='PARAMET', status='old', &
    &  access='sequential', form='formatted')
  call parset
  call rewnml(ifpar, jfpar)
  read (ifpar, nmstdo, iostat=istat)
  call cstnml(jfpar, 'icedcoco', 'nmstdo', istat)
!  write(jfpar, nmstdo)

  lenstd = index(cstdo, ' ')
  write(cstdo(lenstd:lenstd+5), '(a1,i5.5)') '.', myrank
  call rewnml(ifpar, jfpar)
  open(unit=jfpar, file=cstdo, &
    &  access='sequential', form='formatted')
  write(nfstdo, *) 'MESSAGE OUTPUT FOR RANK', myrank
  call rewnml(ifpar, jfpar)
  read (ifpar, nmrun, iostat=istat)
  call cstnml(jfpar, 'icedcoco', 'nmrun', istat)
!  write(jfpar, nmrun)
  write(nfstdo, *) 'Run name :'//crun

  ofinal = .false.
  oinit = .false.
  call tmstup( &
    &           tstrt,   tend,     dt)
  call restrt( &
    &           tstrt, &
    &              ub,     vb,     tb, &
    &              hb,   ubtb,   vbtb, &
    &               w,    amv,    ahv, &
    &              ab,    hib,    uib,    vib,    tib,    hsb, &
    &             tsi, &
    &              ft,  swabs,     fs, &
    &            taux,   tauy,   ptop )
  call ocstup( &
    &              ub,     vb,     tb, &
    &              hb,   ubtb,   vbtb,      w,      r, &
    &             amv,    ahv, &
    &              ft,   ptop, &
    &              dt )
  oinit = .true.
  call ocean ( &
    &              ua,     va,     ta, &
    &              ha,   ubta,   vbta, &
    &              ab,    hib,    uib,    vib,    tib,    hsb, &
    &              ub,     vb,     tb, &
    &              hb,   ubtb,   vbtb, &
    &              aa,    hia,    uia,    via,    tia,    hsa, &
    &               w,      r,    amv,    ahv, &
    &              ft,  swabs,     fs, &
    &            taux,   tauy,   ptop,    tsi, &
    &              nt,     tt,   itst,     ts,    its, &
    &            ntss,    tss, &
    &          oflout, oflstk )
  oinit = .false.
  call forsto( &
    &              ua,     va,     ta, &
    &              ha,   ubta,   vbta, &
    &              ub,     vb,     tb, &
    &              hb,   ubtb,   vbtb )
  call forsti( &
    &              aa,    hia,    uia,    via,    tia,    hsa, &
    &              ab,    hib,    uib,    vib,    tib,    hsb )
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
       &              ab,    hib,    uib,    vib,    tib,    hsb, &
       &              ub,     vb,     tb, &
       &              hb,   ubtb,   vbtb, &
       &              aa,    hia,    uia,    via,    tia,    hsa, &
       &               w,      r,    amv,    ahv, &
       &              ft,  swabs,     fs, &
       &            taux,   tauy,   ptop,    tsi, &
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
          &              aa,    hia,    uia,    via,    tia,    hsa, &
          &             tsi, &
          &              ft,  swabs,     fs, &
          &            taux,   tauy,   ptop, &
          &          orsout, orsrwd )
     else
        call finout( &
          &              tt,     nt, &
          &              ub,     vb,     tb, &
          &              hb,   ubtb,   vbtb, &
          &               w,    amv,    ahv, &
          &              ab,    hib,    uib,    vib,    tib,    hsb, &
          &             tsi, &
          &              ft,  swabs,     fs, &
          &            taux,   tauy,   ptop, &
          &          orsout, orsrwd )
     end if
     if (orsout) then
        ofinal = .true.
        call ocean ( &
          &           ua,     va,     ta, &
          &           ha,   ubta,   vbta, &
          &           ab,    hib,    uib,    vib,    tib,    hsb, &
          &           ub,     vb,     tb, &
          &           hb,   ubtb,   vbtb, &
          &           aa,    hia,    uia,    via,    tia,    hsa, &
          &            w,      r,    amv,    ahv, &
          &           ft,  swabs,     fs, &
          &         taux,   tauy,   ptop,    tsi, &
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
       &              aa,    hia,    uia,    via,    tia,    hsa, &
       &             tsi, &
       &              ft,  swabs,     fs, &
       &            taux,   tauy,   ptop, &
       &          orsout, orsrwd )
  else
     call finout( &
       &              tt,     nt, &
       &              ub,     vb,     tb, &
       &              hb,   ubtb,   vbtb, &
       &               w,    amv,    ahv, &
       &              ab,    hib,    uib,    vib,    tib,    hsb, &
       &             tsi, &
       &              ft,  swabs,     fs, &
       &            taux,   tauy,   ptop, &
       &          orsout, orsrwd )
  end if
  call ocean ( &
    &              ua,     va,     ta, &
    &              ha,   ubta,   vbta, &
    &              ab,    hib,    uib,    vib,    tib,    hsb, &
    &              ub,     vb,     tb, &
    &              hb,   ubtb,   vbtb, &
    &              aa,    hia,    uia,    via,    tia,    hsa, &
    &               w,      r,    amv,    ahv, &
    &              ft,  swabs,     fs, &
    &            taux,   tauy,   ptop,    tsi, &
    &              nt,     tt,   itst,     ts,    its, &
    &            ntss,    tss, &
    &          oflout, oflstk )
  call clcend('FINOUT')
  call clcout
  call parfin

end program icedcoco

! *********************************************************************

subroutine parset

! --- information -----------------------------------------------------
!
!  Initialization for inter-node communication using MPI
!
!  HISTORY
!     '00.05.30  H.Hasumi: from COCO2
!     '12.10.16  Y.Komuro: for COCO5.0
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
#ifndef OPT_IO_COCOMPI
  use bgs3d
  use bgsid
#endif
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
#ifndef OPT_IO_COCOMPI
  call gs3dst
  call gsidst
#endif
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
!     '12.10.16  Y.Komuro: for COCO5.0
!
! ---------------------------------------------------------------------

  use zocnod, only: &
    &   ierr   
  use ufile
  implicit none

#include "mpif.h"

#ifdef OPT_IO_COCOMPI
  call mpi_filcls_all
#endif
  call mpi_finalize(ierr)

  return
end subroutine parfin
