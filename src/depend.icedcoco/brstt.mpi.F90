module brstt

! --- information -----------------------------------------------------
!
!  Reading and writing restart file
!
!  HISTORY
!     '02.10.10  H.Hasumi: from COCO3.4
!     '07.04.23  H.Hasumi
!     '07.05.01  H.Hasumi: removal of negative salinity, which is
!                          sometime used for missing values
!     '07.09.25  H.Hasumi: multi-category sea ice
!     '07.09.27  H.Hasumi: 1-layer sea ice thermodynamics
!     '07.11.06  T.Suzuki: for ES
!     '08.11.25  Y.Komuro: default value of TSI changed
!                          eliminate unphysical negative AB(IJ, 0)
!     '10.04.14  M.Kurogi: (COCO4.4 tripolar code by Dr. Suzuki)
!     '12.09.10  H.Tatebe: for COCO5.0 in F90
!     '15.04.07  M.Kurogi: MPI-IO
! ---------------------------------------------------------------------

  use zocfil,   only  :   ncf

  implicit none
#include "mpif.h"
#include "coco.h"
  private

  integer(4),               save  ::  nfinit,      nfrest
  integer(4),               save  ::   idate(1:6)
  character(len=ncf)              ::  cfinit,      cfrest
  character(len=16),        save  ::  chead(1:64)
  data cfinit, cfrest / 'not-specified', 'not-specified' /

  integer, save :: mpi_fh_w, mpi_fh_r
  integer(kind=mpi_offset_kind), save :: disp, dispw=0
  integer :: icread

!  character(len=16),        save  ::  cheadtx, cheadty
!  character(len=16),        save  ::  cheadvx, cheadvy
  character(len=16),        save  ::   chrnum
  real(8), save  ::  dundef = -1.d20
 
  public  ::  restrt,  rstadd,  finadd, finout

  interface  print_stats
     module procedure &
          & print_stats_2d, &
          & print_stats_3d
  end interface print_stats

contains

  subroutine restrt(                                                  &
    &       tstrt,                                                    &
    &          ub,     vb,     tb,                                    &
    &          hb,   ubtb,   vbtb,                                    &
    &           w,    amv,    ahv,                                    & 
    &          ab,    hib,    uib,    vib,    tib,    hsb,            &
    &         tsi,                                                    & 
    &          ft,  swabs,     fs,                                    &
    &        taux,   tauy,   ptop)

    use zocdim,   only  :                                             &
    &      nxdim,  nydim,  nzdim,  ntdim,  nztdim,    nic,            &
    &     nxgdim, nygdim,  igstr,  jgstr,    kstr,                    &
    &        nxg,    nyg,     nz
    use zocnod,   only  :                                             &
    &     myrank,  iroot, mpi_comm_ogcm
    use zocphy,   only  :                                             &
    &       dtds
    use zocfil
    use ufile
    use ucaln
    use bshft
    use mpiio
    use zocgrd,   only  : hic
    use zocout,   only  :                                             &
    &     loglev
    implicit none

    real(8),   intent(in)     ::  tstrt
    real(8),   intent(inout)  ::     ub(nxdim, nydim, nzdim)
    real(8),   intent(inout)  ::     vb(nxdim, nydim, nzdim)
    real(8),   intent(inout)  ::     tb(nxdim, nydim, nzdim, ntdim)
    real(8),   intent(inout)  ::     hb(nxdim, nydim)
    real(8),   intent(inout)  ::   ubtb(nxdim, nydim)
    real(8),   intent(inout)  ::   vbtb(nxdim, nydim)
    real(8),   intent(inout)  ::      w(nxdim, nydim, nzdim)
    real(8),   intent(inout)  ::    amv(nxdim, nydim, nzdim)
    real(8),   intent(inout)  ::    ahv(nxdim, nydim, nzdim)
    real(8),   intent(inout)  ::     ab(nxdim, nydim, 0:nic)
    real(8),   intent(inout)  ::    hib(nxdim, nydim, 0:nic)
    real(8),   intent(inout)  ::    uib(nxdim, nydim)
    real(8),   intent(inout)  ::    vib(nxdim, nydim)
    real(8),   intent(inout)  ::    tib(nxdim, nydim, 0:nic)
    real(8),   intent(inout)  ::    hsb(nxdim, nydim, 0:nic)
    real(8),   intent(inout)  ::    tsi(nxdim, nydim, 0:nic)
    real(8),   intent(inout)  ::     ft(nxdim, nydim, ntdim)
    real(8),   intent(inout)  ::  swabs(nxdim, nydim)
    real(8),   intent(inout)  ::     fs(nxdim, nydim)
    real(8),   intent(inout)  ::   taux(nxdim, nydim)
    real(8),   intent(inout)  ::   tauy(nxdim, nydim)
    real(8),   intent(inout)  ::   ptop(nxdim, nydim)

!---- local varibles
    real(8)           ::     tt,    ttt
    real(8),    save  ::     si
    character(len=16) ::  cdate
    integer(4)        ::  ixdim,  jydim,  kzdim
    integer(4)        ::      i,      j,      k,      l
    integer(4), save  ::  ifpar,  jfpar,  istat
    integer(4), save  :: irstrt,   ierr

    namelist /nmfini/ cfinit, irstrt
    namelist /nmislt/ si
    data     si / 5.d0 /
    data irstrt / 0 /

    READ_NAMELIST( nmfini)
    READ_NAMELIST( nmislt)

    call mpi_filopn(mpi_fh_r, cfinit, 'READ')
    disp=0

    call mpi_read_chead(chead, mpi_fh_r, disp, icread)
    ub=0.d0
    if(icread == 1024) call mpi_read_3d(ub, mpi_fh_r,disp)
    if (loglev > 0) then
!       call print_stats(ub, chead(3), 'V')
       call print_stats(ub, chead(3))
    end if
    !$acc update device(ub)
    
    call mpi_read_chead(chead, mpi_fh_r, disp, icread)
    vb=0.d0
    if(icread == 1024) call mpi_read_3d(vb, mpi_fh_r,disp)
    if (loglev > 0) then
!       call print_stats(vb, chead(3), 'V')
       call print_stats(vb, chead(3))
    end if
    !$acc update device(vb)
    
    call mpi_read_chead(chead, mpi_fh_r, disp, icread)
    tb(:,:,:,1)=0.d0
    if(icread == 1024) call mpi_read_3d(tb, mpi_fh_r,disp)
    if (loglev > 0) then
       call print_stats(tb(:,:,:,1), chead(3))
    end if
    
    call mpi_read_chead(chead, mpi_fh_r, disp, icread)
    tb(:,:,:,2)=0.d0
    if(icread == 1024) call mpi_read_3d(tb(1,1,1,2), mpi_fh_r,disp)
    if (loglev > 0) then
       call print_stats(tb(:,:,:,2), chead(3))
    end if
    !$acc update device(tb)
    !$acc kernels default(present)
    do k=1,nzdim
    do j=1,nydim
    do i=1,nxdim
       tb(i,j,k,2)=abs(tb(i,j,k,2))
    end do
    end do
    end do
    !$acc end kernels
    
    call mpi_read_chead(chead, mpi_fh_r, disp, icread)
    hb=0.d0
    if(icread == 1024) call mpi_read_2d(hb, mpi_fh_r,disp)
    if (loglev > 0) then
       call print_stats(hb, chead(3))
    end if
    !$acc update device(hb)
    
    call mpi_read_chead(chead, mpi_fh_r, disp, icread)
    ubtb=0.d0
    if(icread == 1024) call mpi_read_2d(ubtb, mpi_fh_r,disp)
    !$acc update device(ubtb)
    
    call mpi_read_chead(chead, mpi_fh_r, disp, icread)
    vbtb=0.d0
    if(icread == 1024) call mpi_read_2d(vbtb, mpi_fh_r,disp)
    !$acc update device(vbtb)
    
    call mpi_read_chead(chead, mpi_fh_r, disp, icread)
    w=0.d0
    if(icread == 1024) call mpi_read_3d(w, mpi_fh_r,disp)
    !$acc update device(w)
    
    call mpi_read_chead(chead, mpi_fh_r, disp, icread)
    ab=0.d0
    if(icread == 1024) call mpi_read_id(ab, mpi_fh_r,disp)
    !$acc update device(ab)
    
    call mpi_read_chead(chead, mpi_fh_r, disp, icread)
    hib=0.d0
    do k = 1, nic
       hib(:,:,k) = hic(k)
    end do
    if(icread == 1024) call mpi_read_id(hib, mpi_fh_r,disp)
    !$acc update device(hib)
    
    call mpi_read_chead(chead, mpi_fh_r, disp, icread)
    uib=0.d0
    if(icread == 1024) call mpi_read_2d(uib, mpi_fh_r,disp)
    !$acc update device(uib)
    
    call mpi_read_chead(chead, mpi_fh_r, disp, icread)
    vib=0.d0
    if(icread == 1024) call mpi_read_2d(vib, mpi_fh_r,disp)
    !$acc update device(vib)
    
    call mpi_read_chead(chead, mpi_fh_r, disp, icread)
    tib=-0.1d0
    if(icread == 1024) call mpi_read_id(tib, mpi_fh_r,disp)
    !$acc update device(tib)
    
    call mpi_read_chead(chead, mpi_fh_r, disp, icread)
    hsb=0.d0
    if(icread == 1024) call mpi_read_id(hsb, mpi_fh_r,disp)
    !$acc update device(hsb)
    
    call mpi_read_chead(chead, mpi_fh_r, disp, icread)
    ft(:,:,1)=0.d0
    if(icread == 1024) call mpi_read_2d(ft, mpi_fh_r,disp)
    !$acc update device(ft)
    
    call mpi_read_chead(chead, mpi_fh_r, disp, icread)
    swabs=0.d0
    if(icread == 1024) call mpi_read_2d(swabs, mpi_fh_r,disp)
    !$acc update device(swabs)
    
    call mpi_read_chead(chead, mpi_fh_r, disp, icread)
    ft(:,:,2)=0.d0
    if(icread == 1024) call mpi_read_2d(ft(1,1,2), mpi_fh_r,disp)
    !$acc update device(ft)
    
    call mpi_read_chead(chead, mpi_fh_r, disp, icread)
    fs=0.d0
    if(icread == 1024) call mpi_read_2d(fs, mpi_fh_r,disp)
    !$acc update device(fs)
    
    call mpi_read_chead(chead, mpi_fh_r, disp, icread)
    taux=0.d0
    if(icread == 1024) call mpi_read_2d(taux, mpi_fh_r,disp)
    !$acc update device(taux)
    
    call mpi_read_chead(chead, mpi_fh_r, disp, icread)
    tauy=0.d0
    if(icread == 1024) call mpi_read_2d(tauy, mpi_fh_r,disp)
    !$acc update device(tauy)
    
    call mpi_read_chead(chead, mpi_fh_r, disp, icread)
    amv=0.d0
    if(icread == 1024) call mpi_read_3d(amv, mpi_fh_r,disp)
    !$acc update device(amv)
    
    call mpi_read_chead(chead, mpi_fh_r, disp, icread)
    ahv=0.d0
    if(icread == 1024) call mpi_read_3d(ahv, mpi_fh_r,disp)
    !$acc update device(ahv)
    
    call mpi_read_chead(chead, mpi_fh_r, disp, icread)
    ptop=0.d0
    if(icread == 1024) call mpi_read_2d(ptop, mpi_fh_r,disp)
    !$acc update device(ptop)
    
    call mpi_read_chead(chead, mpi_fh_r, disp, icread)
    tsi=dtds*si
    if(icread == 1024) call mpi_read_id(tsi, mpi_fh_r,disp)
    !$acc update device(tsi)

    do l = 3, ntdim
       call mpi_read_chead(chead, mpi_fh_r, disp, icread)
       tb(:,:,:,l)=0.d0
       if(icread == 1024) call mpi_read_3d(tb(1,1,1,l), mpi_fh_r,disp)
    end do
    !$acc update device(tb)
    do l = 3, ntdim
       call mpi_read_chead(chead, mpi_fh_r, disp, icread)
       ft(:,:,l)=0.d0
       if(icread == 1024) call mpi_read_2d(ft(1,1,l), mpi_fh_r,disp)
    end do
    !$acc update device(ft)

    call shift_pack_begin
    call shift1(   ub, nxdim, nydim,  nzdim, -1.d0, -1, -1)
    call shift1(   vb, nxdim, nydim,  nzdim, -1.d0, -1, -1)
    call shift1(   tb, nxdim, nydim, nztdim,  1.d0,  0,  0)
    call shift1(  amv, nxdim, nydim,  nzdim,  1.d0, -1, -1)
    call shift1(    w, nxdim, nydim,  nzdim,  1.d0,  0,  0)
    call shift1(  ahv, nxdim, nydim,  nzdim,  1.d0,  0,  0)
    call shift1(   hb, nxdim, nydim,      1,  1.d0,  0,  0)
    call shift1( ubtb, nxdim, nydim,      1, -1.d0, -1, -1)
    call shift1( vbtb, nxdim, nydim,      1, -1.d0, -1, -1)
    call shift1(   ab, nxdim, nydim,  nic+1,  1.d0,  0,  0)
    call shift1(  hib, nxdim, nydim,  nic+1,  1.d0,  0,  0)
    call shift1(  hsb, nxdim, nydim,  nic+1,  1.d0,  0,  0)
    call shift1(  tsi, nxdim, nydim,  nic+1,  1.d0,  0,  0)
    call shift1(  tib, nxdim, nydim,  nic+1,  1.d0,  0,  0)
    call shift1(  uib, nxdim, nydim,      1, -1.d0, -1, -1)
    call shift1(  vib, nxdim, nydim,      1, -1.d0, -1, -1)
    call shift1(   ft, nxdim, nydim,  ntdim,  1.d0,  0,  0)
    call shift1(   fs, nxdim, nydim,      1,  1.d0,  0,  0)
    call shift1(swabs, nxdim, nydim,      1,  1.d0,  0,  0)
    call shift1( ptop, nxdim, nydim,      1,  1.d0,  0,  0)
    call shift1( taux, nxdim, nydim,      1, -1.d0, -1, -1)
    call shift1( tauy, nxdim, nydim,      1, -1.d0, -1, -1)
    call shift_pack_end
    call shift_unpack(   ub,  1)
    call shift_unpack(   vb,  2)
    call shift_unpack(   tb,  3)
    call shift_unpack(  amv,  4)
    call shift_unpack(    w,  5)
    call shift_unpack(  ahv,  6)
    call shift_unpack(   hb,  7)
    call shift_unpack( ubtb,  8)
    call shift_unpack( vbtb,  9)
    call shift_unpack(   ab, 10)
    call shift_unpack(  hib, 11)
    call shift_unpack(  hsb, 12)
    call shift_unpack(  tsi, 13)
    call shift_unpack(  tib, 14)
    call shift_unpack(  uib, 15)
    call shift_unpack(  vib, 16)
    call shift_unpack(   ft, 17)
    call shift_unpack(   fs, 18)
    call shift_unpack(swabs, 19)
    call shift_unpack( ptop, 20)
    call shift_unpack( taux, 21)
    call shift_unpack( tauy, 22)

    !$acc kernels default(present)
    ab(1:nxdim,1:nydim,0) = 1.d0
    do k = 1, nic
       do j = 1, nydim
          do i = 1, nxdim
             ab(i, j, 0) = ab(i, j, 0) - ab(i, j, k)
          end do
       end do
    end do

!   because of the finit precision, AB(IJ, 0) could be a small negative
!   value, which causes some problems
    do j = 1, nydim
       do i = 1, nxdim
          ab(i, j, 0) = max(0.d0, ab(i, j ,0))
       end do
    end do
    !$acc end kernels
    
    if ( myrank == iroot ) then
       cdate = chead(50)
       read(cdate, '(i6.6,5i2.2)') idate
       call cyh2ss(  ttt, idate )
       write(jfpar, *) ' start time :', idate
       if ( irstrt /= 0 ) then
          if ( ttt /= tstrt ) then
             write(jfpar, *) '*** start time error ***'
             call mpi_abort(mpi_comm_ogcm, 1, ierr)
          end if
       end if
    end if

  end subroutine restrt
    
! =====================================================================
  subroutine rstadd(                                                  &
    &      additm,   oeof,                                            &
    &       ixdim,  jydim,  kzdim,                                    &
#ifndef OPT_TRIPOLE
    &      ccitem,   clas )
#else
    &      ccitem,   clas,                                            &
    &        fact,   ioff,   joff )
#endif

    use zocdim,   only  :                                             &
    &      nxdim,  nydim,  nzdim,    nic,                             &
    &      igstr,  jgstr,   kstr,                                     &
    &        nxg,    nyg,     nz
    use zocnod,   only  :                                             &
    &     myrank,  iroot
    use bshft
    use mpiio

    implicit none

    integer(4),   intent(in)     ::   ixdim,  jydim,  kzdim
    real(8),      intent(inout)  ::  additm(ixdim, jydim, kzdim)
    logical,      intent(inout)  ::    oeof
    character(*), intent(in)     ::  ccitem,   clas
#ifdef OPT_TRIPOLE
! set 1 for scalar, set -1 for vector
    real(8),      intent(in)    ::    fact  
!(IOFF,JOFF)=(0,0)for T-cell, (-1,-1)for V-cell 
    integer(4),   intent(in)    ::    ioff,   joff 
#endif

!---- local variables
    integer(4)   ::     i,     j,     k
    integer(4)   ::  ierr

    shift_gpu=.false.
    
    if ( clas(1:3) == 'OCN' ) then
       oeof = .false.
       call mpi_read_chead(chead, mpi_fh_r, disp, icread)
       if(icread .ne. 1024) then
          oeof = .true.
       end if
       if(.not. oeof) then
          additm(:,:,:kstr-1)=0.d0
          additm(:,:,kstr+nz:)=0.d0
          call mpi_read_3d(additm, mpi_fh_r,disp)

#ifdef OPT_TRIPOLE
          call shift1(additm, nxdim, nydim, nzdim,                    &
    &                          fact,  ioff,  joff  )
#else
          call shift1(additm, nxdim, nydim, nzdim)
#endif
       end if
    else if ( clas(1:3) == 'ICE' ) then
       oeof = .false.
       call mpi_read_chead(chead, mpi_fh_r, disp, icread)
       if(icread .ne. 1024) then
          oeof = .true.
       end if
       if(.not. oeof) then
          call mpi_read_id(additm, mpi_fh_r,disp)

#ifdef OPT_TRIPOLE
          call shift1(additm, nxdim, nydim, nic+1,                    &
    &                          fact,  ioff,  joff  )
#else
          call shift1(additm, nxdim, nydim, nic+1)
#endif
       end if
    else
       oeof = .false.
       call mpi_read_chead(chead, mpi_fh_r, disp, icread)
       if(icread .ne. 1024) then
          oeof = .true.
       end if
       if(.not. oeof) then
          call mpi_read_2d(additm, mpi_fh_r,disp)
#ifdef OPT_TRIPOLE
          call shift1(additm, nxdim, nydim,     1,                    &
    &                          fact,  ioff,  joff)
#else
          call shift1(additm, nxdim, nydim, 1)
#endif
       end if
    end if

  end subroutine rstadd

! =====================================================================

  subroutine finout(                                                  &
    &         tt, ntstep,                                             &
    &         ub,     vb,     tb,                                     &
    &         hb,   ubtb,   vbtb,                                     &
    &          w,    amv,    ahv,                                     &
    &         ab,    hib,    uib,    vib,    tib,    hsb,             &
    &        tsi,                                                     &
    &         ft,  swabs,     fs,                                     &
    &       taux,   tauy,   ptop,                                     &
    &     orsout, orsrwd )

    use zocdim,   only  :                                             &
    &      nxdim,  nydim,  nzdim,  ntdim,     nic,                    &
    &     nxgdim, nygdim,  igstr,  jgstr,    kstr,                    &
    &        nxg,    nyg,     nz,   nxyg,   nxyzg
    use zocnod,   only  :                                             &
    &     myrank,  iroot, mpi_comm_ogcm
    use ufile
    use ucaln
    use mpiio

    implicit none

    real(8),    intent(in)     ::  tt
    integer(4), intent(in)     ::  ntstep
    real(8),    intent(inout)  ::      ub(nxdim, nydim, nzdim)
    real(8),    intent(inout)  ::      vb(nxdim, nydim, nzdim)
    real(8),    intent(inout)  ::      tb(nxdim, nydim, nzdim, ntdim)
    real(8),    intent(inout)  ::      hb(nxdim, nydim)
    real(8),    intent(inout)  ::    ubtb(nxdim, nydim)
    real(8),    intent(inout)  ::    vbtb(nxdim, nydim)
    real(8),    intent(inout)  ::       w(nxdim, nydim, nzdim)
    real(8),    intent(inout)  ::     amv(nxdim, nydim, nzdim)
    real(8),    intent(inout)  ::     ahv(nxdim, nydim, nzdim)
    real(8),    intent(inout)  ::      ab(nxdim, nydim, 0:nic)
    real(8),    intent(inout)  ::     hib(nxdim, nydim, 0:nic)
    real(8),    intent(inout)  ::     uib(nxdim, nydim)
    real(8),    intent(inout)  ::     vib(nxdim, nydim)
    real(8),    intent(inout)  ::     tib(nxdim, nydim, 0:nic)
    real(8),    intent(inout)  ::     hsb(nxdim, nydim, 0:nic)
    real(8),    intent(inout)  ::     tsi(nxdim, nydim, 0:nic)
    real(8),    intent(inout)  ::      ft(nxdim, nydim, ntdim)
    real(8),    intent(inout)  ::   swabs(nxdim, nydim)
    real(8),    intent(inout)  ::      fs(nxdim, nydim)
    real(8),    intent(inout)  ::    taux(nxdim, nydim)
    real(8),    intent(inout)  ::    tauy(nxdim, nydim)
    real(8),    intent(inout)  ::    ptop(nxdim, nydim)
    logical,    intent(in)     ::  orsout,    orsrwd

!---- local variables
    integer(4)        ::      i,      j,      k,     l
    integer(4), save  ::  ifpar,  jfpar,  istat
    integer(4)        ::   ierr
    character(len=ncf) :: crun = '(RUN NAME WAS NOT SET)'

    character(len=8)  :: hdate
    character(len=10) :: htime
    character(len=5)  :: hzone
    integer :: ivalues(1:8)
    character(len=16) :: citem
    logical,    save  ::  ofirst = .true.

    namelist /nmfrst/ cfrest
    namelist /nmrun/ crun

    if ( ofirst ) then
       READ_NAMELIST( nmfrst )
       READ_NAMELIST( nmrun  )
       call mpi_filopn(mpi_fh_w, cfrest, 'WRITE')
       ofirst = .false.
    end if

    if (crun(1:1) == '(') then
       chrnum = 'COCO stand-alone'
    else
       chrnum = crun(1:16)
    end if

    if ( .not. orsout ) return

    do i = 1, 64
       write(chead(i),'(16x)')
    end do
    if (orsrwd) then
       dispw=0
    end if
    call css2yh(  idate, tt)
    write(chead(1), '(i16)') 9010
    chead(2) = chrnum
    write(chead(25), '(i16)') nint(tt / 3.6d3)
    chead(26) = 'HOUR'
    chead(38) = 'UR8'
    write(chead(27), '(i4.4,2i2.2,1x,3i2.2,1x)') idate(1:6)
    write(chead(48), '(i4.4,2i2.2,1x,3i2.2,1x)') idate(1:6)
    write(chead(49), '(i4.4,2i2.2,1x,3i2.2,1x)') idate(1:6)
    write(chead(50), '(i6.6,5i2.2)') idate(1:6)
    write(chead(30), '(i16)') 1
    write(chead(31), '(i16)') nxg
    write(chead(33), '(i16)') 1
    write(chead(34), '(i16)') nyg
    write(chead(36), '(i16)') 1
    write(chead(39), '(e16.7)') dundef
    chead(40) = chead(39)
    chead(41) = chead(39)
    chead(42) = chead(39)
    chead(43) = chead(39)
    write(chead(44), '(i16)') 1
    write(chead(46), '(i16)') 0
    write(chead(47), '(e16.7)') 0.d0
    chead(61) = 'COCO'
    chead(63) = 'COCO'
    call date_and_time(hdate, htime, hzone, ivalues)
    write(chead(60), '(i4.4,2i2.2,1x,3i2.2,1x)') ivalues(1:3), ivalues(5:7)
    chead(62) = chead(60)
    !$acc update self(ub, vb, tb, hb,ubtb, vbtb, w, ab, hib, uib, vib, tib, hsb, ft, swabs, fs)
    !$acc update self(taux, tauy, amv, ahv, ptop, tsi)
    call edhead('UO', 'ocean zonal velocity', 'cm/s', 'OCLVTV')
    call mpi_write_header(chead, mpi_fh_w, dispw)
    call mpi_write_3d(ub, mpi_fh_w,dispw)    
    
    call edhead('VO', 'ocean meridional velocity', 'cm/s', 'OCLVTV')
    call mpi_write_header(chead, mpi_fh_w, dispw)
    call mpi_write_3d(vb, mpi_fh_w,dispw)    

    call edhead('TO', 'ocean temperature', 'degC', 'OCLVTT')
    call mpi_write_header(chead, mpi_fh_w, dispw)
    call mpi_write_3d(tb, mpi_fh_w,dispw)    

    call edhead('SO', 'ocean salinity', 'psu', 'OCLVTT')
    call mpi_write_header(chead, mpi_fh_w, dispw)
    call mpi_write_3d(tb(1,1,1,2), mpi_fh_w,dispw)    

    call edhead('SHO', 'sea surface height', 'cm', 'OCSFCT' )
    call mpi_write_header(chead, mpi_fh_w, dispw)
    call mpi_write_2d(hb, mpi_fh_w,dispw)    

    call edhead('UBTO', 'ocean zonal transport', 'cm^2/s', 'OCSFCV')
    call mpi_write_header(chead, mpi_fh_w, dispw)
    call mpi_write_2d(ubtb, mpi_fh_w,dispw)    
       
    call edhead('VBTO', 'ocean meridional transport', 'cm^2/s', 'OCSFCV')
    call mpi_write_header(chead, mpi_fh_w, dispw)
    call mpi_write_2d(vbtb, mpi_fh_w,dispw)    

    call edhead('WO', 'ocean vertical velocity', 'cm/s', 'OCLVMT')
    call mpi_write_header(chead, mpi_fh_w, dispw)
    call mpi_write_3d(w, mpi_fh_w,dispw)    

    call edhead('AI', 'ice concentration', 'N.D.', 'OCICET')
    call mpi_write_header(chead, mpi_fh_w, dispw)
    call mpi_write_id(ab, mpi_fh_w,dispw)    
       
    call edhead('HI', 'ice thickness', 'cm', 'OCICET')
    call mpi_write_header(chead, mpi_fh_w, dispw)
    call mpi_write_id(hib, mpi_fh_w,dispw)    
   
    call edhead('UI', 'ice zonal velocity', 'cm/s', 'OCSFCV')
    call mpi_write_header(chead, mpi_fh_w, dispw)
    call mpi_write_2d(uib, mpi_fh_w,dispw)    
   
    call edhead('VI', 'ice meridional velocity', 'cm/s', 'OCSFCV')
    call mpi_write_header(chead, mpi_fh_w, dispw)
    call mpi_write_2d(vib, mpi_fh_w,dispw)    

    call edhead('TI', 'ice temperature', 'degC', 'OCICET')
    call mpi_write_header(chead, mpi_fh_w, dispw)
    call mpi_write_id(tib, mpi_fh_w,dispw)    
    
    call edhead('HS', 'snow thickness', 'cm', 'OCICET')
    call mpi_write_header(chead, mpi_fh_w, dispw)
    call mpi_write_id(hsb, mpi_fh_w,dispw)    

    call edhead('FT', 'sea surface temperature flux', 'K cm/s', 'OCSFCT')
    call mpi_write_header(chead, mpi_fh_w, dispw)
    call mpi_write_2d(ft, mpi_fh_w,dispw)    

    call edhead('SWABS', 'absorbed shortwave', 'erg/cm^2/s', 'OCSFCT')
    call mpi_write_header(chead, mpi_fh_w, dispw)
    call mpi_write_2d(swabs, mpi_fh_w,dispw)    

    call edhead('FW', 'sea surface freshwater flux', 'cm/s', 'OCSFCT')
    call mpi_write_header(chead, mpi_fh_w, dispw)
    call mpi_write_2d(ft(1,1,2), mpi_fh_w,dispw)    
    
    call edhead('FS', 'sea surface salinity flux', 'psu cm/s', 'OCSFCT')
    call mpi_write_header(chead, mpi_fh_w, dispw)
    call mpi_write_2d(fs, mpi_fh_w,dispw)    

    call edhead('TAUX', 'zonal wind stress', 'dyn/cm^2', 'OCSFCV')
    call mpi_write_header(chead, mpi_fh_w, dispw)
    call mpi_write_2d(taux, mpi_fh_w,dispw)    

    call edhead('TAUY', 'meridional wind stress', 'dyn/cm^2', 'OCSFCV')
    call mpi_write_header(chead, mpi_fh_w, dispw)
    call mpi_write_2d(tauy, mpi_fh_w,dispw)    
    
    call edhead('AMV', 'ocean vertical viscosity', 'cm^2/s', 'OCLVMV')
    call mpi_write_header(chead, mpi_fh_w, dispw)
    call mpi_write_3d(amv, mpi_fh_w,dispw)    

    call edhead('AHV', 'ocean vertical diffusivity', 'cm^2/s', 'OCLVMT')
    call mpi_write_header(chead, mpi_fh_w, dispw)
    call mpi_write_3d(ahv, mpi_fh_w,dispw)    

    call edhead('PTOP', 'sea surface presssure (inc. sea ice)', 'dyn/cm^2', 'OCSFCT')
    call mpi_write_header(chead, mpi_fh_w, dispw)
    call mpi_write_2d(ptop, mpi_fh_w,dispw)    

    call edhead('TSI', 'sea ice surface temperature', 'degC', 'OCICET')
    call mpi_write_header(chead, mpi_fh_w, dispw)
    call mpi_write_id(tsi, mpi_fh_w,dispw)    
    
    do l = 3, ntdim
       write(citem, '(a,i2.2)') 'TRACER', l
       call edhead(citem, '', '', 'OCLVTT')
       call mpi_write_header(chead, mpi_fh_w, dispw)
       call mpi_write_3d(tb(1,1,1,l), mpi_fh_w,dispw)    
    end do
       
    do l = 3, ntdim
       write(citem, '(a,i2.2)') 'TRCFLX', l
       call edhead(citem, '', '', 'OCSFCT')
       call mpi_write_header(chead, mpi_fh_w, dispw)
       call mpi_write_2d(ft(1,1,l), mpi_fh_w,dispw)    
    end do
    
    if ( myrank == iroot ) then
       write(jfpar, *) '*** Write restart file ***'
       write(jfpar, *) ' time :', idate
       write(jfpar, *) ' step :', ntstep
    end if
    call mpi_barrier(mpi_comm_ogcm, ierr)

  end subroutine finout

! =====================================================================

  subroutine finadd(                                                  &
    &             additm,  ixdim,  jydim,  kzdim,                     &
    &             ccitem,   clas  ) 

    use zocdim,   only  :                                             &
    &      igstr,  jgstr,   kstr,                                     &
    &        nxg,    nyg,     nz,   nxyg,  nxyzg, nic
    use zocnod,   only  :                                             &
    &     myrank,  iroot
    use mpiio
    implicit none

    integer(4),   intent(in)     ::   ixdim,  jydim,  kzdim
    real(8),      intent(inout)  ::  additm(ixdim, jydim, kzdim)
    character(*), intent(in)     ::  ccitem,   clas
    
    !$acc update self(additm)
    chead(3) = ccitem
    write(chead(14), '(16x)')
    write(chead(15), '(16x)')
    write(chead(16), '(16x)')

#ifdef OPT_TRIPOLE
    write(chead(29), '(a,i0)') 'OCLONTPT', nxg
    write(chead(32), '(a,i0)') 'OCLATTPT', nyg
#else
    write(chead(29), '(i15,a1)') nxg, 'X'
    write(chead(32), '(i15,a1)') nyg, 'Y'
#endif

    if ( clas(1:3) == 'OCN' ) then
       write(chead(35), '(a,i0)') 'OCDEPT', nz
       write(chead(37), '(i16)') nz
       write(chead(64), '(i16)') nxyzg
       call mpi_write_header(chead, mpi_fh_w, dispw)
       call mpi_write_3d(additm, mpi_fh_w,dispw)    
    else if ( clas(1:3) == 'ICE' ) then
       write(chead(35), '(a)') 'NUMBER1000'
       write(chead(37), '(i16)') nic
       write(chead(64), '(i16)') nxyg*nic
       call mpi_write_header(chead, mpi_fh_w, dispw)
       call mpi_write_id(additm, mpi_fh_w,dispw)    
    else
       write(chead(35), '(a)') 'SFC1'
       write(chead(37), '(i16)') 1
       write(chead(64), '(i16)') nxyg
       call mpi_write_header(chead, mpi_fh_w, dispw)
       call mpi_write_2d(additm, mpi_fh_w,dispw)
    end if

  end subroutine finadd

! =====================================================================

  subroutine edhead(                 &
       &             ccitem,         &
       &              htitl,  hunit, &
       &              cclas)
    use zocdim, only : nxg, nyg, nxyg, nxyzg, nic, nz

    character(*), intent(in) :: ccitem,  cclas,  htitl,  hunit
    character(len=32) ::  ctitl

    ctitl = htitl

    chead(3) = ccitem
    chead(14) = ctitl(1:16)
    chead(15) = ctitl(17:32)
    chead(16) = hunit

#ifdef OPT_TRIPOLE
    if (cclas(6:6) == 'V') then
       write(chead(29), '(a,i0)') 'OCLONTPV', nxg
       write(chead(32), '(a,i0)') 'OCLATTPV', nyg
    else
       write(chead(29), '(a,i0)') 'OCLONTPT', nxg
       write(chead(32), '(a,i0)') 'OCLATTPT', nyg
    end if
#else
    write(chead(29), '(i15,a1)') nxg, 'X'
    write(chead(32), '(i15,a1)') nyg, 'Y'
#endif

    if (cclas(3:5) == 'SFC') then
       write(chead(35), '(a)') 'SFC1'
       write(chead(37), '(i16)') 1
       write(chead(64), '(i16)') nxyg
    else if (cclas(3:5) == 'ICE') then
       write(chead(35), '(a)') 'NUMBER1000'
       write(chead(37), '(i16)') nic
       write(chead(64), '(i16)') nxyg*nic
    else if (cclas(3:5) == 'LVT' .or. cclas(3:5) == 'LVM') then
       write(chead(35), '(2a,i0)') 'OCDEP', cclas(5:5), nz
       write(chead(37), '(i16)') nz
       write(chead(64), '(i16)') nxyzg
    else
       write(chead(35), '(a,i0)') 'OCDEPT', nz
       write(chead(37), '(i16)') nz
       write(chead(64), '(i16)') nxyzg
    end if

  end subroutine edhead

! =====================================================================
  subroutine print_stats_2d(data, cname, cpos)

    use ufile
    use zocdim,   only  :                                             &
    &     istr, iend, jstr, jend, kstr, nxdim
    use zocmsk,   only  :                                             &
    &      amskt,  amskv
    use zocnod,   only  :                                             &
    &     mpi_comm_ogcm
!#ifdef OPT_TOUZA
!    use TOUZA_Std_log, only: msg
!#endif
    real(8),          intent(in)           :: data(:,:)
    character(len=*), intent(in)           :: cname
    character(1),     intent(in), optional :: cpos
    
    logical, save              :: ofirst = .true.
    logical, save, allocatable :: omask(:, :)
    integer, save              :: dnumg

    real(8) ::  dmax,  dmin, dmaxg, dming
    real(8) ::  dsum,  dave, dsumg, daveg
    real(8) ::  dvar, dvarg, dstdg
    integer ::  dnum
    integer :: ifpar, jfpar,  ierr
    integer ::    ij,     i,     j
    logical :: oposv
    
    logical, save ::  omsk = .false.  !! if true, does not work correctly

    if (ofirst) then
       allocate(omask(size(data,1),size(data,2)))
       omask(:,:) = .false.
       if (omsk) then
          oposv = .false.
          if (present(cpos)) then
             if (cpos == 'V') then
                oposv = .true.
             end if
          end if
          if (oposv) then
             do j = jstr, jend
             do i = istr, iend
                ij = nxdim*(j-1) + i
                if (amskv(ij, kstr) > 0.D0) then
                   omask(i, j) = .true.
                end if
             end do
             end do
          else
             do j = jstr, jend
             do i = istr, iend
                ij = nxdim*(j-1) + i
                if (amskt(ij, kstr) > 1.D0) then
                   omask(i, j) = .true.
                end if
             end do
             end do
          end if
       else
          do j = jstr, jend
          do i = istr, iend
             omask(i, j) = .true.
          end do
          end do
       end if
       dnum=count(omask)
       call mpi_allreduce( &
            &  dnum, dnumg, 1, mpi_integer4, &
            &  mpi_sum, mpi_comm_ogcm, ierr)
       ofirst = .false.
    end if
    
    dmax=maxval(data, omask)
    dmin=minval(data, omask)
    call mpi_allreduce( &
         &  dmax, dmaxg, 1, mpi_real8, &
         &  mpi_max, mpi_comm_ogcm, ierr)
    call mpi_allreduce( &
         &  dmin, dming, 1, mpi_real8, &
         &  mpi_min, mpi_comm_ogcm, ierr)
    dsum=sum(data, omask)
    call mpi_allreduce( &
         &  dsum, dsumg, 1, mpi_real8, &
         &  mpi_sum, mpi_comm_ogcm, ierr)
    daveg=dsumg/dble(dnumg)
    dvar=sum(data**2.d0,omask)
    call mpi_allreduce( &
         &  dvar, dvarg, 1, mpi_real8, &
         &  mpi_sum, mpi_comm_ogcm, ierr)
    dvarg=dvarg/dble(dnumg)-daveg**2.d0
    dstdg=sqrt(dvarg)

    call rewnml(ifpar, jfpar)
!#ifdef OPT_TOUZA
!    call msg('("MAX, MIN, AVE, SD, NUM of '//trim(cname)//' : ")', dmaxg, dming, daveg, dstdg, dnumg, jfpar)
!#else
    write(jfpar,*) 'MAX, MIN, AVE, SD, NUM of '//trim(cname)//' : ', dmaxg, ',', dming, ',', daveg, ',', dstdg, ',', dnumg
!    call flush(jfpar)
!#endif
    call flush(jfpar)
    
  end subroutine print_stats_2d

  subroutine print_stats_3d(data, cname, cpos)

    use ufile
    use zocdim,   only  :                                             &
    &     istr, iend, jstr, jend, kstr, kend, nxdim, nxydim
    use zocmsk,   only  :                                             &
    &      amskt,  amskv
    use zocnod,   only  :                                             &
    &     mpi_comm_ogcm
!#ifdef OPT_TOUZA
!    use TOUZA_Std_log, only: msg
!#endif
    
    real(8),          intent(in)           :: data(:,:,:)
    character(len=*), intent(in)           :: cname
    character(1),     intent(in), optional :: cpos
    
    logical, save              :: ofirst = .true.
    logical, save, allocatable :: omask(:,:,:)
    integer, save              :: dnumg

    real(8) ::  dmax,  dmin, dmaxg, dming
    real(8) ::  dsum,  dave, dsumg, daveg
    real(8) ::  dvar, dvarg, dstdg
    integer ::  dnum
    integer :: ifpar, jfpar,  ierr
    integer ::    ij,     i,     j,     k
    logical :: oposv
    
    logical, save ::  omsk = .false.  !! if true, does not work correctly

    if (ofirst) then
       allocate(omask(size(data,1),size(data,2),size(data,3)))
       omask(:,:,:) = .false.
       oposv = .false.
       if (omsk) then
          if (present(cpos)) then
             if (cpos == 'V') then
                oposv = .true.
             end if
          end if
          if (oposv) then
             do j = jstr, jend
             do i = istr, iend
                ij = nxdim*(j-1) + i
                do k = kstr, kend
                   if (amskv(ij, k) > 0.D0) then
                      omask(i, j, k) = .true.
                   end if
                end do
             end do
             end do
          else
             do j = jstr, jend
             do i = istr, iend
                ij = nxdim*(j-1) + i
                do k = kstr, kend
                   if (amskt(ij, k) > 0.D0) then
                      omask(i, j, k) = .true.
                   end if
                end do
             end do
             end do
          end if
       else
          do k = kstr, kend
          do j = jstr, jend
          do i = istr, iend
             omask(i, j, k) = .true.
          end do
          end do
          end do
       end if
       dnum=count(omask)
       call mpi_allreduce( &
            &  dnum, dnumg, 1, mpi_integer4, &
            &  mpi_sum, mpi_comm_ogcm, ierr)
       ofirst = .false.
    end if
    
    dmax=maxval(data, omask)
    dmin=minval(data, omask)
    call mpi_allreduce( &
         &  dmax, dmaxg, 1, mpi_real8, &
         &  mpi_max, mpi_comm_ogcm, ierr)
    call mpi_allreduce( &
         &  dmin, dming, 1, mpi_real8, &
         &  mpi_min, mpi_comm_ogcm, ierr)
    dsum=sum(data, omask)
    call mpi_allreduce( &
         &  dsum, dsumg, 1, mpi_real8, &
         &  mpi_sum, mpi_comm_ogcm, ierr)
    daveg=dsumg/dble(dnumg)
    dvar=sum(data**2.d0,omask)
    call mpi_allreduce( &
         &  dvar, dvarg, 1, mpi_real8, &
         &  mpi_sum, mpi_comm_ogcm, ierr)
    dvarg=dvarg/dble(dnumg)-daveg**2.d0
    dstdg=sqrt(dvarg)

    call rewnml(ifpar, jfpar)
!#ifdef OPT_TOUZA
!    call msg('("MAX, MIN, AVE, SD, NUM of '//trim(cname)//' : ")', dmaxg, dming, daveg, dstdg, dnumg, jfpar)
!#else
    write(jfpar,*) 'MAX, MIN, AVE, SD, NUM of '//trim(cname)//' : ', dmaxg, ',', dming, ',', daveg, ',', dstdg, ',', dnumg
!#endif
    call flush(jfpar)

  end subroutine print_stats_3d

end module brstt
