module brstt

  use zocfil,   only  :   ncf
  use zocdim,   only  :   nxg, nyg, nz, nxgdim, nygdim, nzdim, ntdim

  implicit none

#include "mpif.h"

  private

  integer, save :: mpi_fh_w, mpi_fh_r
  integer(kind=mpi_offset_kind), save :: disp, dispw=0
  integer :: icread
  character(len=16), save :: chrnum
  real(8), save :: dundef = -1.d20
  
  integer(4) :: idate(6)
  character(len=16) :: chead(1:64)
  character(len=16) :: cdate
  integer(4), save ::  ifpar,  jfpar
  
  integer(4), save ::   nfinit, nfrest
   
  character(len=ncf) :: cfinit, cfrest
  integer(4), save ::  irstrt
  namelist /nmfini/ cfinit, irstrt
  namelist /nmfrst/ cfrest
  
  data cfinit, cfrest / 'not-specified', 'not-specified' /
  data irstrt / 0 /
  public :: restrt, finout, rstadd, finadd

  character(len=16) :: ctrnam(ntdim), cftnam(ntdim)
  data ctrnam(1:2) / 'TO', 'SO' /
  data cftnam(1:2) / 'FT', 'FW' / 
#ifdef OPT_OECO2
  data ctrnam(3:ntdim) /                                &
       &         'NO3',  'PHY',  'ZOO',  'DET',  'CA',  &
       &       'CACO3', 'TCO2',  'ALK',   'O2', 'FED',  &
       &       'DETFE', 'DIAZ',  'PO4',  'N2O', 'AGE' /
#endif  

  interface  print_stats
     module procedure &
          & print_stats_2d, &
          & print_stats_3d
  end interface print_stats

contains
  
  subroutine restrt(          &
#ifdef OPT_NPZD
    &   pco2o,     ff,        &
#endif /* OPT_NPZD */
    &   tstrt,     t,    ft )

    use zocdim,   only  :                                             &
    &      nxdim,  nydim,  nzdim,  ntdim,  nztdim,    nic,            &
    &     nxgdim, nygdim,  igstr,  jgstr,    kstr,                    &
    &        nxg,    nyg,     nz
    use zocnod,   only  :                                             &
    &     myrank,  iroot, mpi_comm_ogcm
    use zocfil
    use zocout,   only  : loglev
    use ufile
    use ucaln
    use bshft
    use mpiio

    implicit none
   
#ifdef OPT_NPZD
    real(8), intent(inout) :: pco2o(nxdim,nydim), ff(nxdim,nydim)
#endif
    real(8), intent(inout) ::     t(nxdim, nydim, nzdim, ntdim)
    real(8), intent(inout) ::    ft(nxdim, nydim, ntdim)
    real(8), intent(in)    :: tstrt
    real(8)    ::    ttt
    real(8)    ::    dmaxg
    integer(4) ::      l
    integer(4), save :: istat
    integer(4) :: ierr

    call rewnml(ifpar, jfpar)
    read(ifpar, nmfini, iostat=istat)
    call cstnml( jfpar, 'restrt', 'nfini', istat )
    write(jfpar, nmfini)
    
    call mpi_filopn(mpi_fh_r, cfinit, 'READ')
    disp = 0
    if (myrank == iroot) then
       write(jfpar, *) '*** read initial condition file ***'
    end if
      
    t  = 0.d0
    ft = 0.d0
#ifdef OPT_NPZD
    pco2o = 0.d0
    ff    = 0.d0
#endif

    do l = 3, ntdim
       call mpi_read_chead(chead, mpi_fh_r, disp, icread)
       if (icread == 1024) then
          call mpi_read_3d(t(1, 1, 1, l), mpi_fh_r, disp)
       end if
#ifdef OPT_OECO2
       if ( chead(3)(1:3) == 'NO3' .or. l == 3 ) then
          call print_stats(t(:, :, :, l), chead(3), 'T', dmaxg)
          if ( dmaxg < 1.d-5 ) then
             write(jfpar,*) 'I.C. is abonormal: max. of NO3 < 1.d-5'
             call flush(jfpar)
             stop
          end if
       endif
#else
       if (loglev > 0) then
          call print_stats(t(:, :, :, l), chead(3))
       end if
#endif
    end do
    !$acc update device(t)
    
    do l = 3, ntdim
       call mpi_read_chead(chead, mpi_fh_r, disp, icread)
       if (icread == 1024) then
          call mpi_read_2d(ft(1, 1, l), mpi_fh_r, disp)
       end if
       if (loglev > 0) then
          call print_stats(ft(:, :, l), chead(3))
       end if
    end do
    !$acc update device(ft)

#ifdef OPT_NPZD
    call mpi_read_chead(chead, mpi_fh_r, disp, icread)
    if (icread == 1024) then
       call mpi_read_2d(pco2o, mpi_fh_r, disp)
       if (loglev > 0) then
          call print_stats(pco2o, chead(3))
       end if
    end if
    call mpi_read_chead(chead, mpi_fh_r, disp, icread)
    if (icread == 1024) then
       call mpi_read_2d(ff, mpi_fh_r, disp)
       if (loglev > 0) then
          call print_stats(ff, chead(3))
       end if
    end if
    !$acc update device(pco2o, ff)
#endif /* OPT_NPZD */

#ifdef OPT_TRIPOLE
    call shift1(     t,                     &
                 nxdim,   nydim, nztdim,    &
    &             1.d0,       0,      0 )
    call shift1(    ft,                     &
    &            nxdim,   nydim,  ntdim,    &
    &             1.d0,       0,      0 )
#ifdef OPT_NPZD
    call shift2( pco2o,      ff,            &
    &            nxdim,   nydim,      1,    &
    &             1.d0,       0,      0 )
#endif
#else
    call shift1(      t,                    &
     &            nxdim,  nydim, nztdim)
    call shift1(     ft,                    &
     &            nxdim,  nydim,  ntdim)
#ifdef OPT_NPZD
    call shift1( pco2o,                     &
    &            nxdim,   nydim,      1 )
    call shift1(    ff,                     &
    &            nxdim,   nydim,      1 )
#endif
#endif

    if (myrank == iroot) then
       cdate = chead(50)
       read(cdate, '(i6.6,5i2.2)') idate
       call cyh2ss(  ttt, idate )
       write(jfpar, *) ' start time :', idate
       if (irstrt /= 0) then
          if (ttt /= tstrt) then
             write(jfpar, *) '*** start time error ***'
             call mpi_abort(mpi_comm_ogcm, 1, ierr)
          end if
       end if
    end if

  end subroutine restrt
   
! =====================================================================

  subroutine finout(                &
#ifdef OPT_NPZD
   &         pco2o,     ff,         &
#endif /* OPT_NPZD */
   &            tt, ntstep,         &
   &             t,     ft,     hb, &
   &        orsout, orsrwd)

    use zocdim,   only  :                          &
   &      nxdim,  nydim,  nzdim,  ntdim,     nic,  &
   &     nxgdim, nygdim,  igstr,  jgstr,    kstr,  &
   &        nxg,    nyg,     nz,   nxyg,   nxyzg
    use zocnod,   only  :                          &
   &     myrank,  iroot, mpi_comm_ogcm
    use ufile
    use ucaln
    use mpiio
    
    implicit none

#ifdef OPT_NPZD
    real(8),   intent(inout) ::   pco2o(nxdim, nydim), ff(nxdim, nydim)
#endif 
    real(8),   intent(inout) ::       t(nxdim, nydim, nzdim, ntdim)
    real(8),   intent(inout) ::      ft(nxdim, nydim, ntdim)
    real(8),   intent(in)    ::      hb(nxdim, nydim)
    real(8),   intent(in)    ::      tt
    integer(4),intent(in)    ::  ntstep
    logical,   intent(in)    ::  orsout, orsrwd
!---- local
    integer(4)         ::      i,      l
    integer(4), save   :: istat
    integer(4)         :: ierr
    character(len=ncf) :: crun = '(RUN NAME WAS NOT SET)'
    character(len=8)   :: hdate
    character(len=10)  :: htime
    character(len=5)   :: hzone
    integer            :: ivalues(1:8)
    character(len=16)  :: citem
    logical,    save   :: ofirst = .true.

    namelist /nmfrst/ cfrest
    namelist /nmrun/ crun

    if (ofirst) then
       call rewnml(ifpar, jfpar)
       read(ifpar, nmfrst, iostat=istat)
       call cstnml( jfpar, 'finout', 'nmfrst', istat )
       write(jfpar, nmfrst)
       call rewnml(ifpar, jfpar)
       read(ifpar, nmrun, iostat=istat)
       call cstnml( jfpar, 'finout', 'nmrun', istat )
       call mpi_filopn(mpi_fh_w, cfrest, 'WRITE')
       ofirst = .false.
    end if

    if (crun(1:1) == '(') then
       chrnum = 'offline OECO2'
    else
       chrnum = crun(1:16)
    end if
    
    if (.not. orsout) return

    do i = 1, 64
       write(chead(i), '(16x)')
    end do
    if (orsrwd) then
       dispw = 0
    end if

    call css2yh( idate, tt )
    write(chead(1) , '(i16)' ) 9010
    chead(2) = chrnum
    write(chead(25), '(i16)'  ) nint(tt / 3.6d3)
    chead(26) = 'HOUR'
    chead(38) = 'UR8'
    write(chead(27), '(i4.4,2i2.2,1x,3i2.2,1x)') idate(1:6)
    write(chead(48), '(i4.4,2i2.2,1x,3i2.2,1x)') idate(1:6)
    write(chead(49), '(i4.4,2i2.2,1x,3i2.2,1x)') idate(1:6)
    write(chead(50), '(i6.6,5i2.2)') idate(1:6)
    write(chead(30), '(i16)'  ) 1
    write(chead(31), '(i16)'  ) nxg
    write(chead(33), '(i16)'  ) 1
    write(chead(34), '(i16)'  ) nyg
    write(chead(36), '(i16)'  ) 1
    write(chead(39), '(e16.7)') dundef
    chead(40) = chead(39)
    chead(41) = chead(39)
    chead(42) = chead(39)
    chead(43) = chead(39)
    write(chead(44), '(i16)'  ) 1
    write(chead(46), '(i16)'  ) 0
    write(chead(47), '(e16.7)') 0.d0
    chead(61) = 'offline OECO2'
    chead(63) = 'offline OECO2'
    call date_and_time(hdate, htime, hzone, ivalues)
    write(chead(60), '(i4.4,2i2.2,1x,3i2.2,1x)') ivalues(1:3), ivalues(5:7)
    chead(62) = chead(60)

    !$acc update self(t, ft)
#ifdef OPT_NPZD
    !$acc update self(pco2o, ff)
#endif

    do l = 3, ntdim
       call edhead(trim(ctrnam(l)), '', '', 'OCLVTT')
       call mpi_write_header(chead, mpi_fh_w, dispw)
       call mpi_write_3d(t(1, 1, 1, l), mpi_fh_w, dispw)
    end do
    
    do l = 3, ntdim
       write(citem, '(a,i2.2)') 'TRCFLX', l
       call edhead(citem, '', '', 'OCSFCT')
       call mpi_write_header(chead, mpi_fh_w, dispw)
       call mpi_write_2d(ft(1, 1, l), mpi_fh_w, dispw)
    end do

#ifdef OPT_NPZD
    call edhead('pCO2', '', '', 'OCSFCT')
    call mpi_write_header(chead, mpi_fh_w, dispw)
    call mpi_write_2d(pco2o, mpi_fh_w, dispw)

    call edhead('FF', '', '', 'OCSFCT')
    call mpi_write_header(chead, mpi_fh_w, dispw)
    call mpi_write_2d(ff, mpi_fh_w, dispw)
#endif

    if (myrank == iroot) then
       write(jfpar, *) '*** Write restart file ***'
       write(jfpar, *) ' time :', idate
       write(jfpar, *) ' step :', ntstep
    end if
    call mpi_barrier(mpi_comm_ogcm, ierr)
    
  end subroutine finout

! =====================================================================

  subroutine edhead(                 &
       &             ccitem,         &
       &              htitl,  hunit, &
       &              cclas)
    use zocdim, only : nxg, nyg, nxyg, nxyzg, nic, nz

    implicit none

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

  subroutine print_stats_3d(data, cname, cpos, vmax )

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
    
    real(8),          intent(in)              :: data(:,:,:)
    character(len=*), intent(in)              :: cname
    character(1),     intent(in),    optional :: cpos
    real(8),          intent(inout), optional :: vmax
    
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

    vmax = dmaxg

  end subroutine print_stats_3d

! =====================================================================

  subroutine rstadd(                          &
    &      additm,   oeof,                    &
    &       ixdim,  jydim,  kzdim,            &
#ifndef OPT_TRIPOLE                           
    &      ccitem,   clas )                   
#else                                         
    &      ccitem,   clas,                    &
    &        fact,   ioff,   joff )
#endif
    
    use zocdim,   only  :                     &
    &      nxdim,  nydim,  nzdim,    nic,     &
    &      igstr,  jgstr,   kstr,             &
    &        nxg,    nyg,     nz              
    use zocnod,   only  :                     &
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
    real(8),      intent(in)     ::    fact  
!(IOFF,JOFF)=(0,0)for T-cell, (-1,-1)for V-cell 
    integer(4),   intent(in)     ::    ioff,   joff 
#endif

!---- local variables
    integer(4)   ::     i,     j,     k
    integer(4)   ::  ierr

    if ( clas(1:3) == 'OCN' ) then
       oeof = .false.
       call mpi_read_chead(chead, mpi_fh_r, disp, icread)
       if (icread .ne. 1024) then
          oeof = .true.
       end if
       if (.not. oeof) then
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
       if (icread .ne. 1024) then
          oeof = .true.
       end if
       if (.not. oeof) then
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
       if (icread .ne. 1024) then
          oeof = .true.
       end if
       if (.not. oeof) then
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

  subroutine finadd(                                       &
    &             additm,  ixdim,  jydim,  kzdim,          &
    &             ccitem,   clas  ) 

    use zocdim,   only  :                                  &
    &      igstr,  jgstr,   kstr,                          &
    &        nxg,    nyg,     nz,   nxyg,  nxyzg, nic
    use zocnod,   only  :                                  &
    &     myrank,  iroot
    use mpiio
    implicit none

    integer(4),   intent(in)     ::   ixdim,  jydim,  kzdim
    real(8),      intent(inout)  ::  additm(ixdim, jydim, kzdim)
    character(*), intent(in)     ::  ccitem,   clas

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

end module brstt
