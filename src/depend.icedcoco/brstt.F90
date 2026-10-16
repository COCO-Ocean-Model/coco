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
! ---------------------------------------------------------------------

  use zocfil,   only  :   ncf

  implicit none
#include "mpif.h"
#include "coco.h"
  private

  real(8),    allocatable,  save  ::    buf2(:,:),   buf3(:,:,:),   bufi(:,:,:)
  real(8),    allocatable,  save  ::     g2d(:,:),    g3d(:,:,:),    gid(:,:,:)
  integer(4),               save  ::  nfinit,      nfrest
  integer(4),               save  ::   idate(1:6)
  logical,                  save  ::  ofirst
  character(len=ncf)              ::  cfinit,      cfrest
  character(len=16),        save  ::   chead(64) 
  data ofirst / .true. /
  data cfinit, cfrest / 'not-specified', 'not-specified' /

  character(len=16),        save  ::   chrnum
  real(8), save  ::  dundef = -1.d20

  public  ::  restrt,  rstadd,  finadd, finout

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
    use bgsid
    use bgs2d
    use bgs3d
    use ucaln
    use bshft

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
    real(8)           ::    ttt
    real(8),    save  ::     si
    character(len=16) ::  cdate
    integer(4)        ::      i,      j,      k,      l
    integer(4), save  ::  ifpar,  jfpar,  istat
    integer(4), save  :: irstrt,   ierr
    logical           :: oeof
    
    namelist /nmfini/ cfinit, irstrt
    namelist /nmislt/ si
    data     si / 5.d0 /
    data irstrt / 0 /

    READ_NAMELIST( nmfini )
    READ_NAMELIST( nmislt )

    if ( myrank == iroot ) then
       call filopn(nfinit, cfinit, 'READ')
       write(jfpar, *) '*** read initial condition file ***'
    end if

    if (myrank /= iroot ) then
       allocate ( buf2(1,1) )
       allocate ( buf3(1,1,1) )
       allocate ( bufi(1,1,1) )
       allocate (  g2d(1,1) )
       allocate (  g3d(1,1,1) )
       allocate (  gid(1,1,1) )
    else        
       allocate ( buf2(nxg,nyg) )
       allocate ( buf3(nxg,nyg,nz) )
       allocate ( bufi(nxg,nyg,nic) )
       allocate ( g2d(nxgdim,nygdim) )
       allocate ( g3d(nxgdim,nygdim,nzdim) )
       allocate ( gid(nxgdim,nygdim,0:nic) )
    end if

    call rstadd(ub,          oeof, nxdim, nydim, nzdim, 'U',     'OCN', -1.d0,  -1,  -1, vinit=0.d0 )
    call rstadd(vb,          oeof, nxdim, nydim, nzdim, 'V',     'OCN', -1.d0,  -1,  -1, vinit=0.d0 )
    call rstadd(tb(1,1,1,1), oeof, nxdim, nydim, nzdim, 'T',     'OCN',  1.d0,   0,   0, vinit=0.d0 )
    call rstadd(tb(1,1,1,2), oeof, nxdim, nydim, nzdim, 'S',     'OCN',  1.d0,   0,   0, vinit=0.d0 )
    call rstadd(hb,          oeof, nxdim, nydim,     1, 'SH',    'SFC',  1.d0,   0,   0, vinit=0.d0 )
    call rstadd(ubtb,        oeof, nxdim, nydim,     1, 'UBT',   'SFC', -1.d0,  -1,  -1, vinit=0.d0 )
    call rstadd(vbtb,        oeof, nxdim, nydim,     1, 'VBT',   'SFC', -1.d0,  -1,  -1, vinit=0.d0 )
    call rstadd(w,           oeof, nxdim, nydim, nzdim, 'W',     'OCN',  1.d0,   0,   0, vinit=0.d0 )
    call rstadd(ab,          oeof, nxdim, nydim, nic+1, 'AI',    'ICE',  1.d0,   0,   0, vinit=0.d0 )
    call rstadd(hib,         oeof, nxdim, nydim, nic+1, 'HI',    'ICE',  1.d0,   0,   0, vinit=0.d0 )
    call rstadd(uib,         oeof, nxdim, nydim,     1, 'UI',    'SFC', -1.d0,  -1,  -1, vinit=0.d0 )
    call rstadd(vib,         oeof, nxdim, nydim,     1, 'VI',    'SFC', -1.d0,  -1,  -1, vinit=0.d0 )
    call rstadd(tib,         oeof, nxdim, nydim, nic+1, 'TI',    'ICE',  1.d0,   0,   0, vinit=-0.1d0)
    call rstadd(hsb,         oeof, nxdim, nydim, nic+1, 'HS',    'ICE',  1.d0,   0,   0, vinit=0.d0 )
    call rstadd(ft(1,1,1),   oeof, nxdim, nydim,     1, 'FT',    'SFC',  1.d0,   0,   0, vinit=0.d0 )
    call rstadd(swabs,       oeof, nxdim, nydim,     1, 'SWABS', 'SFC',  1.d0,   0,   0, vinit=0.d0 )
    call rstadd(ft(1,1,2),   oeof, nxdim, nydim,     1, 'FW',    'SFC',  1.d0,   0,   0, vinit=0.d0 )
    call rstadd(fs,          oeof, nxdim, nydim,     1, 'FS',    'SFC',  1.d0,   0,   0, vinit=0.d0 )
    call rstadd(taux,        oeof, nxdim, nydim,     1, 'TAUX',  'SFC', -1.d0,  -1,  -1, vinit=0.d0 )
    call rstadd(tauy,        oeof, nxdim, nydim,     1, 'TAUY',  'SFC', -1.d0,  -1,  -1, vinit=0.d0 )
    call rstadd(amv,         oeof, nxdim, nydim, nzdim, 'AMV',   'OCN',  1.d0,  -1,  -1, vinit=0.d0 )
    call rstadd(ahv,         oeof, nxdim, nydim, nzdim, 'AHV',   'OCN',  1.d0,   0,   0, vinit=0.d0 )
    call rstadd(ptop,        oeof, nxdim, nydim,     1, 'PTOP',  'SFC',  1.d0,   0,   0, vinit=0.d0 )
    call rstadd(tsi,         oeof, nxdim, nydim, nic+1, 'TSI',   'ICE',  1.d0,   0,   0, vinit=dtds * si)

    do l = 3, ntdim
       call rstadd(tb(1,1,1,l), oeof, nxdim, nydim, nzdim, 'TB', 'OCN',  1.d0,   0,   0, vinit=0.d0 )
    end do
    do l = 3, ntdim
           call rstadd(ft(1,1,  l), oeof, nxdim, nydim,     1, 'FT', 'SFC',  1.d0,   0,   0, vinit=0.d0 )
    end do
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
    &      ccitem,   clas,                                            &
    &       fact0,  ioff0,  joff0,  vinit)

    use zocdim,   only  :                                             &
    &      nxdim,  nydim,  nzdim,    nic,                             &
    &      igstr,  jgstr,   kstr,                                     &
    &        nxg,    nyg,     nz
    use zocnod,   only  :                                             &
    &     myrank,  iroot, mpi_comm_ogcm
    use bgs2d
    use bgs3d
    use bgsid
    use bshft

    implicit none

    integer(4),   intent(in)     ::   ixdim,  jydim,  kzdim
    real(8),      intent(inout)  ::  additm(ixdim, jydim, kzdim)
    logical,      intent(inout)  ::    oeof
    character(*), intent(in)     ::  ccitem,   clas

    real(8),      intent(in), optional :: fact0        ! 1 for scalar, -1 for vector
    integer(4),   intent(in), optional :: ioff0, joff0 !((0,0)for T-cell, (-1,-1)for V-cell 
    real(8),      intent(in), optional :: vinit

!---- local variables
    integer(4)   :: i, j, k, ierr, ioff, joff
    real(8)      :: fact

    if(present(fact0)) then
       fact=fact0
       ioff=ioff0
       joff=joff0
    end if
    if(present(vinit)) then
       additm=vinit
    end if

    if ( clas(1:3) == 'OCN' ) then
       if ( myrank == iroot ) then
          oeof = .true.
          read(nfinit, end=309) chead
          read(nfinit) buf3
          do k = 1, nz
             do j = 1, nyg
                do i = 1, nxg
                   g3d(igstr+i-1, jgstr+j-1, kstr+k-1) = buf3(i, j, k)
                end do
             end do
          end do
          oeof = .false.
309       continue
       end if
       call mpi_bcast(oeof, 1, mpi_logical, iroot, mpi_comm_ogcm, ierr)
       if ( .not. oeof ) then
          call scatter_3d(additm, g3d)
          additm(:,:,:kstr-1)=0.d0
          additm(:,:,kstr+nz:)=0.d0
          call shift1(additm, nxdim, nydim, nzdim, fact,  ioff,  joff  )
       end if
    else if ( clas(1:3) == 'ICE' ) then
       if ( myrank == iroot ) then
          oeof = .true.
          read(nfinit, end=314) chead
          read(nfinit) bufi
          do k = 1, nic
             do j = 1, nyg
                do i = 1, nxg
                   gid(igstr+i-1, jgstr+j-1, k) = bufi(i, j, k)
                end do
             end do
          end do
          oeof = .false.
314       continue
       end if
       call mpi_bcast(oeof, 1, mpi_logical, iroot, mpi_comm_ogcm, ierr)
       if ( .not. oeof ) then
          call scatter_id(additm, gid)
          call shift1(additm, nxdim, nydim, nic+1, fact,  ioff,  joff  )
       end if
    else
       if ( myrank == iroot ) then
          oeof = .true.
          read(nfinit, end=319) chead
          read(nfinit) buf2
          do j = 1, nyg
             do i = 1, nxg
                g2d(igstr+i-1, jgstr+j-1) = buf2(i, j)
             end do
          end do
          oeof = .false.
319       continue
       end if
       call mpi_bcast(oeof, 1, mpi_logical, iroot, mpi_comm_ogcm, ierr)
       if ( .not. oeof ) then
          call scatter_2d(additm, g2d)
          call shift1(additm, nxdim, nydim, 1, fact,  ioff,  joff)
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
    use bgsid
    use bgs2d
    use bgs3d
    use ucaln

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

    namelist /nmfrst/ cfrest
    namelist /nmrun/ crun

    if ( ofirst ) then
       READ_NAMELIST( nmfrst )
       READ_NAMELIST( nmrun  )
       if ( myrank == iroot ) then
          call filopn(nfrest, cfrest, 'WRITE')
       end if
       ofirst = .false.
    end if

    if (crun(1:1) == '(') then
       chrnum = 'COCO stand-alone'
    else
       chrnum = crun(1:16)
    end if

    if ( .not. orsout ) return

    if ( myrank == iroot ) then
       do i = 1, 64
          write(chead(i),'(16x)')
       end do
       if (orsrwd) then
          rewind(nfrest)
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
    end if

    call finadd(ub,          nxdim, nydim, nzdim, 'UO',   'OCN', 'ocean zonal velocity',                 'cm/s',       'OCLVTV')
    call finadd(vb,          nxdim, nydim, nzdim, 'VO',   'OCN', 'ocean meridional velocity',            'cm/s',       'OCLVTV')
    call finadd(tb,          nxdim, nydim, nzdim, 'TO',   'OCN', 'ocean temperature',                    'degC',       'OCLVTT')
    call finadd(tb(1,1,1,2), nxdim, nydim, nzdim, 'SO',   'OCN', 'ocean salinity',                       'psu',        'OCLVTT')
    call finadd(hb,          nxdim, nydim,     1, 'SHO',  'SFC', 'sea surface height',                   'cm',         'OCSFCT')
    call finadd(ubtb,        nxdim, nydim,     1, 'UBTO', 'SFC', 'ocean zonal transport',                'cm^2/s',     'OCSFCV')
    call finadd(vbtb,        nxdim, nydim,     1, 'VBTO', 'SFC', 'ocean meridional transport',           'cm^2/s',     'OCSFCV')
    call finadd(w,           nxdim, nydim, nzdim, 'WO',   'OCN', 'ocean vertical velocity',              'cm/s',       'OCLVMT')
    call finadd(ab,          nxdim, nydim, nic+1, 'AI',   'ICE', 'ice concentration',                    'N.D.',       'OCICET')
    call finadd(hib,         nxdim, nydim, nic+1, 'HI',   'ICE', 'ice thickness',                        'cm',         'OCICET')
    call finadd(uib,         nxdim, nydim,     1, 'UI',   'SFC', 'ice zonal velocity',                   'cm/s',       'OCSFCV')
    call finadd(vib,         nxdim, nydim,     1, 'VI',   'SFC', 'ice meridional velocity',              'cm/s',       'OCSFCV')
    call finadd(tib,         nxdim, nydim, nic+1, 'TI',   'ICE', 'ice temperature',                      'degC',       'OCICET')
    call finadd(hsb,         nxdim, nydim, nic+1, 'HS',   'ICE', 'snow thickness',                       'cm',         'OCICET')
    call finadd(ft,          nxdim, nydim,     1, 'FT',   'SFC', 'sea surface temperature flux',         'K cm/s',     'OCSFCT')
    call finadd(swabs,       nxdim, nydim,     1, 'SWABS','SFC', 'absorbed shortwave',                   'erg/cm^2/s', 'OCSFCT')
    call finadd(ft(1,1,2),   nxdim, nydim,     1, 'FW',   'SFC', 'sea surface freshwater flux',          'cm/s',       'OCSFCT')
    call finadd(fs,          nxdim, nydim,     1, 'FS',   'SFC', 'sea surface salinity flux',            'psu cm/s',   'OCSFCT')
    call finadd(taux,        nxdim, nydim,     1, 'TAUX', 'SFC', 'zonal wind stress',                    'dyn/cm^2',   'OCSFCV')
    call finadd(tauy,        nxdim, nydim,     1, 'TAUY', 'SFC', 'meridional wind stress',               'dyn/cm^2',   'OCSFCV')
    call finadd(amv,         nxdim, nydim, nzdim, 'AMV',  'OCN', 'ocean vertical viscosity',             'cm^2/s',     'OCLVMV')
    call finadd(ahv,         nxdim, nydim, nzdim, 'AHV',  'OCN', 'ocean vertical diffusivity',           'cm^2/s',     'OCLVMT')
    call finadd(ptop,        nxdim, nydim,     1, 'PTOP', 'SFC', 'sea surface presssure (inc. sea ice)', 'dyn/cm^2',   'OCSFCT')
    call finadd(tsi,         nxdim, nydim, nic+1, 'TSI',  'ICE', 'sea ice surface temperature',          'degC',       'OCICET')
    
    do l = 3, ntdim
       write(citem, '(a,i2.2)') 'TRACER', l
       call finadd(tb(1,1,1,l), nxdim, nydim, nzdim, citem, 'OCN', '', '', 'OCLVTT')
    end do
      
    do l = 3, ntdim
       write(citem, '(a,i2.2)') 'TRCFLX', l
       call finadd(ft(1,1,l),   nxdim, nydim,     1, citem, 'SFC', '', '', 'OCSFCT')
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
    &             ccitem,   clas,  htitl,  hunit,  cclas) 

    use zocdim,   only  :                                             &
    &      igstr,  jgstr,   kstr,                                     &
    &        nxg,    nyg,     nz,   nxyg,  nxyzg, nic
    use zocnod,   only  :                                             &
    &     myrank,  iroot
    use bgs2d
    use bgs3d
    use bgsid

    implicit none

    integer(4),   intent(in)     ::   ixdim,  jydim,  kzdim
    real(8),      intent(inout)  ::  additm(ixdim, jydim, kzdim)
    character(*), intent(in)     ::  ccitem,   clas
    character(*), intent(in), optional ::  htitl,  hunit, cclas

!---- local variables
    integer(4)   ::     i,     j,     k
    

    if(present(htitl)) then
       call edhead(ccitem, htitl, hunit, cclas)
    else
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
       else if ( clas(1:3) == 'ICE' ) then
          write(chead(35), '(a)') 'NUMBER1000'
          write(chead(37), '(i16)') nic
          write(chead(64), '(i16)') nxyg*nic
       else
          write(chead(35), '(a)') 'SFC1'
          write(chead(37), '(i16)') 1
          write(chead(64), '(i16)') nxyg
       end if
    end if

    if ( clas(1:3) == 'OCN' ) then
       call gather_3d(g3d, additm)
       if ( myrank == iroot ) then
          do k = 1, nz
             do j = 1, nyg
                do i = 1, nxg
                   buf3(i, j, k) = g3d(igstr+i-1, jgstr+j-1, kstr+k-1)
                end do
             end do
          end do
          write(nfrest) chead
          write(nfrest) buf3
       end if
    else if ( clas(1:3) == 'ICE' ) then
       call gather_id(gid, additm)
       if ( myrank == iroot ) then
          do k = 1, nic
             do j = 1, nyg
                do i = 1, nxg
                   bufi(i, j, k) = gid(igstr+i-1, jgstr+j-1, k)
                end do
             end do
          end do
          write(nfrest) chead
          write(nfrest) bufi
       end if
    else
       call gather_2d(g2d, additm)
       if ( myrank == iroot ) then
          do j = 1, nyg
             do i = 1, nxg
                buf2(i, j) = g2d(igstr+i-1, jgstr+j-1)
             end do
          end do
          write(nfrest) chead
          write(nfrest) buf2
       end if
    end if
  end subroutine finadd

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

end module brstt
