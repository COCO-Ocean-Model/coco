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

    call rewnml(ifpar, jfpar)
    read(ifpar, nmfini, iostat=istat)
    call cstnml( jfpar, 'restrt', 'nfini', istat )
    write(jfpar, nmfini)
    
    call rewnml(ifpar, jfpar)
    read(ifpar, nmislt, iostat=istat)
    call cstnml( jfpar, 'restrt', 'nmislt', istat )
    write(jfpar, nmislt)

    if ( myrank == iroot ) then
       call filopn(nfinit, cfinit, 'READ')
       write(jfpar, *) '*** read initial condition file ***'
    end if

    tb    = 0.d0
    ft    = 0.d0
    swabs = 0.d0
    fs    = 0.d0
    taux  = 0.d0
    tauy  = 0.d0
    amv   = 0.d0
    ahv   = 0.d0
    ptop  = 0.d0
    tsi   = dtds * si
    hib   = 0.d0
    hsb   = 0.d0
    tib   = -0.1d0

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

    if ( myrank == iroot ) then
       rewind(nfinit)
       read(nfinit, end=109) chead
       read(nfinit) buf3
       do k = 1, nz
          do j = 1, nyg
             do i = 1, nxg
                g3d(igstr+i-1, jgstr+j-1, kstr+k-1) = buf3(i, j, k)
             end do
          end do
       end do
109    continue
    end if
    call scatter_3d(ub, g3d)

    if ( myrank == iroot ) then
       read(nfinit, end=119) chead
       read(nfinit) buf3
       do k = 1, nz
          do j = 1, nyg
             do i = 1, nxg
                g3d(igstr+i-1, jgstr+j-1, kstr+k-1) = buf3(i, j, k)
             end do
          end do
       end do
119    continue
    end if
    call scatter_3d(vb, g3d)
    
    if ( myrank == iroot ) then
       read(nfinit, end=129) chead
       read(nfinit) buf3
       do k = 1, nz
          do j = 1, nyg
             do i = 1, nxg
                g3d(igstr+i-1, jgstr+j-1, kstr+k-1) = buf3(i, j, k)
             end do
          end do
       end do
 129     continue
    end if
    call scatter_3d(tb, g3d)
      
    if ( myrank == iroot ) then
       read(nfinit, end=139) chead
       read(nfinit) buf3
       do k = 1, nz
          do j = 1, nyg
             do i = 1, nxg
                g3d(igstr+i-1, jgstr+j-1, kstr+k-1) = abs( buf3(i, j, k) )
             end do
          end do
       end do
139    continue
    end if
    call scatter_3d(tb(1, 1, 1, 2), g3d)

    if ( myrank == iroot ) then
       read(nfinit, end=141) chead
       read(nfinit) buf2
       do j = 1, nyg
          do i = 1, nxg
             g2d(igstr+i-1, jgstr+j-1) = buf2(i, j)
          end do
       end do
141    continue
    end if
    call scatter_2d(hb, g2d)
    
    if ( myrank == iroot ) then
       read(nfinit, end=143) chead
       read(nfinit) buf2
       do j = 1, nyg
          do i = 1, nxg
             g2d(igstr+i-1, jgstr+j-1) = buf2(i, j)
          end do
       end do
143    continue
    end if
    call scatter_2d(ubtb, g2d)

    if ( myrank == iroot ) then
       read(nfinit, end=145) chead
       read(nfinit) buf2
       do j = 1, nyg
          do i = 1, nxg
             g2d(igstr+i-1, jgstr+j-1) = buf2(i, j)
          end do
       end do
145    continue
    end if
    call scatter_2d(vbtb, g2d)

    if ( myrank == iroot ) then
       read(nfinit, end=159) chead
       read(nfinit) buf3
       do k = 1, nz
          do j = 1, nyg
             do i = 1, nxg
                g3d(igstr+i-1, jgstr+j-1, kstr+k-1) = buf3(i, j, k)
             end do
          end do
       end do
159    continue
    end if
    call scatter_3d(w, g3d)

    if ( myrank == iroot ) then
       read(nfinit, end=301) chead
       read(nfinit) bufi
       do k = 1, nic
          do j = 1, nyg
             do i = 1, nxg
                gid(igstr+i-1, jgstr+j-1, k) = bufi(i, j, k)
             end do
          end do
       end do
301    continue
    end if
    call scatter_id(ab, gid)
    
    if ( myrank == iroot ) then
       read(nfinit, end=311) chead
       read(nfinit) bufi
       do k = 1, nic
          do j = 1, nyg
             do i = 1, nxg
                gid(igstr+i-1, jgstr+j-1, k) = bufi(i, j, k)
             end do
          end do
       end do
311    continue
    end if
    call scatter_id(hib, gid)

    if ( myrank == iroot ) then
       read(nfinit, end=321) chead
       read(nfinit) buf2
       do j = 1, nyg
          do i = 1, nxg
             g2d(igstr+i-1, jgstr+j-1) = buf2(i, j)
          end do
       end do
321    continue
    end if
    call scatter_2d(uib, g2d)

    if ( myrank == iroot ) then
       read(nfinit, end=331) chead
       read(nfinit) buf2
       do j = 1, nyg
          do i = 1, nxg
             g2d(igstr+i-1, jgstr+j-1) = buf2(i, j)
          end do
       end do
331    continue
    end if
    call scatter_2d(vib, g2d)

    if ( myrank == iroot ) then
       read(nfinit, end=336) chead
       read(nfinit) bufi
       do k = 1, nic
          do j = 1, nyg
             do i = 1, nxg
                gid(igstr+i-1, jgstr+j-1, k) = bufi(i, j, k)
             end do
          end do
       end do
336    continue
    end if
    call scatter_id(tib, gid)

    if ( myrank == iroot ) then
       read(nfinit, end=341) chead
       read(nfinit) bufi
       do k = 1, nic
          do j = 1, nyg
             do i = 1, nxg
                gid(igstr+i-1, jgstr+j-1, k) = bufi(i, j, k)
             end do
          end do
       end do
341    continue
    end if
    call scatter_id(hsb, gid)
    
    if ( myrank == iroot ) then
       read(nfinit, end=189) chead
       read(nfinit) buf2
       do j = 1, nyg
          do i = 1, nxg
             g2d(igstr+i-1, jgstr+j-1) = buf2(i, j)
          end do
       end do
189    continue
    end if
    call scatter_2d(ft, g2d)
    
    if ( myrank == iroot ) then
       read(nfinit, end=401) chead
       read(nfinit) buf2
       do j = 1, nyg
          do i = 1, nxg
             g2d(igstr+i-1, jgstr+j-1) = buf2(i, j)
          end do
       end do
401    continue
    end if
    call scatter_2d(swabs, g2d)

    if ( myrank == iroot ) then
       read(nfinit, end=199) chead
       read(nfinit) buf2
       do j = 1, nyg
          do i = 1, nxg
             g2d(igstr+i-1, jgstr+j-1) = buf2(i, j)
          end do
       end do
199    continue
    end if
    call scatter_2d(ft(1, 1, 2), g2d)
    
    if ( myrank == iroot ) then
       read(nfinit, end=411) chead
       read(nfinit) buf2
       do j = 1, nyg
          do i = 1, nxg
             g2d(igstr+i-1, jgstr+j-1) = buf2(i, j)
          end do
       end do
411    continue
    end if
    call scatter_2d(fs, g2d)

    if ( myrank == iroot ) then
       read(nfinit, end=421) chead
       read(nfinit) buf2
       do j = 1, nyg
          do i = 1, nxg
             g2d(igstr+i-1, jgstr+j-1) = buf2(i, j)
          end do
       end do
421    continue
    end if
    call scatter_2d(taux, g2d)

    if ( myrank == iroot ) then
       read(nfinit, end=431) chead
       read(nfinit) buf2
       do j = 1, nyg
          do i = 1, nxg
             g2d(igstr+i-1, jgstr+j-1) = buf2(i, j)
          end do
       end do
 431     continue
    end if
    call scatter_2d(tauy, g2d)

    if ( myrank == iroot ) then
       read(nfinit, end=169) chead
       read(nfinit) buf3
       do k = 1, nz
          do j = 1, nyg
             do i = 1, nxg
                g3d(igstr+i-1, jgstr+j-1, kstr+k-1) = buf3(i, j, k)
             end do
          end do
       end do
169    continue
    end if
    call scatter_3d(amv, g3d)

    if ( myrank == iroot ) then
       read(nfinit, end=179) chead
       read(nfinit) buf3
       do k = 1, nz
          do j = 1, nyg
             do i = 1, nxg
                g3d(igstr+i-1, jgstr+j-1, kstr+k-1) = buf3(i, j, k)
             end do
          end do
       end do
179    continue
    end if
    call scatter_3d(ahv, g3d)

    if ( myrank == iroot ) then
       read(nfinit, end=441) chead
       read(nfinit) buf2
       do j = 1, nyg
          do i = 1, nxg
             g2d(igstr+i-1, jgstr+j-1) = buf2(i, j)
          end do
       end do
441    continue
    end if
    call scatter_2d(ptop, g2d)

    if ( myrank == iroot ) then
       read(nfinit, end=451) chead
       read(nfinit) bufi
       do k = 1, nic
          do j = 1, nyg
             do i = 1, nxg
                gid(igstr+i-1, jgstr+j-1, k) = bufi(i, j, k)
             end do
          end do
       end do
451    continue
    end if
    call scatter_id(tsi, gid)

    do l = 3, ntdim
       if ( myrank == iroot ) then
          read(nfinit, end=209) chead
          read(nfinit) buf3
          do k = 1, nz
             do j = 1, nyg
                do i = 1, nxg
                   g3d(igstr+i-1, jgstr+j-1, kstr+k-1) = buf3(i, j, k)
                end do
             end do
          end do
209       continue
       end if
       call scatter_3d(tb(1, 1, 1, l), g3d)
    end do

    do l = 3, ntdim
       if ( myrank == iroot ) then
          read(nfinit, end=219) chead
          read(nfinit) buf2
          do j = 1, nyg
             do i = 1, nxg
                g2d(igstr+i-1, jgstr+j-1) = buf2(i, j)
             end do
          end do
219       continue
       end if
       call scatter_2d(ft(1, 1, l), g2d)
    end do

#ifdef OPT_TRIPOLE
    call shift2(    ub,      vb,                                      &
    &             nxdim,   nydim,  nzdim,                             &
    &             -1.d0,      -1,     -1 ) 
    call shift1(   tb,                                                &
    &             nxdim,   nydim, nztdim,                             &
    &              1.d0,       0,      0)
    call shift1(   amv,                                               &
    &            nxdim,   nydim,  nzdim,                              &
    &             1.d0,      -1,     -1 )
    call shift2(     w,     ahv,                                      &
    &            nxdim,   nydim,  nzdim,                              &
    &             1.d0,       0,      0 )
    call shift1(    hb,                                               &
    &            nxdim,   nydim,      1,                              &
    &             1.d0,       0,      0 )
    call shift2(  ubtb,    vbtb,                                      &
    &            nxdim,   nydim,      1,                              &
    &            -1.d0,      -1,     -1 )
    call shift3(    ab,     hib,    hsb,                              &
    &            nxdim,   nydim,  nic+1,                              &
    &             1.d0,       0,      0 )
    call shift2(   tsi,     tib,                                      &
    &            nxdim,   nydim,  nic+1,                              &
    &             1.d0,       0,      0 )
    call shift2(   uib,     vib,                                      &
    &            nxdim,   nydim,      1,                              &
    &            -1.d0,      -1,     -1 )
    call shift1(    ft,                                               &
    &            nxdim,   nydim,  ntdim,                              &
    &             1.d0,       0,      0 )
    call shift3(    fs,   swabs,   ptop,                              &
    &            nxdim,   nydim,      1,                              &
    &             1.d0,       0,      0 )
    call shift2(  taux,    tauy,                                      &
    &            nxdim,   nydim,      1,                              &
    &            -1.d0,      -1,     -1 )
#else
    call shift3(    ub,     vb,      w,  nxdim,  nydim,  nzdim )
    call shift1(    tb,                  nxdim,  nydim, nztdim )
    call shift2(   amv,    ahv,          nxdim,  nydim,  nzdim )
    call shift3(    hb,   ubtb,   vbtb,  nxdim,  nydim,      1 )
    call shift3(    ab,    hib,    hsb,  nxdim,  nydim,  nic+1 )
    call shift2(   tsi,    tib,          nxdim,  nydim,  nic+1 )
    call shift2(   uib,    vib,          nxdim,  nydim,      1 )
    call shift1(    ft,                  nxdim,  nydim,  ntdim )
    call shift2( swabs,     fs,          nxdim,  nydim,      1 )
    call shift3(  taux,   tauy,   ptop,  nxdim,  nydim,      1 )
#endif

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
#ifdef OPT_TRIPOLE
! set 1 for scalar, set -1 for vector
    real(8),      intent(in)    ::    fact  
!(IOFF,JOFF)=(0,0)for T-cell, (-1,-1)for V-cell 
    integer(4),   intent(in)    ::    ioff,   joff 
#endif

!---- local variables
    integer(4)   ::     i,     j,     k
    integer(4)   ::  ierr

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
       call mpi_bcast                                                 &
    &        (oeof, 1, mpi_logical, iroot, mpi_comm_ogcm, ierr)
       if ( .not. oeof ) then
          call scatter_3d(additm, g3d)
#ifdef OPT_TRIPOLE
          call shift1(additm, nxdim, nydim, nzdim,                    &
    &                          fact,  ioff,  joff  )
#else
          call shift1(additm, nxdim, nydim, nzdim)
#endif
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
       call mpi_bcast                                                 &
    &        (oeof, 1, mpi_logical, iroot, mpi_comm_ogcm, ierr)
       if ( .not. oeof ) then
          call scatter_id(additm, gid)
#ifdef OPT_TRIPOLE
          call shift1(additm, nxdim, nydim, nic+1,                    &
    &                          fact,  ioff,  joff  )
#else
          call shift1(additm, nxdim, nydim, nic+1)
#endif
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
       call mpi_bcast                                                 &
    &        (oeof, 1, mpi_logical, iroot, mpi_comm_ogcm, ierr)
       if ( .not. oeof ) then
          call scatter_2d(additm, g2d)
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
       call rewnml(ifpar, jfpar)
       read(ifpar, nmfrst, iostat=istat)
       call cstnml( jfpar, 'finout', 'nmfrst', istat )
       write(jfpar, nmfrst)
       call rewnml(ifpar, jfpar)
       read(ifpar, nmrun, iostat=istat)
       call cstnml( jfpar, 'finout', 'nmrun', istat )
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

    call gather_3d(g3d, ub)
    if ( myrank == iroot ) then
       do k = 1, nz
          do j = 1, nyg
             do i = 1, nxg
                buf3(i, j, k) = g3d(igstr+i-1, jgstr+j-1, kstr+k-1)
             end do
          end do
       end do
       call edhead('UO', 'ocean zonal velocity', 'cm/s', 'OCLVTV')
       write(nfrest) chead
       write(nfrest) buf3
    end if
    
    call gather_3d(g3d, vb)
    if ( myrank == iroot ) then
       do k = 1, nz
          do j = 1, nyg
             do i = 1, nxg
                buf3(i, j, k) = g3d(igstr+i-1, jgstr+j-1, kstr+k-1)
             end do
          end do
       end do
       call edhead('VO', 'ocean meridional velocity', 'cm/s', 'OCLVTV')
       write(nfrest) chead
       write(nfrest) buf3
    end if

    call gather_3d(g3d, tb)
    if ( myrank == iroot ) then
       do k = 1, nz
          do j = 1, nyg
             do i = 1, nxg
                buf3(i, j, k) = g3d(igstr+i-1, jgstr+j-1, kstr+k-1)
             end do
          end do
       end do
       call edhead('TO', 'ocean temperature', 'degC', 'OCLVTT')
       write(nfrest) chead
       write(nfrest) buf3
    end if
    
    call gather_3d(g3d, tb(1, 1, 1, 2))
    if ( myrank == iroot ) then
       do k = 1, nz
          do j = 1, nyg
             do i = 1, nxg
                buf3(i, j, k) = g3d(igstr+i-1, jgstr+j-1, kstr+k-1)
             end do
          end do
       end do
       call edhead('SO', 'ocean salinity', 'psu', 'OCLVTT')
       write(nfrest) chead
       write(nfrest) buf3
    end if
    
    call gather_2d(g2d, hb)
    if ( myrank == iroot ) then
       do j = 1, nyg
          do i = 1, nxg
             buf2(i, j) = g2d(igstr+i-1, jgstr+j-1)
          end do
       end do
       call edhead('SHO', 'sea surface height', 'cm', 'OCSFCT' )
       write(nfrest) chead
       write(nfrest) buf2
    end if

    call gather_2d(g2d, ubtb)
    if ( myrank == iroot ) then
       do j = 1, nyg
          do i = 1, nxg
             buf2(i, j) = g2d(igstr+i-1, jgstr+j-1)
          end do
       end do
       call edhead('UBTO', 'ocean zonal transport', 'cm^2/s', 'OCSFCV')
       write(nfrest) chead
       write(nfrest) buf2
    end if
    
    call gather_2d(g2d, vbtb)
    if ( myrank == iroot ) then
       do j = 1, nyg
          do i = 1, nxg
             buf2(i, j) = g2d(igstr+i-1, jgstr+j-1)
          end do
       end do
       call edhead('VBTO', 'ocean meridional transport', 'cm^2/s', 'OCSFCV')
       write(nfrest) chead
       write(nfrest) buf2
    end if
    
    call gather_3d(g3d, w)
    if ( myrank == iroot ) then
       do k = 1, nz
          do j = 1, nyg
             do i = 1, nxg
                buf3(i, j, k) = g3d(igstr+i-1, jgstr+j-1, kstr+k-1)
             end do
          end do
       end do
       call edhead('WO', 'ocean vertical velocity', 'cm/s', 'OCLVMT')
       write(nfrest) chead
       write(nfrest) buf3
    end if
    
    call gather_id(gid, ab)
    if ( myrank == iroot ) then
       do k = 1, nic
          do j = 1, nyg
             do i = 1, nxg
                bufi(i, j, k) = gid(igstr+i-1, jgstr+j-1, k)
             end do
          end do
       end do
       call edhead('AI', 'ice concentration', 'N.D.', 'OCICET')
       write(nfrest) chead
       write(nfrest) bufi
    end if
    
    call gather_id(gid, hib)
    if ( myrank == iroot ) then
       do k = 1, nic
          do j = 1, nyg
             do i = 1, nxg
                bufi(i, j, k) = gid(igstr+i-1, jgstr+j-1, k)
             end do
          end do
       end do
       call edhead('HI', 'ice thickness', 'cm', 'OCICET')
       write(nfrest) chead
       write(nfrest) bufi
    end if
    
    call gather_2d(g2d, uib)
    if ( myrank == iroot ) then
       do j = 1, nyg
          do i = 1, nxg
             buf2(i, j) = g2d(igstr+i-1, jgstr+j-1)
          end do
       end do
       call edhead('UI', 'ice zonal velocity', 'cm/s', 'OCSFCV')
       write(nfrest) chead
       write(nfrest) buf2
    end if
    
    call gather_2d(g2d, vib)
    if ( myrank == iroot ) then
       do j = 1, nyg
          do i = 1, nxg
             buf2(i, j) = g2d(igstr+i-1, jgstr+j-1)
          end do
       end do
       call edhead('VI', 'ice meridional velocity', 'cm/s', 'OCSFCV')
       write(nfrest) chead
       write(nfrest) buf2
    end if
    
    call gather_id(gid, tib)
    if ( myrank == iroot ) then
       do k = 1, nic
          do j = 1, nyg
             do i = 1, nxg
                bufi(i, j, k) = gid(igstr+i-1, jgstr+j-1, k)
             end do
          end do
       end do
       call edhead('TI', 'ice temperature', 'degC', 'OCICET')
       write(nfrest) chead
       write(nfrest) bufi
    end if
    
    call gather_id(gid, hsb)
    if ( myrank == iroot ) then
       do k = 1, nic
          do j = 1, nyg
             do i = 1, nxg
                bufi(i, j, k) = gid(igstr+i-1, jgstr+j-1, k)
             end do
          end do
       end do
       call edhead('HS', 'snow thickness', 'cm', 'OCICET')
       write(nfrest) chead
       write(nfrest) bufi
    end if
    
    call gather_2d(g2d, ft)
    if ( myrank == iroot ) then
       do j = 1, nyg
          do i = 1, nxg
             buf2(i, j) = g2d(igstr+i-1, jgstr+j-1)
          end do
       end do
       call edhead('FT', 'sea surface temperature flux', 'K cm/s', 'OCSFCT')
       write(nfrest) chead
       write(nfrest) buf2
    end if
    
    call gather_2d(g2d, swabs)
    if ( myrank == iroot ) then
       do j = 1, nyg
          do i = 1, nxg
             buf2(i, j) = g2d(igstr+i-1, jgstr+j-1)
          end do
       end do
       call edhead('SWABS', 'absorbed shortwave', 'erg/cm^2/s', 'OCSFCT')
       write(nfrest) chead
       write(nfrest) buf2
    end if
    
    call gather_2d(g2d, ft(1, 1, 2))
    if ( myrank == iroot ) then
       do j = 1, nyg
          do i = 1, nxg
             buf2(i, j) = g2d(igstr+i-1, jgstr+j-1)
          end do
       end do
       call edhead('FW', 'sea surface freshwater flux', 'cm/s', 'OCSFCT')
       write(nfrest) chead
       write(nfrest) buf2
    end if
    
    call gather_2d(g2d, fs)
    if ( myrank == iroot ) then
       do j = 1, nyg
          do i = 1, nxg
             buf2(i, j) = g2d(igstr+i-1, jgstr+j-1)
          end do
       end do
       call edhead('FS', 'sea surface salinity flux', 'psu cm/s', 'OCSFCT')
       write(nfrest) chead
       write(nfrest) buf2
    end if

    call gather_2d(g2d, taux)
    if ( myrank == iroot ) then
       do j = 1, nyg
          do i = 1, nxg
             buf2(i, j) = g2d(igstr+i-1, jgstr+j-1)
          end do
       end do
       call edhead('TAUX', 'zonal wind stress', 'dyn/cm^2', 'OCSFCV')
       write(nfrest) chead
       write(nfrest) buf2
    end if

    call gather_2d(g2d, tauy)
    if ( myrank == iroot ) then
       do j = 1, nyg
          do i = 1, nxg
             buf2(i, j) = g2d(igstr+i-1, jgstr+j-1)
          end do
       end do
       call edhead('TAUY', 'meridional wind stress', 'dyn/cm^2', 'OCSFCV')
       write(nfrest) chead
       write(nfrest) buf2
    end if
    
    call gather_3d(g3d, amv)
    if ( myrank == iroot ) then
       do k = 1, nz
          do j = 1, nyg
             do i = 1, nxg
                buf3(i, j, k) = g3d(igstr+i-1, jgstr+j-1, kstr+k-1)
             end do
          end do
       end do
       call edhead('AMV', 'ocean vertical viscosity', 'cm^2/s', 'OCLVMV')
       write(nfrest) chead
       write(nfrest) buf3
    end if
    
    call gather_3d(g3d, ahv)
    if ( myrank == iroot ) then
       do k = 1, nz
          do j = 1, nyg
             do i = 1, nxg
                buf3(i, j, k) = g3d(igstr+i-1, jgstr+j-1, kstr+k-1)
             end do
          end do
       end do
       call edhead('AHV', 'ocean vertical diffusivity', 'cm^2/s', 'OCLVMT')
       write(nfrest) chead
       write(nfrest) buf3
    end if
    
    call gather_2d(g2d, ptop)
    if ( myrank == iroot ) then
       do j = 1, nyg
          do i = 1, nxg
             buf2(i, j) = g2d(igstr+i-1, jgstr+j-1)
          end do
       end do
       call edhead('PTOP', 'sea surface presssure (inc. sea ice)', 'dyn/cm^2', 'OCSFCT')
       write(nfrest) chead
       write(nfrest) buf2
    end if
    
    call gather_id(gid, tsi)
    if ( myrank == iroot ) then
       do k = 1, nic
          do j = 1, nyg
             do i = 1, nxg
                bufi(i, j, k) = gid(igstr+i-1, jgstr+j-1, k)
             end do
          end do
       end do
       call edhead('TSI', 'sea ice surface temperature', 'degC', 'OCICET')
       write(nfrest) chead
       write(nfrest) bufi
    end if
    
    do l = 3, ntdim
       call gather_3d(g3d, tb(1, 1, 1, l))
       if ( myrank == iroot ) then
          do k = 1, nz
             do j = 1, nyg
                do i = 1, nxg
                   buf3(i, j, k) = g3d(igstr+i-1, jgstr+j-1, kstr+k-1)
                end do
             end do
          end do
          write(citem, '(a,i2.2)') 'TRACER', l
          call edhead(citem, '', '', 'OCLVTT')
          write(nfrest) chead
          write(nfrest) buf3
       end if
    end do
      
    do l = 3, ntdim
       call gather_2d(g2d, ft(1, 1, l))
       if ( myrank == iroot ) then
          do j = 1, nyg
             do i = 1, nxg
                buf2(i, j) = g2d(igstr+i-1, jgstr+j-1)
             end do
          end do
          write(citem, '(a,i2.2)') 'TRCFLX', l
          call edhead(citem, '', '', 'OCSFCT')
          write(nfrest) chead
          write(nfrest) buf2
       end if
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
    use bgs2d
    use bgs3d
    use bgsid

    implicit none

    integer(4),   intent(in)     ::   ixdim,  jydim,  kzdim
    real(8),      intent(inout)  ::  additm(ixdim, jydim, kzdim)
    character(*), intent(in)     ::  ccitem,   clas

!---- local variables
    integer(4)   ::     i,     j,     k
    
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
       call gather_3d(g3d, additm)
       if ( myrank == iroot ) then
          do k = 1, nz
             do j = 1, nyg
                do i = 1, nxg
                   buf3(i, j, k) = g3d(igstr+i-1, jgstr+j-1, kstr+k-1)
                end do
             end do
          end do
          write(chead(35), '(a,i0)') 'OCDEPT', nz
          write(chead(37), '(i16)') nz
          write(chead(64), '(i16)') nxyzg
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
          write(chead(35), '(a)') 'NUMBER1000'
          write(chead(37), '(i16)') nic
          write(chead(64), '(i16)') nxyg*nic
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
          write(chead(35), '(a)') 'SFC1'
          write(chead(37), '(i16)') 1
          write(chead(64), '(i16)') nxyg
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
