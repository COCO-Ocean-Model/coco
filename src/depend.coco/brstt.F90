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
!     '12.09.20  H.Tatebe: for COCO5.0 in F90
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
  character(len=16)               ::   chead(64) 
  data ofirst / .true. /
  data cfinit, cfrest / 'not-specified', 'not-specified' /

  public  ::  restrt,  rstadd,  finadd, finout

contains

  subroutine restrt(                                                  &
    &       tstrt,                                                    &
    &          ub,     vb,     tb,                                    &
    &          hb,   ubtb,   vbtb,                                    &
    &           w,    amv,    ahv,                                    & 
    &          ft )

    use zocdim,   only  :                                             &
    &      nxdim,  nydim,  nzdim,  ntdim,  nztdim,                    &
    &     nxgdim, nygdim,  igstr,  jgstr,    kstr,                    &
    &        nxg,    nyg,     nz
    use zocnod,   only  :                                             &
    &     myrank,  iroot
    use zocfil
    use ufile

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
    real(8),   intent(inout)  ::     ft(nxdim, nydim, ntdim)

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

    if (myrank /= iroot ) then
       allocate ( buf2(1,1) )
       allocate ( buf3(1,1,1) )
       allocate (  g2d(1,1) )
       allocate (  g3d(1,1,1) )
    else        
       allocate ( buf2(nxg,nyg) )
       allocate ( buf3(nxg,nyg,nz) )
       allocate ( g2d(nxgdim,nygdim) )
       allocate ( g3d(nxgdim,nygdim,nzdim) )
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

!    if ( myrank == iroot ) then
!       read(nfinit, end=169) chead
!       read(nfinit) buf3
!       do k = 1, nz
!          do j = 1, nyg
!             do i = 1, nxg
!                g3d(igstr+i-1, jgstr+j-1, kstr+k-1) = buf3(i, j, k)
!             end do
!          end do
!       end do
!169    continue
!    end if
!    call scatter_3d(amv, g3d)
!
!    if ( myrank == iroot ) then
!       read(nfinit, end=179) chead
!       read(nfinit) buf3
!       do k = 1, nz
!          do j = 1, nyg
!             do i = 1, nxg
!                g3d(igstr+i-1, jgstr+j-1, kstr+k-1) = buf3(i, j, k)
!             end do
!          end do
!       end do
!179    continue
!    end if
!    call scatter_3d(ahv, g3d)
!   
!    if ( myrank == iroot ) then
!       read(nfinit, end=189) chead
!       read(nfinit) buf2
!       do j = 1, nyg
!          do i = 1, nxg
!             g2d(igstr+i-1, jgstr+j-1) = buf2(i, j)
!          end do
!       end do
!189    continue
!    end if
!    call scatter_2d(ft, g2d)
!    
!    if ( myrank == iroot ) then
!       read(nfinit, end=199) chead
!       read(nfinit) buf2
!       do j = 1, nyg
!          do i = 1, nxg
!             g2d(igstr+i-1, jgstr+j-1) = buf2(i, j)
!          end do
!       end do
!199    continue
!    end if
!    call scatter_2d(ft(1, 1, 2), g2d)
    
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
    call shift1(    ft,                                               &
    &            nxdim,   nydim,  ntdim,                              &
    &             1.d0,       0,      0 )
#else
    call shift3(    ub,     vb,      w,  nxdim,  nydim,  nzdim )
    call shift1(    tb,                  nxdim,  nydim, nztdim )
    call shift2(   amv,    ahv,          nxdim,  nydim,  nzdim )
    call shift3(    hb,   ubtb,   vbtb,  nxdim,  nydim,      1 )
    call shift1(    ft,                  nxdim,  nydim,  ntdim )
#endif

    if ( myrank == iroot ) then
       cdate = chead(50)
       read(cdate, '(i6.6,5i2.2)') idate
       call cyh2ss(  ttt, idate )
       write(jfpar, *) ' start time :', idate
       if ( irstrt /= 0 ) then
          if ( ttt /= tstrt ) then
             write(jfpar, *) '*** start time error ***'
             call mpi_abort(mpi_comm_world, 1, ierr)
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
    &      nxdim,  nydim,  nzdim,                                     &
    &      igstr,  jgstr,   kstr,                                     &
    &        nxg,    nyg,     nz
    use zocnod,   only  :                                             &
    &     myrank,  iroot

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

    return

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
    &        (oeof, 1, mpi_logical, iroot, mpi_comm_world, ierr)
       if ( .not. oeof ) then
          call scatter_3d(additm, g3d)
#ifdef OPT_TRIPOLE
          call shift1(additm, nxdim, nydim, nzdim,                    &
    &                          fact,  ioff,  joff  )
#else
          call shift1(additm, nxdim, nydim, nzdim)
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
    &        (oeof, 1, mpi_logical, iroot, mpi_comm_world, ierr)
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
    &         ft,                                                     &
    &     orsout, orsrwd )

    use zocdim,   only  :                                             &
    &      nxdim,  nydim,  nzdim,  ntdim,                             &
    &     nxgdim, nygdim,  igstr,  jgstr,    kstr,                    &
    &        nxg,    nyg,     nz,   nxyg,   nxyzg
    use zocnod,   only  :                                             &
    &     myrank,  iroot
    use ufile

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
    real(8),    intent(inout)  ::      ft(nxdim, nydim, ntdim)
    logical,    intent(in)     ::  orsout,    orsrwd

!---- local variables
    integer(4)        ::      i,      j,      k,     l
    integer(4), save  ::  ifpar,  jfpar,  istat
    integer(4)        ::   ierr
    namelist /nmfrst/ cfrest

    if ( ofirst ) then
       call rewnml(ifpar, jfpar)
       read(ifpar, nmfrst, iostat=istat)
       call cstnml( jfpar, 'finout', 'nmfrst', istat )
       write(jfpar, nmfrst)
       if ( myrank == iroot ) then
          call filopn(nfrest, cfrest, 'WRITE')
       end if
       ofirst = .false.
    end if

    if ( .not. orsout ) return

    if ( myrank == iroot ) then
       if (orsrwd) then
          rewind(nfrest)
       end if
       call css2yh(  idate, tt)
       write(chead(50), '(i6.6,5i2.2)') idate
       write(chead(27), '(i4.4,2i2.2,a1,3i2.2,a1)')                   &
    &        idate(1), idate(2), idate(3), ' ',                       &
    &        idate(4), idate(5), idate(6), ' '
       write(chead(29), '(i15,a1)') nxg, 'X'
       write(chead(30), '(i16)') 1
       write(chead(31), '(i16)') nxg
       write(chead(32), '(i15,a1)') nyg, 'Y'
       write(chead(33), '(i16)') 1
       write(chead(34), '(i16)') nyg
       write(chead(36), '(i16)') 1
       chead(38) = 'REAL8'
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
       chead(3) = 'U'
       write(chead(35), '(i15,a1)') nz, 'Z'
       write(chead(37), '(i16)') nz
       write(chead(64), '(i16)') nxyzg
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
       chead(3) = 'V'
       write(chead(35), '(i15,a1)') nz, 'Z'
       write(chead(37), '(i16)') nz
       write(chead(64), '(i16)') nxyzg
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
       chead(3) = 'T'
       write(chead(35), '(i15,a1)') nz, 'Z'
       write(chead(37), '(i16)') nz
       write(chead(64), '(i16)') nxyzg
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
       chead(3) = 'S'
       write(chead(35), '(i15,a1)') nz, 'Z'
       write(chead(37), '(i16)') nz
       write(chead(64), '(i16)') nxyzg
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
       chead(3) = 'SH'
       write(chead(35), '(i15,a1)') 1, 'Z'
       write(chead(37), '(i16)') 1
       write(chead(64), '(i16)') nxyg
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
       chead(3) = 'UBT'
       write(chead(35), '(i15,a1)') 1, 'Z'
       write(chead(37), '(i16)') 1
       write(chead(64), '(i16)') nxyg
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
       chead(3) = 'VBT'
       write(chead(35), '(i15,a1)') 1, 'Z'
       write(chead(37), '(i16)') 1
       write(chead(64), '(i16)') nxyg
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
       chead(3) = 'W'
       write(chead(35), '(i15,a1)') nz, 'Z'
       write(chead(37), '(i16)') nz
       write(chead(64), '(i16)') nxyzg
       write(nfrest) chead
       write(nfrest) buf3
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
       chead(3) = 'AMV'
       write(chead(35), '(i15,a1)') nz, 'Z'
       write(chead(37), '(i16)') nz
       write(chead(64), '(i16)') nxyzg
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
       chead(3) = 'AHV'
       write(chead(35), '(i15,a1)') nz, 'Z'
       write(chead(37), '(i16)') nz
       write(chead(64), '(i16)') nxyzg
       write(nfrest) chead
       write(nfrest) buf3
    end if
    
    
    call gather_2d(g2d, ft)
    if ( myrank == iroot ) then
       do j = 1, nyg
          do i = 1, nxg
             buf2(i, j) = g2d(igstr+i-1, jgstr+j-1)
          end do
       end do
       chead(3) = 'FT'
       write(chead(35), '(i15,a1)') 1, 'Z'
       write(chead(37), '(i16)') 1
       write(chead(64), '(i16)') nxyg
       write(nfrest) chead
       write(nfrest) buf2
    end if
    
    call gather_2d(g2d, ft(1,1,2))
    if ( myrank == iroot ) then
       do j = 1, nyg
          do i = 1, nxg
             buf2(i, j) = g2d(igstr+i-1, jgstr+j-1)
          end do
       end do
       chead(3) = 'FW'
       write(chead(35), '(i15,a1)') 1, 'Z'
       write(chead(37), '(i16)') 1
       write(chead(64), '(i16)') nxyg
       write(nfrest) chead
       write(nfrest) buf2
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
          write(chead(3), '(a6,i2.2)') 'TRACER', l
          write(chead(35), '(i15,a1)') nz, 'Z'
          write(chead(37), '(i16)') nz
          write(chead(64), '(i16)') nxyzg
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
          write(chead(3), '(a6,i2.2)') 'TRCFLX', l
          write(chead(35), '(i15,a1)') 1, 'Z'
          write(chead(37), '(i16)') 1
          write(chead(64), '(i16)') nxyg
          write(nfrest) chead
          write(nfrest) buf2
       end if
    end do
    
    if ( myrank == iroot ) then
       write(jfpar, *) '*** Write restart file ***'
       write(jfpar, *) ' time :', idate
       write(jfpar, *) ' step :', ntstep
    end if
    call mpi_barrier(mpi_comm_world, ierr)

  end subroutine finout

! =====================================================================

  subroutine finadd(                                                  &
    &             additm,  ixdim,  jydim,  kzdim,                     &
    &             ccitem,   clas  ) 

    use zocdim,   only  :                                             &
    &      igstr,  jgstr,   kstr,                                     &
    &        nxg,    nyg,     nz,   nxyg,  nxyzg
    use zocnod,   only  :                                             &
    &     myrank,  iroot

    implicit none

    integer(4),   intent(in)     ::   ixdim,  jydim,  kzdim
    real(8),      intent(inout)  ::  additm(ixdim, jydim, kzdim)
    character(*), intent(in)     ::  ccitem,   clas

!---- local variables
    integer(4)   ::     i,     j,     k

    
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
          chead(3) = ccitem
          write(chead(35), '(i15,a1)') nz, 'Z'
          write(chead(37), '(i16)') nz
          write(chead(64), '(i16)') nxyzg
          write(nfrest) chead
          write(nfrest) buf3
       end if
    else
       call gather_2d(g2d, additm)
       if ( myrank == iroot ) then
          do j = 1, nyg
             do i = 1, nxg
                buf2(i, j) = g2d(igstr+i-1, jgstr+j-1)
             end do
          end do
          chead(3) = ccitem
          write(chead(35), '(i15,a1)') 1, 'Z'
          write(chead(37), '(i16)') 1
          write(chead(64), '(i16)') nxyg
          write(nfrest) chead
          write(nfrest) buf2
       END IF
    END IF
    
  end subroutine finadd

end module brstt

