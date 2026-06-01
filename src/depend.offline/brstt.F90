module brstt

  use zocfil,   only  :   ncf
  use zocdim,   only  :   nxg, nyg, nz, nxgdim, nygdim, nzdim, ntdim

  implicit none

#include "mpif.h"

  private

  real(8), save ::  g2d(nxgdim, nygdim)
  real(8), save ::  g3d(nxgdim, nygdim, nzdim)
  real(8), save :: buf2(nxg, nyg) 
  real(8), save :: buf3(nxg, nyg, nz)  

  integer, save :: mpi_fh_w, mpi_fh_r
  integer(kind=mpi_offset_kind), save :: disp, dispw=0
  integer :: icread
  
  integer(4) :: idate(6)
  character(len=16) :: chead(1:64)
  character(len=16) :: cdate
  integer(4), save ::  ifpar,  jfpar
  
  integer(4), save ::   nfinit, nfrest
   
  character(len=ncf) :: cfinit, cfrest
  namelist /nmfini/ cfinit, irstrt
  namelist /nmfrst/ cfrest
  
  data cfinit, cfrest / 'not-specified', 'not-specified' /
  integer(4), save ::  irstrt
  data irstrt / 0 /
  public :: restrt, finout

  character(len=16) :: ctrnam(ntdim), cftnam(ntdim)
  data ctrnam(1:2) / 'TO', 'SO' /
  data cftnam(1:2) / 'FT', 'FW' / 
#ifdef OPT_OECO2
  data ctrnam(3:ntdim) /                                &
       &         'NO3',  'PHY',  'ZOO',  'DET',  'CA',  &
       &       'CACO3', 'TCO2',  'ALK',   'O2', 'FED',  &
       &       'DETFE', 'DIAZ',  'PO4',  'N2O', 'AGE' /
#endif  

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
    &     myrank,  iroot, mpi_comm_ogcm, ierr
    use zocphy,   only  :                                             &
    &       dtds
    use zocfil
    use ufile
    use ucaln
    use bgs2d
    use bgs3d
    use bshft
    use mpiio
    use zocgrd,   only  : hic
    use zocout,   only  : loglev

    implicit none
   
#ifdef OPT_NPZD
    real(8), intent(inout) :: pco2o(nxdim,nydim), ff(nxdim,nydim)
#endif
    real(8), intent(inout) ::     t(nxdim, nydim, nzdim, ntdim)
    real(8), intent(inout) ::    ft(nxdim, nydim, ntdim)
    real(8), intent(in)    :: tstrt
    real(8)    ::    tt,    ttt
    integer(4) ::     i,      j,      k,      l
    integer(4), save :: istat

    call rewnml(ifpar, jfpar)
    read(ifpar, nmfini, iostat=istat)
    call cstnml( jfpar, 'restrt', 'nfini', istat )
    write(jfpar, nmfini)
    
    call filopn(nfinit, cfinit, 'READ')
    write(jfpar, *) '*** read initial condition file ***'
      
    t  = 0.d0
    ft = 0.d0

    do l = 3, ntdim
       if (myrank == iroot) then
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
       call scatter_3d(t(1, 1, 1, l), g3d)
    end do
    
    do l = 3, ntdim
       if (myrank == iroot) then
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

#ifdef OPT_NPZD
    if (myrank == iroot) then
       read(nfinit, end=229) chead
       read(nfinit) buf2
       do j = 1, nyg
          do i = 1, nxg
             g2d(igstr+i-1, jgstr+j-1) = buf2(i, j)
          end do
       end do
229    continue
    end if
    call scatter_2d(pco2o, g2d)
    if (myrank == iroot) then
       read(nfinit, end=239) chead
       read(nfinit) buf2
       do j = 1, nyg
          do i = 1, nxg
             g2d(igstr+i-1, jgstr+j-1) = buf2(i, j)
          end do
       end do
239    continue
    end if
    call scatter_2d(ff, g2d)
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
     &            nxdim,  nydim,      1)
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
             call mpi_abort(mpi_comm_world, 1, ierr)
          end if
       end if
    end if
    
    call filcls(nfinit)

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
   &     myrank,  iroot, mpi_comm_ogcm, ierr
    use bgs2d
    use bgs3d
    use ufile
    use ucaln
    use mpiio
    
    implicit none

#ifdef OPT_NPZD
    real(8),   intent(in) ::   pco2o(nxdim, nydim), ff(nxdim, nydim)
#endif 
    real(8),   intent(in) ::       t(nxdim, nydim, nzdim, ntdim)
    real(8),   intent(in) ::      ft(nxdim, nydim, ntdim)
    real(8),   intent(in) ::      hb(nxdim, nydim)
    real(8),   intent(in) ::      tt
    integer(4),intent(in) ::  ntstep
    logical,   intent(in) ::  orsout, orsrwd
!---- local
    integer(4)         ::     i,      j,      k,      l
    integer(4), save   :: istat
    logical,    save   :: ofirst = .true.

    if (ofirst) then
       call rewnml(ifpar, jfpar)
       read(ifpar, nmfrst, iostat=istat)
       call cstnml( jfpar, 'finout', 'nmfrst', istat )
       write(jfpar, nmfrst)
       call filopn(nfrest, cfrest, 'WRITE')
       ofirst = .false.
    end if
    
    if (.not. orsout) return

    if (myrank ==iroot) then
       if (orsrwd) then
          rewind(nfrest)
       end if
       call css2yh( idate,    tt )

       write(chead(1) , '(i16)' ) 9010
       write(chead(25), '(i16)'  ) 0
       chead(26) = 'HOUR'
       write(chead(27), '(i4.4,2i2.2,a1,3i2.2,a1)')      &
         &        idate(1), idate(2), idate(3), ' ',     &
         &        idate(4), idate(5), idate(6), ' '
       write(chead(29), '(a8,i0)') 'OCLONTPT', nxg
       write(chead(30), '(i16)'  ) 1
       write(chead(31), '(i16)'  ) nxg
       write(chead(32), '(a8,i0)') 'OCLATTPT', nyg
       write(chead(33), '(i16)'  ) 1
       write(chead(34), '(i16)'  ) nyg
       write(chead(36), '(i16)'  ) 1
       chead(38) = 'UR8'
       write(chead(50), '(i6.6,5i2.2)') idate
    end if

    do l = 3, ntdim
       call gather_3d(g3d, t(1, 1, 1, l))
       if (myrank == iroot) then
          do k = 1, nz
             do j = 1, nyg
                do i = 1, nxg
                   buf3(i, j, k) = g3d(igstr+i-1, jgstr+j-1, kstr+k-1)
                end do
             end do
          end do
          write(chead(3) , '(a) '     ) trim(ctrnam(l))
          write(chead(35), '(a6,i2)'  ) 'OCDEPT', nz
          write(chead(37), '(i16)'    ) nz
          write(chead(64), '(i16)'    ) nxyzg
          write(nfrest) chead
          write(nfrest) buf3
       end if
    end do
    
    do l = 3, ntdim
       call gather_2d(g2d, ft(1, 1, l))
       if (myrank == iroot) then
          do j = 1, nyg
             do i = 1, nxg
                buf2(i, j) = g2d(igstr+i-1, jgstr+j-1)
             end do
          end do
          write(chead(3) , '(a6,i2.2)') 'TRCFLX', l
          write(chead(35), '(a4)'     ) 'SFC1'
          write(chead(37), '(i16)'    ) 1
          write(chead(64), '(i16)'    ) nxyg
          write(nfrest) chead
          write(nfrest) buf2
       end if
    end do

#ifdef OPT_NPZD
    call gather_2d(g2d, pco2o)
    if (myrank == iroot) then
       do j = 1, nyg
          do i = 1, nxg
             buf2(i, j) = g2d(igstr+i-1, jgstr+j-1)
          end do
       end do
       write(chead(3) , '(a4)' ) 'pCO2'
       write(chead(35), '(a4)' ) 'SFC1'
       write(chead(37), '(i16)') 1
       write(chead(64), '(i16)') nxyg
       write(nfrest) chead
       write(nfrest) buf2
    end if
    call gather_2d(g2d, ff)
    if (myrank == iroot) then
       do j = 1, nyg
          do i = 1, nxg
             buf2(i, j) = g2d(igstr+i-1, jgstr+j-1)
          end do
       end do
       write(chead(3) , '(a2)' ) 'FF'
       write(chead(35), '(a4)' ) 'SFC1'
       write(chead(37), '(i16)') 1
       write(chead(64), '(i16)') nxyg
       write(nfrest) chead
       write(nfrest) buf2
    end if
#endif

    if (myrank == iroot) then
       write(jfpar, *) '*** Write restart file ***'
       write(jfpar, *) ' time :', idate
       write(jfpar, *) ' step :', ntstep
    end if
    call mpi_barrier(mpi_comm_world, ierr)
    
  end subroutine finout

end module brstt
