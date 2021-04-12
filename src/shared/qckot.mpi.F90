module qckot
! --- information -----------------------------------------------------
!
!  Output the data
!
!  HISTORY
!     '99.04.14  H.Hasumi: from CCSR2
!     '01.05.09  H.Hasumi
!     '01.12.07  H.Hasumi
!     '02.05.29  H.Hasumi: combine parallel and nonparallel
!     '03.06.03  H.Hasumi: consistency with the change in atmct.F
!     '07.04.23  H.Hasumi
!     '07.09.25  H.Hasumi: multi-category sea ice
!     '08.??.??  Y.Komuro: SNGOUG/DBLOUG are allocated in all the nodes                           
!     '12.10.09  M.kurogi: for COCO5.0
!     '15.04.07  M.Kurogi: MPI-IO
!     '21.03.05  Y.Komuro: Sigma coordinate output
!
! ---------------------------------------------------------------------
  use zocdim, only : nxydim, nzdim
  use zocfil, only : nfomax, ncf
  use zocout

  implicit none
#include "mpif.h"
  private
  public chkset, chkstk, chkout, chekin, putsig

  integer, save   ::  nxitm(nfomax),  nyitm(nfomax),  nzitm(nfomax)
  integer, save   :: nfunit(nfomax)
  integer, save   :: itopad(nfomax), ktopad(nfomax), nszitm(nfomax)
  integer, save   :: itopat(nfomax), itopas(nfomax)
  real(8), save   ::  fctavr(nfomax), cntavr(nfomax),  ttold(nfomax)
  logical, save   :: osingl(nfomax),   oadd(nfomax)
  logical, save   :: ofctav(nfomax)
  integer, save   :: iflout(nfomax), ioavrg(nfomax)
  integer, save   ::  ixstr(nfomax),  ixend(nfomax)
  integer, save   :: jystr(nfomax),  jyend(nfomax)
  integer, save   :: kzstr(nfomax),  kzend(nfomax)
  integer, save   :: nohitm
  character, save :: citem(nfomax)*16, cfitem(nfomax)*(ncf)
  character, save ::  clas(nfomax)*6
  character, save :: ctitl(nfomax)*32, cunit(nfomax)*16
  character, save :: cvcord(nfomax)*16
  integer, save :: nvcord(nfomax), nhcord(nfomax), isvint(nfomax)
  data citem  / nfomax*'                ' /
  data cntavr / nfomax*0.d0 /
  data iflout / nfomax*0 /
  data ioavrg / nfomax*0 /

  integer, save :: nbtnb(nxydim, nchmax)
  real*8, save :: dzmnb(nxydim, nzdim, nchmax)
  real*8, allocatable, save :: sigma(:,:,:,:)
  integer, allocatable, save :: korg(:,:,:,:), ksdst(:,:,:,:)
  real*8, allocatable, save :: dkrep(:,:,:,:), thick(:,:,:,:)
  logical, save :: oscvtb(nchmax)

  real*8, allocatable, save :: &
   &      c0(:), c1(:), c2(:), c3(:), c4(:), c5(:), c6(:), &
   &      d0(:), d1(:), d2(:), d3(:), d4(:),               &
   &      d5(:), d6(:), d7(:), d8(:), d9(:)

  real*8, allocatable, save :: lsig(:,:), lsigp(:,:), dsig(:,:)
  real*8, save :: zref(0:nncmax) = 0.0d0  !! unit [cm]
  integer, save :: nsig(0:nncmax), nsnzmx
  integer, save :: nnc, nsigmx , nlist
  character(len=16), save :: cname(nncmax) = 'not-a-coordinate'

  integer, save ::  ifpar, jfpar

  character(len=16), save :: chrnum
  real(8), save  ::  dundef = -1.d20

contains
  subroutine chkset
    use zocdim, only : nxg, nyg, nx, ny, nz, myrank, iroot
    use zocfil, only : nfomax
    use ufile
    use ucaln
    implicit none
    integer :: iitem, iohitm
    real(8) ::   tstrt, tostrt

    integer :: itstrt(6) = (/ 0, 0, 0, 0, 0, 0 /)
    integer ::  itend(6) = (/ 0, 0, 0, 0, 0, 0 /)
    integer ::   iutstp  =   4 
    real(8) ::    tmstp  = 1.d10
    integer ::   ntsplt  = 0 
    integer :: iostrt(6, nfomax), ioend(6, nfomax)
    integer :: iointv(nfomax), iuintv(nfomax)
    integer :: iodstr(6) = (/ 0, 0, 0, 0, 0, 0 /)
    integer :: iodend(6) = (/ 0, 0, 0, 0, 0, 0 /)
    integer :: iodint =1, iudint=1, iodavr=0, iodsng=1, iodsvi=0
    integer :: iohstr(6), iohend(6), iohint, iuhint, iohavr, iohsng
    integer :: ioxstr, ioxend, ioystr, ioyend, iozstr, iozend, iosvin
    character ::  cohfil*(ncf), cohitm*16, cohvco*16
    character(len=ncf) ::  crun = '(RUN NAME WAS NOT SET)'

    namelist /nmtime/ itstrt, itend, tmstp, iutstp, ntsplt
    namelist /nmdout/ iodstr, iodend, iodint, iudint, iodavr, iodsng, iodsvi
    namelist /nmhist/ cohitm, cohfil, cohvco,                           &
         &            iohstr, iohend, iohint, iuhint, iohavr, iohsng,   &
         &            ioxstr, ioxend, ioystr, ioyend, iozstr, iozend, iosvin
    namelist /nmrun/ crun
    
    integer :: i, n
    integer :: istat
    integer, save   :: isingl(nfomax)
    data isingl / nfomax*1 /

    nsig(0) = nz
    nsig(1:nncmax) = 0
    call csgset
    nsnzmx = maxval(nsig)
    nworks = max(cnwrks*maxval(nsig(1:nncmax))*nxydim, 1)
    allocate(owrksg(nworks))
        
    call rewnml(ifpar, jfpar)
    read(ifpar, nmtime, iostat=istat)
    call cstnml(jfpar, 'chkset', 'nmtime', istat)
    
    call rewnml(ifpar, jfpar)
    read(ifpar, nmdout, iostat=istat)
    call cstnml(jfpar, 'chkset', 'nmdout', istat)

    call rewnml(ifpar, jfpar)
    read (ifpar, nmrun, iostat=istat)
    call cstnml(jfpar, 'chkset', 'nmrun', istat)

    if (crun(1:1) == '(') then
       chrnum = 'COCO stand-alone'
    else
       chrnum = crun(1:16)
    end if
 
    do iitem = 1, nfomax
       do i = 1, 6
          iostrt(i, iitem) = iodstr(i)
          ioend (i, iitem) = iodend(i)
       end do
       iointv(iitem) = iodint
       iuintv(iitem) = iudint
       ioavrg(iitem) = iodavr
       isingl(iitem) = iodsng
       isvint(iitem) = iodsvi
       ixstr(iitem) = 1
       jystr(iitem) = 1
       ixend(iitem) = nxg
       jyend(iitem) = nyg
       kzstr(iitem) = 1
       kzend(iitem) = nz
    end do
            
    call rewnml(ifpar, jfpar)
    nohitm = 0
    
    do
       cohitm = 'not-specified'
       cohfil = 'not-specified'
       cohvco = 'not-specified'
       iohstr(1) = -1
       iohend(1) = -1
       iohint = -1
       iuhint = -1
       iohavr = -1
       iohsng = -1
       ioxstr = -1
       ioxend = -1
       ioystr = -1
       ioyend = -1
       iozstr = -1
       iozend = -1
       iosvin = -1

       read(ifpar, nmhist, iostat=istat)
       if(istat < 0) exit

       if (cohitm(1:13) /= 'not-specified') then
          nohitm = nohitm + 1
          iohitm = nohitm
          iflout(iohitm) = 1
          citem(iohitm) = cohitm
          cfitem(iohitm) = cohfil
          cvcord(iohitm) = cohvco
          nvcord(iohitm) = 0
          do n=1, nnc
             if ( trim(cname(n)) == trim(cvcord(iohitm)) ) then
                nvcord(iohitm) = n
                exit
             end if
          end do
!         if (cohitm(1:5) == 'PDENV') then  !! special treatment for pden*
!            write(citem(iohitm)(6:7), '(i2.2)') nvcord(iohitm)
!            nvcord(iohitm) = 0
!         else if (cohitm(1:4) == 'PDEN') then
          if (cohitm(1:4) == 'PDEN') then    !! special treatment for pden
             if (nvcord(iohitm) > 0) then
                write(citem(iohitm)(5:6), '(i2.2)') nvcord(iohitm)
                nvcord(iohitm) = 0
             end if
          end if
          if (cohitm(1:5) == 'DZSIG') then
             if (nvcord(iohitm) == 0) then  !! cancel unless sigma output
                iflout(iohitm) = 0
                citem(iohitm) = '                '
                cfitem(iohitm) = '                '
                cvcord(iohitm) = '                '
                nohitm = nohitm - 1
                write(jfpar, *) &
                  &  '### DZSIG on z-coordinate will not be output ###'
                cycle
             end if
          end if      
          if (nvcord(iohitm) > 0) then        !! sigma output
             kzend(iohitm) = nsig(nvcord(iohitm))
          end if
 
          call mpi_filopn(nfunit(iohitm), cohfil, 'WRITE')

          if (iohstr(1) >= 0) then
             do i = 1, 6
                iostrt(i, iohitm) = iohstr(i)
             end do
          end if
          if (iohend(1) >= 0) then
             do i = 1, 6
                ioend(i, iohitm) = iohend(i)
             end do
          end if
          if (iohint > 0) then
             iointv(iohitm) = iohint
          end if
          if (iuhint > 0) then
             iuintv(iohitm) = iuhint
          end if
          if (iohavr >= 0) then
             ioavrg(iohitm) = iohavr
          end if
          if (iohsng >= 0) then
             isingl(iohitm) = iohsng
          end if
          if (iosvin >= 0) then
             isvint(iohitm) = iosvin
          end if
          if ((cohitm(1:5) == 'DZSIG')  & !! ignore iosvin when output DZSIG
           &  .and.(isvint(iohitm) == 0)) then
             isvint(iohitm) = 1
          end if

          if ((ioxstr > 0) .and. (ioxstr <= nxg)) then
             ixstr(iohitm) = ioxstr
          end if
          if ((ioxend >= ixstr(iohitm)) .and. (ioxend <= nxg)) then
             ixend(iohitm) = ioxend
          end if
          
          if ((ioystr > 0) .and. (ioystr <= nyg)) then
             jystr(iohitm) = ioystr
          end if
          if ((ioyend >= jystr(iohitm)) .and. (ioyend <= nyg)) then
             jyend(iohitm) = ioyend
          end if
          if ((iozstr > 0) .and. &
           &  (iozstr <= nsig(nvcord(iohitm)))) then  !! <=NZ if non-sigma
             kzstr(iohitm) = iozstr
          end if
          if ((iozend >= kzstr(iohitm)) .and. &
           &  (iozend <= nsig(nvcord(iohitm)))) then  !! <=NZ if non-sigma
             kzend(iohitm) = iozend
          end if
          if (nvcord(iohitm) > 0) then        !! sigma output
            write(jfpar, *) '*** The data will be output on ', &
              &             trim(cvcord(iohitm)), ' coordinate. ***'
          end if
       end if
    end do

    allocate(dbleou(nx*ny*nsnzmx))
    allocate(snglou(nx*ny*nsnzmx))
    
    do i = 1, nwork
       wrkout(i) = 0.d0
    end do
    do i = 1, nworks
       owrksg(i) = .false.
    end do
    do iitem = 1, nohitm
       if (ioavrg(iitem) == 1) then
          fctavr(iitem) = 1.d0
          ofctav(iitem) = .true.
       else
          fctavr(iitem) = 0.d0
          ofctav(iitem) = .false.
       end if
       if (isingl(iitem) == 1) then
          osingl(iitem) = .true.
       else
          osingl(iitem) = .false.
       end if
    end do

    call cyh2ss(          &
         &             tstrt, &
         &            itstrt)
    do iitem = 1, nohitm
       call cyh2ss(                  &
            &            tostrt,            &
            &            iostrt(1, iitem))
       ttold(iitem) = max(tostrt,  tstrt)
    end do

    return
  end subroutine chkset
! =====================================================================

  subroutine chkstk(        &
       &             oflstk)
    use zocfil, only : nfomax
    implicit none
    logical, intent(in) ::  oflstk(nfomax)
    integer :: iitem

    do iitem = 1, nohitm
       if (oflstk(iitem)) then
          cntavr(iitem) = cntavr(iitem) * fctavr(iitem) + 1.d0
          oadd(iitem) = .true.
       else
          oadd(iitem) = .false.
       end if
    end do
    
    return
  end subroutine chkstk
! =====================================================================

  subroutine chkout(        &
       &             oflout)
    use zocdim, only : nx, ny, nxy, nxg, nyg, nxdim, nxydim, &
         &  istr, jstr, myrank, ijnode, iroot, irank, jrank, ierr, &
         &  nic, nz
    use zocfil, only : nfomax
    use zocgrd, only : tt, nt
    use zocout, only : dbleou, wrkout
    use ucaln
    use mpiio
    implicit none

    logical, intent(in) :: oflout(nfomax)
    integer :: iitem
    integer ::  ixdim,  jydim,  kzdim
    character :: chead(64)*16
    data chead  / 64*'                ' /
    real(8) :: tout, tdur
    integer ::  idate(6) 
    integer :: ijk, ijkm
    integer :: i, j, k
    
    integer (kind=mpi_offset_kind), save :: disp(nfomax)=0
    integer :: ifile
    integer :: istart(3), igsize(3), isize(3)
    integer :: nsize
    integer(8) :: int1, int2, int3, int4, nsize2

    character(len=8)  :: hdate
    character(len=10) :: htime
    character(len=5)  :: hzone
    integer :: ivalues(1:8)

    do iitem = 1, nohitm
       if (oflout(iitem)) then
          if (iflout(iitem) == 1) then
             
             if (ioavrg(iitem) == 1) then
                tout = (tt + ttold(iitem)) * 0.5d0
                tdur = tt - ttold(iitem)
             else
                tout = tt
                tdur = 0.d0
             end if
             call css2yh(idate, tout)
             ixdim = ixend(iitem) - ixstr(iitem) + 1
             jydim = jyend(iitem) - jystr(iitem) + 1
             kzdim = kzend(iitem) - kzstr(iitem) + 1
             write(chead(1), '(i16)') 9010
             chead(2) = chrnum
             chead(3) = citem(iitem)
             chead(14) = ctitl(iitem)(1:16)
             chead(15) = ctitl(iitem)(17:32)
             chead(16) = cunit(iitem)
             write(chead(25), '(i16)') nint(tout / 3.6d3)
             chead(26) = 'HOUR'
             write(chead(28), '(i16)') nint(tdur / 3.6d3)
             write(chead(27), '(i4.4,2i2.2,1x,3i2.2,1x)') idate(1:6)
             write(chead(50), '(i6.6,5i2.2)') idate(1:6)
#ifdef OPT_TRIPOLE
             if (clas(iitem)(6:6) == 'V') then
                write(chead(29), '(a,i0)') 'OCLONTPV', nxitm(iitem)
                write(chead(32), '(a,i0)') 'OCLATTPV', nyitm(iitem)
             else
                write(chead(29), '(a,i0)') 'OCLONTPT', nxitm(iitem)
                write(chead(32), '(a,i0)') 'OCLATTPT', nyitm(iitem)
             end if
#else
             write(chead(29), '(i15,a1)') nxitm(iitem), 'X'
             write(chead(32), '(i15,a1)') nyitm(iitem), 'Y'
#endif
             write(chead(30), '(i16)') ixstr(iitem)
             write(chead(31), '(i16)') ixend(iitem)
             write(chead(33), '(i16)') jystr(iitem)
             write(chead(34), '(i16)') jyend(iitem)
             if ( nzitm(iitem) == 1 ) then
                write(chead(35), '(a4)' ) 'SFC1'
             else if ( clas(iitem)(3:5) == 'ICE' ) then
                write(chead(35), '(a10)' ) 'NUMBER1000'
             else if ( nvcord(iitem) == 0 ) then  !! non-sigma
                if ( clas(iitem)(3:5) == 'LVM' ) then
                   write(chead(35), '(a,i0)') 'OCDEPM',  nzitm(iitem)
                else
                   write(chead(35), '(a,i0)') 'OCDEPT',  nzitm(iitem)
                end if
             else if ( nvcord(iitem) > 0 ) then  !! sigma
                write(chead(35), '(a,i0)') 'OCDEPSIG',  nzitm(iitem)
             end if
             write(chead(36), '(i16)') kzstr(iitem)
             write(chead(37), '(i16)') kzend(iitem)
             write(chead(64), '(i16)') ixdim*jydim*kzdim
             
             write(chead(39), '(e16.7)') dundef
             chead(40) = chead(39)
             chead(41) = chead(39)
             chead(42) = chead(39)
             chead(43) = chead(39)
             write(chead(44), '(i16)') 1
             write(chead(46), '(i16)') 0
             write(chead(47), '(e16.7)') 0.d0
             
             if (ioavrg(iitem) == 1) then
                tout = ttold(iitem)
                call css2yh( idate, tout )
                write(chead(48), '(i4.4,2i2.2,1x,3i2.2,1x)') idate(1:6)
                tout = tt
                call css2yh( idate, tout )
                write(chead(49), '(i4.4,2i2.2,1x,3i2.2,1x)') idate(1:6)
                tout = (tt + ttold(iitem)) * 0.5d0
                call css2yh( idate, tout )
             else
                write(chead(48), '(i4.4,2i2.2,1x,3i2.2,1x)') idate(1:6)
                write(chead(49), '(i4.4,2i2.2,1x,3i2.2,1x)') idate(1:6)
             end if
             
             chead(61) = 'COCO'
             chead(63) = 'COCO'
             call date_and_time(hdate, htime, hzone, ivalues)
             write(chead(60), '(i4.4,2i2.2,1x,3i2.2,1x)') &
               &                             ivalues(1:3), ivalues(5:7)
             chead(62) = chead(60)
             
             if (myrank < ijnode) then
                do k = kzstr(iitem), kzend(iitem)
                   do j = 1, ny
                      do i = 1, nx
                         ijk = (k - kzstr(iitem)) * nxy &
                              &  + (j - 1) * nx + i
                         ijkm = (k - kzstr(iitem)) * nxydim &
                              &  + (jstr + j - 2) * nxdim &
                              &  +  istr + i - 1
                         if ( nvcord(iitem) > 0 ) then  !! if sigma
                            if (.not.owrksg(itopas(iitem)+ijkm-1)) then
                               dbleou(ijk) = -999.d0   !! missing
                            else
                               if ( isvint(iitem) > 0 ) then
                                  dbleou(ijk) &
                                 &  = wrkout(itopad(iitem) + ijkm - 1) &
                                 &  / cntavr(iitem)
                               else
                                  dbleou(ijk) &
                                 &  = wrkout(itopad(iitem) + ijkm - 1) &
                                 &  / wrkout(itopat(iitem)+ijkm-1)
                               end if
                            end if
                         else
                            dbleou(ijk) &
                           &  = wrkout(itopad(iitem) + ijkm - 1) &
                           &  / cntavr(iitem)
                         end if
                      end do
                   end do
                end do
             end if

             istart=(/irank*nx, jrank*ny, 0/)
             igsize=(/nxg, nyg, kzdim/)
             isize =(/nx , ny , kzdim/)

             if (osingl(iitem)) then
                int1=nxg
                int2=nyg
                int3=kzdim
                int4=4
                nsize2=int1*int2*int3*int4
                nsize=nxg*nyg*kzdim*4
#if defined(OPT_IO_SEQUENTIAL) || defined(OPT_IO_SEQUENTIAL_H8)
                chead(38) = 'UR4'
                call mpi_write_header(chead, nfunit(iitem), disp(iitem))
#ifdef OPT_IO_SEQUENTIAL
                call info_seq(nfunit(iitem), disp(iitem), nsize)
#elif defined(OPT_IO_SEQUENTIAL_H8)
                call info_seq8(nfunit(iitem), disp(iitem), nsize2)
#endif
#endif
                ! swap endian
                do k = kzstr(iitem), kzend(iitem)
                   do j = 1, ny
                      do i = 1, nx
                         ijk = (k - kzstr(iitem)) * nxy + (j - 1) * nx + i
                         snglou(ijk) = dbleou(ijk)
                         call reverse_real4(snglou(ijk)) 
                      end do
                   end do
                end do
                
                call mpi_type_create_subarray( &
                     & 3, igsize, isize, istart, &
                     & mpi_order_fortran,        &
                     &  mpi_real4, ifile, ierr)

                call mpi_type_commit(ifile, ierr)

                call mpi_file_set_view(         &
                     & nfunit(iitem),disp(iitem), &
                     & mpi_real4,ifile,"native",  &
                     & mpi_info_null,ierr)

                call mpi_file_write_all(                 &
                     & nfunit(iitem), snglou, nx*ny*kzdim, &
                     & mpi_real4, mpi_status_ignore, ierr)
                
                disp(iitem) = disp(iitem) + nsize2
#ifdef OPT_IO_SEQUENTIAL
                call info_seq(nfunit(iitem), disp(iitem), nsize)
#elif defined(OPT_IO_SEQUENTIAL_H8)
                call info_seq8(nfunit(iitem), disp(iitem), nsize2)
#endif

             else

                int1=nxg
                int2=nyg
                int3=kzdim
                int4=8
                nsize2=int1*int2*int3*int4
                nsize=nxg*nyg*kzdim*8
#if defined(OPT_IO_SEQUENTIAL) || defined(OPT_IO_SEQUENTIAL_H8)
                chead(38) = 'UR8'
                call mpi_write_header(chead, nfunit(iitem), disp(iitem))                    
#ifdef OPT_IO_SEQUENTIAL
                call info_seq(nfunit(iitem), disp(iitem), nsize)
#elif defined(OPT_IO_SEQUENTIAL_H8)
                call info_seq8(nfunit(iitem), disp(iitem), nsize2)
#endif
#endif
                ! swap endian
                do k = kzstr(iitem), kzend(iitem)
                   do j = 1, ny
                      do i = 1, nx
                         ijk = (k - kzstr(iitem)) * nxy  + (j - 1) * nx + i
                         call reverse_real8(dbleou(ijk)) 
                      end do
                   end do
                end do

                call mpi_type_create_subarray( &
                     & 3, igsize, isize, istart, &
                     & mpi_order_fortran,        &
                     & mpi_real8, ifile, ierr)
                
                call mpi_type_commit(ifile, ierr)

                call mpi_file_set_view(         &
                     & nfunit(iitem),disp(iitem), &
                     & mpi_real8,ifile,"native",  &
                     & mpi_info_null,ierr)

                call mpi_file_write_all(                 &
                     & nfunit(iitem), dbleou, nx*ny*kzdim, &
                     & mpi_real8, mpi_status_ignore, ierr)
                disp(iitem) = disp(iitem) + nsize2
#ifdef OPT_IO_SEQUENTIAL
                call info_seq(nfunit(iitem), disp(iitem), nsize)
#elif defined(OPT_IO_SEQUENTIAL_H8)
                call info_seq8(nfunit(iitem), disp(iitem), nsize2)
#endif
             end if

             if (myrank < ijnode) then
                do i = 1, nszitm(iitem)
                   wrkout(itopad(iitem) + i - 1) = 0.d0
                end do
                if (itopas(iitem) > 0) then
                   do i = 1, nszitm(iitem)
                      owrksg(itopas(iitem) + i - 1) = .false.
                   end do
                end if
                if (itopat(iitem) > 0) then
                   do i = 1, nszitm(iitem)
                      wrkout(itopat(iitem) + i - 1) = 0.d0
                   end do
                end if
             end if

          end if
          cntavr(iitem) = 0.d0
          ttold(iitem) = tt

          if (myrank == iroot) then
             write(jfpar, *) '*** File output ***'
             write(jfpar, *) ' Item :', citem(iitem),           &
                  &          ' Time :', idate, '  Step :', nt
          end if
       end if
    end do

    return
  end subroutine chkout
! =====================================================================

  subroutine chekin(                                        &
       &              ditem, ccitem,                        &
       &              htitl,  hunit,                        &
       &             nxitem, nyitem, nzitem, nditem,   cclas)
    use zocdim, only : nxydim, nxg, nyg, kstr, nic
    use zocout, only : nwork, wrkout 
    use zocfil, only : nfstdo
    use ufile
    integer,      intent(in) :: nxitem, nyitem, nzitem, nditem 
    real(8),      intent(in) :: ditem(nditem)
    character(*), intent(in) :: ccitem,  cclas,  htitl,  hunit

    real(8), allocatable :: sigitm(:), sigthk(:)
    logical, allocatable :: osigex(:)
    integer :: sdim

    integer :: item, iohitm
    integer, save :: jtopad = 1, jtopas = 1
    character ::  ctitem*16
    integer :: i
    logical, save   ::  ofirst(nfomax)
    data ofirst / nfomax*.true. /

    do iohitm = 1, nohitm
       ctitem = ccitem
       if (citem(iohitm) /= ccitem) then
          cycle
       end if
       item = iohitm

       if (ofirst(item)) then
          call rewnml(ifpar, jfpar)
          if (cclas(3:5) == 'SFC') then
             kzstr(item) = 1
             kzend(item) = 1
          else if (cclas(3:5) == 'ICE') then
             kzstr(item) = 1
             kzend(item) = nic
          end if

          nxitm(item) = nxg
          nyitm(item) = nyg
          clas(item)(1:6) = cclas(1:6)
          if ((cclas(3:5) == 'LVT' .or. cclas(3:5) == 'LVM').and. &
             & (nvcord(item) > 0)) then   !! sigma output
             nzitm(item) = nsig(nvcord(item))
          else
             nzitm(item) = nzitem
          end if
          nszitm(item) = nxydim * (kzend(item) - kzstr(item) + 1)
          ctitl(item) = htitl
          cunit(item) = hunit

!         Horizontal position: 1:T, 2:V
          if (cclas(6:6) == 'T') then
             nhcord(item) = 1
          elseif (cclas(6:6) == 'V') then
             nhcord(item) = 2
          else
             nhcord(item) = 0
          end if

          itopad(item) = 0
          itopas(item) = 0
          itopat(item) = 0
          if (iflout(item) == 1) then
             itopad(item) = jtopad
             jtopad = jtopad + nszitm(item)
             if (jtopad > nwork+1) then
                write(jfpar, *) '### WORK AREA SHORTAGE ###'
                write(jfpar, *) '### THE ITEM NUMBERED', item, &
                      &          'WILL NOT BE OUTPUT ###'
                itopad(item) = 0
                iflout(item) = 0
                jtopad = jtopad - nszitm(item)
             end if
             if ((nvcord(item) > 0).and.(nhcord(item) > 0)) then
                itopas(item) = jtopas
                jtopas = jtopas + nszitm(item)
                if (jtopas > nworks+1) then
                   write(jfpar, *) &
                   & '### WORK AREA SHORTAGE (SIGMA FLAG) ###'
                   write(jfpar, *) '### THE ITEM NUMBERED', item, &
                   &             'WILL NOT BE OUTPUT ###'
                   itopad(item) = 0
                   itopas(item) = 0
                   iflout(item) = 0
                   jtopad = jtopad - nszitm(item)
                   jtopas = jtopas - nszitm(item)
                end if
                if (isvint(item) == 0) then
                   itopat(item) = jtopad
                   jtopad = jtopad + nszitm(item)
                   if (jtopad > nwork+1) then
                      write(jfpar, *) &
                      & '### WORK AREA SHORTAGE (SIGMA THICKNESS) ###'
                      write(jfpar, *) '### THE ITEM NUMBERED', item, &
                      &             'WILL NOT BE OUTPUT ###'
                      itopad(item) = 0
                      itopas(item) = 0
                      itopat(item) = 0
                      iflout(item) = 0
                      jtopad = jtopad - 2 * nszitm(item)
                      jtopas = jtopas - nszitm(item)
                   end if
                end if
             end if

             if (cclas(3:5) == 'SFC') then
                ktopad(item) = 1
             else if (cclas(3:5) == 'ICE') then
                ktopad(item) = nxydim + 1
             else if (cclas(3:5) == 'LVT' .or. cclas(3:5) == 'LVM') then
                if (nvcord(item) > 0) then   !! sigma output
                   ktopad(item) = (kzstr(item) - 1) * nxydim + 1
                else
                   ktopad(item) = (kstr + kzstr(item) - 2) * nxydim + 1
                end if
             else
                 write(nfstdo, *) '### NO SUCH OUTPUT CLASS AS', &
                      &           cclas, '###'
                ktopad(item) = 1
             end if
          end if
          
          ofirst(item) = .false.
       end if
       
       if (oadd(item)) then
          if (iflout(item) == 1) then
             if ((nvcord(item) > 0) &
                & .and.(nhcord(item) > 0)) then   !! sigma output
                sdim = nsig(nvcord(item))
                allocate(sigitm(nxydim*sdim),sigthk(nxydim*sdim), &
                   &      osigex(nxydim*sdim))
                call cvsigc( &
                   &         sigitm, sigthk, osigex, &
                   &          ditem,                 &
                   &           sdim,  nvcord(item), nhcord(item))
                do i = 1, nszitm(item)
                   wrkout(itopad(item) + i - 1) = &
             &    wrkout(itopad(item) + i - 1) * fctavr(item) &
             &  + sigitm(ktopad(item) + i - 1)
                   owrksg(itopas(item) + i - 1) = &
             &     (owrksg(itopas(item) + i - 1) .and. ofctav(item)) &
             &     .or. osigex(ktopad(item) + i - 1)
                end do
                if (itopat(item) > 0) then
                   do i = 1, nszitm(item)
                      wrkout(itopat(item) + i - 1) = &
                &    wrkout(itopat(item) + i - 1) * fctavr(item) &
                &  + sigthk(ktopad(item) + i - 1)
                   end do
                end if
                deallocate(sigitm,sigthk,osigex)
             else
                do i = 1, nszitm(item)
                   wrkout(itopad(item) + i - 1) = &
             &    wrkout(itopad(item) + i - 1) * fctavr(item) &
             &  + ditem(ktopad(item) + i - 1)
                end do
             end if
          end if
       end if
    end do
    return
  end subroutine chekin

  subroutine csgset
    use zocdim, only: &
      &   nxdim,  nydim, nxydim,  nzdim,  ntdim,   istr,   iend, &
      &    jstr,   jend,   kstr,     le,     ln,    lne, &
      &  myrank,  iroot,   ierr
    use zocgrd, only: &
      &      dz,    dzv
    use zocmsk, only: &
#ifdef OPT_BBL
      &  amsktb, amskvb,  nbotv
#endif
      &  nbot
    use zocfil, only: &
      &     ncf
    use ufile
    use mpiio
    use bshfi

    implicit none

    integer :: mpi_fh_t
    integer(kind=mpi_offset_kind) :: disp_t
    integer :: i, j, ij, k, n, nh
    integer :: ijl, ijle, ijln, ijlne
    integer :: istat, nfsgco
    logical :: oexist

    character :: cfsgco*(ncf) = 'not-specified'
    namelist /nmsgco/ cfsgco

    call rewnml(ifpar, jfpar)
    read(ifpar, nmsgco, iostat=istat)
    call cstnml(jfpar, 'csgset', 'nmsgco', istat)

    if ( myrank == iroot ) then
       inquire(file=trim(cfsgco), exist=oexist)
    end if
    call mpi_bcast( oexist, 1, mpi_logical, &
      &             iroot, mpi_comm_world, ierr )
    if (oexist) then
       call mpi_filopn(mpi_fh_t, cfsgco, 'READ')
       disp_t = 0
       call mpi_read_root_int_sgl(nnc, mpi_fh_t, disp_t)
    else
       nnc = 0
       write(jfpar, *) &
         &  '*** Sigma-coordinate setting not found.', &
         &  ' output sigma0 only. ***'
    end if
    call mpi_bcast( nnc, 1, mpi_integer, &
      &             iroot, mpi_comm_world, ierr )
    if (nnc > 0) then
       call mpi_read_root_char(cname(1:nnc), 16, nnc, mpi_fh_t, disp_t)
       call mpi_read_root(zref(1:nnc), nnc, mpi_fh_t, disp_t)
       call mpi_read_root_int(nsig(1:nnc), nnc, mpi_fh_t, disp_t)
       call mpi_bcast( cname, 16*nncmax, mpi_character, &
         &             iroot, mpi_comm_world, ierr )
       call mpi_bcast( zref, nncmax+1, mpi_real8, &
         &             iroot, mpi_comm_world, ierr )
       call mpi_bcast( nsig, nncmax+1, mpi_integer, &
         &             iroot, mpi_comm_world, ierr )
       nsigmx=maxval(nsig)
       allocate(lsig(nsigmx,nnc),lsigp(1:nsigmx+1,nnc),dsig(nsigmx,nnc))
       call mpi_read_root(lsig, nsigmx*nnc, mpi_fh_t, disp_t)
       call mpi_read_root(lsigp, (nsigmx+1)*nnc, mpi_fh_t, disp_t)
       call mpi_read_root(dsig, nsigmx*nnc, mpi_fh_t, disp_t)
       call mpi_filcls(mpi_fh_t)
       call mpi_bcast( lsig, nsigmx*nnc, mpi_real8, &
         &             iroot, mpi_comm_world, ierr )
       call mpi_bcast( lsigp, (nsigmx+1)*nnc, mpi_real8, &
         &             iroot, mpi_comm_world, ierr )
       call mpi_bcast( dsig, nsigmx*nnc, mpi_real8, &
         &             iroot, mpi_comm_world, ierr )
       nlist = nsigmx + int(real(nzdim)*1.25)
       allocate(korg(nxydim,nlist,nnc,nchmax),  &
         &      ksdst(nxydim,nlist,nnc,nchmax), &
         &      dkrep(nxydim,nlist,nnc,nchmax), &
         &      thick(nxydim,nlist,nnc,nchmax))
       write(jfpar, *) &
         &    '*** Sigma-coordinate setting is found. ***'
       do n = 1, nnc
          write(jfpar, *) ' ID:', n, ', Name: ',trim(cname(n)), &
            &            ', Reference depth:', zref(n),         &
            &            ', # of layers:', nsig(n)
       end do
    end if

    allocate(sigma(nxydim,nzdim,0:nnc,nchmax))

    if (nnc == 0) then
       write(jfpar,*) &
         & '*** Varibles on sigma-coordinate will not be output. ***'
    end if

    call secofx

!   nbtnb(ij,1): nbot without BBL at T-point 
#ifdef OPT_BBL
    do ij = 1, nxydim
       nbtnb(ij,1) = nbot(ij) - int(amsktb(ij))
    end do
#else
    do ij = 1, nxydim
       nbtnb(ij,1) = nbot(ij)
    end do
#endif

!   nbtnb(ij,2): nbot without BBL at V-point 
!   At the last grid of each row, force toset NBTNB to KSTR-1:
!      calculating there might cause trouble in diagnosing density
!      (these grids are never reffered to for data output). 
#ifdef OPT_BBL
    do ij = 1, nxydim
       nbtnb(ij,2) = nbotv(ij) - int(amskvb(ij))
       if (mod(ij,nxdim) == 0) then
          nbtnb(ij,2) = kstr-1
       end if
    end do
#else
    do ij = 1, nxydim
       nbtnb(ij,2) = kstr - 1
    end do
    do j = jstr-1, jend+1
       do i = istr-1, iend+1
          ijl = (j-1)*nxdim + i
          ijle = ijl + le
          ijln = ijl + ln
          ijlne = ijl + lne
          if (      (nbot(ijl  ) > (kstr-1)) &
            & .and. (nbot(ijle ) > (kstr-1)) &
            & .and. (nbot(ijln ) > (kstr-1)) &
            & .and. (nbot(ijlne) > (kstr-1))) then
             nbtnb(ijl,2) = min(nbot(ijl  ), nbot(ijle ), &
               &                nbot(ijln ), nbot(ijlne))
          end if
       end do
    end do
#ifdef OPT_TRIPOLE
    call shft1i( &
      &          nbtnb(1, 2), &
      &          nxdim,  nydim)
    call shftinv( &
      &          nbtnb(1, 2), &
      &          nxdim,  nydim)
    call shft1i( &
      &          nbtnb(1, 2), &
      &          nxdim,  nydim)
#else
    call shft1i( &
      &          nbtnb(1, 2), &
      &          nxdim,  nydim)
#endif
#endif

    do nh = 1, nchmax
       do k = 1, nzdim
          do ij = 1, nxydim
             dzmnb(ij, k, nh) = 0.0d0
          end do
       end do
    end do

!   dzmnb(ij,1): dzm without BBL at T-point 
    do ij = 1, nxydim
       if (nbtnb(ij, 1) >= kstr) then
          dzmnb(ij, kstr, 1) = 0.5d0 * dz(ij, kstr)
          do k = kstr+1, nbtnb(ij, 1)
             dzmnb(ij, k, 1) = 0.5d0 * (dz(ij, k-1) + dz(ij, k))
          end do
          dzmnb(ij, nbtnb(ij, 1)+1, 1) = 0.5d0 * dz(ij, nbtnb(ij, 1))
       end if
    end do

!   dzmnb(ij,2): dzm without BBL at V-point 
    do ij = 1, nxydim
       if (nbtnb(ij, 2) >= kstr) then
          dzmnb(ij, kstr, 2) = 0.5d0 * dzv(ij, kstr)
          do k = kstr+1, nbtnb(ij, 2)
             dzmnb(ij, k, 2) = 0.5d0 * (dzv(ij, k-1) + dzv(ij, k))
          end do
          dzmnb(ij, nbtnb(ij, 2)+1, 2) = 0.5d0 * dzv(ij, nbtnb(ij, 2))
       end if
    end do

    return
  end subroutine csgset

! ======================================================================
  subroutine putsig( &
    &                     t )
    use zocdim, only: &
      &      nx,     ny,     nz, nxydim,  nzdim, nxyzdm,  ntdim, &
      &    kstr,   kend, ijtstr, ijtend,  oinit
    use zocmsk, only: &
#ifdef OPT_BBL
      &  amsktb,
#endif
      &   amskt,   nbot

    implicit none

    real(8), intent(in) :: t(nxydim, nzdim, ntdim)

    real(8) :: sigout(nxydim,nzdim)
    real(8) :: tl, sl, p1, p2
    integer :: ij, k, n, nh, nl
    character(len=7) :: cvnam
    character(len=32) :: cvmes

    do nh = 1, nchmax
       do n = 0, nnc
          do k = 1, nzdim
             do ij = 1, nxydim
                sigma(ij, k, n, nh) = 0.0d0
             end do
          end do
       end do
       do n = 1, nnc
          do nl = 1, nlist
             do ij = 1, nxydim
                korg(ij, nl, n, nh) = -1
                ksdst(ij, nl, n, nh) = -1
                dkrep(ij, nl, n, nh) = 0.0d0
                thick(ij, nl, n, nh) = 0.0d0
             end do
          end do
       end do
       oscvtb(nh) = .false.
    end do

    do n = 0, nnc
       do ij = 1, nxydim
          do k = kstr, nbot(ij)
!             tl = t(ij, k, 1) * amskt(ij, k)
!             sl = t(ij, k, 2) * amskt(ij, k)
             tl = t(ij, k, 1)
             sl = t(ij, k, 2)
             p1 = c0(n)                                               &
         &  + (c1(n) + (c2(n) + c3(n) * tl) * tl) * tl                &
         &  + (c4(n) + c5(n) * tl + c6(n) * sl) * sl
             p2 = d0(n)                                               &
         &  + (d1(n) + (d2(n) + (d3(n) + d4(n) * tl) * tl) * tl) * tl &
         &  + (d5(n) + (d6(n) + d7(n) * tl * tl) * tl                 &
                     + (d8(n) + d9(n) * tl * tl) * sqrt(sl)) * sl
!             sigma(ij, k, n, 1) = (p1 / p2 - 1.0d3) * amskt(ij, k)
             sigma(ij, k, n, 1) = p1 / p2 - 1.0d3
          end do
#ifdef OPT_BBL
          sigma(ij, kend, n, 1) =                                     &
            &       sigma(ij, nbot(ij), n, 1) * amsktb(ij)            &
            &       + sigma(ij, kend, n, 1) * (1.0d0 - amsktb(ij))
          sigma(ij, nbot(ij), n, 1) = sigma(ij, nbot(ij), n, 1)       &
            &                  * (1.0d0-amsktb(ij))
#endif
       end do
    end do

    if (oinit) then
       return
    end if

    do n=0, nnc
       do k=1, nzdim
          do ij=1, nxydim
             sigout(ij, k) = sigma(ij, k, n, 1)
          end do
       end do
#ifdef OPT_BBL
       do ij=ijtstr, ijtend
          k = nbot(ij)
          sigout(ij, k) = sigout(ij, kend) * amsktb(ij) &
            &       + sigout(ij, k) * (1.0d0 - amsktb(ij))
       end do
#endif
       cvmes = '                                '
       if (n == 0) then
          cvnam = 'PDEN   '
          cvmes = 'potential density at T point    '
       else
          write(cvnam, '(a4,i2.2,1x)') 'PDEN', n
          cvmes = 'pden. at T-pt (' // trim(cname(n)) // ')' 
       end if
       call chekin(sigout, cvnam, &
         &      cvmes, 'kg/m^3',  &
         &      nx, ny, nz, nxyzdm, 'OCLVTT')
!         do k=1, nzdim
!            do ij=1, nxydim
!               sigout(ij, k) = sigma(ij, k, n, 2)
!            end do
!         end do
!         do ij=ijvstr, ijvend
!            k = nbotv(ij)
!            sigout(ij, k) = sigout(ij, kend) * amskvb(ij) &
!              &       + sigout(ij, k) * (1.0d0 - amskvb(ij))
!         end do
!         cvmes = '                                '
!         if (n == 0) then
!            cvnam = 'PDENV  '
!            cvmes = 'potential density at V point    '
!         else
!            write(cvnam, '(a5,i2.2)') 'PDENV', n
!            cvmes = 'pden. at V-pt (' // trim(cname(n)) // ')' 
!         end if
!         call chekin(sigout, cvnam, &
!           &      'potential density at V point', 'kg/m^3', &
!           &      nx, ny, nz, nxyzdm, 'OCLVTV')
       if (n == 0) then
          cycle
       end if

       do k=1, nzdim
          do ij=1, nxydim
             sigout(ij, k) = 1.0d0
          end do
       end do
       cvmes = 'Sigma thickness(' // trim(cname(n)) // ')' 
       call chekin(sigout, 'DZSIG', &
         &      cvmes, 'cm',        &
         &      nx, ny, nz, nxyzdm, 'OCLVTT')
    end do

    return
  end subroutine putsig
! ======================================================================
  subroutine cvsigc( &
    &                sigitm, sigthk, osigex, &
    &                zitm3d,                 &
    &                  sdim,  ncsig,    nch )

    use zocdim, only: &
      &   nxdim,  nzdim, nxydim,   kstr,   kend, &
      &      le,     ln,    lne
#ifdef OPT_BBL
    use zocmsk, only: &
      &  amskvb
#endif
    implicit none

    integer, intent(in)  ::  sdim,  ncsig,    nch
    real(8), intent(out) :: sigitm( nxydim,   sdim), &
      &                     sigthk( nxydim,   sdim)
    logical, intent(out) :: osigex( nxydim,   sdim)
    real(8), intent(in)  :: zitm3d( nxydim,  nzdim)

    REAL(8) :: zitmi
    integer :: ij, k, n, nl

    if (.not.oscvtb(nch)) then  !! first conversion in each step
       if (nch == 2) then  !! calculate sigma for V-point
          do n = 1, nnc
             do ij = 1, nxydim - nxdim - 1
                do k = kstr, nbtnb(ij,nch)
                   sigma(ij,k,n,nch) = 0.25d0 * &
              &    ( sigma(ij    ,k,n,1) + sigma(ij+le ,k,n,1)   &
              &    + sigma(ij+ln ,k,n,1) + sigma(ij+lne,k,n,1) )
                end do
#ifdef OPT_BBL
                sigma(ij,kend,n,nch) = 0.25d0 *                      &
           &    ( sigma(ij    ,kend,n,1) + sigma(ij+le ,kend,n,1)    &
           &    + sigma(ij+ln ,kend,n,1) + sigma(ij+lne,kend,n,1) )  &
           &      * amskvb(ij)                                       &
           &    + sigma(ij,kend,n,nch) * (1.0d0 - amskvb(ij))
#endif
             end do
          end do
       end if
       call mkcvtb( nch )
       oscvtb(nch) = .true.
    end if

    sigitm(:,:) = 0.0d0
    sigthk(:,:) = 0.0d0
    osigex(:,:) = .false.
    do ij=1, nxydim
       nl = 1
       do while (korg(ij,nl,ncsig,nch) > 0)
          zitmi =                                       &
            &   (1.0d0-dkrep(ij,nl,ncsig,nch))          &
            &     * zitm3d(ij,korg(ij,nl,ncsig,nch)-1)  &
            &   + dkrep(ij,nl,ncsig,nch)                &
            &     * zitm3d(ij,korg(ij,nl,ncsig,nch))
          sigitm(ij,ksdst(ij,nl,ncsig,nch)) =           &
            &    sigitm(ij,ksdst(ij,nl,ncsig,nch))      &
            &    + thick(ij,nl,ncsig,nch) * zitmi
          sigthk(ij,ksdst(ij,nl,ncsig,nch)) =           &
            &    sigthk(ij,ksdst(ij,nl,ncsig,nch))      &
            &    + thick(ij,nl,ncsig,nch)
          osigex(ij,ksdst(ij,nl,ncsig,nch)) = .true.
          nl = nl + 1
       end do
    end do

    return
  end subroutine cvsigc
! **********************************************************************
  subroutine mkcvtb( &
      &                nch )
    use zocdim, only: &
      &   ijstr,  ijend,   kstr,   kend
    use zocgrd, only: &
      &     dz0
#ifdef OPT_BBL
    use zocmsk, only: &
      &  amskvb
#endif
    implicit none
    integer, intent(in) :: nch

    real(8) :: csgn, dk0, dkn
    integer :: ij, k, ks, ksb, n, nl, nh

    do n = 1, nnc
       do ij = ijstr, ijend
          if (nbtnb(ij,nch) < kstr) then
             cycle
          end if
          nl = 1
!         Search a bin to start with 
          do k=1, nsig(n)
             if (lsigp(k+1,n) > sigma(ij,kstr,n,nch)) then
                ks = k
                korg(ij,nl,n,nch) = kstr
                ksdst(ij,nl,n,nch) = ks
                dkrep(ij,nl,n,nch) = 1.0d0
                thick(ij,nl,n,nch) = dzmnb(ij,kstr,nch)
                exit
             end if
          end do
          if (korg(ij,1,n,nch) < 0) then
             write(jfpar,*) &
               &  ' ### PUTSIG: FAIL TO FIND THE STARTING BIN. ###'
!             write(jfpar,*) ij,sigma(ij,kstr,n,nch)
             stop
          end if
 
          do k = kstr+1, nbtnb(ij,nch)
!            csgn: +1 for normal stratification, -1 for reverse
             csgn = sign(1.0d0, sigma(ij,k,n,nch)-sigma(ij,k-1,n,nch))
             dk0 = 0.0d0
             do
                nl = nl + 1
                ksb = ks + (1+int(csgn))/2
                if (((csgn > 0.0d0)                                 &
                  &  .and.(lsigp(ksb,n) > sigma(ij,k,n,nch))).or.  &
                  &  ((csgn < 0.0d0)                                &
                  &  .and.(lsigp(ksb,n) <= sigma(ij,k,n,nch)))) then
                   korg(ij,nl,n,nch) = k
                   ksdst(ij,nl,n,nch) = ks
                   dkrep(ij,nl,n,nch) = 0.5d0*(dk0+1.0d0)
                   thick(ij,nl,n,nch) = dzmnb(ij,k,nch)*(1.0d0-dk0)
                   exit
                end if
                dkn = (lsigp(ksb,n)-sigma(ij,k-1,n,nch))           &
                  & / (sigma(ij,k,n,nch)-sigma(ij,k-1,n,nch))
                korg(ij,nl,n,nch) = k
                ksdst(ij,nl,n,nch) = ks
                dkrep(ij,nl,n,nch) = 0.5d0*(dk0+dkn)
                thick(ij,nl,n,nch) = dzmnb(ij,k,nch)*(dkn-dk0)
                ks = ks + int(csgn)
                dk0 = dkn
             end do
          end do

          nl = nl + 1
          k = nbtnb(ij,nch)+1
          korg(ij,nl,n,nch) = k
          ksdst(ij,nl,n,nch) = ks
          dkrep(ij,nl,n,nch) = 0.0d0
          thick(ij,nl,n,nch) = dzmnb(ij,k,nch)
             
#ifdef OPT_BBL
          if (amskvb(ij) == 1.0d0) then
             nl = nl + 1
!            Search a bin for BBL.
             do k=1, nsig(n)
                if (lsigp(k+1,n) > sigma(ij,kend,n,nch)) then
                   korg(ij,nl,n,nch) = kend
                   ksdst(ij,nl,n,nch) = k
                   dkrep(ij,nl,n,nch) = 1.0d0
                   thick(ij,nl,n,nch) = dz0(kend)
                   exit
                end if
             end do
             if (korg(ij,nl,n,nch) < 0) then
                write(jfpar,*) &
                  &  ' ### PUTSIG: FAIL TO FIND THE BIN FOR BBL. ###'
!                write(jfpar,*) i,j,sigma(ij,kend,n,nch)
                stop
             end if
          end if
#endif
       end do
    end do

    return
  end subroutine mkcvtb
! **********************************************************************
  subroutine secofx

! --- information -----------------------------------------------------
!
!  Coefficients for the approximated equation of state by McDougall et
! al. (2003, JAOT).
!
!  HISTORY
!     '19.06.04  Y.Komuro: from xprst.F
!
! ---------------------------------------------------------------------

    implicit none
    real*8, parameter :: &
      &      p10     =  9.99843699d+2,  &
      &      p1t     =  7.35212840d+0,  &
      &      p1tt    = -5.45928211d-2,  &
      &      p1ttt   =  3.98476704d-4,  &
      &      p1s     =  2.96938239d+0,  &
      &      p1st    = -7.23268813d-3,  &
      &      p1ss    =  2.12382341d-3,  &
      &      p1p     =  1.04004591d-2,  &
      &      p1ptt   =  1.03970529d-7,  &
      &      p1ps    =  5.18761880d-6,  &
      &      p1pp    = -3.24041825d-8,  &
      &      p1pptt  = -1.23869360d-11, &
      &      p20     =  1.d0,           &
      &      p2t     =  7.28606739d-3,  &
      &      p2tt    = -4.60835542d-5,  &
      &      p2ttt   =  3.68390573d-7,  &
      &      p2tttt  =  1.80809186d-10, &
      &      p2s     =  2.14691708d-3,  &
      &      p2st    = -9.27062484d-6,  &
      &      p2sttt  = -1.78343643d-10, &
      &      p2ss    =  4.76534122d-6,  &
      &      p2sstt  =  1.63410736d-9,  &
      &      p2p     =  5.30848875d-6,  &
      &      p2ppttt = -3.03175128d-16, &
      &      p2pppt  = -1.27934137d-17

    real*8 :: z    !! unit [m]
    integer :: n

    allocate( &
      &   c0(0:nnc), c1(0:nnc), c2(0:nnc), c3(0:nnc),            &
      &   c4(0:nnc), c5(0:nnc), c6(0:nnc),                       &
      &   d0(0:nnc), d1(0:nnc), d2(0:nnc), d3(0:nnc), d4(0:nnc), &
      &   d5(0:nnc), d6(0:nnc), d7(0:nnc), d8(0:nnc), d9(0:nnc))

    do n = 0, nnc
       z = zref(n) / 1.0d2
       c0(n) = p10 + (p1p + p1pp * z) * z
       c1(n) = p1t
       c2(n) = p1tt + (p1ptt + p1pptt * z) * z
       c3(n) = p1ttt
       c4(n) = p1s + p1ps * z
       c5(n) = p1st
       c6(n) = p1ss

       d0(n) = p20 + p2p * z
       d1(n) = p2t + p2pppt * z * z * z
       d2(n) = p2tt
       d3(n) = p2ttt + p2ppttt * z * z
       d4(n) = p2tttt
       d5(n) = p2s
       d6(n) = p2st
       d7(n) = p2sttt
       d8(n) = p2ss
       d9(n) = p2sstt
    end do

    return
  end subroutine secofx
end module qckot
