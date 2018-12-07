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
!     '12.10.09  M.Kurogi: for COCO5.0
!     '15.04.07  M.Kurogi: MPI-IO
! ---------------------------------------------------------------------
 use zocfil, only : nfomax, ncf
 implicit none
#include "mpif.h"
 private
 public chkset, chkstk, chkout, chekin


 integer, save   ::  nxitm(nfomax),  nyitm(nfomax),  nzitm(nfomax)
 integer, save   :: nfunit(nfomax)
 integer, save   :: itopad(nfomax), ktopad(nfomax), nszitm(nfomax)
 real(8), save   ::  fctavr(nfomax), cntavr(nfomax),  ttold(nfomax)
 logical, save   :: osingl(nfomax),   oadd(nfomax)
 integer, save   :: iflout(nfomax), ioavrg(nfomax)
 integer, save   ::  ixstr(nfomax),  ixend(nfomax)
 integer, save   :: jystr(nfomax),  jyend(nfomax)
 integer, save   :: kzstr(nfomax),  kzend(nfomax)
 integer, save   :: nohitm
 character, save :: citem(nfomax)*16, cfitem(nfomax)*(ncf)
 data citem  / nfomax*'                ' /
 data cntavr / nfomax*0.d0 /
 data iflout / nfomax*0 /
 data ioavrg / nfomax*0 /


 integer, save ::  ifpar, jfpar



contains
 subroutine chkset
 use zocdim, only : nxg, nyg, nz, myrank, iroot
 use zocfil, only : nfomax
 use zocout, only : nwork, wrkout
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
 integer :: iodint =1, iudint=1, iodavr=0, iodsng=1
 integer :: iohstr(6), iohend(6), iohint, iuhint, iohavr, iohsng
 integer :: ioxstr, ioxend, ioystr, ioyend, iozstr, iozend
 character ::  cohfil*(ncf), cohitm*16

 namelist /nmtime/ itstrt, itend, tmstp, iutstp, ntsplt
 namelist /nmdout/ iodstr, iodend, iodint, iudint, iodavr, iodsng
 namelist /nmhist/ cohitm, cohfil,                                   &
  &                  iohstr, iohend, iohint, iuhint, iohavr, iohsng, &
  &                  ioxstr, ioxend, ioystr, ioyend, iozstr, iozend

 integer :: i
 integer :: istat 
 integer, save   :: isingl(nfomax)
 data isingl / nfomax*1 /
 
 call rewnml(ifpar, jfpar)
 read(ifpar, nmtime, iostat=istat)
 call cstnml(jfpar, 'chkset', 'nmtime', istat)

 call rewnml(ifpar, jfpar)
 read(ifpar, nmdout, iostat=istat)
 call cstnml(jfpar, 'chkset', 'nmdout', istat)

 do iitem = 1, nfomax
    do i = 1, 6
       iostrt(i, iitem) = iodstr(i)
       ioend (i, iitem) = iodend(i)
    end do
    iointv(iitem) = iodint
    iuintv(iitem) = iudint
    ioavrg(iitem) = iodavr
    isingl(iitem) = iodsng
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

    read(ifpar, nmhist, iostat=istat)
    if(istat < 0) exit

    if (cohitm(1:13) /= 'not-specified') then
       nohitm = nohitm + 1
       iohitm = nohitm
       iflout(iohitm) = 1
       citem(iohitm) = cohitm
       cfitem(iohitm) = cohfil


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
       if ((iozstr > 0) .and. (iozstr <= nz)) then
          kzstr(iohitm) = iozstr
       end if
       if ((iozend >= kzstr(iohitm)) .and. (iozend <= nz)) then
          kzend(iohitm) = iozend
       end if
    end if
 end do

 do i = 1, nwork
    wrkout(i) = 0.d0
 end do
 do iitem = 1, nohitm
    if (ioavrg(iitem) == 1) then
       fctavr(iitem) = 1.d0
    else
       fctavr(iitem) = 0.d0
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
 use zocout, only :  dbleou, snglou, wrkout
 use ucaln
 use mpiio
 implicit none

 logical, intent(in) :: oflout(nfomax)
 integer :: iitem
 integer ::  ixdim,  jydim,  kzdim
 character :: chead(64)*16
 data chead  / 64*'                ' /
 real(8) :: tout
 integer ::  idate(6) 
 integer :: ijk, ijkm
 integer :: i, j, k

 integer (kind=mpi_offset_kind), save :: disp(nfomax)=0
 integer :: ifile
 integer :: istart(3), igsize(3), isize(3)
 integer :: nsize
 integer(8) :: int1, int2, int3, int4, nsize2

 character(16) :: cnx, cny, cnz

 do iitem = 1, nohitm
    if (oflout(iitem)) then
       if (iflout(iitem) == 1) then

          if (ioavrg(iitem) == 1) then
             tout = (tt + ttold(iitem)) * 0.5d0
          else
             tout = tt
          end if
          call css2yh(            &
  &                       idate,  &
  &                        tout)
          ixdim = ixend(iitem) - ixstr(iitem) + 1
          jydim = jyend(iitem) - jystr(iitem) + 1
          kzdim = kzend(iitem) - kzstr(iitem) + 1

! for GTOOL3 header
          write(chead(1), '(i16)') 9010
          chead(3) = citem(iitem)
          write(chead(14), '(16x)') ! title
          write(chead(15), '(16x)') ! title (cont.)
          write(chead(16), '(16x)') ! unit
          write(chead(25), '(i16)') nint(tout / 3.6d3)
          chead(26) = 'HOUR'
          write(chead(27), '(i4.4,2i2.2,a1,3i2.2,a1)')   &
  &              idate(1), idate(2), idate(3), ' ',      &
  &              idate(4), idate(5), idate(6), ' '
          write(chead(50), '(i6.6,5i2.2)') idate
#ifdef OPT_TRIPOLE
          write(cnx, '(i16)') nxitm(iitem)
          write(cny, '(i16)') nyitm(iitem)
          write(chead(29), '(a)') 'OCLONTPT'//trim(adjustl(cnx))
          write(chead(32), '(a)') 'OCLATTPT'//trim(adjustl(cny))
#else
          write(chead(29), '(i15,a1)') nxitm(iitem), 'X'
          write(chead(32), '(i15,a1)') nyitm(iitem), 'Y'
#endif
          write(chead(30), '(i16)') ixstr(iitem)
          write(chead(31), '(i16)') ixend(iitem)
          write(chead(33), '(i16)') jystr(iitem)
          write(chead(34), '(i16)') jyend(iitem)
          if ( nzitm(iitem) == 1 ) then
             write(chead(35), '(a)' ) 'SFC1'
          else if ( nzitm(iitem) == nic ) then
             write(chead(35), '(a)' ) 'NUMBER1000'
          else
!             write(cnz, '(i16)') nzitm(iitem)
             write(cnz, '(i16)') nz
             write(chead(35), '(a)') 'OCDEPT'//trim(adjustl(cnz))
          end if
          write(chead(36), '(i16)') kzstr(iitem)
          write(chead(37), '(i16)') kzend(iitem)
          write(chead(64), '(i16)') ixdim*jydim*kzdim

          if (ioavrg(iitem) == 1) then
             tout = ttold(iitem)
             call css2yh( idate, tout )
             write(chead(48), '(i4.4,2i2.2,a1,3i2.2,a1)')   &
                  &  idate(1), idate(2), idate(3), ' ',     &
                  &  idate(4), idate(5), idate(6), ' '
             tout = tt
             call css2yh( idate, tout )
             write(chead(49), '(i4.4,2i2.2,a1,3i2.2,a1)')   &
                  &  idate(1), idate(2), idate(3), ' ',     &
                  &  idate(4), idate(5), idate(6), ' '
             tout = (tt + ttold(iitem)) * 0.5d0
             write(chead(28), '(i16)') nint((tt-ttold(iitem)) / 3.6d3)
             call css2yh( idate, tout )
          else
             write(chead(48), '(i4.4,2i2.2,a1,3i2.2,a1)')   &
                  &  idate(1), idate(2), idate(3), ' ',     &
                  &  idate(4), idate(5), idate(6), ' '
             write(chead(49), '(i4.4,2i2.2,a1,3i2.2,a1)')   &
                  &  idate(1), idate(2), idate(3), ' ',     &
                  &  idate(4), idate(5), idate(6), ' '
          end if

          if (myrank < ijnode) then
             do k = kzstr(iitem), kzend(iitem)
                do j = 1, ny
                   do i = 1, nx
                      ijk = (k - kzstr(iitem)) * nxy     &
  &                        + (j - 1) * nx + i
                      ijkm = (k - kzstr(iitem)) * nxydim &
  &                        + (jstr + j - 2) * nxdim      &
  &                        + istr + i - 1
                      dbleou(ijk)                        &
  &                 = wrkout(itopad(iitem) + ijkm - 1)   &
  &                   / cntavr(iitem)
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
#ifdef OPT_IO_SEQUENTIAL
                  chead(38) = 'UR4'
                  call mpi_write_header(chead, nfunit(iitem), disp(iitem))
                  call info_seq(nfunit(iitem), disp(iitem), nsize)
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
     &                 3, igsize, isize, istart, &
     &                 mpi_order_fortran,        &
     &                 mpi_real4, ifile, ierr)

                  call mpi_type_commit(ifile, ierr)

                  call mpi_file_set_view(         &
     &                 nfunit(iitem),disp(iitem), &
     &                 mpi_real4,ifile,"native",  &
     &                 mpi_info_null,ierr)

                  call mpi_file_write_all(                 &
     &                 nfunit(iitem), snglou, nx*ny*kzdim, &
     &                 mpi_real4, mpi_status_ignore, ierr)

                  disp(iitem) = disp(iitem) + nsize2
#ifdef OPT_IO_SEQUENTIAL
                  call info_seq(nfunit(iitem), disp(iitem), nsize)
#endif

             else

                  int1=nxg
                  int2=nyg
                  int3=kzdim
                  int4=8
                  nsize2=int1*int2*int3*int4
                  nsize=nxg*nyg*kzdim*8
#ifdef OPT_IO_SEQUENTIAL
                  chead(38) = 'UR8'
                  call mpi_write_header(chead, nfunit(iitem), disp(iitem))                    
                  call info_seq(nfunit(iitem), disp(iitem), nsize)
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
     &                 3, igsize, isize, istart, &
     &                 mpi_order_fortran,        &
     &                 mpi_real8, ifile, ierr)

                  call mpi_type_commit(ifile, ierr)

                  call mpi_file_set_view(         &
     &                 nfunit(iitem),disp(iitem), &
     &                 mpi_real8,ifile,"native",  &
     &                 mpi_info_null,ierr)

                  call mpi_file_write_all(                 &
     &                 nfunit(iitem), dbleou, nx*ny*kzdim, &
     &                 mpi_real8, mpi_status_ignore, ierr)
                  disp(iitem) = disp(iitem) + nsize2
#ifdef OPT_IO_SEQUENTIAL
                  call info_seq(nfunit(iitem), disp(iitem), nsize)
#endif
          end if



          if (myrank < ijnode) then
          do i = 1, nszitm(iitem)
             wrkout(itopad(iitem) + i - 1) = 0.d0
          end do
          end if


       end if
       cntavr(iitem) = 0.d0
       ttold(iitem) = tt

       if (myrank == iroot) then
       write(jfpar, *) '*** File output ***'
       write(jfpar, *) ' Item :', citem(iitem),           &
  &                     '  Time :', idate, '  Step :', nt
       end if
    end if
 end do

 return
 end subroutine chkout
! =====================================================================

 subroutine chekin(                                        &
     &              ditem, ccitem,                         &
     &             nxitem, nyitem, nzitem, nditem,   clas)
 use zocdim, only : nxydim, nxg, nyg, kstr, nic
 use zocout, only : nwork, wrkout 
 use zocfil, only : nfstdo
 use ufile
 integer,      intent(in) :: nxitem, nyitem, nzitem, nditem 
 real(8),      intent(in) :: ditem(nditem)
 character(*), intent(in) :: ccitem, clas

 integer :: item, iohitm
 integer ::  jtopad=1
 character ::  ctitem*16
 integer :: i
 logical, save   ::  ofirst(nfomax)
 data ofirst / nfomax*.true. /

 do iohitm = 1, nohitm
    ctitem = ccitem
    if (citem(iohitm) == ccitem) then
       item = iohitm
       go to 123
    end if
 end do
 return
 123  continue

 if (ofirst(item)) then
    call rewnml(ifpar, jfpar)
    if (clas(1:3) == 'SFC') then
       kzstr(item) = 1
       kzend(item) = 1
    else if (clas(1:3) == 'ICE') then
       kzstr(item) = 1
       kzend(item) = nic
    end if

    nxitm(item) = nxg
    nyitm(item) = nyg
    nzitm(item) = nzitem
    nszitm(item) = nxydim * (kzend(item) - kzstr(item) + 1)

    if (iflout(item) == 1) then
       itopad(item) = jtopad
       jtopad = jtopad + nszitm(item)
       if (jtopad > nwork+1) then
          write(jfpar, *) '### WORK AREA SHORTAGE ###'
          write(jfpar, *) '### THE ITEM NUMBERED', item, &
  &                        'WILL NOT BE OUTPUT ###'
          itopad(item) = 0
          iflout(item) = 0
          jtopad = jtopad - nszitm(item)
       end if
       if (clas(1:3) == 'SFC') then
          ktopad(item) = 1
       else if (clas(1:3) == 'ICE') then
          ktopad(item) = nxydim + 1
       else if (clas(1:3) == 'OCN') then
          ktopad(item) = (kstr + kzstr(item) - 2) * nxydim + 1
       else
          write(nfstdo, *) '### NO SUCH OUTPUT CLASS AS', &
  &                        clas, '###'
          ktopad(item) = 1
       end if
    else
       itopad(item) = 0
    end if

    ofirst(item) = .false.
 end if

 if (oadd(item)) then
    if (iflout(item) == 1) then
       do i = 1, nszitm(item)
          wrkout(itopad(item) + i - 1) =               &
  &       wrkout(itopad(item) + i - 1) * fctavr(item)  &
  &     + ditem(ktopad(item) + i - 1)
       end do
    end if
 end if

 return
 end subroutine chekin
end module qckot
