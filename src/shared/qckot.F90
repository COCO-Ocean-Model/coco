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
!
! ---------------------------------------------------------------------
  use zocfil, only : nfomax, ncf

  implicit none
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
  character, save ::  clas(nfomax)*6
  character, save :: ctitl(nfomax)*32, cunit(nfomax)*16
  data citem  / nfomax*'                ' /
  data cntavr / nfomax*0.d0 /
  data iflout / nfomax*0 /
  data ioavrg / nfomax*0 /

  integer, save ::  ifpar, jfpar

  character(len=16), save :: chrnum
  real(8), save  ::  dundef = -1.d20

contains
  subroutine chkset
    use zocdim, only : nxg, nyg, nz, myrank, iroot
    use zocfil, only : nfomax
    use zocout, only : nwork, wrkout, dbloug, sngoug
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
    character(len=ncf) ::  crun = '(RUN NAME WAS NOT SET)'

    namelist /nmtime/ itstrt, itend, tmstp, iutstp, ntsplt
    namelist /nmdout/ iodstr, iodend, iodint, iudint, iodavr, iodsng
    namelist /nmhist/ cohitm, cohfil,                                   &
         &                  iohstr, iohend, iohint, iuhint, iohavr, iohsng, &
         &                  ioxstr, ioxend, ioystr, ioyend, iozstr, iozend
    namelist /nmrun/ crun
    
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

          if (myrank == iroot) then
             call filopn(nfunit(iohitm), cohfil, 'WRITE')
          end if
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

    allocate(dbloug(1, 1, 1))
    allocate(sngoug(1, 1, 1))
    if (myrank == iroot) then
       deallocate(dbloug, sngoug)
       allocate(dbloug(nxg, nyg, nz))
       allocate(sngoug(nxg, nyg, nz))
    end if
    
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
         &  istr, jstr, myrank, ijnode, iroot, &
         &  nic, nz
    use zocfil, only : nfomax
    use zocgrd, only : tt, nt
    use zocout, only : dbleou, dbloug, sngoug, wrkout
    use bgs3d
    use ucaln
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
             else
                if ( clas(iitem)(3:5) == 'LVM' ) then
                   write(chead(35), '(a,i0)') 'OCDEPM',  nzitm(iitem)
                else
                   write(chead(35), '(a,i0)') 'OCDEPT',  nzitm(iitem)
                end if
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
             write(chead(60), '(i4.4,2i2.2,1x,3i2.2,1x)') ivalues(1:3), ivalues(5:7)
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
                         dbleou(ijk) = &
                              & wrkout(itopad(iitem) + ijkm - 1) &
                              & / cntavr(iitem)
                      end do
                   end do
                end do
             end if
             call gather_chk(dbloug, dbleou)
             if (myrank == iroot) then
                if (osingl(iitem)) then
                   do k = 1, kzdim
                      do j = 1, nyg
                         do i = 1, nxg
                            sngoug(i, j, k) = dbloug(i, j, k)
                         end do
                      end do
                   end do
                   chead(38) = 'UR4'
                   write(nfunit(iitem)) chead
                   write(nfunit(iitem)) &
                        & (((sngoug(i, j, k), &
                        &  i = ixstr(iitem), ixend(iitem)), &
                        &  j = jystr(iitem), jyend(iitem)), &
                        &  k = 1, kzdim)
                else
                   chead(38) = 'UR8'
                   write(nfunit(iitem)) chead
                   write(nfunit(iitem)) &
                        & (((dbloug(i, j, k), &
                        &  i = ixstr(iitem), ixend(iitem)), &
                        &  j = jystr(iitem), jyend(iitem)), &
                        &  k = 1, kzdim)
                end if
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
123 continue

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
       nzitm(item) = nzitem
       nszitm(item) = nxydim * (kzend(item) - kzstr(item) + 1)
       ctitl(item) = htitl
       cunit(item) = hunit

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
          if (cclas(3:5) == 'SFC') then
             ktopad(item) = 1
          else if (cclas(3:5) == 'ICE') then
             ktopad(item) = nxydim + 1
          else if (cclas(3:5) == 'LVT' .or. cclas(3:5) == 'LVM') then
             ktopad(item) = (kstr + kzstr(item) - 2) * nxydim + 1
          else
             write(nfstdo, *) '### NO SUCH OUTPUT CLASS AS', &
                  &           cclas, '###'
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
             wrkout(itopad(item) + i - 1) = &
                  &  wrkout(itopad(item) + i - 1) * fctavr(item) &
                  & + ditem(ktopad(item) + i - 1)
          end do
       end if
    end if

    return
  end subroutine chekin
end module qckot
