module utint

  implicit none
  private

  public  ::  tmintp
#ifdef OPT_BODY
  public  ::  tmintb
#endif

contains

! --- information -----------------------------------------------------
!
!  Time interpolation of the boundary condition data
!
!  HISTORY
!     '02.10.10  H.Hasumi: from COCO3.4
!     '07.04.23  H.Hasumi
!     '08.07.08  H.Hasumi
!     '12.09.02  H.Tatebe: rewrite in F95 format
!     '15.04.07  M.Kurogi: MPI-IO
!
! ---------------------------------------------------------------------

  subroutine tmintp(  ditem, iitem  )

    use zocdim,  only  :                      &
         nxdim,  nydim, ntdim,                &
           nxg,    nyg,    nx,    ny,         &
          istr,   jstr
    use zocgrd,  only  :     tt
    use zocfil,  only  :    ncf
    use zocnod,  only  :  iroot,  myrank, mpi_comm_ogcm
    use ufile
    use ucaln
    use mpiio
    use ncfio
    use mpi
    
    implicit none

    real(8),    intent(inout)  ::  ditem(nxdim, nydim)
    integer(4), intent(in)     ::  iitem

!---- local
    integer(4), parameter      ::  nitem = 1
    integer(4), save           ::  idate1(6,nitem), idate2(6,nitem)
    real(8),    save           ::   time1(nitem),    time2(nitem)
    real(8),    save           ::   data1(nx,ny,nitem)
    real(8),    save           ::   data2(nx,ny,nitem)
    integer(4), save           ::  iyskiq(nitem)
    logical,    save           ::  osngld(nitem),   ofirst(nitem)
    logical,    save           ::      of,   oeof

    integer(4)          ::  idates(6),    idatet(6)
    integer(4)          ::       iytt,    iydiff
    real(8)             ::      times,     timet
    real(8)             ::      tintv,     tintp,    tintq
    integer(4)          ::          i,         j
    integer(4)          ::      ifpar,     jfpar,    istat,   ierr
    character(len=ncf)  ::     cfitem
    character(len=16)   ::     cdate
    character(len=16)   ::     chead(1:64)
    integer(4)          ::    iyskip(nitem)
    namelist /nmskip/ iyskip

    integer(4), save    ::  nfitem(nitem)
    character(len=ncf) :: cfsh
    namelist /nmsfbc/ cfsh
    data cfsh / 'not-specified' /

    data osngld / nitem*.false. /
    data ofirst / nitem*.true. /
    data of / .true. /
    data iyskip / nitem*1 /

    integer, save :: mpi_fh(nitem)
    integer (kind=mpi_offset_kind), save :: disp(nitem)
    integer :: icread

    logical, save :: is_ncf(nitem) = .false.
    integer :: n

    real(8) :: undef
    real(8), parameter :: undef_rtol = 1.d-2
    real(8), parameter :: undef_repl(nitem) = 0.d0

    if ( of ) then
       call rewnml( ifpar, jfpar )
       read(ifpar, nmsfbc, iostat = istat )
       call cstnml( jfpar, 'tmintp', 'nmsfbc', istat )
       of = .false.
    end if

    if (iitem == 1) then
       cfitem = cfsh
    else
       write(jfpar, *) '*** TMINTP: NO SUCH ITEM ***'
       call mpi_abort(mpi_comm_ogcm, 1, ierr)
    end if

    if ( ofirst(iitem) ) then

#ifdef OPT_IO_NCF
       n = len_trim(cfitem)
       is_ncf(iitem) = (n>=3 .and. cfitem(n-2:n) == '.nc')
#endif
       if (is_ncf(iitem)) then
          call nc_filopn(mpi_fh(iitem), cfitem, 'READ')
       else
          call mpi_filopn(mpi_fh(iitem), cfitem, 'READ')
       end if
       disp(iitem)=0

       call rewnml(ifpar, jfpar)
       read(ifpar, nmskip, iostat = istat )
       call cstnml( jfpar, 'tmintp', 'nmskip', istat )

       if ( iyskip(iitem) < 1 ) then
          iyskiq(iitem) = 1
       else
          iyskiq(iitem) = iyskip(iitem)
       end if
       call css2yh( idatet, tt )
       iytt = idatet(1)

!------------------------
       disp(iitem)=0
       if (is_ncf(iitem)) then
          call nc_read_chead(chead, mpi_fh(iitem), disp(iitem), icread)
          call nc_read_sfc(data1(1,1,iitem), mpi_fh(iitem), disp(iitem))
       else
          call mpi_read_chead(chead, mpi_fh(iitem), disp(iitem), icread)
          call mpi_read_sfc(data1(1,1,iitem), mpi_fh(iitem), disp(iitem))
       end if
       read(chead(39), *, iostat = istat) undef
       if (istat == 0) then
          where (abs(data1(:,:,iitem) - undef) <= undef_rtol * abs(undef))
             data1(:,:,iitem) = undef_repl(iitem)
          end where
       end if
!------------------------

       cdate = chead(50)
       read(cdate, '(i6.6,5i2.2)') (idate1(i, iitem), i = 1, 6)
       idatet(1) = 0
       idates(1) = 0
       idates(2:6) = idate1(2:6,iitem)
       call cyh2ss( timet, idatet )
       call cyh2ss( times, idates )
       if ( timet >= times ) then
          idates(1) = iytt
       else
          idates(1) = iytt-1
       end if
       call cyh2ss( time1(iitem), idates )

!------------------------
       if (is_ncf(iitem)) then
          call nc_read_chead(chead, mpi_fh(iitem), disp(iitem), icread)
       else
          call mpi_read_chead(chead, mpi_fh(iitem), disp(iitem), icread)
       end if
       if (icread .ne. 1024) go to 98
       if (is_ncf(iitem)) then
          call nc_read_sfc(data2(1,1,iitem), mpi_fh(iitem), disp(iitem))
       else
          call mpi_read_sfc(data2(1,1,iitem), mpi_fh(iitem), disp(iitem))
       end if
       read(chead(39), *, iostat = istat) undef
       if (istat == 0) then
          where (abs(data2(:,:,iitem) - undef) <= undef_rtol * abs(undef))
             data2(:,:,iitem) = undef_repl(iitem)
          end where
       end if
!------------------------

       cdate = chead(50)
       read(cdate, '(i6.6,5i2.2)') (idate2(i, iitem), i = 1, 6)
       iydiff = idate2(1, iitem) - idate1(1, iitem)
       idatet(1) = idates(1) + iydiff
       idatet(2:6) = idate2(2:6, iitem)
       call cyh2ss( time2(iitem), idatet )
       if (       (tt >=  time1(iitem))                               &
    &       .and. (tt <=  time2(iitem))) go to 99
90     continue
       time1(iitem) = time2(iitem)
       idate1(1:6, iitem) = idate2(1:6, iitem)
       data1(1:nx,1:ny,iitem) = data2(1:nx,1:ny,iitem)

!------------------------
       if (is_ncf(iitem)) then
          call nc_read_chead(chead, mpi_fh(iitem), disp(iitem), icread)
       else
          call mpi_read_chead(chead, mpi_fh(iitem), disp(iitem), icread)
       end if
       if (icread .ne. 1024) go to 97
       if (is_ncf(iitem)) then
          call nc_read_sfc(data2(1,1,iitem), mpi_fh(iitem), disp(iitem))
       else
          call mpi_read_sfc(data2(1,1,iitem), mpi_fh(iitem), disp(iitem))
       end if
       read(chead(39), *, iostat = istat) undef
       if (istat == 0) then
          where (abs(data2(:,:,iitem) - undef) <= undef_rtol * abs(undef))
             data2(:,:,iitem) = undef_repl(iitem)
          end where
       end if
!------------------------

       cdate = chead(50)
       read(cdate, '(i6.6,5i2.2)') (idate2(i, iitem), i = 1, 6)
       iydiff = idate2(1, iitem) - idate1(1, iitem)
       idatet(1)   = idates(1) + iydiff
       idatet(2:6) = idate2(2:6, iitem)
       call cyh2ss(  time2(iitem),  idatet )
       if (       (tt >= time1(iitem))                                &
    &       .and. (tt <= time2(iitem))) go to 99
       go to 90
97     continue
       time1(iitem) = time2(iitem)
       idate1(1:6,iitem) = idate2(1:6, iitem)
       data1(1:nx,1:ny,iitem) = data2(1:nx,1:ny,iitem)


!------------------------
       disp(iitem)=0
       if (is_ncf(iitem)) then
          call nc_read_chead(chead, mpi_fh(iitem), disp(iitem), icread)
          call nc_read_sfc(data2(1,1,iitem), mpi_fh(iitem), disp(iitem))
       else
          call mpi_read_chead(chead, mpi_fh(iitem), disp(iitem), icread)
          call mpi_read_sfc(data2(1,1,iitem), mpi_fh(iitem), disp(iitem))
       end if
       read(chead(39), *, iostat = istat) undef
       if (istat == 0) then
          where (abs(data2(:,:,iitem) - undef) <= undef_rtol * abs(undef))
             data2(:,:,iitem) = undef_repl(iitem)
          end where
       end if
!------------------------

       cdate = chead(50)
       read(cdate, '(i6.6,5i2.2)') (idate2(i, iitem), i = 1, 6)
       call css2yh(   idates,  time1(iitem)  )
       idatet(1) = idates(1)
       idatet(2:6) = idate2(2:6, iitem)
       call cyh2ss(  times,  idates  )
       call cyh2ss(  timet,  idatet  )
       if (timet <= times) then
          iydiff = iyskiq(iitem)
       else
          iydiff = iyskiq(iitem) - 1
       end if
       idatet(1) = idates(1) + iydiff
       call cyh2ss(  time2(iitem),  idatet  )
       if (       (tt >= time1(iitem))                                &
    &       .and. (tt <= time2(iitem))) go to 99

96     write(jfpar, *) '*** TMINTP: UNEXPECTED ERROR ***'
       write(jfpar, *) '-> Please inform the developer of the situation.'
       call mpi_abort(mpi_comm_ogcm, 1, ierr)
98     osngld(iitem) = .true.
99     continue
       ofirst(iitem) = .false.
    end if

    if ( osngld(iitem) ) then
       do j = 1, ny
          do i = 1, nx
             ditem(istr+i-1, jstr+j-1) = data1(i, j, iitem)
          end do
       end do
    else
       if (tt > time2(iitem)) then
          time1(iitem) = time2(iitem)
          idate1(1:6, iitem) = idate2(1:6, iitem)
          data1(1:nx,1:ny,iitem) = data2(1:nx,1:ny,iitem)

!------------------------
       if (is_ncf(iitem)) then
          call nc_read_chead(chead, mpi_fh(iitem), disp(iitem), icread)
       else
          call mpi_read_chead(chead, mpi_fh(iitem), disp(iitem), icread)
       end if
       if(icread .ne. 1024) go to 997
       if (is_ncf(iitem)) then
          call nc_read_sfc(data2(1,1,iitem), mpi_fh(iitem), disp(iitem))
       else
          call mpi_read_sfc(data2(1,1,iitem), mpi_fh(iitem), disp(iitem))
       end if
       read(chead(39), *, iostat = istat) undef
       if (istat == 0) then
          where (abs(data2(:,:,iitem) - undef) <= undef_rtol * abs(undef))
             data2(:,:,iitem) = undef_repl(iitem)
          end where
       end if
!------------------------

          cdate = chead(50)
          read(cdate, '(i6.6,5i2.2)') (idate2(i, iitem), i = 1, 6)
          iydiff = idate2(1, iitem) - idate1(1, iitem)
          call css2yh(  idates,  time1(iitem) )
          idatet(1)   = idates(1) + iydiff
          idatet(2:6) = idate2(2:6,iitem)
          call cyh2ss(  time2(iitem),  idatet )
          go to 999
997       continue

!------------------------
       disp(iitem)=0
       if (is_ncf(iitem)) then
          call nc_read_chead(chead, mpi_fh(iitem), disp(iitem), icread)
       else
          call mpi_read_chead(chead, mpi_fh(iitem), disp(iitem), icread)
       endif
       if(icread .ne. 1024) go to 997
       if (is_ncf(iitem)) then
          call nc_read_sfc(data2(1,1,iitem), mpi_fh(iitem), disp(iitem))
       else
          call mpi_read_sfc(data2(1,1,iitem), mpi_fh(iitem), disp(iitem))
       end if
       read(chead(39), *, iostat = istat) undef
       if (istat == 0) then
          where (abs(data2(:,:,iitem) - undef) <= undef_rtol * abs(undef))
             data2(:,:,iitem) = undef_repl(iitem)
          end where
       end if
!------------------------

          cdate = chead(50)
          read(cdate, '(i6.6,5i2.2)') (idate2(i, iitem), i = 1, 6)
          call css2yh(  idates,  time1(iitem) )
          idatet(1)   = idates(1)
          idatet(2:6) = idate2(2:6, iitem)
          call cyh2ss(  times,  idates  )
          call cyh2ss(  timet,  idatet  )
          if ( timet <= times ) then
             iydiff = iyskiq(iitem)
          else
             iydiff = iyskiq(iitem) - 1
          end if
          idatet(1) = idates(1) + iydiff
          call cyh2ss(  time2(iitem), idatet  )

999       continue
       end if
       tintv = time2(iitem) - time1(iitem)
       tintp = (tt - time1(iitem)) / tintv
       tintq = (time2(iitem) - tt) / tintv
       do j = 1, ny
          do i = 1, nx
             ditem(istr+i-1, jstr+j-1) = tintq * data1(i, j, iitem)   &
    &                                  + tintp * data2(i, j, iitem)
          end do
       end do

    end if

  end subroutine tmintp

#ifdef OPT_BODY

  subroutine tmintb(  ditem,  iitem  )


! --- information -----------------------------------------------------
!
!  Time interpolation of the boundary condition data
!
!  HISTORY
!     '99.10.04  H.Hasumi
!     '00.12.07  H.Hasumi: combine par and non-par routines
!     '01.05.10  H.Hasumi
!     '01.12.07  H.Hasumi
!     '02.06.02  H.Hasumi: tracer dimension
!     '12.09.02  H.Tatebe: rewrite in F95 format
!
! ---------------------------------------------------------------------

    use zocdim,  only  :                                              &
         nxdim,  nydim, nzdim, ntdim,                                 &
           nxg,    nyg,    nx,    ny,    nz,                          &
          istr,   jstr,  kstr
    use zocgrd,  only  :     tt
    use zocfil,  only  :    ncf
    use zocnod,  only  :  iroot,  myrank, mpi_comm_ogcm
    use ufile
    use ucaln
    use mpiio
    use ncfio
    use mpi

    implicit none

    real(8),    intent(inout)  ::  ditem(nxdim, nydim, nzdim)
    integer(4), intent(in)     ::  iitem

!---- local variables
    integer(4), parameter      ::  nitem = 5
    integer(4), save           ::  idate1(6,nitem), idate2(6,nitem)
    real(8),    save           ::   time1(nitem),    time2(nitem)
    real(8),    save           ::   data1(nx,ny,nz,nitem)
    real(8),    save           ::   data2(nx,ny,nz,nitem)
    integer(4), save           ::  iyskiq(nitem)
    logical,    save           ::  osngld(nitem),   ofirst(nitem)
    logical,    save           ::      of,   oeof

    integer(4)          ::  idates(6),    idatet(6)
    integer(4)          ::       iytt,    iydiff
    real(8)             ::      times,     timet
    real(8)             ::      tintv,     tintp,    tintq
    integer(4)          ::          i,         j,        k
    integer(4)          ::      ifpar,     jfpar,    istat,   ierr
    character(len=ncf)  ::     cfitem
    character(len=16)   ::     cdate
    character(len=16)   ::     chead(1:64)
    integer(4)          ::    iyskib(nitem)
    namelist /nmskib/ iyskib

    integer(4), save    ::  nfitem(nitem)
    character(len=ncf) :: cft, cfs, cfu, cfv, cfahv
    namelist /nmbody/ cft, cfs, cfu, cfv, cfahv
    data cft   / 'not-specified' /
    data cfs   / 'not-specified' /
    data cfu   / 'not-specified' /
    data cfv   / 'not-specified' /
    data cfahv / 'not-specified' /

    data osngld / nitem*.false. /
    data ofirst / nitem*.true. /
    data of / .true. /
    data iyskib / nitem*1 /

    integer, save :: mpi_fh(nitem)
    integer (kind=mpi_offset_kind), save :: disp(nitem)
    integer :: icread

    integer, save :: is_ncf(nitem) = .false.
    integer :: n
    
    real(8) :: undef
    real(8), parameter :: undef_rtol = 1.d-2
    real(8), parameter :: undef_repl(nitem) = (/10.d0, 30.d0, 0.d0, 0.d0, 0.d0/)
    
    if ( of ) then
       call rewnml( ifpar, jfpar )
       read(ifpar, nmbody, iostat = istat )
       call cstnml( jfpar, 'tmintb', 'nmbody', istat )
       of = .false.
    end if

    if (iitem == 1) then
       cfitem = cft
    else if (iitem == 2) then
       cfitem = cfs
    else if (iitem == 3) then
       cfitem = cfu
    else if (iitem == 4) then
       cfitem = cfv
    else if (iitem == 5) then
       cfitem = cfahv
    else
       write(jfpar, *) '*** TMINTB: NO SUCH ITEM ***'
       call mpi_abort(mpi_comm_ogcm, 1, ierr)
    end if

    if ( ofirst(iitem) ) then

#ifdef OPT_IO_NCF
       n = len_trim(cfitem)
       is_ncf(iitem) = (n>=3 .and. cfitem(n-2:n) == '.nc')
#endif
       if (is_ncf(iitem)) then
          call nc_filopn(mpi_fh(iitem), cfitem, 'READ')
       else
          call mpi_filopn(mpi_fh(iitem), cfitem, 'READ')
       end if
       disp(iitem)=0

       call rewnml(ifpar, jfpar)
       read(ifpar, nmskib, iostat = istat )
       call cstnml( jfpar, 'tmintb', 'nmskib', istat )
       if (iyskib(iitem) < 1) then
          iyskiq(iitem) = 1
       else
          iyskiq(iitem) = iyskib(iitem)
       end if
       call css2yh(  idatet,  tt  )
       iytt = idatet(1)

!------------------------
       disp(iitem)=0
       if (is_ncf(iitem)) then
          call nc_read_chead(chead, mpi_fh(iitem), disp(iitem), icread)
          call nc_read_bdy(data1(1,1,1,iitem), mpi_fh(iitem), disp(iitem))
       else
          call mpi_read_chead(chead, mpi_fh(iitem), disp(iitem), icread)
          call mpi_read_bdy(data1(1,1,1,iitem), mpi_fh(iitem), disp(iitem))
       end if
       read(chead(39), *, iostat = istat) undef
       if (istat == 0) then
          where (abs(data1(:,:,:,iitem) - undef) <= undef_rtol * abs(undef))
             data1(:,:,:,iitem) = undef_repl(iitem)
          end where
       end if
!------------------------

       cdate = chead(50)
       read(cdate, '(i6.6,5i2.2)') (idate1(i, iitem), i = 1, 6)
       idatet(1) = 0
       idates(1) = 0
       idates(2:6) = idate1(2:6, iitem)
       call cyh2ss(  timet,  idatet  )
       call cyh2ss(  times,  idates  )
       if ( timet >= times ) then
          idates(1) = iytt
       else
          idates(1) = iytt-1
       end if
       call cyh2ss(  time1(iitem),  idates  )

!------------------------
       if (is_ncf(iitem)) then
          call nc_read_chead(chead, mpi_fh(iitem), disp(iitem), icread)
       else
          call mpi_read_chead(chead, mpi_fh(iitem), disp(iitem), icread)
       end if
       if(icread .ne. 1024) go to 98
       if (is_ncf(iitem)) then
          call nc_read_bdy(data2(1,1,1,iitem), mpi_fh(iitem), disp(iitem))
       else
          call mpi_read_bdy(data2(1,1,1,iitem), mpi_fh(iitem), disp(iitem))
       end if
       read(chead(39), *, iostat = istat) undef
       if (istat == 0) then
          where (abs(data2(:,:,:,iitem) - undef) <= undef_rtol * abs(undef))
             data2(:,:,:,iitem) = undef_repl(iitem)
          end where
       end if
!------------------------

       cdate = chead(50)
       read(cdate, '(i6.6,5i2.2)') (idate2(i, iitem), i = 1, 6)
       iydiff = idate2(1, iitem) - idate1(1, iitem)
       idatet(1)   = idates(1) + iydiff
       idatet(2:6) = idate2(2:6, iitem)
       call cyh2ss(  time2(iitem),  idatet  )
       if (      (tt >= time1(iitem))                               &
    &      .and. (tt <= time2(iitem))) go to 99
90     continue
       time1(iitem) = time2(iitem)
       idate1(1:6, iitem) = idate2(1:6, iitem)
       data1(1:nx,1:ny,1:nz,iitem) = data2(1:nx,1:ny,1:nz,iitem)

!------------------------
       if (is_ncf(iitem)) then
          call nc_read_chead(chead, mpi_fh(iitem), disp(iitem), icread)
       else
          call mpi_read_chead(chead, mpi_fh(iitem), disp(iitem), icread)
       end if
       if(icread .ne. 1024) go to 97
       if (is_ncf(iitem)) then
          call nc_read_bdy(data2(1,1,1,iitem), mpi_fh(iitem), disp(iitem))
       else
          call mpi_read_bdy(data2(1,1,1,iitem), mpi_fh(iitem), disp(iitem))
       end if
       read(chead(39), *, iostat = istat) undef
       if (istat == 0) then
          where (abs(data2(:,:,:,iitem) - undef) <= undef_rtol * abs(undef))
             data2(:,:,:,iitem) = undef_repl(iitem)
          end where
       end if
!------------------------

       cdate = chead(50)
       read(cdate, '(i6.6,5i2.2)') (idate2(i, iitem), i = 1, 6)
       iydiff = idate2(1, iitem) - idate1(1, iitem)
       idatet(1)   = idates(1) + iydiff
       idatet(2:6) = idate2(2:6, iitem)
       call cyh2ss(  time2(iitem),  idatet  )
       if (      (tt >= time1(iitem))                                 &
    &      .and. (tt <= time2(iitem))) go to 99
       go to 90
97     continue
       time1(iitem) = time2(iitem)
       idate1(1:6,iitem) = idate2(1:6,iitem)
       data1(1:nx,1:ny,1:nz,iitem) = data2(1:nx,1:ny,1:nz,iitem)

!------------------------
       disp(iitem)=0
       if (is_ncf(iitem)) then
          call nc_read_chead(chead, mpi_fh(iitem), disp(iitem), icread)
          call nc_read_bdy(data2(1,1,1,iitem), mpi_fh(iitem), disp(iitem))
       else
          call mpi_read_chead(chead, mpi_fh(iitem), disp(iitem), icread)
          call mpi_read_bdy(data2(1,1,1,iitem), mpi_fh(iitem), disp(iitem))
       end if
       read(chead(39), *, iostat = istat) undef
       if (istat == 0) then
          where (abs(data2(:,:,:,iitem) - undef) <= undef_rtol * abs(undef))
             data2(:,:,:,iitem) = undef_repl(iitem)
          end where
       end if
!------------------------

       cdate = chead(50)
       read(cdate, '(i6.6,5i2.2)') (idate2(i, iitem), i = 1, 6)
       call css2yh(  idates,  time1(iitem) )
       idatet(1) = idates(1)
       idatet(2:6) = idate2(2:6, iitem)
       call cyh2ss( times, idates )
       call cyh2ss( timet, idatet )
       if (timet <= times) then
          iydiff = iyskiq(iitem)
       else
          iydiff = iyskiq(iitem) - 1
       end if
       idatet(1) = idates(1) + iydiff
       call cyh2ss( time2(iitem), idatet )
       if (      (tt >= time1(iitem))                                 &
    &      .and. (tt <= time2(iitem))) go to 99
96     write(jfpar, *) '*** TMINTP: UNEXPECTED ERROR ***'
       write(jfpar, *) '-> Please inform the developer of the situation.'
       call mpi_abort(mpi_comm_ogcm, 1, ierr)
98     osngld(iitem) = .true.
99     continue
       ofirst(iitem) = .false.
    end if

    if (osngld(iitem)) then
       do k = 1, nz
          do j = 1, ny
             do i = 1, nx
                ditem(istr+i-1, jstr+j-1, kstr+k-1)                   &
    &         = data1(i, j, k, iitem)
             end do
          end do
       end do
    else
       if ( tt > time2(iitem) ) then
          time1(iitem) = time2(iitem)
          idate1(1:6, iitem) = idate2(1:6, iitem)
          data1(1:nx,1:ny,1:nz,iitem) = data2(1:nx,1:ny,1:nz,iitem)

!------------------------
       if (is_ncf(iitem)) then
          call nc_read_chead(chead, mpi_fh(iitem), disp(iitem), icread)
       else
          call mpi_read_chead(chead, mpi_fh(iitem), disp(iitem), icread)
       end if
       if (icread .ne. 1024) goto 997
       if (is_ncf(iitem)) then
          call nc_read_bdy(data2(1,1,1,iitem), mpi_fh(iitem), disp(iitem))
       else
          call mpi_read_bdy(data2(1,1,1,iitem), mpi_fh(iitem), disp(iitem))
       end if
       read(chead(39), *, iostat = istat) undef
       if (istat == 0) then
          where (abs(data2(:,:,:,iitem) - undef) <= undef_rtol * abs(undef))
             data2(:,:,:,iitem) = undef_repl(iitem)
          end where
       end if
!------------------------

          cdate = chead(50)
          read(cdate, '(i6.6,5i2.2)') (idate2(i, iitem), i = 1, 6)
          iydiff = idate2(1, iitem) - idate1(1, iitem)
          call css2yh(  idates, time1(iitem) )
          idatet(1)   = idates(1) + iydiff
          idatet(2:6) = idate2(2:6, iitem)
          call cyh2ss(  time2(iitem),  idatet  )
          go to 999
997       continue

!------------------------
       disp(iitem)=0
       if (is_ncf(iitem)) then
          call nc_read_chead(chead, mpi_fh(iitem), disp(iitem), icread)
       else
          call mpi_read_chead(chead, mpi_fh(iitem), disp(iitem), icread)
       end if
       if (icread .ne. 1024) goto 997
       if (is_ncf(iitem)) then
          call nc_read_bdy(data2(1,1,1,iitem), mpi_fh(iitem), disp(iitem))
       else
          call mpi_read_bdy(data2(1,1,1,iitem), mpi_fh(iitem), disp(iitem))
       end if
       read(chead(39), *, iostat = istat) undef
       if (istat == 0) then
          where (abs(data2(:,:,:,iitem) - undef) <= undef_rtol * abs(undef))
             data2(:,:,:,iitem) = undef_repl(iitem)
          end where
       end if
!------------------------

          cdate = chead(50)
          read(cdate, '(i6.6,5i2.2)') (idate2(i, iitem), i = 1, 6)
          call css2yh(  idates,  time1(iitem) )
          idatet(1)   = idates(1)
          idatet(2:6) = idate2(2:6, iitem)
          call cyh2ss(  times,  idates )
          call cyh2ss(  timet,  idatet )
          if (timet <= times) then
             iydiff = iyskiq(iitem)
          else
             iydiff = iyskiq(iitem) - 1
          end if
          idatet(1) = idates(1) + iydiff
          call cyh2ss(  time2(iitem),  idatet )

999       continue
       end if
       tintv = time2(iitem) - time1(iitem)
       tintp = (tt - time1(iitem)) / tintv
       tintq = (time2(iitem) - tt) / tintv
       do k = 1, nz
          do j = 1, ny
             do i = 1, nx
                ditem(istr+i-1, jstr+j-1, kstr+k-1)                   &
    &         = tintq * data1(i, j, k, iitem)                         &
    &         + tintp * data2(i, j, k, iitem)
             end do
          end do
       end do

    end if

  end subroutine tmintb
#endif

end module utint


