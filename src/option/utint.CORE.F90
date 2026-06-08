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
!     '08.08.06  Y.Komuro: for CORE data
!     '08.08.26  Y.Komuro: diagnosing evaporation in execution
!     '08.09.02  Y.Komuro: diagnosing wind stress
!     '12.10.06  M.Kurogi: rewrite in F95 format
!
! ---------------------------------------------------------------------
  
  subroutine tmintp(  ditem, iitem  )

    use zocdim,  only  :                                              &
         nxdim,  nydim, ntdim,                                        &
           nxg,    nyg,    nx,    ny,                                 &
          istr,   jstr
    use zocgrd,  only  :     tt
    use zocfil,  only  :    ncf
    use zocnod,  only  :  iroot,  myrank, mpi_comm_ogcm
    use ufile
    use bgs2d
    use ucaln

    implicit none
#include "mpif.h"
#include "coco.h"
    real(8),    intent(inout)  ::  ditem(nxdim, nydim)
    integer(4), intent(in)     ::  iitem

!---- local variables
    integer(4), parameter      ::  nitem = ( ntdim - 2 ) * 2 + 10 
    integer(4), save           ::  idate1(6,nitem), idate2(6,nitem)
    real(8),    save           ::   time1(nitem),    time2(nitem)
    real(8),    save           ::   data1(nx,ny,nitem)
    real(8),    save           ::   data2(nx,ny,nitem)
    integer(4), save           ::  iyskiq(nitem)
    logical,    save           ::  osngld(nitem),   ofirst(nitem)
    logical,    save           ::      of,   oeof

    real(8)             ::   datag(nxg,nyg)
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
    character(len=ncf)  ::  cfusfc,   cfvsfc
    character(len=ncf)  ::  cftsfc,   cfqsfc,   cfpplr,   cfsflx
    character(len=ncf)  ::  cfswdw,   cflwdw,   cfpsfc,   cfssfc
    character(len=ncf)  ::  cftref(ntdim), cftdmp(ntdim)

    namelist /nmsfbc/  cfusfc, cfvsfc, cftsfc, cfqsfc,                  &
    &                  cfpplr, cfsflx, cfswdw, cflwdw, cfpsfc, cfssfc,  &
    &                  cftref, cftdmp


    data cfusfc, cfvsfc / 'not-specified', 'not-specified' /
    data cftsfc, cfqsfc / 'not-specified', 'not-specified' /
    data cfpplr, cfsflx / 'not-specified', 'not-specified' /
    data cfswdw, cflwdw / 'not-specified', 'not-specified' /
    data cfpsfc, cfssfc / 'not-specified', 'not-specified' /
    data cftref / ntdim*'not-specified' /
    data cftdmp / ntdim*'not-specified' /

    data osngld / nitem*.false. /
    data ofirst / nitem*.true. /
    data of / .true. /
    data iyskip / nitem*1 /

    if ( of ) then
       READ_NAMELIST( nmsfbc )
       of = .false.
    end if

    if (iitem == 1) then
       cfitem = cfusfc
    else if (iitem == 2) then
       cfitem = cfvsfc
    else if (iitem == 3) then
       cfitem = cftsfc
    else if (iitem == 4) then
       cfitem = cfqsfc
    else if (iitem == 5) then
       cfitem = cfpplr
    else if (iitem == 6) then
       cfitem = cfsflx
    else if (iitem == 7) then
       cfitem = cfswdw
    else if (iitem == 8) then
       cfitem = cflwdw
    else if (iitem == 9) then
       cfitem = cfpsfc
    else if (iitem == 10) then
       cfitem = cfssfc
    else if (iitem > nitem ) then
       write(jfpar, *) '*** TMINTP: NO SUCH ITEM ***'
       call mpi_abort(mpi_comm_ogcm, 1, ierr)
    else if ( mod(iitem, 2) == 1 ) then
       i = (iitem - 9) / 2
       cfitem = cftref(i)
    else
       i = (iitem - 10) / 2
       cfitem = cftdmp(i)
    end if

    if ( ofirst(iitem) ) then
       if ( myrank == iroot ) then
          call filopn(nfitem(iitem), cfitem, 'READ')
       end if
       READ_NAMELIST( nmskip )

       if ( iyskip(iitem) < 1 ) then
          iyskiq(iitem) = 1
       else
          iyskiq(iitem) = iyskip(iitem)
       end if
       call css2yh( idatet, tt )
       iytt = idatet(1)

       if ( myrank == iroot ) then
          rewind(nfitem(iitem))
          read(nfitem(iitem)) chead
          read(nfitem(iitem)) ((datag(i, j), i = 1, nxg), j = 1, nyg)
       end if
       call mpi_bcast(chead, 1024, mpi_character,                     &
    &                 iroot, mpi_comm_ogcm, ierr)
       call scatter_sfc(data1(1, 1, iitem), datag)
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

       if ( myrank == iroot ) then
          oeof = .false.
          read(nfitem(iitem), end=198) chead
          read(nfitem(iitem)) ((datag(i, j), i = 1, nxg), j = 1, nyg)
          go to 298
198       oeof = .true.
298       continue
       end if
       call mpi_bcast                                                 &
    &        (oeof, 1, mpi_logical, iroot, mpi_comm_ogcm, ierr)
       if (oeof) go to 98
       call mpi_bcast(chead, 1024, mpi_character,                     &
    &                 iroot, mpi_comm_ogcm, ierr)
       call scatter_sfc(data2(1, 1, iitem), datag)
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
       if ( myrank == iroot ) then
          oeof = .false.
          read(nfitem(iitem), end=197) chead
          read(nfitem(iitem)) ((datag(i, j), i = 1, nxg), j = 1, nyg)
          go to 297
197       oeof = .true.
297       continue
       end if
       call mpi_bcast                                                 &
    &        (oeof, 1, mpi_logical, iroot, mpi_comm_ogcm, ierr)
       if (oeof) go to 97
       call mpi_bcast(chead, 1024, mpi_character,                     &
    &                 iroot, mpi_comm_ogcm, ierr)
       call scatter_sfc(data2(1, 1, iitem), datag)
       cdate = chead(50)
       read(cdate, '(i6.6,5i2.2)') (idate2(i, iitem), i = 1, 6)
       iydiff = idate2(1, iitem) - idate1(1, iitem)
       call css2yh( idates, time1(iitem))
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
       if ( myrank == iroot ) then
          rewind(nfitem(iitem))
          read(nfitem(iitem)) chead
          read(nfitem(iitem)) ((datag(i, j), i = 1, nxg), j = 1, nyg)
       end if
       call mpi_bcast(chead, 1024, mpi_character,                     &
    &                 iroot, mpi_comm_ogcm, ierr)
       call scatter_sfc(data2(1, 1, iitem), datag)
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
          if ( myrank == iroot ) then
             oeof = .false.
             read(nfitem(iitem), end=897) chead
             read(nfitem(iitem)) ((datag(i, j), i = 1, nxg), j = 1, nyg)
             go to 797
897          oeof = .true.
797          continue
          end if
          call mpi_bcast                                              &
    &           (oeof, 1, mpi_logical, iroot, mpi_comm_ogcm, ierr)
          if (oeof) go to 997
          call mpi_bcast(chead, 1024, mpi_character,                  &
    &                    iroot, mpi_comm_ogcm, ierr)
          call scatter_sfc( data2(1, 1, iitem), datag )
          cdate = chead(50)
          read(cdate, '(i6.6,5i2.2)') (idate2(i, iitem), i = 1, 6)
          iydiff = idate2(1, iitem) - idate1(1, iitem)
          call css2yh(  idates,  time1(iitem) )
          idatet(1)   = idates(1) + iydiff
          idatet(2:6) = idate2(2:6,iitem)
          call cyh2ss(  time2(iitem),  idatet )
          go to 999
997       continue
          if ( myrank == iroot ) then
             oeof = .false.
             rewind(nfitem(iitem))
             read(nfitem(iitem), end=697) chead
             read(nfitem(iitem)) ((datag(i, j), i = 1, nxg), j = 1, nyg)
             go to 597
697          oeof = .true.
597          continue
          end if
          call mpi_bcast                                              &
    &           (oeof, 1, mpi_logical, iroot, mpi_comm_ogcm, ierr)
          if ( oeof ) go to 997
          call mpi_bcast(chead, 1024, mpi_character,                  &
    &                    iroot, mpi_comm_ogcm, ierr)
          call scatter_sfc(data2(1, 1, iitem), datag)
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
    use bgs3d
    use ucaln

    implicit none

#include "mpif.h"

    real(8),    intent(inout)  ::  ditem(nxdim, nydim, nzdim)
    integer(4), intent(in)     ::  iitem

!---- local variables
    integer(4), parameter      ::  nitem = 4
    integer(4), save           ::  idate1(6,nitem), idate2(6,nitem)
    real(8),    save           ::   time1(nitem),    time2(nitem)
    real(8),    save           ::   data1(nx,ny,nz,nitem)
    real(8),    save           ::   data2(nx,ny,nz,nitem)
    integer(4), save           ::  iyskiq(nitem)
    logical,    save           ::  osngld(nitem),   ofirst(nitem)
    logical,    save           ::      of,   oeof

    real(8)             ::   datag(nxg,nyg,nz)
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
    character(len=ncf)  ::  cftbdy(ntdim), cftdmb(ntdim)
    namelist /nmbody/ cftbdy, cftdmb
    data cftbdy / ntdim*'not-specified' /
    data cftdmb / ntdim*'not-specified' /

    data osngld / nitem*.false. /
    data ofirst / nitem*.true. /
    data of / .true. /
    data iyskib / nitem*1 /

    if ( of ) then
       READ_NAMELIST( nmbody )
       of = .false.
    end if

    if (iitem == 1) then
       cfitem = cftbdy(1)
    else if (iitem == 2) then
       cfitem = cftbdy(2)
    else if (iitem == 3) then
       cfitem = cftdmb(1)
    else if (iitem == 4) then
       cfitem = cftdmb(2)
    else
       write(jfpar, *) '*** TMINTB: NO SUCH ITEM ***'
       call mpi_abort(mpi_comm_ogcm, 1, ierr)
    end if

    if ( ofirst(iitem) ) then
       if ( myrank == iroot ) then
          call filopn(nfitem(iitem), cfitem, 'READ')
       end if
       READ_NAMELIST( nmskib )
       if (iyskib(iitem) < 1) then
          iyskiq(iitem) = 1
       else
          iyskiq(iitem) = iyskib(iitem)
       end if
       call css2yh(  idatet,  tt  )
       iytt = idatet(1)
       
       if ( myrank == iroot ) then
          rewind(nfitem(iitem))
          read(nfitem(iitem)) chead
          read(nfitem(iitem)) (((datag(i, j, k), i = 1, nxg),         &
    &                                            j = 1, nyg),         &
    &                                            k = 1, nz)
       end if
       call mpi_bcast(chead, 1024, mpi_character,                     &
    &                  iroot, mpi_comm_ogcm, ierr)
       call scatter_bdy(data1(1, 1, 1, iitem), datag)
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
       if ( myrank == iroot ) then
          oeof = .false.
          read(nfitem(iitem), end=198) chead
          read(nfitem(iitem)) (((datag(i, j, k), i = 1, nxg),         &
    &                                            j = 1, nyg),         &
    &                                            k = 1, nz)
          go to 298
198       oeof = .true.
298       continue
       end if
       call mpi_bcast                                                 &
    &        (oeof, 1, mpi_logical, iroot, mpi_comm_ogcm, ierr)
       if (oeof) go to 98
       call mpi_bcast(chead, 1024, mpi_character,                     &
    &                 iroot, mpi_comm_ogcm, ierr)
       call scatter_bdy(data2(1, 1, 1, iitem), datag)
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
       if ( myrank == iroot ) then
          oeof = .false.
          read(nfitem(iitem), end=197) chead
          read(nfitem(iitem)) (((datag(i, j, k), i = 1, nxg),         &
    &                                            j = 1, nyg),         &
    &                                            k = 1, nz)
          go to 297
197       oeof = .true.
297       continue
       end if
       call mpi_bcast                                                 &
    &        (oeof, 1, mpi_logical, iroot, mpi_comm_ogcm, ierr)
       if ( oeof ) go to 97
       call mpi_bcast(chead, 1024, mpi_character,                     &
    &                 iroot, mpi_comm_ogcm, ierr)
       call scatter_bdy(data2(1, 1, 1, iitem), datag)
       cdate = chead(50)
       read(cdate, '(i6.6,5i2.2)') (idate2(i, iitem), i = 1, 6)
       iydiff = idate2(1, iitem) - idate1(1, iitem)
       call css2yh( idates, time1(iitem))
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
       if ( myrank == iroot ) then
          rewind(nfitem(iitem))
          read(nfitem(iitem)) chead
          read(nfitem(iitem)) (((datag(i, j, k), i = 1, nxg),         &
    &                                            j = 1, nyg),         &
    &                                            k = 1, nz)
       end if
       call mpi_bcast(chead, 1024, mpi_character,                     &
    &                 iroot, mpi_comm_ogcm, ierr)
       call scatter_bdy(data2(1, 1, 1, iitem), datag)
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
          if ( myrank == iroot ) then
             oeof = .false.
             read(nfitem(iitem), end=897) chead
             read(nfitem(iitem)) (((datag(i, j, k), i = 1, nxg),      &
    &                                               j = 1, nyg),      &
    &                                               k = 1, nz)
             go to 797
897          oeof = .true.
797          continue
          end if
          call mpi_bcast                                              &
    &           (oeof, 1, mpi_logical, iroot, mpi_comm_ogcm, ierr)
          if (oeof) go to 997
          call mpi_bcast(chead, 1024, mpi_character,                  &
   &                     iroot, mpi_comm_ogcm, ierr)
          call scatter_bdy(data2(1, 1, 1, iitem), datag)
          cdate = chead(50)
          read(cdate, '(i6.6,5i2.2)') (idate2(i, iitem), i = 1, 6)
          iydiff = idate2(1, iitem) - idate1(1, iitem)
          call css2yh(  idates, time1(iitem) )
          idatet(1)   = idates(1) + iydiff
          idatet(2:6) = idate2(2:6, iitem)
          call cyh2ss(  time2(iitem),  idatet  )
          go to 999
997       continue
          if ( myrank == iroot ) then
             oeof = .false.
             rewind(nfitem(iitem))
             read(nfitem(iitem), end=697) chead
             read(nfitem(iitem)) (((datag(i, j, k), i = 1, nxg),      &
    &                                               j = 1, nyg),      &
    &                                               k = 1, nz)
             go to 597
697          oeof = .true.
597          continue
          end if
          call mpi_bcast                                              &
    &           (oeof, 1, mpi_logical, iroot, mpi_comm_ogcm, ierr)
          if (oeof) go to 997
          call mpi_bcast(chead, 1024, mpi_character,                  &
    &                    iroot, mpi_comm_ogcm, ierr)
          call scatter_bdy(data2(1, 1, 1, iitem), datag)
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
