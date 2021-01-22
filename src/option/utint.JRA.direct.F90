module utint
    use zocdim,  only  : nxy
  implicit none

  private

  public  ::  tmintp, tmintp_direct
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
  
  subroutine tmintp_direct(  ditem, iitem  )

    use zocdim,  only  :                                              &
         nxdim,  nydim, ntdim,                                        &
           nxg,    nyg,    nx,    ny,                                 &
          istr,   jstr
    use zocgrd,  only  :     tt, dx,dy,hxt,hyt
    use zocfil,  only  :    ncf
    use zocnod,  only  :  iroot,  myrank, irank, jrank
    use zocphy,  only  :  dwatr
    use ufile
    use ucaln
#ifdef OPT_IO_COCOMPI
    use mpiio
#else
    use bgs2d
#endif

    implicit none

#include "mpif.h"

    real(8),    intent(inout)  ::  ditem(nxdim, nydim)
    integer(4), intent(in)     ::  iitem

!---- local variables
    integer(4), parameter      ::  nitem = ( ntdim - 2 ) * 2 + 11 
    integer(4), save           ::  idate1(6,nitem), idate2(6,nitem)
    real(8),    save           ::   time1(nitem),    time2(nitem)
    real(8),    save           ::   data1(nx,ny,nitem), data2(nx,ny,nitem)
    logical,    save           ::  of=.true., ofirst(nitem)=.true.

    real(8)             ::   data(nx,ny), time
    integer(4)          ::  idates(6),    idatet(6)
    real(8)             ::      tintv,     tintp,    tintq
    integer(4)          ::          i,         j,        ij,      n,      nn
    integer(4)          ::      ifpar,     jfpar,    istat,    ierr
    character(len=ncf)  ::     cfitem
    character(len=16)   ::     cdate
    character(len=16)   ::     chead(1:64)

    integer, save ::  nfitem(nitem)
#ifdef OPT_IO_COCOMPI
    integer, save :: mpi_fh(nitem)
    integer :: mpi_fh_t
    integer(kind=mpi_offset_kind) :: disp_t
    integer, allocatable :: ibuf(:)
#endif
    integer :: nf
    character(len=ncf)  ::  cfusfc,   cfvsfc
    character(len=ncf)  ::  cftsfc,   cfqsfc,   cfprec,   cfsflx, cfroff
    character(len=ncf)  ::  cfswdw,   cflwdw,   cfpsfc,   cfssfc
    character(len=ncf)  ::  cftref(ntdim), cftdmp(ntdim)
    character(len=ncf)  ::  grid_jra, roff_map
    integer :: imax
    integer, allocatable :: ip(:), jp(:), iq(:), jq(:)
    real(8), allocatable :: wt(:)
    integer, save :: imaxn
    integer, save, allocatable :: ipn(:), jpn(:), iqn(:), jqn(:)
    real(8), save, allocatable :: wtn(:)
    logical, save :: odirect=.false.

    namelist /nmsfbc/  cfusfc, cfvsfc, cftsfc, cfqsfc,                  &
    &                  cfprec, cfsflx, cfroff, cfswdw, cflwdw, cfpsfc, cfssfc,  &
    &                  cftref, cftdmp, grid_jra, roff_map

    data grid_jra /'not-specified'/
    data roff_map /'not-specified'/
    data cfusfc, cfvsfc / 'not-specified', 'not-specified' /
    data cftsfc, cfqsfc / 'not-specified', 'not-specified' /
    data cfprec, cfsflx / 'not-specified', 'not-specified' /
    data cfswdw, cflwdw / 'not-specified', 'not-specified' /
    data cfpsfc, cfssfc / 'not-specified', 'not-specified' /
    data cfroff / 'not-specified' /
    data cftref / ntdim*'not-specified' /
    data cftdmp / ntdim*'not-specified' /


    integer,parameter :: nx0=640, ny0=320
    real(4) :: direct(nx0, ny0)
    real(8), save :: alon(nx0), alat(ny0)
    integer, save :: mask(nx0,ny0)=1

    if ( of ) then
       call rewnml( ifpar, jfpar )
       read(ifpar, nmsfbc, iostat = istat )
       call cstnml( jfpar, 'tmintp', 'nmsfbc', istat )
       of = .false.

#ifdef OPT_IO_COCOMPI
       call mpi_filopn(mpi_fh_t, grid_jra, 'READ')
       disp_t = 0
       call mpi_read_root(alon, nx0, mpi_fh_t, disp_t)
       call mpi_read_root(alat, ny0, mpi_fh_t, disp_t)
       if (.not. mpi_iseof(mpi_fh_t, disp_t)) then
          call mpi_read_root_int(mask, nx0*ny0, mpi_fh_t, disp_t)
       end if
       call mpi_filcls(mpi_fh_t)
#else
       if ( myrank == iroot ) then
          call filopn(nf, grid_jra, 'READ')
          read(nf) alon
          read(nf) alat
          read(nf, end=901) mask
901       continue
          call filcls(nf)
       end if
#endif
       call mpi_bcast(alon,   nx0,     mpi_real8, iroot, mpi_comm_world, ierr)
       call mpi_bcast(alat,   ny0,     mpi_real8, iroot, mpi_comm_world, ierr)
       call mpi_bcast(mask, nx0*ny0, mpi_integer, iroot, mpi_comm_world, ierr)


       if ( trim(roff_map) /= 'not-specified') odirect=.true.
       call rewnml( ifpar, jfpar )
       write(jfpar,*)'runoff odirect=', odirect

       if (odirect) then ! runoff mapping array
#ifdef OPT_IO_COCOMPI
          call mpi_filopn(mpi_fh_t, roff_map, 'READ')
          disp_t = 0
          allocate(ibuf(1))
          call mpi_read_root_int(ibuf, 1, mpi_fh_t, disp_t)
          imax = ibuf(1)
          deallocate(ibuf)
#else
          if ( myrank == iroot ) then
             call filopn(nf, roff_map, 'READ')       
             read(nf)imax
          end if
#endif
          call mpi_bcast(imax,    1, mpi_integer, iroot, mpi_comm_world, ierr) 
          allocate(ip(imax), jp(imax), iq(imax), jq(imax), wt(imax))
          wt(:) = 1.d0
#ifdef OPT_IO_COCOMPI
          call mpi_read_root_int(ip, imax, mpi_fh_t, disp_t)
          call mpi_read_root_int(jp, imax, mpi_fh_t, disp_t)
          call mpi_read_root_int(iq, imax, mpi_fh_t, disp_t)
          call mpi_read_root_int(jq, imax, mpi_fh_t, disp_t)
          if (.not. mpi_iseof(mpi_fh_t, disp_t)) then
             call mpi_read_root(wt, imax, mpi_fh_t, disp_t)
             write(jfpar, *) '*** Weighted distribution of runoff applied. ***'
          end if
          call mpi_filcls(mpi_fh_t)
#else
          if ( myrank == iroot ) then
             read(nf)ip ! runoff points (JRA55-do grid)
             read(nf)jp
             read(nf)iq ! runoff points (model grid)
             read(nf)jq
             read(nf, end=902) wt  ! weight
             write(jfpar, *) '*** Weighted distribution of runoff applied. ***'
902          continue
             call filcls(nf)
          end if
#endif
          call mpi_bcast(  ip, imax, mpi_integer, iroot, mpi_comm_world, ierr) 
          call mpi_bcast(  jp, imax, mpi_integer, iroot, mpi_comm_world, ierr) 
          call mpi_bcast(  iq, imax, mpi_integer, iroot, mpi_comm_world, ierr) 
          call mpi_bcast(  jq, imax, mpi_integer, iroot, mpi_comm_world, ierr) 
          call mpi_bcast(  wt, imax, mpi_real8,   iroot, mpi_comm_world, ierr) 
  
! Make list vector of distribution array for each node
          imaxn = 0
          do n = 1, imax
             i =iq(n) -(nx*irank)
             j =jq(n) -(ny*jrank)
             ij=i+istr-1 + nxdim * (j+jstr-1 -1)
             
             if ((i < 1).or.(i > nx).or.(j < 1).or.(j > ny)) then  ! outside the node
                wt(n) = 0.0d0
             else
                iq(n) = i
                jq(n) = j
!                wt(n) = wt(n)/(dx*dy(ij)*hxt(ij)*hyt(ij) *1.d-4) ! m^3/s -> m/s
                wt(n) = wt(n)/(dx*dy(ij)*hxt(ij)*hyt(ij) *1.d-4)/dwatr ! kg/s -> m/s
                if (wt(n) > 0.0d0) then
                   imaxn = imaxn + 1
                end if
             end if
          end do

          if (imaxn == 0) then  ! no distination grid in the node
             allocate(ipn(1),jpn(1),iqn(1),jqn(1),wtn(1))  ! dummy
             ipn(1) = 1
             jpn(1) = 1
             iqn(1) = 1
             jqn(1) = 1
             wtn(1) = 0.0d0
          else
             allocate(ipn(imaxn),jpn(imaxn),iqn(imaxn),jqn(imaxn),wtn(imaxn))
             nn = 0
             do n = 1, imax
                if (wt(n) > 0.0d0) then
                   nn = nn + 1
                   ipn(nn) = ip(n)
                   jpn(nn) = jp(n)
                   iqn(nn) = iq(n)
                   jqn(nn) = jq(n)
                   wtn(nn) = wt(n)
                end if
             end do
             if (nn /= imaxn) then
                write(jfpar,*) '### failed to make list of runoff vector. ###'
                call mpi_abort(mpi_comm_world, 1, ierr)
             end if
          end if

          deallocate(ip, jp, iq, jq, wt)
       end if

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
       cfitem = cfprec
    else if (iitem == 6) then
       cfitem = cfsflx
    else if (iitem == 7) then
       cfitem = cfswdw
    else if (iitem == 8) then
       cfitem = cflwdw
    else if (iitem == 9) then
       cfitem = cfpsfc
    else if (iitem == 10) then
       cfitem = cfroff
    else if (iitem == 11) then
       cfitem = cfssfc
    else if (iitem > nitem ) then
       write(jfpar, *) '*** TMINTP: NO SUCH ITEM ***'
       call mpi_abort(mpi_comm_world, 1, ierr)
    else if ( mod(iitem, 2) == 0 ) then
       i = (iitem - 10) / 2
       cfitem = cftref(i)
    else
       i = (iitem - 11) / 2
       cfitem = cftdmp(i)
    end if


   if ( ofirst(iitem) ) then
#ifdef OPT_IO_COCOMPI
      call mpi_filopn(mpi_fh(iitem), cfitem, 'READ')
#else
      if ( myrank == iroot ) call filopn(nfitem(iitem), cfitem, 'READ')
#endif
      if (iitem == 10) then
         call read_runoff(time, data)
      else
         call read_dat(time, data)
      end if
      time2(    iitem)=time
      data2(:,:,iitem)=data(:,:)

      if (time >= tt ) then
         call rewnml( ifpar, jfpar )
         write(jfpar,*)'error: 1st time of data > initial date'
         stop
      end if
      ofirst(iitem)=.false.
   end if

   if (tt > time2(iitem) ) then
   do
          time1(          iitem) = time2(          iitem)
          data1(1:nx,1:ny,iitem) = data2(1:nx,1:ny,iitem)

          if (iitem == 10) then
             call read_runoff(time, data)
          else
             call read_dat(time, data)
          end if

          time2(    iitem)=time
          data2(:,:,iitem)=data(:,:)
          if ( time2(iitem) >= tt ) exit
   end do
   end if

!--------- time interpolation
       tintv = time2(iitem) - time1(iitem)
       tintp = (tt - time1(iitem)) / tintv
       tintq = (time2(iitem) - tt) / tintv
       do j = 1, ny
          do i = 1, nx
             ditem(istr+i-1, jstr+j-1) = tintq * data1(i, j, iitem)   &
    &                                  + tintp * data2(i, j, iitem)
          end do
       end do

    contains
      !=========================================================================
      subroutine read_dat(time, dout)
        real(8) :: dout(nx,ny)
        real(8) :: time
        integer :: icread
#ifdef OPT_IO_COCOMPI
        integer (kind=mpi_offset_kind), save :: disp(nitem)
        integer, save :: nrec(nitem)=0
        integer(8) :: isize, iskip

        if ( ofirst(iitem) ) disp(iitem)=0
        call mpi_read_direct(chead, direct, mpi_fh(iitem), disp(iitem), icread)
        read(chead(50), '(i6.6,5i2.2)') (idatet(i), i = 1, 6)
        nrec(iitem)=nrec(iitem)+1
        call rewnml( ifpar, jfpar )
        if (nrec(iitem) == 1) write(jfpar,'(a, 6i6)')'    first record   :',idatet
        if (nrec(iitem) == 2) write(jfpar,'(a, 6i6)')'    skipped to here:',idatet

        if (ofirst(iitem)) then ! skip records
           call cyh2ss( time, idatet )
           iskip=int( (tt-time)/(3.d0*3600.d0) ) ! 3h interval
           isize=disp(iitem)
           disp(iitem)=0 ! rewind
           disp(iitem)=disp(iitem) + iskip*isize
        end if
#else
        if ( myrank == iroot ) then
           read(nfitem(iitem)) chead
           read(nfitem(iitem)) direct
        end if
        call mpi_bcast(chead, 1024, mpi_character, iroot, mpi_comm_world, ierr)     
        call mpi_bcast(direct, nx0*ny0, mpi_real4, iroot, mpi_comm_world, ierr)     
#endif

        call intpsfc(iitem, direct, alon, alat, dout, mask)

        cdate = chead(50)
        read(cdate, '(i6.6,5i2.2)') (idatet(i), i = 1, 6)
        call cyh2ss( time, idatet )
      end subroutine read_dat

      !=========================================================================
      subroutine read_runoff(time,dout)
        real(8), intent(out) :: dout(nx,ny), time
        integer, parameter :: nx0=1440, ny0=720
        character(len=16)  :: chead(1:64)
        real(8) ::   datag(nxg,nyg)
        real(4) :: direct(nx0,ny0)
        integer :: n, i0,j0, i1,j1, ij
#ifdef OPT_IO_COCOMPI
        integer (kind=mpi_offset_kind), save :: disp
        integer :: icread
        integer, save :: nrec=0
        integer(8) :: isize, iskip
#endif
        if ( ofirst(iitem) ) disp=0
        if (odirect) then !=====================================================
#ifdef OPT_IO_COCOMPI
           call mpi_read_direct(chead, direct, mpi_fh(iitem), disp, icread)
           read(chead(50), '(i6.6,5i2.2)') (idatet(i), i = 1, 6)

           nrec=nrec+1
           call rewnml( ifpar, jfpar )
           if (nrec == 1) write(jfpar,'(a, 6i6)')'    first record   :',idatet
           if (nrec == 2) write(jfpar,'(a, 6i6)')'    skipped to here:',idatet

           if (ofirst(iitem)) then ! skip records
              call cyh2ss( time, idatet )
              iskip=int( (tt-time)/(24.d0*3600.d0) ) ! runoff: 24h interval
              isize=disp
              disp=0 ! rewind
              disp=disp + iskip*isize
           end if
#else
           if ( myrank == iroot ) then
              read(nfitem(iitem))chead
              read(nfitem(iitem))direct
           end if
           call mpi_bcast(chead, 1024, mpi_character, iroot, mpi_comm_world, ierr)
           call mpi_bcast(direct, nx0*ny0, mpi_real4, iroot, mpi_comm_world, ierr)
#endif
           dout(:,:)=0.d0
           do n=1,imaxn
              dout(iqn(n),jqn(n)) = dout(iqn(n),jqn(n)) + wtn(n) * direct(ipn(n),jpn(n))
           end do
        else  !==================================================================
#ifdef OPT_IO_COCOMPI
           call mpi_read_chead(chead, mpi_fh(iitem), disp, icread)
           call mpi_read_sfc(dout, mpi_fh(iitem), disp)
#else
           if ( myrank == iroot ) then
              read(nfitem(iitem)) chead
              read(nfitem(iitem)) datag
           end if
           call mpi_bcast(chead, 1024, mpi_character, iroot, mpi_comm_world, ierr)
           call scatter_sfc(dout, datag)
#endif
       end if !==================================================================

       read(chead(50), '(i6.6,5i2.2)') (idatet(i), i = 1, 6)
       call cyh2ss( time, idatet )

      end subroutine read_runoff
  end subroutine tmintp_direct

  subroutine intpsfc(iitem, direct, alon, alat, data1, mask)
  use zocdim
  use zocgrd
  use zocmsk
  use ufile
  implicit none
#include "mpif.h"

  integer :: iitem
  integer, parameter :: nitem = (ntdim-2)*2+10
  real(8) ::  data1(nxy)
  integer,parameter :: nx0=640, ny0=320
  real(4) ::  direct(nx0, ny0)
  real(8) ::  direct0(0:nx0+1, 0:ny0+1)
  real(8)  :: alon(nx0), alat(ny0)

  real(8), save ::  alon0(nitem,0:nx0+1), alat0(nitem,0:ny0+1)

  integer :: ij, ij0,  i, j, i0, j0, ijd
  integer, save ::  ilon0(nitem,nxy), jlat0(nitem,nxy)
  integer :: ilon, jlat

  real(8) :: d1, d2, d3, d4, d12, d34, pi
  integer ::  ifpar,  jfpar
  logical, save :: oerr=.false.

  logical, save :: ofirst(1:nitem)=.true.
  integer :: mask(nx0,ny0)

  real(8), save :: glont0(nitem,nxydim)

  if(ofirst(iitem)) then
     pi = 4.d0 * atan(1.d0)
     do i=1,nx0
        alon0(iitem,i)=alon(i)
     end do

     do j=1,ny0
        alat0(iitem,j)=alat(j)
     end do
     alon0(iitem,0)    =alon(nx0) - 2.0d0*pi
     alon0(iitem,nx0+1)=alon(1)   + 2.0d0*pi

     glont0(iitem, :) = glont(:)
     do ij = 1, nxydim
        if (glont0(iitem, ij) < alon0(iitem,0)) then
           glont0(iitem, ij) = glont0(iitem, ij) + 2.0d0*pi
        end if
        if (glont0(iitem, ij) > alon0(iitem,nx0+1)) then
           glont0(iitem, ij) = glont0(iitem, ij) - 2.0d0*pi
        end if
     end do
     alat0(iitem,0)    = alat0(iitem,1  )-(alat0(iitem,2  )-alat0(iitem,1    ))
     alat0(iitem,ny0+1)= alat0(iitem,ny0)+(alat0(iitem,ny0)-alat0(iitem,ny0-1))
  end if


  do j=1,ny0
     do i=1,nx0
        direct0(i,j)=direct(i,j)
     end do
  end do

  do j = 1, ny0
     direct0(0, j)     = direct0(nx0, j)
     direct0(nx0+1, j) = direct0(1, j)
  end do

  do i = 0, nx0+1
     direct0(i, 0)     = direct0(i, 1)
     direct0(i, ny0+1) = direct0(i, ny0)
  end do


  if(ofirst(iitem)) then 
     do ij0=1,nxy
        ijd = mod(ij0-1,nx)+istr + nxdim*(int((ij0-1)/nx)+jstr-1)
        do i = 1, nx0+1
           if (alon0(iitem,i) .gt. glont0(iitem, ijd)) then
              ilon0(iitem,ij0) = i
              go to 101
           end if
        end do
101     continue
        do j = 1, ny0+1
           if (alat0(iitem,j) .gt. glatt(ijd)) then
              jlat0(iitem,ij0) = j
              go to 102
           end if
        end do
102     continue
     end do
  end if


  do ij0=1,nxy
     ijd = mod(ij0-1,nx)+istr + nxdim*(int((ij0-1)/nx)+jstr-1)
     ilon=ilon0(iitem,ij0)  
     jlat=jlat0(iitem,ij0)  

     d1 = direct0(ilon-1, jlat-1)
     d2 = direct0(ilon, jlat-1)
     d3 = direct0(ilon-1, jlat)
     d4 = direct0(ilon, jlat)

     if( max(d1,d2,d3,d4) .ge. 1.d19) then
!         oerr=.true.
        data1(ij0)=1.d20
     else

        d12 = (  d1 * (alon0(iitem,ilon) - glont0(iitem,ijd))          &
             &       + d2 * (glont0(iitem,ijd) - alon0(iitem,ilon-1))) &
             &      / (alon0(iitem,ilon) - alon0(iitem,ilon-1))

        d34 = (  d3 * (alon0(iitem,ilon) - glont0(iitem,ijd))          &
             &       + d4 * (glont0(iitem,ijd) - alon0(iitem,ilon-1))) &
             &      / (alon0(iitem,ilon) - alon0(iitem,ilon-1))

        data1(ij0)                                                 &
             &      = (  d12 * (alat0(iitem,jlat) - glatt(ijd))    &
             &         + d34 * (glatt(ijd) - alat0(iitem,jlat-1))) &
             &        / (alat0(iitem,jlat) - alat0(iitem,jlat-1))
     end if
  end do


!      if(oerr) then
!         call rewnml(ifpar, jfpar)
!         write(jfpar,*)'error: missing value in surface forcing data'
!         call mpi_abort(mpi_comm_world, 1, ierr)
!      end if

  ofirst(iitem)=.false.

  if(iitem .le. 2) then
     call intp_spline(alon, alat, mask, direct, data1)
  end if

!  return
end subroutine intpsfc

!====================================================================
  subroutine intp_spline(alon, alat, mask, dat4, data1)
  use zocdim,  only  : nx, ny, istr, jstr, nxdim, nxydim
  use zocgrd,  only  : glont, glatt
  implicit none
  integer, parameter :: nx0=640, ny0=320
  real(8) :: alon(nx0),alat(ny0)
  integer :: mask(nx0, ny0)
  real(4) :: dat4(nx0, ny0)
  integer, save :: is(ny0)
  logical, save :: ocycle(ny0)
  logical, save :: of=.true.
  real(8), save :: pi

  real(8) :: dat(0:nx0+1,ny0), alonc(0:nx0+1), y2c(0:nx0+1,1:ny0)
  real(8) :: buf(nx0)
  real(8) :: fy(ny0), fy2(ny0), fyb(ny0,nx), fy2b(ny0,nx), xb(nx)
  real(8) :: data1(nx,ny), tmp(nx,ny)
  real(8) :: x,y,value
  integer :: i,j,j1, ijd

  real(8), save :: glont0(nxydim)

  if (of) then
     ocycle(:)=.true.
     do j=1,ny0
        do i=1,nx0
           if(mask(i,j) ==0 ) ocycle(j)=.false.
        end do
     end do

     is(:)=0
     do j=1,ny0
        if(ocycle(j)) cycle
        if(mask(1,j) == 0) cycle
        do i=nx0,1,-1
           if(mask(i,j) .eq. 0) exit
           is(j)=is(j) -1
        end do
     end do
     pi = 4.d0 * atan(1.d0)

     alonc(1:nx0)=alon(1:nx0)
     alonc(0)=alon(nx0) - 2.d0*pi
     alonc(nx0+1)=2.d0*pi + alon(1)

     glont0(:) = glont(:)
     where(glont0 < alonc(0))     glont0 = glont0 + 2.d0*pi
     where(glont0 > alonc(nx0+1)) glont0 = glont0 - 2.d0*pi

     of=.false.
  end if

  where(mask ==0) dat4=-1.d20 !set missing value



  dat(1:nx0,1:ny0)=dat4(1:nx0,1:ny0)
  dat(   0,:)=dat(nx0,:)
  dat(nx0+1,:)=dat( 1,:)

  do j=1,ny0
     if(ocycle(j)) then
       call spline_cyclic(alonc,dat(:,j),y2c(:,j))
     else   
       !cshift is unnecessary for alon because d(alon) is const.
       buf=cshift(dat(1:nx0,j), is(j))
       call spline(alon, buf, y2c(1:nx0,j))
       y2c(1:nx0,j)=cshift(y2c(1:nx0,j), -is(j))
       y2c(    0,j)=y2c(nx0,j)
       y2c(nx0+1,j)=y2c(  1,j)
     end if
  end do

  do j=1,ny
  do i=1,nx
     ijd = i+istr-1 + nxdim*(j+jstr-2)
     y=glatt(ijd)
     x=glont0(ijd)

     if( j .eq. 1) then ! effective for lat/lon grid
        do j1=1,ny0
           call splint(alonc(0:nx0+1), dat(0:nx0+1,j1), y2c(0:nx0+1,j1), x, fy(j1))
        end do
        call spline(alat, fy, fy2)
          xb(i)=x
         fyb(:,i)= fy(:)
        fy2b(:,i)=fy2(:)

     else
        if( abs(x-xb(i)) .le. 1.e-10 ) then
           fy= fyb(:,i)
          fy2=fy2b(:,i)
        else
           do j1=1,ny0
              call splint(alonc(0:nx0+1), dat(0:nx0+1,j1), y2c(0:nx0+1,j1), x, fy(j1))
           end do
           call spline(alat, fy, fy2)
        end if
     end if     

     call splint(alat, fy, fy2, y, value)

        tmp(i,j)=value
  end do
  end do

  do j=1,ny
  do i=1,nx
     if(tmp(i,j) > -1.d19) data1(i,j)=tmp(i,j)
  end do
  end do


  end subroutine intp_spline
!====================================================================
  subroutine spline(x,y,y2)
  implicit none
  integer :: n, i
  real(8), intent(in ) :: x(:), y(:)
  real(8), intent(out) :: y2(:)
  real(8),allocatable :: a(:), b(:), c(:), r(:)

  n=size(x)
  allocate(a(n), b(n), c(n), r(n))
  do i=2,n-1
    a(i)=(x(i  )-x(i-1))/6.d0
    b(i)=(x(i+1)-x(i-1))/3.d0
    c(i)=(x(i+1)-x(i  ))/6.d0
    r(i)=(y(i+1)-y(i  ))/(x(i+1)-x(i  ))  &
       - (y(i  )-y(i-1))/(x(i  )-x(i-1))

    if( min(y(i-1), y(i+1)) < -1.d19 ) then
     a(i)=0.d0
     c(i)=0.d0
     r(i)=0.d0
    end if
  end do
  y2(1)=0.d0 !natural spline
  y2(n)=0.d0   
  call tridag(a(2:n-1), b(2:n-1), c(2:n-1), r(2:n-1), y2(2:n-1))

  deallocate(a,b,c,r)
  return
  end subroutine spline
!====================================================================
  subroutine spline_cyclic(x,y,y2)
  implicit none
  integer :: n, i
  real(8), intent(in ) :: x(0:), y(0:)
  real(8), intent(out) :: y2(0:)
  real(8) :: alpha, beta
  real(8),allocatable :: a(:), b(:), c(:), r(:)

  n=size(x)-2  
  allocate(a(n), b(n), c(n), r(n))
  do i=1,n
    a(i)=(x(i  )-x(i-1))/6.d0
    b(i)=(x(i+1)-x(i-1))/3.d0
    c(i)=(x(i+1)-x(i  ))/6.d0
    r(i)=(y(i+1)-y(i  ))/(x(i+1)-x(i  ))  &
       - (y(i  )-y(i-1))/(x(i  )-x(i-1))
  end do

  beta =a(1)
  alpha=c(n)
  call tridag_cyclic (a, b, c, alpha, beta, r, y2(1:n))

  y2(0  )=y2(n)
  y2(n+1)=y2(1)
  deallocate(a,b,c,r)
  return
  end subroutine spline_cyclic
!====================================================================
  subroutine splint(x,y,y2, xx, value)
  implicit none
  real(8),intent(in ) :: x(:), y(:), y2(:), xx
  real(8),intent(out) :: value
  real(8) :: a, b, c, d, h
  integer :: nx, i,i0

  nx=size(x)
  i0=-999
  do i=1,nx-1
     if( (xx .ge. x(i)) .and. (xx .le. x(i+1)) ) then
        i0=i
     end if     
  end do
  if (i0 >0) then
     h=(x(i0+1)-x(i0))
     a=(x(i0+1)-xx)/h
     b=1.d0-a
     c=a*(a*a-1)*h*h/6.d0
     d=b*(b*b-1)*h*h/6.d0

     value=a*y(i0)+ b*y(i0+1)+c*y2(i0)+ d*y2(i0+1)

     if ( min(y(i0), y(i0+1)) < -1.d19 ) value=-1.d20 !missing value
  end if

  if(i0 < -1) then
    value= -1.d20 !missing value
  end if


  do i=1,nx
     if( xx .eq. x(i) ) value=y(i)
  end do

  end subroutine splint
!====================================================================
  subroutine tridag(a, b, c, r, u)
  implicit none
  integer :: n, j
  real(8), intent(in ) :: a(:), b(:), c(:), r(:)
  real(8), intent(out) :: u(:)
  real(8) :: bet
  real(8),allocatable :: gam(:)
  n=size(a)

  allocate(gam(n))
  bet=b(1)
  u(1)=r(1)/bet
  do j=2,n
     gam(j)=c(j-1)/bet
     bet=b(j)-a(j)*gam(j)
     if (bet .eq. 0.d0)then
        write(*,*)'tridag failed'
        stop
     end if
        
     u(j)=(r(j)-a(j)*u(j-1))/bet
  end do
  do j=n-1,1,-1
     u(j)=u(j)-gam(j+1)*u(j+1)
  end do
  deallocate(gam)
  return
  end subroutine tridag
!====================================================================
  subroutine tridag_cyclic (a, b, c, alpha, beta, r, x)      
!     From Numerical Recipes p67-68
  implicit none
  real(8), intent(in ) ::  a(:), b(:), c(:), alpha, beta, r(:)
  real(8), intent(out) ::  x(:)
  integer :: n, i
  real(8) :: fact, gamma
  real(8),allocatable :: bb(:), u(:), z(:)

  n=size(a)
  allocate(bb(n))
  allocate( u(n))
  allocate( z(n))

  gamma=-b(1)
  bb(1)=b(1)-gamma
  bb(n)=b(n)-alpha*beta/gamma

  do i=2,n-1
     bb(i)=b(i)
  end do
      
  call tridag(a, bb, c, r, x)

  u(1)=gamma
  u(n)=alpha
  do i=2,n-1
     u(i)=0.d0
  end do
      
  call tridag(a, bb, c, u, z)
  fact=( x(1) + beta*x(n)/gamma ) /(1.d0 +z(1) +beta*z(n)/ gamma )
      
  do i=1,n
     x(i)=x(i)-fact*z(i)
  end do

  deallocate(bb)
  deallocate( u)
  deallocate( z)
  return
  end subroutine tridag_cyclic



#ifndef OPT_IO_COCOMPI
  subroutine tmintp(  ditem, iitem  )

    use zocdim,  only  :                                              &
         nxdim,  nydim, ntdim,                                        &
           nxg,    nyg,    nx,    ny,                                 &
          istr,   jstr
    use zocgrd,  only  :     tt
    use zocfil,  only  :    ncf
    use zocnod,  only  :  iroot,  myrank
    use ufile
    use bgs2d
    use ucaln

    implicit none

#include "mpif.h"

    real(8),    intent(inout)  ::  ditem(nxdim, nydim)
    integer(4), intent(in)     ::  iitem

!---- local variables
    integer(4), parameter      ::  nitem = ( ntdim - 2 ) * 2 + 11 
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
    character(len=ncf)  ::  cftsfc,   cfqsfc,   cfprec,   cfsflx, cfroff
    character(len=ncf)  ::  cfswdw,   cflwdw,   cfpsfc,   cfssfc
    character(len=ncf)  ::  cftref(ntdim), cftdmp(ntdim)
    character(len=ncf)  ::  grid_jra, roff_map

    namelist /nmsfbc/  cfusfc, cfvsfc, cftsfc, cfqsfc,                  &
    &                  cfprec, cfsflx, cfroff, cfswdw, cflwdw, cfpsfc, cfssfc,  &
    &                  cftref, cftdmp, grid_jra, roff_map

    data cfusfc, cfvsfc / 'not-specified', 'not-specified' /
    data cftsfc, cfqsfc / 'not-specified', 'not-specified' /
    data cfprec, cfsflx / 'not-specified', 'not-specified' /
    data cfswdw, cflwdw / 'not-specified', 'not-specified' /
    data cfpsfc, cfssfc / 'not-specified', 'not-specified' /
    data cfroff / 'not-specified' /
    data cftref / ntdim*'not-specified' /
    data cftdmp / ntdim*'not-specified' /

    data osngld / nitem*.false. /
    data ofirst / nitem*.true. /
    data of / .true. /
    data iyskip / nitem*1 /

    if ( of ) then
       call rewnml( ifpar, jfpar )
       read(ifpar, nmsfbc, iostat = istat )
       call cstnml( jfpar, 'tmintp', 'nmsfbc', istat )
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
       cfitem = cfprec
    else if (iitem == 6) then
       cfitem = cfsflx
    else if (iitem == 7) then
       cfitem = cfswdw
    else if (iitem == 8) then
       cfitem = cflwdw
    else if (iitem == 9) then
       cfitem = cfpsfc
    else if (iitem == 10) then
       cfitem = cfroff
    else if (iitem == 11) then
       cfitem = cfssfc
    else if (iitem > nitem ) then
       write(jfpar, *) '*** TMINTP: NO SUCH ITEM ***'
       call mpi_abort(mpi_comm_world, 1, ierr)
    else if ( mod(iitem, 2) == 0 ) then
       i = (iitem - 10) / 2
       cfitem = cftref(i)
    else
       i = (iitem - 11) / 2
       cfitem = cftdmp(i)
    end if

    if ( ofirst(iitem) ) then
       if ( myrank == iroot ) then
          call filopn(nfitem(iitem), cfitem, 'READ')
       end if
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

       if ( myrank == iroot ) then
          rewind(nfitem(iitem))
          read(nfitem(iitem)) chead
          read(nfitem(iitem)) ((datag(i, j), i = 1, nxg), j = 1, nyg)
       end if
       call mpi_bcast(chead, 1024, mpi_character,                     &
    &                 iroot, mpi_comm_world, ierr)
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
    &        (oeof, 1, mpi_logical, iroot, mpi_comm_world, ierr)
       if (oeof) go to 98
       call mpi_bcast(chead, 1024, mpi_character,                     &
    &                 iroot, mpi_comm_world, ierr)
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
    &        (oeof, 1, mpi_logical, iroot, mpi_comm_world, ierr)
       if (oeof) go to 97
       call mpi_bcast(chead, 1024, mpi_character,                     &
    &                 iroot, mpi_comm_world, ierr)
       call scatter_sfc(data2(1, 1, iitem), datag)
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
       if ( myrank == iroot ) then
          rewind(nfitem(iitem))
          read(nfitem(iitem)) chead
          read(nfitem(iitem)) ((datag(i, j), i = 1, nxg), j = 1, nyg)
       end if
       call mpi_bcast(chead, 1024, mpi_character,                     &
    &                 iroot, mpi_comm_world, ierr)
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
       call mpi_abort(mpi_comm_world, 1, ierr)
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
    &           (oeof, 1, mpi_logical, iroot, mpi_comm_world, ierr)
          if (oeof) go to 997
          call mpi_bcast(chead, 1024, mpi_character,                  &
    &                    iroot, mpi_comm_world, ierr)
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
    &           (oeof, 1, mpi_logical, iroot, mpi_comm_world, ierr)
          if ( oeof ) go to 997
          call mpi_bcast(chead, 1024, mpi_character,                  &
    &                    iroot, mpi_comm_world, ierr)
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
    use zocnod,  only  :  iroot,  myrank, newcomm
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
       call rewnml( ifpar, jfpar )
       read(ifpar, nmbody, iostat = istat )
       call cstnml( jfpar, 'tmintb', 'nmbody', istat )
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
       call mpi_abort(mpi_comm_world, 1, ierr)
    end if

    if ( ofirst(iitem) ) then
       if ( myrank == iroot ) then
          call filopn(nfitem(iitem), cfitem, 'READ')
       end if
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
       
       if ( myrank == iroot ) then
          rewind(nfitem(iitem))
          read(nfitem(iitem)) chead
          read(nfitem(iitem)) (((datag(i, j, k), i = 1, nxg),         &
    &                                            j = 1, nyg),         &
    &                                            k = 1, nz)
       end if
       call mpi_bcast(chead, 1024, mpi_character,                     &
    &                  iroot, newcomm, ierr)
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
    &        (oeof, 1, mpi_logical, iroot, newcomm, ierr)
       if (oeof) go to 98
       call mpi_bcast(chead, 1024, mpi_character,                     &
    &                 iroot, newcomm, ierr)
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
    &        (oeof, 1, mpi_logical, iroot, newcomm, ierr)
       if ( oeof ) go to 97
       call mpi_bcast(chead, 1024, mpi_character,                     &
    &                 iroot, newcomm, ierr)
       call scatter_bdy(data2(1, 1, 1, iitem), datag)
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
       if ( myrank == iroot ) then
          rewind(nfitem(iitem))
          read(nfitem(iitem)) chead
          read(nfitem(iitem)) (((datag(i, j, k), i = 1, nxg),         &
    &                                            j = 1, nyg),         &
    &                                            k = 1, nz)
       end if
       call mpi_bcast(chead, 1024, mpi_character,                     &
    &                 iroot, newcomm, ierr)
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
       call mpi_abort(mpi_comm_world, 1, ierr)
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
    &           (oeof, 1, mpi_logical, iroot, newcomm, ierr)
          if (oeof) go to 997
          call mpi_bcast(chead, 1024, mpi_character,                  &
   &                     iroot, newcomm, ierr)
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
    &           (oeof, 1, mpi_logical, iroot, newcomm, ierr)
          if (oeof) go to 997
          call mpi_bcast(chead, 1024, mpi_character,                  &
    &                    iroot, newcomm, ierr)
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

#else
  subroutine tmintp(  ditem, iitem  )

    use zocdim,  only  :                                              &
         nxdim,  nydim, ntdim,                                        &
           nxg,    nyg,    nx,    ny,                                 &
          istr,   jstr
    use zocgrd,  only  :     tt
    use zocfil,  only  :    ncf
    use zocnod,  only  :  iroot,  myrank
    use ufile
    use ucaln
    use mpiio
    implicit none

#include "mpif.h"

    real(8),    intent(inout)  ::  ditem(nxdim, nydim)
    integer(4), intent(in)     ::  iitem

!---- local variables
    integer(4), parameter      ::  nitem = ( ntdim - 2 ) * 2 + 11 
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
    character(len=ncf)  ::  cfusfc,   cfvsfc
    character(len=ncf)  ::  cftsfc,   cfqsfc,   cfprec,   cfsflx, cfroff
    character(len=ncf)  ::  cfswdw,   cflwdw,   cfpsfc,   cfssfc
    character(len=ncf)  ::  cftref(ntdim), cftdmp(ntdim)
    character(len=ncf)  ::  grid_jra, roff_map

    namelist /nmsfbc/  cfusfc, cfvsfc, cftsfc, cfqsfc,                  &
    &                  cfprec, cfsflx, cfroff, cfswdw, cflwdw, cfpsfc, cfssfc,  &
    &                  cftref, cftdmp, grid_jra, roff_map

    data cfusfc, cfvsfc / 'not-specified', 'not-specified' /
    data cftsfc, cfqsfc / 'not-specified', 'not-specified' /
    data cfprec, cfsflx / 'not-specified', 'not-specified' /
    data cfswdw, cflwdw / 'not-specified', 'not-specified' /
    data cfpsfc, cfssfc / 'not-specified', 'not-specified' /
    data cfroff / 'not-specified' /
    data cftref / ntdim*'not-specified' /
    data cftdmp / ntdim*'not-specified' /

    data osngld / nitem*.false. /
    data ofirst / nitem*.true. /
    data of / .true. /
    data iyskip / nitem*1 /

    integer, save :: mpi_fh(nitem)
    integer (kind=mpi_offset_kind), save :: disp(nitem)
    integer :: icread

    if ( of ) then
       call rewnml( ifpar, jfpar )
       read(ifpar, nmsfbc, iostat = istat )
       call cstnml( jfpar, 'tmintp', 'nmsfbc', istat )
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
       cfitem = cfprec
    else if (iitem == 6) then
       cfitem = cfsflx
    else if (iitem == 7) then
       cfitem = cfswdw
    else if (iitem == 8) then
       cfitem = cflwdw
    else if (iitem == 9) then
       cfitem = cfpsfc
    else if (iitem == 10) then
       cfitem = cfroff
    else if (iitem == 11) then
       cfitem = cfssfc
    else if (iitem > nitem ) then
       write(jfpar, *) '*** TMINTP: NO SUCH ITEM ***'
       call mpi_abort(mpi_comm_world, 1, ierr)
    else if ( mod(iitem, 2) == 0 ) then
       i = (iitem - 10) / 2
       cfitem = cftref(i)
    else
       i = (iitem - 11) / 2
       cfitem = cftdmp(i)
    end if

    if ( ofirst(iitem) ) then
       call mpi_filopn(mpi_fh(iitem), cfitem, 'READ')
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
       call mpi_read_chead(chead, mpi_fh(iitem), disp(iitem), icread)
       call mpi_read_sfc(data1(1,1,iitem), mpi_fh(iitem), disp(iitem))
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
       call mpi_read_chead(chead, mpi_fh(iitem), disp(iitem), icread)
       if (icread .ne. 1024) go to 98
       call mpi_read_sfc(data2(1,1,iitem), mpi_fh(iitem), disp(iitem))
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
       call mpi_read_chead(chead, mpi_fh(iitem), disp(iitem), icread)
       if (icread .ne. 1024) go to 97
       call mpi_read_sfc(data2(1,1,iitem), mpi_fh(iitem), disp(iitem))
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
       call mpi_read_chead(chead, mpi_fh(iitem), disp(iitem), icread)
       call mpi_read_sfc(data2(1,1,iitem), mpi_fh(iitem), disp(iitem))
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
       call mpi_abort(mpi_comm_world, 1, ierr)
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
       call mpi_read_chead(chead, mpi_fh(iitem), disp(iitem), icread)
       if(icread .ne. 1024) go to 997
       call mpi_read_sfc(data2(1,1,iitem), mpi_fh(iitem), disp(iitem))
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
       call mpi_read_chead(chead, mpi_fh(iitem), disp(iitem), icread)
       if(icread .ne. 1024) go to 997
       call mpi_read_sfc(data2(1,1,iitem), mpi_fh(iitem), disp(iitem))
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
    use zocnod,  only  :  iroot,  myrank
    use ufile
    use ucaln
    use mpiio
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

    integer, save :: mpi_fh(nitem)
    integer (kind=mpi_offset_kind), save :: disp(nitem)
    integer :: icread

    if ( of ) then
       call rewnml( ifpar, jfpar )
       read(ifpar, nmbody, iostat = istat )
       call cstnml( jfpar, 'tmintb', 'nmbody', istat )
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
       call mpi_abort(mpi_comm_world, 1, ierr)
    end if

    if ( ofirst(iitem) ) then
       call mpi_filopn(mpi_fh(iitem), cfitem, 'READ')
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
       call mpi_read_chead(chead, mpi_fh(iitem), disp(iitem), icread)
       call mpi_read_bdy(data1(1,1,1,iitem), mpi_fh(iitem), disp(iitem))
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
       call mpi_read_chead(chead, mpi_fh(iitem), disp(iitem), icread)
       if(icread .ne. 1024) go to 98
       call mpi_read_bdy(data2(1,1,1,iitem), mpi_fh(iitem), disp(iitem))
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
       call mpi_read_chead(chead, mpi_fh(iitem), disp(iitem), icread)
       if(icread .ne. 1024) go to 97
       call mpi_read_bdy(data2(1,1,1,iitem), mpi_fh(iitem), disp(iitem))
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
       call mpi_read_chead(chead, mpi_fh(iitem), disp(iitem), icread)
       call mpi_read_bdy(data2(1,1,1,iitem), mpi_fh(iitem), disp(iitem))
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
       call mpi_abort(mpi_comm_world, 1, ierr)
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
       call mpi_read_chead(chead, mpi_fh(iitem), disp(iitem), icread)
       if (icread .ne. 1024) goto 997
       call mpi_read_bdy(data2(1,1,1,iitem), mpi_fh(iitem), disp(iitem))
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
       call mpi_read_chead(chead, mpi_fh(iitem), disp(iitem), icread)
       if (icread .ne. 1024) goto 997
       call mpi_read_bdy(data2(1,1,1,iitem), mpi_fh(iitem), disp(iitem))
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
#endif

end module utint


