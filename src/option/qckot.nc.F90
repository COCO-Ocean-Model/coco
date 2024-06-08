module nc_io
  use mpi
  use netcdf
  use zocdim, only : mpi_comm_ogcm
  use zocfil, only : nfomax
  implicit none

  logical, save :: ofirst(nfomax) = .true.
  real(8), save :: time0(nfomax) = 0.d0
  private
  public :: nc_write
contains
  subroutine nc_write(cf, fid, time, &
       & nxg, nyg, nx, ny, nz, irank, jrank, &
       & shuffle, deflate_level, &
       & buf4, cvar, v_unit, v_long_name, fill_value, &
       & x, x_unit, x_long_name, &
       & y, y_unit, y_long_name, &
       & z, z_unit, z_long_name, &
#ifdef OPT_EXMASK
       & x_vert, y_vert, &
#endif
       & z_bnd, &
       & t_unit, time1)

    implicit none
    
    integer, parameter :: nbnd = 2
    integer, parameter :: nvert = 4
    
    character(len=*), intent(in) :: cf
    integer, intent(in) :: fid
    real(8), intent(in) :: time
    real(8), intent(in) :: time1

    integer, intent(in) :: nxg, nyg, nx, ny, nz, irank, jrank

    logical, intent(in) :: shuffle
    integer, intent(in) :: deflate_level

    character(len=*), intent(in) :: cvar, v_unit, v_long_name
    real(4), intent(in) :: buf4(nx*ny*nz)
    real(4), intent(in) :: fill_value

#ifdef OPT_EXMASK
    real(8), intent(in) :: x(nxg, nyg)
#else
    real(8), intent(in) :: x(nxg)
#endif
    character(len=*), intent(in) :: x_unit, x_long_name

#ifdef OPT_EXMASK
    real(8), intent(in) :: y(nxg, nyg)
#else
    real(8), intent(in) :: y(nyg)
#endif
    character(len=*), intent(in) :: y_unit, y_long_name

    real(8), intent(in), optional :: z(nz) ! only 3-D variables
    real(8), intent(in), optional :: z_bnd(nbnd, nz)
    character(len=*), intent(in), optional :: z_unit, z_long_name

    character(len=*), intent(in) :: t_unit

#ifdef OPT_EXMASK
    real(8),intent(in), optional :: x_vert(nvert, nxg, nyg), y_vert(nvert, nxg, nyg) ! only EXMASK
#endif

    real(8) :: t(1), tbnd(nbnd)
    integer :: ncid, ndims, v_id, x_id, y_id, z_id, t_id, tbnd_id, zbnd_id
    integer :: x_dimid, y_dimid, z_dimid, t_dimid, bnd_dimid
    integer :: start_x(1), start_y(1), start_z(1), start_t(1)
    integer :: count_x(1), count_y(1), count_z(1), count_t(1)
    
    integer, allocatable :: v_dimids(:), start_v(:), count_v(:)
    integer ::  zb_dimids(2), start_zb(2),  count_zb(2)
    integer ::  tb_dimids(2), start_tb(2),  count_tb(2)
#ifdef OPT_EXMASK
    integer ::  xy_dimids(2), start_xy(2),  count_xy(2)
    integer :: xyv_dimids(3), start_xyv(3), count_xyv(3)
    integer, parameter :: ntype_index = NF90_INT
    integer :: ix_id, iy_id
    integer, allocatable :: ix(:), iy(:)
    integer :: vert_dimid, x_vert_id, y_vert_id
    integer :: i
#endif
    integer :: nt
    integer, parameter :: ntype_coord = NF90_FLOAT
    integer, parameter :: ntype_out = NF90_FLOAT

    if (ofirst(fid)) then
       ! save initial time
       time0(fid) = time1
       start_t(1)=1
       
       ! create (open) file
       call check(nf90_create(cf, IOR(IOR(NF90_CLOBBER,NF90_NETCDF4),NF90_CLASSIC_MODEL), &
            & ncid, comm = mpi_comm_ogcm, info = MPI_INFO_NULL))

       ! define dimensions
       call check(nf90_def_dim(ncid, "bnds",      nbnd,         bnd_dimid))
#ifdef OPT_EXMASK
       if (present(x_vert)) &
     & call check(nf90_def_dim(ncid, "vertices",  nvert,       vert_dimid))
       call check(nf90_def_dim(ncid, "x",         nxg,            x_dimid))
       call check(nf90_def_dim(ncid, "y",         nyg,            y_dimid))
#else
       call check(nf90_def_dim(ncid, "longitude", nxg,            x_dimid))
       call check(nf90_def_dim(ncid, "latitude",  nyg,            y_dimid))
#endif
       if (present(z)) &
     & call check(nf90_def_dim(ncid, "level",     nz,             z_dimid))
       call check(nf90_def_dim(ncid, "time",      NF90_UNLIMITED, t_dimid))

       ! dimension ID
#ifdef OPT_EXMASK
       xy_dimids = (/ x_dimid, y_dimid /)
       if (present(x_vert)) &
     & xyv_dimids(1:3) = (/ vert_dimid, x_dimid, y_dimid /)
#endif
       if (time /= time1) &
     & tb_dimids(1:2) = (/bnd_dimid, t_dimid/)
       if (present(z)) then
          ndims=4
          allocate( v_dimids(ndims))
          v_dimids(1:ndims) = (/x_dimid, y_dimid, z_dimid, t_dimid/)
          zb_dimids(1:2) = (/bnd_dimid, z_dimid/)
       else
          ndims=3
          allocate( v_dimids(ndims))       
          v_dimids(1:ndims) = (/x_dimid, y_dimid, t_dimid/)
       end if

       ! define variables
#ifdef OPT_EXMASK
       call check(nf90_def_var(ncid, 'longitude',          ntype_coord,  xy_dimids,      x_id ))
       call check(nf90_def_var(ncid, 'latitude',           ntype_coord,  xy_dimids,      y_id ))
       if (present(x_vert)) then
       call check(nf90_def_var(ncid, 'longitude_vertices', ntype_coord, xyv_dimids, x_vert_id ))
       call check(nf90_def_var(ncid, 'latitude_vertices',  ntype_coord, xyv_dimids, y_vert_id ))
       end if
       call check(nf90_def_var(ncid, "x",                  ntype_index,   x_dimid,      ix_id ))
       call check(nf90_def_var(ncid, "y",                  ntype_index,   y_dimid,      iy_id ))
#else
       call check(nf90_def_var(ncid, 'longitude',          ntype_coord,   x_dimid,       x_id ))
       call check(nf90_def_var(ncid, 'latitude',           ntype_coord,   y_dimid,       y_id ))
#endif
       if (present(z)) then
       call check(nf90_def_var(ncid, 'level',              ntype_coord,   z_dimid,       z_id ))
       call check(nf90_def_var(ncid, 'level_bnds',         ntype_coord,  zb_dimids,   zbnd_id ))
       end if
       call check(nf90_def_var(ncid, 'time',               ntype_coord,   t_dimid,       t_id ))
       if (time /= time1) & ! time-averaged output
     & call check(nf90_def_var(ncid, 'time_bnds',          ntype_coord,  tb_dimids,   tbnd_id ))
       call check(nf90_def_var(ncid,  cvar,                ntype_out,     v_dimids,      v_id, &
            & shuffle=shuffle, deflate_level=deflate_level ))

       ! put attributes (unit, bounds)
       call check(nf90_put_att(ncid,      x_id, 'units', x_unit))
       call check(nf90_put_att(ncid,      y_id, 'units', y_unit))
#ifdef OPT_EXMASK
       if (present(x_vert)) then
       call check(nf90_put_att(ncid, x_vert_id,  'units', x_unit))
       call check(nf90_put_att(ncid, y_vert_id,  'units', y_unit))
       call check(nf90_put_att(ncid,      x_id, 'bounds', 'longitude_vertices'))
       call check(nf90_put_att(ncid,      y_id, 'bounds',  'latitude_vertices'))
       end if
       call check(nf90_put_att(ncid,     ix_id,  'units', '1'))
       call check(nf90_put_att(ncid,     iy_id,  'units', '1'))
#endif
       if (present(z)) then
       call check(nf90_put_att(ncid,      z_id,  'units', z_unit))
       call check(nf90_put_att(ncid,   zbnd_id,  'units', z_unit))
       call check(nf90_put_att(ncid,      z_id, 'bounds', 'level_bnds'))
       end if
       call check(nf90_put_att(ncid,      t_id,  'units', t_unit))
       if (time /= time1) then ! time-averaged output
       call check(nf90_put_att(ncid,   tbnd_id,  'units', t_unit))
       call check(nf90_put_att(ncid,      t_id, 'bounds', 'time_bnds'))
       end if
       call check(nf90_put_att(ncid,      v_id,  'units', v_unit))

       ! put attributes (names, axis, etc ...)
#ifdef OPT_EXMASK
       call check(nf90_put_att(ncid, ix_id, 'long_name', 'cell index along first dimension'))
       call check(nf90_put_att(ncid, iy_id, 'long_name', 'cell index along second dimension' ))
       call check(nf90_put_att(ncid, ix_id, 'axis', 'X' ))
       call check(nf90_put_att(ncid, iy_id, 'axis', 'Y' ))
#endif
       call check(nf90_put_att(ncid,  x_id, 'long_name', x_long_name ))
       call check(nf90_put_att(ncid,  y_id, 'long_name', y_long_name ))
       call check(nf90_put_att(ncid,  x_id, 'standard_name', 'longitude'))
       call check(nf90_put_att(ncid,  y_id, 'standard_name', 'latitude' ))
       call check(nf90_put_att(ncid,  x_id, '_CoordinateAxisType', 'Lon'))
       call check(nf90_put_att(ncid,  y_id, '_CoordinateAxisType', 'Lat'))
       if( present(z) ) &
     & call check(nf90_put_att(ncid,  z_id, 'long_name', z_long_name ))
       call check(nf90_put_att(ncid,  t_id, 'long_name', 'time'      ))
       call check(nf90_put_att(ncid,  v_id, 'long_name', v_long_name ))
       if( present(z) ) then
          call check(nf90_put_att(ncid,  v_id, 'coordinates', 'time, level, latitude, longitude' )) ! for CDO remapping
       else
          call check(nf90_put_att(ncid,  v_id, 'coordinates', 'time, latitude, longitude' )) ! for CDO remapping
       end if
       call check(nf90_put_att(ncid,  v_id, '_FillValue', fill_value ))
       
       call check(nf90_enddef(ncid))

       ! put variables
       if (irank+jrank == 0) then ! coordinate info is written only root node
          ! define start and count
          start_x(1)=  1;  start_y(1)= 1; start_z(1)= 1
          count_x(1)=nxg; count_y(1)=nyg; count_z(1)=nz
          start_zb  = (/     1, start_z(1) /)
          count_zb  = (/  nbnd, count_z(1) /)
#ifdef OPT_EXMASK
          start_xy  = (/        start_x(1), start_y(1) /)
          count_xy  = (/        count_x(1), count_y(1) /)
          start_xyv = (/     1, start_x(1), start_y(1) /)
          count_xyv = (/ nvert, count_x(1), count_y(1) /)
          allocate(ix(1:nxg))
          allocate(iy(1:nyg))
          do i = 1, nxg
             ix(i) = i
          end do
          do i = 1, nyg
             iy(i) = i
          end do
          
          call check(nf90_put_var(ncid,      x_id,  x,      start=start_xy,  count=count_xy ))
          call check(nf90_put_var(ncid,      y_id,  y,      start=start_xy,  count=count_xy ))
          call check(nf90_put_var(ncid,     ix_id, ix,      start=start_x,   count=count_x  ))
          call check(nf90_put_var(ncid,     iy_id, iy,      start=start_y,   count=count_y  ))
          if (present(x_vert)) then
          call check(nf90_put_var(ncid, x_vert_id,  x_vert, start=start_xyv, count=count_xyv))
          call check(nf90_put_var(ncid, y_vert_id,  y_vert, start=start_xyv, count=count_xyv))
          end if
          deallocate(ix, iy)
#else
          call check(nf90_put_var(ncid,      x_id,  x,      start=start_x,   count=count_x  ))
          call check(nf90_put_var(ncid,      y_id,  y,      start=start_y,   count=count_y  ))
#endif
          if( present(z) ) then
          call check(nf90_put_var(ncid,      z_id,  z,      start=start_z,   count=count_z  ))
          call check(nf90_put_var(ncid,   zbnd_id,  z_bnd,  start=start_zb,  count=count_zb ))
          end if
       end if
       
       deallocate(v_dimids)
       ofirst(fid) = .false.
    else
       ! open file
       call check(nf90_open_par(cf, nf90_write, comm = mpi_comm_ogcm, info = MPI_INFO_NULL, ncid = ncid))

       ! get current record length
       call check(nf90_inq_dimid(ncid, 'time', t_dimid))
       call check(nf90_inquire_dimension(ncid, t_dimid, len=nt))
       start_t(1) = nt+1

       ! get veriable ID
       call check(nf90_inq_varid(ncid, 'time',         t_id))
       if( time /= time1 ) &  ! time-averaged output
    &  call check(nf90_inq_varid(ncid, 'time_bnds', tbnd_id))
       call check(nf90_inq_varid(ncid,   cvar,         v_id))
    end if

    ! write time
    count_t(1) = 1
    if (t_unit(1:1) == 'D' .or. t_unit(1:1) == 'd') then
       t(1) = (time-time0(fid)) / 86400.d0 ! day
    else
       t(1) = (time-time0(fid)) / 3600.d0 ! hour
    end if
    call check(nf90_var_par_access(ncid, t_id, nf90_collective))
    call check(nf90_put_var(ncid, t_id, t, start=start_t, count=count_t))
    if( time /= time1 )then  ! time-averaged output
       tbnd(1) = (time1-time0(fid))/3600.d0
       tbnd(2) = (time*2-time1-time0(fid))/3600.d0
       start_tb(1:2) = (/   1, start_t(1)/)
       count_tb(1:2) = (/nbnd, count_t(1)/)
       call check(nf90_var_par_access(ncid, tbnd_id, nf90_collective))
       call check(nf90_put_var(ncid, tbnd_id, tbnd, start=start_tb, count=count_tb))
    end if

    ! write data body
    if( present(z) ) then
       ndims=4
       allocate(start_v(ndims))
       allocate(count_v(ndims))
       start_v(1:ndims)= (/ irank*nx+1, jrank*ny+1, 1, start_t(1)/)
       count_v(1:ndims)= (/nx, ny, nz, 1 /)
    else
       ndims=3
       allocate(start_v(ndims))
       allocate(count_v(ndims))
       start_v(1:ndims)= (/ irank*nx+1, jrank*ny+1, start_t(1)/)
       count_v(1:ndims)= (/nx, ny, 1 /)
    end if
    call check(nf90_var_par_access(ncid, v_id, nf90_collective))
    call check(nf90_put_var(ncid, v_id, buf4, start=start_v, count=count_v))
    deallocate(start_v)
    deallocate(count_v)

    ! close file
    call check(nf90_close(ncid))
    
  end subroutine nc_write

  subroutine check(status)
    implicit none
    integer, intent(in)           :: status
    if(status /= nf90_noerr) then 
       write(*,*) trim(nf90_strerror(status))
       stop
    end if
  end subroutine check
end module nc_io


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
!     '21.03.05  Y.Komuro: Sigma coordinate output
!     '21.06.05  H.Tatebe: MR8/MR4 support
!
! ---------------------------------------------------------------------
  use zocdim, only : nxydim, nzdim, mpi_comm_ogcm
  use zocfil, only : nfomax, ncf
  use zocout

  implicit none
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
  character, save ::  cdfmt(nfomax)*16
  integer, save :: nvcord(nfomax), nhcord(nfomax), isvint(nfomax)
  data citem  / nfomax*'                ' /
  data cntavr / nfomax*0.d0 /
  data iflout / nfomax*0 /
  data ioavrg / nfomax*0 /

  integer, save :: nbtnb(nxydim, nchmax)
  real*8, save :: dzmnb(nxydim, nzdim, nchmax)

  real(8), allocatable, save ::  dmskt (:,:,:), dmsktl(:,:)
  real(8), allocatable, save ::  dmskv (:,:,:), dmskvl(:,:) 
  real(8), allocatable, save ::  dmskxl(:,:), dmskyl(:,:)
  real(8), allocatable, save ::    buf3(:,:,:)
  real(8), allocatable, save ::    buf2(:,:)
  real(8), allocatable, save ::  dmskl (:,:)
  
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

  integer, save   :: isingl(nfomax)
  data isingl / nfomax*1 /

contains
  subroutine chkset
    use zocdim, only : nxg, nyg, nx, ny, nz, myrank, iroot, &
                       nxgdim, nygdim, nzdim, &
                       nxdim,  nydim,         &
                       igstr, jgstr, kstr,    &
                       istr,  jstr,  kend
    use zocfil, only : nfomax
    use zocmsk, only : amskt, amskv, amftx, amfty, &
#ifdef OPT_BBL
         & amsktb, amskt1, &
#endif
         & nbot
    use bgs3d
    use bshft
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
    character(16) :: ddfmt = 'not-specified'
    character(16) ::  dfmt

    namelist /nmtime/ itstrt, itend, tmstp, iutstp, ntsplt
    namelist /nmdout/ iodstr, iodend, iodint, iudint, iodavr, iodsng, iodsvi, &
         &             ddfmt
    namelist /nmhist/ cohitm, cohfil, cohvco,                                 &
         &            iohstr, iohend, iohint, iuhint, iohavr, iohsng,         &
         &            ioxstr, ioxend, ioystr, ioyend, iozstr, iozend, iosvin, &
         &              dfmt
    namelist /nmrun/ crun
    
    integer :: i, j, k, n, ij
    integer :: istat

    nsig(0) = nz
    nsig(1:nncmax) = 0
    call csgset
    nsnzmx = maxval(nsig)
    call gs3dst_sig( nsnzmx )
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
       isvint(iitem) = iodsvi
       isingl(iitem) = iodsng
       if ( iodsng == 0 ) cdfmt(iitem) = 'UR8'
       if ( iodsng == 1 ) cdfmt(iitem) = 'UR4'
       if ( ddfmt(1:13) /= 'not-specified' ) cdfmt(iitem) = ddfmt
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
       dfmt   = 'not-specified'

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
             if ( isingl(iohitm) == 0 ) cdfmt(iohitm) = 'UR8'
             if ( isingl(iohitm) == 1 ) cdfmt(iohitm) = 'UR4'
          end if
          if ( dfmt(1:13) /= 'not-specified' ) then
             cdfmt(iohitm) = dfmt
          end if
          if ( cdfmt(iohitm)(1:2) /= 'NC' ) then
             if (myrank == iroot) then
                call filopn(nfunit(iohitm), cohfil, 'WRITE')
             end if
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
    allocate(dbloug(1, 1, 1))
    allocate(sngoug(1, 1, 1))
    if (myrank == iroot) then
       deallocate(dbloug, sngoug)
       allocate(dbloug(nxg, nyg, nsnzmx))
       allocate(sngoug(nxg, nyg, nsnzmx))
    end if

!---- following is from brdge.oms.F of MIROC
!---- setting local mask
    allocate(dmsktl(nxydim, nzdim))
    allocate(dmskvl(nxydim, nzdim))
    allocate(dmskxl(nxydim, nzdim))
    allocate(dmskyl(nxydim, nzdim))
    allocate(dmskl (nxydim, nzdim))
    dmsktl(:,:) = amskt(:,:)
#ifdef OPT_BBL
    do ij = 1, nxydim
       dmsktl(ij, kend    ) = amsktb(ij)
       dmsktl(ij, nbot(ij)) = amskt1(ij)
    end do
#endif
#ifdef OPT_TRIPOLE
    call shift1( dmsktl, nxdim, nydim, nzdim, 1.d0, 0, 0 )
#else
    call shift1( dmsktl, nxdim, nydim, nzdim )
#endif
    dmskvl = 1.d0
    do k = 1, nzdim
       do j = 1, ny
          do i = 1, nx
             ij = ( j + jstr - 2 ) * nxdim + i + istr - 1
             dmskvl(ij, k) =                                  &
                  ( 1.d0 - ( 1.d0 - dmsktl(ij        , k) )   &
                          *( 1.d0 - dmsktl(ij+1      , k) )   &
                          *( 1.d0 - dmsktl(ij+nxdim  , k) )   &
                          *( 1.d0 - dmsktl(ij+nxdim+1, k) ) )
          end do
       end do
    end do
#ifdef OPT_TRIPOLE
    call shift1( dmskvl, nxdim, nydim, nzdim, 1.d0, -1, -1 )
#else
    call shift1( dmskvl, nxdim, nydim, nzdim )
#endif
    dmskxl = 0.d0
    dmskyl = 0.d0
    do k = 1, nzdim
       do ij = 1, nxydim
          if (amftx(ij, k) /= 0) then
             dmskxl(ij,k) = 1.d0
          end if
          if (amfty(ij, k) /= 0) then
             dmskyl(ij,k) = 1.d0
          end if
       end do
    end do
    
!---- setting global mask
    allocate(dmskt(1, 1, 1))
    allocate(dmskv(1, 1, 1))
    allocate( buf3(1, 1, 1))
    allocate( buf2(1, 1))
    if (myrank == iroot) then
       deallocate(dmskt, dmskv, buf3, buf2)
       allocate(dmskt(nxg, nyg, nz))
       allocate(dmskv(nxg, nyg, nz))
       allocate(buf3(nxgdim, nygdim, nzdim))
       allocate(buf2(nx, nz))
    end if
    call gather_3d( buf3, dmsktl )
    if ( myrank == iroot ) then
       do k = 1, nz
          do j = 1, nyg
             do i = 1, nxg
                dmskt(i, j, k) = buf3(i+igstr-1, j+jgstr-1, k+kstr-1)
             end do
          end do
       end do
    end if
    call gather_3d( buf3, dmskvl )
    if ( myrank == iroot ) then
       do k = 1, nz
          do j = 1, nyg
             do i = 1, nxg
                dmskv(i, j, k) = buf3(i+igstr-1, j+jgstr-1, k+kstr-1)
             end do
          end do
       end do
    end if
    deallocate ( buf3 )

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
    use zocdim, only : nx, ny, nxy, nxg, nyg, nxdim, nxydim, nxgdim, &
         &  istr, jstr, myrank, ijnode, iroot, &
         &  nic, nz, kstr
    use zocfil, only : nfomax
    use zocgrd, only : tt, nt
    use zocout, only : dbleou, dbloug, sngoug, wrkout
    use bgs3d
    use ucaln
    use ufile
    implicit none

    logical, intent(in) :: oflout(nfomax)
    integer :: iitem
    integer ::  ixdim,  jydim,  kzdim
    integer ::  nsiz
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
    integer :: ij

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
             nsiz = ixdim*jydim*kzdim
             write(chead(64), '(i16)') nsiz
             
             write(chead(39), '(e16.7)') dundef
             chead(40) = chead(39)
             chead(41) = chead(39)
             chead(42) = chead(39)
             chead(43) = chead(39)
             write(chead(44), '(i16)') 1
             write(chead(46), '(i16)') 0
             write(chead(47), '(e16.7)') 0.d0
             
             if (ioavrg(iitem) == 1) then
                call css2yh( idate, ttold(iitem) )
                write(chead(48), '(i4.4,2i2.2,1x,3i2.2,1x)') idate(1:6)
             else
                call css2yh( idate, tt )
                write(chead(48), '(i4.4,2i2.2,1x,3i2.2,1x)') idate(1:6)
             end if
             call css2yh( idate, tt )
             write(chead(49), '(i4.4,2i2.2,1x,3i2.2,1x)') idate(1:6)
             
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
                               dbleou(ijk) = dundef   !! missing
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
             if ( cdfmt(iitem)(1:2) /= 'NC' ) then
             call gather_chk_sig(dbloug, dbleou, nsnzmx)
             if (myrank == iroot) then
                if ( cdfmt(iitem)(1:3) == 'UR8' ) then
                   chead(38) = 'UR8'
                   write(nfunit(iitem)) chead
                   call gfwrt8( dbloug(ixstr(iitem):ixend(iitem),  &
                                       jystr(iitem):jyend(iitem),  &
                                       kzstr(iitem):kzend(iitem)), & 
                                nfunit(iitem),                     &
                                nsiz )              
                end if
                if ( cdfmt(iitem)(1:3) == 'UR4' ) then
                   chead(38) = 'UR4'
                   write(nfunit(iitem)) chead
                   call gfwrt4( dbloug(ixstr(iitem):ixend(iitem),  &
                                       jystr(iitem):jyend(iitem),  &
                                       kzstr(iitem):kzend(iitem)), & 
                                nfunit(iitem),                     &
                                nsiz )
                end if
!----- buffering fty at the southern most part
                if ( chead(32)(1:9) == 'OCLATTPVS' ) then
                   do k = 1, nz
                      do i = 1, nxg
                         buf2(i,k) = dbloug(i,1,k)
                      end do
                   end do
                end if
                if ( nvcord(iitem) == 0 ) then
                   if ( clas(iitem)(1:6) == 'OCLVTT' .or. &
                        clas(iitem)(1:6) == 'OCLVMT' .or. &
                        clas(iitem)(1:6) == 'OCSFCT'       ) then
                      do k = kzstr(iitem), kzend(iitem)
                         do j = jystr(iitem), jyend(iitem)
                            do i = ixstr(iitem), ixend(iitem)
                               if ( dmskt(i,j,k) == 0.d0 ) then
                                  dbloug(i,j,k) = dundef
                               end if
                            end do
                         end do
                      end do
                   end if
                   if ( clas(iitem)(1:6) == 'OCLVTV' .or. &
                        clas(iitem)(1:6) == 'OCLVMV' .or. &
                        clas(iitem)(1:6) == 'OCSFCV'       ) then
                      do k = kzstr(iitem), kzend(iitem)
                         do j = jystr(iitem), jyend(iitem)
                            do i = ixstr(iitem), ixend(iitem)
                               if ( dmskv(i,j,k) == 0.d0 ) then
                                  dbloug(i,j,k) = dundef
                               end if
                            end do
                         end do
                      end do
                   end if
                   if ( clas(iitem)(1:6) == 'OCICET' ) then
                      do k = kzstr(iitem), kzend(iitem)
                         do j = jystr(iitem), jyend(iitem)
                            do i = ixstr(iitem), ixend(iitem)
                               if ( dmskt(i,j,1) == 0.d0 ) then
                                  dbloug(i,j,k) = dundef
                               end if
                            end do
                         end do
                      end do
                   end if
                end if
                if ( chead(32)(1:9) == 'OCLATTPVS' ) then
                   do k = 1, nz
                      do i = 1, nxg
                         dbloug(i,1,k) = buf2(i,k)
                      end do
                   end do
                end if
                if ( cdfmt(iitem)(1:3) == 'MR8' ) then
                   chead(38) = 'MR8'
                   write(nfunit(iitem)) chead
                   call gfwrm8( dbloug(ixstr(iitem):ixend(iitem),  &
                                       jystr(iitem):jyend(iitem),  &
                                       kzstr(iitem):kzend(iitem)), & 
                                nfunit(iitem),                     &
                                nsiz, dundef )
                end if
                if ( cdfmt(iitem)(1:3) == 'MR4' ) then
                   chead(38) = 'MR4'
                   write(nfunit(iitem)) chead
                   call gfwrm4( dbloug(ixstr(iitem):ixend(iitem),  &
                                       jystr(iitem):jyend(iitem),  &
                                       kzstr(iitem):kzend(iitem)), & 
                                nfunit(iitem),                     &
                                nsiz, dundef )
                end if
             end if
             else ! NetCDF
                if ( nvcord(iitem) == 0 ) then
                   select case (clas(iitem)(6:6))
                   case ('T')
                      dmskl = dmsktl
                   case ('V')
                      dmskl = dmskvl
                   case ('X')
                      dmskl = dmskxl
                   case ('Y')
                      dmskl = dmskyl
                   end select
                   if (clas(iitem)(3:5) == 'ICE') then
                      do k = kzstr(iitem), kzend(iitem)
                         dmskl(:, k) = dmskl(:, kstr)
                      end do
                   end if
                   do k = kzstr(iitem), kzend(iitem)
                   do j = 1, ny
                   do i = 1, nx
                      ijk = (k - kzstr(iitem)) * nxy + (j - 1) * nx + i
                      ij = i+istr-1+(j+jstr-2)*nxdim
                      if ( dmskl(ij,k+kstr-1) == 0.d0 ) then
                         dbleou(ijk) = dundef
                      end if
                   end do
                   end do
                   end do
                end if
                if (osingl(iitem)) then
                   do k = kzstr(iitem), kzend(iitem)
                   do j = 1, ny
                   do i = 1, nx
                      ijk = (k - kzstr(iitem)) * nxy + (j - 1) * nx + i
                      snglou(ijk) = dbleou(ijk)
                   end do
                   end do
                   end do
                   call write_ncfile(iitem, chead, nx, ny, kzdim, snglou)
                else
                   call rewnml(ifpar, jfpar)
                   write(jfpar,*)'Err. Only for REAL4'
                   stop
                end if
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
#ifdef OPT_TRIPOLE
    use zocdim, only : nxydim, nxdim, nydim, nxg, nyg, kstr, nic
    use bshft
#else
    use zocdim, only : nxydim, nxg, nyg, kstr, nic
#endif
    use zocout, only : nwork, wrkout 
    use zocfil, only : nfstdo
    use ufile
    
    integer,      intent(in) :: nxitem, nyitem, nzitem, nditem 
    real(8),      intent(in) :: ditem(nditem)
    character(*), intent(in) :: ccitem,  cclas,  htitl,  hunit

    real(8), allocatable :: sigitm(:), sigthk(:)
    logical, allocatable :: osigex(:)
    integer :: sdim
#ifdef OPT_TRIPOLE
    real(8), allocatable :: flxitm(:)
    integer :: kdim
#endif

    integer :: item, iohitm
    integer, save :: jtopad = 1, jtopas = 1
    character ::  ctitem*16
    integer :: i
    logical, save   ::  ofirst(nfomax)
    data ofirst / nfomax*.true. /

    do iohitm = 1, nohitm
       ctitem = ccitem
       if (citem(iohitm) .ne. ccitem) then
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
          if ( isvint(item) > 0 .and. nvcord(item) > 0 .and. citem(item)(1:5) /= 'DZSIG') then
             cunit(item) = TRIM(hunit)//'*cm'
          end if
!         Horizontal position: 1:T, 2:V, 3:X(flux), 4:Y(flux)
          if (cclas(6:6) == 'T') then
             nhcord(item) = 1
          elseif (cclas(6:6) == 'V') then
             nhcord(item) = 2
          elseif (cclas(6:6) == 'X') then
             nhcord(item) = 3
          elseif (cclas(6:6) == 'Y') then
             nhcord(item) = 4
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
!!              '23.08.29:
!!              cyclic output of fy on sigma co. is not implemented yet.
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
#ifdef OPT_TRIPOLE
                if (nhcord(item) == 4) then  !! y-flux
                   allocate(flxitm(nditem))
                   flxitm(:) = ditem(:)
                   if (cclas(3:5) == 'SFC') then
                      kdim = 1
                   else if (cclas(3:5) == 'ICE') then
                      kdim = nic + 1
                   else if (cclas(3:5) == 'LVT' .or. cclas(3:5) == 'LVM') then
                      kdim = nzdim
                   end if
                   call shiftf1( &
                     &          flxitm, &
                     &           nxdim,  nydim,   kdim)
                   do i = 1, nszitm(item)
                      wrkout(itopad(item) + i - 1) = &
                        &    wrkout(itopad(item) + i - 1) * fctavr(item) &
                        &  + flxitm(ktopad(item) + i - 1)
                   end do
                   deallocate(flxitm)
                else
                   do i = 1, nszitm(item)
                      wrkout(itopad(item) + i - 1) = &
                        &    wrkout(itopad(item) + i - 1) * fctavr(item) &
                        &  + ditem(ktopad(item) + i - 1)
                   end do
                endif
#else
                do i = 1, nszitm(item)
                   wrkout(itopad(item) + i - 1) = &
                     &    wrkout(itopad(item) + i - 1) * fctavr(item) &
                     &  + ditem(ktopad(item) + i - 1)
                end do
#endif
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
      &  amsktb, amskvb,  nbotv, &
#endif
      &    nbot
    use zocfil, only: &
      &     ncf
    use ufile
    implicit none

#include "mpif.h"

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
       if (oexist) then
          call filopn( nfsgco, cfsgco, 'read' )
          rewind( nfsgco )
          read( nfsgco ) nnc
       else
          nnc = 0
          write(jfpar, *) &
           &  '*** Sigma-coordinate setting not found.', &
           &  ' output sigma0 only. ***'
       end if
    end if
    call mpi_bcast( nnc, 1, mpi_integer, &
      &             iroot, mpi_comm_ogcm, ierr )
    if (nnc > 0) then
       if ( myrank == iroot ) then
          read( nfsgco ) cname(1:nnc)
          read( nfsgco ) zref(1:nnc)
          read( nfsgco ) nsig(1:nnc)
       end if
       call mpi_bcast( cname, 16*nncmax, mpi_character, &
         &             iroot, mpi_comm_ogcm, ierr )
       call mpi_bcast( zref, nncmax+1, mpi_real8, &
         &             iroot, mpi_comm_ogcm, ierr )
       call mpi_bcast( nsig, nncmax+1, mpi_integer, &
         &             iroot, mpi_comm_ogcm, ierr )
       nsigmx=maxval(nsig)
       allocate(lsig(nsigmx,nnc),lsigp(1:nsigmx+1,nnc),dsig(nsigmx,nnc))
       if ( myrank == iroot ) then
          read( nfsgco ) lsig
          read( nfsgco ) lsigp
          read( nfsgco ) dsig
          call filcls( nfsgco )
       end if
       call mpi_bcast( lsig, nsigmx*nnc, mpi_real8, &
         &             iroot, mpi_comm_ogcm, ierr )
       call mpi_bcast( lsigp, (nsigmx+1)*nnc, mpi_real8, &
         &             iroot, mpi_comm_ogcm, ierr )
       call mpi_bcast( dsig, nsigmx*nnc, mpi_real8, &
         &             iroot, mpi_comm_ogcm, ierr )
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

! '23.08.29: will be updated for output of flux vars on sigma co..
!   nbtnb(ij,3): nbot without BBL at FX-point, ***temporary*** 
    nbtnb(:,3) = nbtnb(:,1)

!   nbtnb(ij,4): nbot without BBL at FY-point, ***temporary*** 
    nbtnb(:,4) = nbtnb(:,1)

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

! '23.08.29: will be updated for output of flux vars on sigma co..
!   dzmnb(ij,3): dzm without BBL at FX-point, ***temporary*** 
    dzmnb(:,:,3) = dzmnb(:,:,1)

!   dzmnb(ij,4): dzm without BBL at FY-point, ***temporary*** 
    dzmnb(:,:,4) = dzmnb(:,:,1)

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
      &  amsktb, &
#endif
      &    nbot,  amskt

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
! '23.08.29: will be updated for output of flux vars on sigma co..
       else if (nch == 3) then  !! sigma for FX-point, ***temporary***
          sigma(:,:,:,3) = sigma(:,:,:,1)
       else if (nch == 4) then  !! sigma for FY-point, ***temporary***
          sigma(:,:,:,4) = sigma(:,:,:,1)
       endif
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
!             write(0,*) ij,sigma(ij,kstr,n,nch)
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
!                write(0,*) i,j,sigma(ij,kend,n,nch)
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
  
  subroutine write_ncfile(iitem, chead, nx, ny, kzdim, dat)
    use zocdim, only : nxg, nyg, nz, myrank, iroot, irank, jrank, ierr, mpi_comm_ogcm, &
         & kstr, nic
    use nc_io
    use ucaln

    implicit none
#include "mpif.h"
    integer,            intent(in) :: iitem
    character (len=16), intent(in) :: chead(1:64)
    integer,            intent(in) :: nx, ny, kzdim
    real(4),            intent(in) :: dat(nx*ny*kzdim)
    character (len=4)  :: cyr
    character (len=2)  :: cmon, cday, chr, cmin, csec
    character (len=5)  :: cz_unit
    character (len=64) :: cz_long_name
    logical, save :: of = .true.
#ifdef OPT_EXMASK
    real(8), save :: xt(nxg,nyg), yt(nxg,nyg)
    real(8), save :: xv(nxg,nyg), yv(nxg,nyg)
    real(8), save :: xx(nxg,nyg), yx(nxg,nyg)
    real(8), save :: xy(nxg,nyg), yy(nxg,nyg)
    real(8) :: x(nxg,nyg), y(nxg,nyg)
    real(8) :: lon0
    integer, parameter :: nvert = 4
    real(8), save :: xt_vert(nvert,nxg,nyg), yt_vert(nvert,nxg,nyg)
    real(8), save :: xv_vert(nvert,nxg,nyg), yv_vert(nvert,nxg,nyg)
    real(8), save :: xx_vert(nvert,nxg,nyg), yx_vert(nvert,nxg,nyg)
    real(8), save :: xy_vert(nvert,nxg,nyg), yy_vert(nvert,nxg,nyg)
    real(8) :: x_vert(nvert,nxg,nyg), y_vert(nvert,nxg,nyg)
#else
    real(8), save :: xt(nxg), yt(nyg)
    real(8), save :: xv(nxg), yv(nyg)
    real(8), save :: xx(nxg), yx(nyg)
    real(8), save :: xy(nxg), yy(nyg)
    real(8) :: x(nxg), y(nyg)
#endif
    integer, parameter :: nbnd = 2
    real(8), save :: zt(nz), zw(nz), zi(nic)
    real(8), save :: zt_bnd(nbnd,nz), zw_bnd(nbnd,nz), zi_bnd(nbnd,nic)
    real(8), allocatable :: z(:), z_bnd(:,:)
    real(8) :: tout, time1
    integer :: idate(6)

    if (of) then
       call def_coord
       of=.false.
    end if

    select case (clas(iitem)(6:6))
    case ('V')
       x=xv
       y=yv
#ifdef OPT_EXMASK
       x_vert=xv_vert
       y_vert=yv_vert
#endif
    case ('X')
       x=xx
       y=yx
#ifdef OPT_EXMASK
       x_vert=xx_vert
       y_vert=yx_vert
#endif
    case ('Y')
       x=xy
       y=yy
#ifdef OPT_EXMASK
       x_vert=xy_vert
       y_vert=yy_vert
#endif
    case default
       x=xt
       y=yt
#ifdef OPT_EXMASK
       x_vert=xt_vert
       y_vert=yt_vert
#endif
    end select

    allocate(z(1:kzdim))
    allocate(z_bnd(1:nbnd,1:kzdim))
    z(1:kzdim)=zt(1:kzdim)
    z_bnd(1:nbnd,1:kzdim)=zt_bnd(1:nbnd,1:kzdim)
    cz_unit = 'm'
    cz_long_name = 'depth'
    if (nvcord(iitem) > 0) then
       z(1:kzdim)=lsig(1:kzdim,nvcord(iitem))
       z_bnd(1,1:kzdim)=lsigp(1:kzdim,nvcord(iitem))
       z_bnd(2,1:kzdim)=lsigp(2:kzdim+1,nvcord(iitem))
       cz_unit = 'kg/m3'
       write(cz_long_name,'("sigma density (reference depth: "i0" m)")') nint(zref(nvcord(iitem))*1.d-2)
    else
       select case (clas(iitem)(3:5))
       case ('LVM')
          z(1:kzdim)=zw(1:kzdim)
          z_bnd(1:nbnd,1:kzdim)=zw_bnd(1:nbnd,1:kzdim)
       case ('ICE')
          z(1:kzdim)=zi(1:kzdim)
          z_bnd(1:nbnd,1:kzdim)=zi_bnd(1:nbnd,1:kzdim)
          cz_unit = ''
          cz_long_name = 'ice thickness category'
       end select
    end if

    cyr =chead(48)( 1:4)
    cmon=chead(48)( 5:6)
    cday=chead(48)( 7:8)
    chr =chead(48)(10:11)
    cmin=chead(48)(12:13)
    csec=chead(48)(14:15)
    read(chead(48),'(i4.4,2i2.2,1x,3i2.2,1x)') idate
    call cyh2ss(time1, idate)
    read(chead(50),'(i6.6,5i2.2)') idate
    call cyh2ss(tout, idate)

    if (clas(iitem)(3:5) == 'SFC') then
       call nc_write(cf=trim(adjustl(cfitem(iitem)))//'.nc', fid=iitem, time=tout, &
            & nxg=nxg,nyg=nyg, nx=nx,ny=ny,nz=kzdim, irank=irank, jrank=jrank, &
            & shuffle=.true., deflate_level=1, &
            & buf4=dat, cvar=adjustl(citem(iitem)), v_unit=adjustl(cunit(iitem)), v_long_name=adjustl(ctitl(iitem)), &
            & fill_value=-1.e20, &
            & x=x, x_unit='degrees_east',  x_long_name='longitude', &
            & y=y, y_unit='degrees_north', y_long_name='latitude', &
#ifdef OPT_EXMASK
            & x_vert=x_vert, y_vert=y_vert, &
#endif
            & t_unit='hours since '//cyr//'-'//cmon//'-'//cday//' '//chr//':'//cmin//':'//csec, time1=time1)
    else
       call nc_write(cf=trim(adjustl(cfitem(iitem)))//'.nc', fid=iitem, time=tout, &
            & nxg=nxg,nyg=nyg, nx=nx,ny=ny,nz=kzdim, irank=irank, jrank=jrank, &
            & shuffle=.true., deflate_level=1, &
            & buf4=dat, cvar=adjustl(citem(iitem)), v_unit=adjustl(cunit(iitem)), v_long_name=adjustl(ctitl(iitem)), &
            & fill_value=-1.e20, &
            & x=x, x_unit='degrees_east',  x_long_name='longitude', &
            & y=y, y_unit='degrees_north', y_long_name='latitude', &
            & z=z, z_unit=trim(cz_unit),   z_long_name=cz_long_name, &
#ifdef OPT_EXMASK
            & x_vert=x_vert, y_vert=y_vert, &
#endif
            & z_bnd = z_bnd, &
            & t_unit='hours since '//cyr//'-'//cmon//'-'//cday//' '//chr//':'//cmin//':'//csec, time1=time1)
    end if
    deallocate(z)
    deallocate(z_bnd)

  contains

#ifdef OPT_EXMASK
    subroutine mod_lon(a,b)
      implicit none
      real(8), intent(inout) :: a, b
      do while(b - lon0 >= 360.d0)
         b=b-360.d0
      end do
      do while(b - lon0 < 0.d0)
         b=b+360.d0
      end do
      do while(a-b > 180.d0)
         a=a-360.d0
      end do
      do while(a-b < -180.d0)
         a=a+360.d0
      end do
      return
    end subroutine mod_lon
#endif
    
    subroutine def_coord

#ifdef OPT_EXMASK
      use zocdim, only : nxyg, igstr, jgstr, nxgdim, igend, jgend, nxygdm
      use bgs2d
#endif
      use zocgrd, only : &
#ifdef OPT_EXMASK
           & glont, glatt, &
#endif
           & dz0
      use ufile
      implicit none

#ifdef OPT_EXMASK
      real(8) :: glon(nxygdm), glat(nxygdm)
      real(8) :: drad2deg
      real(8) :: glon_ij, glon_ije, glon_ijn, glon_ijw, glon_ijs, glon_ijne, glon_ijse, glon_ijnw, glon_ijsw
      integer :: ij, ijn, ije, ijs, ijw, ijne, ijnw, ijse, ijsw
#endif
      integer :: i, j, k

      integer :: istat 
      integer :: nf_ncg
      character(len=ncf) :: c_ncg
      namelist /nm_ncg/ c_ncg
      data c_ncg /'not-specified'/

#ifdef OPT_EXMASK
      call gather_2d( glon, glont)
      call gather_2d( glat, glatt)
#endif
      if (irank + jrank == 0) then ! coordinate info is defined and written only root node
         call rewnml(ifpar, jfpar)
         read(ifpar, nm_ncg, iostat=istat)
         if (trim(c_ncg) /= 'not-specified') then
            call filopn(nf_ncg, c_ncg, 'READ')
            read(nf_ncg) xt
            read(nf_ncg) yt
            read(nf_ncg) xv
            read(nf_ncg) yv
            read(nf_ncg) zt
            read(nf_ncg) zw
            xx = xt
            yx = yt
            xy = xt
            yy = yt
         else
#ifdef OPT_EXMASK
            drad2deg = 45.d0 / atan(1.d0)
            glon = glon * drad2deg
            glat = glat * drad2deg
            lon0 = (3*glon(igstr+(jgstr-1)*nxgdim)-glon(igstr+1+(jgstr-1)*nxgdim))*0.5d0
            do j = jgstr, jgend
               glon(igstr-1+(j-1)*nxgdim) = glon(igend+(j-1)*nxgdim)
               glon(igend+1+(j-1)*nxgdim) = glon(igstr+(j-1)*nxgdim)
               glat(igstr-1+(j-1)*nxgdim) = glat(igend+(j-1)*nxgdim)
               glat(igend+1+(j-1)*nxgdim) = glat(igstr+(j-1)*nxgdim)
            end do
            do i = igstr-1, igend+1
               glon(i+(jgstr-2)*nxgdim) = 2*glon(i+(jgstr-1)*nxgdim)-glon(i+jgstr*nxgdim)
               glat(i+(jgstr-2)*nxgdim) = 2*glat(i+(jgstr-1)*nxgdim)-glat(i+jgstr*nxgdim)
#ifdef OPT_TRIPOLE
               glon(i+jgend*nxgdim) = glon(nxgdim-i+1+(jgend-1)*nxgdim)
               glat(i+jgend*nxgdim) = glat(nxgdim-i+1+(jgend-1)*nxgdim)
#else
               glon(i+jgend*nxgdim) = 2*glon(i+(jgend-1)*nxgdim)-glon(i+(jgend-2)*nxgdim)
               glat(i+jgend*nxgdim) = 2*glat(i+(jgend-1)*nxgdim)-glat(i+(jgend-2)*nxgdim)
#endif
            end do
            do j = 1, nyg
               do i = 1, nxg
                  ij   = i+igstr-1 + (j+jgstr-2)*nxgdim
                  ijn  = i+igstr-1 + (j+jgstr-1)*nxgdim
                  ijs  = i+igstr-1 + (j+jgstr-3)*nxgdim
                  ije  = i+igstr   + (j+jgstr-2)*nxgdim
                  ijw  = i+igstr-2 + (j+jgstr-2)*nxgdim
                  ijne = i+igstr   + (j+jgstr-1)*nxgdim
                  ijse = i+igstr   + (j+jgstr-3)*nxgdim
                  ijnw = i+igstr-2 + (j+jgstr-1)*nxgdim
                  ijsw = i+igstr-2 + (j+jgstr-3)*nxgdim
                  glon_ij = glon(ij)
                  glon_ijn = glon(ijn)
                  glon_ijs = glon(ijs)
                  glon_ije = glon(ije)
                  glon_ijw = glon(ijw)
                  glon_ijne = glon(ijne)
                  glon_ijse = glon(ijse)
                  glon_ijnw = glon(ijnw)
                  glon_ijsw = glon(ijsw)
                  call mod_lon(glon_ijw, glon_ij)
                  call mod_lon(glon_ijs, glon_ij)
                  call mod_lon(glon_ije, glon_ij)
                  call mod_lon(glon_ijn, glon_ij)
                  call mod_lon(glon_ijne, glon_ij)
                  call mod_lon(glon_ijnw, glon_ij)
                  call mod_lon(glon_ijse, glon_ij)
                  call mod_lon(glon_ijsw, glon_ij)
                  xt(i,j) = glon_ij
                  yt(i,j) = glat(ij)
                  xv(i,j) = (glon_ij  + glon_ije  + glon_ijn  + glon_ijne ) * 0.25d0
                  yv(i,j) = (glat(ij) + glat(ije) + glat(ijn) + glat(ijne)) * 0.25d0
                  yx(i,j) = glat(ij)
                  xy(i,j) = glon_ij
                  xx(i,j) = (glon_ij  + glon_ijw ) * 0.5d0
                  yy(i,j) = (glat(ij) + glat(ijs)) * 0.5d0
                  
                  xt_vert(1,i,j) = (glon_ij  + glon_ijw  + glon_ijs  + glon_ijsw ) * 0.25d0
                  xt_vert(2,i,j) = (glon_ij  + glon_ije  + glon_ijs  + glon_ijse ) * 0.25d0
                  xt_vert(3,i,j) = (glon_ij  + glon_ije  + glon_ijn  + glon_ijne ) * 0.25d0
                  xt_vert(4,i,j) = (glon_ij  + glon_ijw  + glon_ijn  + glon_ijnw ) * 0.25d0
                  yt_vert(1,i,j) = (glat(ij) + glat(ijw) + glat(ijs) + glat(ijsw)) * 0.25d0
                  yt_vert(2,i,j) = (glat(ij) + glat(ije) + glat(ijs) + glat(ijse)) * 0.25d0
                  yt_vert(3,i,j) = (glat(ij) + glat(ije) + glat(ijn) + glat(ijne)) * 0.25d0
                  yt_vert(4,i,j) = (glat(ij) + glat(ijw) + glat(ijn) + glat(ijnw)) * 0.25d0
                  
                  xv_vert(1,i,j) = glon_ij
                  xv_vert(2,i,j) = glon_ije
                  xv_vert(3,i,j) = glon_ijne
                  xv_vert(4,i,j) = glon_ijn
                  yv_vert(1,i,j) = glat(ij)
                  yv_vert(2,i,j) = glat(ije)
                  yv_vert(3,i,j) = glat(ijne)
                  yv_vert(4,i,j) = glat(ijn)

                  xx_vert(1,i,j) = (glon_ijw  + glon_ijsw ) * 0.5d0
                  xx_vert(2,i,j) = (glon_ij   + glon_ijs  ) * 0.5d0
                  xx_vert(3,i,j) = (glon_ij   + glon_ijn  ) * 0.5d0
                  xx_vert(4,i,j) = (glon_ijw  + glon_ijnw ) * 0.5d0
                  yx_vert(1,i,j) = (glat(ijw) + glat(ijsw)) * 0.5d0
                  yx_vert(2,i,j) = (glat(ij)  + glat(ijs) ) * 0.5d0
                  yx_vert(3,i,j) = (glat(ij)  + glat(ijn) ) * 0.5d0
                  yx_vert(4,i,j) = (glat(ijw) + glat(ijnw)) * 0.5d0

                  xy_vert(1,i,j) = (glon_ijs  + glon_ijsw ) * 0.5d0
                  xy_vert(2,i,j) = (glon_ijs  + glon_ijse ) * 0.5d0
                  xy_vert(3,i,j) = (glon_ij   + glon_ije  ) * 0.5d0
                  xy_vert(4,i,j) = (glon_ij   + glon_ijw  ) * 0.5d0
                  yy_vert(1,i,j) = (glat(ijs) + glat(ijsw)) * 0.5d0
                  yy_vert(2,i,j) = (glat(ijs) + glat(ijse)) * 0.5d0
                  yy_vert(3,i,j) = (glat(ij)  + glat(ije) ) * 0.5d0
                  yy_vert(4,i,j) = (glat(ij)  + glat(ijw) ) * 0.5d0

               end do
            end do
#else
            do i = 1, nxg
               xt(i) = dble(i) - 0.5d0
               xv(i) = dble(i)
               xx(i) = dble(i) - 1.d0
               xy(i) = dble(i) - 0.5d0
            end do
            do j = 1, nyg
               yt(j) = dble(j) - 0.5d0
               yv(j) = dble(j)
               yx(j) = dble(j) - 1.d0
               yy(j) = dble(j) - 0.5d0
            end do
#endif
            zt(1) = dz0(kstr) * 0.5d-2
            zw(1) = 0.d0
            zt_bnd(1,1) = zw(1)
            zt_bnd(2,1) = dz0(kstr) * 1.d-2
            zw_bnd(1,1) = - dz0(kstr) * 0.5d-2
            zw_bnd(2,1) = zt(1)
            do k = 2, nz
               zt(k) = zt(k-1) + (dz0(k+kstr-1) + dz0(k+kstr-2)) * 0.5d-2
               zw(k) = zw(k-1) + dz0(k+kstr-2) * 1.d-2
               zt_bnd(1,k) = zt_bnd(2,k-1)
               zt_bnd(2,k) = zt_bnd(1,k) + dz0(k+kstr-1) * 1.d-2
               zw_bnd(1,k) = zw_bnd(2,k-1)
               zw_bnd(2,k) = zt(k)
            end do
         end if
         do k = 1, nic
            zi(k) = dble(k)
            zi_bnd(1,k) = dble(k) - 0.5d0
            zi_bnd(2,k) = dble(k) + 0.5d0
         end do
      end if

    end subroutine def_coord

  end subroutine write_ncfile

end module qckot

