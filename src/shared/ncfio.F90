#ifdef OPT_IO_NCF
module ncfio
  use mpi
  use netcdf
  use zocdim, only : mpi_comm_ogcm
  use zocfil, only : nfomax
  use ufile
  implicit none

  logical, save :: ofirst(nfomax) = .true.
  real(8), save :: time0(nfomax) = 0.d0
  integer :: ifpar, jfpar
  
  private
  public :: nc_write, nc_filopn, nc_read_chead, nc_read_sfc, nc_read_bdy
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

  subroutine nc_read_chead( &
    chead, ncid, disp, icread)
    
    character(16),                 intent(out)   :: chead(1:64)
    integer,                       intent(out)   :: icread
    integer(kind=mpi_offset_kind), intent(inout) :: disp
    integer,                       intent(in)    :: ncid

    integer :: varid, dimid
    integer(4) :: xtime_len
    real(8), allocatable :: time_values(:)
    character(256) :: time_units

    call check(nf90_inq_dimid(ncid, "time", dimid))
    call check(nf90_inquire_dimension(ncid, dimid, len = xtime_len))

    if (disp >= xtime_len) then
       icread = 0
       return
    else
       icread = 1024
       disp = disp + 1
    end if
    
    call check(nf90_inq_varid(ncid, "time", varid))
    allocate(time_values(xtime_len))

    call check(nf90_get_var(ncid, varid, time_values))
    call check(nf90_get_att(ncid, varid, "units", time_units))

    chead(50) = conv_to_chead50(time_values(disp),time_units)
    deallocate(time_values)
    
  end subroutine nc_read_chead

  subroutine nc_read_sfc( &
    data, ncid, disp)
    use zocdim, only : nx, ny, irank, jrank
    real(8),                       intent(out) :: data(nx, ny)
    integer,                       intent(in)  :: ncid
    integer(kind=mpi_offset_kind), intent(in)  :: disp

    real(8) :: data3d(nx, ny, 1)
    real(8) :: fill_val
    integer :: varid, xtype, ndims, dimids(10), natts, nvars
    character(len=NF90_MAX_NAME) :: tmp_name
    integer :: starts(1:3), counts(1:3)
    integer :: i
    integer(4) :: disp4
    logical :: is_coord

    disp4 = disp
    starts(1:3)= (/ irank*nx+1, jrank*ny+1, disp4/)
    counts(1:3)= (/nx, ny, 1 /)

    call check(nf90_inquire(ncid, nVariables = nvars))
    do varid = 1, nvars
       call check(nf90_inquire_variable(ncid, varid, tmp_name, xtype, ndims, dimids, natts))
       if (ndims == 3) then
          is_coord = .false.
          do i = 1, ndims
             call check(nf90_inquire_dimension(ncid, dimids(i), tmp_name))
             if ( index(tmp_name, 'vertices') > 0 .or. &
                  index(tmp_name, 'bnds') > 0 .or. &
                  index(tmp_name, 'vertex') > 0 ) then
                is_coord = .true.
                exit
             end if
          end do
          if (is_coord) then
             cycle
          else
             exit
          end if
       end if
    end do


    call check(nf90_var_par_access(ncid, varid, NF90_COLLECTIVE))
    call check(nf90_get_var(ncid, varid, data3d, start = starts, count = counts))
    data(:,:) = data3d(:,:,1)

    call check(nf90_get_att(ncid, varid, "_FillValue", fill_val))
    where(data == fill_val) data = 0.0d0
    
  end subroutine nc_read_sfc

  subroutine nc_read_bdy( &
    data, ncid, disp)
    use zocdim, only : nx, ny, nz, irank, jrank
    real(8),                       intent(out) :: data(nx, ny, nz)
    integer,                       intent(in)  :: ncid
    integer(kind=mpi_offset_kind), intent(in)  :: disp

    real(8) :: data4d(nx, ny, nz, 1)
    real(8) :: fill_val
    integer :: varid, xtype, ndims, dimids(10), natts, nvars
    character(80) :: tmp_name
    integer :: starts(1:4), counts(1:4)
    integer(4) :: disp4

    disp4 = disp
    
    starts(1:4)= (/ irank*nx+1, jrank*ny+1, 1, disp4/)
    counts(1:4)= (/nx, ny, nz, 1/)

    call check(nf90_inquire(ncid, nVariables = nvars))
    do varid = 1, nvars
       call check(nf90_inquire_variable(ncid, varid, tmp_name, xtype, ndims, dimids, natts))
       if (ndims == 4) then
          exit
       end if
    end do

    call check(nf90_var_par_access(ncid, varid, NF90_COLLECTIVE))
    call check(nf90_get_var(ncid, varid, data4d, start = starts, count = counts))
    data(:,:,:) = data4d(:,:,:,1)

    call check(nf90_get_att(ncid, varid, "_FillValue", fill_val))
    where(data == fill_val) data = 0.0d0
 
  end subroutine nc_read_bdy

  subroutine nc_filopn( &
    ncid, cf, cact)
    character(*), intent(in) :: cf
    character(*), intent(in) :: cact
    integer, intent(out) :: ncid

    if (cact(1:1) == 'R' .or. cact(1:1) == 'r') then
       call check(nf90_open_par(cf, nf90_nowrite,  comm = mpi_comm_ogcm, info = MPI_INFO_NULL, ncid = ncid))
    else
       call check(nf90_open_par(cf, nf90_write, comm = mpi_comm_ogcm, info = MPI_INFO_NULL, ncid = ncid))
    end if
    
  end subroutine nc_filopn

  function conv_to_chead50(time_val, time_units)
    use ucaln
    implicit none
    character(16) :: conv_to_chead50  ! chead(50), YYYYYYMMDDHHMMSS
    real(8), intent(in) :: time_val      ! NetCDFから読み込んだ数値
    character(*), intent(in) :: time_units ! "e.g., seconds since 1970-01-01..."

    integer :: pos0, pos
    character(16) :: unit_part
    real(8) :: total_seconds, tt1, tt2, dtime_val
    integer :: idate(6), idate1(6), idate2(6)

    pos0 = index(time_units, "since")
    if (pos0 == 0) return
    unit_part = adjustl(time_units(1:pos0-1))

    ! get base time
    idate(1:6) = (/0,1,1,0,0,0/)
    pos = index(time_units, '-') 
    if (pos > 0) then
       read(time_units(pos0+5:pos-1), '(I)') idate(1)
       read(time_units(pos+1:pos+2), '(I2)') idate(2)
       read(time_units(pos+4:pos+5), '(I2)') idate(3)
    end if
    pos = index(time_units, ':') 
    if (pos > 0) then
       read(time_units(pos-2:pos-1), '(I2)') idate(4)
       read(time_units(pos+1:pos+2), '(I2)') idate(5)
       read(time_units(pos+4:pos+5), '(I2)') idate(6)
    end if

    if (index(unit_part, 'year') > 0 .or. index(unit_part, 'month') > 0) then
       idate1 = idate
       idate2 = idate
       if (index(unit_part, 'year') > 0) then
          idate1(1) = idate1(1) + int(time_val)
          idate2(1) = idate2(1) + int(time_val) + 1
       else
          idate1(2) = idate1(2) + int(time_val)
          idate2(2) = idate2(2) + int(time_val) + 1
       end if
       dtime_val = time_val - dble(int(time_val))
       call cyh2ss(tt1,idate1)
       call cyh2ss(tt2,idate2)
       total_seconds = tt1 * (1.d0 - dtime_val) + tt2 * dtime_val
    else
       call cyh2ss(total_seconds, idate)
       if (index(unit_part, 'day') > 0) then
          total_seconds = total_seconds + time_val * 8.64d4
       else if (index(unit_part, 'hour') > 0) then
          total_seconds = total_seconds + time_val * 3.6d3
       else if (index(unit_part, 'minute') > 0) then
          total_seconds = total_seconds + time_val * 6.d1
       else if (index(unit_part, 'second') > 0) then
          total_seconds = total_seconds + time_val
       end if
    end if
    
    call css2yh(idate,total_seconds)
    write(conv_to_chead50, '(I6.6, I2.2, I2.2, I2.2, I2.2, I2.2)') idate

  end function conv_to_chead50
  
end module ncfio
#else
subroutine ncf_io
  return
end subroutine ncf_io
#endif
