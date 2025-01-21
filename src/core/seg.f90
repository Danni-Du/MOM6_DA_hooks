program seg
  use gsw_mod_toolbox
  implicit none

  ! Declare variables
  character(len=255), dimension(6) :: flux_filenames, flux_varnames
  character(len=255), dimension(2) :: tracer_filenames, tracer_varnames
  character(len=255), dimension(2) :: current_filenames, current_varnames
  character(len=255) :: grid_filename
  integer :: i,  ix, jy, zz, tt, ix_left, ix_right, jy_south, jy_north
  real(8) :: SA, PT, CT, PRHO
  real(8) :: hflso, hfsso, rlntds, rsntds, taux_left, taux_right, tauy_south, tauy_north, tauamp
  real(8) ::  thetao, so, uo_left, uo_right, vo_south, vo_north, div, thetao_top, thetao_bottom, so_top,so_bottom,so_zgrad,thetao_zgrad, PRHO_top, PRHO_bottom, PRHO_zgrad
  real(8) :: dyCu_left, dyCu_right, dxCv_south, dxCv_north, areacello
  real(8), dimension(75) :: z_l
  real(8), dimension(76) :: z_i
  integer :: zl_index_mld, zl_index10m, zl_index_3mld, right_index
  real(8), dimension(:), allocatable :: PRHO_profile
  real(8) :: value_for_control_PRHO, value_for_control_depth, reference_depth, PRHO_change, PRHO_mld, PRHO_10m,value_for_control_oceanzvars, dummy_var 
  real(8) :: mld_depth
  real(8), dimension(:), allocatable :: zl_to_sigma, zi_to_sigma
  real(8), dimension(15) :: target_sigmas, thetao_zgrad_sigma, so_zgrad_sigma, PRHO_zgrad_sigma, div_sigma
  real(8), dimension(:), allocatable :: thetao_zgrad_profile, so_zgrad_profile, div_profile, PRHO_zgrad_profile
  character(len=255)  :: danni_ANN_name
  real(8), dimension(16,51)  :: l1_weight
  real(8), dimension(16,16)  :: l2_weight, l3_weight
  real(8), dimension(16) :: l1_bias, l2_bias, l3_bias, l1_output, l2_output, l3_output
  real(8) :: ReLU_zero
  real(8), dimension(51) :: ANN_input
  
  ReLU_zero = 0
  target_sigmas = (/0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9/)
  !print *, target_sigmas

  PRHO_change = 0.03
  reference_depth = 10
  value_for_control_depth = 1E6
  value_for_control_PRHO = 999
  value_for_control_oceanzvars = 1E5
    
  grid_filename = '/scratch/cimes/dd7201/DA_increments/DiffDA/M_prod_inc/ocean_z.static.nc'
  
  ! given horizontal location and time
  ix = 201
  jy = 151
  tt = 5
  ! for taux, tauy, uo and vo
  ix_left = ix - 1
  ix_right = ix
  jy_south = jy -1 
  jy_north = jy
 
  ! To read data (don't need this part in MOM6 online implementation)
  ! Initialize file names
  flux_filenames(1) = '/scratch/cimes/dd7201/DA_increments/DiffDA/M_prod_inc/ocean_daily.20060101-20061231.hflso.nc'
  flux_filenames(2) = '/scratch/cimes/dd7201/DA_increments/DiffDA/M_prod_inc/ocean_daily.20060101-20061231.hfsso.nc'
  flux_filenames(3) = '/scratch/cimes/dd7201/DA_increments/DiffDA/M_prod_inc/ocean_daily.20060101-20061231.rlntds.nc'
  flux_filenames(4) = '/scratch/cimes/dd7201/DA_increments/DiffDA/M_prod_inc/ocean_daily.20060101-20061231.rsntds.nc'
  flux_filenames(5) = '/scratch/cimes/dd7201/DA_increments/DiffDA/M_prod_inc/ocean_daily.20060101-20061231.taux.nc'
  flux_filenames(6) = '/scratch/cimes/dd7201/DA_increments/DiffDA/M_prod_inc/ocean_daily.20060101-20061231.tauy.nc'
  tracer_filenames(1) = '/scratch/cimes/dd7201/DA_increments/DiffDA/M_prod_inc/ocean_z_daily.20060101-20061231.thetao.nc'
  tracer_filenames(2) = '/scratch/cimes/dd7201/DA_increments/DiffDA/M_prod_inc/ocean_z_daily.20060101-20061231.so.nc'
  current_filenames(1) = '/scratch/cimes/dd7201/DA_increments/DiffDA/M_prod_inc/ocean_z_daily.20060101-20061231.uo.nc'
  current_filenames(2) = '/scratch/cimes/dd7201/DA_increments/DiffDA/M_prod_inc/ocean_z_daily.20060101-20061231.vo.nc'
  ! Initialize varible names
  flux_varnames(1) = 'hflso'
  flux_varnames(2) = 'hfsso'
  flux_varnames(3) = 'rlntds'
  flux_varnames(4) = 'rsntds'
  flux_varnames(5) = 'taux'
  flux_varnames(6) = 'tauy'
  tracer_varnames(1) = 'thetao'
  tracer_varnames(2) = 'so'
  current_varnames(1) = 'uo'
  current_varnames(2) = 'vo'
  ! Loop over/ choose each file and call the subroutine to read it
  call read_netcdf_file_3d(flux_filenames(1), flux_varnames(1), ix, jy, tt, hflso)
  call read_netcdf_file_3d(flux_filenames(2), flux_varnames(2), ix, jy, tt, hfsso)
  call read_netcdf_file_3d(flux_filenames(3), flux_varnames(3), ix, jy, tt, rlntds)
  call read_netcdf_file_3d(flux_filenames(4), flux_varnames(4), ix, jy, tt, rsntds)
  call read_netcdf_file_3d(flux_filenames(5), flux_varnames(5), ix_left, jy, tt, taux_left)
  call read_netcdf_file_3d(flux_filenames(5), flux_varnames(5), ix_right, jy, tt, taux_right)
  call read_netcdf_file_3d(flux_filenames(6), flux_varnames(6), ix, jy_south, tt, tauy_south)
  call read_netcdf_file_3d(flux_filenames(6), flux_varnames(6), ix, jy_north, tt, tauy_north)
  tauamp = sqrt(((taux_left+taux_right)/2)**2+((tauy_south+tauy_north)/2)**2)
  ! grid info for computing divergence (depth independent)
  call read_grid_file(grid_filename,'dxCv',ix,jy_south,dxCv_south)
  call read_grid_file(grid_filename,'dxCv',ix,jy_north,dxCv_north)
  call read_grid_file(grid_filename,'dyCu',ix_left,jy,dyCu_left)
  call read_grid_file(grid_filename,'dyCu',ix_right,jy,dyCu_right)
  call read_grid_file(grid_filename, 'areacello', ix, jy, areacello)
  call read_netcdf_file_zlzi(tracer_filenames(1), 'z_l', 'z_i',z_l,z_i)
  
  ! To compute the potential density profile using GSW_Fortran 
  ! Absolute Salinity in g/kg
  ! Conservative Temperature in degrees Celsius (converted from potential temprature)
  do zz  = 1, 75
    call read_netcdf_file_4d(tracer_filenames(1), tracer_varnames(1), ix, jy, zz,  tt, PT)
    call read_netcdf_file_4d(tracer_filenames(2), tracer_varnames(2), ix, jy, zz, tt, SA)
    CT = gsw_ct_from_pt(SA,PT)
    PRHO = gsw_sigma0(SA,CT)
    if (zz == 1) then
      allocate(PRHO_profile(1))
      PRHO_profile(1) = PRHO
    else
      call append_value(PRHO_profile, PRHO)
    end if
  end do

  ! 3MLD 
  ! find the first index below 10m
  call find_right_index(z_l, reference_depth, value_for_control_depth, zl_index10m)
  ! the 10m potential density
  call interpolate(z_l(zl_index10m-1),z_l(zl_index10m),PRHO_profile(zl_index10m-1), PRHO_profile(zl_index10m),reference_depth,PRHO_10m)
  ! the MLD potential density
  PRHO_mld = PRHO_10m + PRHO_change
  ! the first z_l index below MLD
  call find_right_index(PRHO_profile, PRHO_mld,value_for_control_PRHO, zl_index_mld)
  ! the MLD depth
  call interpolate(PRHO_profile(zl_index_mld-1),PRHO_profile(zl_index_mld),z_l(zl_index_mld-1),z_l(zl_index_mld),PRHO_mld,mld_depth)
  ! the MLD must be below 10m
  if (mld_depth < 10) then
    mld_depth = 10
  end if
  ! the first z_l index below 3MLD
  call find_right_index(z_l, 3*mld_depth,value_for_control_depth, zl_index_3mld)
  
  ! to confirm that  at (zl_index_3mld + 1) for zgrad and at (zl_index_3mld) for div, the input values are non nan
  call read_netcdf_file_4d(tracer_filenames(1), tracer_varnames(1), ix, jy, zl_index_3mld+1,  tt, thetao)
  call read_netcdf_file_4d(tracer_filenames(2), tracer_varnames(2), ix, jy, zl_index_3mld+1, tt, so)
  call read_netcdf_file_4d(current_filenames(1), current_varnames(1), ix_left, jy, zl_index_3mld, tt, uo_left)
  call read_netcdf_file_4d(current_filenames(1), current_varnames(1), ix_right, jy, zl_index_3mld, tt, uo_right)
  call read_netcdf_file_4d(current_filenames(2), current_varnames(2), ix, jy_south, zl_index_3mld, tt, vo_south)
  call read_netcdf_file_4d(current_filenames(2), current_varnames(2), ix, jy_north, zl_index_3mld, tt, vo_north)
  dummy_var = (thetao + so + uo_left + uo_right + vo_south + vo_north)/6 
  ! if not nan, then get the vertical profiles
  if (abs(dummy_var) < value_for_control_oceanzvars) then
    zi_to_sigma = z_i(2:zl_index_3mld + 1)/mld_depth
    do zz = 1, zl_index_3mld
      call read_netcdf_file_4d(tracer_filenames(1), tracer_varnames(1), ix, jy, zz,  tt, thetao)
      call read_netcdf_file_4d(tracer_filenames(2), tracer_varnames(2), ix, jy, zz, tt, so)
      thetao_top = thetao
      so_top = so
      CT = gsw_ct_from_pt(so_top,thetao_top)
      PRHO_top = gsw_sigma0(so_top,CT)
      
      call read_netcdf_file_4d(tracer_filenames(1), tracer_varnames(1), ix, jy, zz+1,  tt, thetao)
      call read_netcdf_file_4d(tracer_filenames(2), tracer_varnames(2), ix, jy, zz+1, tt, so)
      thetao_bottom = thetao
      so_bottom = so
      CT = gsw_ct_from_pt(so_bottom,thetao_bottom)
      PRHO_bottom = gsw_sigma0(so_bottom,CT)


      thetao_zgrad = (thetao_top - thetao_bottom)/(z_l(zz+1) - z_l(zz))
      so_zgrad = (so_top - so_bottom)/(z_l(zz+1) - z_l(zz))
      PRHO_zgrad = (PRHO_top - PRHO_bottom)/(z_l(zz+1) - z_l(zz))
      ! append them to 1D array
      if (zz == 1) then
        allocate(thetao_zgrad_profile(1))
        thetao_zgrad_profile(1) = thetao_zgrad
        allocate(so_zgrad_profile(1))
        so_zgrad_profile(1) = so_zgrad
        allocate(PRHO_zgrad_profile(1))
        PRHO_zgrad_profile(1) = PRHO_zgrad
      else
        call append_value(thetao_zgrad_profile, thetao_zgrad)
        call append_value(so_zgrad_profile, so_zgrad)
        call append_value(PRHO_zgrad_profile, PRHO_zgrad)
      end if
    end do
    zl_to_sigma = z_l(1:zl_index_3mld)/mld_depth
    do zz = 1, zl_index_3mld
      call read_netcdf_file_4d(current_filenames(1), current_varnames(1), ix_left, jy, zz, tt, uo_left)
      call read_netcdf_file_4d(current_filenames(1), current_varnames(1), ix_right, jy, zz, tt, uo_right)
      call read_netcdf_file_4d(current_filenames(2), current_varnames(2), ix, jy_south, zz, tt, vo_south)
      call read_netcdf_file_4d(current_filenames(2), current_varnames(2), ix, jy_north, zz, tt, vo_north) 
      call compute_current_divergence(uo_left*dyCu_left, uo_right*dyCu_right, vo_south*dxCv_south, vo_north*dxCv_north, areacello, div)
      if (zz == 1) then
        allocate(div_profile(1))
        div_profile(1) = div
      else
        call append_value(div_profile, div)
      end if
     
    end do 
  
  end if
  
  ! interpolate values to target_sigmas
  do i = 1, 15
    call find_right_index(zi_to_sigma, target_sigmas(i),value_for_control_depth, right_index)
    if (right_index == 1) then
      thetao_zgrad_sigma(i) = thetao_zgrad_profile(1)
      so_zgrad_sigma(i) = so_zgrad_profile(1)
      PRHO_zgrad_sigma(i) = PRHO_zgrad_profile(1)
    else
      call interpolate(zi_to_sigma(right_index-1),zi_to_sigma(right_index),thetao_zgrad_profile(right_index-1),thetao_zgrad_profile(right_index),target_sigmas(i),thetao_zgrad_sigma(i))
      call interpolate(zi_to_sigma(right_index-1),zi_to_sigma(right_index),so_zgrad_profile(right_index-1),so_zgrad_profile(right_index),target_sigmas(i),so_zgrad_sigma(i))
      call interpolate(zi_to_sigma(right_index-1),zi_to_sigma(right_index),PRHO_zgrad_profile(right_index-1),PRHO_zgrad_profile(right_index),target_sigmas(i),PRHO_zgrad_sigma(i))
    end if

    call find_right_index(zl_to_sigma, target_sigmas(i),value_for_control_depth, right_index)
    if (right_index == 1) then
      div_sigma(i) = div_profile(1)
    else
      call interpolate(zl_to_sigma(right_index-1),zl_to_sigma(right_index),div_profile(right_index-1),div_profile(right_index),target_sigmas(i),div_sigma(i))
    end if
  end do
  
  !print *, PRHO_zgrad_sigma
  !print *, zl_to_sigma
  !print *, hflso
  !print *, hfsso
  !print *, rlntds
  !print *, rsntds
  !print *, tauamp
  !print *, thetao
  !print *, so
  !print *, uo_left
  !print *, uo_right
  !print *, vo_south
  !print *, vo_north
  !print *, dyCu_left
  !print *, dyCu_right
  !print *, dxCv_south
  !print *, dxCv_north
  !print *, areacello
  !print *, div
  !print *, z_l
  !print *, z_i

  ! load the NN weights and biases
  ! subroutine(input,DA tendency)
  danni_ANN_name = '/scratch/cimes/dd7201/pp_DA_increments/danni_data_20250120/argo_only_clean/networks/danni_ANN_weights.nc'
  call read_ANN_file(danni_ANN_name, l1_weight, l2_weight, l3_weight, l1_bias, l2_bias, l3_bias)
   
  ANN_input(1:15) = thetao_zgrad_sigma*100
  ANN_input(16:30) = PRHO_zgrad_sigma*100
  ANN_input(31:45) = div_sigma*1E7
  ANN_input(46) = mld_depth*0.1
  ANN_input(47) = tauamp*100
  ANN_input(48) = hflso*0.1
  ANN_input(49) = hfsso*0.1
  ANN_input(50) = rlntds*0.1
  ANN_input(51) = rsntds*0.1

  l1_output = max(ReLU_zero, matmul(l1_weight, ANN_input) + l1_bias)
  l2_output = max(ReLU_zero, matmul(l2_weight, l1_output) + l2_bias)
  l3_output = matmul(l3_weight, l2_output) + l3_bias
  !print *, l3_output
  ! l3_output is the predicted flux
  print *, (l3_output(1:15)-l3_output(2:16))/(0.2*mld_depth)
contains
  Subroutine read_ANN_file(filename, l1_weight, l2_weight, l3_weight, l1_bias, l2_bias, l3_bias)
    use netcdf
    implicit none

    integer :: ncid, varid, retval
    character(len = 255) :: varname
    character(len=*), intent(in) :: filename
    real(8), dimension(51, 16)  :: l1_weight_temp
    real(8), dimension(16, 51), intent(out) :: l1_weight
    real(8), dimension(16,16), intent(out) :: l2_weight, l3_weight
    real(8), dimension(16), intent(out) :: l1_bias, l2_bias, l3_bias

    ! Open the NetCDF file
    retval = nf90_open(filename, nf90_nowrite, ncid)
    if (retval /= nf90_noerr) then
      print *, 'Error: Unable to open file'
      stop
    endif

    ! Get the variable ID
    varname = 'l1_weight'
    retval = nf90_inq_varid(ncid, varname, varid)
    if (retval == nf90_noerr) then
      ! Read the dimension values
      retval = nf90_get_var(ncid, varid, l1_weight_temp)
      l1_weight = transpose(l1_weight_temp)
      if (retval /= nf90_noerr) then
        print *, 'Error: Unable to get l1 weight values'
        stop
      endif
    else
      print *, 'Error: l1 weight variable not found'
    endif

    varname = 'l2_weight'
    retval = nf90_inq_varid(ncid, varname, varid)
    if (retval == nf90_noerr) then
      ! Read the dimension values
      retval = nf90_get_var(ncid, varid, l2_weight)
      l2_weight = transpose(l2_weight)
      if (retval /= nf90_noerr) then
        print *, 'Error: Unable to get l2 weight values'
        stop
      endif
    else
      print *, 'Error: l2 weight variable not found'
    endif

    varname = 'l3_weight'
    retval = nf90_inq_varid(ncid, varname, varid)
    if (retval == nf90_noerr) then
      ! Read the dimension values
      retval = nf90_get_var(ncid, varid, l3_weight)
      l3_weight = transpose(l3_weight)
      if (retval /= nf90_noerr) then
        print *, 'Error: Unable to get l3 weight values'
        stop
      endif
    else
      print *, 'Error: l3 weight variable not found'
    endif

    varname = 'l1_bias'
    retval = nf90_inq_varid(ncid, varname, varid)
    if (retval == nf90_noerr) then
      ! Read the dimension values
      retval = nf90_get_var(ncid, varid, l1_bias)
      if (retval /= nf90_noerr) then
        print *, 'Error: Unable to get l1 bias values'
        stop
      endif
    else
      print *, 'Error: l1 bias variable not found'
    endif

    varname = 'l2_bias'
    retval = nf90_inq_varid(ncid, varname, varid)
    if (retval == nf90_noerr) then
      ! Read the dimension values
      retval = nf90_get_var(ncid, varid, l2_bias)
      if (retval /= nf90_noerr) then
        print *, 'Error: Unable to get l2 bias values'
        stop
      endif
    else
      print *, 'Error: l2 bias variable not found'
    endif

    varname = 'l3_bias'
    retval = nf90_inq_varid(ncid, varname, varid)
    if (retval == nf90_noerr) then
      ! Read the dimension values
      retval = nf90_get_var(ncid, varid, l3_bias)
      if (retval /= nf90_noerr) then
        print *, 'Error: Unable to get l3 bias values'
        stop
      endif
    else
      print *, 'Error: l3 bias variable not found'
    endif

    ! Close the NetCDF file
    retval = nf90_close(ncid)
    if (retval /= nf90_noerr) then
      print *, 'Error: Unable to close file'
      stop
    endif
  end subroutine read_ANN_file

  ! 1D linear interpolation; it is guaranteed that x1 <= thisx < x2
  subroutine interpolate(x1,x2,y1,y2,thisx,thisy)
    implicit none
    real(8), intent(in) :: x1, x2, y1, y2, thisx
    real(8), intent(out) :: thisy
    thisy = (thisx-x1)/(x2-x1)*(y2-y1) + y1
  end subroutine interpolate

  ! to get an 1D array with variable size, appending values
  subroutine append_value(old_array, new_value)
    implicit none
    real(8), dimension(:), allocatable, intent(inout) :: old_array
    real(8), intent(in) :: new_value
    integer :: old_size, new_size
    real(8), dimension(:), allocatable :: temp_array
    old_size = size(old_array)
    new_size = old_size + 1

    ! Temporarily allocate a new array to hold the combined values
    allocate(temp_array(new_size))

    ! Copy old values to the new array
    temp_array(1:old_size) = old_array

    ! Append the new value to the new array
    temp_array(new_size) = new_value

    ! Deallocate the old array and reallocate it with the new size
    deallocate(old_array)
    allocate(old_array(new_size))

    ! Copy the combined values back to the original array
    old_array = temp_array

    ! Deallocate the temporary array
    deallocate(temp_array)
  end subroutine append_value

  ! Subroutine to find the index. 1D array (index) is greater than the given value
  Subroutine find_right_index(array1d_for_indexing, value_for_indexing,value_for_control, right_index)
    implicit none
    real(8), intent(in) :: array1d_for_indexing(:)
    integer :: array1d_i
    integer, intent(out) :: right_index
    real(8), intent(in) :: value_for_indexing, value_for_control
    right_index = 0 ! if there's so such right index, it will be 0 
    do array1d_i = 1, size(array1d_for_indexing)
      if (abs(array1d_for_indexing(array1d_i)) > value_for_control) then
        return
      else if (array1d_for_indexing(array1d_i) > value_for_indexing) then
        right_index = array1d_i
        return
      end if
    end do
    


  end subroutine find_right_index
  ! Subroutine to compute divergence
  Subroutine compute_current_divergence(uy_left, uy_right, vx_south, vy_north, area, div)
    implicit none
    real(8), intent(in) :: uy_left, uy_right, vx_south, vy_north, area
    real(8), intent(out) :: div

    div = (uy_right - uy_left + vy_north - vx_south) / area
  
  end subroutine compute_current_divergence

  ! Subroutine to read a grid NetCDF file
  Subroutine read_grid_file(filename, varname, xdimindex, ydimindex, varvalue)
    use netcdf
    implicit none

    integer :: ncid, varid, retval
    integer, intent(in) :: xdimindex, ydimindex
    character(len=*), intent(in) :: filename, varname
    real(8), intent(out) :: varvalue

    ! Open the NetCDF file
    retval = nf90_open(filename, nf90_nowrite, ncid)
    if (retval /= nf90_noerr) then
      print *, 'Error: Unable to open file'
      stop
    endif
    
    ! Get the variable ID
    retval = nf90_inq_varid(ncid, varname, varid)
    if (retval == nf90_noerr) then
      ! Read the values
      retval = nf90_get_var(ncid, varid, varvalue, start = (/xdimindex,ydimindex/))
      if (retval /= nf90_noerr) then
        print *, 'Error: Unable to get the grid info'
        stop
      endif
    else
      print *, 'Error: this grid info not found'
    endif
    
    ! Close the NetCDF file
    retval = nf90_close(ncid)
    if (retval /= nf90_noerr) then
      print *, 'Error: Unable to close file'
      stop
    endif

  end subroutine read_grid_file
  

  ! Subroutine to read a single (time, lev, lat, lon) NetCDF file
  subroutine read_netcdf_file_zlzi(filename, dimname1, dimname2,dimvalues1,dimvalues2)
    use netcdf
    implicit none

    ! Declare variables
    integer :: ncid, varid1, varid2, retval
    real(8), dimension(75), intent(out) :: dimvalues1
    real(8), dimension(76), intent(out) :: dimvalues2
    character(len=*), intent(in) :: filename, dimname1, dimname2


    ! Open the NetCDF file
    retval = nf90_open(filename, nf90_nowrite, ncid)
    if (retval /= nf90_noerr) then
      print *, 'Error: Unable to open file'
      stop
    endif

    ! Get the variable ID for varname
    retval = nf90_inq_varid(ncid, dimname1, varid1)
    if (retval /= nf90_noerr) then
      print *, 'Error: Unable to get variable ID'
      stop
    endif
    retval = nf90_inq_varid(ncid, dimname2, varid2)
    if (retval /= nf90_noerr) then
      print *, 'Error: Unable to get variable ID'
      stop
    endif


    ! Read the data
    retval = nf90_get_var(ncid, varid1, dimvalues1)

    if (retval /= nf90_noerr) then
      print *, 'Error: Unable to read data'
      stop
    endif

    retval = nf90_get_var(ncid, varid2, dimvalues2)

    if (retval /= nf90_noerr) then
      print *, 'Error: Unable to read data'
      stop
    endif

    ! Close the NetCDF file
    retval = nf90_close(ncid)
    if (retval /= nf90_noerr) then
      print *, 'Error: Unable to close file'
      stop
    endif

  end subroutine read_netcdf_file_zlzi

  ! Subroutine to read a single (time, lev, lat, lon) NetCDF file
  subroutine read_netcdf_file_4d(filename, varname, xdimindex, ydimindex, zdimindex, timedimindex, varvalue)
    use netcdf 
    implicit none

    ! Declare variables
    integer :: ncid, varid, retval
    integer, intent(in) :: xdimindex, ydimindex, timedimindex, zdimindex
    real(8), intent(out) :: varvalue
    character(len=*), intent(in) :: filename, varname

 
    ! Open the NetCDF file
    retval = nf90_open(filename, nf90_nowrite, ncid)
    if (retval /= nf90_noerr) then
      print *, 'Error: Unable to open file'
      stop
    endif

    ! Get the variable ID for varname
    retval = nf90_inq_varid(ncid, varname, varid)
    if (retval /= nf90_noerr) then
      print *, 'Error: Unable to get variable ID'
      stop
    endif


    ! Read the data
    retval = nf90_get_var(ncid, varid, varvalue, start = (/xdimindex, ydimindex,zdimindex, timedimindex/))

    if (retval /= nf90_noerr) then
      print *, 'Error: Unable to read data'
      stop
    endif

    ! Close the NetCDF file
    retval = nf90_close(ncid)
    if (retval /= nf90_noerr) then
      print *, 'Error: Unable to close file'
      stop
    endif

  end subroutine read_netcdf_file_4d

  ! Subroutine to read a single (time, lat, lon) NetCDF file
  subroutine read_netcdf_file_3d(filename, varname,xdimindex, ydimindex, timedimindex, varvalue)
    use netcdf 
    implicit none

    ! Declare variables
    integer :: ncid, varid, retval
    integer, intent(in) :: xdimindex, ydimindex, timedimindex
    real(8) :: varvalue
    character(len=*), intent(in) :: filename, varname

 
    ! Open the NetCDF file
    retval = nf90_open(filename, nf90_nowrite, ncid)
    if (retval /= nf90_noerr) then
      print *, 'Error: Unable to open file'
      stop
    endif

    ! Get the variable ID for varname
    retval = nf90_inq_varid(ncid, varname, varid)
    if (retval /= nf90_noerr) then
      print *, 'Error: Unable to get variable ID'
      stop
    endif


    ! Read the data
    retval = nf90_get_var(ncid, varid, varvalue, start = (/xdimindex, ydimindex, timedimindex/))
    if (retval /= nf90_noerr) then
      print *, 'Error: Unable to read data'
      stop
    endif

    ! Close the NetCDF file
    retval = nf90_close(ncid)
    if (retval /= nf90_noerr) then
      print *, 'Error: Unable to close file'
      stop
    endif

  end subroutine read_netcdf_file_3d
end program seg
