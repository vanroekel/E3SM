module agi_read

  !-----------------------------------------------------------------!
  ! Module to read silver iodide from an input file in E3SM simulations
  !
  !
  !-----------------Code history------------------------------------!
  ! Authors: Luke Van Roekel
  !
  !-----------------------------------------------------------------!
  use cam_abortutils,   only: endrun
  use pio
  use cam_pio_utils, only : cam_pio_openfile
  use infnan,       only: nan, assignment(=)
  use physics_types,          only: physics_state, physics_ptend
  use physics_buffer, only: physics_buffer_desc
  use shr_kind_mod,  only: r8=>shr_kind_r8, cx => SHR_KIND_CX, &
                           cs =>SHR_KIND_CS, cxx =>SHR_KIND_CXX
  use shr_log_mod,   only: errMsg=>shr_log_errMsg
  use shr_const_mod, only: shr_const_rdair, shr_const_g
  use phys_control,  only: phys_getopts
  use physconst,     only: pi
  use cam_logfile,   only: iulog
  use spmd_utils,    only: masterproc
  use perf_mod,      only: t_startf, t_stopf
  use mpishorthand
  use cam_history_support, only: fillvalue
  use ppgrid,         only: pver, pcols
  use input_data_utils, only: time_coordinate

  implicit none

  private
  save

  ! ----------------- !
  ! Public interfaces !
  ! ----------------- !

  public :: agi_read_data_init,  & ! initializes reading of file
            agi_read_data_adv  ! advances position in agi file 
       
  integer :: ntimes
  integer :: last_index = 1
  integer :: agi_vid

  logical :: initialized = .false.

  type :: agi_emis_native_grid
     !------------------------------------------------------------------
     !"agi_emis_native_grid" is forcing from files which has to be on the
     !native grid (only in horizontal,vertical resolution may be different
     !from the model's grid resolution)
     !That is, forcing files has to be on same grid as the grid used for
     !the model run
     !------------------------------------------------------------------

     !Number of levels in the 3D forcing file
     integer                               :: lev_frc

     !Data structure to store two time samples from a file to do time interpolation in the next step
     !(pcols,lev_frc,begchunk:endchunk,2)
     real(r8), pointer, dimension(:,:,:,:) :: native_grid_flds_tslices

     !Data structure to store data after time interpolation from two time samples
     !(pcols,lev_frc,begchunk:endchunk)
     real(r8), pointer, dimension(:,:,:)   :: native_grid_flds

     !Data structure to keep track of time
     type(time_coordinate) :: time_coord

     !agi emission name
     character( len = cx)  :: agi_emis_name_ngrd


     !Level bounds read from input file
     real(r8), pointer, dimension(:,:) :: lev_bnds

     !Forcing file name
     character( len = cx)  :: input_file

     !Units of forcing data
     character( len = cs)  :: units

     !logical to control first data read
     logical               :: initialized

     !pbuf index to store read in data in pbuf
     integer               :: pbuf_ndx = -1
  end type agi_emis_native_grid
  type(agi_emis_native_grid) :: native_grid_agi

  logical :: readiv = .false.
  logical :: has_fixed_ubs = .false.
  logical :: cam_outfld = .false.

  logical :: rmv_file = .false.
  logical :: dimnames_set = .false.

  integer :: number_flds
  character(len=8) :: dim1name, dim2name
  integer, parameter :: huge_int = huge(1)
contains

  subroutine agi_read_data_init(state, pbuf2d, agi_emis_data_file,agiEmis_idx)

    use cam_history,      only: addfld, add_default
    use physics_types,    only: physics_state
    use ppgrid,           only: pcols, begchunk, endchunk
    use physics_buffer,   only: physics_buffer_desc, pbuf_get_index
    use physics_types,    only: physics_state
    use physics_buffer,   only: physics_buffer_desc
    use cam_grid_support, only: cam_grid_id, cam_grid_check
    use cam_grid_support, only: cam_grid_get_dim_names
    use dyn_grid,         only: get_horiz_grid_dim_d
    use dycore,           only: dycore_is
    use cam_pio_utils,    only: cam_pio_openfile
    use pio,              only: file_desc_t, pio_nowrite, pio_closefile, pio_inq_dimid, pio_bcast_error, &
         pio_seterrorhandling, pio_noerr, pio_inquire_dimension, pio_get_att, pio_inq_varid, pio_get_var


    implicit none

    !arguments
    type(physics_state), intent(in)    :: state(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    !local vars
    type(file_desc_t)   :: fh
    character(len=16)   :: spc_name
    character(len=cxx)  :: err_str
    character(len=cx),intent(in)    :: agi_emis_data_file
    character(len=cx) :: input_file
    integer :: ndx, istat, i, astat, m, n, mm, c
    integer :: grid_id, agiEmis_idx
    integer :: dimlevid, var_id, errcode, dim1id, dim2id, dim1len, dim2len
    integer :: dimbndid, nbnd
    integer :: hdim1_d, hdim2_d    ! model grid size

    real(r8) :: dtime

    if (.not. dimnames_set) then
       grid_id = cam_grid_id('physgrid')
       if (.not. cam_grid_check(grid_id)) then
          call endrun('no "physgrid" grid:'//errmsg(__FILE__,__LINE__))
       endif
       !dim1name and dim2name are populated here with the grid dimension the model is running on (e.g. ne30, lat, lon etc.)
       !For SE grid, dim1name = dim2name = "ncol"
       !For FV grid, dim1name = lon, dim2name = lat
       call cam_grid_get_dim_names(grid_id, dim1name, dim2name)
       if(dim1name .ne. dim2name) then
          call endrun('Must run on SE grid, dim1 and dim2 for agi_read suggest lat/lon or FV grid')
       endif
       dimnames_set = .true.
    end if

    !--------------------------------------------------------------------------------
    ! allocate forcings type array for native grid forcing files
    !--------------------------------------------------------------------------------
!    allocate( native_grid_agi, stat=astat )
!    if( astat /= 0 ) then
!       write(err_str,*) 'failed to allocate native_grid_agi array; error = ',astat,',',errmsg(__FILE__, __LINE__)
!       call endrun(err_str)
!    end if

    !-----------------------------------------------------------------------
    !       initialize variables for native grid forcing files          
    !-----------------------------------------------------------------------

    input_file     = agi_emis_data_file 
    native_grid_agi%input_file = input_file

    native_grid_agi%initialized    = .false.
    dtime = 0.0_r8 
    call native_grid_agi%time_coord%initialize(trim(adjustl(input_file)), &
               force_time_interp=.true., delta_days=dtime)

    !-----------------------------------------------------------------------
    !       Open file
    !-----------------------------------------------------------------------
    call cam_pio_openfile(fh, trim(adjustl(input_file)), PIO_NOWRITE)

   !ask PIO to return the control if it experiences an error so that we can
    !handle it explicitly in the code
    call pio_seterrorhandling(fh, pio_bcast_error)
   
    !-----------------------------------------------------------------------
    !       Sanity checks for the native grid
    !-----------------------------------------------------------------------

    !if forcing file is on a different grid than the model grid
    !(e.g. model is running on an FV grid and forcing netcdf file is on an SE grid), exit with an error
    if(pio_inq_dimid(fh, trim(adjustl(dim2name)), dim1id) /= pio_noerr) then
    !pio_inq_dimid function tries to find dim1name in file with id "fh"
    !if it can't find dim1name, it means there is a mismacth in model and netcdf
    !file grid
       call endrun('grid mismatch, failed to find '//dim1name//' dimension in file:'&
            &' '//trim(input_file)//' '&
            &' '//errmsg(__FILE__,__LINE__))
    endif

    !find if the model and netcdf file has same grid resolution
    call get_horiz_grid_dim_d(hdim1_d,hdim2_d) !get model dim lengths
    if( dycore_is('SE') )  then
       if(pio_inquire_dimension(fh, dim1id, len = dim1len) ==  pio_noerr) then
          if(dim1len /= hdim1_d ) then !compare model grid length with file's
             write(err_str,*)'Netcdf file grid size(',dim1len,') should be same as model grid size(',&
                  hdim1_d,'), netcdf file is:'//trim(input_file)
             call endrun(err_str//errmsg(__FILE__,__LINE__))
          endif
       else
          call endrun('failed while inquiring dimensions of file:'//trim(adjustl(native_grid_agi%input_file))//' '&
               &' '//errmsg(__FILE__,__LINE__))
       endif
    elseif( dycore_is('LR')) then
       if(pio_inq_dimid(fh, trim(adjustl(dim2name)), dim2id) .ne. pio_noerr) then !obtain lat dimension of model
          call endrun('failed while inquiring dimension'//trim(adjustl(dim2name))//' from file:'&
               &' '//trim(adjustl(native_grid_agi%input_file))//' '//errmsg(__FILE__,__LINE__))
       endif
       if(pio_inquire_dimension(fh, dim1id, len = dim1len) ==  pio_noerr .and. &
          pio_inquire_dimension(fh, dim2id, len = dim2len) ==  pio_noerr) then !compare grid and model's dims
          if(dim1len /= hdim1_d .or. dim2len /= hdim2_d)then
             write(err_str,*)'Netcdf file grid size(',dim1len,' x ',dim2len,') should be same as model grid size(',&
                   hdim1_d,' x ',hdim2_d,'), netcdf file is:'//trim(adjustl(input_file))
             call endrun(err_str//errmsg(__FILE__,__LINE__))
          endif
       else
          call endrun('failed while inquiring dimensions of file:'//trim(adjustl(input_file))//' '&
                       &' '//errmsg(__FILE__,__LINE__))
       endif
    else
       call endrun('Only SE or LR(FV) grids are supported currently:'//errmsg(__FILE__,__LINE__))
    endif

    !Find the value of vertical levels in the forcing file
    if( pio_inq_dimid(fh, 'lev', dimlevid) ==  pio_noerr ) then
       if ( pio_inquire_dimension(fh, dimlevid, len =  native_grid_agi%lev_frc) /=  pio_noerr ) then
          write(err_str,*)'failed to obtain value of "lev" dimension from file:',&
                           trim(adjustl(input_file)),',',errmsg(__FILE__, __LINE__)
          call endrun(err_str)
       endif
       !obtain level bounds needed for vertical interpolation
       if( pio_inq_varid(fh, 'lev_bnds', var_id) ==  pio_noerr ) then
          !get dimension "bound"
          if( pio_inq_dimid(fh, 'bound', dimbndid) ==  pio_noerr ) then
             if ( pio_inquire_dimension(fh, dimbndid, len = nbnd) ==  pio_noerr ) then
                !"nbnd" has to be 2 (it is obvious but adding a check here doesn't hurt)
                if(nbnd /= 2) then
                   write(err_str,*)'"bound" should be equal to 2, bound=',nbnd,' in file:', &
                                    trim(adjustl(input_file)),',',errmsg(__FILE__, __LINE__)
                   call endrun(err_str)
                endif
                allocate(native_grid_agi%lev_bnds(nbnd,native_grid_agi%lev_frc))
                if (pio_get_var(fh, var_id,native_grid_agi%lev_bnds) /=  pio_noerr ) then
                   write(err_str,*)'failed to read "lev_bnds" variable from file:',&
                                    trim(adjustl(input_file)),',',errmsg(__FILE__, __LINE__)
                   call endrun(err_str)
                endif
             else
                write(err_str,*)'failed to obtain value of "bound" dimension from file:',&
                                 trim(adjustl(input_file)),',',errmsg(__FILE__, __LINE__)
                call endrun(err_str)
             endif
          else
             write(err_str,*)'failed to inquire "bound" dimension from file:',&
                              trim(adjustl(input_file)),',',errmsg(__FILE__, __LINE__)
             call endrun(err_str)
          endif
       else
          write(err_str,*)'failed to obtain "lev_bnds" variable from file:',&
                           trim(adjustl(input_file)),',',errmsg(__FILE__, __LINE__)
          call endrun(err_str)
       endif
    else
       write(err_str,*)'Dimension "lev" is not found in:',&
                        trim(adjustl(input_file)),',',errmsg(__FILE__, __LINE__)
       call endrun(err_str)
    endif

    !close file
    call pio_closefile(fh)

    !allocate arrays to stroe data for interpolation in time
    allocate(native_grid_agi%native_grid_flds_tslices(pcols, native_grid_agi%lev_frc,begchunk:endchunk,2), stat=astat )

    if( astat/= 0 ) then
       write(err_str,*) 'failed to allocate native_grid_agi%native_grid_flds_tslices array; '&
                       &'error = ',astat,',',errmsg(__FILE__, __LINE__)
       call endrun(err_str)
    endif
    !allocate arrays to hold data before the vertical interpolation
    allocate(native_grid_agi%native_grid_flds(pcols, native_grid_agi%lev_frc,begchunk:endchunk), stat=astat )
    if( astat/= 0 ) then
       write(err_str,*) 'failed to allocate native_grid_agi%native_grid_flds array; error = ',&
                         astat,',',errmsg(__FILE__, __LINE__)
       call endrun(err_str)
    endif

    !get pbuf index to store the field in pbuf
    native_grid_agi%pbuf_ndx = agiEmis_idx 
    if(errcode < 0 ) then
       write(err_str,*)'failed to get pbuf index for AgI errorcode is:',errcode,',',errmsg(__FILE__, __LINE__)
       call endrun(err_str)
    endif

!    call addfld( 'AgI_emis', (/ 'lev' /), 'A',  'kg/s',     &
!            'Silver Iodide Emission rate ' )
!    call add_default( 'AgI_emis', 1, ' ' )

    number_flds = 0
    if (associated(native_grid_agi%native_grid_flds_tslices))  number_flds = 1
    !read the forcing file once to initialize variable including time cordinate
    call advance_native_grid_data_agi( native_grid_agi )
    native_grid_agi%initialized = .true.

    if( number_flds < 1 ) then
       if ( masterproc ) then
          write(err_str,*) 'There are no AgI emissions ',errmsg(__FILE__, __LINE__)
          call endrun(err_str)
       endif
    end if

   end subroutine agi_read_data_init

   subroutine agi_read_data_adv(state, pbuf2d, agitend)

    use perf_mod,     only: t_startf, t_stopf
    use tracer_data,  only: advance_trcdata
    use physics_types,only: physics_state
    use ppgrid,       only: pcols, pver, begchunk, endchunk
    use string_utils, only: to_lower, GLC
    use cam_history,  only: outfld
    use physconst,    only: mwdry       ! molecular weight dry air ~ kg/kmole
    use physconst,    only: boltz                ! J/K/molecule
    ! C.-C. Chen
    use physics_buffer, only : physics_buffer_desc, pbuf_get_field, pbuf_get_chunk

    implicit none

    type(physics_state), intent(in)    :: state(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    real(r8), dimension(pcols,pver,begchunk:endchunk), intent(out) :: agitend
    integer  :: ind, c, ncol, i, caseid, m, pbuf_ndx
    real(r8) :: to_mmr(pcols,pver)

    real(r8),pointer :: tmpptr(:,:)
    character(len = cs) :: units_spc

    call t_startf('agi_read_data_adv')

    pbuf_ndx = native_grid_agi%pbuf_ndx

    call advance_native_grid_data_agi( native_grid_agi )

    !    call vert_interp_agi( state, pbuf_ndx, native_grid_agi, pbuf2d)

	call t_stopf('agi_read_data_adv')

    agitend(:,:,:) = native_grid_agi%native_grid_flds(:,:,:)

   end subroutine agi_read_data_adv

   subroutine advance_native_grid_data_agi( native_grid_strct )

     use ppgrid,         only: pcols, begchunk, endchunk
    use ncdio_atm,      only: infld
    use cam_pio_utils,  only: cam_pio_openfile
    use pio,            only: file_desc_t, pio_nowrite, pio_closefile

    implicit none

    !args
    type(agi_emis_native_grid), intent (inout) :: native_grid_strct

    !local vars
    type(file_desc_t) :: fh
    character(len=cs) :: spc_name

    logical  :: read_data
    integer  :: indx2_pre_adv
    logical  :: found

    !Decide whether to read new data or not (e.g. data may needs to be read on month boundaries )
    read_data = native_grid_strct%time_coord%read_more() .or. .not. native_grid_strct%initialized

    !Find time index to decide whether to read new data or recycle previously read data
    indx2_pre_adv = native_grid_strct%time_coord%indxs(2)

    !compute weights for time interpolation (time_coord%wghts) by advancing in time
    call native_grid_strct%time_coord%advance()

    if ( read_data ) then

       !open file
       call cam_pio_openfile(fh, trim(adjustl(native_grid_strct%input_file)), PIO_NOWRITE)

       ! read time-level 1
       if (native_grid_strct%initialized .and. native_grid_strct%time_coord%indxs(1) == indx2_pre_adv) then
          ! skip the read if the needed vals for time level 1 are present in time-level 2
          native_grid_strct%native_grid_flds_tslices(:,:,:,1) = native_grid_strct%native_grid_flds_tslices(:,:,:,2)
       else
          !NOTE: infld call doesn't do any interpolation in space, it just reads in the data
          call infld('AgI_emis', fh, 'ncol', 'ncol', 'lev',&
            1, pcols, 1, native_grid_strct%lev_frc, begchunk,endchunk, &
               native_grid_strct%native_grid_flds_tslices(:,:,:,1), found, &
               gridname='physgrid', timelevel=native_grid_strct%time_coord%indxs(1))
          if (.not. found) then
             call endrun('AgI_emis not found '//errmsg(__FILE__,__LINE__))
          endif
       endif

       ! read time level 2
       call infld('AgI_emis', fh, 'ncol', 'ncol','lev',&
         1, pcols, 1, native_grid_strct%lev_frc, begchunk, endchunk, &
            native_grid_strct%native_grid_flds_tslices(:,:,:,2), found, &
            gridname='physgrid', timelevel=native_grid_strct%time_coord%indxs(2))

       if (.not. found) then
          call endrun('AgI_emis' // ' not found '//errmsg(__FILE__,__LINE__))
       endif

       !close file
       call pio_closefile(fh)
    endif

    ! interpolate between time-levels
    ! If time:bounds is in the dataset, and the dataset calendar is compatible with EAM's,
    ! then the time_coordinate class will produce time_coord%wghts(2) == 0.0,
    ! generating fluxes that are piecewise constant in time.

    print *, 'agi1 = ',maxval(native_grid_strct%native_grid_flds_tslices(:,:,:,1)), &
                       minval(native_grid_strct%native_grid_flds_tslices(:,:,:,1))
    print *, 'agi2 = ',maxval(native_grid_strct%native_grid_flds_tslices(:,:,:,2)), &
                       minval(native_grid_strct%native_grid_flds_tslices(:,:,:,2))
    print *, 'agiloc = ',maxloc(native_grid_strct%native_grid_flds_tslices(:,:,:,1))
    print *, 'weights = ', native_grid_strct%time_coord%wghts(:)
    print *, 'index = ',native_grid_strct%time_coord%indxs(2), native_grid_strct%time_coord%indxs(1)
    if (native_grid_strct%time_coord%wghts(2) == 0.0_r8) then
       native_grid_strct%native_grid_flds(:,:,:) = native_grid_strct%native_grid_flds_tslices(:,:,:,1)
    else
       native_grid_strct%native_grid_flds(:,:,:) = native_grid_strct%native_grid_flds_tslices(:,:,:,1) + &
            native_grid_strct%time_coord%wghts(2) * (native_grid_strct%native_grid_flds_tslices(:,:,:,2) - &
            native_grid_strct%native_grid_flds_tslices(:,:,:,1))
    endif

   end subroutine advance_native_grid_data_agi

   subroutine vert_interp_agi( state, pbuf_ndx, native_grid_strct, pbuf2d )

    !-------------------------------------------------------------------
    !    This subroutine interpolates in vertical direction. Finally the
    !    interpolated data is stored  in pbuf
    ! called by: aircraft_emit_adv (local)
    !-------------------------------------------------------------------

    use physics_types,  only: physics_state
    use physics_buffer, only: physics_buffer_desc, pbuf_get_field, pbuf_get_chunk
    use ppgrid,         only: pver, pcols, begchunk, endchunk
    use cam_history,    only: outfld
    use constituents,   only: cnst_name
    use co2_cycle,      only: c_i


    !args
    type(physics_state), intent(in)        :: state(begchunk:endchunk)
    type(physics_buffer_desc), pointer     :: pbuf2d(:,:)
    type(agi_emis_native_grid), intent(in) :: native_grid_strct
    integer, intent(in)                    :: pbuf_ndx

    !local vars
    type(physics_buffer_desc), pointer :: pbuf_chnk(:)

    integer  :: ic, ncol, icol
    real(r8) :: vrt_interp_field(pcols, pver,begchunk:endchunk)
    real(r8),pointer :: tmpptr_native_grid(:,:)

    real(r8) :: nzil, nzir, mzil, mzir                       ! level bounds
    real(r8) :: ovrl, ovrr, ovrf                             ! overlap bounds
    real(r8) :: ovrmat(pver, native_grid_strct%lev_frc) ! overlap fractions matrix
    integer  :: kinp, kmdl                              ! vertical indexes
    integer  :: lchnk, k
    real(r8) :: ftem(pcols, begchunk:endchunk)

    ! Vertical interpolation follows Joeckel (2006) ACP method for conservation
    do ic = begchunk,endchunk
       ncol = state(ic)%ncol
       do icol = 1,ncol
       do kinp = 1, native_grid_strct%lev_frc
          nzil = native_grid_strct%lev_bnds(1,kinp)
          nzir = native_grid_strct%lev_bnds(2,kinp)
          do kmdl = 1, pver
             mzil               = state(ic)%zi(icol, kmdl)
             mzir               = state(ic)%zi(icol, kmdl+1)
             ovrl               = MAX( MIN(nzil, nzir), MIN(mzil, mzir) )
             ovrr               = MIN( MAX(nzil, nzir), MAX(mzil, mzir) )
             ovrf               = INT( SIGN(0.5_r8, ovrr-ovrl) + 1._r8 )
             ovrmat(kmdl,kinp) = ABS(ovrr - ovrl) * ovrf / ABS(nzir - nzil)
          end do
       end do
       vrt_interp_field(icol,:,ic)  = MATMUL( ovrmat, native_grid_strct%native_grid_flds(icol,:,ic) )
       end do

    ftem(:ncol, ic) = vrt_interp_field(:ncol, 1, ic)
    do k = 2, pver
       ftem(:ncol, ic) = ftem(:ncol, ic) + vrt_interp_field(:ncol, k, ic)
    end do
    end do
    !future work: these two (above abd below) do loops should be combined into one.

    !add field to pbuf
    !$OMP PARALLEL DO PRIVATE (IC, NCOL, tmpptr_native_grid, pbuf_chnk, vrt_interp_field)
    do ic = begchunk, endchunk
    ncol = state(ic)%ncol
    pbuf_chnk => pbuf_get_chunk(pbuf2d, ic)
    call pbuf_get_field(pbuf_chnk, pbuf_ndx, tmpptr_native_grid )
    tmpptr_native_grid(:ncol,:) = vrt_interp_field(:ncol,:,ic)
    end do
  end subroutine vert_interp_agi

end module agi_read
