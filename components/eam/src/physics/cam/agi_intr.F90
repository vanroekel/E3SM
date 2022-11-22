module agi_intr

  !-----------------------------------------------------------------!
  ! Module to include silver iodide in E3SM simulations
  !
  !
  !-----------------Code history------------------------------------!
  ! Authors: Luke Van Roekel
  !
  !-----------------------------------------------------------------!

  use shr_kind_mod, only: r8=>shr_kind_r8
  use shr_log_mod,  only: errMsg=>shr_log_errMsg
  use phys_control, only: phys_getopts
  use physconst,    only: pi
  use cam_logfile,  only: iulog
  use spmd_utils,   only: masterproc
  use perf_mod,     only: t_startf, t_stopf
  use mpishorthand
  use cam_history_support, only: fillvalue
  use ppgrid        only: pver, pcols

  implicit none

  private
  save

  ! ----------------- !
  ! Public interfaces !
  ! ----------------- !

  public :: agi_register,  & ! add silver iodide to pbuf
            agi_readnl,    & ! read namelist related to silver iodide
            agi_init_cnst, & ! initialize a constant point source of silver iodide
            agi_emiss_data, &   ! read emission from data with times.
            agi_pointsource

  integer :: agi_idx

  !namelist variables
  logical :: agi_enable = .false. ! AgI is off by default
  logical :: agi_point_source_emis = .true. ! if true use a fixed constant in time point source
  logical :: agi_data_emis = .false. ! if true do emission by data file
  logical :: agi_microphysics = .false. ! if true add impact of agi to MG2
  real(r8) :: agi_lat = huge(1.0_r8)
  real(r8) :: agi_lon = huge(1.0_r8)
  real(r8) :: agi_emis_rate = huge(1.0_r8) ! emission of AgI at a point source in kg/s
  character(len=*) :: agi_emis_data_file = 'agi_emission.bin'
  real(r8),parameter :: rad_to_deg = 180.0/pi

  contains

  subroutine agi_readnl(nlfile)

    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit

    character(len=*), intent(in) :: nlfile
    character(len=*), parameter :: subname = 'agi_read_nl'

    integer :: unitn, ierr

    namelist /agi_nl/ agi_enable, agi_point_source_emis, agi_data_emis,  &
                      agi_microphysics, agi_lat, agi_lon, agi_height, &
                      agi_thickness, agi_emis_rate, agi_emis_data_file

    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'agi_nl', status='ierr' )
       if (ierr == 0) then
          read(unitn, agi_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)

       if (agi_enable .and. agi_point_source_emis .and.  agi_data_emis) then
         call endrun(subname // ':: ERROR cannot enable both point source and data emissions')
       end if

       if (agi_enable .and. .not. agi_point_source_emis .and. .not. agi_data_emis) then
         call endrun(subname // ':: ERROR must enable either point source or data emissions for AgI')
       end if

    end if

#ifdef SPMD
    call mpibcast(agi_enable,            1, mpilog, 0, mpicom)
    call mpibcast(agi_point_source_emis, 1, mpilog, 0, mpicom)
    call mpibcast(agi_data_emis,         1, mpilog, 0, mpicom)
    call mpibcast(agi_microphysics,      1, mpilog, 0, mpicom)
    call mpibcast(agi_lat,               1, mpir8,  0, mpicom)
    call mpibcast(agi_lon,               1, mpir8,  0, mpicom)
    call mpibcast(agi_height,            1, mpir8,  0, mpicom)
    call mpibcast(agi_thickness,         1, mpir8,  0, mpicom)
    call mpibcast(agi_emis_rate,         1, mpir8,  0, mpicom)
    call mpibcast(agi_emis_data_file,    128, mpichar, 0, mpicom)
#endif

  end subroutine agi_readnl

!=============================================================================================

  subroutine agi_register

    ! Register Agi in the physics buffer

    call pbuf_add_field('AgI', 'global', dtype_r8, (/pcols,pver/), agi_idx)

  end subroutine agi_register

!==============================================================================================

  subroutine agi_pointsource(state)

    use ppgrid only: pver, pcols

    type(physics_state), intent(in)    :: state
    real(r8), pointer, dimension(:,:) :: agi
    real(r8), parameter :: rtd = 180.0_r8 / pi
    integer :: ncol, lchnk, agi_indx, itim_old
    real(r8) :: agiTemp
    use rgrid,          only: nlon
  real(r8) :: pdel(plon,plev)           ! Pressure depth of layer
  real(r8) :: pint(plon,plevp)          ! Pressure at interfaces
  real(r8) :: pmid(plon,plev)           ! Pressure at midpoint

    ncol = state%ncol
    lchnk = state%ncol

    itim_old = pbuf_old_tim_idx() !gets time index 

    do k = 1, pver
       z3(:ncol,k) = state%zm(:ncol,k) + state%phis(:ncol)*rga
    end do

    call plevs0(nlon(lat), plon, plev, ps(1,lat,n3), pint, pmid, pdel)
	agi_indx = pbuf_get_index('AgI')
    agi = pbug_get_field(pbuf, agi_indx, agi, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
    !find points to release AgI
    do icol = 1, ncol
       do k = 1, pver
          agi(icol,k) = 0.0_r8
          if(abs(state%lon(icol)*rtd - agi_lon) < 0.5_r8 .and. abs(state%lat(icol)*rtd - agi_lat) < 0.5_r8 &
             .and. abs(stat%zm(icol, k) - agi_height) < agi_thickness) then
             agiMass = agi_emis_rate* !need conversion of mass rate to number concentration
! use mass of dry air and then ideal gass law at a given pressure to get number conc


  end subroutine agi_pointsource

