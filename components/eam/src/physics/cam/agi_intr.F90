module agi_intr

  !-----------------------------------------------------------------!
  ! Module to include silver iodide in E3SM simulations
  !
  !
  !-----------------Code history------------------------------------!
  ! Authors: Luke Van Roekel
  !
  !-----------------------------------------------------------------!

  use physics_types,          only: physics_state, physics_ptend
  use physics_buffer, only: physics_buffer_desc
  use shr_kind_mod,  only: r8=>shr_kind_r8
  use shr_log_mod,   only: errMsg=>shr_log_errMsg
  use shr_const_mod, only: shr_const_rdair, shr_const_g
  use phys_control,  only: phys_getopts
  use physconst,     only: pi
  use cam_logfile,   only: iulog
  use constituents,   only: pcnst
  use spmd_utils,    only: masterproc
  use perf_mod,      only: t_startf, t_stopf
  use mpishorthand
  use cam_history_support, only: fillvalue
  use ppgrid,         only: pver, pcols

  implicit none

  private
  save

  ! ----------------- !
  ! Public interfaces !
  ! ----------------- !

  public :: agi_register,  & ! add silver iodide to pbuf
            agi_readnl,    & ! read namelist related to silver iodide
            agi_emiss_data, &   ! read emission from data with times.
            agi_pointsource, &
            agi_tend,  &
            agi_ini

  !namelist variables
  logical :: agi_enable = .false. ! AgI is off by default
  logical :: agi_point_source_emis = .true. ! if true use a fixed constant in time point source
  logical :: agi_data_emis = .false. ! if true do emission by data file
  logical :: agi_microphysics = .false. ! if true add impact of agi to MG2
  logical :: lq(pcnst)
  real(r8) :: agi_lat = huge(1.0_r8)
  real(r8) :: agi_lon = huge(1.0_r8)
  real(r8) :: agi_emis_rate = huge(1.0_r8) ! emission of AgI at a point source in kg/s
  real(r8) :: agi_height = 0.0_r8
  real(r8) :: agi_thickness = 100.0_r8
  character(len=256) :: agi_emis_data_file
  real(r8),parameter :: rad_to_deg = 180.0/pi
  real(r8),dimension(:,:),allocatable :: massAtmos, densAtmos

  integer :: agi_idx, agiq_idx

  integer, public :: &
    ix_qagi = 0,         &
    ix_agi = 0

  contains

  subroutine agi_emiss_data()

  end subroutine agi_emiss_data
  subroutine agi_readnl(nlfile)

    use cam_abortutils,  only: endrun
    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit

    character(len=*), intent(in) :: nlfile
    character(len=*), parameter :: subname = 'agi_read_nl'

    integer :: unitn, ierr

    namelist /agi_ctl_nl/ agi_enable, agi_point_source_emis, agi_data_emis,  &
                      agi_microphysics, agi_lat, agi_lon, agi_height, &
                      agi_thickness, agi_emis_rate, agi_emis_data_file

    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'agi_ctl_nl', status=ierr )
       if (ierr == 0) then
          read(unitn, agi_ctl_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       else
          call endrun(subname // ':: ERROR namelist group not found')
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

    use cam_history,            only: addfld
    use physics_buffer,         only: pbuf_add_field, dtype_r8
    use constituents,   only: cnst_add
    ! Register Agi in the physics buffer

    if (.not. agi_enable) return

    call pbuf_add_field('AgI_nadv', 'global', dtype_r8, (/pcols,pver/), agi_idx)
    call pbuf_add_field('qAgI_nadv', 'global', dtype_r8, (/pcols,pver/), agiq_idx)
    call cnst_add('qAgI',0._r8,0._r8,0._r8,ix_qagi,longname='mixing ratio of AgI',cam_outfld=.false.)
    call cnst_add('AgI',0._r8,0._r8,0._r8,ix_agi,longname='number concentration of AgI', cam_outfld=.false.)

  end subroutine agi_register

!=============================================================================================

  subroutine agi_ini(pbuf2d)

    use cam_history,            only: addfld, add_default
    use physics_buffer,         only: pbuf_get_index, pbuf_add_field, pbuf_set_field, physics_buffer_desc
    implicit none
    !  Input Variables
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    integer :: agi_idx, agiq_idx

    if (.not. agi_enable) return

    ! Register Agi in the physics buffer

    agi_idx = pbuf_get_index('AgI_nadv')
    agiq_idx = pbuf_get_index('qAgI_nadv')

    lq(1:pcnst) = .true.

    call pbuf_set_field(pbuf2d, agi_idx, 0.0_r8)
    call pbuf_set_field(pbuf2d, agiq_idx, 0.0_r8)
    call addfld('agi_Nc', (/'ilev'/), 'A', 'molec/m3', 'number concentration of silver iodide')
    call addfld('agi_q', (/'ilev'/), 'A', 'kg/kg', 'mixing ratio of silver iodide')

    call add_default('agi_Nc',1,' ')
    call add_default('agi_q',1,' ')


  end subroutine agi_ini

!==============================================================================================
  subroutine agi_tend(state1, ptend1, pbuf1, ztodt)

    use ppgrid, only: pver, pcols
    use rgrid, only: nlon
    use physics_types,          only: physics_state, physics_ptend
    use physics_buffer, only: physics_buffer_desc
    use dyn_grid,         only: get_horiz_grid_dim_d

    type(physics_state), intent(in) :: state1
    type(physics_ptend), intent(inout) :: ptend1
    type(physics_buffer_desc), pointer, intent(in) :: pbuf1(:)
    real(r8), intent(in) :: ztodt
    integer :: plon, plat

    if(.not. agi_enable) return

    call get_horiz_grid_dim_d( plon, plat )
    if(agi_point_source_emis) then
      call agi_pointsource(state1, ptend1, pbuf1, ztodt, plon, plat)
    elseif(agi_data_emis) then
      !stub non-call to agi data emis
    else
      !do nothing
    endif

    ! this is probably where the microphysics part comes to add to ice nuclei, still need a PBUF field (FIXME)

  end subroutine agi_tend

!==============================================================================================

  subroutine agi_pointsource(state, ptend, pbuf, dt, plon, plat)

    use ppgrid, only: pver, pcols
    use cam_history, only: outfld
    use rgrid,          only: nlon
    use physconst,          only: gravit, rga, rair, cpair, latvap, rearth, pi, cappa
    use physics_types,          only: physics_state, physics_ptend, physics_ptend_init
    use physics_buffer, only: pbuf_old_tim_idx, physics_buffer_desc, pbuf_get_field, pbuf_get_index
    use pmgrid,           only: plev
    use phys_grid,        only : get_area_all_p

    type(physics_buffer_desc), pointer :: pbuf(:)
    type(physics_state), intent(in)    :: state
    type(physics_ptend), intent(inout) :: ptend
    real(r8), intent(in) :: dt
    integer, intent(in) :: plon,plat
    real(r8), pointer, dimension(:,:) :: agi,agiq
    real(r8), parameter :: rtd = 180.0_r8 / pi
    integer :: k, icol, ncol, lchnk, itim_old, agiq_indx
    real(r8) :: nConcTend, mixingratioTend, agiTemp, massAtmos, densAtmos, volAtmos
    real(r8) :: pdel(plon,plev)           ! Pressure depth of layer
    real(r8) :: pint(plon,plev)          ! Pressure at interfaces
    real(r8) :: pmid(plon,plev)           ! Pressure at midpoint
    real(r8) :: z3(pcols,pver)
    real(r8), allocatable, dimension(:) :: area
    real(r8), parameter :: molec_agi = 0.23477 !kg/mol
    ncol = state%ncol
    lchnk = state%lchnk 

    call physics_ptend_init(ptend, state%psetcols, 'agi', lq=lq)

    itim_old = pbuf_old_tim_idx() !gets time index 

    allocate(area(ncol))
    do k = 1, pver
       z3(:ncol,k) = state%zm(:ncol,k) + state%phis(:ncol)*rga
    end do

    call get_area_all_p(lchnk, ncol, area)
    area = area * rearth**2

    !use delp to get mass of the air in a layer should be dp/g*dA
!    call plevs0(nlon(lat), plon, plev, ps(1,lat,n3), pint, pmid, pdel)
!	agi_indx = pbuf_get_index('AgI')
    call pbuf_get_field(pbuf, agi_idx, agi, start=(/1,1/), kount=(/pcols,pver/))
!    agiq_indx = pbuf_get_index('qAgI')
    call pbuf_get_field(pbuf, agiq_idx, agiq, start=(/1,1/), kount=(/pcols,pver/))
    !find points to release AgI
    
    do icol = 1, ncol
       do k = 1, pver
          agi(icol,k) = state%q(icol,k,ix_agi)
          agiq(icol,k) = state%q(icol,k,ix_qagi)
 !         agi(icol,k) = 0.0_r8
 !         agiq(icol,k) = 0.0_r8
          if(abs(state%lon(icol)*rtd - agi_lon) < 2.5_r8 .and. abs(state%lat(icol)*rtd - agi_lat) < 2.5_r8 &
             .and. abs(state%zm(icol, k) - agi_height) < agi_thickness) then

             !compute mass of atmosphere in a grid cell
             massAtmos = 1.0_r8/shr_const_g*state%pdel(icol,k)
             densAtmos = state%pmid(icol,k)/shr_const_rdair/state%t(icol,k)
             volAtmos = massAtmos/densAtmos
             mixingratioTend = agi_emis_rate/massAtmos
             nconcTend = agi_emis_rate/molec_agi*6.022e23_r8/volAtmos

             agi(icol,k) = agi(icol,k) + dt*nconcTend
             agiq(icol,k) = agiq(icol,k) + dt*mixingratioTend
 ptend%q(icol,k,ix_agi) = nconcTend
ptend%q(icol,k,ix_qagi) = mixingratioTend

           end if
       end do
    end do

    ! calls to outfield to save data
    call outfld('agi_Nc', agi,pcols,lchnk)
    call outfld('agi_q', agiq,pcols,lchnk)

    deallocate(area)
  end subroutine agi_pointsource

end module agi_intr
