!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vmix_cvmix
!
!> \brief MPAS ocean vertical mixing interface to new CVMix that runs on GPU
!> \author Luke Van Roekel
!> \date   16 January 2019
!> \details
!>  This module contains the routines for calls into CVMix-Kokkos
!>
!
!-----------------------------------------------------------------------

module ocn_vmix_cvmix2_kokkos

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timer
   use mpas_constants
   use mpas_log

   use ocn_constants
   use ocn_config
   use ocn_mesh
   use ocn_diagnostics_variables

   use ocn_constants

   use,intrinsic :: iso_c_binding

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_vmix_coefs_cvmix2_kokkos_build, &
             ocn_vmix_cvmix2_kokkos_init,        &
             ocn_vmix_cvmix2_kokkos_finalize

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   type, bind(c) :: CVMix_mixingOpts
     real (C_DOUBLE) :: bl1, bl2, bl3, bl4, prandtl
     real (C_DOUBLE) :: PP_nu_zero, PP_alpha, PP_exp, PP_nu_b, PP_kappa_b, KPP_nu_zero, KPP_Ri_zero, KPP_exp
     real (c_double) :: convect_diff, convect_visc, BVsqr_convect
     logical (C_BOOL) :: bryanLewis, lBruntVaisala, lShearKPP, lShearPP, &
                         lEkman, lMonOb, lKPP, computelangmuir
     logical (C_BOOL) :: LANGMUIR_LWF16, LANGMUIR_lifk17, parabolicNonLocal
     real (c_double) :: Ri_crit, surf_layer_ext, c_s, minVtsqr, rho_o, minOSBLUnderIce
     real (c_double) :: a_s, c_m, a_m, non_local_coeff, numRiSmoothPasses
   end type CVMix_mixingOpts

   INTERFACE
     SUBROUTINE f_cvmix2_dummy(nCells,nVertLevels,layerThickness,mixingOpts) &
                     BIND(c, NAME='c_compute_cvmix2_mixing_dummy')

             USE mpas_kind_types
             USE, INTRINSIC :: ISO_C_BINDING
             import :: CVMix_mixingOpts
             IMPLICIT none
             integer ( c_int ), value :: nCells,nVertLevels
             real(kind=c_double) :: layerThickness(nVertLevels,nCells)
             type (CVMix_mixingOpts) :: mixingOpts
     END SUBROUTINE f_cvmix2_dummy
   END INTERFACE


   INTERFACE
     SUBROUTINE f_kokkos_initialize(nVertLevels,nCols) &
       BIND(c, NAME='c_kokkos_initialize')
       USE, INTRINSIC :: ISO_C_BINDING
       IMPLICIT none
       integer(c_int),value :: nCols,nVertLevels
     END SUBROUTINE f_kokkos_initialize
   END INTERFACE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   INTERFACE
     SUBROUTINE f_kokkos_finalize() &
       BIND(c, NAME='c_kokkos_finalize')
       USE, INTRINSIC :: ISO_C_BINDING
       IMPLICIT NONE
     END SUBROUTINE f_kokkos_finalize
   END INTERFACE

   interface
     subroutine f_cvmix2_compute_mixing(nCells,nVertLevels,visc,diff, layerThickness, U, V, &
       NLT, N2, RI_grad, rho, fCor, sfcBuoyancyForcing, sfcFrictionVelocity, maxLevel, &
       OSBL, efactor, lasl, wind10, mixingOpts)                                             &
       BIND(c, NAME='c_compute_cvmix2_mixing') 

        use mpas_kind_types
        use,intrinsic :: iso_c_binding
        import :: CVMix_mixingOpts
        implicit none
        integer ( c_int ), value :: nCells,nVertLevels
        real(kind=RKIND) :: visc(nVertLevels+1,nCells), diff(nVertLevels+1,nCells),      &   
           layerThickness(nVertLevels,nCells), NLT(nVertLevels+1,nCells),       &
           N2(nVertLevels+1,nCells),                 &
           Ri_grad(nVertLevels+1,nCells), sfcBuoyancyForcing(nCells),           &
           sfcFrictionVelocity(nCells), U(nVertLevels, nCells), V(nVertLevels,nCells),  &
           rho(nVertLevels,nCells), fCor(nCells),                  &
           OSBL(nCells), efactor(nCells),lasl(nCells),   &   
           wind10(nCells)
        integer (c_int) :: maxLevel(nCells)
        type (CVMix_mixingOpts) :: mixingOpts
     end subroutine f_cvmix2_compute_mixing 
   end interface

   logical :: cvmixOn, cvmixBackgroundOn, cvmixConvectionOn, cvmixKPPOn
   real (kind=RKIND) :: backgroundVisc, backgroundDiff

   type (CVMix_mixingOpts) :: mixingOpts

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vmix_coefs_cmvix_build
!
!> \brief   Computes mixing coefficients using CVMix
!> \author  Todd Ringler
!> \date    04 February 2013
!> \details
!>  This routine computes the vertical mixing coefficients for momentum
!>  and tracers by calling CVMix routines.
!
!-----------------------------------------------------------------------

   subroutine ocn_vmix_coefs_cvmix2_kokkos_build(statePool, forcingPool, err, timeLevelIn)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

! Things needed
! nCells, nVertLevs, viscosity array, diff array, layer thickness
! velocityZonal, velocityMeridional, nonlocal T, nonLocal S, BruntVaisala, gradient rich, p_rho, 
! fcor, sfcBuoy, sfcFric, maxLevelCell, OSBLdepth, langmuir efactor, lasl, 10m wind, mixingOpts
! efactor and lasl sent in doesn't matter what's inn them.

      integer, intent(in), optional :: timeLevelIn !< Input: time level for state pool

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: &
         statePool         !< Input/Output: state information

      type (mpas_pool_type), intent(inout) :: &
         forcingPool   !< Input/Output: forcing information

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), pointer :: &
        iceFraction, windSpeed10m, lasl

      real (kind=RKIND), dimension(:,:), pointer :: &
         layerThickness

      real (kind=RKIND), dimension(:,:,:), pointer :: vertNonLocalFlux
      integer, pointer :: index_vertNonLocalFluxTemp

      integer :: k, i, iCell, timeLevel, nCells

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing mixing-related fields
      ! note that the user can choose multiple options and the
      !   mixing fields have to be added/merged together
      !
      !-----------------------------------------------------------------

      !
      ! assume no errors during initialization and set to 1 when error is encountered
      !
      err=0

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      !
      ! only build up viscosity/diffusivity if CVMix is turned on
      !
      if ( .not. cvmixOn ) return

      !
      ! set parameters
      !
      !
      ! set pointers for fields related to position on sphere
      !

      !
      ! set pointers for fields related to vertical mesh
      !
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)
      !
      ! set pointers for fields related to ocean forcing state
      !
      call mpas_pool_get_array(forcingPool, 'iceFraction', iceFraction)
      call mpas_pool_get_array(forcingPool, 'windSpeed10m', windSpeed10m)

      nCells = nCellsAll

      !$omp do schedule(runtime)
      do iCell = 1, nCells
         vertViscTopOfCell(:, iCell) = 0.0_RKIND
         vertDiffTopOfCell(:, iCell) = 0.0_RKIND
         vertNonLocalFlux(:, :, iCell) = 0.0_RKIND
         boundaryLayerDepth(iCell) = layerThickness(1,iCell)
      end do
      !$omp end do

      call mpas_timer_start('cvmix2.0 compute coefficients', .false.)
     
      call f_cvmix2_compute_mixing(nCells, nVertLevels, vertViscTopOfCell,      &
                  vertDiffTopOfCell, layerThickness, velocityZonal,       &
                  velocityMeridional, vertNonLocalFlux,            &
                  BruntVaisalaFreqTop,          &
                  RiTopOfCell, potentialDensity, fCell,                   &
                  surfaceBuoyancyForcing, surfaceFrictionVelocity, &
                  maxLevelCell, boundaryLayerDepth, iceFraction, lasl,    &
                  windSpeed10m, mixingOpts)

      call mpas_timer_stop('cvmix2.0 compute coefficients')

   !--------------------------------------------------------------------

   end subroutine ocn_vmix_coefs_cvmix2_kokkos_build!}}}

   subroutine ocn_vmix_cvmix2_kokkos_finalize()!{{{

           call f_kokkos_finalize()

   end subroutine ocn_vmix_cvmix2_kokkos_finalize!}}}

!***********************************************************************
!
!  routine ocn_vmix_cvmix2_kokkos_init
!
!> \brief   Initializes ocean vertical mixing quantities for cvmix 2.0
!> \author  Luke Van Roekel 
!> \date    19 February 2019
!> \details
!>  This routine initializes a variety of quantities related to
!>  vertical mixing in the ocean for cvmix2
!
!-----------------------------------------------------------------------


   subroutine ocn_vmix_cvmix2_kokkos_init(domain,err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      implicit none

      type (domain_type), intent(inout) :: domain !< Input/Output: domain information
      type (mpas_pool_type), pointer :: meshPool

      integer, intent(out) :: err !< Output: error flag

      integer :: nCells
      integer, dimension(:), pointer :: nCellsArray
      integer, pointer :: nVertLevels, config_cvmix_num_ri_smooth_loops
      type (block_type), pointer :: block


      cvmixOn = config_use_cvmix
      if (.not.config_use_cvmix) return

      ! init backround mixing
      backgroundDiff = 0.0_RKIND
      backgroundVisc = 0.0_RKIND

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
         call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

         nCells = nCellsArray( size(nCellsArray) )
         call f_kokkos_initialize(nVertLevels, nCells)
         block => block % next
      end do

      if(config_cvmix2_background_mixing_scheme == 'bryanLewis') then
      ! set Bryan and Lewis background profile specific values
         mixingOpts%bl1 = config_cvmix_BryanLewis_bl1
         mixingOpts%bl2 = config_cvmix_BryanLewis_bl2
         mixingOpts%bl3 = config_cvmix_BryanLewis_transitionDepth
         mixingOpts%bl4 = config_cvmix_BryanLewis_transitionWidth
         mixingOpts%bryanLewis = .true.
         mixingOpts%prandtl = config_cvmix_prandtl_number
      elseif(config_cvmix2_background_mixing_scheme == 'constant') then
         backgroundVisc = config_cvmix_background_viscosity
         backgroundDiff = config_cvmix_background_diffusion
      else
         call mpas_log_write("unknown option for config_cvmix2_background_mixing_scheme, allowed values are"  &
                   // " 'bryanLewis' or 'constant'", MPAS_LOG_CRIT)
      endif
      cvmixConvectionOn = config_use_cvmix_convection
      cvmixKPPOn = config_use_cvmix_kpp

      
      !
      ! When CVMix is turned on, all other vertical mixing schemes should be off
      ! Test to make sure this is the case.
      !
      ! test here, err=1 if a problem

      !
      ! pull nVertLevels out of the mesh structure
      !
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nVertLevels', nVertLevels)

      !
      ! initialize shear-based mixing
      !
      if (config_use_cvmix_shear) then
            mixingOpts%lShearKPP = .false.
            mixingOpts%lShearPP = .false.
            if( trim(config_cvmix_shear_mixing_scheme) == 'KPP' ) then
                mixingOpts%lShearKPP = .true.
                mixingOpts%KPP_nu_zero = config_cvmix_shear_KPP_nu_zero
                mixingOpts%KPP_Ri_zero = config_cvmix_shear_KPP_Ri_zero
                mixingOpts%KPP_exp = config_cvmix_shear_KPP_exp
            elseif ( trim(config_cvmix_shear_mixing_scheme) == 'PP' ) then
                mixingOpts%lShearPP = .true.
                mixingopts%PP_alpha = config_cvmix_shear_PP_alpha
                mixingOpts%PP_exp = config_cvmix_shear_PP_exp
                mixingOpts%PP_nu_zero = config_cvmix_shear_PP_nu_zero
                mixingOpts%PP_nu_b = config_cvmix_background_viscosity 
                mixingOpts%PP_kappa_b = config_cvmix_background_diffusion
            else 
                call mpas_log_write("unknown option for config_cvmix2_shear_mixing_scheme, allowed values are"  &
                        // " 'KPP' or 'PP'", MPAS_LOG_CRIT)
            endif

            mixingOpts%numRiSmoothPasses = config_cvmix_num_ri_smooth_loops
      endif

      !
      ! initialize convective mixing
      !
      if (config_use_cvmix_convection) then

        ! config_cvmix_convective_basedOnBVF is not supported at this time
        if (.not.config_cvmix_convective_basedOnBVF) then
            call mpas_log_write("config_cvmix_convective_basedOnBVF = .false. is not supported. Change to true.", MPAS_LOG_CRIT)
            err = 1
            return
        endif

        mixingOpts%lBruntVaisala = config_use_cvmix_convection
        mixingOpts%BVsqr_convect = config_cvmix_convective_triggerBVF
        mixingOpts%convect_diff = config_cvmix_convective_diffusion
        mixingOpts%convect_visc = config_cvmix_convective_viscosity

      endif

      !
      ! FUTURE: initialize tidal mixing
      !

      !
      ! FUTURE: initialize double diffusion
      !

      !
      ! FUTURE: initialize KPP boundary layer scheme
      !


      if (config_use_cvmix_kpp) then
         mixingOpts%Ri_crit = config_cvmix_kpp_criticalBulkRichardsonNumber
         mixingOpts%lEkman = config_cvmix_kpp_EkmanOBL
         mixingOpts%lMonOb = config_cvmix_kpp_MonObOBL
         mixingOpts%lKPP = config_use_cvmix_kpp
         mixingOpts%surf_layer_ext = config_cvmix_kpp_surface_layer_extent
         mixingOpts%c_s = config_cvmix2_cs_value
         mixingOpts%minVtsqr = config_cvmix2_kpp_minimumVt2 
         mixingOpts%rho_o = config_density0
         mixingOpts%minOSBLUnderIce = configure_cvmix_kpp_minimum_OBL_under_sea_ice
         mixingOpts%a_s = config_cvmix2_as_value
         mixingOpts%c_m = config_cvmix2_cm_value 
         mixingOpts%a_m = config_cvmix2_am_value
         mixingOpts%non_local_coeff = config_cvmix2_nonLocalTermStrength
         mixingOpts%parabolicNonLocal = config_cvmix2_useParabolicNonlocal

        !langmuir stuff
         mixingOpts%computeLangmuir = config_cvmix2_computeLangmuir
         if(mixingOpts%computeLangmuir) then
            mixingOpts%LANGMUIR_LWF16 = .true.     
            if( config_cvmix2_addLFK17_langmuir ) then
               mixingOpts%LANGMUIR_lifk17 = config_cvmix2_addLFK17_langmuir
            endif
         else
            mixingOpts%LANGMUIR_LWF16 = .false.
         endif

      endif

   !--------------------------------------------------------------------

   end subroutine ocn_vmix_cvmix2_kokkos_init!}}}

!***********************************************************************

end module ocn_vmix_cvmix2_kokkos

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
