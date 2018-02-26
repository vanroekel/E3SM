module elm_interface_dataType

!=================================================================================================
! ALM Thermal(T)-Hydrology (H) & BioGeoChemistry (BGC) Interface: Data Type (Variables)
! created: 8/25/2015
! update: 9/16/2016, 2/2/2017, May-2017, June-2017
! update: 5/13/2019 (all are IMPLICITLY column-based coupling,
!        esp. after ELM v2 data-structure modification @3/12/2019, commit 1bf22e32d)
!=================================================================================================
  ! USES:
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use shr_kind_mod          , only : r8 => shr_kind_r8
  use shr_infnan_mod        , only : nan => shr_infnan_nan, assignment(=)

  use elm_interface_thType
  use elm_interface_bgcType

  implicit none

  private

!-------------------------------------------------------------------------------------------------
  type, public :: elm_interface_data_type

     ! coupling level (grid/column/patch):
     integer :: cpl_level = 2                                ! coupling level: 1 - grid, 2 - column (default), 3 - patch

     ! col dimension
     real(r8), pointer :: z                                     (:,:) =>null()  ! layer depth (m) (-nlevsno+1:nlevgrnd)
     real(r8), pointer :: zi                                    (:,:) =>null()  ! layer depth (m) (-nlevsno+1:nlevgrnd)
     real(r8), pointer :: dz                                    (:,:) =>null()  ! layer thickness (m)  (-nlevsno+1:nlevgrnd)
     integer,  pointer :: snl                                   (:)   =>null()  ! number of snow layers (negative)

     ! soilstate_vars:
     real(r8), pointer :: bd                                    (:,:) =>null()  ! bulk density of dry soil material [kg/m^3] (CN)
     real(r8), pointer :: hksat                                 (:,:) =>null()  ! hydraulic conductivity at saturation (mm H2O /s)
     real(r8), pointer :: bsw                                   (:,:) =>null()  ! Clapp and Hornberger "b" (nlevgrnd)
     real(r8), pointer :: watsat                                (:,:) =>null()  ! volumetric soil water at saturation (porosity)
     real(r8), pointer :: watmin                                (:,:) =>null()  ! minimum volumetric soil water (nlevsoi)
     real(r8), pointer :: sucsat                                (:,:) =>null()  ! minimum soil suction (mm) (nlevgrnd)
     real(r8), pointer :: sucmin                                (:,:) =>null()  ! minimum allowable soil liquid suction pressure (mm) [Note: sucmin is a negative value, while sucsat is a positive quantity]
     real(r8), pointer :: watfc                                 (:,:) =>null()  ! volumetric soil water at field capacity (nlevsoi)
     real(r8), pointer :: porosity                              (:,:) =>null()  ! soil porisity (1-bulk_density/soil_density) (VIC)=>null()
     real(r8), pointer :: eff_porosity                          (:,:) =>null()  ! effective porosity = porosity - vol_ice (nlevgrnd)
     real(r8), pointer :: cellorg                               (:,:) =>null()  ! organic matter for gridcell containing column (1:nlevsoi)
     real(r8), pointer :: rootfr                                (:,:) =>null()  ! fraction of roots in each soil layer (nlevgrnd)

     real(r8), pointer :: tkfrz                                 (:,:) =>null()  ! thermal conductivity, frozen soil  [W/m-K] (new) (nlevgrnd)
     real(r8), pointer :: tkdry                                 (:,:) =>null()  ! thermal conductivity, dry soil (W/m/Kelvin) (nlevgrnd)
     real(r8), pointer :: tkwet                                 (:,:) =>null()  ! thermal conductivity, saturated soil [W/m-K] (new) (nlevgrnd)
     real(r8), pointer :: csol                                  (:,:) =>null()  ! heat capacity, soil solids (J/m**3/Kelvin) (nlevgrnd)

     ! thermal-hydrology:
     type(elm_interface_th_datatype) :: th

     ! biogeochemistry:
     type(elm_interface_bgc_datatype):: bgc

     !
  contains
     procedure , public  :: Init
     procedure , private :: InitAllocate
  end type elm_interface_data_type
!-------------------------------------------------------------------------------------------------
  

contains

!-------------------------------------------------------------------------------------------------
  subroutine Init(this, bounds)
     use decompMod               , only : bounds_type
     class(elm_interface_data_type) :: this
     type(bounds_type), intent(in)  :: bounds

     call this%InitAllocate (bounds)

     call this%th%Init(bounds)

     call this%bgc%Init(bounds)

  end subroutine Init
!-------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------

  subroutine InitAllocate(this, bounds)
    ! USES
    use elm_varpar            , only : nlevsno, nlevgrnd
    use elm_varcon            , only : spval
    use decompMod             , only : bounds_type

    ! ARGUMENTS:
    real(r8) :: ival  = 0.0_r8  ! initial value
    class(elm_interface_data_type)      :: this
    type(bounds_type), intent(in)       :: bounds

    ! LOCAL VARIABLES:
    integer  :: begc, endc
    !------------------------------------------------------------------------
    begc = bounds%begc; endc= bounds%endc
    if (this%cpl_level==1) then
      begc = bounds%begg; endc= bounds%endg
    elseif (this%cpl_level==3) then
      begc = bounds%begp; endc= bounds%endp
    endif

    ! col:
    allocate(this%z                     (begc:endc,-nlevsno+1:nlevgrnd))      ; this%z                    (:,:) = spval
    allocate(this%zi                    (begc:endc,-nlevsno+0:nlevgrnd))      ; this%zi                   (:,:) = spval
    allocate(this%dz                    (begc:endc,-nlevsno+1:nlevgrnd))      ; this%dz                   (:,:) = spval
    allocate(this%snl                   (begc:endc))                          ; this%snl                  (:)   = 0

    ! soilstate_vars:
    allocate(this%bd                (begc:endc,1:nlevgrnd))               ; this%bd               (:,:) = spval
    allocate(this%hksat             (begc:endc,1:nlevgrnd))               ; this%hksat            (:,:) = spval
    allocate(this%bsw               (begc:endc,1:nlevgrnd))               ; this%bsw              (:,:) = spval
    allocate(this%watsat            (begc:endc,1:nlevgrnd))               ; this%watsat           (:,:) = spval
    allocate(this%watmin            (begc:endc,1:nlevgrnd))               ; this%watmin           (:,:) = spval
    allocate(this%sucsat            (begc:endc,1:nlevgrnd))               ; this%sucsat           (:,:) = spval
    allocate(this%sucmin            (begc:endc,1:nlevgrnd))               ; this%sucmin           (:,:) = spval
    allocate(this%watfc             (begc:endc,1:nlevgrnd))               ; this%watfc            (:,:) = spval
    allocate(this%porosity          (begc:endc,1:nlevgrnd))               ; this%porosity         (:,:) = spval
    allocate(this%eff_porosity      (begc:endc,1:nlevgrnd))               ; this%eff_porosity     (:,:) = spval
    allocate(this%cellorg           (begc:endc,1:nlevgrnd))               ; this%cellorg          (:,:) = spval
    allocate(this%rootfr            (begc:endc,1:nlevgrnd))               ; this%rootfr           (:,:) = spval

    allocate(this%tkwet             (begc:endc,1:nlevgrnd))               ; this%tkwet            (:,:) = spval
    allocate(this%tkdry             (begc:endc,1:nlevgrnd))               ; this%tkdry            (:,:) = spval
    allocate(this%tkfrz             (begc:endc,1:nlevgrnd))               ; this%tkfrz            (:,:) = spval
    allocate(this%csol              (begc:endc,1:nlevgrnd))               ; this%csol             (:,:) = spval

  end subroutine InitAllocate
!-------------------------------------------------------------------------------------------------

end module elm_interface_dataType
