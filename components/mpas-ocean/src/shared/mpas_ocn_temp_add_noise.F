! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_noise
!
!> \brief adds noise to temperature 
!> \author Luke Van Roekel
!> \date   20 March 2023
!> \details
!>  This module contains the routine for adding white  
!>  noise to LES configurations
!
!-----------------------------------------------------------------------

module ocn_temp_noise

   use mpas_timer
   use mpas_derived_types
   use mpas_pool_routines
   use ocn_constants
   use ocn_config
   use ocn_mesh

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

   public :: ocn_temperature_add_noise, &
             ocn_temperature_noise_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: noiseOn

   real(kind=RKIND) :: maxTime,   & !time to end the noise
                       curTime

   integer :: minNoiseLevel, maxNoiseLevel

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_temp_add_noise
!
!> \brief   Computes a white noise perturbation on temperature 
!> \author  Luke Van Roekel
!> \date    7 June 2022
!> \details
!>  This routine computes the a white noise tendency to apply to temperature
!
!-----------------------------------------------------------------------

   subroutine ocn_temperature_add_noise(tracer, index_temperature, dt, err)
      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real(kind=RKIND),intent(in) :: dt

      integer,intent(in) :: index_temperature

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
         tracer    !< Input/Output: velocity tendency

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: Error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iCell, k, nCells

      real(kind=RKIND) :: randomPerturbation, finalPert

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing tendencies
      ! note that the user can choose multiple options and the
      !   tendencies will be added together
      !
      !-----------------------------------------------------------------

      err = 0

      if(.not.noiseOn) return


      if(curTime > maxTime) return
      curTime = curTime + dt

      call mpas_timer_start("noise tend")

      nCells = nCellsAll

      !$omp parallel
      !$omp do schedule(runtime) private(k,randomPerturbation,finalPert,minNoiseLevel,maxNoiseLevel)
      do iCell = 1, nCells
         minNoiseLevel = max(minLevelCell(iCell),config_LES_noise_min_level)
         maxNoiseLevel = min(maxLevelCell(iCell),config_LES_noise_max_level)

         do k = minNoiseLevel, maxNoiseLevel
           call random_number(randomPerturbation)
           randomPerturbation = 2.0_RKIND * randomPerturbation - 1.0_RKIND
           finalPert = config_LES_temp_noise_perturbation_magnitude * randomPerturbation
           tracer(index_temperature,k,iCell) = tracer(index_temperature,k,iCell) + dt*finalPert 
         end do
      end do
      !$omp end do
      !$omp end parallel

      call mpas_timer_stop("noise tend")

   !--------------------------------------------------------------------

   end subroutine ocn_temperature_add_noise

!***********************************************************************
!
!  routine ocn_temperature_noise_init
!
!> \brief   Initializes ocean temperature noise tendency 
!> \author  Luke Van Roekel
!> \date    7 June 2022
!> \details
!>  This routine initializes quantities related to white noise
!>  tendency of temperature
!
!-----------------------------------------------------------------------

   subroutine ocn_temperature_noise_init(err)

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      
      noiseOn = .false.

      if(config_LES_temp_noise_enable) then
         noiseOn = .true.
         maxTime = config_LES_temp_noise_max_time
         curTime = 0.0_RKIND

      end if

      err = 0

   !--------------------------------------------------------------------

   end subroutine ocn_temperature_noise_init

!***********************************************************************

end module ocn_temp_noise

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
