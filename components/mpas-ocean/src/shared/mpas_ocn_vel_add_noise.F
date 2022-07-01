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
!> \brief adds noise to normal velocity
!> \author Luke Van Roekel
!> \date   7 June 2022
!> \details
!>  This module contains the routine for adding white  
!>  noise to LES configurations
!
!-----------------------------------------------------------------------

module ocn_vel_noise

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

   public :: ocn_vel_add_noise, &
             ocn_vel_noise_init

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
!  routine ocn_vel_add_noise
!
!> \brief   Computes a white noise perturbation on normal velocity 
!> \author  Luke Van Roekel
!> \date    7 June 2022
!> \details
!>  This routine computes the a white noise tendency to apply to normal velocity 
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_add_noise(normalVelocity, dt, err)!{{{
      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real(kind=RKIND),intent(in) :: dt

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         normalVelocity    !< Input/Output: velocity tendency

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

      integer :: iEdge, k, nEdges

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

      nEdges = nEdgesOwned

      !$omp parallel
      !$omp do schedule(runtime) private(k,randomPerturbation,finalPert,minNoiseLevel,maxNoiseLevel)
      do iEdge = 1, nEdges
         call random_number(randomPerturbation)
         randomPerturbation = 2.0_RKIND * randomPerturbation - 1.0_RKIND
         finalPert = config_LES_noise_perturbation_magnitude * randomPerturbation
         minNoiseLevel = max(minLevelEdgeTop(iEdge),config_LES_noise_min_level)
         maxNoiseLevel = min(maxLevelEdgeTop(iEdge),config_LES_noise_max_level)

         do k = minNoiseLevel, maxNoiseLevel
            normalVelocity(k,iEdge) = normalVelocity(k,iEdge) + finalPert
         end do
      end do
      !$omp end do
      !$omp end parallel

      call mpas_timer_stop("noise tend")

   !--------------------------------------------------------------------

   end subroutine ocn_vel_add_noise!}}}

!***********************************************************************
!
!  routine ocn_vel_noise_init
!
!> \brief   Initializes ocean normal velocity noise tendency 
!> \author  Luke Van Roekel
!> \date    7 June 2022
!> \details
!>  This routine initializes quantities related to white noise
!>  tendency of normal velocity
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_noise_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      
      noiseOn = .false.

      if(config_LES_noise_enable) then
         noiseOn = .true.
         maxTime = config_LES_noise_max_time
         curTime = 0.0_RKIND
      end if

      err = 0

   !--------------------------------------------------------------------

   end subroutine ocn_vel_noise_init!}}}

!***********************************************************************

end module ocn_vel_noise

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
