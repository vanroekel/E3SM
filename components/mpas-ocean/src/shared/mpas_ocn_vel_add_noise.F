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

   logical :: first, noiseOn

   real(kind=RKIND) :: maxTime,   & !time to end the noise
                       curTime

   integer :: minNoiseLevel, maxNoiseLevel

   real(kind=RKIND),dimension(:,:),allocatable :: smooth1, smooth2 

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
      !----------------------------------------------------------------

      integer, intent(out) :: err !< Output: Error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iEdge, k, nEdges, iCell, nCells, c1, c2, iNeighbor, edgeCount

      real(kind=RKIND) :: quant, cmask, areaSum
      real(kind=RKIND) :: uav,vav,randomPerturbation(2), finalPert(2), cosAng, sinAng
      real(kind=RKIND),dimension(:,:),allocatable :: unoise, vnoise
      real(kind=RKIND),dimension(nVertLevels) :: noise1, noise2

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

      !steps for a v2 noise
      ! create a noise on cells and only allow 2 to nvertLevels-1
      ! smooth in vertical with a 1-2-1 filter
      ! interpolate from centers to edges

      nEdges = nEdgesHalo(1)
      nCells = nCellsAll
      if(first) then
        first = .false.
      allocate(unoise(nVertLevels,nCells),vnoise(nVertLevels,nCells))

      !$omp parallel
      !$omp do schedule(runtime) private(noise1,noise2,k,randomPerturbation,finalPert,minNoiseLevel,maxNoiseLevel)
      do iCell = 1, nCells
         minNoiseLevel = max(minLevelCell(iCell),config_LES_noise_min_level)
         maxNoiseLevel = min(maxLevelCell(iCell),config_LES_noise_max_level)

         unoise(:,iCell) = 0.0_RKIND
         vnoise(:,iCell) = 0.0_RKIND
!         cosAng = cos(angleEdge(iEdge))
!         sinAng = sin(angleEdge(iEdge))
         do k = minNoiseLevel, maxNoiseLevel
           call random_number(randomPerturbation)
           randomPerturbation = 2.0_RKIND * randomPerturbation - 1.0_RKIND
           finalPert = config_LES_noise_perturbation_magnitude * randomPerturbation
           unoise(k,iCell) = finalPert(1)
           vnoise(k,iCell) = finalPert(2)
!           normalVelocity(k,iEdge) = normalVelocity(k,iEdge) + cosAng*finalPert(1) + sinAng*finalPert(2)
         end do

         noise1(:) = unoise(:,iCell)
         noise2(:) = vnoise(:,iCell)
         do k = minNoiseLevel+1,maxNoiseLevel-1
            noise1(k) = 0.25_RKIND*(unoise(k-1,iCell) + 2.0_RKIND*unoise(k,iCell) + unoise(k+1,iCell))
            noise2(k) = 0.25_RKIND*(vnoise(k-1,iCell) + 2.0_RKIND*vnoise(k,iCell) + vnoise(k+1,iCell))
         end do

         unoise(:,iCell) = noise1(:)
         vnoise(:,iCell) = noise2(:)
      end do
      !$omp end do
      !$omp end parallel

      nCells = nCellsHalo(1)
      !next smooth in the horizontal with basic laplacian filter
      smooth1(:,:) = 0.0_RKIND
      smooth2(:,:) = 0.0_RKIND
      !$omp parallel
      !$omp do schedule(runtime) private(quant, cmask, k, nEdges, areaSum, edgeCount, iEdge, iNeighbor)
      do iCell=1,nCells
         nEdges = nEdgesOnCell(iCell)
         areaSum = 0.0_RKIND
         edgeCount = 0

         do iEdge = 1,nEdges
            iNeighbor = cellsOnCell(iEdge,iCell)
            cmask = cellMask(minLevelCell(iNeighbor),iNeighbor)
            quant = cmask*areaCell(iNeighbor)
            areaSum = areaSum + 2.0_RKIND*quant
            edgeCount = edgeCount + cmask

            do k = 1,maxLevelCell(iCell)
               smooth1(k,iCell) = smooth1(k,iCell) + 2.0_RKIND*quant*unoise(k,iNeighbor)
               smooth2(k,iCell) = smooth2(k,iCell) + 2.0_RKIND*quant*vnoise(k,iNeighbor)
            end do

            areaSum = areaSum + edgeCount + areaCell(iCell)
            do k=1,maxLevelCell(iCell)
               smooth1(k,iCell) = smooth1(k,iCell) + unoise(k,iCell)*edgeCount*areaCell(iCell)
               smooth1(k,iCell) = smooth1(k,iCell) / areaSum
               smooth2(k,iCell) = smooth2(k,iCell) + vnoise(k,iCell)*edgeCount*areaCell(iCell)
               smooth2(k,iCell) = smooth2(k,iCell) / areaSum
            end do
        end do

    end do
    !$omp end do
    deallocate(unoise,vnoise)
    end if


      ! Now interpolate to edges
      !$omp parallel
      !$omp do schedule(runtime) private(c1, c2, uav, vav, cosAng, sinAng)
      do iEdge = 1, nEdgesOwned
         c1 = cellsOnEdge(1,iEdge)
         c2 = cellsOnEdge(2,iEdge)
         cosAng = cos(angleEdge(iEdge))
         sinAng = sin(angleEdge(iEdge))
         do k = minLevelEdgeTop(iEdge),maxLevelEdgeTop(iEdge)
            uav = 0.5_RKIND*(smooth1(k,c1) + smooth1(k,c2))
            vav = 0.5_RKIND*(smooth2(k,c1) + smooth2(k,c2))
            normalVelocity(k,iEdge) = normalVelocity(k,iEdge) + dt*(uav*cosAng + vav*sinAng)
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

      allocate(smooth1(nVertLevels,nCellsAll),smooth2(nVertLevels,nCellsAll))

      if(config_LES_noise_enable) then
         first = .true.
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
