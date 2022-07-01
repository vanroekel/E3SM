! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vvel_vert_advection
!
!> \brief MPAS vertical advection for (non-hydrostatic) vertical velocity
!> \author Luke Van Roekel, Darren Engwirda
!> \date   October 2021
!> \details
!>  This module contains routines for advection of (non-hydrostatic) vertical 
!>  velocity.
!
!-------------------------------------------------------------------------------

module ocn_vvel_vert_advection

   ! module includes
   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_io_units
   use mpas_threading
   use ocn_mesh
   use mpas_tracer_advection_helpers
   use ocn_tracer_advection_shared

   implicit none
   private
   save

   ! public method interfaces
   public :: ocn_vvel_vert_advection_tend, &
             ocn_vvel_vert_advection_init

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_vvel_vert_advection_tend
!
!> \brief MPAS vertical advection tendency for (non-hydrostatic) vertical velocity
!> \author Luke Van Roekel, Darren Engwirda
!> \date   27/10/21
!> \details
!>  NOTE this is only active for non-hydrostatic mode
!
!-----------------------------------------------------------------------
   subroutine ocn_vvel_vert_advection_tend(verticalVelocity, &!{{{
                                             layerThickness, dt, &
                                             vertAleTransportTop, tend)

      implicit none

      real (kind=RKIND), dimension(:,:), intent(in) :: verticalVelocity !< Input: current vertVel values
      real (kind=RKIND), dimension(:,:), intent(in) :: layerThickness !< Input: Thickness
      real (kind=RKIND), intent(in) :: dt !< Input: Timestep
      real (kind=RKIND), dimension(:,:), intent(inout) :: tend !< Input/Output: VertVel tendency
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         vertAleTransportTop  !< [in] velocity through layer top levels
      integer :: i, iCell, iEdge, k, cell1, cell2
      integer :: nCells, nVertLevels
      integer :: nEdges, maxEdges
      real (kind=RKIND), parameter :: eps = 1.e-10_RKIND
      real (kind=RKIND) :: levelThickness, vertVelocityUp, vertVelocityDn


      nCells = nCellsAll
      ! -w_ale * dw/dz
      !$omp parallel
      !$omp do schedule(runtime)
      !$omp private(k, levelThickness, vertVelocityUp, vertVelocityDn)
      do iCell = 1, nCells
         !w_ale = 0.0 at k==1
         do k = 2, maxLevelCell(iCell)
            levelThickness = 0.5_RKIND * (layerThickness(k-1,iCell) + &
                                          layerThickness(k  ,iCell))

            vertVelocityUp = 0.5_RKIND * (verticalVelocity(k-1,iCell) + &
                                          verticalVelocity(k  ,iCell))
            vertVelocityDn = 0.5_RKIND * (verticalVelocity(k+1,iCell) + &
                                          verticalVelocity(k  ,iCell))
            
            tend(k,iCell) = tend(k,iCell) - vertAleTransportTop(k,iCell) * &
                                (vertVelocityUp - vertVelocityDn) / levelThickness
         end do
      end do
      !$omp end do
      !$omp end parallel


end subroutine ocn_vvel_vert_advection_tend!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_vvel_vert_advection_init
!
!> \brief MPAS initialize vertical velocity advection tendency.
!> \author Luke Van Roekel, Darren Engwirda
!> \date   27/10/2021
!> \details
!>  This routine initializes constants and choices for the (non-hydrostatic) 
!>  vertical velocity advection tendency
!
!-------------------------------------------------------------------------------

   subroutine ocn_vvel_vert_advection_init(err) !{{{

      implicit none

      !*** input parameters

      !*** output parameters

      integer, intent(out) :: err !< [out] Error Flag

      err = 0 ! set error code to success

      ! don't have actual choices for vvel advection here, unlike tracers   

   end subroutine ocn_vvel_vert_advection_init!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

end module ocn_vvel_vert_advection

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
