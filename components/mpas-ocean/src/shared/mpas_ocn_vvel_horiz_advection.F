! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vvel_horiz_advection
!
!> \brief MPAS horizontal advection for (non-hydrostatic) vertical velocity
!> \author Luke Van Roekel, Darren Engwirda
!> \date   October 2021
!> \details
!>  This module contains routines for horizontal advection of vertical velocity
!>  via a simple TRiSK-type scheme at vertical levels (ie. half-layers).
!
!-------------------------------------------------------------------------------

module ocn_vvel_horiz_advection

   ! module includes
   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_io_units
   use mpas_threading
   use ocn_mesh
   use ocn_config
   use mpas_tracer_advection_helpers
   use ocn_tracer_advection_shared
   use ocn_diagnostics_variables

   implicit none
   private
   save

   ! public method interfaces
   public :: ocn_vvel_horiz_advection_tend, &
             ocn_vvel_horiz_advection_init

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_vvel_horiz_advection_tend
!
!> \brief MPAS horizontal advection tendency for (non-hydrostatic) vertical velocity
!> \author Luke Van Roekel, Darren Engwirda
!> \date   27/10/21
!> \details
!>  NOTE this is only active for non-hydrostatic mode
!
!-----------------------------------------------------------------------
   subroutine ocn_vvel_horiz_advection_tend(verticalVelocity, normalThicknessFlux, &!{{{
                                             layerThickness, dt, tend)

      implicit none

      real (kind=RKIND), dimension(:,:), intent(in) :: normalThicknessFLux !< Input: Thichness weighted velocitiy
      real (kind=RKIND), dimension(:,:), intent(in) :: layerThickness !< Input: Thickness
      real (kind=RKIND), dimension(:,:), intent(in) :: verticalVelocity
      real (kind=RKIND), intent(in) :: dt !< Input: Timestep
      real (kind=RKIND), dimension(:,:), intent(inout) :: tend !< Input/Output: Vertical Velocity tendency

      integer :: i, iCell, jEdge, iEdge, k, cell1, cell2
      integer :: nVertLevels, nCells, nEdges

      real (kind=RKIND) :: vertVelocityFlux, vvel_weight, invAreaCell1

      real (kind=RKIND), dimension(:,:), allocatable :: horizFlux
      real (kind=RKIND), dimension(:), allocatable :: invVolCell

      real (kind=RKIND), parameter :: eps = 1.e-10_RKIND

      ! Get dimensions
     nVertLevels = size(verticalVelocity,dim=1)

      ! Initialize pointers
      nCells = nCellsHalo( 2 )
      nEdges = nEdgesHalo( 3 )
      allocate(horizFlux(nVertLevels, nEdges))
      allocate(invVolCell(nVertLevels))

      !$omp parallel
       !$omp do schedule(runtime)
        do iEdge = 1, nEdges
         horizFlux(:, iEdge) = 0.0_RKIND
        end do
        !$omp end do
      !$omp end parallel

      ! compute edge fluxes for div(u*h*w)
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp private(k, cell1, cell2, vertVelocityFlux)
      do iEdge = 1, nEdges
         cell1 = cellsOnEdge(1, iEdge)
         cell2 = cellsOnEdge(2, iEdge)
         !omp do schedule(runtime)
         do k = 1, maxLevelEdgeTop(iEdge)!
            if (k .gt. 1) then
               vertVelocityFlux = &
                  0.5_RKIND * (normalThicknessFlux(k-1,iEdge) + &
                               normalThicknessFlux(k-0,iEdge)) * &
                  0.5_RKIND * (verticalVelocity(k,cell1) + verticalVelocity(k,cell2))
            else
               vertVelocityFlux = & ! surface
                  1.0_RKIND * (normalThicknessFlux(k-0,iEdge)) * &
                  0.5_RKIND * (verticalVelocity(k,cell1) + verticalVelocity(k,cell2))
            endif

            horizFlux(k,iEdge) = &
               edgeMask(k,iEdge) * dvEdge(iEdge) * vertVelocityFlux
         end do
      end do
      !$omp end do
      !$omp end parallel

      ! eval. tendencies as div(u*h*w) / (A*h)
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp private(k, iEdge, jEdge, invVolCell)
        do iCell = 1, nCells
         invVolCell(:) =-1.0e34_RKIND
         invVolCell(1) = 1.0_RKIND / (areaCell(iCell) * layerThickness(1,iCell))
         do k = 2, maxLevelCell(iCell)
            invVolCell(k) = 2.0_RKIND / ( &
               areaCell(iCell) * (layerThickness(k-1,iCell) + layerThickness(k-0,iCell)))
         end do

         do jEdge = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(jEdge,iCell)
            do k = 1, maxLevelEdgeTop(iEdge)
               tend(k,iCell) = tend(k,iCell) + &
                  edgeSignOnCell(jEdge,iCell) * invVolCell(k) * horizFlux(k,iEdge) 
            end do
          end do
        end do
        !$omp end do
      !$omp end parallel

      deallocate(horizFlux)
      deallocate(invVolCell)

   end subroutine ocn_vvel_horiz_advection_tend!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_vvel_horiz_advection_init
!
!> \brief MPAS initialize vertical velocity advection tendency.
!> \author Luke Van Roekel, Darren Engwirda
!> \date   27/10/2021
!> \details
!>  This routine initializes constants and choices for the (non-hydrostatic) 
!>  vertical velocity advection tendency
!
!-------------------------------------------------------------------------------

   subroutine ocn_vvel_horiz_advection_init(err) !{{{

      implicit none

      !*** input parameters

      !*** output parameters

      integer, intent(out) :: err !< [out] Error Flag

      err = 0 ! set error code to success

      ! don't have actual choices for vvel advection here, unlike tracers   

   end subroutine ocn_vvel_horiz_advection_init!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

end module ocn_vvel_horiz_advection

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
