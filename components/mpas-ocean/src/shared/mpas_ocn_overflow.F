! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_overflow
!
!> \brief MPAS ocean overflow parameterization 
!> \author Luke Van Roekel, Steven Brus, Mat Maltrud 
!> \date  June 2024 
!> \details
!>  This module contains the main driver routine for computing
!>  the impacts of dense overflows in z-coordinate configurations at 
!>  low resolution.  The scheme is based on Beckmann and Doscher 1997 and 
!>  https://doi.org/10.1175/1520-0485(1997)027%3C0581:AMFIRO%3E2.0.CO;2 
!>  Doscher and Beckmann (2000) https://doi.org/10.1175/1520-0426(2000)017%3C0698:EOABBL%3E2.0.CO;2
!>  tendencies due to advection and diffusion in a bottom sigma layer

!>  It provides an init and a tend function. Each are described below.
!
!-----------------------------------------------------------------------

Module ocn_overflow

   use mpas_timer
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_threading

   use ocn_constants
   use ocn_config
   use ocn_mesh
   use ocn_diagnostics_variables

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

   public :: ocn_overflow_tend, &
             ocn_overflow_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: overflowOn, useOverflowAdvection
   real (kind=RKIND) :: eddyDiff2, overflowDiffMin, overflowDiffMax

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_overflow_tend
!
!> \brief   Computes tendency of dense overflows on tracer 
!> \author  Luke Van Roekel 
!> \date    June 2024 
!> \details
!>  This routine computes the impact of dense overflows on tracers via 
!>  advection and diffusion of tracers in a bottom sigma layer 
!
!-----------------------------------------------------------------------

   subroutine ocn_overflow_tend(normalVelocity, layerThickness, tracers, tend, err)
      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
        tracers !< Input: tracer quantities

      real (kind=RKIND), dimension(:,:), intent(in) :: &
        normalVelocity,  &
        layerThickness

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
         tend          !< Input/Output: velocity tendency

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

      integer :: iCell, iEdge, cell1, cell2
      integer :: nEdges, i, k, iTracer, num_tracers, nCells

      real (kind=RKIND) :: tmp1, tmp2, normalThicknessFlux,invAreaCell
      real (kind=RKIND) :: tracer_turb_flux, flux, r_tmp, ltEdge, div_hu, vertBot
      real (kind=RKIND),dimension(:), allocatable :: tempTend, eddyDiff
      err = 0

      if (.not.overflowOn) return

      call mpas_timer_start("dense overflow")

      nEdges = nEdgesHalo(1)
      allocate(eddyDiff(nEdges))
      do iEdge = 1,nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         if (maxLevelEdgeTop(iEdge) > 0) then
             tmp1 = (density(maxLevelCell(cell2),cell2) - density(maxLevelCell(cell1),cell1)) / dcEdge(iEdge) * &
                           (bottomDepth(cell2) - bottomDepth(cell1)) / dcEdge(iEdge)
             tmp2 = normalVelocity(maxLevelEdgeTop(iEdge),iEdge)*  &
                           (bottomDepth(cell2) - bottomDepth(cell1)) / dcEdge(iEdge)
             alphaH2(iEdge) = tmp1
             if(tmp1 < -3e-6) then
                alphaH(iEdge) = 1.00_RKIND
                eddyDiff(iEdge) = dcEdge(iEdge)*abs(normalVelocity(maxLevelEdgeTop(iEdge),iEdge))
                eddyDiff(iEdge) = max(overflowDiffMin,min(overflowDiffMax,eddyDiff2))
             else
                alphaH(iEdge) = 0.0_RKIND
                eddyDiff(iEdge) = overflowDiffMin
             end if
         end if
      end do

      num_tracers = size(tracers, dim=1)

      nCells = nCellsOwned

      !alphaH, eddyDiff gets computed in a separate loop to not compute each edge twice
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp private(invAreaCell, i, iEdge, cell1, cell2, r_tmp, k, iTracer, &
      !$omp         tracer_turb_flux, flux)
      do iCell = 1, nCells
        alphaZ(iCell) = 0.0_RKIND ! only set to 1 if alphaH is 1 on one of the edges
        invAreaCell = 1.0_RKIND / areaCell(iCell)
        do i = 1, nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i, iCell)
          cell1 = cellsOnEdge(1,iEdge)
          cell2 = cellsOnEdge(2,iEdge)

          if (maxLevelEdgeTop(iEdge) > 0) then
             if (alphaH(iEdge) > 0.0_RKIND) then
                 alphaZ(iCell) = 1.0_RKIND
             end if
             r_tmp = eddyDiff(iEdge) * dvEdge(iEdge) / dcEdge(iEdge)

             do iTracer = 1, num_tracers
                ! \kappa_2 \nabla \phi on edge
                tracer_turb_flux = tracers(iTracer, maxLevelCell(cell2), cell2) - &
                                          tracers(iTracer, maxLevelCell(cell1), cell1)

                ! div(h \kappa_2 \nabla \phi) at cell center
                flux = 0.5_RKIND*(layerThickness(maxLevelCell(cell1),cell1) + &
                                  layerThickness(maxLevelCell(cell2),cell2)) &
                                  * tracer_turb_flux * r_tmp

                tend(iTracer, maxLevelCell(iCell), iCell) = tend(iTracer, maxLevelCell(iCell), iCell) -  &
                                                            edgeSignOnCell(i, iCell) * flux * invAreaCell
               end do
           end if
        end do
      end do
      !$omp end do
      !$omp end parallel

! upstream advection next

      if (useOverflowAdvection) then

        ! need to compute vertical velocity from the sigma layer

        allocate(tempTend(num_tracers))
        alphaH(:) = 1.0_RKIND
        alphaZ(:) = 1.0_RKIND
        do iCell = 1, nCells
          tempTend(:) = 0.0_RKIND
           invAreaCell = 1.0_RKIND / areaCell(icell)
           div_hu = 0.0_RKIND
           do i = 1, nEdgesOnCell(iCell)
              iEdge = edgesOnCell(i, iCell)
              cell1 = cellsOnEdge(1, iEdge)
              cell2 = cellsOnEdge(2, iEdge)
              if (maxLevelEdgeTop(iEdge) > 0) then
                 ltEdge =min(layerThickness(maxLevelCell(cell1),cell1), &
                                                layerThickness(maxLevelCell(cell2),cell2))
                 normalThicknessFlux = ltEdge*normalVelocity(maxLevelEdgeTop(iEdge),iEdge)

                 div_hu = div_hu - ltEdge*edgeSignOnCell(i, iCell)*dvEdge(iEdge)* &
                             alphaH(iEdge)*normalVelocity(maxLevelEdgeTop(iEdge),iEdge)*invAreaCell

                 do iTracer = 1, num_tracers
                    tempTend(iTracer) = tempTend(iTracer) +      &
                                    alphaH(iEdge)*edgeSignOnCell(i, iCell)*invAreaCell*dvEdge(iEdge)*      &
                                    (max(0.0_RKIND,normalThicknessFlux)*tracers(iTracer, maxLevelCell(cell1),cell1) &
                                  +  min(0.0_RKIND,normalThicknessFlux)*tracers(iTracer, maxLevelCell(cell2),cell2))
                 end do
            end if
          end do

          if (maxLevelCell(iCell) > 10000) then
             vertBot = div_hu
             do iTracer = 1, num_tracers
                tempTend(iTracer) = tempTend(iTracer) + &
                                              alphaZ(iCell)*(min(0.0_RKIND,vertBot)*tracers(iTracer,maxLevelCell(iCell)-1,iCell) + &
                                               max(0.0_RKIND,vertBot)*tracers(iTracer,maxLevelCell(iCell),iCell))
             end do
          end if
          do iTracer = 1, num_tracers
             tend(iTracer,maxLevelCell(iCell),iCell) = tend(iTracer,maxLevelCell(iCell),iCell) - &
                                                       tempTend(iTracer)
          end do
     !          alphaZ(iCell) = 0.0_RKIND
       end do
       alphaH(:) = 0.0_RKIND
       alphaZ(:) = 0.0_RKIND
       deallocate(tempTend)
      else
       alphaH(:) = 0.0_RKIND
       alphaZ(:) = 0.0_RKIND
      end if
      call mpas_timer_stop("dense overflow") 
   !--------------------------------------------------------------------

   end subroutine ocn_overflow_tend!}}}

!***********************************************************************
!
!  routine ocn_tracer_hmix_del2_init
!
!> \brief   Initializes ocean tracer horizontal mixing quantities
!> \author  Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date    September 2011
!> \details
!>  This routine initializes a variety of quantities related to
!>  Laplacian horizontal velocity mixing in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_overflow_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      err = 0

      overflowOn = .false.

      if ( config_use_overflow ) then
         overflowOn = .true.
         useOverflowAdvection = config_use_overflow_advection
         overflowDiffMin = config_overflow_eddy_diffusivity_minimum
         overflowDiffMax = config_overflow_eddy_diffusivity_maximum
      endif

   !--------------------------------------------------------------------

   end subroutine ocn_overflow_init!}}}

!***********************************************************************

end module ocn_overflow

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
