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

   use ocn_tracer_hmix_del2
   use ocn_constants
   use ocn_config

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

   logical :: del2On
   real (kind=RKIND) :: eddyDiff2


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

   subroutine ocn_overflow_tend(meshPool, layerThickEdgeMean, tracers, tend, err)
      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickEdgeMean !< Input: mean thickness at edges

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
        tracers !< Input: tracer quantities

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
      integer :: i, k, iTracer, num_tracers, nCells
      integer, dimension(:), pointer :: nCellsArray

      integer, dimension(:), pointer :: minLevelEdgeBot, maxLevelEdgeTop, nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge, edgesOnCell, edgeSignOnCell

      real (kind=RKIND) :: invAreaCell
      real (kind=RKIND) :: tracer_turb_flux, flux, r_tmp

      real (kind=RKIND), dimension(:), pointer :: areaCell, dvEdge, dcEdge
      real (kind=RKIND), dimension(:), pointer :: meshScalingDel2

      err = 0

      if (.not.overflowOn) return

      call mpas_timer_start("tracer del2")

      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      num_tracers = size(tracers, dim=1)

      call mpas_pool_get_array(meshPool, 'minLevelEdgeBot', minLevelEdgeBot)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'meshScalingDel2', meshScalingDel2)

      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)

      nCells = nCellsArray( 1 )

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

          r_tmp = meshScalingDel2(iEdge) * eddyDiff2 * dvEdge(iEdge) / dcEdge(iEdge)

          tmp1 = (density(maxLevelCell(cell2),cell2) - density(maxLevelCell(cell1),cell1)) / dcEdge(iEdge) * &
                        (bottomDepth(cell2) - bottomDepth(cell1)) / dcEdge(iEdge)
          tmp2 = edgeSignOnCell(i, iCell)*velBot(1,iEdge)*(bottomDepth(cell2) - bottomDepth(cell1)) / dcEdge(iEdge)
          if(tmp1 < 0 and tmp2 > 0) then
             alphaH(iEdge) = 1.0_RKIND
             alphaZ(iCell) = 1.0_RKIND
             eddyDiff2 = dcEdge(iEdge)*abs(velBot(1,iEdge))
             eddyDiff2 = max(overflowDiffMin,min(overflowDiffMax,eddyDiff2))
          else
             alphaH(iEdge) = 0.0_RKIND
             eddyDiff2 = overflowDiffMin
          end if
          r_tmp = meshScalingDel2(iEdge) * eddyDiff2 * dvEdge(iEdge) / dcEdge(iEdge)

          do iTracer = 1, num_tracers
             ! \kappa_2 \nabla \phi on edge
             tracer_turb_flux = tracers(iTracer, maxLevelCell(cell2), cell2) - &
                                       tracers(iTracer, maxLevelCell(cell1), cell1)

             ! div(h \kappa_2 \nabla \phi) at cell center
             flux = layerThickEdgeBot(1, iEdge) * tracer_turb_flux * r_tmp

             tend(iTracer, maxLevelCell(iCell), iCell) = tend(iTracer, maxLevelCell(iCell), iCell) -  &
                                                         edgeSignOnCell(i, iCell) * flux * invAreaCell
            end do

        end do
      end do
      !$omp end do
      !$omp end parallel

! upstream advection next

      if (useOverflowAdvection) then
        ! need to compute vertical velocity from the sigma layer

        do iCell = 1, nCells
           invAreaCell = 1.0_RKIND / areaCell(icell)
           do i = 1, nEdgesOnCell(iCell)
              iEdge = edgesOnCell(i, iCell)
              cell1 = cellsOnEdge(1, iEdge)
              cell2 = cellsOnEdge(2, iEdge)
              layerThickEdgeBot = 0.5_RKIND*(layerThickness(maxLevelCell(cell1),cell1) + &
                                             layerThickness(maxLevelCell(cell2),cell2))
              normalThicknessFlux = layerThickEdgeBot*normalVelocity(maxLevelEdgeBot(iEdge),iEdge)

              tend(iTracer, maxLevelCell(iCell), iCell) = tend(iTracer, maxLevelCell(iCell), iCell) -      &
                                    alphaH(iEdge)*edgeSignOnCell(i, iCell)*invAreaCell*dvEdge(iEdge)*      &
                                    (max(0.0_RKIND,normalThicknessFlux)*tracers(maxLevelCell(cell1),cell1) &
                                  +  min(0.0_RKIND,normalThicknessFlux)*tracers(maxLevelCell(cell2),cell2))

          end do
          ! need vertical advection too
          tend(iTracer, maxLevelCell(iCell), iCell) = tend(iTracer
       end do

       ! vertical advection
   !--------------------------------------------------------------------

   end subroutine ocn_tracer_hmix_del2_tend!}}}

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

   subroutine ocn_tracer_hmix_del2_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      err = 0

      del2on = .false.

      if ( config_use_tracer_del2 ) then
         if ( config_tracer_del2 > 0.0_RKIND ) then
            del2On = .true.
            eddyDiff2 = config_tracer_del2
         endif
      endif

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_hmix_del2_init!}}}

!***********************************************************************

end module ocn_tracer_hmix_del2

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
