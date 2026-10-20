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
!> \brief MPAS ocean bottom boundary layer (dense downslope flow) parameterization
!> \author Luke Van Roekel, Steven Brus, Mat Maltrud
!> \date   June 2024, revised October 2026
!> \details
!>  This module parameterizes the downslope flow of dense bottom water
!>  across topographic steps in z-level (partial bottom cell) configurations.
!>  It acts at every step edge, i.e. an edge whose two cells have different
!>  bottom levels. The cell with the shallower bottom level kS is the shallow
!>  cell, the other (bottom level kD > kS) the deep cell. The scheme is active
!>  where the bottom water of the shallow cell is denser than the bottom water
!>  of the deep cell, compared at a common pressure. Two independent pieces:
!>
!>  Diffusive (Beckmann and Doscher 1997, Doscher and Beckmann 2000)
!>    https://doi.org/10.1175/1520-0485(1997)027<0581:AMFIRO>2.0.CO;2
!>    https://doi.org/10.1175/1520-0426(2000)017<0698:EOABBL>2.0.CO;2
!>    a down-gradient flux between the two bottom cells, shallow(kS) and
!>    deep(kD), with diffusivity dcEdge*|u| bounded by a minimum and maximum.
!>
!>  Advective (Campin and Goosse 1999, as in the MITgcm down_slope package)
!>    https://doi.org/10.3402/tellusa.v51i3.13468
!>    a closed loop of volume transport Q carried by tracers only:
!>    shallow(kS) -> deep(kD), up the deep column from kD to kS, and back
!>    across the edge to shallow(kS). The loop moves no volume, so thickness,
!>    momentum and the vertical velocity are untouched, and every leg is
!>    upwind, so the scheme is conservative and monotone as long as Q is
!>    bounded by the volume of the cells on the loop.
!>
!>  Edge transports are computed once per tracer tendency call
!>  (ocn_overflow_transport) and applied to every tracer group
!>  (ocn_overflow_tend). CPU only: there are no OpenACC directives.
!
!-----------------------------------------------------------------------

module ocn_overflow

   use mpas_timer
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_log
   use mpas_constants

   use ocn_constants
   use ocn_config
   use ocn_mesh
   use ocn_diagnostics_variables

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_overflow_transport, &
             ocn_overflow_tend, &
             ocn_overflow_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: overflowOn, overflowDiffusionOn, overflowAdvectionOn, &
              overflowGravityVelocity

   real (kind=RKIND) :: overflowDiffMin, overflowDiffMax, &
                        overflowGravityCoeff, overflowMaxVolumeFraction, &
                        overflowDensityThreshold

   ! sound speed used to shift density between the small pressure
   ! difference of two cells at the same level (m/s)
   real (kind=RKIND), parameter :: soundSpeed = 1500.0_RKIND

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_overflow_transport
!
!> \brief   Computes the bottom boundary layer transports at step edges
!> \author  Luke Van Roekel
!> \date    October 2026
!> \details
!>  For each step edge, computes the density difference between the
!>  shallow and deep bottom water at the deep bottom pressure, the
!>  advective loop transport overflowAdvectiveTransport (m3/s) and the
!>  diffusive exchange overflowDiffusiveTransport (m3/s, the flux per
!>  unit tracer difference). Both are zero where the scheme is inactive.
!>
!>  Density difference: density at the same level kS in both cells
!>  (corrected for their small pressure difference), plus the stability
!>  of the deep column from kS down to kD, using displacedDensity(k)
!>  (layer k displaced to the depth of layer k+1):
!>    drho = rho_s(kS) - rho_d(kS)
!>         + sum_{k=kS}^{kD-1} [displacedDensity_d(k) - density_d(k+1)]
!>  drho > 0 means the shallow bottom water would sink below the deep
!>  bottom cell.
!
!-----------------------------------------------------------------------

   subroutine ocn_overflow_transport(layerThickness, dt, err)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness  !< [in] layer thickness

      real (kind=RKIND), intent(in) :: &
         dt              !< [in] time step (s), for the volume limit

      ! these variables input from shared diagnostics module:
      !    density, displacedDensity, pressure, normalTransportVelocity

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< [out] error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: iEdge, k, kS, kD, cellS, cellD, nEdges

      real (kind=RKIND) :: drho, uDown, slope, eddyDiff, hEdge, &
                           volMin, transportMax

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      err = 0

      if (.not. overflowOn) return

      call mpas_timer_start("overflow transport")

      nEdges = nEdgesHalo(1)

      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp private(k, kS, kD, cellS, cellD, drho, uDown, slope, &
      !$omp         eddyDiff, hEdge, volMin, transportMax)
      do iEdge = 1, nEdges
         overflowDensityDifference(iEdge) = 0.0_RKIND
         overflowAdvectiveTransport(iEdge) = 0.0_RKIND
         overflowDiffusiveTransport(iEdge) = 0.0_RKIND

         if (maxLevelCell(cellsOnEdge(1,iEdge)) < 1 .or. &
             maxLevelCell(cellsOnEdge(2,iEdge)) < 1) cycle
         if (maxLevelCell(cellsOnEdge(1,iEdge)) == &
             maxLevelCell(cellsOnEdge(2,iEdge))) cycle

         if (maxLevelCell(cellsOnEdge(1,iEdge)) < &
             maxLevelCell(cellsOnEdge(2,iEdge))) then
            cellS = cellsOnEdge(1,iEdge)
            cellD = cellsOnEdge(2,iEdge)
            ! normalTransportVelocity is positive from cell 1 to cell 2
            uDown = normalTransportVelocity(maxLevelCell(cellS),iEdge)
         else
            cellS = cellsOnEdge(2,iEdge)
            cellD = cellsOnEdge(1,iEdge)
            uDown = -normalTransportVelocity(maxLevelCell(cellS),iEdge)
         end if
         kS = maxLevelCell(cellS)
         kD = maxLevelCell(cellD)

         ! the edge must be open at kS and the deep column must reach up to kS
         if (maxLevelEdgeTop(iEdge) < kS .or. kS < minLevelCell(cellD)) cycle

         drho = density(kS,cellS) - density(kS,cellD) &
              + (pressure(kS,cellD) - pressure(kS,cellS)) / soundSpeed**2
         do k = kS, kD - 1
            drho = drho + displacedDensity(k,cellD) - density(k+1,cellD)
         end do
         overflowDensityDifference(iEdge) = drho

         if (drho <= overflowDensityThreshold) cycle

         ! no single leg may move more than a fraction of the smallest
         ! cell volume on the loop per time step
         volMin = areaCell(cellS) * layerThickness(kS,cellS)
         do k = kS, kD
            volMin = min(volMin, areaCell(cellD) * layerThickness(k,cellD))
         end do
         transportMax = overflowMaxVolumeFraction * volMin / dt

         if (overflowAdvectionOn) then
            if (overflowGravityVelocity) then
               slope = (bottomDepth(cellD) - bottomDepth(cellS)) / dcEdge(iEdge)
               uDown = overflowGravityCoeff * gravity * drho / rho_sw * slope
            end if
            if (uDown > 0.0_RKIND) then
               overflowAdvectiveTransport(iEdge) = min(transportMax, &
                  uDown * layerThickness(kS,cellS) * dvEdge(iEdge))
            end if
         end if

         if (overflowDiffusionOn) then
            eddyDiff = dcEdge(iEdge) * &
                       abs(normalTransportVelocity(kS,iEdge))
            eddyDiff = max(overflowDiffMin, min(overflowDiffMax, eddyDiff))
            hEdge = 0.5_RKIND * (layerThickness(kS,cellS) + &
                                 layerThickness(kD,cellD))
            overflowDiffusiveTransport(iEdge) = min(transportMax, &
               eddyDiff * hEdge * dvEdge(iEdge) / dcEdge(iEdge))
         end if
      end do
      !$omp end do
      !$omp end parallel

      call mpas_timer_stop("overflow transport")

   !--------------------------------------------------------------------

   end subroutine ocn_overflow_transport!}}}

!***********************************************************************
!
!  routine ocn_overflow_tend
!
!> \brief   Adds the bottom boundary layer tendency for a tracer group
!> \author  Luke Van Roekel
!> \date    October 2026
!> \details
!>  Applies the edge transports from ocn_overflow_transport to the
!>  thickness-weighted tracer tendency (tracer * m/s) of owned cells.
!>  Each cell gathers the contributions of its own edges, so there are
!>  no write conflicts between threads.
!>    advective loop, transport Q > 0 at an edge:
!>      shallow(kS):          Q * (T_d(kS)   - T_s(kS))
!>      deep(k), kS<=k<kD:    Q * (T_d(k+1)  - T_d(k))
!>      deep(kD):             Q * (T_s(kS)   - T_d(kD))
!>    diffusive exchange, transport Qd:
!>      shallow(kS):          Qd * (T_d(kD) - T_s(kS))
!>      deep(kD):             Qd * (T_s(kS) - T_d(kD))
!>  all divided by the cell area.
!
!-----------------------------------------------------------------------

   subroutine ocn_overflow_tend(tracers, tend, err)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         tracers  !< [in] tracers in the current group

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
         tend     !< [inout] tracer tendency for the group

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< [out] error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: iCell, iEdge, i, k, n, kS, kD, cellS, cellD, &
                 nCells, nTracers

      real (kind=RKIND) :: invAreaCell, qAdv, qDif

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      err = 0

      if (.not. overflowOn) return

      call mpas_timer_start("overflow tend")

      nCells = nCellsOwned
      nTracers = size(tracers, dim=1)

      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp private(invAreaCell, i, iEdge, k, n, kS, kD, cellS, cellD, &
      !$omp         qAdv, qDif)
      do iCell = 1, nCells
         invAreaCell = 1.0_RKIND / areaCell(iCell)
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i,iCell)
            qAdv = overflowAdvectiveTransport(iEdge) * invAreaCell
            qDif = overflowDiffusiveTransport(iEdge) * invAreaCell
            if (qAdv <= 0.0_RKIND .and. qDif <= 0.0_RKIND) cycle

            ! transports are only nonzero at step edges, where the
            ! bottom levels of the two cells differ
            if (maxLevelCell(cellsOnEdge(1,iEdge)) < &
                maxLevelCell(cellsOnEdge(2,iEdge))) then
               cellS = cellsOnEdge(1,iEdge)
               cellD = cellsOnEdge(2,iEdge)
            else
               cellS = cellsOnEdge(2,iEdge)
               cellD = cellsOnEdge(1,iEdge)
            end if
            kS = maxLevelCell(cellS)
            kD = maxLevelCell(cellD)

            if (iCell == cellS) then
               do n = 1, nTracers
                  tend(n,kS,iCell) = tend(n,kS,iCell) &
                     + qAdv * (tracers(n,kS,cellD) - tracers(n,kS,cellS)) &
                     + qDif * (tracers(n,kD,cellD) - tracers(n,kS,cellS))
               end do
            else
               do k = kS, kD - 1
                  do n = 1, nTracers
                     tend(n,k,iCell) = tend(n,k,iCell) &
                        + qAdv * (tracers(n,k+1,cellD) - tracers(n,k,cellD))
                  end do
               end do
               do n = 1, nTracers
                  tend(n,kD,iCell) = tend(n,kD,iCell) &
                     + (qAdv + qDif) * (tracers(n,kS,cellS) - tracers(n,kD,cellD))
               end do
            end if
         end do
      end do
      !$omp end do
      !$omp end parallel

      call mpas_timer_stop("overflow tend")

   !--------------------------------------------------------------------

   end subroutine ocn_overflow_tend!}}}

!***********************************************************************
!
!  routine ocn_overflow_init
!
!> \brief   Initializes the bottom boundary layer parameterization
!> \author  Luke Van Roekel
!> \date    October 2026
!> \details
!>  Reads the overflow_parameterization namelist options.
!
!-----------------------------------------------------------------------

   subroutine ocn_overflow_init(err)!{{{

      integer, intent(out) :: err !< [out] error flag

      err = 0

      overflowOn = .false.
      overflowDiffusionOn = .false.
      overflowAdvectionOn = .false.

      if (.not. config_use_overflow) return

      overflowDiffusionOn = config_use_overflow_diffusion
      overflowAdvectionOn = config_use_overflow_advection
      overflowOn = overflowDiffusionOn .or. overflowAdvectionOn

      overflowDiffMin = config_overflow_eddy_diffusivity_minimum
      overflowDiffMax = config_overflow_eddy_diffusivity_maximum
      overflowGravityCoeff = config_overflow_gravity_coefficient
      overflowMaxVolumeFraction = config_overflow_max_volume_fraction
      overflowDensityThreshold = config_overflow_density_threshold

      select case (trim(config_overflow_advection_velocity))
      case ('resolved')
         overflowGravityVelocity = .false.
      case ('gravity')
         overflowGravityVelocity = .true.
      case default
         call mpas_log_write( &
            'ocn_overflow_init: unknown config_overflow_advection_velocity ' // &
            trim(config_overflow_advection_velocity), MPAS_LOG_CRIT)
         err = 1
      end select

      if (overflowMaxVolumeFraction <= 0.0_RKIND .or. &
          overflowMaxVolumeFraction > 0.5_RKIND) then
         call mpas_log_write( &
            'ocn_overflow_init: config_overflow_max_volume_fraction must be in (0, 0.5]', &
            MPAS_LOG_CRIT)
         err = 1
      end if

   end subroutine ocn_overflow_init!}}}

!***********************************************************************

end module ocn_overflow

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
