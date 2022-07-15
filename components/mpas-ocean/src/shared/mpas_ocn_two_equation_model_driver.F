! Copyright (c) 2013-2018,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_two_equation_model
!
!> \brief MPAS ocean two equation model driver
!> \author Luke Van Roekel
!> \date   June 2021
!> \details
!>    This module contains all necessary routines to compute 
!>    tendencies for the two equation turbulence closure 
!>    this is primarily meant to be for nonhydrostatic mode, but could be
!>    extended to the hydrostatic code as a GOTM replacement.
!!
!-----------------------------------------------------------------------

module ocn_two_equation_model

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_kind_types
   use mpas_dmpar
   use mpas_log
   use ocn_constants
   use ocn_mesh
   use ocn_tke_advection_mono
   use ocn_glsPsi_advection_mono
   use ocn_tke_glsPsi_transport
   use ocn_shear_production
   use ocn_buoyancy_production
   use ocn_diagnostics_variables
   use ocn_two_equation_structure_functions
   use ocn_config
   use ocn_tracer_advection_shared

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
   ! Public Member Functions
   !
   !--------------------------------------------------------------------

   public :: ocn_two_equation_init, &
             ocn_compute_two_equation_tendency, &
             ocn_compute_strain, &
             ocn_compute_two_equation_dissipation, &
             ocn_compute_two_equation_diff, &
             ocn_compute_two_equation_molecular_viscosity

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   real(kind=RKIND) :: &
      c_mu, c_psi1, c_psi2, c_psi3, sigma_k, sigma_psi, c_p, c_m, c_n, &
      fkin, feps

   logical :: twoEquationMixingOn
!***********************************************************************

contains

!***********************************************************************
!
!   routine ocn_compute_two_equation_dissipation
!>  \brief    Compute the tendency terms for the k and psi terms 
!>  \author   Luke Van Roekel
!>  \date     June 2021
!>  \details
!>            
!>            
!>            
!>            
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_two_equation_dissipation(tke, glsPsi, tkeTend, glsPsiTend)

     real(kind=RKIND), dimension(:,:), intent(in) :: tke,        & ! tke at current time level
                                                     glsPsi        ! generalized variable from 
                                                                   ! Umlauf and Burchard 2003

     real(kind=RKIND), dimension(:,:), intent(inout) :: tkeTend,    & ! Tendency of TKE
                                                        glsPsiTend    ! tendency of Psi variable

     real(kind=RKIND) :: dissipation, wshear, c2pans

     integer :: iCell, k

     if(config_LES_mode) then
        do iCell=1,nCellsOwned
           do k=minLevelCell(iCell),maxLevelCell(iCell)
              tkeTend(k,iCell) = tkeTend(k,iCell) - 8.0E-6_RKIND*tke(k,iCell)
              glsPsiTend(k,iCell) = 0.0_RKIND
           end do
        end do
     else
     do iCell = 1, nCellsOwned
        do k = 1, maxLevelCell(iCell)
           !
           if ( (abs(shearProduction(k,iCell)) + abs(buoyancyProduction(k,iCell)) ) .ge.  1.0E-10_RKIND) THEN   !FSP_PANS
              wshear = abs(shearProduction(k,iCell)) / ( abs(shearProduction(k,iCell)) + abs(buoyancyProduction(k,iCell)))
           else
              wshear = 0.50_RKIND
           endif
           c2pans =  c_psi2 * fkin/feps &
                  + (c_psi1*wshear + c_psi3*(1.0_RKIND-wshear)) * (1.0_RKIND - fkin/feps)
           dissipation = c_mu**(3.0+c_p/c_n)*tke(k,iCell)**(1.5+c_m/c_n)* &
                         glsPsi(k,iCell)**(-1.0/c_n)
           tkeTend(k,iCell) = tkeTend(k,iCell) - dissipation
           !glsPsiTend(k,iCell) = glsPsiTend(k,iCell) - c_psi2*glsPsi(k,iCell) / &
           glsPsiTend(k,iCell) = glsPsiTend(k,iCell) - c2pans*glsPsi(k,iCell) / &
                                 (tke(k,iCell) + 1.0E-15_RKIND)*dissipation
        end do
     end do
     end if
   end subroutine ocn_compute_two_equation_dissipation

!***********************************************************************
!
!   routine ocn_compute_two_equation_tendency
!>  \brief    Compute the tendency terms for the k and psi terms 
!>  \author   Luke Van Roekel
!>  \date     June 2021
!>  \details
!>            
!>            
!>            
!>            
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_two_equation_tendency(tkeTend, glsPsiTend, tke, glsPsi, velocityZonal, &
                    velocityMeridional, normalVelocity, verticalVelocity, layerThickness,        &
                    activeTracers, dt, err)

     real(kind=RKIND),dimension(:,:), intent(in) :: velocityZonal, velocityMeridional, &
                                                    normalVelocity, layerThickness,    &
                                                    verticalVelocity, tke, glsPsi

     real(kind=RKIND),dimension(:,:), intent(inout) :: tkeTend, glsPsiTend

     real (kind=RKIND), dimension(:,:,:), intent(in) :: activeTracers
     real(kind=RKIND),dimension(:,:), allocatable :: normalThicknessFlux

     integer,intent(inout) :: err

     real(kind=RKIND), intent(in) :: dt

     integer :: iCell, k, nCells, nEdges, iEdge, cell1, cell2, i

     real(kind=RKIND) :: dbdz
     real(kind=RKIND),dimension(:),allocatable :: dbdy

     nCells = nCellsHalo(2)

     do iCell=1,nCells
        do k=1,nVertLevels
           tkeTend(k,iCell) = 0.0_RKIND
           glsPsiTend(k,iCell) = 0.0_RKIND
        end do
     end do

     allocate(dbdy(nVertLevels))
     nEdges = nEdgesAll
     allocate(normalThicknessFlux(nVertLevels,nEdges))
     do iEdge=1,nEdges
        do k=1,nVertLevels
           normalThicknessFlux(k,iEdge) = normalVelocity(k,iEdge)*layerThickEdgeFlux(k,iEdge)
        end do
     end do

     !Compute the strain for the remaining terms in the tendency
     call ocn_compute_strain(velocityZonal, velocityMeridional,  &
                verticalVelocity, activeTracers, layerThickness, err)

     if(.not. config_LES_Mode) call ocn_compute_two_equation_molecular_viscosity(activeTracers)

     tkeAdvTend(:,:) = tkeTend(:,:)

     call ocn_tke_advection_mono_tend(tkeTend, tke, layerThickness,    &
                                      normalThicknessFlux, verticalVelocity, dt,       &
                                      advCoefs, advCoefs3rd,            &
                                      nAdvCellsForEdge, advCellsForEdge,&
                                      advMaskHighOrder)

     tkeAdvTend(:,:) = tkeTend(:,:) - tkeAdvTend(:,:)

     tkeTransTend(:,:) = tkeTend(:,:)
     call ocn_tke_transport_horizontal_compute(tke, tkeTend, layerThickness, sigma_k, fkin, feps, err)
     tkeTransTend(:,:) = tkeTend(:,:) - tkeTransTend(:,:)

     if (.not. config_LES_mode) then

        call ocn_glsPsi_advection_mono_tend(glsPsiTend, glsPsi, layerThickness,    &
                                            normalThicknessFlux, verticalVelocity, dt,       &
                                            advCoefs, advCoefs3rd,            &
                                            nAdvCellsForEdge, advCellsForEdge,&
                                            advMaskHighOrder)

         call ocn_glsPsi_transport_horizontal_compute(glsPsi, glsPsiTend, layerThickness, &
                                        c_psi2, sigma_psi, fkin, feps, err)

     endif

     tkeDissTend(:,:) = tkeTend(:,:)
     call ocn_compute_two_equation_dissipation(tke, glsPsi, tkeTend, glsPsiTend)
     tkeDissTend(:,:) = tkeTend(:,:) - tkeDissTend(:,:)
     call ocn_buoyancy_production_compute(tkeTend, glsPsiTend, glsPsi, tke, c_psi3, &
                                         surfaceBuoyancyForcing, err)
     call ocn_shear_production_compute(tkeTend, glsPsiTend, tke, glsPsi, c_psi1, &
                                       surfaceFrictionVelocity, c_mu, err)

     if (config_two_equation_mesoscale) then
     !add source term to tke tendency equivalent to kappaGM*(M4/N2)
     !   over cells then loop around edges and use edgeareafraction to interpolate to center
        do iCell=1,nCellsAll
              dbdy(:) = 0.0_RKIND
              do i=1,nEdgesOnCell(iCell)
                 iEdge = edgesOnCell(i,iCell)
                 !dbdy = dbdy + edgeAreaFactionOfCell(i,iCell)*( 
                 cell1 = cellsOnEdge(1,iEdge)
                 cell2 = cellsOnEdge(2,iEdge)
                 do k=minLevelEdgeTop(iEdge),maxLevelEdgeTop(iEdge)
                    dbdy(k) = dbdy(k) + edgeAreaFractionOfCell(i,iCell)*(density(k,cell2) - density(k,cell1)) / &
                               dcEdge(iEdge)*gravity/rho_sw
                 end do
              end do

           do k=minLevelCell(iCell),maxLevelCell(iCell)
              if(k==minLevelCell(iCell)) then
                 dbdz = max(BruntVaisalaFreqTop(k+1,iCell),1.0E-20_RKIND)
              elseif(k==maxLevelCell(iCell)) then
                 dbdz = max(BruntVaisalaFreqTop(k,iCell),1.0E-20_RKIND)
              else
                 dbdz = max(0.5_RKIND*(BruntVaisalaFreqTop(k,iCell) + BruntVaisalaFreqTop(k+1,iCell)),1.0E-20_RKIND)
              end if
              if(dbdz > 1e-10_RKIND) then
                 tkeTend(k,iCell) = tkeTend(k,iCell) + config_GM_constant_kappa* & 
                                              max(0.0_RKIND,min(abs(dbdy(k))/dbdz,0.005_RKIND)) * &
                                              abs(dbdy(k))
              end if
           end do
        end do
     end if

     deallocate(dbdy)
     deallocate(normalThicknessFlux)
  end subroutine ocn_compute_two_equation_tendency

!***********************************************************************
!
!   routine ocn_compute_two_equation_diff
!>  \brief    Computes the diffusivity for the two equation model
!>  \author   Luke Van Roekel
!>  \date     May 2021
!>  \details
!>            Computes the diffusivity/viscosity for the two equation model
!>            this is based on the tke and dissipation and diagnosed length 
!>            scale
!>
!>            The diffusivity and viscosity are placed in the vertDiffTopOfCell
!>            and vertViscTopOfCell arrays respectively
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_two_equation_diff(tke, glsPsi, layerThickness, dt, err)

     real(kind=RKIND),dimension(:,:), intent(inout) :: glsPsi, tke

     real(kind=RKIND),dimension(:,:), intent(in) :: layerThickness

     integer,intent(inout) :: err

     real(kind=RKIND), intent(in) :: dt

     real(kind=RKIND),dimension(:),allocatable :: dbdy,N2edge,intValue
     real(kind=RKIND) :: tkeAvg, glsPsiAvg, dbdz !FSP
     real(kind=RKIND) :: length, dissipation, length_grid, length_stable, dcAv
     real(kind=RKIND) :: Lrhines, valEdge, gradH, sumE, sumD, N2ref
     integer :: cell1, cell2, j, i, iCell, k, nCells, nEdges, iEdge, spot(2)

     !compute new diffusivity based on older values
     call ocn_canuto_structure_function_compute(Sm, Sh, tke, glsPsi, c_p, c_m, c_n, c_mu)
     nCells = nCellsAll

     allocate(intValue(nCells),N2edge(nVertLevels),dbdy(nVertLevels))
     if (config_two_equation_mesoscale) then
         do iCell=1,nCellsAll
              dbdy(:) = 0.0_RKIND
              do i=1,nEdgesOnCell(iCell)
                 iEdge = edgesOnCell(i,iCell)
                 !dbdy = dbdy + edgeAreaFactionOfCell(i,iCell)*( 
                 cell1 = cellsOnEdge(1,iEdge)
                 cell2 = cellsOnEdge(2,iEdge)
                 do k=minLevelEdgeTop(iEdge),maxLevelEdgeTop(iEdge)
                    dbdy(k) = dbdy(k) + edgeAreaFractionOfCell(i,iCell)*(density(k,cell2) - density(k,cell1)) / &
                               dcEdge(iEdge)*gravity/rho_sw
                 end do
              end do

           sumE = 0.0_RKIND
           N2ref = max(1.0E-10_RKIND,maxval(BruntVaisalaFreqTop(2:maxLevelCell(iCell)-1,iCell)))
           sumD = 0.0_RKIND
           do k=minLevelCell(iCell),maxLevelCell(iCell)
              dbdz = max(1.0E-10_RKIND,BruntVaisalaFreqTop(k,iCell))
              sumE = sumE + tke(k,iCell)*layerThickness(k,iCell)
              sumD = sumD + dbdz/N2ref*min(abs(dbdy(k))/dbdz,0.01_RKIND)
           end do
           intValue(iCell) = sumE/sumD*0.04_RKIND
        end do
       
        !add in the length scale from Jansen 2015
        do iEdge=1,nEdgesAll
           cell1 = cellsOnEdge(1,iEdge)
           cell2 = cellsOnEdge(2,iEdge)
           gradH = (bottomDepth(cell2) - bottomDepth(cell1))/dcEdge(iEdge)
           gradH = fEdge(iEdge)*gradH / (0.5_RKIND*(bottomDepth(cell1) + bottomDepth(cell2)))
           N2edge(:) = 0.0_RKIND
           do k=minLevelEdgeTop(iEdge),maxLevelEdgeTop(iEdge)
              N2edge(k) = max(1.0E-10_RKIND,0.5_RKIND*(BruntVaisalaFreqTop(k,cell1) + &
                                                       BruntVaisalaFreqTop(k,cell2)))
           enddo
           N2ref = maxval(N2edge)
           do k=minLevelEdgeTOp(iEdge),maxLevelEdgeTop(iEdge)
              valEdge = 0.5_RKIND*(intValue(cell1)+intValue(cell2))*N2edge(k)/N2ref
              mesoEddyDiffusivity(k,iEdge) = max(50.0_RKIND,min(15000.0_RKIND,valEdge))
         end do
         end do
     else
     if (config_LES_mode) then
        vertViscTopOfCell(:,:) = 0.0_RKIND
        vertDiffTopOfCell(:,:) = 0.0_RKIND
        !$omp parallel
        !$omp do schedule(runtime) private(k,iEdge,dcAv,length_grid,length_stable,length)
        do iCell=1,nCellsAll
           dcAv = 0.0_RKIND
           do j=1,nEdgesOnCell(iCell)
              iEdge = edgesOnCell(j, iCell)
              dcAv = dcAv + dcEdge(iEdge)
           end do
           dcAv = dcAv / (1.0E-15_RKIND + nEdgesOnCell(iCell))
           do k=2,maxLevelCell(iCell)
              !compute length scale
              !define the viscosity
              tkeAvg = 0.5_RKIND*(tke(k,iCell) + tke(k-1,iCell))

              length_grid = sqrt(dcAv * 0.5_RKIND*(layerThickness(k-1,iCell) + layerThickness(k,iCell)))
              length_stable = 0.76_RKIND*sqrt(tkeAvg)/ &
                              (sqrt(max(BruntVaisalaFreqTop(k,iCell),0.0_RKIND)) + &
                                                             1.0E-15_RKIND)
              length = min(length_grid, length_stable)
              vertViscTopOfCell(k,iCell) = smagorinsky_coefficient*length*sqrt(tkeAvg)
              vertDiffTopOfCell(k,iCell) = vertViscTopOfCell(k,iCell) !Assume turbulent Pr of 1
           end do
        end do
        !$omp end do
        !$omp end parallel
     else
        !$omp parallel
        !$omp do schedule(runtime) private(tkeAvg,glsPsiAvg,k,length,dissipation)
        do iCell=1,nCells
           vertViscTopOfCell(1,iCell) = 0.0_RKIND
           vertViscTopOfCell(maxLevelCell(iCell)+1:nVertLevels,iCell) = 0.0_RKIND
           do k=2,maxLevelCell(iCell)
              tkeAvg   =0.5_RKIND*( tke(k,iCell)+tke(k-1,iCell) ) !FSP
              glsPsiAvg=0.5_RKIND*( glsPsi(k,iCell)+glsPsi(k-1,iCell) ) !FSP

              !!!dissipation = c_mu**(3.0+c_p/c_n)*tke(k,iCell)**(1.5+c_m/c_n)* &
              !!!              glsPsi(k,iCell)**(-1.0/c_n)
              !!!length = c_mu**3*tke(k,iCell)**1.5/(1.0E-15_RKIND + dissipation)
              !!!vertViscTopOfCell(k,iCell) = length*Sm(k,iCell)*sqrt(tke(k,iCell)) + &
              dissipation = c_mu**(3.0+c_p/c_n)*tkeAvg**(1.5+c_m/c_n)* &
                            glsPsiAvg**(-1.0/c_n)
              length = c_mu**3*tkeAvg**1.5/(1.0E-15_RKIND + dissipation)
              !vertViscTopOfCell(k,iCell) = length*Sm(k,iCell)*sqrt(tkeAvg) + &
              vertViscTopOfCell(k,iCell) = length*sqrt(tkeAvg) + &
                                           dynamicViscosity(k,iCell)/rho_sw
           end do
           !FSP
           vertViscTopOfCell(1,iCell) = vertViscTopOfCell(2,icell)
           vertViscTopOfCell(maxLevelCell(iCell)+1,iCell) = vertViscTopOfCell(maxLevelCell(iCell),icell)
        end do
        !$omp end do 
        !$omp end parallel
     end if

     call ocn_tke_transport_vertical_compute(tke,layerThickness,sigma_k,fkin,feps,dt)
     if(.not. config_LES_mode) then
        call ocn_glsPsi_transport_vertical_compute(glsPsi,layerThickness,sigma_psi,fkin,feps,dt)
     end if

     if (config_LES_mode) then
        !$omp parallel
        !$omp do schedule(runtime) private(k,dcAv,iEdge,length_grid,length_stable,length)
        do iCell=1,nCells
           dcAv = 0.0_RKIND
           do j=1,nEdgesOnCell(iCell)
              iEdge = edgesOnCell(j, iCell)
              dcAv = dcAv + dcEdge(iEdge)
           end do
           dcAv = dcAv / (1.0E-15_RKIND + nEdgesOnCell(iCell))
           do k=2,maxLevelCell(iCell)
              !compute length scale
              !define the viscosity
              length_grid = sqrt(dcAv * 0.5_RKIND*(layerThickness(k-1,iCell) + layerThickness(k,iCell)))
              length_stable = 0.76_RKIND*sqrt(max(tke(k,iCell),0.0_RKIND))/(sqrt(max(BruntVaisalaFreqTop(k,iCell),0.0_RKIND)) + &
                                                             1.0E-15_RKIND)
              length = min(length_grid, length_stable)
              vertViscTopOfCell(k,iCell) = smagorinsky_coefficient*length*sqrt(tke(k,iCell)) !+ &
                     !         dynamicViscosity(k,iCell)/rho_sw
              vertDiffTopOfCell(k,iCell) = vertViscTopOfCell(k,iCell) !Assume turbulent Pr of 1
           end do
        end do
        !$omp end do
        !$omp end parallel
     else
       ! Recompute now that everything is fully up to date
        call ocn_canuto_structure_function_compute(Sm, Sh, tke, glsPsi, c_p, c_m, c_n, c_mu)
        !$omp parallel
        !$omp do schedule(runtime) private(tkeAvg,glsPsiAvg,k,dissipation,length)
        do iCell=1,nCells
           vertDiffTopOfCell(1,iCell) = 0.0_RKIND
           vertViscTopOfCell(1,iCell) = 0.0_RKIND
           vertDiffTopOfCell(maxLevelCell(iCell)+1:nVertLevels,iCell) = 0.0_RKIND
           vertViscTopOfCell(maxLevelCell(iCell)+1:nVertLevels,iCell) = 0.0_RKIND
           do k=2,maxLevelCell(iCell)
              tkeAvg   =0.5_RKIND*( tke(k,iCell)+tke(k-1,iCell) ) !FSP
              glsPsiAvg=0.5_RKIND*( glsPsi(k,iCell)+glsPsi(k-1,iCell) ) !FSP
              dissipation = c_mu**(3.0+c_p/c_n)*tkeAvg**(1.5+c_m/c_n)* &
                            glsPsiAvg**(-1.0/c_n)
              length = c_mu**3*tkeAvg**1.5/(1.0E-15_RKIND + dissipation)
              !vertViscTopOfCell(k,iCell) = length*Sm(k,iCell)*sqrt(tke(k,iCell)) 
              !vertViscTopOfCell(k,iCell) = length*Sm(k,iCell)*sqrt(tkeAvg) + &
              vertViscTopOfCell(k,iCell) = length*sqrt(tkeAvg) + &
                                           dynamicViscosity(k,iCell)/rho_sw
              !!!vertDiffTopOfCell(k,iCell) = length*Sh(k,iCell)*sqrt(tkeAvg)
              vertDiffTopOfCell(k,iCell) = length*sqrt(tkeAvg)


              !!!dissipation = c_mu**(3.0+c_p/c_n)*tke(k,iCell)**(1.5+c_m/c_n)* &
              !!!              glsPsi(k,iCell)**(-1.0/c_n)
              !!!length = c_mu**3*tke(k,iCell)**1.5/(1.0E-15_RKIND + dissipation)
              !!!!vertViscTopOfCell(k,iCell) = length*Sm(k,iCell)*sqrt(tke(k,iCell)) 
              !!!vertViscTopOfCell(k,iCell) = length*Sm(k,iCell)*sqrt(tke(k,iCell)) + &
              !!!                             dynamicViscosity(k,iCell)/rho_sw
              !!!vertDiffTopOfCell(k,iCell) = length*Sh(k,iCell)*sqrt(tke(k,iCell))
           end do
           !FSP
           vertViscTopOfCell(1,iCell) = vertViscTopOfCell(2,icell)
           vertViscTopOfCell(maxLevelCell(iCell)+1,iCell) = vertViscTopOfCell(maxLevelCell(iCell),icell)
           vertDiffTopOfCell(1,iCell) = vertDiffTopOfCell(2,icell)
           vertDiffTopOfCell(maxLevelCell(iCell)+1,iCell) = vertDiffTopOfCell(maxLevelCell(iCell),icell)
        end do
        !$omp end do 
        !$omp end parallel
     end if

   end if !meso_eddy closure if

   deallocate(N2edge, intValue, dbdy)
   end subroutine ocn_compute_two_equation_diff

!***********************************************************************
!
!   routine ocn_compute_strain
!>  \brief    Computes the strain rate tensor magnitude for smagorinsky 
!>  \author   Luke Van Roekel
!>  \date     May 2021
!>  \details
!>            Computes the strain rate tensor at cell centers.  Uses
!>            the vector reconstruct coefficients to get x and y derivs
!>            For 2D need dudx, dudy, dvdx, dvdy
!>            For 3D adding dwdx and dwdy
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_strain(velocityZonal, velocityMeridional,  &
                verticalVelocity, activeTracers, layerThickness, err)

     real(kind=RKIND),dimension(:,:),intent(in) :: &
         velocityZonal, &  !Input: Zonal Velocity for gradients
         velocityMeridional, & !Input: Meridional Velocity for gradients
         layerThickness

     real(kind=RKIND),dimension(:,:),intent(in),optional :: &
         verticalVelocity !Input: optionally for nonhydrostaticMode include w

     real(kind=RKIND),dimension(:,:,:),intent(in) :: &
         activeTracers

     integer, intent(in) :: err

     integer :: iCell, iEdge, j, k, cell1, cell2

     real(kind=RKIND) :: eSign, coeff1, coeff2, factor

     ! initialize gradients

     do iCell=1,nCellsAll
        ugrad(:,:,iCell) = 0.0_RKIND
        vgrad(:,:,iCell) = 0.0_RKIND
        tgrad(:,:,iCell) = 0.0_RKIND

        do j = 1,nEdgesOnCell(iCell) !for device may be better to move next lines inside k loop
           iEdge = edgesOnCell(j,iCell)
           cell1 = cellsOnEdge(1,iEdge)
           cell2 = cellsOnEdge(2,iEdge)
           eSign = edgeSignOnCell(j,iCell)
           coeff1 = coeffs_reconstruct(1,j,iCell)
           coeff2 = coeffs_reconstruct(2,j,iCell)

           do k = 1,maxLevelEdgeTop(iEdge)
              tgrad(1,k,iCell) = tgrad(1,k,iCell) + coeff1 * eSign * & 
                            (activeTracers(1,k,cell2) - activeTracers(1,k,cell1)) / dcEdge(iEdge)
              tgrad(2,k,iCell) = tgrad(2,k,iCell) + coeff2 * eSign * & 
                            (activeTracers(1,k,cell2) - activeTracers(1,k,cell1)) / dcEdge(iEdge)
              sgrad(1,k,iCell) = sgrad(1,k,iCell) + coeff1 * eSign * & 
                            (activeTracers(2,k,cell2) - activeTracers(2,k,cell1)) / dcEdge(iEdge)
              sgrad(2,k,iCell) = sgrad(2,k,iCell) + coeff2 * eSign * &
                            (activeTracers(2,k,cell2) - activeTracers(2,k,cell1)) / dcEdge(iEdge)
              ugrad(1,k,iCell) = ugrad(1,k,iCell) + coeff1 * eSign * &
                           (velocityZonal(k,cell2) - velocityZonal(k,cell1)) / dcEdge(iEdge)
              ugrad(2,k,iCell) = ugrad(2,k,iCell) + coeff2 * eSign * &
                            (velocityZonal(k,cell2) - velocityZonal(k,cell1)) / dcEdge(iEdge)
              vgrad(1,k,iCell) = vgrad(1,k,iCell) + coeff1 * eSign * &
                            (velocityMeridional(k,cell2) - velocityMeridional(k,cell1)) /         &
                            dcEdge(iEdge)
              vgrad(2,k,iCell) = vgrad(2,k,iCell) + coeff2 * eSign * &
                            (velocityMeridional(k,cell2) - velocityMeridional(k,cell1)) /         &
                            dcEdge(iEdge)

          end do
          if(maxLevelEdgeTop(iEdge) == 0) then
             if (cell1 == iCell) then
                factor = -1.0_RKIND
             else
                factor = 1.0_RKIND
             end if

             do k=1,maxLevelCell(iCell)
!               ugrad(1,k,iCell) = ugrad(1,k,iCell) + abs(eSign*coeff1*velocityZonal(k,iCell)/dcEdge(iEdge))
               ugrad(2,k,iCell) = ugrad(2,k,iCell) + abs(velocityZonal(k,iCell)/dcEdge(iEdge))
!               vgrad(1,k,iCell) = vgrad(1,k,iCell) + abs(eSign*coeff1*velocityMeridional(k,iCell)/dcEdge(iEdge))
!               vgrad(2,k,iCell) = vgrad(2,k,iCell) + abs(eSign*coeff2*velocityMeridional(k,iCell)/dcEdge(iEdge))
            end do
          end if
      end do !end n edges loop
       ! instead of doing this put in surface BCs for generation
       ! compute dudz dvdz
       ugrad(3,1,iCell) = (velocityZonal(1,iCell) - velocityZonal(2,iCell)) /   &
                            (0.5_RKIND*(layerThickness(1,iCell) + layerThickness(2,iCell)))
       vgrad(3,1,iCell) = (velocityMeridional(1,iCell) - velocityMeridional(2,iCell)) /   &
                            (0.5_RKIND*(layerThickness(1,iCell) + layerThickness(2,iCell)))
       tgrad(3,1,iCell) = (activeTracers(1,1,iCell) - activeTracers(1,2,iCell)) /   &
                            (0.5_RKIND*(layerThickness(1,iCell) + layerThickness(2,iCell)))
       sgrad(3,1,iCell) = (activeTracers(2,1,iCell) - activeTracers(2,2,iCell)) /   &
                            (0.5_RKIND*(layerThickness(1,iCell) + layerThickness(2,iCell)))

       do k=2,maxLevelCell(iCell)
          tgrad(3,k,iCell) = 0.5_RKIND*(tgrad(3,k-1,iCell) + (activeTracers(1,k-1,iCell) - activeTracers(1,k,iCell)) / &
                        (0.5_RKIND*(layerThickness(k-1,iCell) + layerThickness(k,iCell))))
          sgrad(3,k,iCell) = 0.5_RKIND*(sgrad(3,k-1,iCell) + (activeTracers(2,k-1,iCell) - activeTracers(2,k,iCell)) / &
                        (0.5_RKIND*(layerThickness(k-1,iCell) + layerThickness(k,iCell))))
          ugrad(3,k,iCell) = 0.5_RKIND*(ugrad(3,k-1,iCell) + (velocityZonal(k-1,iCell) - velocityZonal(k,iCell)) / &
                        (0.5_RKIND*(layerThickness(k-1,iCell) + layerThickness(k,iCell))))
          vgrad(3,k,iCell) = 0.5_RKIND*(vgrad(3,k-1,iCell) + (velocityMeridional(k-1,iCell) - velocityMeridional(k,iCell)) / &
                        (0.5_RKIND*(layerThickness(k-1,iCell) + layerThickness(k,iCell))))
       end do

    end do ! end cell loop

    ! add w part if requested 
    if(present(verticalVelocity)) then
       do iCell=1,nCellsOwned
          wgrad(:,:,iCell) = 0.0_RKIND

          do j = 1,nEdgesOnCell(iCell) !for device may be better to move next lines inside k loop
             iEdge = edgesOnCell(j,iCell)
             cell1 = cellsOnEdge(1,iEdge)
             cell2 = cellsOnEdge(2,iEdge)
         
             coeff1 = coeffs_reconstruct(1,j,iCell)
             coeff2 = coeffs_reconstruct(2,j,iCell)
             eSign = edgeSignOnCell(j,iCell)
             do k=1,maxLevelEdgeTop(iEdge)
                wgrad(1,k,iCell) = wgrad(1,k,iCell) + coeff1 * eSign * & 
                             (verticalVelocity(k,cell2) - verticalVelocity(k,cell1)) / dcEdge(iEdge)
                wgrad(2,k,iCell) = wgrad(2,k,iCell) + coeff2 * eSign * & 
                             (verticalVelocity(k,cell2) - verticalVelocity(k,cell1)) / dcEdge(iEdge)
             end do

             if(maxLevelEdgeTop(iEdge) == 0) then
                do k=1,maxLevelCell(iCell)
                  wgrad(1,k,iCell) = wgrad(1,k,iCell) - coeff1*eSign*verticalVelocity(k,iCell)/dcEdge(iEdge)
                  wgrad(2,k,iCell) = wgrad(2,k,iCell) - coeff2*eSign*verticalVelocity(k,iCell)/dcEdge(iEdge)
                end do
             end if
          end do !end n edges loop

          do k=2,maxLevelCell(iCell)
             wgrad(3,k,iCell) = (verticalVelocity(k-1,iCell) - verticalVelocity(k,iCell)) / &
                         ((layerThickness(k,iCell)))
          end do

       end do !end icell loop
    end if !if present w 

 end subroutine ocn_compute_strain

!***********************************************************************
!
!   routine ocn_two_equation_init
!>  \brief    Configures options for the two equation model 
!>  \author   Luke Van Roekel
!>  \date     May 2021
!>  \details
!>            Initializes variables and allocates internal variables for 
!>            computing strain tensor and k and psi 
!
!----------------------------------------------------------------------

   subroutine ocn_two_equation_init( domain, err )

      type (domain_type), intent(inout) :: domain
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool
      integer, intent(out) :: err
      real(kind=RKIND) ,dimension(:,:), pointer :: glsPsi, tke         !FSP
      real(kind=RKIND) ,dimension(:,:), pointer :: glsPsiNew, tkeNew   !FSP
      real(kind=RKIND), dimension(:,:), pointer :: layerThickness
      real(kind=RKIND), dimension(:,:), pointer :: distanceToBoundary
      real(kind=RKIND), dimension(:,:), pointer :: boundaryNormalAngle

      integer :: cell1, cell2, iEdge, k

      call mpas_pool_get_subpool(domain%blocklist%structs,'mesh',meshPool)
      call mpas_pool_get_subpool(domain%blocklist%structs,'state',statePool)
 
      twoEquationMixingOn = .false.
      if(config_use_two_equation_turbulence_model) then

         twoEquationMixingOn = .true.

         !PANS paremeters to set range of resolved scales
         fkin=min(max(fkin_value,0.05_RKIND),1.0_RKIND)
         feps=min(max(feps_value,0.05_RKIND),1.0_RKIND)

         if (config_two_equation_model_choice == 'kl') then !MY 1982
            c_mu = 0.5544_RKIND
            sigma_k = 1.96_RKIND
            sigma_psi = 1.96_RKIND
            c_psi1 = 0.9_RKIND
            c_psi2 = 0.5_RKIND
            c_psi3 = 0.9_RKIND
            c_p = 0.0_RKIND
            c_m = 1.0_RKIND
            c_n = 1.0_RKIND
         elseif (config_two_equation_model_choice == 'epsilon') then ! Rotta 87
            c_mu = 0.5477_RKIND
            sigma_k = 1.0_RKIND
            sigma_psi = 1.3_RKIND
            c_psi1 = 1.44_RKIND
            c_psi2 = 1.92_RKIND
            c_psi3 = -0.629_RKIND
            c_p = 3.0_RKIND
            c_m = 3.0_RKIND / 2.0_RKIND
            c_n = -1.0_RKIND
        elseif (config_two_equation_model_choice == 'omega') then !Wilcox 88
            c_mu = 0.5477_RKIND
            sigma_k = 2.0_RKIND
            sigma_psi = 2.0_RKIND
            c_psi1 = 0.555_RKIND
            c_psi2 = 0.833_RKIND
            c_psi3 = -0.642_RKIND
            c_p = -1.0_RKIND
            c_m = 0.5_RKIND
            c_n = -1.0_RKIND
         elseif (config_two_equation_model_choice == 'ktau') then !Zilinkivech Wilcox 88
            c_mu = 0.5477_RKIND
            sigma_k = 1.46_RKIND
            sigma_psi = 10.8_RKIND
            c_psi1 = 0.173_RKIND
            c_psi2 = 0.225_RKIND
            c_psi3 = 0.0_RKIND
            c_p = -3.0_RKIND
            c_m = 0.5_RKIND
            c_n = 1.0_RKIND
         elseif (config_two_equation_model_choice == 'user-defined') then
            c_mu = config_user_c_mu
            sigma_k = config_user_sigma_k
            sigma_psi = config_user_sigma_psi
            c_psi1 = config_user_c_psi1
            c_psi2 = config_user_c_psi2
            c_psi3 = config_user_c_psi3
            c_p = config_user_c_p
            c_m = config_user_c_m
            c_n = config_user_c_n
         else !invalid choice 
            call mpas_log_write("Unknown option for config_two_equation_model_choice, valid values are: " // &
                                  " 'epsilon', 'kl', 'omega', 'ktau', 'user-defined'", MPAS_LOG_CRIT)
         endif

         if(config_use_cvmix .and. .not. config_two_equation_mesoscale) then
            call mpas_log_write("Cannot use two equation turbulence closure and cvmix at the same time: " // &
                    "set either config_use_cvmix or config_use_two_equation_model to .false.", MPAS_LOG_CRIT)
         end if

         !Initial conditions for tke and glsPsi
         call mpas_pool_get_array(statePool,'tke',tke,1)             !FSP
         call mpas_pool_get_array(statePool,'tke',tkeNew,2)       !FSP
         call mpas_pool_get_array(statePool,'glsPsi',glsPsi,1)       !FSP
         call mpas_pool_get_array(statePool,'glsPsi',glsPsiNew,2) !FSP
         tke(:,:)       = config_tke_ini                                                          !FSP
         tkeNew(:,:)    = tke(:,:)                                                                !FSP
         glsPsi(:,:)    = 1.0E-30_RKIND!( config_nut_ini**c_n) / ( c_mu**(-c_p)*config_tke_ini**(c_n/2-c_m) )   !FSP 
         glsPsiNew(:,:) = glsPsi(:,:)                                                             !FSP 


         !init advection of tke and glsPsi

         call ocn_tke_advection_mono_init(config_num_halos, config_two_equation_hor_adv_order,        &
                                                config_two_equation_vert_adv_order,              &
                                                config_two_equation_coef_3rd_order,              &
                                                config_two_equation_check_monotonicity, err)

         call ocn_glsPsi_advection_mono_init(config_num_halos, config_two_equation_hor_adv_order,        &
                                                config_two_equation_vert_adv_order,                 &
                                                config_two_equation_coef_3rd_order,                 &
                                                config_two_equation_check_monotonicity, err)

         call mpas_pool_get_array(meshPool, 'distanceToBoundary', distanceToBoundary)
         call mpas_pool_get_array(meshPool, 'boundaryNormalAngle', boundaryNormalAngle)
         if(.not. associated(distanceToBoundary)) then
            distanceToBoundary(:,:) = 1.0E50_RKIND
            ! make an assumption that the distance to boundary is just dcEdge/2.0 when edgeMask==0
            do iEdge=1,nEdgesAll
               do k=1,nVertLevels
                  if(edgeMask(k,iEdge) == 0) then
                     cell1 = cellsOnEdge(1,iEdge)
                     cell2 = cellsOnedge(2,iEdge)
                     distanceToBoundary(k,cell1) = dcEdge(iEdge)/2.0_RKIND
                     distanceToBoundary(k,cell2) = dcEdge(iEdge)/2.0_RKIND
                  endif
               end do
            end do
         end if

         if(.not. associated(boundaryNormalAngle)) then
            if(config_two_equation_model_choice == 'epsilon') then
               call mpas_log_write("Unable to use k-epsilon configuration without specifying the " // &
                        "boundaryNormalAngle array", MPAS_LOG_CRIT)
            end if
         end if
      end if

   end subroutine ocn_two_equation_init

   subroutine ocn_compute_two_equation_molecular_viscosity(activeTracers)

     !computes the molecular viscosity of water following Sharqawy et al (2010)
     real(kind=RKIND),dimension(:,:,:), intent(in) :: activeTracers
     real(kind=RKIND) :: A, B, C, dcEdgeAverage
     real(kind=RKIND) :: tAvg, sAvg, viscAvg, A1, B1, freshwaterMu

     real(kind=RKIND) :: boundaryTangentialVel
     integer :: iCell, j, k

     A = 2.414E-5_RKIND
     B = 247.8_RKIND
     C = 140.0_RKIND

     do iCell=1,nCellsAll
        !compute dcEdgeAverage
        dcEdgeAverage = 0.0_RKIND
        do j=1,nEdgesOnCell(iCell)
           dcEdgeAverage = dcEdgeAverage + dcEdge(edgesOnCell(j,iCell))
        end do
        dcEdgeAverage = dcEdgeAverage / nEdgesOnCell(iCell) + 1.0E-15_RKIND
        dynamicViscosity(1,iCell) = 0.0_RKIND
        do k = 1,maxLevelCell(iCell)-1
           tAvg = 0.5_RKIND*(activeTracers(1,k,iCell) + activeTracers(1,k+1,iCell))
           sAvg = 0.5_RKIND*(activeTracers(2,k,iCell) + activeTracers(2,k+1,iCell)) / 1000.0_RKIND
!!!           viscAvg = 1.0E-15_RKIND + 0.5_RKIND*(dynamicViscosity(k,iCell) + dynamicViscosity(k+1,iCell))
!           freshwaterMu = A*10.0_RKIND**(B/(tAvg + 273.15_RKIND - C))
           freshwaterMu = 4.2844E-5 + (0.157_RKIND*(tAvg + 64.993_RKIND)**2.0 - 91.296_RKIND)**(-1.0)
           A1 = 1.541_RKIND + 0.01998*tAvg - 9.52E-5_RKIND*tAvg**2.0
           B1 = 7.974_RKIND - 0.07561*tAvg + 4.724E-4_RKIND*tAvg**2.0
           dynamicViscosity(k+1,iCell) = freshwaterMu*(1.0_RKIND + &
                                        A1*sAvg + &
                                        B1*sAvg**2)
          if(dynamicViscosity(k+1,iCell) < 0) then
            print *, activeTracers(1,:,iCell)
            print *, 'hi',k, tavg,savg, viscavg, freshwaterMu, A1, B1, dynamicViscosity(k+1,iCell)
            stop
          endif
!!!           boundaryTangentialVel = velocityZonal(k,iCell)*sin(boundaryNormalAngle(k,iCell)) + &
!!!                                  velocityMeridional(k,iCell)*cos(boundaryNormalAngle(k,iCell))
!!!           !need to find the edge in question and then use that dcEdge -- maybe just do average for performance reasons?
!!!
!!!           uStar(k,iCell) = sqrt(boundaryCellMask(k,iCell)*abs(boundaryTangentialVel/(0.5_RKIND*dcEdgeAverage)) &
!!!                          * viscAvg/rho_sw)
!!!           yPlus(k,iCell) = uStar(k,iCell)*distanceToBoundary(k,iCell) / (viscAvg/rho_sw)
        end do
        dynamicViscosity(1,iCell) = dynamicViscosity(2,iCell)
        k=maxLevelCell(iCell)
        dynamicViscosity(maxLevelCell(iCell)+1:nVertLevels+1,iCell) = dynamicViscosity(maxLevelCell(iCell),iCell)
!
        do k = 1,maxLevelCell(iCell)-1
           viscAvg = 0.5_RKIND*(dynamicViscosity(k,iCell) + dynamicViscosity(k+1,iCell))
           boundaryTangentialVel = velocityZonal(k,iCell)*sin(boundaryNormalAngle(k,iCell)) + &
                                   velocityMeridional(k,iCell)*cos(boundaryNormalAngle(k,iCell))
           !need to find the edge in question and then use that dcEdge -- maybe just do average for performance reasons?

           uStar(k,iCell) = sqrt(boundaryCellMask(k,iCell)*abs(boundaryTangentialVel/(0.5_RKIND*dcEdgeAverage)) &
                          * viscAvg/rho_sw)
           yPlus(k,iCell) = uStar(k,iCell)*distanceToBoundary(k,iCell) / (viscAvg/rho_sw)
        end do
!
        k=maxLevelCell(iCell)
        viscAvg = 0.5_RKIND*(dynamicViscosity(k,iCell) + dynamicViscosity(k+1,iCell))
        boundaryTangentialVel = velocityZonal(k,iCell)*sin(boundaryNormalAngle(k,iCell)) + &
                                velocityMeridional(k,iCell)*cos(boundaryNormalAngle(k,iCell))
        !need to find the edge in question and then use that dcEdge -- maybe just do average for performance reasons?

        uStar(k,iCell) = sqrt(boundaryCellMask(k,iCell)*abs(boundaryTangentialVel/(0.5_RKIND*dcEdgeAverage)) &
                          * viscAvg/rho_sw)
        yPlus(k,iCell) = uStar(k,iCell)*distanceToBoundary(k,iCell) / (viscAvg/rho_sw)

     end do

     !compute ustar

   end subroutine ocn_compute_two_equation_molecular_viscosity
end module ocn_two_equation_model

