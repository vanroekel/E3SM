!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vmix_tke
!
!> \brief MPAS ocean vertical mixing TKE scheme
!> \author Todd Ringler
!> \date   29 June 2022
!> \details
!>  This module contains the routines to cmpute diffusivity and 
!>   viscosity based on the Gaspar 1990 TKE scheme
!
!-----------------------------------------------------------------------

module ocn_vmix_tke

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timer
   use mpas_constants
   use mpas_log

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

   public :: ocn_vmix_tke_compute, &
             ocn_vmix_tke_init, &
			 ocn_vmix_tke_finalize

   private :: ocn_vmix_tke_compute_length, &
              ocn_vmix_tke_compute_coeffs

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   real(kind=RKIND),dimension(:,:),allocatable :: lup, ldo, niwFactor
   real(kind=RKIND) :: ck, ceps, emin0, emin, alpha, lmin0, lmin
   real(kind=RKIND) :: nu0, kappa0, kappaC, gam, alpha_conv, tkeEnhance

   contains
!***********************************************************************
!
!  routine ocn_vmix_tke_compute_coeffs
!
!> \brief   Computes mixing coefficients using Gaspar 1990
!> \author  Luke Van Roekel
!> \date    29 June 2022
!> \details
!> 
!>  
!
!-----------------------------------------------------------------------

   subroutine ocn_vmix_tke_compute_coeffs(tke,compute_convect)

     real(kind=RKIND),dimension(:,:),intent(in) :: tke

     logical,intent(in) :: compute_convect

     integer :: iCell, k, nCells
     real(kind=RKIND) :: lk, Pr, tkeVal

     nCells = nCellsHalo(2)

     do iCell=1,nCells
        do k=1,nVertLevels+1
           vertDiffTopOfCell(k,iCell) = 0.0_RKIND
           vertViscTopOfCell(k,iCell) = 0.0_RKIND
        end do
     end do

     do iCell = 1,nCells
        do k=minLevelCell(iCell)+1,maxLevelCell(iCell)
           lk = min(lup(k,iCell),ldo(k,iCell))
           if(RiTopOfCell(k,iCell) <= 0.2_RKIND) then
              Pr = 1.0_RKIND
           elseif(RiTopOfCell(k,iCell) > 0.2_RKIND .and. RiTopOfCell(k,iCell) <= 2.0_RKIND) then
              Pr = 5.0_RKIND*RiTopOfCell(k,iCell)
           else
              Pr = 10.0_RKIND
           end if

          ! if(BruntVaisalaFreqTop(k,iCell) < 0.0_RKIND) then
          !    vertViscTopOfCell(k,iCell) = kappaC
          !    vertDiffTopOfCell(k,iCell) = kappaC
          ! else
              tkeVal = max(tke(k,iCell),emin)
              vertViscTopOfCell(k,iCell) = ck*lk*sqrt(tkeVal)
              vertDiffTopOfCell(k,iCell) = vertViscTopOfCell(k,iCell) / Pr
          ! end if
        end do
     end do

     if(compute_convect) then
     do iCell=1,nCells
        do k=minLevelCell(iCell)+1,maxLevelCell(iCell)
           if(BruntVaisalaFreqTop(k,iCell)<0.0_RKIND) then
              vertViscTopOfCell(k,iCell) = kappaC
              vertDiffTopOfCell(k,iCell) = kappaC
           end if
        end do
     end do

   end if
   end subroutine ocn_vmix_tke_compute_coeffs

!***********************************************************************
!
!  routine ocn_vmix_tke_compute_length
!
!> \brief   Computes length scale using Gaspar 1990
!> \author  Luke Van Roekel
!> \date    29 June 2022
!> \details
!> 
!>  
!
!-----------------------------------------------------------------------

   subroutine ocn_vmix_tke_compute_length(tke,layerThickness,ssh,windStressZonal,windStressMeridional)

      real(kind=RKIND),dimension(:,:),intent(in) :: tke,layerThickness

      real(kind=RKIND),dimension(:),intent(in) :: ssh,windStressZonal,windStressMeridional

      integer :: nCells, iCell, k

      nCells = nCellsHalo(2)

      do iCell = 1,nCells
         do k=minLevelCell(iCell)+1,maxLevelCell(iCell)
            if(BruntVaisalaFreqTop(k,iCell) > 0.0_RKIND) then
               lup(k,iCell) = max(lmin,sqrt(2.0_RKIND*tke(k,iCell)/BruntVaisalaFreqTop(k,iCell)))
               ldo(k,iCell) = lup(k,iCell)
            else
               lup(k,iCell) = ssh(iCell) - zTop(k,iCell)
               ldo(k,iCell) = zTop(k,iCell) + bottomDepth(iCell)
            end if
         end do
      end do

      ! impose surface and bottom boundary conditions
      do iCell=1,nCells
         k = minLevelCell(iCell)
         lup(k,iCell) = max(0.41_RKIND*2.0E5_RKIND*sqrt(windStressZonal(iCell)**2 + &
                                                    windStressMeridional(iCell)**2) / &
                                               (rho_sw*gravity),lmin0)
         ldo(k,iCell) = lup(k,iCell)

         k = maxLevelCell(iCell)
         lup(k,iCell) = lmin
         ldo(k,iCell) = lmin
      end do

      ! limit the length scale based on derivative. start from top and do ldown and bottom up with lup

      do iCell=1,nCells
         do k=minLevelCell(iCell)+1,maxLevelCell(iCell)
            ldo(k,iCell) = min(ldo(k,iCell),ldo(k-1,iCell)+layerThickness(k-1,iCell))
         end do
      end do

      do iCell=1,nCells
         do k=maxLevelCell(iCell),minLevelCell(iCell),-1
            lup(k,iCell) = min(lup(k,iCell),lup(k+1,iCell)+layerThickness(k,iCell))
         end do
      end do

   end subroutine ocn_vmix_tke_compute_length

!
!  routine ocn_vmix_tke_init
!
!> \brief   Computes length scale using Gaspar 1990
!> \author  Luke Van Roekel
!> \date    29 June 2022
!> \details
!> 
!>  
!
!-----------------------------------------------------------------------

   subroutine ocn_vmix_tke_init(err)
      integer,intent(inout) :: err

      integer :: k, iCell

      real(kind=RKIND) :: efold

      err = 0

      if(config_use_vmix_tke .and. config_use_cvmix) then
         call mpas_log_write("Both TKE and CVMIX schemes are enabled, must chose one or the other",MPAS_LOG_CRIT)
      end if

      allocate(lup(nVertLevels+1,nCellsAll),ldo(nVertLevels+1,nCellsAll))
      allocate(niwFactor(nVertLevels+1,nCellsOwned))

      ck = config_vmix_tke_ck
      ceps = config_vmix_tke_ce
      emin0 = config_vmix_tke_emin_sfc
      emin = config_vmix_tke_emin_int
      lmin0 = config_vmix_tke_lenmin_sfc
      lmin = config_vmix_tke_lenmin_int
      alpha = config_vmix_tke_wb_alpha
      gam = config_vmix_tke_niw_gamma
      nu0 = config_vmix_tke_nu_bkg
      kappa0 = config_vmix_tke_kappa_bkg
      kappaC = config_vmix_tke_kappa_conv
      tkeEnhance = config_vmix_tke_transport_enhancement_factor
      alpha_conv = config_vmix_tke_sfc_convective_scaling
      efold = config_vmix_tke_niw_efold_depth

      !For now simply have a fixed in space efolding depth for simplicity 
      do iCell=1,nCellsOwned
!         if(latCell(iCell) < -0.6981317_RKIND) then
!            efold = 30.0_RKIND
!         else if(latCell(icell) >= -0.6981317_RKIND .and. latCell(iCell) <= 0.0_RKIND) then
!            efold = -25.0_RKIND/0.6981317_RKIND*latCell(iCell) + 5.0_RKIND
!         else if(latCell(iCell) > 0.0_RKIND .and. latCell(iCell) <= 0.2617993878_RKIND) then
!            efold = 5.0_RKIND/0.2617993878_RKIND*latCell(iCell) + 5.0_RKIND
!         else
!            efold = 10.0_RKIND
!         end if
         do k=minLevelCell(iCell),maxLevelCell(iCell)
            niwFactor(k,iCell) = gam*exp(zTop(k,iCell) / efold)
         end do
         niwFactor(1,iCell) = 0.0_RKIND
      end do

   end subroutine ocn_vmix_tke_init

!  routine ocn_vmix_tke_finalize
!
!> \brief   Computes length scale using Gaspar 1990
!> \author  Luke Van Roekel
!> \date    29 June 2022
!> \details
!>  
!>  
!
!-----------------------------------------------------------------------

   subroutine ocn_vmix_tke_finalize()

      deallocate(lup,ldo,niwFactor)

   end subroutine ocn_vmix_tke_finalize

!  routine ocn_vmix_tke_compute
!
!> \brief   Computes length scale using Gaspar 1990
!> \author  Luke Van Roekel
!> \date    29 June 2022
!> \details
!> 
!>  
!
!-----------------------------------------------------------------------

   subroutine ocn_vmix_tke_compute(statePool, forcingPool, dt, err)

      type (mpas_pool_type), intent(inout) :: statePool
      type (mpas_pool_type), intent(inout) :: forcingPool

      integer, intent(out) :: err


      real(kind=RKIND),dimension(:,:),pointer :: tkeNew, tkeOld, layerThickness

      real(kind=RKIND),dimension(:),pointer :: windStressZonal, windStressMeridional, ssh

      real(kind=RKIND),intent(in) :: dt

      real(kind=RKIND),dimension(:), allocatable :: bTemp, C, rTemp
      real(kind=RKIND) :: stressMag, sfcTKE, shearProduction,  &
                          dissipation, buoyancyProduction, A, m, length, nuMid

      integer :: nCells, iCell, k, N, Nsurf

      if(.not. config_use_vmix_tke) return


      call mpas_timer_start('vmix tke compute loop', .false.)

      err = 0
      call mpas_pool_get_array(statePool, 'tke', tkeNew, 2)
      call mpas_pool_get_array(statePool, 'tke', tkeOld, 1)
      call mpas_pool_get_array(statePool, 'ssh', ssh, 2)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 2)
      call mpas_pool_get_array(forcingPool, 'windStressZonal', windStressZonal)
      call mpas_pool_get_array(forcingPool, 'windStressMeridional', windStressMeridional)

      allocate(bTemp(nVertLevels+1), C(nVertLevels+1), rTemp(nVertLevels+1))

      call ocn_vmix_tke_compute_length(tkeOld,layerThickness,ssh,windStressZonal,windStressMeridional)
      call ocn_vmix_tke_compute_coeffs(tkeOld,.false.)

      sfcTKE = (15.8_RKIND*alpha)**(2.0_RKIND/3.0_RKIND)/(rho_sw*2.0_RKIND)

      nCells = nCellsOwned
      do iCell=1,nCells
         k=minLevelCell(iCell)
         sfcTKE = sfcTKE + alpha_conv*(0.4_RKIND*abs(surfaceBuoyancyForcing(iCell))*dThreshMLD(iCell))**(2.0/3.0)
         stressMag = sqrt(windStressZonal(iCell)**2.0 + windStressMeridional(iCell)**2.0)
         tkeNew(k,iCell) = max(emin0, sfcTKE*stressMag)

         k=maxLevelCell(iCell)+1
         tkeNew(k,iCell) = emin
      end do

      ! update RHS due to local terms (no transport)
      do iCell=1,nCells
         do k=minLevelCell(iCell)+1,maxLevelCell(iCell)
            shearProduction = vertViscTopOfCell(k,iCell)*shearSquared(k,iCell)
            buoyancyProduction = -vertDiffTopOfCell(k,iCell)*BruntvaisalaFreqTop(k,iCell)
            length = sqrt(ldo(k,iCell)*lup(k,iCell))
            dissipation = -ceps*tkeOld(k,iCell)**1.5_RKIND / length
            tkeNew(k,iCell) = max(emin, tkeOld(k,iCell) + dt*(shearProduction + buoyancyProduction + &
                                                              dissipation))
         end do
      end do

      !now transport of TKE

      do iCell=1,nCells
         Nsurf = minLevelCell(iCell)
         N = maxLevelCell(iCell)+1

         bTemp(Nsurf) = 1.0_RKIND
         rTemp(Nsurf) = tkeNew(Nsurf,iCell)

         vertViscTopOfCell(Nsurf,iCell) = ck*min(lup(Nsurf,iCell),ldo(Nsurf,iCell))*sqrt(tkeNew(Nsurf,iCell))
         do k=Nsurf+1,N-1
            nuMid = tkeEnhance*0.5_RKIND*(vertViscTopOfCell(k-1,iCell) + vertViscTopOfCell(k,iCell))
            A = -dt*nuMid/(0.5_RKIND*layerThickness(k-1,iCell)*(layerThickness(k-1,iCell)+layerThickness(k,iCell)))
            m = A/bTemp(k-1)
            nuMid = tkeEnhance*0.5_RKIND*(vertViscTopOfCell(k,iCell) + vertViscTopOfCell(k+1,iCell))
            C(k) = -dt*nuMid / (0.5_RKIND*layerThickness(k-1,iCell)*(layerThickness(k-1,iCell) + &
                                                                  layerThickness(k,iCell)))
            bTemp(k) = 1.0_RKIND - A - C(k) - m*C(k-1)

            rTemp(k) = tkeNew(k,iCell) - m*rTemp(k-1)
         end do

         vertViscTopOfCell(Nsurf,iCell) = 0.0_RKIND
         do k= N-1, Nsurf+1,-1
            tkeNew(k,iCell) = (rTemp(k) - C(k)*tkeNew(k+1,iCell)) / bTemp(k)
         end do
      end do

      do iCell=1,nCells
         sfcTKE = tkeNew(minLevelCell(iCell),iCell)
         do k=minLevelCell(iCell)+1,maxLevelCell(iCell)
            tkeNew(k,iCell) = tkeNew(k,iCell) + sfcTKE*niwFactor(k,iCell)
         end do
      end do

      call ocn_vmix_tke_compute_length(tkeNew,layerThickness,ssh,windStressZonal,windStressMeridional)
      call ocn_vmix_tke_compute_coeffs(tkeNew,.true.)

      call mpas_timer_stop('vmix tke compute loop')
   end subroutine ocn_vmix_tke_compute

end module ocn_vmix_tke


