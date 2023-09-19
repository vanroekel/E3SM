module agi_physics

  !-----------------------------------------------------------------!
  ! Module to compute the fraction of silver iodide available for  
  ! ice nucleation
  !
  !-----------------Code history------------------------------------!
  ! Authors: Luke Van Roekel
  !
  !-----------------------------------------------------------------!

  use physics_types,          only: physics_state, physics_ptend
  use physics_buffer, only: physics_buffer_desc
  use shr_kind_mod,  only: r8=>shr_kind_r8
  use shr_log_mod,   only: errMsg=>shr_log_errMsg
  use shr_const_mod, only: shr_const_rdair, shr_const_g
  use phys_control,  only: phys_getopts
  use physconst,     only: pi
  use cam_logfile,   only: iulog
  use constituents,   only: pcnst
  use spmd_utils,    only: masterproc
  use perf_mod,      only: t_startf, t_stopf
  use mpishorthand
  use cam_history_support, only: fillvalue
  use ppgrid,         only: pver, pcols

  implicit none

  private
  save

  ! ----------------- !
  ! Public interfaces !
  ! ----------------- !

  public :: agi_deposition, &
            agi_total_avail, &
            agi_condensation_freezing, &
            agi_contact_freezing, &
            agi_immersion_freezing, &
            agi_scavenging, &
            agi_immersed, & ! amount scavenged and not activated, done by tracking tendencies -- oThis one is confusing
            ! Meyers et al 1995 says this one is smaller than others and can maybe be ignored. Will implement last.
            agi_ccn

  contains

  subroutine agi_total_avail(state, pbuf, dt)

      use ppgrid, only: pver, pcols, begchunk, endchunk
      use cam_history, only: outfld
	  use constituents, only: cnst_get_ind
      use rgrid,          only: nlon
      use physconst,          only: tmelt, gravit, rga, rair, cpair, latvap, rearth, pi, cappa, epsilo
      use physics_types,          only: physics_state, physics_ptend, physics_ptend_init
      use physics_buffer, only: pbuf_old_tim_idx, physics_buffer_desc, pbuf_get_field, pbuf_get_index
      use pmgrid,           only: plev
      use wv_sat_methods, only: &
           qsat_water => wv_sat_qsat_water, &
           qsat_ice => wv_sat_qsat_ice

      type(physics_buffer_desc), pointer :: pbuf(:)
      type(physics_state), intent(inout) :: state
      real(r8), pointer, dimension(:,:) :: fAgi, agiAvail, agi, dissipation
      real(r8), dimension(pcols,pver) :: qvl, esl, sat_i, sat_w
      real(r8) :: qvi, esi, ftemp
	  real(r8), intent(in) :: dt
      real(r8), dimension(7,pcols,pver) :: fractionAgi ! seven is the number of independent processes for activation
      integer :: diss_idx, agi_indx, ix_agi, agiAvail_indx, fagi_indx
      integer :: ix_numliq, icol, k, ix_q
      fagi_indx = pbuf_get_index('fAgi')
	  call pbuf_get_field(pbuf, fagi_indx, fAgi)

      diss_idx = pbuf_get_index('dissipation')
      call pbuf_get_field(pbuf, diss_idx, dissipation)

      call cnst_get_ind('AgI', ix_agi)
	  call cnst_get_ind('availAgI', agiAvail_indx)
	  
	  ! Retrieve state%q index for water vapor specific humidity Q
	  !  Below, qvl is saturation specific humidity over water (kg/kg) and
	  !  qvi is saturation specific humidity over ice (kg/kg).
	  call cnst_get_ind('Q', ix_q)

      !loads the total available agi particles
      do icol = 1,pcols
         do k = 1,pver
            call qsat_water(state%t(icol,k), state%pmid(icol, k), esl(icol,k), qvl(icol,k))
            if(state%t(icol,k) > tmelt) then
               qvi = qvl(icol,k)
               esi = esl(icol,k)
            else
               call qsat_ice(state%t(icol,k), state%pmid(icol,k), esi, qvi)
            endif
            wv_mixing_ratio = state%q(icol,k,ix_q) / (1.0_r8 - state%q(icol,k,ix_q))  ! water vapor mixing ratio = Q/(1-Q), where Q is water vapor specific humidity (Bohren & Albrecht book 1998 edition Eq. 5.10)
            wv_pressure = ( wv_mixing_ratio / (wv_mixing_ratio + epsilo) ) * state%pmid(icol, k)  ! water vapor pressure (Wallace & Hobbs book 1977 edition Eq. 2.61)
            sat_w(icol,k) = wv_pressure/esl(icol,k) ! water vapor saturation ratio wrt water ("water saturation ratio"):  Bohren & Albrecht book 1998 edition Eq. 5.11
            sat_i(icol,k) = wv_pressure/esi         ! water vapor saturation ratio wrt ice ("ice saturation ratio"):  :  Bohren & Albrecht book 1998 edition Eq. 5.11

            agi(icol,k) = state%q(icol,k,ix_agi)
         end do
      end do

      call cnst_get_ind('NUMLIQ', ix_numliq)

      call agi_scavenging(state%t, state%pmid, dissipation, sat_w, state%q(:,:,ix_numliq), &
	                      agi, esl, qvl, fractionAgi, 6, dt) 
      call agi_deposition(state%t,sat_i,fractionAgi,1)
      call agi_condensation_freezing(state%t,sat_w,fractionAgi,2)
      call agi_contact_freezing(sat_i,state%t,fractionAgi,3)
      call agi_immersion_freezing(state%t,fractionAgi,4)
      call agi_ccn(sat_w,fractionAgi,5)
  !    call agi_ctf(sat_i, state%t, fractionAgi, 7)

      do icol = 1,pcols
	     do k = 1, pver
		    ftemp = fractionAgi(1,icol,k) + fractionAgi(2,icol,k) + fractionAgi(3,icol,k) + &
			        fractionAgi(4,icol,k) + fractionAgi(5,icol,k) !+ fractionAgi(7,icol,k)
			fAgi(icol,k) = min(1.0_r8, ftemp)
			state%q(icol,k,agiAvail_indx) = fAgi(icol,k)*agi(icol,k)
	     end do
	  end do

  end subroutine agi_total_avail

  subroutine agi_contact_freezing(sat_i, T, fAgi, spot)

    use ppgrid, only: pver, pcols
	integer, intent(in) :: spot
    real(r8),dimension(:,:,:), intent(out) :: fAgi
	real(r8),dimension(pcols,pver), intent(in) :: sat_i, T

	real(r8), parameter:: a = 0.0878_r8, &
	                      b = -3.7947_r8, &
						  c = 52.3167_r8, &
						  d = -255.4484_r8, &
						  e = 568.3257_r8, &
						  f = -460.4234_r8, &
						  g = -63.1248_r8

    integer :: icol, k

    real(r8) :: sm1

    do icol = 1, pcols
       do k = 1, pver
	      fAgi(spot,icol,k) = 0.0_r8
	      if(sat_i(icol,k) .gt. 1.058_r8 .and. T(icol,k) .lt. 269.2_r8) then
	 	     sm1 = sat_i(icol,k) - 1.0_r8
	         fAgi(spot,icol,k) = fAgi(6,icol,k)*(a+b*sm1 + c*sm1**2._r8 + d*sm1**3._r8 + e*sm1**4._r8 + &
		  	 					 				 f*sm1**5._r8 + g*sm1**6._r8)
          end if
	   end do
    end do

  end subroutine agi_contact_freezing

  subroutine agi_ccn(sat_w, fAgi, spot)
    use ppgrid, only: pver, pcols

    real(r8),dimension(pcols,pver) :: sat_w 
    integer, intent(in) :: spot
    real(r8),dimension(:,:,:), intent(out) :: fAgi
    integer :: icol,k

    do icol = 1,pcols
       do k=1,pver
          if(sat_w(icol,k) <1.05_r8 .and. sat_w(icol,k) .ge. 1.0_r8) then
             fAgi(spot,icol,k) = 5.0_r8*(sat_w(icol,k) - 1._r8)**1.5_r8
          else
             fAgi(spot,icol,k) = 0.0_r8
          end if
       end do
    end do

  end subroutine agi_ccn

  subroutine agi_scavenging(T, P, dissipation, sat_w, Nc, Nagi, esl, &
                            qsat, fscavAgi, spot, dt)

    use physconst, only: rair, rh2o, latvap, rhoh2o, boltz, pi, tmelt
    use ppgrid, only: pver, pcols

	real(r8),dimension(pcols,pver) :: qsat, sat_w, esl
    real(r8), dimension(:,:), intent(in) :: T, P, dissipation, Nc, Nagi
    real(r8), intent(in) :: dt
    real(r8), dimension(:,:,:), intent(out) :: fscavAgi

    integer :: icol, k, m, i, spot
    !need a bunch of USE statements
    real(r8), parameter :: &
      lambdaAO = 6.6E-8_r8,   & !AgI mean free path [m]
      Dv = 2.2E-5_r8,         & !vapor diffusivity (estimate) [m^2/s]
      Ka = 0.023_r8,          & !thermal conductivity of air [W/(mK)]
      Fv = 1.0_r8,            & !ventilation of vapor from scavenging param:  Fv is unitless, assumes fall velocity of cloud droplets = 0, and comes from Vuckovic et al. 2022 Sec. 3.1.1.c (thermophoresis)
      Fh = 1.0_r8,            & !ventilation of heat from scavenging param:  Fh is unitless, assumes fall velocity of cloud droplets = 0, and comes from Vuckovic et al. 2022 Sec. 3.1.1.d (diffusiophoresis)
      Asl = 1.257_r8,         & !factor in alpha equation
      Bsl = 0.4_r8,           & !factor in alpha equation
      Csl = 1.1_r8,           & !factor in alpha equation
      Po = 1013.25_r8,        & !reference pressure in param of final lambda [hPa]
      Tref = 293.15_r8,       & !reference temperature in param of final lambda [K]
      Dagi = 0.04E-6_r8,      & !geometric mean diameter of Agi -- taken from Xue et al 2013
      sigAgi = 2.0_r8           !geometric standard deviation of AgI -- from Xue et al 2013

      real(r8) :: No, lam, mu, mu_a, nu_a, tdegC, rho, Cb0, Cb3, Ctb0, Ctb3
      real(r8) :: Cthdf, Nconctend, factor, lama, Nkn, alpha, phiTH, A3, P_hPa, cscf

	  real(r8),dimension(9) :: Map, Mc  !moments of aerosol distributions

      do icol = 1, pcols
	     do k = 1, pver
		    rho = P(icol,k) / (rair*T(icol,k))
		    !air dynamic viscosity -- taken from Vuckovic et al 2022 in Atmospheric research (bottom page 3)
			tdegC = T(icol,k) - tmelt
		    if(tdegC .ge. 0.0_r8) then
			   mu_a = (1.718_r8 + 0.0049_r8*tdegC)*1e-3_r8 ! it is 1e-3 due to conversion from Poise to SI units
			else
			   mu_a = (1.718_r8 + 0.0049_r8*tdegC - 1.2E-5_r8 * tdegC**2.0_r8) * 1e-3_r8 !using 1e-3 due to conversion from poise to SI units
			endif
			nu_a = mu_a / rho ! kinematic viscosity is just dynamic visc / rho

            i = 1
            do m=-2,6
               !Moments of aerosol particles (Map) are from Caro et al (2004) Eq C.3
			   Map(i) = Nagi(icol,k) * Dagi**m * exp(m**2._r8 * log(sigAgi)**2._r8 / 2.0_r8)
               !Mass of cloud particle is from Morrison and Gettleman 2008
               !mu is acutually eta MG08 Eq2)
			   mu = 0.0005714_r8*Nc(icol,k) + 0.2714_r8
               !lam is MG08 Eq 3
			   lam = ((pi * rhoh2o * Nc(icol,k) * gamma(mu+4._r8))/(6.0_r8*qsat(icol,k) * gamma(mu+1._r8)))**(1.0_r8/3.0_r8)
               !No is Eq 4 of MG08
			   No = Nc(icol,k) * lam**(mu+1._r8)/gamma(mu+1._r8)
               !Mc is integration of these terms -- consult your favorite table of integrals
			   Mc(i) = No * gamma(m+mu+1._r8) / lam**(m+mu)
			   i = i+1
            end do

            ! Now compute three coefficients:  Cb*, Ctb*, and Cthdf following Caro et al. equation C.9, where
            !   the * is a placeholder for m=0 (number concentration) or m=3 (mass) to indicate the moment.
            
            !this term is Caro et al 2004 Equation c.9, but since it is a function of 'm' create two 
            !constants of it
            factor = 2.0_r8 * boltz * T(icol,k) / (3.0_r8 * mu_a)
			Cb0 = factor/Map(3)
			Cb3 = factor/Map(6)
            !Ctb is also equation c.9 of Caro et al 2004
			!note that in 'factor' below, Caro et al 2004 suggests eps_c (dissipation) could be set constant
			!should be around 46.2 cm2/s-3 for convective clouds 
			factor = 3.0_r8 * pi / (2.0_r8 * sqrt(15.0_r8)) * (dissipation(icol,k) / nu_a)**0.5_r8
			Ctb0 = factor/Map(3)
			Ctb3 = factor/Map(6)

            ! Now compute constants C associated with thermophoresis and diffusiophoresis (Caro et al. equation C.9).
            !   To compute Cthdf, first compute phiTH (Caro et al. equation 33) and A3 (Caro et al. equation 10), plus a few dependencies.
            !Lama is defined in the first sentence of the first paragraph of Vuckovic et al 2022 page 4
            !just above eq 6
            P_hPa = P(icol,k) * 0.01_r8   ! compute air pressure in hPa (P is in units Pa)
            lama = lambdaAO * (Po / P_hPa) * (T(icol,k) / Tref)
            !Nkn is defined in the same place as lama
			Nkn = 2.0_r8 * lama / Dagi
            !alpha also in same place as lama
			alpha = Asl + Bsl * exp(-1._r8 * Csl / Nkn)
            !phiTH is equation 33 of Caro et al 2004
			phiTH = ((1.0_r8 + alpha * Nkn) * Nkn) / ((1.0_r8 + 3.0_r8 * Nkn) + (1.0_r8 + 5.0_r8 * Nkn))   ! unitless

            !A3 is equation 10 of Caro et al 2004
            A3 = (rhoh2o*rh2o*T(icol,k)/(esl(icol,k)*Dv) + rhoh2o*latvap/(Ka*T(icol,k))*(latvap/(rh2o*T(icol,k)) - 1.0_r8))**(-1.0_r8)
            !Equation C.9 of Caro et al 2004
            Cthdf = (-(2.0_r8*pi*phiTH*Fh*latvap*rhoh2o)/(P(icol,k)) + (2.4_r8*pi*rhoh2o*Fv)/rho)*sat_w(icol,k)*A3  ! Use P in units Pa here;  Cthdf is unitless

            !Equation C.8 of Caro et al 2004 -- Note I'm applying this as a temnporary tendency as we need a fraction
            !  Here, we assume a single value of j (aerosol mode), likely j=0 since AgI is mostly in the nuclei mode according to Xue et al. (2013) page 1436.
			! First, explicitly define Cunningham slip correction factor (cscf).  Based on Vuckvic et al. (2022),
			!   cscf = 1 + alpha*Nkn    ! (using variable names in this code)
			!   Note: Caro et al. (2004) uses "alpha" to denote the entire Cunningham slip correction factor (not the -part- of the CSCF as in Vuckovic et al. 2022),
			!         so I will explicitly define cscf to avoid more confusion.
			cscf = 1._r8 + alpha*Nkn    ! From Vuckovic et al. 2022 page 4
			Nconctend = Nagi(icol,k) + dt*(Cb0*((Map(3) + 2.0_r8*cscf*lama*Map(2))*Mc(3) + Map(3)*(Mc(3) + 2.0_r8*cscf*lama*Mc(2)) + &
			                    Mc(4)*(Map(2) + 2.0_r8*cscf*lama*Map(1)) + Map(4)*(Mc(2) + 2.0_r8*cscf*lama*Mc(1))) + &
								Ctb0*(Map(6)*Mc(3) + 3.0_r8*Map(5)*Mc(4) + 3.0_r8*Map(4)*Mc(5) + Map(3)*Mc(6)) + &
								Cthdf*Map(3)*Mc(4))

!			Qrattend(icol,k) =  Cb3*((Map(6) + 2.0_r8*lama*Map(5))*Mc(3) + Map(6)*(Mc(3) + 2.0_r8*lama*Mc(2)) + &
!			                    Mc(4)*(Map(5) + 2.0_r8*lama*Map(4)) + Map(7)*(Mc(2) + 2.0_r8*lama*Mc(1))) + &
!								Ctb0*(Map(9)*Mc(3) + 3.0_r8*Map(8)*Mc(4) + 3.0_r8*Map(7)*Mc(5) + Map(6)*Mc(6))

            !fscav is the change in Agi normalized by original is the scavenging fractionation
            fscavAgi(spot,icol,k) = (Nagi(icol,k) - Nconctend) / Nagi(icol,k)
         end do
      end do
  end subroutine agi_scavenging

  subroutine agi_immersion_freezing(temp, fAgi, spot)

    use ppgrid, only: pver, pcols

    real(r8),dimension(pcols,pver) :: temp
    integer, intent(in) :: spot
    real(r8),dimension(:,:,:), intent(out) :: fAgi

    fAgi(spot,:,:) = 0.0_r8
     ! Meyers et al 1995 suggests that this term is much smaller than the other modes so we leave it out for now
  end subroutine agi_immersion_freezing

  subroutine agi_immersed()
     ! this will be a helper routine to calculate the number of drops immersed in drops that are unactivated
  end subroutine agi_immersed

  subroutine agi_condensation_freezing(temp, sat_w, fAgi, spot)

    use ppgrid, only: pver, pcols

    real(r8),dimension(pcols,pver) :: temp, sat_w
    integer, intent(in) :: spot
    real(r8),dimension(:,:,:), intent(out) :: fAgi

    real(r8), parameter :: &
       To = 10.0_r8,    &
       a = 900.0_r8

    integer :: icol, k

    do icol = 1,pcols
       do k = 1,pver
          if(temp(icol,k) < 268.66_r8 .and. sat_w(icol,k) > 1.0_r8) then
             fAgi(spot,icol,k) = a*((268.66_r8 - temp(icol,k))/To)**3.0_r8*(sat_w(icol,k) - 1.0_r8)**2.0_r8
          else
             fAgi(spot,icol,k) = 0.0_r8
          end if

       end do
    end do

  end subroutine agi_condensation_freezing


  subroutine agi_deposition(temp, sat_i, fAgi, spot)

    use ppgrid, only: pver, pcols

    real(r8),dimension(pcols,pver) :: temp, sat_i
    integer, intent(in) :: spot
    real(r8),dimension(:,:,:), intent(out) :: fAgi

    real(r8), parameter :: &
       To = 10.0_r8,    &
       a = -3.25E-3_r8, &
       b = 5.39E-5_r8,  &
       c = 4.35E-2_r8,  &
       d = 1.55E-4_r8,  &
       e = -0.07_r8

    integer :: icol, k

    do icol = 1,pcols
       do k = 1,pver

          if(temp(icol,k) < 268.2_r8 .and. sat_i(icol,k) > 1.04_r8) then
             fAgi(spot,icol,k) = a*(sat_i(icol,k) - 1) + b*((273.16_r8 - temp(icol,k)) / To) + &
                                 c*(sat_i(icol,k) - 1)**2.0 + d*((273.16_r8 - temp(icol,k)) / To) + &
                                 e*(sat_i(icol,k) - 1)**3.0
          else
            fAgi(spot,icol,k) = 0.0_r8
          end if
       end do
    end do

   end subroutine agi_deposition

end module agi_physics
