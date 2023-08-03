#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <cstddef>
#include <Kokkos_Core.hpp>
#include "cvmix_global_const.h"
#include "cvmix_kokkos_defs.hpp"
#ifdef GPU
    #include <cuda.h>
#endif

using namespace std;


extern"C"{
struct CVMix_mixingOpts{
    double bl1, bl2, bl3, bl4, prandtl;
    double PP_nu_zero, PP_alpha, PP_exp, PP_nu_b, PP_kappa_b, KPP_nu_zero, KPP_Ri_zero, KPP_exp;
    double convect_diff, convect_visc, BVsqr_convect;
    bool bryanLewis, lBruntVaisala, lShearKPP, lShearPP,
         lEkman, lMonOb, lKPP, computeLangmuir;
    bool LANGMUIR_LWF16, LANGMUIR_lifk17, parabolicNonLocal;
    double Ri_crit, surf_layer_ext, c_s, minVtsqr, rho_o, minOSBLUnderIce;
    double a_s, c_m, a_m, non_local_coeff, numRiSmoothPasses;
} ;

}
void compute_mixing(const int, const int, double *, double *, double *, double *,
            double *, double *, double *, double *, double *, double *,
            double *, double *, int *, double *, double *, double *, double *, CVMix_mixingOpts *);

double compute_kpp_ustokes_SL_model(double, double);

void compute_kpp_EFactor_model(double *, double *, double *, double *, int);
extern"C"{

void c_kokkos_finalize( ) { 
   delete(KPPvars);
   Kokkos::finalize();
}


void c_kokkos_initialize(int nVertLevels, int nCols) {

    #ifdef GPU
    // Initialize Host mirror device
        Kokkos::HostSpace::execution_space::initialize(1);
        const unsigned device_count = Kokkos::Cuda::detect_device_count();
    //Use the last device:
        Kokkos::Cuda::initialize( Kokkos::Cuda::SelectDevice(device_count-1) );
    #else
        Kokkos::initialize();
    #endif


    KPPvars = new KPPvariables();

    KPPvars->kv = ViewDoubleType("viscosity", nVertLevels+1, nCols);
    KPPvars->kh = ViewDoubleType("diffusivity", nVertLevels+1, nCols);
    KPPvars->nlt = ViewDoubleType("non-local Temperature", nVertLevels+1, nCols);
    KPPvars->lt =  ViewDoubleType("layer Thickness", nVertLevels, nCols);
    KPPvars->u =   ViewDoubleType("U-velocity", nVertLevels, nCols);
    KPPvars->v =   ViewDoubleType("V-velocity", nVertLevels, nCols);
    KPPvars->f =   ViewColType("Coriolis Parameter", nCols);
    KPPvars->n2 =  ViewDoubleType("brunt Vaisala frequency", nVertLevels+1, nCols);
    KPPvars->rho = ViewDoubleType("density", nVertLevels, nCols);
    KPPvars->ri =  ViewDoubleType("gradient richardson number", nVertLevels+1, nCols);
    KPPvars->sfc_buoy = ViewColType("surface buoyancy forcing", nCols);
    KPPvars->us =  ViewColType("surface friction velocity", nCols);
    KPPvars->maxLevs = ViewIntColType("max active level", nCols);
    KPPvars->zw =  ViewDoubleType("interface depths",nVertLevels+1,nCols);
    KPPvars->zt =  ViewDoubleType("mid cell depths",nVertLevels,nCols);
    KPPvars->dzw=  ViewDoubleType("distance between layer interfaces",nVertLevels+1,nCols);
    KPPvars->maskt = ViewDoubleType("cell center mask",nVertLevels,nCols);
    KPPvars->maskw = ViewDoubleType("cell interface mask",nVertLevels+1,nCols);
    KPPvars->h = ViewColType("boundary layer depth", nCols);
    KPPvars->interfaceForcings = ViewDoubleType("interfaceForcing for bld calculation", nVertLevels,nCols);
    KPPvars->bldArray = ViewDoubleType("holding for boundary layer depth temporary", nVertLevels,nCols);
    KPPvars->bulkRi = ViewDoubleType("bulk richardson number for boundary layer depth", nVertLevels,nCols);
    KPPvars->ws = ViewDoubleType("turbulent scalar velocity scale", nVertLevels,nCols);
    KPPvars->wm = ViewDoubleType("turbulent momentum velocity scale", nVertLevels,nCols);
    KPPvars->sfcAverageIndex = ViewDoubleType("surface area average index (vertical)", nVertLevels,nCols);
    KPPvars->uVelocitySum = ViewDoubleType("U velocity sum in vertical (KPP)", nVertLevels,nCols);
    KPPvars->vVelocitySum = ViewDoubleType("V veloctiy sum in vertical (KPP)", nVertLevels,nCols);
    KPPvars->potentialDensitySum = ViewDoubleType("pot. den. sum in vertical for KPP", nVertLevels,nCols);
    KPPvars->bulkRichardsonNumberBuoy = ViewDoubleType("Ri_bulk buoy component for KPP", nVertLevels,nCols);
    KPPvars->bulkRichardsonNumberShear = ViewDoubleType("Ri_bulk shear component for KPP", nVertLevels,nCols);
    KPPvars->N_cntr = ViewDoubleType("sqrt of bruntVaisala", nVertLevels,nCols);
    KPPvars->LaSL = ViewColType("surface layer average Langmuir niumber", nCols);
    KPPvars->Langmuir_EFactor = ViewColType("langmuir enhancement factor",nCols);
    KPPvars->Minv = View3DType("array for solving coefficients", 3,3,nCols);
    KPPvars->Vt2 = ViewDoubleType("turbulent velocity shear", nVertLevels,nCols);
    KPPvars->rhs = ViewDoubleType("rhs for coefficient solve", 3,nCols);
    KPPvars->coeffs = ViewDoubleType("coefficient list for polynomial fit", 3,nCols);
    KPPvars->OBL_Mdiff = ViewDoubleType("KPP momentum diffusivities",nVertLevels,nCols);
    KPPvars->OBL_Tdiff = ViewDoubleType("KPP scalar diffusivities",nVertLevels,nCols);
    KPPvars->intOBL = ViewColType("integer position of OBL depth", nCols);
    KPPvars->riSmoothed = ViewDoubleType("smoothed gradient richardson number", nVertLevels, nCols);
    KPPvars->iceFrac = ViewColType("ice fraction from sea ice model", nCols);

    KPPvars->h_iceFrac = Kokkos::create_mirror_view(KPPvars->iceFrac);
    KPPvars->h_efactor = Kokkos::create_mirror_view(KPPvars->Langmuir_EFactor);
    KPPvars->h_lasl = Kokkos::create_mirror_view(KPPvars->LaSL);
    KPPvars->h_OSBL = Kokkos::create_mirror_view(KPPvars->h);
    KPPvars->h_kv = Kokkos::create_mirror_view(KPPvars->kv);
    KPPvars->h_kh = Kokkos::create_mirror_view(KPPvars->kh);
    KPPvars->h_nlt = Kokkos::create_mirror_view(KPPvars->nlt);
    KPPvars->h_rib = Kokkos::create_mirror_view(KPPvars->bulkRi);
    KPPvars->h_lt = Kokkos::create_mirror_view(KPPvars->lt);
    KPPvars->h_u = Kokkos::create_mirror_view(KPPvars->u);
    KPPvars->h_v = Kokkos::create_mirror_view(KPPvars->v);
    KPPvars->h_n2 = Kokkos::create_mirror_view(KPPvars->n2);
    KPPvars->h_ri = Kokkos::create_mirror_view(KPPvars->ri);
    KPPvars->h_rho = Kokkos::create_mirror_view(KPPvars->rho);
    KPPvars->h_f = Kokkos::create_mirror_view(KPPvars->f);
    KPPvars->h_sfc_buoy = Kokkos::create_mirror_view(KPPvars->sfc_buoy);
    KPPvars->h_us = Kokkos::create_mirror_view(KPPvars->us);
    KPPvars->h_h = Kokkos::create_mirror_view(KPPvars->h);
    KPPvars->h_maxLevs = Kokkos::create_mirror_view(KPPvars->maxLevs);

}

void c_compute_cvmix2_mixing(const int nCols, const int nVertLevels,  double *visc, double *diff,
    double *layerThickness, double *U, double *V, double *NLT, double *N2, double *RI_grad,
    double *density, double *fCor, double *sfcBuoyancyForcing, double *sfcFrictionVelocity,
    int *maxLevel, double *OSBL, double *iceFraction, double *lasl, double *wind10,
    CVMix_mixingOpts* mixingOpts)
{

  double efactor[nCols];

// from a logical compute LaSL and La_efactor
if(mixingOpts->computeLangmuir){
    for(int iCol = 0; iCol < nCols; iCol++){
        double us = compute_kpp_ustokes_SL_model(wind10[iCol],OSBL[iCol]);
        lasl[iCol] = sqrt(sfcFrictionVelocity[iCol]/us);
    }
    compute_kpp_EFactor_model(wind10, sfcFrictionVelocity,
                                OSBL, efactor, nCols);
} else{
    for(int iCol = 0; iCol < nCols; iCol++){
      lasl[iCol] = 1.0E10;
      efactor[iCol] = 1.0;
    }
}

compute_mixing(nCols,nVertLevels,visc,diff,NLT,layerThickness,U,V,N2,density,
            RI_grad, fCor, sfcBuoyancyForcing,sfcFrictionVelocity,maxLevel,OSBL,
            iceFraction, lasl, efactor,mixingOpts);

}
}

void compute_mixing(const int nCols, const int nVertLevels, double *visc, double *diff,
            double *NLT, double *layerThickness, double *U, double *V,
            double *N2, double *density, double *RI_grad, double *fCor, double *sfcBuoyancyForcing,
            double *sfcFrictionVelocity, int *maxLevel, double *OSBL, double *iceFraction, double *lasl,
            double *efactor, CVMix_mixingOpts* mixingOpts)
{

    auto& uVelocitySum = KPPvars->uVelocitySum;
    auto& vVelocitySum = KPPvars->vVelocitySum;
    auto& potentialDensitySum = KPPvars->potentialDensitySum;
    auto& bulkRichardsonNumberBuoy = KPPvars->bulkRichardsonNumberBuoy;
    auto& bulkRichardsonNumberShear = KPPvars->bulkRichardsonNumberShear;
    auto& N_cntr = KPPvars->N_cntr;
    auto& stokesDrift = KPPvars->stokesDrift;
    auto& LaSL = KPPvars->LaSL;
    auto& bulkRi = KPPvars->bulkRi;
    auto& Vt2 = KPPvars->Vt2;
    auto& Minv = KPPvars->Minv;
    auto& h = KPPvars->h;
    auto& rhs = KPPvars->rhs;
    auto& coeffs = KPPvars->coeffs;
    auto& intOBL = KPPvars->intOBL;
    auto& surfaceAverageIndex = KPPvars->sfcAverageIndex;
    auto& kh = KPPvars->kh;
    auto& zt = KPPvars->zt;
    auto& zw = KPPvars->zw;
    auto& kv = KPPvars->kv;
    auto& dzw = KPPvars->dzw;
    auto& n2 = KPPvars->n2;
    auto& ri = KPPvars->ri;
    auto& lt = KPPvars->lt;
    auto& maskw = KPPvars->maskw;
    auto& maskt = KPPvars->maskt;
    auto& rho = KPPvars->rho;
    auto& maxLevs = KPPvars->maxLevs;
    auto& nlt = KPPvars->nlt;
    auto& ws = KPPvars->ws;
    auto& wm = KPPvars->wm;
    auto& OBL_Tdiff = KPPvars->OBL_Tdiff;
    auto& OBL_Mdiff = KPPvars->OBL_Mdiff;
    auto& bldArray = KPPvars->bldArray;
    auto& sfc_buoy = KPPvars->sfc_buoy;
    auto& us = KPPvars->us;
    auto& u = KPPvars->u;
    auto& v = KPPvars->v;
    auto& f = KPPvars->f;
    auto& interfaceForcings = KPPvars->interfaceForcings;
    auto& Langmuir_EFactor = KPPvars->Langmuir_EFactor;
    auto& riSmoothed = KPPvars->riSmoothed;
    auto& iceFrac = KPPvars->iceFrac;

    auto& h_nlt = KPPvars->h_nlt;
    auto& h_rib = KPPvars->h_rib;
    auto& h_lt = KPPvars->h_lt;
    auto& h_u = KPPvars->h_u;
    auto& h_v = KPPvars->h_v;
    auto& h_rho = KPPvars->h_rho;
    auto& h_ri = KPPvars->h_ri;
    auto& h_n2 = KPPvars->h_n2;
    auto& h_kv = KPPvars->h_kv;
    auto& h_kh = KPPvars->h_kh;
    auto& h_OSBL = KPPvars->h_OSBL;
    auto& h_sfc_buoy = KPPvars->h_sfc_buoy;
    auto& h_us = KPPvars->h_us;
    auto& h_maxLevs = KPPvars->h_maxLevs;
    auto& h_f = KPPvars->h_f;
    auto& h_lasl = KPPvars->h_lasl;
    auto& h_efactor = KPPvars->h_efactor;
    auto& h_iceFrac = KPPvars->h_iceFrac;

    for(int iCol=0; iCol < nCols; iCol++){
        for(int i=0; i<nVertLevels; i++){
            int ii = i+iCol*nVertLevels;
            int k = ii / nVertLevels;
            int j = ii - k*nVertLevels;

            h_lt(j,k) = (layerThickness)[ii];
            h_u(j,k) = (U)[ii];
            h_v(j,k) = (V)[ii];
            h_rho(j,k) = (density)[ii];
        }

        for(int i=0; i<(nVertLevels+1); i++){
            int ii = i+iCol*(nVertLevels+1);
            int k = ii / (nVertLevels+1);
            int j = ii - k*(nVertLevels+1);

            h_ri(j,k) = max(cvmix_zero,(RI_grad)[ii]);
            h_n2(j,k) = (N2)[ii];
            h_kv(j,k) = visc[ii];
            h_kh(j,k) = diff[ii];
        }

        h_OSBL(iCol) = (OSBL)[iCol];
        h_sfc_buoy(iCol) = (sfcBuoyancyForcing)[iCol];
        h_us(iCol) = (sfcFrictionVelocity)[iCol];
        h_maxLevs(iCol) = (maxLevel)[iCol];
        h_f(iCol) = fCor[iCol];
        h_lasl(iCol) = lasl[iCol];
        h_efactor(iCol) = efactor[iCol];
        h_iceFrac(iCol) = iceFraction[iCol];
//host views of eFactor lasl wind10

    }

    Kokkos::deep_copy(iceFrac,h_iceFrac);
    Kokkos::deep_copy(lt,h_lt);
    Kokkos::deep_copy(kv,h_kv);
    Kokkos::deep_copy(kh,h_kh);
    Kokkos::deep_copy(u,h_u);
    Kokkos::deep_copy(v,h_v);
    Kokkos::deep_copy(n2,h_n2);
    Kokkos::deep_copy(ri,h_ri);
    Kokkos::deep_copy(rho,h_rho);
    Kokkos::deep_copy(h,h_OSBL);
    Kokkos::deep_copy(f,h_f);
    Kokkos::deep_copy(sfc_buoy,h_sfc_buoy);
    Kokkos::deep_copy(us,h_us);
    Kokkos::deep_copy(maxLevs,h_maxLevs);
    Kokkos::deep_copy(LaSL,h_lasl);
    Kokkos::deep_copy(Langmuir_EFactor,h_efactor);

    auto bl1 = mixingOpts->bl1;
    auto bl2 = mixingOpts->bl2;
    auto bl3 = mixingOpts->bl3;
    auto bl4 = mixingOpts->bl4;
    auto prandtl = mixingOpts->prandtl;

    auto nu_zero = mixingOpts->PP_nu_zero;
    auto loc_exp = mixingOpts->PP_exp;
    auto PP_alpha = mixingOpts->PP_alpha;
    auto PP_nu_b = mixingOpts->PP_nu_b;
    auto PP_kappa_b = mixingOpts->PP_kappa_b;

    auto kpp_nu_zero = mixingOpts->KPP_nu_zero;
    auto kpp_loc_exp = mixingOpts->KPP_exp;
    auto KPP_Ri_zero = mixingOpts->KPP_Ri_zero;
    auto bryanLewis = mixingOpts->bryanLewis;
    auto lShearPP = mixingOpts->lShearPP;
    auto lShearKPP = mixingOpts->lShearKPP;
    auto lBruntVaisala = mixingOpts->lBruntVaisala;
    auto lEkman = mixingOpts->lEkman;
    auto lMonOb = mixingOpts->lMonOb;
    auto lKPP = mixingOpts->lKPP;
    auto LANGMUIR_lifk17 = mixingOpts->LANGMUIR_lifk17;
    auto LANGMUIR_LWF16 = mixingOpts->LANGMUIR_LWF16;

    auto convect_diff = mixingOpts->convect_diff;
    auto convect_visc = mixingOpts->convect_visc;
    auto cvmix_zeroD = 0.0;
    auto cvmix_oneD = 1.0;
    auto cvmix_p5D = 0.5;
    auto cvmix_piD = 3.14159265358979323846;
    auto kappa = 0.4;
    auto gravity = 9.8061;
    auto beta = 0.2;
    auto CVmin = 1.7;
    auto numRiSmoothPasses = mixingOpts->numRiSmoothPasses;
    auto minOSBLUnderIce = mixingOpts->minOSBLUnderIce;

    // Set values for shape functions
// Parabolic non-local will require scaling of non-local still -- TODO

    double Mshape1 = cvmix_zeroD;
    double Mshape2 = cvmix_oneD;
    double Mshape3 = -2.0;
    double Mshape4 = cvmix_oneD;
    double Tshape1 = Mshape1;
    double Tshape2 = Mshape2;
    double Tshape3 = Mshape3;
    double Tshape4 = Mshape4;
    double TshapeNL1 = Tshape1;
    double TshapeNL2 = Tshape2;
    double TshapeNL3 = Tshape3;
    double TshapeNL4 = Tshape4;

    if(mixingOpts->parabolicNonLocal){
      TshapeNL1 = cvmix_oneD;
      TshapeNL2 = -2.0;
      TshapeNL3 = cvmix_oneD;
      TshapeNL4 = cvmix_zeroD;
    }

    auto c_s = mixingOpts->c_s;
    auto surf_layer_ext = mixingOpts->surf_layer_ext;
    auto a_s = mixingOpts->a_s;
    auto a_m = mixingOpts->a_m;
    auto c_m = mixingOpts->c_m;
    auto minVtsqr = mixingOpts->minVtsqr;
    auto rho_o = mixingOpts->rho_o;
    auto Ri_crit = mixingOpts->Ri_crit;
    auto non_local_coeff = mixingOpts->non_local_coeff*kappa*pow(kappa*surf_layer_ext*c_s,1./3.);

    Kokkos::parallel_for(nCols, KOKKOS_LAMBDA(int iCol){

        for(int k=0; k < maxLevs(iCol); k++){
            maskt(k,iCol) = cvmix_oneD;
            maskw(k,iCol) = cvmix_oneD;
        }

        for(int k=maxLevs(iCol); k < nVertLevels; k++){
            maskt(k,iCol) = cvmix_zeroD;
            maskw(k,iCol) = cvmix_zeroD;
        }
        maskw(nVertLevels,iCol) = cvmix_zeroD;

// Build depth arrays
        zw(0,iCol) = cvmix_zeroD;
        dzw(0,iCol) = cvmix_p5D*lt(0,iCol);
        zt(0,iCol) = -cvmix_p5D*lt(0,iCol);

        for(int k=1; k<nVertLevels; k++){
            zw (k,iCol) = zw(k-1,iCol) - lt(k-1,iCol);
            zt (k,iCol) = zw(k  ,iCol) - cvmix_p5D*lt(k,iCol);
            dzw(k,iCol) = zt(k-1,iCol) - zt(k,iCol);
        }

        double zbot = zw(nVertLevels-1,iCol) - lt(nVertLevels-1,iCol);
        zw(nVertLevels,iCol) = zbot;
        dzw(nVertLevels,iCol) = zt(nVertLevels-1,iCol) - zw(nVertLevels,iCol);

        if(lShearPP || lShearKPP){
              for(int iSmooth=0; iSmooth < numRiSmoothPasses; iSmooth++){
                for(int k=2; k < nVertLevels; k++){
                  riSmoothed(k,iCol) = (ri(k-1,iCol) + 2.0*ri(k,iCol) + ri(k+1,iCol)) / 4.0;
               }
                for(int k=2; k< nVertLevels; k++){
                  ri(k,iCol) = riSmoothed(k,iCol);
                }
              }
        }

        // Add Background Mixing Bryan Lewis
        if(bryanLewis){
            for(int k=1; k < nVertLevels+1; k++){
                kh(k,iCol) = bl1 + (bl2/cvmix_piD)*atan(bl3*(zw(k,iCol)-bl4));
                kv(k,iCol) = prandtl*kh(k,iCol);
            }
        }

        // Add PP81 shear instability driven mixing
        if(lShearPP){
            for(int k=1; k<nVertLevels+1; k++){
                double denom = cvmix_oneD + PP_alpha*ri(k,iCol);
                double khTemp = (nu_zero / pow(denom,loc_exp) + PP_nu_b)*maskw(k,iCol);
                kv(k,iCol) = kv(k,iCol)+khTemp;
                kh(k,iCol) = kh(k,iCol) + khTemp / denom + PP_kappa_b;

            }
        }

        // Add Large et al. 1994 Shear instability driven mixing
        if(lShearKPP){
            for(int k=1; k<nVertLevels+1; k++){
                double khTemp = kpp_nu_zero*pow(cvmix_oneD -
                        pow((min(ri(k,iCol),KPP_Ri_zero) /KPP_Ri_zero ),2),kpp_loc_exp);
                kh(k,iCol) = kh(k,iCol) + khTemp;
                kv(k,iCol) = kv(k,iCol) + khTemp;
            }
        }

        // Tidal


        // double diffusion

        if(lKPP){


        int topIndex = 0;
        for(int k=0; k < nVertLevels; k++){
          bldArray(k,iCol) = abs(zw(k+1,iCol));
          interfaceForcings(k,iCol) = sfc_buoy(iCol);

          double sfc_layer_depth = bldArray(k,iCol)*surf_layer_ext;

          surfaceAverageIndex(k, iCol) = 0;
          for(int kIndexOBL=topIndex; kIndexOBL < k; kIndexOBL++){
            if(zw(kIndexOBL+1,iCol) < -sfc_layer_depth){
                surfaceAverageIndex(k,iCol) = kIndexOBL;
                break;
            }
          }
          topIndex = (int) max(cvmix_zeroD,surfaceAverageIndex(k,iCol) - 1);
        }

        // Compute turbulent scales ws here
        //
        if(sfc_buoy(iCol) >= 0) {
          for(int k=0; k<nVertLevels; k++){
            double L = pow(us(iCol),3) / (kappa * sfc_buoy(iCol));
            double sigma = -zt(k,iCol) / bldArray(k,iCol);
            double zeta = sigma * bldArray(k,iCol) / (L + 1.0E-10);

            ws(k,iCol) = kappa*us(iCol) / (1.+5.*zeta);
          }
        }
        else{
          for(int k=0; k < nVertLevels; k++){
            double L = pow(us(iCol),3) / (kappa * sfc_buoy(iCol) + 1.0E-10);
            double sigma = -zt(k,iCol) / bldArray(k,iCol);
            double sigma1 = min(sigma,surf_layer_ext);
            double zeta = sigma1 * bldArray(k,iCol) / (L + 1.0E-10);

            if(zeta < -1.0 || us(iCol) == cvmix_zeroD){
              ws(k,iCol) = kappa*pow((a_s*pow(us(iCol),3) - c_s*kappa*sigma1*bldArray(k,iCol)*sfc_buoy(iCol)),1./3.);
            }
            else{
              ws(k,iCol) = kappa*us(iCol)*pow((1.0 - 16.0*zeta),1.0/2.0);
            }
          }
        }
        // average quantities
        uVelocitySum(0,iCol) = u(0,iCol);
        vVelocitySum(0,iCol) = v(0,iCol);
        potentialDensitySum(0,iCol) = rho(0,iCol);
        for(int k=1; k < nVertLevels; k++){
          uVelocitySum(k,iCol) = uVelocitySum(k-1,iCol) + u(k,iCol);
          vVelocitySum(k,iCol) = vVelocitySum(k-1,iCol) + v(k,iCol);
          potentialDensitySum(k,iCol) = potentialDensitySum(k-1,iCol) + rho(k,iCol);
        }

        for(int k=0; k < nVertLevels; k++){
          int averageLocation = surfaceAverageIndex(k,iCol);
          double uAv = uVelocitySum(averageLocation,iCol) / (surfaceAverageIndex(k,iCol) + 1. );
          double vAv = vVelocitySum(averageLocation,iCol) / (surfaceAverageIndex(k,iCol) + 1.);
          double deltaVelocitySquared = pow((uAv - u(k,iCol)),2) + pow((vAv - v(k,iCol)),2);
          double rhoAv = potentialDensitySum(averageLocation,iCol) / (surfaceAverageIndex(k,iCol) + 1.);

          bulkRichardsonNumberShear(k,iCol) = max(deltaVelocitySquared, 1.0E-15);
          bulkRichardsonNumberBuoy(k,iCol) = gravity * (rho(k,iCol) - rhoAv) / rho_o;
        }

        // compute unresolved shear
        double Vtc = sqrt(beta / (c_s*surf_layer_ext)) / pow(kappa,2);
        for(int k = 0; k < nVertLevels; k++){
          N_cntr(k,iCol) = sqrt(max(n2(k+1,iCol),cvmix_zeroD));
          double CvT = 2.1 - 200.0*N_cntr(k,iCol);

         if( LANGMUIR_lifk17 ){
          double c_ST = 0.17;
          double c_CT = 0.15;
          double p_LT = 2.0;
          double c_LT = 0.083;
          double VTC = sqrt((c_CT*sfc_buoy(iCol)*zt(k,iCol) + c_ST*pow(us(iCol),3) +
                     c_LT*pow(us(iCol),3)*pow(LaSL(iCol),(-1.*p_LT)))/ws(k,iCol));
          Vt2(k,iCol) = max(-max(CvT,CVmin)*VTC*zt(k,iCol)*N_cntr(k,iCol)/Ri_crit,minVtsqr);
         } else{
          Vt2(k,iCol) = max(-max(CvT,CVmin)*Vtc*zt(k,iCol)*N_cntr(k,iCol)*ws(k,iCol) / Ri_crit, minVtsqr);
         }

    }

    // bulk richardson number computation // could probably fuse with above loop
    double scaling = cvmix_oneD - cvmix_p5D*surf_layer_ext;

    for(int k=0; k < nVertLevels; k++){
      double num = -scaling*zt(k,iCol)*bulkRichardsonNumberBuoy(k,iCol);
      double denom = bulkRichardsonNumberShear(k,iCol) + Vt2(k,iCol);
      bulkRi(k,iCol) = num / (denom + 1.0E-10);
    }
    // bld
    double OBL_limit = abs(zt(maxLevs(iCol)-1,iCol));

    if(lEkman){
      double Ekman = abs(zt(maxLevs(iCol)-1,iCol));
      if(f(iCol) != cvmix_zeroD) Ekman = 0.7*us(iCol) / abs(f(iCol));
      OBL_limit = min(OBL_limit, Ekman);
    }

    if(lMonOb){
      double moninObukhov = abs(zt(maxLevs(iCol)-1,iCol));
      if(sfc_buoy(iCol) > cvmix_zeroD) moninObukhov = pow(us(iCol),3) / (sfc_buoy(iCol)*kappa);
      OBL_limit = min(OBL_limit, moninObukhov);
    }

    int kindex = 0;
    for(int k=0; k < nVertLevels; k++){
      if(bulkRi(k,iCol) > Ri_crit){
        kindex= k;
        break;
      }
    }

    if(kindex == nVertLevels-1){
      h(iCol) = abs(OBL_limit);
    } else if(kindex == 0){
      h(iCol) = abs(zt(kindex,iCol));
    } else {
      // find OBL depth using quadratic interpolation
      double zt0 = (zt(kindex,iCol));
      double ztm1 = (zt(kindex-1,iCol));
      double det = -(pow(zt0 - ztm1,2));

      Minv(0,0,iCol) = -cvmix_oneD / det;
      Minv(0,1,iCol) = cvmix_oneD / det;
      Minv(0,2,iCol) = -cvmix_oneD / (zt0 - ztm1);
      Minv(1,0,iCol) = 2.0*ztm1 / det;
      Minv(1,1,iCol) = -2.0*ztm1 / det;
      Minv(1,2,iCol) = (zt0 + ztm1) / (zt0 - ztm1);
      Minv(2,0,iCol) = -(pow(ztm1,2)) / det;
      Minv(2,1,iCol) = zt0*(2.0*ztm1 - zt0) / det;
      Minv(2,2,iCol) = -zt0*ztm1 / (zt0 - ztm1);

      rhs(0,iCol) = bulkRi(kindex,iCol);
      rhs(1,iCol) = bulkRi(kindex-1,iCol);
      if(kindex==1){
        rhs(2,iCol) = cvmix_zeroD;
      } else{
        rhs(2,iCol) = (bulkRi(kindex-1,iCol) - bulkRi(kindex-2,iCol)) / (ztm1 - (zt(kindex-2,iCol)));
      }

      for(int k=0; k < 3; k++){
            coeffs(k,iCol) = 0.0;
        }
      for(int k2=0; k2 < 3; k2++){
        for(int k=0; k < 3; k++){
          coeffs(k2,iCol) = coeffs(k2,iCol) + Minv(2-k2,k,iCol)*rhs(k,iCol);
        }
      }

      // solve quadratic fit
      double a = coeffs(2,iCol);
      double b = coeffs(1,iCol);
      double c = coeffs(0,iCol) - Ri_crit;

      h(iCol) = (-b - sqrt(pow(b,2) - 4.0*a*c)) / (2.0*a);

      h(iCol) = min(-h(iCol), OBL_limit);
      h(iCol) = max(h(iCol),abs(zt(0,iCol))/2.0);
      if(iceFrac(iCol) > 0.05){
        h(iCol) = max(h(iCol), minOSBLUnderIce);
      }

    }
    //
    // compute intOBL
    //
    for(int kw = 0; kw < nVertLevels; kw++){
      if(h(iCol) < abs(zw(kw+1,iCol))){
        if(h(iCol) < abs(zt(kw,iCol))){
          intOBL(iCol) = (float)kw + 0.25;
        } else{
          intOBL(iCol) = (float)kw + 0.75;
        }
        break;
      }
    }

    // define a variable for stability and non-stability
    int kwup = floor(intOBL(iCol)) ;
    int ktup = round(intOBL(iCol)) - 1;
    // 10/23/18 -- Start with simple shapes only as it's easier

    double MshapeAtS;
    double TshapeAtS;
    double GAtS;

    double delta = (h(iCol)+zt(ktup,iCol))/(zt(ktup,iCol)-zt(ktup+1,iCol));
    if(ktup == maxLevs(iCol)){
        delta = (h(iCol)+zt(ktup,iCol))/(zt(ktup,iCol)-zw(ktup+1,iCol));
    }

    // compute wm and ws now and sigma
    if(sfc_buoy(iCol) >= 0) {
      for(int k=1; k<nVertLevels; k++){
        double L = pow(us(iCol),3) / (kappa * sfc_buoy(iCol));
        double sigma = -zw(k,iCol) / h(iCol);
        double zeta = sigma * h(iCol) / (L + 1.0E-10);

        ws(k,iCol) = kappa*us(iCol) / (1.+5.*zeta);
        wm(k,iCol) = ws(k,iCol);

        MshapeAtS = Mshape1 + Mshape2*sigma + Mshape3*pow(sigma,2.0) +
                    Mshape4*pow(sigma,3.0);

        TshapeAtS = Tshape1 + Tshape2*sigma + Tshape3*pow(sigma,2.0) +
                    Tshape4*pow(sigma,3.0);

        GAtS = 0.0;

       double MixingCoefEnhancement = cvmix_oneD;

       if( LANGMUIR_lifk17 ){
        MixingCoefEnhancement = Langmuir_EFactor(iCol);
       }

        OBL_Mdiff(k,iCol) = h(iCol)*wm(k,iCol)*MshapeAtS*maskw(k,iCol)*MixingCoefEnhancement;
        OBL_Tdiff(k,iCol) = h(iCol)*ws(k,iCol)*TshapeAtS*maskw(k,iCol)*MixingCoefEnhancement;
        nlt(k,iCol) = 0.0;
      }
    }
    else{
      for(int k=1; k < nVertLevels; k++){
        double L = pow(us(iCol),3) / (kappa * sfc_buoy(iCol) + 1.0E-10);
        double sigma = -zw(k,iCol) / h(iCol);
        double zeta = sigma * h(iCol) / (L + 1.0E-10);

        if(zeta < -0.2 || us(iCol) == cvmix_zeroD){
          double sigma1 = min(surf_layer_ext, sigma);
          wm(k,iCol) = kappa*pow((a_m*pow(us(iCol),3) - c_m*kappa*sigma1*h(iCol)*sfc_buoy(iCol)),1./3.);
        } else{
          wm(k,iCol) = kappa*us(iCol)*pow((1.0 - 16.0*zeta),0.25);
        }

        if(zeta < -1.0 || us(iCol) == cvmix_zeroD){
          double sigma1 = min(surf_layer_ext, sigma);
          ws(k,iCol) = kappa*pow((a_s*pow(us(iCol),3) - c_s*kappa*sigma1*h(iCol)*sfc_buoy(iCol)),1./3.);
        }
        else{
          ws(k,iCol) = kappa*us(iCol)*pow((1.0 - 16.0*zeta),1.0/2.0);
        }

        MshapeAtS = Mshape1 + Mshape2*sigma + Mshape3*pow(sigma,2.0) +
                    Mshape4*pow(sigma,3.0);

        TshapeAtS = Tshape1 + Tshape2*sigma + Tshape3*pow(sigma,2.0) +
                    Tshape4*pow(sigma,3.0);

        GAtS = TshapeNL1 + TshapeNL2*sigma + TshapeNL3*pow(sigma,2.0) +
                TshapeNL4*pow(sigma,3.0);

       double MixingCoefEnhancement = cvmix_oneD;

       if( LANGMUIR_LWF16 ){
        MixingCoefEnhancement = Langmuir_EFactor(iCol);
       }

        OBL_Mdiff(k,iCol) = h(iCol)*wm(k,iCol)*MshapeAtS*MixingCoefEnhancement;
        OBL_Tdiff(k,iCol) = h(iCol)*ws(k,iCol)*TshapeAtS*MixingCoefEnhancement;
        nlt(k,iCol) = non_local_coeff*GAtS;
      }
    }

    double sigma_ktup = -zt(ktup,iCol) / h(iCol);
    MshapeAtS = Mshape1 + Mshape2*sigma_ktup + Mshape3*pow(sigma_ktup,2.0) +
                Mshape4*pow(sigma_ktup,3.0);
    TshapeAtS = Tshape1 + Tshape2*sigma_ktup + Tshape3*pow(sigma_ktup,2.0) +
                Tshape4*pow(sigma_ktup,3.0);

    double ws_ktup, wm_ktup;

    if(sfc_buoy(iCol) >= 0){
        double L = pow(us(iCol),3) / (kappa * sfc_buoy(iCol));
        double zeta = sigma_ktup * h(iCol) / (L + 1.0E-10);

        ws_ktup = kappa*us(iCol) / (1. + 5.*zeta);
        wm_ktup = ws_ktup;
    } else{
        double L = pow(us(iCol),3) / (kappa * sfc_buoy(iCol));
        double sigma1 = min(sigma_ktup,surf_layer_ext);
        double zeta = sigma1 / (L + 1.0E-10);

        if(zeta < -0.2 || us(iCol) == 0){
            wm_ktup = kappa*pow((a_m*pow(us(iCol),3) - c_m*kappa*sigma1*h(iCol)*sfc_buoy(iCol)),1./3.);
        } else{
            wm_ktup = kappa*us(iCol)*pow((1.0-16.0*zeta),0.25);
        }

        if(zeta < -1.0 || us(iCol) == 0){
            ws_ktup = kappa*pow((a_s*pow(us(iCol),3.0) - c_s*kappa*sigma1*h(iCol)*sfc_buoy(iCol)),1./3.);
        } else{
            ws_ktup = kappa*us(iCol)*pow((1.0 - 16.0*zeta),1.0/2.0);
        }
    }

    double Mdiff_ktup = h(iCol) * wm_ktup * Langmuir_EFactor(iCol) * MshapeAtS;
    double Tdiff_ktup = h(iCol) * ws_ktup * Langmuir_EFactor(iCol) * TshapeAtS;

    double enh_Mdiff, enh_Tdiff;
    double omd = cvmix_oneD - delta;
    double old_Tdiff = OBL_Tdiff(ktup+1,iCol);

    // Enhanced diffusivity parameterization
    if(ktup == kwup){

        enh_Mdiff = pow(omd,2)*Mdiff_ktup + pow(delta,2)*kv(ktup+1,iCol);
        enh_Tdiff = pow(omd,2)*Tdiff_ktup + pow(delta,2)*kh(ktup+1,iCol);

        kv(ktup+1,iCol) = omd*kv(ktup+1,iCol) + delta*enh_Mdiff;
        kh(ktup+1,iCol) = omd*kh(ktup+1,iCol) + delta*enh_Tdiff;

        OBL_Mdiff(ktup+1,iCol) = kv(ktup+1,iCol);
        OBL_Tdiff(ktup+1,iCol) = kh(ktup+1,iCol);
    } else{
        enh_Mdiff = pow(omd,2)*Mdiff_ktup + pow(delta,2)*OBL_Mdiff(ktup+1,iCol);
        enh_Tdiff = pow(omd,2)*Tdiff_ktup + pow(delta,2)*OBL_Tdiff(ktup+1,iCol);

        OBL_Mdiff(ktup+1,iCol) = omd*kv(ktup+1,iCol) + delta*enh_Mdiff;
        OBL_Tdiff(ktup+1,iCol) = omd*kh(ktup+1,iCol) + delta*enh_Tdiff;

        if(old_Tdiff != cvmix_zeroD){
            nlt(ktup+1,iCol) = nlt(ktup+1,iCol)*OBL_Tdiff(ktup+1,iCol) / old_Tdiff;
        } else{
            nlt(ktup+1,iCol) = cvmix_zeroD;
        }
    }

    //Simply overwrite diffusivities in OBL -- ADD option to add diffusivity later
    for(int k=1; k< nVertLevels+1; k++){
        if(k <= kwup){
            kh(k,iCol) = OBL_Tdiff(k,iCol);
            kv(k,iCol) = OBL_Mdiff(k,iCol);
        } else{
            nlt(k,iCol) = cvmix_zeroD;
        }
    }

// TODO for now overwrite boundary layer diffusivity -- possibly add later

        // Convection -- TODO: Fix to add the other option, just step for now
        if(lBruntVaisala){
            for(int k=1; k< nVertLevels+1; k++){
                if(n2(k,iCol) <= cvmix_zeroD and k > kwup){
                    kh(k,iCol) = kh(k,iCol) + convect_diff;
                    kv(k,iCol) = kv(k,iCol) + convect_visc;
                }
            }
        }
    }

    if(lBruntVaisala and !lKPP){
        for(int k=1; k<= nVertLevels+1; k++){
            if(n2(k,iCol) <= cvmix_zeroD){
                kh(k,iCol) = kh(k,iCol) + convect_diff;
                kv(k,iCol) = kv(k,iCol) + convect_visc;
            }
        }
    }

    for(int k=1; k < nVertLevels+1; k++){
        kh(k,iCol) = kh(k,iCol)*maskw(k,iCol);
        kv(k,iCol) = kv(k,iCol)*maskw(k,iCol);
    }
    });

    Kokkos::deep_copy(h_kh,kh);
    Kokkos::deep_copy(h_kv,kv);
    Kokkos::deep_copy(h_OSBL,h);
    Kokkos::deep_copy(h_nlt,nlt);
    Kokkos::deep_copy(h_rib,bulkRi);

    for(int iCol=0; iCol < nCols; iCol++){
        // Pack data up for return to calling model
        for(int i=0; i<nVertLevels+1; i++){
            int ii = i+iCol*(nVertLevels+1);
            int k = ii / (nVertLevels+1);
            int j = ii - k*(nVertLevels+1);

            (diff)[ii] = h_kh(j,k);
            (visc)[ii] = h_kv(j,k);
            (NLT)[ii] = h_nlt(j,k);
        }

        (OSBL)[iCol] = h_OSBL(iCol);


    }

}
//need to read in langmuirefactor and return it
void compute_kpp_EFactor_model(double* u10, double* ustar, double* hbl,
                               double* laE, int nCols){
/*
! This function returns the enhancement factor, given the 10-meter
! wind (m/s), friction velocity (m/s) and the boundary layer depth (m).
!
! Based on CVMix implementation by Qing Li, 160606
*/

  double cvmix_zero = 0.0;
  double cvmix_one = 1.0;
// Local variables
  double us_sl, lasl_sqr_i;

  for(int i=0; i < nCols; i++){
    if (u10[i] > cvmix_zero && ustar[i] > cvmix_zero){
      // surface layer averaged Stokes drift
      us_sl = compute_kpp_ustokes_SL_model(u10[i], hbl[i]);

      // LaSL^{-2}
      lasl_sqr_i = us_sl/ustar[i];

      // enhancement factor (Li et al., 2016)
      laE[i] = sqrt(cvmix_one + cvmix_one/pow(1.5,2)*lasl_sqr_i
                 +cvmix_one/pow(5.4,4)*pow(lasl_sqr_i,2));
    } else{
      laE[i] = cvmix_one;
    }
  }
}

double compute_kpp_ustokes_SL_model(double u10, double hbl){
/*
! This function returns the surface layer averaged Stokes drift, given
! the 10-meter wind (m/s) and the boundary layer depth (m).
!
! Based on CVMix implementation by Qing Li, 20180130
*/
    // parameters
    // ratio of U19.5 to U10 (Holthuijsen, 2007)
    double u19p5_to_u10 = 1.075;
    // ratio of mean frequency to peak frequency for
    // Pierson-Moskowitz spectrum (Webb, 2011)
    double fm_to_fp = 1.296;
    // ratio of surface Stokes drift to U10
    double us_to_u10 = 0.0162;
    // loss ratio of Stokes transport
    double r_loss = 0.667;
		double cvmix_zero = 0.0;
		double cvmix_one = 1.0;
		double gravity = 9.8061;
		double cvmix_PI = 3.14159265358979323846;
    double kpp_ustokes_SL_model;

    if (u10 > cvmix_zero) {
      // surface Stokes drift
      double us = us_to_u10*u10;

      // significant wave height from Pierson-Moskowitz
      // spectrum (Bouws, 1998)
      double hm0 = 0.0246*pow(u10,2);

      // peak frequency (PM, Bouws, 1998)
      double tmp = 2.0*cvmix_PI*u19p5_to_u10*u10;
      double fp = 0.877*gravity/tmp;

      // mean frequency
      double fm = fm_to_fp*fp;

      /* total Stokes transport (a factor r_loss is applied to account
        for the effect of directional spreading, multidirectional waves
        and the use of PM peak frequency and PM significant wave height
        on estimating the Stokes transport) */
      double vstokes = 0.125*cvmix_PI*r_loss*fm*pow(hm0,2);

      /* the general peak wavenumber for Phillips' spectrum
       (Breivik et al., 2016) with correction of directional spreading */
      double kphil = 0.176*us/vstokes;

      /* surface layer averaged Stokes dirft with Stokes drift profile
       estimated from Phillips' spectrum (Breivik et al., 2016)
       the directional spreading effect from Webb and Fox-Kemper, 2015
       is also included */

      double kstar = kphil*2.56;
      // surface layer
      double z0 = 0.2*abs(hbl);
      double z0i = 1.0 / z0;
      // term 1 to 4
			double r1, r2, r3, r4;
      r1 = (0.151/kphil*z0i-0.84)*(cvmix_one-exp(-2.0*kphil*z0));
      r2 = -(0.84+0.0591/kphil*z0i)*sqrt(2.0*cvmix_PI*kphil*z0)
             *erfc(sqrt(2.0*kphil*z0));
      r3 = (0.0632/kstar*z0i+0.125)
            *(cvmix_one-exp(-2.0*kstar*z0));
      r4 = (0.125+0.0946/kstar*z0i)*sqrt(2.0*cvmix_PI*kstar*z0)
             *erfc(sqrt(2.0*kstar*z0));
      kpp_ustokes_SL_model = us*(0.715+r1+r2+r3+r4);
    } else {
      kpp_ustokes_SL_model = cvmix_zero;
    }

		return kpp_ustokes_SL_model;
} // extern"C"

