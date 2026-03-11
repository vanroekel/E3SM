//===-- ocn/KPPMix.cpp - KPP Boundary Layer Mixing Implementation --*- C++
//-*-===//
//
/// \file
/// \brief Implementation of KPP boundary layer mixing scheme
///
/// This file implements the KPPMix class for computing ocean boundary layer
/// mixing following Large et al. (1994) with optional Langmuir enhancement.
//
//===----------------------------------------------------------------------===//

#include "KPPMix.h"
#include "DataTypes.h"
#include "Error.h"
#include "GlobalConstants.h"
#include "KPPConstants.h"
#include "Logging.h"
#include "OmegaKokkos.h"
#include "VertCoord.h"
#include <limits>

namespace OMEGA {

// Singleton instance
KPPMix *KPPMix::Instance = nullptr;

/// Constructor for KPPMix
KPPMix::KPPMix(const std::string &Name_in, const HorzMesh *Mesh_in,
               const VertCoord *VCoord_in)
    : Name(Name_in), Mesh(Mesh_in), VCoord(VCoord_in) {

   // Allocate output arrays
   VertDiff = Array2DReal("VertDiff", Mesh->NCellsAll, VCoord->NVertLayers + 1);
   VertVisc = Array2DReal("VertVisc", Mesh->NCellsAll, VCoord->NVertLayers + 1);
   BoundaryLayerDepth = Array1DReal("BoundaryLayerDepth", Mesh->NCellsAll);
   IndexBoundaryLayerDepth =
       Array1DI4("IndexBoundaryLayerDepth", Mesh->NCellsAll);
   VertNonLocalFlux = Array2DReal("VertNonLocalFlux", Mesh->NCellsAll,
                                  VCoord->NVertLayers + 1);
      BulkRichardsonNumber = Array2DReal("BulkRichardsonNumber", Mesh->NCellsAll,
                                VCoord->NVertLayers + 1);
      TurbulentVelocityScale = Array2DReal("TurbulentVelocityScale",
                                 Mesh->NCellsAll,
                                 VCoord->NVertLayers + 1);

   // Set field names
   VertDiffFldName     = "VertDiff";
   VertViscFldName     = "VertVisc";
   OBLDepthFldName     = "BoundaryLayerDepth";
   NonLocalFluxFldName = "VertNonLocalFlux";
   BulkRichardsonFldName = "BulkRichardsonNumber";
   TurbulentVelScaleFldName = "TurbulentVelocityScale";

   if (Name != "Default") {
      VertDiffFldName.append(Name);
      VertViscFldName.append(Name);
      OBLDepthFldName.append(Name);
      NonLocalFluxFldName.append(Name);
      BulkRichardsonFldName.append(Name);
      TurbulentVelScaleFldName.append(Name);
   }

   defineFields();
}

/// Destructor for KPPMix
KPPMix::~KPPMix() {}

/// Get singleton instance
KPPMix *KPPMix::getInstance() { return Instance; }

/// Destroy singleton instance
void KPPMix::destroyInstance() {
   delete Instance;
   Instance = nullptr;
}

/// Initialize KPPMix from configuration
void KPPMix::init() {
   if (!Instance) {
      Instance = new KPPMix("Default", HorzMesh::getDefault(),
                            VertCoord::getDefault());
   }

   Error Err;
   KPPMix *DefKPPMix   = KPPMix::getInstance();
   Config *OmegaConfig = Config::getOmegaConfig();

   // Get VertMix config group
   Config VertMixConfig("VertMix");
   Err += OmegaConfig->get(VertMixConfig);
   CHECK_ERROR_ABORT(Err, "KPPMix::init: VertMix group not found in Config");

   // Get KPP config subgroup
   Config KPPConfig("KPP");
   Err += VertMixConfig.get(KPPConfig);
   if (Err.isFail()) {
      LOG_INFO("KPPMix::init: KPP subgroup not found, using defaults");
      return; // Continue with defaults
   }

   // Read KPP parameters
   bool enable = true;
   Err += KPPConfig.get("Enable", enable);
   DefKPPMix->Enabled = enable;

   Err += KPPConfig.get("CriticalBulkRichardsonNumber",
                        DefKPPMix->CriticalRichardson);
   Err += KPPConfig.get("StopOBLSearch", DefKPPMix->StopOBLSearchMult);
   Err += KPPConfig.get("SurfaceLayerExtent", DefKPPMix->SurfaceLayerExtent);

   // OBL depth scheme - store as string
   std::string scheme_str = "SimpleShapes";
   Err += KPPConfig.get("BoundaryLayerDepthScheme", scheme_str);
   DefKPPMix->OBLDepthSchemeStr = scheme_str;

   // Wave and flux options
   Err += KPPConfig.get("UseLangmuirCirculation",
                        DefKPPMix->UseLangmuirCirculation);
   Err += KPPConfig.get("UseNonLocalFlux", DefKPPMix->UseNonLocalFlux);
   Error DebugErr = KPPConfig.get("DebugDiagnostics", DefKPPMix->DebugDiagnostics);
   if (!DebugErr.isSuccess()) {
      DebugErr.reset();
      DefKPPMix->DebugDiagnostics = false;
   }

   // Background mixing
   Err += KPPConfig.get("BackgroundViscosity", DefKPPMix->BackgroundVisc);
   Err += KPPConfig.get("BackgroundDiffusivity", DefKPPMix->BackgroundDiff);

   LOG_INFO("KPPMix::init: KPP initialized");
}

/// Main computation routine
void KPPMix::computeKPPMix(const Array2DReal &PotentialDensity,
                           const Array2DReal &NormalVelocity,
                           const Array2DReal &TangentialVelocity,
                           const Array1DReal &SurfaceFrictionVelocity,
                           const Array1DReal &SurfaceBuoyancyFlux,
                           const Array2DReal &BruntVaisalaFreqSq,
                           const Array1DReal &IceFraction,
                           const Array1DReal &WindSpeed10m) {

   if (!Enabled) {
      return; // Skip if disabled
   }

   // =======================================================================
   // Stage 1: Compute OBL Depth
   // =======================================================================
   computeOBLDepth(PotentialDensity, NormalVelocity, TangentialVelocity,
                   SurfaceFrictionVelocity, SurfaceBuoyancyFlux,
                   BruntVaisalaFreqSq, IceFraction, WindSpeed10m);

   // =======================================================================
   // Stage 2: Compute Mixing Coefficients
   // =======================================================================
   computeMixingCoefficients(PotentialDensity, NormalVelocity,
                             TangentialVelocity, SurfaceFrictionVelocity,
                             SurfaceBuoyancyFlux);

   if (DebugDiagnostics) {
      logDiagnostics(PotentialDensity, NormalVelocity, TangentialVelocity,
                     SurfaceFrictionVelocity, SurfaceBuoyancyFlux,
                     WindSpeed10m);
   }
}

void KPPMix::logDiagnostics(const Array2DReal &PotentialDensity,
                            const Array2DReal &NormalVelocity,
                            const Array2DReal &TangentialVelocity,
                            const Array1DReal &SurfaceFrictionVelocity,
                            const Array1DReal &SurfaceBuoyancyFlux,
                            const Array1DReal &WindSpeed10m) {

   using namespace KPP;

   const auto MinLayerCellH = createHostMirrorCopy(VCoord->MinLayerCell);
   const auto MaxLayerCellH = createHostMirrorCopy(VCoord->MaxLayerCell);
   const auto ZInterfaceH   = createHostMirrorCopy(VCoord->ZInterface);
   const auto DensityH      = createHostMirrorCopy(PotentialDensity);
   const auto UVelH         = createHostMirrorCopy(NormalVelocity);
   const auto VVelH         = createHostMirrorCopy(TangentialVelocity);
   const auto UStarH        = createHostMirrorCopy(SurfaceFrictionVelocity);
   const auto B0H           = createHostMirrorCopy(SurfaceBuoyancyFlux);
   const auto OBLDepthH     = createHostMirrorCopy(BoundaryLayerDepth);
   const auto OBLIndexH     = createHostMirrorCopy(IndexBoundaryLayerDepth);

   const int NCellsAll = Mesh->NCellsAll;
   if (NCellsAll <= 0) {
      return;
   }

   const int ICell = 0;
   const int KMin  = MinLayerCellH(ICell);
   const int KMax  = MaxLayerCellH(ICell);
   const int NVertLayers = VCoord->NVertLayers;

   if (KMin > KMax) {
      return;
   }

   const Real rho_ref = 1025.0_Real;
   const int KSurf = Kokkos::min(KMin, NVertLayers - 1);
   const Real rho_surf = DensityH(ICell, KSurf);
   const Real u_surf   = UVelH(ICell, KSurf);
   const Real v_surf   = VVelH(ICell, KSurf);
   const Real u_star   = UStarH(ICell);
   const Real u_star_eff = Kokkos::fmax(KPP::MIN_USTAR, u_star);
   const Real b0       = B0H(ICell);
   Real u10 = 0.0_Real;
   if (WindSpeed10m.extent(0) > 0) {
      const auto Wind10mH = createHostMirrorCopy(WindSpeed10m);
      u10 = Wind10mH(ICell);
   }
      const Real langmuir_factor =
         UseLangmuirCirculation ? ComputeEnhancementFactor(u10, u_star_eff, 50.0)
                              : 1.0_Real;
   const Real b0_eff = b0 * langmuir_factor;

   LOG_INFO(
       "KPP debug: cell={} h_obl={} m k_obl={} u*={} b0={} b0_eff={} langmuir={}",
       ICell, OBLDepthH(ICell), OBLIndexH(ICell), u_star, b0, b0_eff,
       langmuir_factor);

   const int KTop = Kokkos::min(KMax, KMin + 3);
   const int k_obl = OBLIndexH(ICell);
   const Real h_obl = OBLDepthH(ICell);

   for (int K = KMin; K <= KTop; ++K) {
      const int kCell = Kokkos::min(K, NVertLayers - 1);
      const int kInt  = Kokkos::min(K + 1, NVertLayers);
      const Real z_depth = Kokkos::abs(ZInterfaceH(ICell, kInt));

      const Real rho_k = DensityH(ICell, kCell);
      const Real u_k   = UVelH(ICell, kCell);
      const Real v_k   = VVelH(ICell, kCell);

      const Real delta_rho = rho_k - rho_surf;
      const Real delta_b   = Gravity * delta_rho / rho_ref;
      const Real du        = u_k - u_surf;
      const Real dv        = v_k - v_surf;
      const Real shear2    = du * du + dv * dv;
        const Real w_turb =
           ComputeTurbulentVelocityScale(u_star_eff, b0_eff, z_depth);
      const Real vel_scale2 = shear2 + w_turb * w_turb + 1.0e-12_Real;
      const Real ri_b       = delta_b * z_depth / vel_scale2;

      Real sigma = 0.0_Real;
      if (K <= k_obl) {
         sigma = -1.0_Real * static_cast<Real>(K - KMin) /
                 static_cast<Real>(k_obl - KMin + 1);
         sigma = Kokkos::fmax(-1.0_Real, Kokkos::fmin(0.0_Real, sigma));
      }

      const Real z_local = -sigma * h_obl;
      Real zeta = 0.0_Real;
      const Real denom = VonKar * b0;
      if (Kokkos::abs(denom) > 1.0e-16_Real) {
         const Real l_mo = (u_star_eff * u_star_eff * u_star_eff) / denom;
         if (Kokkos::abs(l_mo) > 1.0e-16_Real) {
            zeta = z_local / l_mo;
         }
      }

      const Real phi_m = KPP::KPPProfileM2(zeta);
      const Real phi_s = KPP::KPPProfileS2(zeta);

      LOG_INFO(
          "KPP debug top: cell={} k={} z={} ri_b={} zeta={} phi_m={} phi_s={}",
          ICell, K, z_depth, ri_b, zeta, phi_m, phi_s);
   }
}

/// Stage 1: Compute OBL depth using bulk Richardson search  
void KPPMix::computeOBLDepth(const Array2DReal &PotentialDensity,
                             const Array2DReal &NormalVelocity,
                             const Array2DReal &TangentialVelocity,
                             const Array1DReal &SurfaceFrictionVelocity,
                             const Array1DReal &SurfaceBuoyancyFlux,
                             const Array2DReal &BruntVaisalaFreqSq,
                             const Array1DReal &IceFraction,
                             const Array1DReal &WindSpeed10m) {

   using namespace KPP;

   I4 NVertLayers = VCoord->NVertLayers;
   (void)BruntVaisalaFreqSq;

   // =======================================================================
   // Compute Langmuir enhancement factors if wind speed is available
   // =======================================================================
   Array1DReal LangmuirFactor("LangmuirFactor", Mesh->NCellsAll);
   parallelFor(
       "KPP-Langmuir", {Mesh->NCellsAll}, KOKKOS_LAMBDA(I4 ICell) {
          if (UseLangmuirCirculation) {
             const Real uStar = SurfaceFrictionVelocity(ICell);
             const Real u10   = (WindSpeed10m.extent(0) > 0)
                                    ? WindSpeed10m(ICell)
                                    : 0.0_Real;
             LangmuirFactor(ICell) = ComputeEnhancementFactor(u10, uStar, 50.0);
          } else {
             LangmuirFactor(ICell) = 1.0;
          }
       });

   // =======================================================================
   // Stage 1: Compute OBL depth - simplified implementation
   // =======================================================================

   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);
   OMEGA_SCOPE(ZInterface, VCoord->ZInterface);
   OMEGA_SCOPE(ZMid, VCoord->ZMid);
   OMEGA_SCOPE(LocPotentialDensity, PotentialDensity);
   OMEGA_SCOPE(LocNormalVelocity, NormalVelocity);
   OMEGA_SCOPE(LocTangentialVelocity, TangentialVelocity);
   OMEGA_SCOPE(LocIceFraction, IceFraction);
   OMEGA_SCOPE(LocLangmuirFactor, LangmuirFactor);
   OMEGA_SCOPE(LocBoundaryLayerDepth, BoundaryLayerDepth);
   OMEGA_SCOPE(LocIndexBoundaryLayerDepth, IndexBoundaryLayerDepth);
   OMEGA_SCOPE(LocBulkRichardson, BulkRichardsonNumber);

   deepCopy(BulkRichardsonNumber, 0.0_Real);

   parallelFor(
       "KPP-OBLDepth", {Mesh->NCellsAll}, KOKKOS_LAMBDA(I4 ICell) {
          using namespace KPP;

          Real u_star   = SurfaceFrictionVelocity(ICell);
          Real b0       = SurfaceBuoyancyFlux(ICell);

          const I4 KMin = MinLayerCell(ICell);
          const I4 KMax = MaxLayerCell(ICell);
          const I4 KIntTop  = Kokkos::min(KMin + 1, NVertLayers);
          const I4 KIntDeep = Kokkos::min(KMax + 1, NVertLayers);

          // Suppression under heavy ice
          const Real iceFrac = LocIceFraction(ICell);
          if (ShouldSuppressOBL(iceFrac, 0)) {
             LocBoundaryLayerDepth(ICell)      = MIN_OBL_UNDER_ICE;
             LocIndexBoundaryLayerDepth(ICell) = KMin;
             return;
          }

          // Full bulk-Richardson search from actual fields.
          const Real rho_ref  = 1025.0_Real;
          const I4 KSurf      = Kokkos::min(KMin, NVertLayers - 1);
          const Real rho_surf = LocPotentialDensity(ICell, KSurf);
          const Real u_surf   = LocNormalVelocity(ICell, KSurf);
          const Real v_surf   = LocTangentialVelocity(ICell, KSurf);

          I4 k_obl          = KMax;
          Real obl_depth    = Kokkos::abs(ZInterface(ICell, KIntDeep));
          bool found_obl    = false;
          const Real ri_crit = CriticalRichardson;
          const Real b0_eff  = b0 * LocLangmuirFactor(ICell);

          for (I4 k = KMin; k <= KMax; ++k) {
             const I4 kCell = Kokkos::min(k, NVertLayers - 1);
             const I4 kInt = Kokkos::min(k + 1, NVertLayers);
             const Real z_depth = Kokkos::abs(ZInterface(ICell, kInt));
             if (z_depth < 1.0e-12)
                continue;

             const Real rho_k = LocPotentialDensity(ICell, kCell);
             const Real u_k   = LocNormalVelocity(ICell, kCell);
             const Real v_k   = LocTangentialVelocity(ICell, kCell);

             const Real delta_rho = rho_k - rho_surf;
             const Real delta_b   = Gravity * delta_rho / rho_ref;

             const Real du = u_k - u_surf;
             const Real dv = v_k - v_surf;
             const Real shear2 = du * du + dv * dv;

             const Real w_turb =
                 ComputeTurbulentVelocityScale(u_star, b0_eff, z_depth);
             const Real vel_scale2 = shear2 + w_turb * w_turb + 1.0e-12;

             const Real ri_b = delta_b * z_depth / vel_scale2;
             LocBulkRichardson(ICell, kInt) = ri_b;

             if (ri_b >= ri_crit) {
                k_obl     = k;
                obl_depth = z_depth;
                found_obl = true;
                break;
             }
          }

          if (!found_obl) {
             k_obl     = KMax;
             obl_depth = Kokkos::abs(ZInterface(ICell, KIntDeep));
          }

          if (OBLDepthSchemeStr == "MatchBoth") {
             const Real surface_depth = Kokkos::abs(ZInterface(ICell, KIntTop));
             obl_depth                = 0.3 * surface_depth + 0.7 * obl_depth;
          }

          Real surface_thickness = 1.0;
          if (KMin + 1 <= KMax) {
             surface_thickness = Kokkos::abs(ZMid(ICell, KMin + 1) - ZMid(ICell, KMin));
          } else {
             surface_thickness =
                 Kokkos::abs(ZInterface(ICell, KIntTop) - ZInterface(ICell, KMin));
          }
          const Real water_depth = Kokkos::abs(ZInterface(ICell, KIntDeep));
          obl_depth = ConstrainOBLDepth(obl_depth, surface_thickness, water_depth,
                                        iceFrac);

          I4 k_final = KMax;
          for (I4 k = KMin; k < KMax; ++k) {
             const Real z_above = Kokkos::abs(ZInterface(ICell, k));
             const Real z_below = Kokkos::abs(ZInterface(ICell, k + 1));
             if (obl_depth >= z_above && obl_depth <= z_below) {
                k_final = k;
                break;
             }
          }

          LocBoundaryLayerDepth(ICell)      = obl_depth;
          LocIndexBoundaryLayerDepth(ICell) = k_final;
       });

   LOG_INFO("KPPMix::computeOBLDepth: OBL depth computed");
}

/// Stage 2: Compute mixing coefficients within and below OBL
void KPPMix::computeMixingCoefficients(const Array2DReal &PotentialDensity,
                                       const Array2DReal &NormalVelocity,
                                       const Array2DReal &TangentialVelocity,
                                       const Array1DReal &SurfaceFrictionVelocity,
                                       const Array1DReal &SurfaceBuoyancyFlux) {

   using namespace KPP;

   (void)PotentialDensity;
   (void)NormalVelocity;
   (void)TangentialVelocity;

   I4 NVertLayers = VCoord->NVertLayers;

   // =======================================================================
   // Capture data for use in lambda
   // =======================================================================
   OMEGA_SCOPE(LocBoundaryLayerDepth, BoundaryLayerDepth);
   OMEGA_SCOPE(LocIndexBoundaryLayerDepth, IndexBoundaryLayerDepth);
   OMEGA_SCOPE(LocVertDiff, VertDiff);
   OMEGA_SCOPE(LocVertVisc, VertVisc);
   OMEGA_SCOPE(LocVertNonLocalFlux, VertNonLocalFlux);
   OMEGA_SCOPE(LocTurbulentVelocityScale, TurbulentVelocityScale);
   OMEGA_SCOPE(LocSurfaceFrictionVelocity, SurfaceFrictionVelocity);
   OMEGA_SCOPE(LocSurfaceBuoyancyFlux, SurfaceBuoyancyFlux);
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);

      // Capture member variables for use in lambda
   Real LocBackgroundDiff = BackgroundDiff;
   Real LocBackgroundVisc = BackgroundVisc;
   bool LocUseNonLocalFlux = UseNonLocalFlux;
   bool LocParabolicNonLocal = (OBLDepthSchemeStr == "ParabolicNonLocal");
   const Real LocKappa = VonKar;

   // =======================================================================
   // Initialize with background mixing
   // =======================================================================
   parallelFor(
       "KPP-Coeffs-Init", {Mesh->NCellsAll, NVertLayers + 1},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          LocVertDiff(ICell, K)         = LocBackgroundDiff;
          LocVertVisc(ICell, K)         = LocBackgroundVisc;
          LocVertNonLocalFlux(ICell, K) = 0.0;
          LocTurbulentVelocityScale(ICell, K) = 0.0;
       });

   // =======================================================================
   // Stage 2: Compute KPP profile-based mixing coefficients
   // =======================================================================

   parallelFor(
       "KPP-MixingCoeffs", {Mesh->NCellsAll}, KOKKOS_LAMBDA(I4 ICell) {
          // Get OBL information for this cell
          Real h_obl = LocBoundaryLayerDepth(ICell);
          I4 k_obl   = LocIndexBoundaryLayerDepth(ICell);

          const I4 KMin = MinLayerCell(ICell);
          const I4 KMax = MaxLayerCell(ICell);

          // =============================================================
          // Compute turbulent velocity scale w_s
          // =============================================================
          Real u_star = LocSurfaceFrictionVelocity(ICell);
          Real b0     = LocSurfaceBuoyancyFlux(ICell);

          // =============================================================
          // Compute mixing coefficients at each interface
          // =============================================================
          for (I4 k = KMin; k <= KMax + 1; ++k) {
             // Check if within OBL
             if (k <= k_obl) {
                // Normalized depth: σ = -z/h_OBL
                Real sigma = -1.0 * static_cast<Real>(k - KMin) /
                             static_cast<Real>(k_obl - KMin + 1);
                sigma = Kokkos::fmax(-1.0, Kokkos::fmin(0.0, sigma));

               // CVMix-style turbulent scales: w = kappa*u*/phi in general,
               // with explicit free-convection limits when u*=0.
               const Real sigma_coord = -sigma; // [0,1]
               const Real sigma_loc =
                  Kokkos::fmin(SurfaceLayerExtent, Kokkos::fmax(0.0_Real, sigma_coord));

               Real zeta = 0.0_Real;
               Real w_m_turb = 0.0_Real;
               Real w_s_turb = 0.0_Real;

               if (u_star > 0.0_Real) {
                  const Real u3 = u_star * u_star * u_star;
                  zeta = sigma_loc * h_obl * b0 * LocKappa /
                       Kokkos::max(u3, 1.0e-20_Real);

                  const Real phi_m = KPP::KPPProfileM2(zeta);
                  const Real phi_s = KPP::KPPProfileS2(zeta);
                  const Real phi_inv_m =
                     1.0_Real / Kokkos::max(phi_m, 1.0e-12_Real);
                  const Real phi_inv_s =
                     1.0_Real / Kokkos::max(phi_s, 1.0e-12_Real);

                  w_m_turb = LocKappa * u_star * phi_inv_m;
                  w_s_turb = LocKappa * u_star * phi_inv_s;
               } else if (b0 < 0.0_Real) {
                  // Free-convection edge case (u*=0, unstable forcing).
                  const Real c_m = 16.0_Real;
                  const Real c_s = 16.0_Real;
                  const Real wm3 = -c_m * sigma_loc * h_obl * LocKappa * b0;
                  const Real ws3 = -c_s * sigma_loc * h_obl * LocKappa * b0;
                  w_m_turb = LocKappa * Kokkos::pow(Kokkos::max(0.0_Real, wm3),
                                           1.0_Real / 3.0_Real);
                  w_s_turb = LocKappa * Kokkos::pow(Kokkos::max(0.0_Real, ws3),
                                           1.0_Real / 3.0_Real);
               }

                // ========================================================
                // Momentum mixing: K_m = u_star * w_s * M1(σ) * M2(zeta)
                // ========================================================
                Real m1 = KPP::KPPProfileM1(sigma);
                LocVertVisc(ICell, k) = Kokkos::fmax(
                  KPP::MIN_COEFFICIENT, h_obl * w_m_turb * m1);

                // ========================================================
                // Tracer mixing: K_s = u_star * w_s * S1(σ) * S2(zeta)
                // ========================================================
                Real s1 = KPP::KPPProfileS1(sigma);
                LocVertDiff(ICell, k) = Kokkos::fmax(
                  KPP::MIN_COEFFICIENT, h_obl * w_s_turb * s1);
                        LocTurbulentVelocityScale(ICell, k) = w_s_turb;

                // ========================================================
                // Non-local flux: G(σ)
                // ========================================================
                if (LocUseNonLocalFlux) {
                   Real g_sigma = LocParabolicNonLocal
                                      ? KPP::KPPProfileGParabolicNonLocal(sigma)
                                      : KPP::KPPProfileG(sigma);
                   LocVertNonLocalFlux(ICell, k) = g_sigma;
                } else {
                   LocVertNonLocalFlux(ICell, k) = 0.0;
                }

             } else {
                // Below OBL: use background values
                LocVertDiff(ICell, k)         = LocBackgroundDiff;
                LocVertVisc(ICell, k)         = LocBackgroundVisc;
                LocVertNonLocalFlux(ICell, k) = 0.0;
               LocTurbulentVelocityScale(ICell, k) = 0.0;
             }
          }
       });

   LOG_INFO("KPPMix::computeMixingCoefficients: Phase 2 mixing coefficients "
            "computed");
}

/// Register fields with I/O system
void KPPMix::defineFields() {
   const Real FillValue = -9.99e30;

   // BoundaryLayerDepth on cells
   std::vector<std::string> CellDims(1);
   CellDims[0] = "NCells";
   auto OBLDepthField =
       Field::create(OBLDepthFldName,                // field name
                     "ocean boundary layer depth",  // long name
                     "m",                           // units
                     "",                            // CF standard name
                     0.0,                            // min valid value
                     std::numeric_limits<Real>::max(), // max valid value
                     FillValue,                      // fill value
                     1,                              // number of dims
                     CellDims);

   // KPP non-local tracer flux profile on cell-layer interfaces
   std::vector<std::string> FluxDims(2);
   FluxDims[0] = "NCells";
   FluxDims[1] = "NVertLayersP1";
   auto NonLocalFluxField =
       Field::create(NonLocalFluxFldName,            // field name
                     "KPP non-local tracer flux profile", // long name
                     "1",                           // units
                     "",                            // CF standard name
                     std::numeric_limits<Real>::lowest(), // min valid value
                     std::numeric_limits<Real>::max(),    // max valid value
                     FillValue,                      // fill value
                     2,                              // number of dims
                     FluxDims);

         auto BulkRichardsonField =
            Field::create(BulkRichardsonFldName,          // field name
                      "bulk Richardson number",      // long name
                      "1",                           // units
                      "",                            // CF standard name
                      std::numeric_limits<Real>::lowest(), // min valid value
                      std::numeric_limits<Real>::max(),    // max valid value
                      FillValue,                      // fill value
                      2,                              // number of dims
                      FluxDims);

         auto TurbulentVelScaleField =
            Field::create(TurbulentVelScaleFldName,       // field name
                      "KPP turbulent velocity scale", // long name
                      "m s-1",                       // units
                      "",                            // CF standard name
                      0.0,                            // min valid value
                      std::numeric_limits<Real>::max(), // max valid value
                      FillValue,                      // fill value
                      2,                              // number of dims
                      FluxDims);

   // Group KPP-specific outputs for convenient stream selection.
   auto KPPGroup = FieldGroup::create("KPPMix");
   KPPGroup->addField(OBLDepthFldName);
   KPPGroup->addField(NonLocalFluxFldName);
      KPPGroup->addField(BulkRichardsonFldName);
      KPPGroup->addField(TurbulentVelScaleFldName);

   OBLDepthField->attachData<Array1DReal>(BoundaryLayerDepth);
   NonLocalFluxField->attachData<Array2DReal>(VertNonLocalFlux);
   BulkRichardsonField->attachData<Array2DReal>(BulkRichardsonNumber);
   TurbulentVelScaleField->attachData<Array2DReal>(TurbulentVelocityScale);

   LOG_INFO("KPPMix::defineFields: registered {}, {}, {}, {}",
            OBLDepthFldName, NonLocalFluxFldName,
            BulkRichardsonFldName, TurbulentVelScaleFldName);
}

} // namespace OMEGA
