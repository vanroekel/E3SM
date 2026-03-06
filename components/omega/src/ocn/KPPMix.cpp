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
#include "VertMix.h"
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

   // Set field names
   VertDiffFldName     = "VertDiff";
   VertViscFldName     = "VertVisc";
   OBLDepthFldName     = "BoundaryLayerDepth";
   NonLocalFluxFldName = "VertNonLocalFlux";

   if (Name != "Default") {
      VertDiffFldName.append(Name);
      VertViscFldName.append(Name);
      OBLDepthFldName.append(Name);
      NonLocalFluxFldName.append(Name);
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

   // Start from VertMix background defaults when available
   if (auto *BaseVertMix = VertMix::getInstance(); BaseVertMix != nullptr) {
      DefKPPMix->BackgroundVisc = BaseVertMix->BackVisc;
      DefKPPMix->BackgroundDiff = BaseVertMix->BackDiff;
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

   // Optional KPP-specific background overrides
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
   const Real LocRhoRef = VCoord->Rho0;

   parallelFor(
       "KPP-OBLDepth", {Mesh->NCellsAll}, KOKKOS_LAMBDA(I4 ICell) {
          using namespace KPP;

          Real u_star   = SurfaceFrictionVelocity(ICell);
          Real b0       = SurfaceBuoyancyFlux(ICell);

          const I4 KMin = MinLayerCell(ICell);
          const I4 KMax = MaxLayerCell(ICell);
          const I4 KIntTop  = Kokkos::min(KMin + 1, NVertLayers);
          const I4 KIntDeep = Kokkos::min(KMax + 1, NVertLayers);

          // Ice treatment: disable Langmuir under ice (>= 5%), but always compute OBL
          const Real iceFrac = LocIceFraction(ICell);
          Real lang_factor = LocLangmuirFactor(ICell);
          if (iceFrac >= 0.05_Real) {
             lang_factor = 1.0_Real; // No wave effects under ice
          }

          // Full bulk-Richardson search from actual fields
          const I4 KSurf      = Kokkos::min(KMin, NVertLayers - 1);

          I4 k_obl          = KMax;
          Real obl_depth    = Kokkos::abs(ZInterface(ICell, KIntDeep));
          bool found_obl    = false;
          const Real ri_crit = CriticalRichardson;
          const Real b0_eff  = b0 * lang_factor;

          for (I4 k = KMin; k <= KMax; ++k) {
             const I4 kCell = Kokkos::min(k, NVertLayers - 1);
             const I4 kInt = Kokkos::min(k + 1, NVertLayers);
             const Real z_depth = Kokkos::abs(ZInterface(ICell, kInt));
             if (z_depth < 1.0e-12)
                continue;

             const Real z_top = Kokkos::abs(ZInterface(ICell, KMin));
             const Real top_layer_thickness =
                 Kokkos::fmax(1.0e-12_Real,
                              Kokkos::abs(ZInterface(ICell, KIntTop) -
                                          ZInterface(ICell, KMin)));
             const Real avg_depth =
                 Kokkos::fmin(Kokkos::fmax(SurfaceLayerExtent * z_depth,
                                           top_layer_thickness),
                              Kokkos::abs(ZInterface(ICell, KIntDeep)));

             Real rho_surf_num = 0.0_Real;
             Real u_surf_num   = 0.0_Real;
             Real v_surf_num   = 0.0_Real;
             Real surf_wgt_sum = 0.0_Real;

             for (I4 ks = KMin; ks <= KMax; ++ks) {
                const I4 ksCell = Kokkos::min(ks, NVertLayers - 1);
                const I4 ksInt  = Kokkos::min(ks + 1, NVertLayers);
                const Real z_layer_top = Kokkos::abs(ZInterface(ICell, ks));
                const Real z_layer_bot = Kokkos::abs(ZInterface(ICell, ksInt));

                if (z_layer_top >= avg_depth) {
                   break;
                }

                const Real z_overlap_top = Kokkos::fmax(z_layer_top, z_top);
                const Real z_overlap_bot = Kokkos::fmin(z_layer_bot, avg_depth);
                const Real dz_overlap = z_overlap_bot - z_overlap_top;
                if (dz_overlap <= 0.0_Real) {
                   continue;
                }

                const Real rho_ks = LocPotentialDensity(ICell, ksCell);
                const Real u_ks   = LocNormalVelocity(ICell, ksCell);
                const Real v_ks   = LocTangentialVelocity(ICell, ksCell);

                rho_surf_num += rho_ks * dz_overlap;
                u_surf_num += u_ks * dz_overlap;
                v_surf_num += v_ks * dz_overlap;
                surf_wgt_sum += dz_overlap;
             }

             const Real rho_surf =
                 (surf_wgt_sum > 0.0_Real)
                     ? rho_surf_num / surf_wgt_sum
                     : LocPotentialDensity(ICell, KSurf);
             const Real u_surf =
                 (surf_wgt_sum > 0.0_Real)
                     ? u_surf_num / surf_wgt_sum
                     : LocNormalVelocity(ICell, KSurf);
             const Real v_surf =
                 (surf_wgt_sum > 0.0_Real)
                     ? v_surf_num / surf_wgt_sum
                     : LocTangentialVelocity(ICell, KSurf);

             const Real rho_k = LocPotentialDensity(ICell, kCell);
             const Real u_k   = LocNormalVelocity(ICell, kCell);
             const Real v_k   = LocTangentialVelocity(ICell, kCell);

             const Real delta_rho = rho_k - rho_surf;
             const Real delta_b   = Gravity * delta_rho / LocRhoRef;

             const Real du = u_k - u_surf;
             const Real dv = v_k - v_surf;
             const Real shear2 = du * du + dv * dv;

             const Real w_turb =
                 ComputeTurbulentVelocityScale(u_star, b0_eff, z_depth);
             const Real vel_scale2 = shear2 + w_turb * w_turb + 1.0e-12;

             const Real ri_b = delta_b * z_depth / vel_scale2;

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

          // MatchBoth does not directly modify OBL depth in this implementation.

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
   // Pre-compute gradient Richardson numbers for stability correction
   // =======================================================================
   Array2DReal GradientRichardsonNum("GradientRichardsonNum", Mesh->NCellsAll,
                                     NVertLayers + 1);

   // Simplified gradient Richardson computation
   parallelFor(
       "KPP-GradRichardson", {Mesh->NCellsAll, NVertLayers},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          // Placeholder: set to typical stable value
          GradientRichardsonNum(ICell, K) = 0.1;
       });

   // =======================================================================
   // Capture data for use in lambda
   // =======================================================================
   OMEGA_SCOPE(LocBoundaryLayerDepth, BoundaryLayerDepth);
   OMEGA_SCOPE(LocIndexBoundaryLayerDepth, IndexBoundaryLayerDepth);
   OMEGA_SCOPE(LocVertDiff, VertDiff);
   OMEGA_SCOPE(LocVertVisc, VertVisc);
   OMEGA_SCOPE(LocVertNonLocalFlux, VertNonLocalFlux);
   OMEGA_SCOPE(LocGradRichardson, GradientRichardsonNum);
   OMEGA_SCOPE(LocSurfaceFrictionVelocity, SurfaceFrictionVelocity);
   OMEGA_SCOPE(LocSurfaceBuoyancyFlux, SurfaceBuoyancyFlux);
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);

   // Capture member variables for use in lambda
   Real LocBackgroundDiff = BackgroundDiff;
   Real LocBackgroundVisc = BackgroundVisc;
   if (auto *BaseVertMix = VertMix::getInstance(); BaseVertMix != nullptr) {
      LocBackgroundDiff = BaseVertMix->BackDiff;
      LocBackgroundVisc = BaseVertMix->BackVisc;
   }
   bool LocUseNonLocalFlux = UseNonLocalFlux;
   bool LocParabolicNonLocal = (OBLDepthSchemeStr == "ParabolicNonLocal");

   // =======================================================================
   // Initialize with background mixing using shared VertMix-style deepCopy
   // =======================================================================
   deepCopy(LocVertDiff, LocBackgroundDiff);
   deepCopy(LocVertVisc, LocBackgroundVisc);
   deepCopy(LocVertNonLocalFlux, 0.0_Real);

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
          Real w_s    = KPP::ComputeTurbulentVelocityScale(u_star, b0, h_obl);

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

                // Get gradient Richardson at this layer
                Real ri_g = (k < NVertLayers) ? LocGradRichardson(ICell, k)
                                              : 0.0;

                // ========================================================
                // Momentum mixing: K_m = u_star * w_s * M1(σ) * M2(σ, Ri_g)
                // ========================================================
                Real m1 = KPP::KPPProfileM1(sigma);
                Real m2 = KPP::KPPProfileM2(sigma, ri_g);

                LocVertVisc(ICell, k) = Kokkos::fmax(
                    KPP::MIN_COEFFICIENT, u_star * w_s * m1 * m2);

                // ========================================================
                // Tracer mixing: K_s = u_star * w_s * S1(σ) * S2(σ, Ri_g)
                // ========================================================
                Real s1 = KPP::KPPProfileS1(sigma);
                Real s2 = KPP::KPPProfileS2(sigma, ri_g);

                LocVertDiff(ICell, k) = Kokkos::fmax(
                    KPP::MIN_COEFFICIENT, u_star * w_s * s1 * s2);

                // ========================================================
                // Non-local flux: G(σ)
                // ========================================================
                if (LocUseNonLocalFlux) {
                   Real g_sigma = LocParabolicNonLocal
                                      ? KPP::KPPProfileGParabolicNonLocal(sigma)
                                      : KPP::KPPProfileG(sigma);
                   Real s2_flux  = KPP::KPPProfileS2(sigma, ri_g);
                   LocVertNonLocalFlux(ICell, k) = g_sigma * s2_flux;
                } else {
                   LocVertNonLocalFlux(ICell, k) = 0.0;
                }

             } else {
                // Below OBL: use background values
                LocVertDiff(ICell, k)         = LocBackgroundDiff;
                LocVertVisc(ICell, k)         = LocBackgroundVisc;
                LocVertNonLocalFlux(ICell, k) = 0.0;
             }
          }
       });

   LOG_INFO("KPPMix::computeMixingCoefficients: Phase 2 mixing coefficients "
            "computed");
}

/// Register fields with I/O system
void KPPMix::defineFields() {
   LOG_INFO("KPPMix::defineFields: Stub - I/O registration pending");
}

} // namespace OMEGA
