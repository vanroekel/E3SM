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
#include "Field.h"
#include "FieldGroup.h"
#include "KPPComputeOBLDepth.h"
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
   if (Err) {
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

   // OBL depth scheme
   std::string scheme_str;
   Err += KPPConfig.get("BoundaryLayerDepthScheme", scheme_str);
   if (scheme_str == "MatchBoth") {
      DefKPPMix->OBLDepthScheme = KPP::OBLDepthScheme::MatchBoth;
   } else {
      DefKPPMix->OBLDepthScheme = KPP::OBLDepthScheme::SimpleShapes;
   }

   // Wave and flux options
   Err += KPPConfig.get("UseLangmuirCirculation",
                        DefKPPMix->UseLangmuirCirculation);
   Err += KPPConfig.get("UseNonLocalFlux", DefKPPMix->UseNonLocalFlux);

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
   computeOBLDepth(SurfaceFrictionVelocity, SurfaceBuoyancyFlux,
                   BruntVaisalaFreqSq, IceFraction, WindSpeed10m);

   // =======================================================================
   // Stage 2: Compute Mixing Coefficients
   // =======================================================================
   computeMixingCoefficients(PotentialDensity, NormalVelocity,
                             TangentialVelocity);
}

/// Stage 1: Compute OBL depth using bulk Richardson search
void KPPMix::computeOBLDepth(const Array1DReal &SurfaceFrictionVelocity,
                             const Array1DReal &SurfaceBuoyancyFlux,
                             const Array2DReal &BruntVaisalaFreqSq,
                             const Array1DReal &IceFraction,
                             const Array1DReal &WindSpeed10m) {

   // =======================================================================
   // Pre-compute velocity shear squared (|dVel/dz|²)
   // For KPP, shear is computed from vertical differences in velocity
   // =======================================================================
   I4 NVertLayers = VCoord->NVertLayers;
   Array2DReal VelocityShearSquared("VelocityShearSquared", Mesh->NCellsAll,
                                    NVertLayers + 1);

   // Initialize velocity shear to small positive value (placeholder)
   // TODO Phase 2: Compute from actual velocity gradients
   parallelFor(
       "KPP-VelocityShear", {Mesh->NCellsAll, NVertLayers + 1},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          VelocityShearSquared(ICell, K) = 1.0e-8;
       });

   // =======================================================================
   // Compute Langmuir enhancement factors if wind speed is available
   // =======================================================================
   Array1DReal LangmuirFactor("LangmuirFactor", Mesh->NCellsAll);
   parallelFor(
       "KPP-Langmuir", {Mesh->NCellsAll}, KOKKOS_LAMBDA(I4 ICell) {
          if (!WindSpeed10m.is_empty() &&
              SurfaceFrictionVelocity(ICell) > KPP::MIN_USTAR) {
             Real h_bl_estimate = 50.0; // Placeholder
             Real u_stokes =
                 KPP::EstokesSLModel(WindSpeed10m(ICell), h_bl_estimate);
             Real langmuir_num = KPP::ComputeLangmuirNumber(
                 SurfaceFrictionVelocity(ICell), u_stokes);
             LangmuirFactor(ICell) = KPP::ComputeEnhancementFactor(
                 WindSpeed10m(ICell), SurfaceFrictionVelocity(ICell),
                 h_bl_estimate);
          } else {
             LangmuirFactor(ICell) = 1.0;
          }
       });

   // =======================================================================
   // Stage 1: Compute OBL depth using bulk Richardson criterion
   // =======================================================================

   // Capture mesh and coordinate data for use in lambda
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);
   OMEGA_SCOPE(ZInterface, VCoord->ZInterface);
   OMEGA_SCOPE(ZMid, VCoord->ZMid);
   OMEGA_SCOPE(LayerThickness, VCoord->LayerThickness);
   OMEGA_SCOPE(LocBoundaryLayerDepth, BoundaryLayerDepth);
   OMEGA_SCOPE(LocIndexBoundaryLayerDepth, IndexBoundaryLayerDepth);

   parallelFor(
       "KPP-OBLDepth", {Mesh->NCellsAll}, KOKKOS_LAMBDA(I4 ICell) {
          using namespace KPP;

          // Get cell-level data
          Real u_star   = SurfaceFrictionVelocity(ICell);
          Real b0       = SurfaceBuoyancyFlux(ICell);
          Real b0_enh   = b0 * LangmuirFactor(ICell);
          Real ice_frac = IceFraction(ICell);

          const I4 KMin = MinLayerCell(ICell);
          const I4 KMax = MaxLayerCell(ICell);

          // Check suppression condition
          if (ShouldSuppressOBL(ice_frac, 0)) {
             LocBoundaryLayerDepth(ICell)      = MIN_OBL_UNDER_ICE;
             LocIndexBoundaryLayerDepth(ICell) = KMin;
             return;
          }

          // Cumulative sum bulk Richardson search
          Real buoy_sum  = 0.0;
          Real shear_sum = 0.0;
          I4 k_obl       = KMin;

          for (I4 k = KMin; k <= KMax; ++k) {
             Real z_depth = Kokkos::abs(ZInterface(k + 1));

             if (z_depth < 1.0e-10) {
                continue; // Skip surface
             }

             // Accumulate Richardson components
             if (k > KMin) {
                Real dz = Kokkos::abs(ZMid(k) - ZMid(k - 1));
                buoy_sum += BruntVaisalaFreqSq(ICell, k) * dz;
                shear_sum += VelocityShearSquared(ICell, k) * dz;
             }

             // Compute bulk Richardson number
             Real u_star_sq =
                 Kokkos::fmax(MIN_USTAR * MIN_USTAR, u_star * u_star);
             Real denom =
                 u_star_sq * u_star_sq * Kokkos::fmax(1.0e-15, shear_sum);
             Real richardson = (b0_enh * z_depth * buoy_sum) / denom;

             // Check stopping criterion
             if (richardson >= CriticalRichardson) {
                k_obl = k;
                break;
             }
          }

          // Store OBL depth
          LocBoundaryLayerDepth(ICell) = ConstrainOBLDepth(
              Kokkos::abs(ZInterface(k_obl + 1)),
              Kokkos::abs(ZMid(KMin + 1) - ZMid(KMin)), 10000.0, ice_frac);
          LocIndexBoundaryLayerDepth(ICell) = k_obl;
       });

   LOG_INFO("KPPMix::computeOBLDepth: OBL depth computed");
}

/// Stage 2: Compute mixing coefficients within and below OBL
void KPPMix::computeMixingCoefficients(const Array2DReal &PotentialDensity,
                                       const Array2DReal &NormalVelocity,
                                       const Array2DReal &TangentialVelocity) {

   // =======================================================================
   // Initialize with background mixing
   // =======================================================================
   parallelFor(
       "KPP-Coeffs-Init", {Mesh->NCellsAll, VCoord->NVertLayers + 1},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          VertDiff(ICell, K)         = BackgroundDiff;
          VertVisc(ICell, K)         = BackgroundVisc;
          VertNonLocalFlux(ICell, K) = 0.0;
       });

   // =======================================================================
   // TODO Phase 2: Compute mixing coefficients from KPP profiles
   // - Evaluate KPP_M1, KPP_M2, KPP_S1, KPP_S2 profile functions
   // - Apply turbulent velocity scale and stability functions
   // - Compute diffusivity = K_m * w_s * g(σ)
   // - Compute viscosity = K_m * w_s * g_m(σ)
   // - Compute non-local flux if enabled
   // =======================================================================

   LOG_INFO("KPPMix::computeMixingCoefficients: Stub - Phase 2 implementation "
            "pending");
}

/// Register fields with I/O system
void KPPMix::defineFields() {
   // TODO: Register fields for output
   // - Create Field objects for VertDiff, VertVisc, etc.
   // - Create FieldGroup for organizing KPP diagnostics
   // - Attach Kokkos arrays to Field objects

   LOG_INFO("KPPMix::defineFields: Stub - I/O registration pending");
}

} // namespace OMEGA
