# K-Profile Parameterization (KPP) Port to OMEGA

## Executive Summary

This document specifies the porting of the CVMix K-Profile Parameterization (KPP) scheme to the OMEGA ocean model following OMEGA's C++/Kokkos architecture and design patterns. The implementation will **augment** the current VertMix module with a sophisticated KPP boundary layer parameterization, keeping the existing ConvectiveMix and ShearMix schemes as optional alternatives.

---

## 1. Overview of KPP

The K-Profile Parameterization is a non-local mixing scheme that represents subgrid-scale mixing in the ocean boundary layer. Key features include:

### 1.1 Core Components
- **Boundary Layer Depth (OBL)**: Computed from bulk Richardson number criterion
- **Mixed-layer Mixing**: Enhanced (non-local) diffusion within the boundary layer
- **Non-local Transport**: Vertical transport of tracers due to coherent plumes
- **Langmuir Circulation**: Wave-driven enhancement factors (theory-based)
- **Stratification Effects**: Accounts for Brunt-Väisälä frequency and shear

### 1.2 Key CVMix Algorithms to Port
From `cvmix_kpp.F90`:
- `cvmix_kpp_compute_bulk_Richardson()` - Compute bulk Richardson number
- `cvmix_kpp_compute_OBL_depth()` - Find OBL depth from Richardson criterion
- `cvmix_kpp_compute_turbulent_scales()` - Calculate turbulent velocity scales
- `cvmix_coeffs_kpp()` - Main mixing coefficient calculation
- `cvmix_kpp_EFactor_model()` - Langmuir enhancement from wind

---

## 2. Architecture and Design

### 2.1 Design Philosophy
- Follow OMEGA's **functor-based Kokkos** pattern (like VertMix)
- Use **singleton pattern** for instance management
- Implement as **Kokkos functors** for GPU compatibility
- Maintain **field-based I/O** integration
- Support **configuration-driven** parameters
- **Augment VertMix**: Coexist with existing schemes
- **Place in src/ocn/**: Alongside VertMix.h/cpp

### 2.2 Class Structure

```cpp
class KPPMix {
 public:
   // Singleton pattern
   static KPPMix* getInstance();
   static void init();
   static void destroyInstance();

   // Output fields
   Array2DReal VertDiff;              // Vertical diffusivity (m²/s)
   Array2DReal VertVisc;              // Vertical viscosity (m²/s)
   Array1DReal BoundaryLayerDepth;    // OBL depth (m)
   Array1DI4 IndexBoundaryLayerDepth; // OBL as layer index
   Array2DReal VertNonLocalFlux;      // Non-local flux coefficient

   // Main computation routine
   void computeKPPMix(const Array2DReal& PotentialDensity,
                      const Array2DReal& NormalVelocity,
                      const Array2DReal& TangentialVelocity,
                      const Array1DReal& SurfaceFrictionVelocity,
                      const Array1DReal& SurfaceBuoyancyFlux,
                      const Array2DReal& BruntVaisalaFreqSq,
                      const Array1DReal& IceFraction);

 private:
   KPPMix();
   ~KPPMix();
};
```

### 2.3 File Structure

New files in `src/ocn/`:
```
├── KPPMix.h                    (Main class)
├── KPPMix.cpp                  (Implementation)
├── KPPConstants.h              (Parameters and profiles)
├── KPPComputeOBLDepth.h        (OBL depth functor)
├── KPPComputeBulkRichardson.h  (Richardson functor)
├── KPPNonLocalFlux.h           (Non-local flux functor)
└── (VertMix.h/cpp modified)
```

---

## 3. Input Requirements

### 3.1 Required State Variables
- `PotentialDensity[nCells][nLevels]` - From Eos calculation (kg/m³)
- `NormalVelocity[nEdges][nLevels]` - (m/s)
- `TangentialVelocity[nEdges][nLevels]` - (m/s)

### 3.2 Required Forcing/Auxiliary
- `SurfaceFrictionVelocity[nCells]` - u* from wind stress (m/s)
- `SurfaceBuoyancyFlux[nCells]` - Surface buoyancy flux (m²/s³)
- `BruntVaisalaFreqSq[nCells][nLevels+1]` - N² at interfaces (s⁻²)
- `IceFraction[nCells]` - Sea ice coverage (0-1)
  - Used to disable Langmuir enhancement when ≥ 0.05
  - Used to apply minimum OBL depth constraint when ≥ 0.15
- `LandIceMask[nCells]` - Land ice mask (0 or 1, optional)
  - Disables Langmuir enhancement under land ice
- `WindSpeed10m[nCells]` - For Langmuir enhancement (m/s)

### 3.3 Grid/Mesh
- `zw_iface`, `zt_cntr` - Depth arrays (m)
- `LayerThickness[nCells][nLevels]` - (m)
- `Coriolis_F[nCells]` - (1/s)
- Mesh topology via HorzMesh

---

## 4. Output Variables

### 4.1 Main Outputs (required)
1. **VertDiff[nCells][nLevels+1]** - Vertical diffusivity (m²/s)
2. **VertVisc[nCells][nLevels+1]** - Vertical viscosity (m²/s)
3. **BoundaryLayerDepth[nCells]** - OBL depth (m)
4. **IndexBoundaryLayerDepth[nCells]** - OBL layer index

### 4.2 Non-local Flux Terms
- **VertNonLocalFlux[nCells][nLevels+1]** - Profile G(σ) (dimensionless)
  - Applied: `d(T)/dt = -d/dz(G(σ) × Q_surf)`
  - Reference: [mpas_ocn_tracer_nonlocalflux.F](mpas-ocean/src/shared/mpas_ocn_tracer_nonlocalflux.F)
  - Non-zero only within OBL

---

## 5. Algorithm Summary

### 5.1 Two-Stage Computation

**Stage 1: Boundary Layer Depth** (~O(nCells × nLevels))
For each cell:
1. Check ice conditions:
   - If `landIceMask=0 AND iceFraction < 0.05`: Compute Langmuir enhancement
   - Else: Set Langmuir enhancement = 1.0 (no wave effects under ice)

2. Initialize cumulative sums for thickness-weighted surface layer averaging:
   - `densitySum[1] = density[1] × thickness[1]`
   - `thicknessSum[1] = thickness[1]`
   - `velocitySum_edge[1] = velocity_edge[1] × thickness_edge[1]` (for each edge)
   - `thicknessEdgeSum[1] = thickness_edge[1]` (for each edge)

3. **Loop through layers k=1 to nLevels computing bulk Richardson**:
   
   a. **Set test OBL depth** = depth at bottom of layer k
   
   b. **Determine surface layer extent**: `d_surface = surfLayerExtent × OBL_test` (typically 0.1)
   
   c. **Find surface layer index**: Deepest layer within d_surface
   
   d. **Update cumulative sums** for next iteration:
      - `thicknessSum[k] = thicknessSum[k-1] + thickness[k]`
      - `densitySum[k] = densitySum[k-1] + thickness[k] × density[k]`
      - (similarly for velocity on each edge)
   
   e. **Compute surface layer averages** (thickness-weighted):
      - `avgDensity = densitySum[surfaceLayerIndex[k]] / thicknessSum[surfaceLayerIndex[k]]`
      - `avgVelocity = velocitySum[surfaceLayerIndex[k]] / thicknessEdgeSum[surfaceLayerIndex[k]]`
   
   f. **Compute differences from surface layer average**:
      - `ΔB[k] = g × (density[k] - avgDensity) / ρ_ref`
      - `ΔV²[k] = |avgVelocity - velocity[k]|²` (averaged over cell edges)
   
   g. **Compute turbulent velocity scale** w_s for this potential OBL depth:
      - `w_s[k] = turbulent_scale(σ=surfLayerExtent, OBL_test, u*, B_surface)`
   
   h. **Compute unresolved shear**: `Vt²[k] = f(w_s[k], N², Langmuir_factor)`
   
   i. **Compute bulk Richardson number**:
      - `scaling = 1.0 - 0.5 × surfLayerExtent`
      - `Ri_bulk[k] = -scaling × z_center[k] × ΔB[k] / (ΔV²[k] + Vt²[k])`
   
   j. **Check stopping criterion**: If `Ri_bulk[k] > Ri_crit × stopFactor`: break

4. **Interpolate** to find exact OBL depth where Ri_b = Ri_crit (MatchBoth or SimpleShapes)

5. **Apply OBL constraints**:
   - Minimum: `OBL = max(OBL, 0.5 × thickness[1])`
   - If `iceFraction > 0.15`: `OBL = max(OBL, minimumOBLUnderSeaIce)`
   - Maximum: `OBL = min(OBL, abs(z_center[nLevels]))`

6. Store OBL depth and layer index

**Stage 2: Mixing Coefficients** (~O(nCells × nLevels))
For each cell:
1. Compute turbulent velocity scales at all depths within OBL: w_m(σ), w_s(σ)
2. Apply KPP profile within OBL: ν(z) = u* × w_s(σ) × w(σ) + ν_bg
3. Compute non-local flux G(σ) within OBL
4. Store background mixing below OBL

### 5.2 Non-local Flux Details

**Vertical Divergence Form:**
`d(T)/dt = -d/dz(G(σ) × Q_surf)`

**Physical Meaning:**
- Non-local flux represents coherent plume transport within the boundary layer
- Only active within OBL (G=0 below OBL depth)
- Magnitude depends on surface forcing and boundary layer profile shape
- Zero under ice if using non-local flux only for buoyancy forcing regime

**Omega Implementation:**
1. KPP computes `VertNonLocalFlux[nCells][nLevels+1]` with G(σ) within OBL
2. In Tendencies: Loop over tracers, apply: `tend += (G(k) - G(k+1)) × Q_surf`
3. Zero boundary conditions at surface and below OBL
4. G(σ) typically follows parabolic or cubic profile within OBL

See [mpas_ocn_tracer_nonlocalflux.F](mpas-ocean/src/shared/mpas_ocn_tracer_nonlocalflux.F) for reference implementation.

### 5.3 Key Quantities

| Variable | Units | Range | Notes |
|----------|-------|-------|-------|
| u\* | m/s | 0.0001+ | Friction velocity (no upper limit) |
| B₀ | m²/s³ | -0.5 to 0.01 | Buoyancy flux |
| h_OBL | m | 5–500 | Boundary layer depth |
| Ri_b | — | 0–1 | Bulk Richardson = (d-d_r)×ΔB/(ΔV²+Vt²) where ΔB,ΔV are differences from **surface layer average** |
| κ | m²/s | 1e-6 to 1e-2 | Coefficient |

---

## 6. Configuration Parameters

### 6.1 YAML Config Structure

```yaml
VertMix:
  Scheme: KPP              # KPP, SimpleConvectiveShear, or Background

  KPP:
    Enable: true

    # OBL depth method
    BoundaryLayerDepthScheme: MatchBoth    # MatchBoth or SimpleShapes
    # MatchBoth: Interpolate OBL between Ri_crit and surface points
    # SimpleShapes: Traditional - OBL where Ri_crit criterion met

    # OBL parameters
    CriticalBulkRichardsonNumber: 0.3
    StopOBLSearch: 1.0
    SurfaceLayerExtent: 0.1

    # OBL bounds
    MinimumOBLDepth: 0.0                   # (m) typically 0.5 × first layer thickness
    MaximumOBLDepth: 0.0                   # (m) 0 = use bottom depth
    MinimumOBLUnderSeaIce: 10.0            # (m) applied when iceFraction > 0.15

    # Ice thresholds
    IceFractionThresholdForLangmuir: 0.05  # Disable Langmuir when iceFraction ≥ this
    IceFractionThresholdForMinimumOBL: 0.15 # Apply minimum OBL when iceFraction ≥ this

    # Wave enhancement (disabled automatically under ice)
    UseLangmuirCirculation: true           # Theory-based model

    # Background mixing
    BackgroundViscosity: 1.0e-4
    BackgroundDiffusivity: 1.0e-5

    # Non-local flux
    UseNonLocalFlux: true
```

### 6.2 Hard-Coded Constants (KPPConstants.h)

```cpp
namespace OMEGA::KPP {
  const Real PEC_LANGMUIR = 0.5;      // Langmuir parameterization
  const Real ZETA_M_SCALE = 0.4;      // Momentum scale
  const Real ZETA_S_SCALE = 0.16;     // Tracer scale
  const Real HUON = 0.03;             // Surface mixing coefficient
  const Real BD = 1.0;                // Buoyancy parameter
  const Real MIN_USTAR = 1.0e-4;      // Min friction velocity (m/s)
}
```

---

## 7. Integration Points

### 7.1 Code Modifications

**src/ocn/VertMix.h/cpp**
- Add `bool EnableKPP` flag, read from config
- Modify dispatch in `computeVertMix()` to call KPPMix if enabled
- Keep ConvectiveMix/ShearMix intact for backward compatibility

**src/ocn/Tendencies.cpp**
- Add non-local flux divergence: `tend += d/dz(G × Q_surf)`
- Reference mpas_ocn_tracer_nonlocalflux.F for pattern
- Apply within OBL region only

**Dependencies**
- Requires: `Config`, `Field`, `FieldGroup`, `HorzMesh`, `VertCoord`, `Eos`
- Requires: `BruntVaisalaFreq` computed first

---

## 8. Parallel Implementation

### 8.1 Kokkos Pattern

```cpp
parallelForOuter("KPPMix", {Mesh->NCellsAll},
  KOKKOS_LAMBDA(I4 ICell, const TeamMember &Team) {
    // Stage 1: OBL depth (cell-wise)
    KPPComputeOBLDepth obl_functor(...);
    Real h_obl = obl_functor(ICell);

    // Stage 2: Coefficients (cell-K chunked)
    KPPComputeCoeffs coeff_functor(...);
    parallelForInner(Team, vertRangeChunked(...),
      KOKKOS_LAMBDA(int KChunk) {
        coeff_functor(ICell, KChunk, h_obl);
      });
  });
```

### 8.2 Performance Notes on Stage 1

⚠️ **Performance Optimization**: Bulk Richardson loop uses cumulative sum pattern
- **Cumulative sums**: Maintain running totals for density and velocity to avoid O(nLevels²)
- **Surface layer averaging**: Thickness-weighted averages computed incrementally
- **Early exit**: Stop loop when Ri > Ri_critical criterion met
- **Key insight**: Each iteration computes surface layer average for a different test OBL depth
- **Future optimization**: Integrate `parallelSearchInner` when available (20-40% reduction)

---

## 9. Testing and Validation

### 9.1 Unit Tests
- **Surface layer averaging**: Verify thickness-weighted averages correct for varying layer thickness
- **Surface layer index**: Verify correct determination of surface layer extent (0.1×OBL)
- OBL depth against CVMix reference
- Richardson accumulation logic with cumulative sums
- Profile functions (M1, M2, S1, S2, G)
- Non-local flux conservation

### 9.2 Regression Tests
- **Bit-for-bit with MPAS**: Compare bulk Richardson computation with MPAS-Ocean for identical inputs
- **Ice treatment**: 
  - Verify OBL computed for all iceFraction values 0-1
  - Verify Langmuir disabled for iceFraction ≥ 0.05
  - Verify minimum OBL applied for iceFraction ≥ 0.15
  - Test under land ice (landIceMask = 1)
- Tracer budget verification with non-local flux
- Boundary condition checks

### 9.3 Physical Validation
- **Conservation**: Thickness-weighted averaging conserves mass
- **Dimensional analysis**: All Richardson terms dimensionally consistent
- **LMD94 compliance**: Verify matches Large et al. (1994) equations
- **Two-stage timing**: Verify turbulent scales computed before Richardson calculation

---

## 10. Implementation Roadmap (7 phases, 12-17 hours)

### Phase 1: Core Infrastructure (3-4 hours)
- KPPConstants.h with profile parameters
- OBL depth computation with cumulative sum pattern
- MatchBoth and SimpleShapes schemes
- Unit tests

### Phase 2: Mixing Coefficients (2-3 hours)
- Turbulent scales computation
- KPP profile functions (M1, M2, S1, S2)
- Stability corrections
- VertDiff/VertVisc arrays

### Phase 3: Non-local Flux (1-2 hours)
- G(σ) profile computation
- Reference vs mpas_ocn_tracer_nonlocalflux.F
- OBL boundary constraints

### Phase 4: Wave Enhancement (1-2 hours)
- `cvmix_kpp_EFactor_model()` from wind
- Langmuir enhancement application
- Wind speed integration

### Phase 5: Integration (2-3 hours)
- Config dispatch in VertMix
- Tendencies non-local flux divergence
- OceanRun loop integration
- Config validation

### Phase 6: Testing (2-3 hours)
- Unit tests vs CVMix reference
- Regression with MPAS data
- GPU/CPU performance
- Code review

### Phase 7: Optimization (1-2 hours)
- Profile Stage 1 bottlenecks
- Integrate parallelSearchInner if available
- Memory optimization

---

## 11. Performance Considerations

### 11.1 Computational Complexity
- **Stage 1**: O(nCells × nLevels) with variable work (OPTIMIZE: cumulative sums, early exit)
- **Stage 2**: O(nCells × nLevels) with fixed work per pair
- **Target**: 1.5–2.5× current VertMix cost

### 11.2 GPU Memory
- Main arrays: VertDiff, VertVisc, NonLocalFlux (standard size)
- Temporaries: ~10 scalars per cell for Richardson accumulation
- Estimate: +10% vs current

---

## 12. References

- **Non-local Flux**: [mpas_ocn_tracer_nonlocalflux.F](mpas-ocean/src/shared/mpas_ocn_tracer_nonlocalflux.F)
- **MPAS KPP**: [mpas_ocn_vmix_cvmix.F](mpas-ocean/src/shared/mpas_ocn_vmix_cvmix.F)
- **OMEGA VertMix**: [src/ocn/VertMix.h/cpp](src/ocn/VertMix.h)
- **CVMix**: cvmix.github.io
- **Papers**: Large et al. 1994, Li et al. 2016

---

## 13. Confirmed Decisions ✅

1. **Augment VertMix**: Keep old schemes; add KPP option via config
2. **Theory-based waves only**: Wind → Langmuir enhancement (no active wave model)
   - Automatically disabled when iceFraction ≥ 0.05 or under land ice
3. **Main outputs**: VertDiff, VertVisc, BoundaryLayerDepth, IndexBoundaryLayerDepth
4. **Use Omega Eos**: For density calculations
5. **Location**: src/ocn/ alongside VertMix
6. **OBL Matching**: Support MatchBoth and SimpleShapes schemes
7. **Always compute OBL**: Boundary layer computed under all conditions (open water, sea ice, land ice)
   - Langmuir enhancement disabled under ice but OBL still computed
   - Minimum OBL depth constraint applied when iceFraction > 0.15
8. **Surface Layer Averaging**: Bulk Richardson uses thickness-weighted averages over surface layer (0 to 0.1×OBL)

---

## 14. Open Questions

**parallelSearchInner Integration:**
- If available in your workspace: Should I integrate and use for Stage 1 OBL search?
- Expected benefit: 20-40% reduction in Stage 1 computation

## 15. Physical Correctness Notes

### Bulk Richardson Number
The bulk Richardson number measures stratification **relative to a well-mixed surface layer**, not relative to the surface point. This is physically correct because:
- The surface layer (0 to ε×d, typically 0.1×d) is assumed well-mixed
- Comparing properties at depth d to the surface layer average captures the correct physics
- This matches Large et al. (1994) formulation exactly

### Ice Treatment
Boundary layer depth is computed under all ice conditions because:
- Shear-driven mixing continues under ice (wind stress transmitted through ice)
- Convective mixing continues under ice (surface cooling)
- Only wave-driven processes (Langmuir) are disabled under ice
- Minimum OBL depth under sea ice accounts for altered mixing dynamics
- This approach is validated in MPAS-Ocean and matches observations

---

**Status**: Requirements finalized with physical corrections, ready for Phase 1 implementation  
**Last Updated**: 2026-02-27  
**Physical Review**: Completed - Corrected bulk Richardson (surface layer averaging), ice treatment, and algorithm ordering
