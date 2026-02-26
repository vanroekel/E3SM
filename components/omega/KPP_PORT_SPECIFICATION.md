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
1. Compute surface friction u* and buoyancy flux B_0
2. Compute Langmuir enhancement factor from wind speed
3. **Loop through layers computing bulk Richardson**:
   - Use **cumulative sum pattern** to avoid redundant calculation
   - Stop when `Ri_b > Ri_crit` (typically 0.3)
   - ⚠️ **Performance Note**: Early exit critical; inner loop has growing complexity
4. Optional: Apply MatchBoth interpolation scheme
5. Apply min/max OBL constraints
6. Store OBL depth and layer index

**Stage 2: Mixing Coefficients** (~O(nCells × nLevels))
For each cell:
1. Compute turbulent velocity scales: w_s
2. Apply KPP profile within OBL: ν(z) = u* × w_s × w(σ) + ν_bg
3. Compute non-local flux G(σ) within OBL
4. Store background mixing below OBL

### 5.2 Non-local Flux Details

**Vertical Divergence Form:**
`d(T)/dt = -d/dz(G(σ) × Q_surf)`

**Omega Implementation:**
1. KPP computes `VertNonLocalFlux[nCells][nLevels+1]`
2. In Tendencies: Loop over tracers, apply: `tend += (G(k) - G(k+1)) × Q_surf`
3. Zero boundary conditions at surface and below OBL

See [mpas_ocn_tracer_nonlocalflux.F](mpas-ocean/src/shared/mpas_ocn_tracer_nonlocalflux.F) for reference implementation.

### 5.3 Key Quantities

| Variable | Units | Range | Notes |
|----------|-------|-------|-------|
| u\* | m/s | 0.0001+ | Friction velocity (no upper limit) |
| B₀ | m²/s³ | -0.5 to 0.01 | Buoyancy flux |
| h_OBL | m | 5–500 | Boundary layer depth |
| Ri_b | — | 0–1 | Bulk Richardson criterion |
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
    Minimum_OBL_under_sea_ice: 5.0         # (m)

    # Wave enhancement
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

⚠️ **Critical Performance Issue**: OBL bulk Richardson loop has growing inner complexity
- **Problem**: Computing shear and buoyancy accumulation can lead to O(nLevels²) behavior
- **Solution**: Use **cumulative sum pattern** - maintain running totals, don't recompute
- **Early exit**: Stop loop when Ri > Ri_critical criterion met
- **Future optimization**: Integrate `parallelSearchInner` when available (20-40% reduction)

---

## 9. Testing and Validation

### 9.1 Unit Tests
- OBL depth against CVMix reference
- Richardson accumulation logic
- Profile functions (M1, M2, S1, S2, G)
- Non-local flux conservation

### 9.2 Regression Tests
- Compare with MPAS-Ocean KPP for identical inputs
- Tracer budget verification with non-local flux
- Boundary condition checks

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
3. **Main outputs**: VertDiff, VertVisc, BoundaryLayerDepth, IndexBoundaryLayerDepth
4. **Use Omega Eos**: For density calculations
5. **Location**: src/ocn/ alongside VertMix
6. **OBL Matching**: Support MatchBoth and SimpleShapes schemes
7. **No Fixed OBL**: Always compute dynamically (no fixed boundary layer option)

---

## 14. Open Questions

**parallelSearchInner Integration:**
- If available in your workspace: Should I integrate and use for Stage 1 OBL search?
- Expected benefit: 20-40% reduction in Stage 1 computation

---

**Status**: Requirements finalized, ready for Phase 1 implementation
**Last Updated**: 2026-02-26
