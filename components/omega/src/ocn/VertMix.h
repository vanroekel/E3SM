#ifndef OMEGA_VERTMIX_H
#define OMEGA_VERTMIX_H
//===-- ocn/VertMix.cpp - Vertical Mixing Coefficients -----------*- C++
//-*-===//
//
/// \file
/// \brief Contains functors for calculating vertical diffusivity and viscosity
///
/// This header defines functors to be called by the time-stepping scheme
/// to calculate the vertical diffusivity and viscosity
//
//===----------------------------------------------------------------------===//

#include "AuxiliaryState.h"
#include "Config.h"
#include "HorzMesh.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"
#include "TimeMgr.h"
#include "VertCoord.h"
#include <string>

namespace OMEGA {

class ConvectiveMix {
 public:
   bool Enabled = true; ///< Enable convective mixing flag

   // Convective mixing parameters
   Real ConvDiff =
       1.0; ///< Convective vertical viscosity and diffusivity (m^2 s^-1)
   Real ConvTriggerBVF = 0.0; ///< Reference density (kg m^-3) at (T,S)=(0,0)

   /// Constructor for ConvectiveMix
   ConvectiveMix(const VertCoord *VCoord);

   KOKKOS_FUNCTION void
   operator()(Array2DReal VertDiff, Array2DReal VertVisc, I4 ICell, I4 KChunk,
              const Array2DReal &BruntVaisalaFreqSq) const {

      const I4 KStart = chunkStart(KChunk, MinLayerCell(ICell) + 1);
      const I4 KLen   = chunkLength(KChunk, KStart, MaxLayerCell(ICell));

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const I4 K = KStart + KVec;

         if (BruntVaisalaFreqSq(ICell, K) < ConvTriggerBVF) {
            VertDiff(ICell, K) += ConvDiff;
            VertVisc(ICell, K) += ConvDiff;
         }
      }
   }

 private:
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
};

class ShearMix {
 public:
   bool Enabled = true; ///< Enable shear mixing flag

   // Shear mixing parameters
   Real BaseShearValue = 0.005; ///< Base shear vertical viscosity and
                                ///< diffusivity (m^2 s^-1) of LMG94
   Real ShearRiCrit = 0.7;      ///< Critical Richardson number of LMG94
   Real ShearExponent =
       3.0; /// Exponent value used interior shear mixing calculation of LMG94
   I4 RiSmoothLoops = 2; ///< Number of smoothing loops for Richardson number

   /// Constructor for ShearMix
   ShearMix(const VertCoord *VCoord);

   KOKKOS_FUNCTION void
   operator()(Array2DReal VertDiff, Array2DReal VertVisc, I4 ICell, I4 KChunk,
              const Array2DReal &GradRichNumSmoothed) const {

      const I4 KStart = chunkStart(KChunk, MinLayerCell(ICell) + 1);
      const I4 KLen   = chunkLength(KChunk, KStart, MaxLayerCell(ICell));

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const I4 K = KStart + KVec;

         if (GradRichNumSmoothed(ICell, K) <= 0.0_Real) {
            VertDiff(ICell, K) += BaseShearValue;
            VertVisc(ICell, K) += BaseShearValue;
         } else if (GradRichNumSmoothed(ICell, K) > 0.0_Real &&
                    GradRichNumSmoothed(ICell, K) < ShearRiCrit) {
            VertDiff(ICell, K) +=
                Kokkos::pow(
                    1.0_Real -
                        (GradRichNumSmoothed(ICell, K) / ShearRiCrit) *
                            (GradRichNumSmoothed(ICell, K) / ShearRiCrit),
                    ShearExponent) *
                BaseShearValue;
            VertVisc(ICell, K) +=
                Kokkos::pow(
                    1.0_Real -
                        (GradRichNumSmoothed(ICell, K) / ShearRiCrit) *
                            (GradRichNumSmoothed(ICell, K) / ShearRiCrit),
                    ShearExponent) *
                BaseShearValue;
         }
      }
   }

 private:
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
};

/// Class for Gradient Richardson Number calculation
class GradRichardsonNum {
 public:
   /// constructor declaration
   GradRichardsonNum(const HorzMesh *Mesh, const VertCoord *VCoord);

   //   The functor takes the full arrays of Richardson number (inout),
   //   the index ICell, and normal and tangential velocities as inputs,
   //   and outputs the Richardson number.
   KOKKOS_FUNCTION void
   operator()(Array2DReal GradRichNum, I4 ICell, I4 KChunk,
              const Array2DReal &NormalVelocity,
              const Array2DReal &TangentialVelocity,
              const Array2DReal &BruntVaisalaFreqSq) const {

      const I4 KStart = chunkStart(KChunk, MinLayerCell(ICell) + 1);
      const I4 KLen   = chunkLength(KChunk, KStart, MaxLayerCell(ICell));

      Real GradRichNumNorm[VecLength];
      Real GradRichNumTmp[VecLength];

      for (int KVec = 0; KVec < KLen; ++KVec) {
         GradRichNumNorm[KVec] = 1.0e-12_Real;
         GradRichNumTmp[KVec]  = 100.0_Real;
      }

      for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
         I4 JEdge = EdgesOnCell(ICell, J);
         I4 JCell = CellsOnCell(ICell, J);

         for (int KVec = 0; KVec < KLen; ++KVec) {
            const I4 K = KStart + KVec;
            I4 K1      = K - 1;
            I4 K2      = K;

            // Skip this edge contribution if it would access
            // invalid edge velocity levels.
            if (K1 > MaxLayerEdgeBot(JEdge) || K2 > MaxLayerEdgeBot(JEdge))
               continue;

            Real DNormVel =
                NormalVelocity(JEdge, K1) - NormalVelocity(JEdge, K2);
            Real DTanVel =
                TangentialVelocity(JEdge, K1) - TangentialVelocity(JEdge, K2);
            Real DzEdge =
                0.5_Real * (GeomZMid(ICell, K1) + GeomZMid(JCell, K1) -
                            (GeomZMid(ICell, K2) + GeomZMid(JCell, K2)));
            Real ShearSquared =
                (DNormVel * DNormVel + DTanVel * DTanVel) / (DzEdge * DzEdge);
            Real RiEdge =
                Kokkos::max(0.0_Real,
                            0.5_Real * (BruntVaisalaFreqSq(ICell, K2) +
                                        BruntVaisalaFreqSq(JCell, K2))) /
                (ShearSquared + 1.0e-12_Real);

            Real Weight           = 0.25_Real * DcEdge(JEdge) * DvEdge(JEdge);
            GradRichNumNorm[KVec] = GradRichNumNorm[KVec] + Weight;
            GradRichNumTmp[KVec]  = GradRichNumTmp[KVec] + Weight * RiEdge;
         }
      }

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const I4 K            = KStart + KVec;
         GradRichNum(ICell, K) = GradRichNumTmp[KVec] / GradRichNumNorm[KVec];
      }
   }

 private:
   Array2DReal GeomZMid;
   Array2DI4 EdgesOnCell;
   Array2DI4 CellsOnCell;
   Array2DI4 CellsOnEdge;
   Array1DI4 NEdgesOnCell;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
   Array1DI4 MaxLayerEdgeBot;
   Array1DReal DcEdge;
   Array1DReal DvEdge;
   I4 NVertLayers;
   I4 NCellsAll;
};

/// Class for Gradient Richardson Number calculation
class OneTwoOneFilter {
 public:
   /// constructor declaration
   OneTwoOneFilter(const VertCoord *VCoord);
   //   The functor takes the full arrays of Richardson number (inout),
   //   the index ICell, and normal and tangential velocities as inputs,
   //   and outputs the Richardson number.
   KOKKOS_FUNCTION void operator()(Array2DReal VarOut, I4 ICell, I4 KChunk,
                                   const Array2DReal &VarIn) const {

      const I4 KStart = chunkStart(KChunk, MinLayerCell(ICell));
      const I4 KLen   = chunkLength(KChunk, KStart, MaxLayerCell(ICell) + 1);

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const I4 K = KStart + KVec;
         if (K > MinLayerCell(ICell) && K < MaxLayerCell(ICell)) {
            // apply 1-2-1 filter
            VarOut(ICell, K) =
                (VarIn(ICell, K - 1) + 2.0_Real * VarIn(ICell, K) +
                 VarIn(ICell, K + 1)) /
                4.0_Real;
         } else {
            VarOut(ICell, K) = VarIn(ICell, K);
         }
      }
   }

 private:
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
};

/// Class for Vertical Mixing Coefficient (VertMix) calculations
class VertMix {
 public:
   /// Get instance of VertMix
   static VertMix *getInstance();

   /// Destroy instance (frees Kokkos views)
   static void destroyInstance();

   Array2DReal VertDiff;    ///< Vertical diffusivity field (m^2 s^-1)
   Array2DReal VertVisc;    ///< Vertical viscosity field (m^2 s^-1)
   Array2DReal GradRichNum; ///< Gradient Richardson number field
   Array2DReal
       GradRichNumSmoothed; ///< Smoothed Gradient Richardson number field

   std::string VertDiffFldName; ///< Field name for vertical diffusivity
   std::string VertViscFldName; ///< Field name for vertical viscosity
   std::string
       GradRichNumFldName; ///< Field name for gradient Richardson number
   std::string GradRichNumSmoothedFldName; ///< Field name for smoothed gradient
                                           ///< Richardson number
   std::string VertMixGroupName;           ///< VertMix group name (for config)
   std::string Name;                       ///< Name of this VertMix instance

   // Background mixing parameters
   Real BackDiff = 1.0e-5; ///< Background vertical diffusivity (m^2 s^-1)
   Real BackVisc = 1.0e-4; ///< Background vertical viscosity (m^2 s^-1)

   ConvectiveMix
       ComputeVertMixConv;       ///< Functor for Convective VertMix calculation
   ShearMix ComputeVertMixShear; ///< Functor for Shear VertMix calculation
   GradRichardsonNum
       ComputeGradRichardsonNum; ///< Functor for Gradient Richardson Number
                                 ///< calculation
   OneTwoOneFilter ComputeOneTwoOneFilter; ///< Functor for 1-2-1 filtering

   /// Compute vertical diffusivity and viscosity for all cells/layers
   void computeVertMix(const Array2DReal &NormalVelocity,
                       const Array2DReal &TangentialVelocity,
                       const Array2DReal &BruntVaisalaFreqSq);

   /// Initialize VertMix from config and mesh
   static void init();

 private:
   /// Private constructor
   VertMix(const std::string &Name, const HorzMesh *Mesh,
           const VertCoord *VCoord);

   /// Private destructor
   ~VertMix();

   static VertMix *Instance; ///< Instance pointer

   // Delete copy and move constructors and assignment operators
   VertMix(const VertMix &)            = delete;
   VertMix &operator=(const VertMix &) = delete;
   VertMix(VertMix &&)                 = delete;
   VertMix &operator=(VertMix &&)      = delete;

   const HorzMesh *Mesh;    ///< Horizontal mesh
   const VertCoord *VCoord; ///< Vertical coordinate

   // Define fields and metadata
   void defineFields();

}; // End class VertMix

} // namespace OMEGA
#endif
