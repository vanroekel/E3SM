#include "SfcStressForcingVars.h"
#include "DataTypes.h"
#include "Field.h"

#include <limits>

namespace OMEGA {

SfcStressForcingVars::SfcStressForcingVars(const std::string &Suffix,
                                           const HorzMesh *Mesh)
    : NormalStressEdge("NormalStressEdge" + Suffix, Mesh->NEdgesSize),
      ZonalStressCell("SfcStressZonal" + Suffix, Mesh->NCellsSize),
      MeridStressCell("SfcStressMeridional" + Suffix, Mesh->NCellsSize),
      LatentHeatFlux("LatentHeatFlux" + Suffix, Mesh->NCellsSize),
      SensibleHeatFlux("SensibleHeatFlux" + Suffix, Mesh->NCellsSize),
      ShortWaveHeatFlux("ShortWaveHeatFlux" + Suffix, Mesh->NCellsSize),
      EvaporationFlux("EvaporationFlux" + Suffix, Mesh->NCellsSize),
      RainFlux("RainFlux" + Suffix, Mesh->NCellsSize),
      RiverRunoffFlux("RiverRunoffFlux" + Suffix, Mesh->NCellsSize),
      IceRunoffFlux("IceRunoffFlux" + Suffix, Mesh->NCellsSize),
      SubglacialRunoffFlux("SubglacialRunoffFlux" + Suffix, Mesh->NCellsSize),
      IcebergFreshWaterFlux("IcebergFreshWaterFlux" + Suffix, Mesh->NCellsSize),
      CellsOnEdge(Mesh->CellsOnEdge), AngleEdge(Mesh->AngleEdge), Interp(Mesh) {
   deepCopy(NormalStressEdge, 0.0_Real);
   deepCopy(ZonalStressCell, 0.0_Real);
   deepCopy(MeridStressCell, 0.0_Real);
   deepCopy(LatentHeatFlux, 0.0_Real);
   deepCopy(SensibleHeatFlux, 0.0_Real);
   deepCopy(ShortWaveHeatFlux, 0.0_Real);
   deepCopy(EvaporationFlux, 0.0_Real);
   deepCopy(RainFlux, 0.0_Real);
   deepCopy(RiverRunoffFlux, 0.0_Real);
   deepCopy(IceRunoffFlux, 0.0_Real);
   deepCopy(SubglacialRunoffFlux, 0.0_Real);
   deepCopy(IcebergFreshWaterFlux, 0.0_Real);
}

void SfcStressForcingVars::registerFields(
    const std::string &MeshName // name of horizontal mesh
) const {

   int NDims = 1;
   std::vector<std::string> DimNames(NDims);
   std::string DimSuffix;
   if (MeshName == "Default") {
      DimSuffix = "";
   } else {
      DimSuffix = MeshName;
   }

   DimNames[0] = "NCells" + DimSuffix;
   std::vector<std::string> EdgeDimNames(1);
   EdgeDimNames[0] = "NEdges" + DimSuffix;

   auto NormalStressEdgeField =
       Field::create(NormalStressEdge.label(),         // field name
                     "edge-normal surface stress",     // long name/describe
                     "N m^{-2}",                       // units
                     "",                               // CF standard Name
                     std::numeric_limits<Real>::min(), // min valid value
                     std::numeric_limits<Real>::max(), // max valid value
                     1,                                // number of dimensions
                     EdgeDimNames                      // dim names
       );

   auto ZonalStressCellField =
       Field::create(ZonalStressCell.label(),          // field name
                     "zonal surface stress",           // long name/describe
                     "N m^{-2}",                       // units
                     "",                               // CF standard Name
                     std::numeric_limits<Real>::min(), // min valid value
                     std::numeric_limits<Real>::max(), // max valid value
                     NDims,                            // number of dimensions
                     DimNames                          // dim names
       );

   auto MeridStressCellField =
       Field::create(MeridStressCell.label(),     // field name
                     "meridional surface stress", // long Name or description
                     "N m^{-2}",                  // units
                     "",                          // CF standard Name
                     std::numeric_limits<Real>::min(), // min valid value
                     std::numeric_limits<Real>::max(), // max valid value
                     NDims,                            // number of dimensions
                     DimNames                          // dimension names
       );

   auto createCellField = [&](const Array1DReal &FieldData,
                              const std::string &LongName,
                              const std::string &Units) {
      auto CellField =
          Field::create(FieldData.label(),                   // field name
                        LongName,                            // long name
                        Units,                               // units
                        "",                                  // CF standard Name
                        std::numeric_limits<Real>::lowest(), // min valid value
                        std::numeric_limits<Real>::max(),    // max valid value
                        NDims,                               // number of dims
                        DimNames                             // dimension names
          );
      FieldGroup::addFieldToGroup(FieldData.label(), "Forcing");
      CellField->attachData<Array1DReal>(FieldData, false);
   };

   FieldGroup::addFieldToGroup(NormalStressEdge.label(), "Forcing");
   FieldGroup::addFieldToGroup(ZonalStressCell.label(), "Forcing");
   FieldGroup::addFieldToGroup(MeridStressCell.label(), "Forcing");

   NormalStressEdgeField->attachData<Array1DReal>(NormalStressEdge, false);
   ZonalStressCellField->attachData<Array1DReal>(ZonalStressCell, false);
   MeridStressCellField->attachData<Array1DReal>(MeridStressCell, false);

   createCellField(LatentHeatFlux, "latent heat flux", "W m^{-2}");
   createCellField(SensibleHeatFlux, "sensible heat flux", "W m^{-2}");
   createCellField(ShortWaveHeatFlux, "shortwave heat flux", "W m^{-2}");
   createCellField(EvaporationFlux, "evaporation freshwater flux",
                   "kg m^{-2} s^{-1}");
   createCellField(RainFlux, "rain freshwater flux", "kg m^{-2} s^{-1}");
   createCellField(RiverRunoffFlux, "river runoff freshwater flux",
                   "kg m^{-2} s^{-1}");
   createCellField(IceRunoffFlux, "ice runoff freshwater flux",
                   "kg m^{-2} s^{-1}");
   createCellField(SubglacialRunoffFlux, "subglacial runoff freshwater flux",
                   "kg m^{-2} s^{-1}");
   createCellField(IcebergFreshWaterFlux, "iceberg freshwater flux",
                   "kg m^{-2} s^{-1}");
}

void SfcStressForcingVars::unregisterFields() const {
   auto destroyIfExists = [](const std::string &FieldName) {
      if (Field::exists(FieldName)) {
         Field::destroy(FieldName);
      }
   };

   destroyIfExists(NormalStressEdge.label());
   destroyIfExists(ZonalStressCell.label());
   destroyIfExists(MeridStressCell.label());
   destroyIfExists(LatentHeatFlux.label());
   destroyIfExists(SensibleHeatFlux.label());
   destroyIfExists(ShortWaveHeatFlux.label());
   destroyIfExists(EvaporationFlux.label());
   destroyIfExists(RainFlux.label());
   destroyIfExists(RiverRunoffFlux.label());
   destroyIfExists(IceRunoffFlux.label());
   destroyIfExists(SubglacialRunoffFlux.label());
   destroyIfExists(IcebergFreshWaterFlux.label());
}

} // namespace OMEGA
