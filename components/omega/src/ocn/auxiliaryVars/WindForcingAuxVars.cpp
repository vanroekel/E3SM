#include "WindForcingAuxVars.h"
#include "DataTypes.h"
#include "Field.h"

#include <limits>

namespace OMEGA {

WindForcingAuxVars::WindForcingAuxVars(const std::string &AuxStateSuffix,
                              const HorzMesh *Mesh,
                              const VertCoord *VCoord)
    : NormalStressEdge("NormalStressEdge" + AuxStateSuffix, Mesh->NEdgesSize),
      ZonalStressCell("WindStressZonal" + AuxStateSuffix, Mesh->NCellsSize),
      MeridStressCell("WindStressMeridional" + AuxStateSuffix,
                      Mesh->NCellsSize),
     TemperaturePistonVelocity("TemperaturePistonVelocity" + AuxStateSuffix,
                        Mesh->NCellsSize),
     SalinityPistonVelocity("SalinityPistonVelocity" + AuxStateSuffix,
                      Mesh->NCellsSize),
     TemperatureSurfaceRestoringValue(
        "TemperatureSurfaceRestoringValue" + AuxStateSuffix,
        Mesh->NCellsSize),
     SalinitySurfaceRestoringValue(
        "SalinitySurfaceRestoringValue" + AuxStateSuffix,
        Mesh->NCellsSize),
     TemperatureInteriorRestoringRate(
        "TemperatureInteriorRestoringRate" + AuxStateSuffix,
        Mesh->NCellsSize, VCoord ? VCoord->NVertLayers : 1),
     SalinityInteriorRestoringRate("SalinityInteriorRestoringRate" + AuxStateSuffix,
                           Mesh->NCellsSize,
                           VCoord ? VCoord->NVertLayers : 1),
     TemperatureInteriorRestoringValue(
        "TemperatureInteriorRestoringValue" + AuxStateSuffix,
        Mesh->NCellsSize, VCoord ? VCoord->NVertLayers : 1),
     SalinityInteriorRestoringValue(
        "SalinityInteriorRestoringValue" + AuxStateSuffix,
        Mesh->NCellsSize, VCoord ? VCoord->NVertLayers : 1),
     LatentHeatFlux("LatentHeatFlux" + AuxStateSuffix, Mesh->NCellsSize),
     SensibleHeatFlux("SensibleHeatFlux" + AuxStateSuffix, Mesh->NCellsSize),
     ShortWaveHeatFlux("ShortWaveHeatFlux" + AuxStateSuffix,
                  Mesh->NCellsSize),
     EvaporationFlux("EvaporationFlux" + AuxStateSuffix, Mesh->NCellsSize),
     RainFlux("RainFlux" + AuxStateSuffix, Mesh->NCellsSize),
     RiverRunoffFlux("RiverRunoffFlux" + AuxStateSuffix, Mesh->NCellsSize),
     IceRunoffFlux("IceRunoffFlux" + AuxStateSuffix, Mesh->NCellsSize),
     SubglacialRunoffFlux("SubglacialRunoffFlux" + AuxStateSuffix,
                     Mesh->NCellsSize),
     IcebergFreshWaterFlux("IcebergFreshWaterFlux" + AuxStateSuffix,
                     Mesh->NCellsSize),
      CellsOnEdge(Mesh->CellsOnEdge), AngleEdge(Mesh->AngleEdge), Interp(Mesh) {
}

void WindForcingAuxVars::registerFields(
    const std::string &AuxGroupName, // name of Auxiliary field group
    const std::string &MeshName      // name of horizontal mesh
) const {

   // Create fields
   const Real FillValue = -9.99e30;
   int NDims            = 1;
   std::vector<std::string> DimNames(NDims);
   std::string DimSuffix;
   if (MeshName == "Default") {
      DimSuffix = "";
   } else {
      DimSuffix = MeshName;
   }

   // Zonal wind stress
   DimNames[0] = "NCells" + DimSuffix;
   auto ZonalStressCellField =
       Field::create(ZonalStressCell.label(),          // field name
                     "zonal wind stress",              // long name/describe
                     "N m^{-2}",                       // units
                     "",                               // CF standard Name
                     std::numeric_limits<Real>::min(), // min valid value
                     std::numeric_limits<Real>::max(), // max valid value
                     FillValue, // scalar for undefined entries
                     NDims,     // number of dimensions
                     DimNames   // dim names
       );
      // Polaris forcing files use lower-camel MPAS-style variable names.
      ZonalStressCellField->addMetadata("InputName", std::string("windStressZonal"));

   // Meridional wind stress
   auto MeridStressCellField =
       Field::create(MeridStressCell.label(),  // field name
                     "meridional wind stress", // long Name or description
                     "N m^{-2}",               // units
                     "",                       // CF standard Name
                     std::numeric_limits<Real>::min(), // min valid value
                     std::numeric_limits<Real>::max(), // max valid value
                     FillValue, // scalar used for undefined entries
                     NDims,     // number of dimensions
                     DimNames   // dimension names
       );
      MeridStressCellField->addMetadata("InputName",
                               std::string("windStressMeridional"));

      auto LatentHeatFluxField =
         Field::create(LatentHeatFlux.label(), "latent heat flux", "W m^{-2}",
                   "", std::numeric_limits<Real>::min(),
                   std::numeric_limits<Real>::max(), FillValue, NDims,
                   DimNames);
      LatentHeatFluxField->addMetadata("InputName", std::string("latentHeatFlux"));

      auto SensibleHeatFluxField =
         Field::create(SensibleHeatFlux.label(), "sensible heat flux",
                   "W m^{-2}", "", std::numeric_limits<Real>::min(),
                   std::numeric_limits<Real>::max(), FillValue, NDims,
                   DimNames);
      SensibleHeatFluxField->addMetadata("InputName",
                                std::string("sensibleHeatFlux"));

      auto ShortWaveHeatFluxField =
         Field::create(ShortWaveHeatFlux.label(), "short-wave heat flux",
                   "W m^{-2}", "", std::numeric_limits<Real>::min(),
                   std::numeric_limits<Real>::max(), FillValue, NDims,
                   DimNames);
      ShortWaveHeatFluxField->addMetadata("InputName",
                                 std::string("shortWaveHeatFlux"));

      auto EvaporationFluxField =
         Field::create(EvaporationFlux.label(), "evaporation freshwater flux",
                   "kg m^{-2} s^{-1}", "", std::numeric_limits<Real>::min(),
                   std::numeric_limits<Real>::max(), FillValue, NDims,
                   DimNames);
      EvaporationFluxField->addMetadata("InputName",
                               std::string("evaporationFlux"));

      auto RainFluxField =
         Field::create(RainFlux.label(), "rain freshwater flux",
                   "kg m^{-2} s^{-1}", "", std::numeric_limits<Real>::min(),
                   std::numeric_limits<Real>::max(), FillValue, NDims,
                   DimNames);
      RainFluxField->addMetadata("InputName", std::string("rainFlux"));

      auto RiverRunoffFluxField =
         Field::create(RiverRunoffFlux.label(), "river runoff freshwater flux",
                   "kg m^{-2} s^{-1}", "", std::numeric_limits<Real>::min(),
                   std::numeric_limits<Real>::max(), FillValue, NDims,
                   DimNames);
      RiverRunoffFluxField->addMetadata("InputName",
                               std::string("riverRunoffFlux"));

      auto IceRunoffFluxField =
         Field::create(IceRunoffFlux.label(), "ice runoff freshwater flux",
                   "kg m^{-2} s^{-1}", "", std::numeric_limits<Real>::min(),
                   std::numeric_limits<Real>::max(), FillValue, NDims,
                   DimNames);
      IceRunoffFluxField->addMetadata("InputName", std::string("iceRunoffFlux"));

      auto SubglacialRunoffFluxField =
         Field::create(SubglacialRunoffFlux.label(),
                   "subglacial runoff freshwater flux",
                   "kg m^{-2} s^{-1}", "", std::numeric_limits<Real>::min(),
                   std::numeric_limits<Real>::max(), FillValue, NDims,
                   DimNames);
      SubglacialRunoffFluxField->addMetadata("InputName",
                                   std::string("subglacialRunoffFlux"));

      auto IcebergFreshWaterFluxField =
         Field::create(IcebergFreshWaterFlux.label(),
                   "iceberg freshwater flux", "kg m^{-2} s^{-1}", "",
                   std::numeric_limits<Real>::min(),
                   std::numeric_limits<Real>::max(), FillValue, NDims,
                   DimNames);
      IcebergFreshWaterFluxField->addMetadata(
         "InputName", std::string("icebergFreshWaterFlux"));

   // Add fields to FieldGroup
   FieldGroup::addFieldToGroup(ZonalStressCell.label(), AuxGroupName);
   FieldGroup::addFieldToGroup(MeridStressCell.label(), AuxGroupName);
   FieldGroup::addFieldToGroup(LatentHeatFlux.label(), AuxGroupName);
   FieldGroup::addFieldToGroup(SensibleHeatFlux.label(), AuxGroupName);
   FieldGroup::addFieldToGroup(ShortWaveHeatFlux.label(), AuxGroupName);
   FieldGroup::addFieldToGroup(EvaporationFlux.label(), AuxGroupName);
   FieldGroup::addFieldToGroup(RainFlux.label(), AuxGroupName);
   FieldGroup::addFieldToGroup(RiverRunoffFlux.label(), AuxGroupName);
   FieldGroup::addFieldToGroup(IceRunoffFlux.label(), AuxGroupName);
   FieldGroup::addFieldToGroup(SubglacialRunoffFlux.label(), AuxGroupName);
   FieldGroup::addFieldToGroup(IcebergFreshWaterFlux.label(), AuxGroupName);

   // Attach data
   ZonalStressCellField->attachData<Array1DReal>(ZonalStressCell);
   MeridStressCellField->attachData<Array1DReal>(MeridStressCell);
   LatentHeatFluxField->attachData<Array1DReal>(LatentHeatFlux);
   SensibleHeatFluxField->attachData<Array1DReal>(SensibleHeatFlux);
   ShortWaveHeatFluxField->attachData<Array1DReal>(ShortWaveHeatFlux);
   EvaporationFluxField->attachData<Array1DReal>(EvaporationFlux);
   RainFluxField->attachData<Array1DReal>(RainFlux);
   RiverRunoffFluxField->attachData<Array1DReal>(RiverRunoffFlux);
   IceRunoffFluxField->attachData<Array1DReal>(IceRunoffFlux);
   SubglacialRunoffFluxField->attachData<Array1DReal>(SubglacialRunoffFlux);
   IcebergFreshWaterFluxField->attachData<Array1DReal>(IcebergFreshWaterFlux);
}

void WindForcingAuxVars::unregisterFields() const {
   Field::destroy(ZonalStressCell.label());
   Field::destroy(MeridStressCell.label());
   Field::destroy(LatentHeatFlux.label());
   Field::destroy(SensibleHeatFlux.label());
   Field::destroy(ShortWaveHeatFlux.label());
   Field::destroy(EvaporationFlux.label());
   Field::destroy(RainFlux.label());
   Field::destroy(RiverRunoffFlux.label());
   Field::destroy(IceRunoffFlux.label());
   Field::destroy(SubglacialRunoffFlux.label());
   Field::destroy(IcebergFreshWaterFlux.label());
}

} // namespace OMEGA
