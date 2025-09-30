/**
 * @file IRTCVisualization.cpp
 * @brief Source of the Boussinesq rotating thermal convection in a spherical
 * shell (Toroidal/Poloidal formulation)
 */

// System includes
//

// Project includes
//
#include "Model/Boussinesq/Shell/RTC/IRTCVisualization.hpp"
#include "Model/Boussinesq/Shell/RTC/gitHash.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Generator/Visualizers/ScalarFieldVisualizer.hpp"
#include "QuICC/Generator/Visualizers/SphericalVerticalFieldVisualizer.hpp"
#include "QuICC/Generator/Visualizers/VectorFieldVisualizer.hpp"
#include "QuICC/Io/Variable/VisualizationFileWriter.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/PhysicalNames/VelocityZ.hpp"
#include "QuICC/PhysicalNames/VorticityZ.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Shell {

namespace RTC {

VectorFormulation::Id IRTCVisualization::SchemeFormulation()
{
   return VectorFormulation::TORPOL;
}

std::string IRTCVisualization::version() const
{
   return std::string(gitHash);
}

void IRTCVisualization::addVisualizers(SharedVisualizationGenerator spVis)
{
   // Shared pointer to basic field visualizer
   Equations::SharedScalarFieldVisualizer spScalar;
   Equations::SharedVectorFieldVisualizer spVector;
   Equations::SharedSphericalVerticalFieldVisualizer spVertical;

   // Add temperature field visualization
   spScalar =
      spVis->addEquation<Equations::ScalarFieldVisualizer>(this->spBackend());
   spScalar->setFields(true, false);
   spScalar->setIdentity(PhysicalNames::Temperature::id());

   // Add velocity field visualization
   spVector =
      spVis->addEquation<Equations::VectorFieldVisualizer>(this->spBackend());
   spVector->setFields(true, false, false);
   spVector->setIdentity(PhysicalNames::Velocity::id());

   // Add vertical velocity visualization
   spVertical = spVis->addEquation<Equations::SphericalVerticalFieldVisualizer>(
      this->spBackend());
   spVertical->setFieldType(FieldType::VECTOR);
   spVertical->setIdentity(PhysicalNames::VelocityZ::id(),
      PhysicalNames::Velocity::id());

   // Add vertical vorticity visualization
   spVertical = spVis->addEquation<Equations::SphericalVerticalFieldVisualizer>(
      this->spBackend());
   spVertical->setFieldType(FieldType::CURL);
   spVertical->setIdentity(PhysicalNames::VorticityZ::id(),
      PhysicalNames::Velocity::id());

   // Add output file
   auto spOut = std::make_shared<Io::Variable::VisualizationFileWriter>(
      spVis->ss().tag());
   spOut->expect(PhysicalNames::Temperature::id());
   spOut->expect(PhysicalNames::Velocity::id());
   spOut->expect(PhysicalNames::VelocityZ::id());
   spOut->expect(PhysicalNames::VorticityZ::id());
   spVis->addHdf5OutputFile(spOut);
}

} // namespace RTC
} // namespace Shell
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC
