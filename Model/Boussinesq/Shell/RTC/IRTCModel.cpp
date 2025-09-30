/**
 * @file IRTCModel.cpp
 * @brief Source of the Boussinesq rotating thermal convection in a spherical
 * shell (Toroidal/Poloidal formulation)
 */

// System includes
//

// Project includes
//
#include "Model/Boussinesq/Shell/RTC/IRTCModel.hpp"
#include "Model/Boussinesq/Shell/RTC/Momentum.hpp"
#include "Model/Boussinesq/Shell/RTC/Transport.hpp"
#include "Model/Boussinesq/Shell/RTC/gitHash.hpp"
#include "QuICC/Io/Variable/FieldProbeWriter.hpp"
#include "QuICC/Io/Variable/ShellNusseltWriter.hpp"
#include "QuICC/Io/Variable/ShellScalarEnergyWriter.hpp"
#include "QuICC/Io/Variable/ShellScalarLSpectrumWriter.hpp"
#include "QuICC/Io/Variable/ShellScalarMSpectrumWriter.hpp"
#include "QuICC/Io/Variable/ShellTorPolEnergyWriter.hpp"
#include "QuICC/Io/Variable/ShellTorPolLSpectrumWriter.hpp"
#include "QuICC/Io/Variable/ShellTorPolMSpectrumWriter.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Shell {

namespace RTC {

VectorFormulation::Id IRTCModel::SchemeFormulation()
{
   return VectorFormulation::TORPOL;
}

std::string IRTCModel::version() const
{
   return std::string(gitHash);
}

void IRTCModel::addEquations(SharedSimulation spSim)
{
   // Add transport equation
   spSim->addEquation<Equations::Boussinesq::Shell::RTC::Transport>(
      this->spBackend());

   // Add Navier-Stokes equation
   spSim->addEquation<Equations::Boussinesq::Shell::RTC::Momentum>(
      this->spBackend());
}

std::map<std::string, std::map<std::string, int>> IRTCModel::configTags() const
{
   std::map<std::string, int> onOff;
   onOff.emplace("enable", 1);

   std::map<std::string, int> options;
   options.emplace("enable", 0);
   options.emplace("numbered", 0);
   options.emplace("only_every", 1);

   std::map<std::string, std::map<std::string, int>> tags;
   // kinetic
   tags.emplace("kinetic_energy", onOff);
   tags.emplace("kinetic_l_spectrum", options);
   tags.emplace("kinetic_m_spectrum", options);
   // temperature
   tags.emplace("temperature_energy", onOff);
   tags.emplace("temperature_l_spectrum", options);
   tags.emplace("temperature_m_spectrum", options);
   tags.emplace("temperature_nusselt", onOff);

   return tags;
}

void IRTCModel::addAsciiOutputFiles(SharedSimulation spSim)
{
   // Create temperature energy writer
   this->enableAsciiFile<Io::Variable::ShellScalarEnergyWriter>(
      "temperature_energy", "temperature", PhysicalNames::Temperature::id(),
      spSim);

   // Create temperature L energy spectrum writer
   this->enableAsciiFile<Io::Variable::ShellScalarLSpectrumWriter>(
      "temperature_l_spectrum", "temperature", PhysicalNames::Temperature::id(),
      spSim);

   // Create temperature M energy spectrum writer
   this->enableAsciiFile<Io::Variable::ShellScalarMSpectrumWriter>(
      "temperature_m_spectrum", "temperature", PhysicalNames::Temperature::id(),
      spSim);

   // Create kinetic energy writer
   this->enableAsciiFile<Io::Variable::ShellTorPolEnergyWriter>(
      "kinetic_energy", "kinetic", PhysicalNames::Velocity::id(), spSim);

   // Create kinetic L energy spectrum writer
   this->enableAsciiFile<Io::Variable::ShellTorPolLSpectrumWriter>(
      "kinetic_l_spectrum", "kinetic", PhysicalNames::Velocity::id(), spSim);

   // Create kinetic M energy spectrum writer
   this->enableAsciiFile<Io::Variable::ShellTorPolMSpectrumWriter>(
      "kinetic_m_spectrum", "kinetic", PhysicalNames::Velocity::id(), spSim);

   // Create nusselt number writer
   this->enableAsciiFile<Io::Variable::ShellNusseltWriter>(
      "temperature_nusselt", "temperature_", PhysicalNames::Temperature::id(),
      spSim);

   // Examples of physical space field probes
   //

   const bool probeVelocity = false;
   const bool probeTemperature = false;

   // Add Velocity probe
   if (probeVelocity)
   {
      const std::vector<MHDFloat> pos = {1.10, Math::PI / 2.0, 0};
      auto spFile = std::make_shared<Io::Variable::FieldProbeWriter>(
         "velocity_", spSim->ss().tag(), pos);
      spFile->expect(PhysicalNames::Velocity::id());
      spSim->addAsciiOutputFile(spFile);
   }

   // Add Temperature probe
   if (probeTemperature)
   {
      const std::vector<MHDFloat> pos = {1.10, Math::PI / 2.0, 0};
      auto spFile = std::make_shared<Io::Variable::FieldProbeWriter>(
         "temperature_", spSim->ss().tag(), pos);
      spFile->expect(PhysicalNames::Temperature::id());
      spSim->addAsciiOutputFile(spFile);
   }
}

} // namespace RTC
} // namespace Shell
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC
