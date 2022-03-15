/**
 * @file IRTCModel.cpp
 * @brief Source of the Boussinesq rotating thermal convection in a spherical shell (Toroidal/Poloidal formulation)
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Model/Boussinesq/Shell/RTC/IRTCModel.hpp"

// Project includes
//
#include "QuICC/Model/Boussinesq/Shell/RTC/Transport.hpp"
#include "QuICC/Model/Boussinesq/Shell/RTC/Momentum.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/PhysicalNames/VelocityZ.hpp"
#include "QuICC/PhysicalNames/VorticityZ.hpp"
#include "QuICC/NonDimensional/Ekman.hpp"
#include "QuICC/NonDimensional/Rayleigh.hpp"
#include "QuICC/NonDimensional/Prandtl.hpp"
#include "QuICC/NonDimensional/RRatio.hpp"
#include "QuICC/NonDimensional/Heating.hpp"
#include "QuICC/NonDimensional/CflInertial.hpp"
#include "QuICC/NonDimensional/Lower1D.hpp"
#include "QuICC/NonDimensional/Upper1D.hpp"
#include "QuICC/Io/Variable/StateFileReader.hpp"
#include "QuICC/Io/Variable/StateFileWriter.hpp"
#include "QuICC/Io/Variable/VisualizationFileWriter.hpp"
#include "QuICC/Io/Variable/ShellScalarEnergyWriter.hpp"
#include "QuICC/Io/Variable/ShellScalarLSpectrumWriter.hpp"
#include "QuICC/Io/Variable/ShellScalarMSpectrumWriter.hpp"
#include "QuICC/Io/Variable/ShellTorPolEnergyWriter.hpp"
#include "QuICC/Io/Variable/ShellTorPolLSpectrumWriter.hpp"
#include "QuICC/Io/Variable/ShellTorPolMSpectrumWriter.hpp"
#include "QuICC/Generator/States/RandomScalarState.hpp"
#include "QuICC/Generator/States/RandomVectorState.hpp"
#include "QuICC/Generator/States/ShellExactScalarState.hpp"
#include "QuICC/Generator/States/ShellExactVectorState.hpp"
#include "QuICC/Generator/States/Kernels/Shell/BenchmarkTempC1.hpp"
#include "QuICC/Generator/Visualizers/ScalarFieldVisualizer.hpp"
#include "QuICC/Generator/Visualizers/VectorFieldVisualizer.hpp"
#include "QuICC/Generator/Visualizers/SphericalVerticalFieldVisualizer.hpp"
#include "QuICC/SpectralKernels/MakeRandom.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Shell {

namespace RTC {

   VectorFormulation::Id IRTCModel::SchemeFormulation()
   {
      return VectorFormulation::TORPOL;
   }

   void IRTCModel::addEquations(SharedSimulation spSim)
   {
      // Add transport equation
      spSim->addEquation<Equations::Boussinesq::Shell::RTC::Transport>(this->spBackend());

      // Add Navier-Stokes equation
      spSim->addEquation<Equations::Boussinesq::Shell::RTC::Momentum>(this->spBackend());
   }

   void IRTCModel::addStates(SharedStateGenerator spGen)
   {
      // Shared pointer to equation
      Equations::SharedShellExactScalarState spScalar;
      Equations::SharedShellExactVectorState spVector;

      Spectral::Kernel::Complex3DMapType tSH;
      std::pair<Spectral::Kernel::Complex3DMapType::iterator,bool> ptSH;

      // Add temperature initial state generator
      spScalar = spGen->addEquation<Equations::ShellExactScalarState>(this->spBackend());
      spScalar->setIdentity(PhysicalNames::Temperature::id());
      switch(3)
      {
         case 0:
            {
               spScalar->setPhysicalNoise(1e-15);
            }
            break;

         case 1:
            {
               spScalar->setPhysicalConstant(1.0);
            }
            break;

         case 2:
            {
               tSH.clear();
               ptSH = tSH.insert(std::make_pair(std::make_pair(3,3), std::map<int,MHDComplex>()));
               ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0,2.0)));
               spScalar->setSpectralModes(tSH);
            }
            break;

         case 3:
            {
               auto ri = spScalar->eqParams().nd(NonDimensional::Lower1D::id());
               auto ro = spScalar->eqParams().nd(NonDimensional::Upper1D::id());
               auto spKernel = std::make_shared<Physical::Kernel::Shell::BenchmarkTempC1>();
               spKernel->init(ri, ro);
               spScalar->setPhysicalKernel(spKernel);
            }
            break;

         case 4:
            {
               auto spKernel = std::make_shared<Spectral::Kernel::MakeRandom>(spGen->ss().has(SpatialScheme::Feature::ComplexSpectrum));
               std::vector<MHDFloat> ratios = {1e4, 1e4, 1e4};
               spKernel->setRatio(ratios);
               spKernel->init(-1e-4, 1e-4);
               spScalar->setSrcKernel(spKernel);
            }
            break;
      }

      // Add velocity initial state generator
      spVector = spGen->addEquation<Equations::ShellExactVectorState>(this->spBackend());
      spVector->setIdentity(PhysicalNames::Velocity::id());
      switch(3)
      {
         // Toroidal only
         case 0:
            {
               // Toroidal
               tSH.clear();
               ptSH = tSH.insert(std::make_pair(std::make_pair(1,1), std::map<int,MHDComplex>()));
               ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0)));
               spVector->setSpectralModes(FieldComponents::Spectral::TOR, tSH);
               // Poloidal
               tSH.clear();
               spVector->setSpectralModes(FieldComponents::Spectral::POL, tSH);
            }
            break;

         // Poloidal only
         case 1:
            {
               // Toroidal
               tSH.clear();
               spVector->setSpectralModes(FieldComponents::Spectral::TOR, tSH);
               // Poloidal
               tSH.clear();
               ptSH = tSH.insert(std::make_pair(std::make_pair(2,0), std::map<int,MHDComplex>()));
               ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0)));
               spVector->setSpectralModes(FieldComponents::Spectral::POL, tSH);
            }
            break;

         // Toroidal & Poloidal
         case 2:
            {
               // Toroidal
               tSH.clear();
               ptSH = tSH.insert(std::make_pair(std::make_pair(1,1), std::map<int,MHDComplex>()));
               ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0)));
               spVector->setSpectralModes(FieldComponents::Spectral::TOR, tSH);
               // Poloidal
               tSH.clear();
               ptSH = tSH.insert(std::make_pair(std::make_pair(2,0), std::map<int,MHDComplex>()));
               ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0)));
               spVector->setSpectralModes(FieldComponents::Spectral::POL, tSH);
            }
            break;

         case 3:
            {
               auto spKernel = std::make_shared<Spectral::Kernel::MakeRandom>(spGen->ss().has(SpatialScheme::Feature::ComplexSpectrum));
               std::vector<MHDFloat> ratios = {1e2, 1e2, 1e2};
               spKernel->setRatio(ratios);
               spKernel->init(-1e-15, 1e-15);
               spVector->setSrcKernel(FieldComponents::Spectral::TOR, spKernel);
               spVector->setSrcKernel(FieldComponents::Spectral::POL, spKernel);
            }
            break;

         case 4:
            {
               auto spKernel = std::make_shared<Spectral::Kernel::MakeRandom>(spGen->ss().has(SpatialScheme::Feature::ComplexSpectrum));
               std::vector<MHDFloat> ratios = {1e4, 1e4, 1e4};
               spKernel->setRatio(ratios);
               spKernel->init(-1e-4, 1e-4);
               spVector->setSrcKernel(FieldComponents::Spectral::TOR, spKernel);
               spVector->setSrcKernel(FieldComponents::Spectral::POL, spKernel);
            }
            break;
      }

      // Add output file
      auto spOut = std::make_shared<Io::Variable::StateFileWriter>(spGen->ss().tag(), spGen->ss().has(SpatialScheme::Feature::RegularSpectrum));
      spOut->expect(PhysicalNames::Temperature::id());
      spOut->expect(PhysicalNames::Velocity::id());
      spGen->addHdf5OutputFile(spOut);
   }

   void IRTCModel::addVisualizers(SharedVisualizationGenerator spVis)
   {
      // Shared pointer to basic field visualizer
      Equations::SharedScalarFieldVisualizer spScalar;
      Equations::SharedVectorFieldVisualizer spVector;
      Equations::SharedSphericalVerticalFieldVisualizer spVertical;

      // Add temperature field visualization
      spScalar = spVis->addEquation<Equations::ScalarFieldVisualizer>(this->spBackend());
      spScalar->setFields(true, false);
      spScalar->setIdentity(PhysicalNames::Temperature::id());

      // Add velocity field visualization
      spVector = spVis->addEquation<Equations::VectorFieldVisualizer>(this->spBackend());
      spVector->setFields(true, false, false);
      spVector->setIdentity(PhysicalNames::Velocity::id());

      // Add vertical velocity visualization
      spVertical = spVis->addEquation<Equations::SphericalVerticalFieldVisualizer>(this->spBackend());
      spVertical->setFieldType(FieldType::VECTOR);
      spVertical->setIdentity(PhysicalNames::VelocityZ::id(), PhysicalNames::Velocity::id());

      // Add vertical vorticity visualization
      spVertical = spVis->addEquation<Equations::SphericalVerticalFieldVisualizer>(this->spBackend());
      spVertical->setFieldType(FieldType::CURL);
      spVertical->setIdentity(PhysicalNames::VorticityZ::id(), PhysicalNames::Velocity::id());

      // Add output file
      auto spOut = std::make_shared<Io::Variable::VisualizationFileWriter>(spVis->ss().tag());
      spOut->expect(PhysicalNames::Temperature::id());
      spOut->expect(PhysicalNames::Velocity::id());
      spOut->expect(PhysicalNames::VelocityZ::id());
      spOut->expect(PhysicalNames::VorticityZ::id());
      spVis->addHdf5OutputFile(spOut);
   }

   void IRTCModel::addAsciiOutputFiles(SharedSimulation spSim)
   {
      // Create temperature energy writer
      auto spTemp = std::make_shared<Io::Variable::ShellScalarEnergyWriter>("temperature", spSim->ss().tag());
      spTemp->expect(PhysicalNames::Temperature::id());
      spSim->addAsciiOutputFile(spTemp);

#if 1
      // Create temperature L energy spectrum writer
      auto spTempL = std::make_shared<Io::Variable::ShellScalarLSpectrumWriter>("temperature", spSim->ss().tag());
      spTempL->expect(PhysicalNames::Temperature::id());
      spTempL->numberOutput();
      spTempL->onlyEvery(10);
      spSim->addAsciiOutputFile(spTempL);

      // Create temperature M energy spectrum writer
      auto spTempM = std::make_shared<Io::Variable::ShellScalarMSpectrumWriter>("temperature", spSim->ss().tag());
      spTempM->expect(PhysicalNames::Temperature::id());
      spTempM->numberOutput();
      spTempM->onlyEvery(10);
      spSim->addAsciiOutputFile(spTempM);
#endif

      // Create kinetic energy writer
      auto spKinetic = std::make_shared<Io::Variable::ShellTorPolEnergyWriter>("kinetic", spSim->ss().tag());
      spKinetic->expect(PhysicalNames::Velocity::id());
      spSim->addAsciiOutputFile(spKinetic);

#if 1
      // Create kinetic L energy spectrum writer
      auto spKineticL = std::make_shared<Io::Variable::ShellTorPolLSpectrumWriter>("kinetic", spSim->ss().tag());
      spKineticL->expect(PhysicalNames::Velocity::id());
      spKineticL->numberOutput();
      spKineticL->onlyEvery(10);
      spSim->addAsciiOutputFile(spKineticL);

      // Create kinetic M energy spectrum writer
      auto spKineticM = std::make_shared<Io::Variable::ShellTorPolMSpectrumWriter>("kinetic", spSim->ss().tag());
      spKineticM->expect(PhysicalNames::Velocity::id());
      spKineticM->numberOutput();
      spKineticM->onlyEvery(10);
      spSim->addAsciiOutputFile(spKineticM);
#endif
   }

}
}
}
}
}
