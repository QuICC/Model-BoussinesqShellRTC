/**
 * @file IRTCBackend.cpp
 * @brief Source of the interface for model backend
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "QuICC/Model/Boussinesq/Shell/RTC/IRTCBackend.hpp"
#include "QuICC/ModelOperator/Time.hpp"
#include "QuICC/ModelOperator/ImplicitLinear.hpp"
#include "QuICC/ModelOperator/ExplicitLinear.hpp"
#include "QuICC/ModelOperator/ExplicitNonlinear.hpp"
#include "QuICC/ModelOperator/ExplicitNextstep.hpp"
#include "QuICC/ModelOperator/Stencil.hpp"
#include "QuICC/ModelOperator/Boundary.hpp"
#include "QuICC/ModelOperatorBoundary/FieldToRhs.hpp"
#include "QuICC/ModelOperatorBoundary/SolverHasBc.hpp"
#include "QuICC/ModelOperatorBoundary/SolverNoTau.hpp"
#include "QuICC/ModelOperatorBoundary/Stencil.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Bc/Name/FixedTemperature.hpp"
#include "QuICC/Bc/Name/FixedFlux.hpp"
#include "QuICC/Bc/Name/StressFree.hpp"
#include "QuICC/Bc/Name/NoSlip.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/NonDimensional/Prandtl.hpp"
#include "QuICC/NonDimensional/Rayleigh.hpp"
#include "QuICC/NonDimensional/Ekman.hpp"
#include "QuICC/NonDimensional/Heating.hpp"
#include "QuICC/NonDimensional/RRatio.hpp"
#include "QuICC/NonDimensional/Lower1d.hpp"
#include "QuICC/NonDimensional/Upper1d.hpp"
#include "QuICC/NonDimensional/CflInertial.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/ICondition.hpp"
#include "QuICC/Tools/IdToHuman.hpp"
#include "QuICC/Resolutions/Tools/IndexCounter.hpp"
#include "QuICC/Equations/CouplingIndexType.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Id.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2Y2.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2Y3.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2Y2SphLapl.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2Y3SphLapl.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I4Y4.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I4Y4SphLapl.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I4Y4SphLapl2.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/Operator.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/Value.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/D1.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/D2.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/R1D1DivR1.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/InsulatingShell.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Stencil/Value.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Stencil/D1.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Stencil/R1D1DivR1.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Stencil/ValueD1.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Stencil/ValueD2.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Shell {

namespace RTC {

   IRTCBackend::IRTCBackend()
      : IModelBackend(), mUseGalerkin(false)
   {
   }

   std::vector<std::string> IRTCBackend::fieldNames() const
   {
      std::vector<std::string> names = {
         PhysicalNames::Velocity().tag(),
         PhysicalNames::Temperature().tag()
      };

      return names;
   }

   std::vector<std::string> IRTCBackend::paramNames() const
   {
      std::vector<std::string> names = {
         NonDimensional::Prandtl().tag(),
         NonDimensional::Rayleigh().tag(),
         NonDimensional::Ekman().tag(),
         NonDimensional::Heating().tag(),
         NonDimensional::RRatio().tag()
      };

      return names;
   }

   std::vector<bool> IRTCBackend::isPeriodicBox() const
   {
      std::vector<bool> periodic = {false, false, false};

      return periodic;
   }

   void IRTCBackend::enableGalerkin(const bool flag)
   {
      this->mUseGalerkin = flag;
   }

   bool IRTCBackend::useGalerkin() const
   {
      return this->mUseGalerkin;
   }

   std::map<std::string,MHDFloat> IRTCBackend::automaticParameters(const std::map<std::string,MHDFloat>& cfg) const
   {
      auto E = cfg.find(NonDimensional::Ekman().tag())->second;
      auto rratio = cfg.find(NonDimensional::RRatio().tag())->second;

      std::map<std::string,MHDFloat> params = {
         {NonDimensional::CflInertial().tag(), 0.1*E}
      };

      bool useGapWidth = true;
      if(useGapWidth)
      {
         params.emplace(NonDimensional::Lower1d().tag(), rratio/(1.0 - rratio));
         params.emplace(NonDimensional::Upper1d().tag(), 1.0/(1.0 - rratio));
      }
      else
      {
         params.emplace(NonDimensional::Lower1d().tag(), rratio);
         params.emplace(NonDimensional::Upper1d().tag(), 1.0);
      }

      return params;
   }

   MHDFloat IRTCBackend::effectiveRa(const NonDimensional::NdMap& nds) const
   {
      auto effRa = nds.find(NonDimensional::Rayleigh::id())->second->value();
      auto ro = nds.find(NonDimensional::Upper1d::id())->second->value();

      // Scaled on gap width
      if(ro != 1.0)
      {
         auto T = 1.0/nds.find(NonDimensional::Ekman::id())->second->value();
         effRa *= T/ro;
      }

      return effRa;
   }

   MHDFloat IRTCBackend::effectiveBg(const NonDimensional::NdMap& nds) const
   {
      MHDFloat effBg = 1.0;
      auto ro = nds.find(NonDimensional::Upper1d::id())->second->value();
      auto rratio = nds.find(NonDimensional::RRatio::id())->second->value();
      auto heatingMode = nds.find(NonDimensional::Heating::id())->second->value();

      if(ro == 1.0)
      {
         // Nothing
         //
      }
      // gap width and internal heating
      else if(heatingMode == 0)
      {
         effBg = 2.0/(ro*(1.0 + rratio));
      }
      // gap width and differential heating
      else if(heatingMode == 1)
      {
         effBg = ro*ro*rratio;
      }

      return effBg;
   }


} // RTC
} // Shell
} // Boussinesq
} // Model
} // QuICC
