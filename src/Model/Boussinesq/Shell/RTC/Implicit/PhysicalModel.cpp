/**
 * @file PhysicalModel.cpp
 * @brief Source of the Boussinesq rotating thermal convection in a spherical shell (Toroidal/Poloidal formulation)
 */

// System includes
//

// Project include
//
#include "QuICC/Model/Boussinesq/Shell/RTC/Implicit/PhysicalModel.hpp"
#include "QuICC/Model/Boussinesq/Shell/RTC/Implicit/ModelBackend.hpp"
#include "QuICC/Model/PyModelBackend.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Shell {

namespace RTC {

namespace Implicit {

   std::string PhysicalModel::PYMODULE()
   {
      return "boussinesq.shell.rtc.implicit.physical_model";
   }

   void PhysicalModel::init()
   {
#ifdef QUICC_MODEL_BOUSSINESQSHELLRTC_IMPLICIT_BACKEND_CPP
      IPhysicalModel<Simulation,StateGenerator,VisualizationGenerator>::init();

      this->mpBackend = std::make_shared<ModelBackend>();
#else
      IPhysicalPyModel<Simulation,StateGenerator,VisualizationGenerator>::init();

      this->mpBackend = std::make_shared<PyModelBackend>(this->PYMODULE(), this->PYCLASS());
#endif
   }

} // Implicit
} // RTC
} // Shell
} // Boussinesq
} // Model
} // QuICC
