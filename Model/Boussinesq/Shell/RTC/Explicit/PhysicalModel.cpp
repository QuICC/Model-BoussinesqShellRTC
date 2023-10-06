/**
 * @file PhysicalModel.cpp
 * @brief Source of the Boussinesq rotating thermal convection in a spherical shell (Toroidal/Poloidal formulation) without coupled solve (standard implementation)
 */

// System includes
//

// Project includes
//
#include "Model/Boussinesq/Shell/RTC/Explicit/PhysicalModel.hpp"
#include "Model/Boussinesq/Shell/RTC/Explicit/ModelBackend.hpp"
#include "QuICC/Model/PyModelBackend.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Shell {

namespace RTC {

namespace Explicit {

   std::string PhysicalModel::PYMODULE()
   {
      return "boussinesq.shell.rtc.explicit.physical_model";
   }

   void PhysicalModel::init()
   {
#ifdef QUICC_MODEL_BOUSSINESQSHELLRTC_EXPLICIT_BACKEND_CPP
      IPhysicalModel<Simulation,StateGenerator,VisualizationGenerator>::init();

      this->mpBackend = std::make_shared<ModelBackend>();
#else
      IPhysicalPyModel<Simulation,StateGenerator,VisualizationGenerator>::init();

      this->mpBackend = std::make_shared<PyModelBackend>(this->PYMODULE(), this->PYCLASS());
#endif
   }

} // Explicit
} // RTC
} // Shell
} // Boussinesq
} // Model
} // QuICC
