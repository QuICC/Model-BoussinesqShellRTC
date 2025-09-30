/**
 * @file PhysicalModel.hpp
 * @brief Implementation of the Boussinesq rotating thermal convection in a
 * spherical shell (Toroidal/Poloidal formulation) without coupled solve
 * (standard implementation)
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SHELL_RTC_EXPLICIT_PHYSICALMODEL_HPP
#define QUICC_MODEL_BOUSSINESQ_SHELL_RTC_EXPLICIT_PHYSICALMODEL_HPP

// System includes
//
#include <string>

// Project includes
//
#include "QuICC/SpatialScheme/3D/SLFl.hpp"
#include "QuICC/Model/PyModelBackend.hpp"
#include "Model/Boussinesq/Shell/RTC/Explicit/ModelBackend.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Shell {

namespace RTC {

namespace Explicit {

/**
 * @brief Implementation of the Boussinesq rotating thermal convection spherical
 * shell model (Toroidal/Poloidal formulation) without coupled solve (standard
 * implementation)
 */
template <typename TBuilder> class PhysicalModel : public TBuilder
{
public:
   /// Typedef for the spatial scheme used
   typedef SpatialScheme::SLFl SchemeType;

   /**
    * @brief Constructor
    */
   PhysicalModel() = default;

   /**
    * @brief Destructor
    */
   virtual ~PhysicalModel() = default;

   /**
    * @brief Initialize specialized backend
    */
   void init() final;

protected:
private:
};

template <typename TBuilder> void PhysicalModel<TBuilder>::init()
{
   TBuilder::init();
#ifdef QUICC_MODEL_BOUSSINESQSHELLRTC_EXPLICIT_BACKEND_CPP

   this->mpBackend = std::make_shared<ModelBackend>();
#else
   std::string pyModule = "boussinesq.shell.rtc.explicit.physical_model";
   std::string pyClass = "PhysicalModel";

   this->mpBackend =
      std::make_shared<PyModelBackend>(pyModule, pyClass);
#endif
}

} // namespace Explicit
} // namespace RTC
} // namespace Shell
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC

#endif // QUICC_MODEL_BOUSSINESQ_SHELL_RTC_EXPLICIT_PHYSICALMODEL_HPP
