/**
 * @file PhysicalModel.hpp
 * @brief Implementation of the Boussinesq rotating thermal convection in a
 * spherical shell (Toroidal/Poloidal formulation)
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SHELL_RTC_IMPLICIT_PHYSICALMODEL_HPP
#define QUICC_MODEL_BOUSSINESQ_SHELL_RTC_IMPLICIT_PHYSICALMODEL_HPP

// System includes
//
#include <string>

// Project includes
//
#include "QuICC/SpatialScheme/3D/SLFm.hpp"
#include "QuICC/Model/PyModelBackend.hpp"
#include "Model/Boussinesq/Shell/RTC/Implicit/ModelBackend.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Shell {

namespace RTC {

namespace Implicit {

/**
 * @brief Implementation of the Boussinesq rotating thermal convection spherical
 * shell model (Toroidal/Poloidal formulation)
 */
template <typename TBuilder> class PhysicalModel : public TBuilder
{
public:
   /// Typedef for the spatial scheme used
   typedef SpatialScheme::SLFm SchemeType;

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
#ifdef QUICC_MODEL_BOUSSINESQSHELLRTC_IMPLICIT_BACKEND_CPP
   this->mpBackend = std::make_shared<ModelBackend>();
#else
   std::string pyModule = "boussinesq.shell.rtc.implicit.physical_model";
   std::string pyClass = "PhysicalModel";

   this->mpBackend =
      std::make_shared<PyModelBackend>(pyModule, pyClass);
#endif
}

} // namespace Implicit
} // namespace RTC
} // namespace Shell
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC

#endif // QUICC_MODEL_BOUSSINESQ_SHELL_RTC_IMPLICIT_PHYSICALMODEL_HPP
