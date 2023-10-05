/**
 * @file PhysicalModel.hpp
 * @brief Implementation of the Boussinesq rotating thermal convection in a spherical shell (Toroidal/Poloidal formulation) without coupled solve (standard implementation)
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SHELL_RTC_EXPLICIT_PHYSICALMODEL_HPP
#define QUICC_MODEL_BOUSSINESQ_SHELL_RTC_EXPLICIT_PHYSICALMODEL_HPP

// System includes
//
#include <string>

// Project includes
//
#include "Model/Boussinesq//Shell/RTC/IRTCModel.hpp"
#include "QuICC/SpatialScheme/3D/SLFl.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Shell {

namespace RTC {

namespace Explicit {

   /**
    * @brief Implementation of the Boussinesq rotating thermal convection spherical shell model (Toroidal/Poloidal formulation) without coupled solve (standard implementation)
    */
   class PhysicalModel: public IRTCModel
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

         /// Python script/module name
         virtual std::string PYMODULE() final;

         /**
          * @brief Initialize specialized backend
          */
         void init() final;

      protected:

      private:
   };

} // Explicit
} // RTC
} // Shell
} // Boussinesq
} // Model
} // QuICC

#endif // QUICC_MODEL_BOUSSINESQ_SHELL_RTC_EXPLICIT_PHYSICALMODEL_HPP
