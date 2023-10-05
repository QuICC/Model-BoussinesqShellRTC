/**
 * @file Momentum.hpp
 * @brief Implementation of the vector Navier-Stokes equation for the Boussinesq rotating thermal convection spherical shell
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SHELL_RTC_MOMENTUM_HPP
#define QUICC_MODEL_BOUSSINESQ_SHELL_RTC_MOMENTUM_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Equations/IVectorEquation.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Shell {

namespace RTC {

   /**
    * @brief Implementation of the vector Navier-Stokes equation for the Boussinesq rotating thermal convection in a spherical shell
    */
   class Momentum: public IVectorEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         Momentum(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend);

         /**
          * @brief Simple empty destructor
          */
         virtual ~Momentum();

         /**
          * @brief Initialize nonlinear interaction kernel
          */
         virtual void initNLKernel(const bool force = false) override;

      protected:
         /**
          * @brief Set variable requirements
          */
         virtual void setRequirements() override;

         /**
          * @brief Set the equation coupling information
          */
         virtual void setCoupling() override;

         /**
          * @brief Set the nonlinear integration components
          */
         virtual void setNLComponents() override;

      private:
   };

}
}
}
}
}

#endif // QUICC_MODEL_BOUSSINESQ_SHELL_RTC_MOMENTUM_HPP
