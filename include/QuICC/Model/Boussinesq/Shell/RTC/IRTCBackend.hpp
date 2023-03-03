/**
 * @file IRTCBackend.hpp
 * @brief Base model backend for RTC model
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SHELL_RTC_IRTCBACKEND_HPP
#define QUICC_MODEL_BOUSSINESQ_SHELL_RTC_IRTCBACKEND_HPP

// System includes
//
#include <string>
#include <vector>
#include <map>
#include <memory>

// Project includes
//
#include "QuICC/Model/IModelBackend.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Shell {

namespace RTC {

   /**
    * @brief Base model backed for RTC model
    */
   class IRTCBackend: public IModelBackend
   {
      public:
         /**
          * @brief Constructor
          */
         IRTCBackend();

         /**
          * @brief Destructor
          */
         virtual ~IRTCBackend() = default;

         /**
          * @brief Get vector of names for the physical fields
          */
         virtual std::vector<std::string> fieldNames() const override;

         /**
          * @brief Get vector of names for the nondimensional parameters
          */
         virtual std::vector<std::string> paramNames() const override;

         /**
          * @brief Get vector of bools about periodic box
          */
         virtual std::vector<bool> isPeriodicBox() const override;

         /**
          * @brief Enable galerkin basis
          */
         virtual void enableGalerkin(const bool flag) override;

         /**
          * @brief Get auotmatically computed parameters based on input parameters
          *
          * @param cfg  Input parameters
          */
         virtual std::map<std::string,MHDFloat> automaticParameters(const std::map<std::string,MHDFloat>& cfg) const override;

      protected:
         /**
          * @brief Use Galerkin basis?
          */
         bool useGalerkin() const;

         /**
          * @brief Compute effective Rayleigh number
          */
         MHDFloat effectiveRa(const NonDimensional::NdMap& nds) const;

         /**
          * @brief Compute effective thermal backgroundn
          */
         MHDFloat effectiveBg(const NonDimensional::NdMap& nds) const;

      private:
         /**
          * @brief Use Galerkin basis?
          */
         bool mUseGalerkin;
   };

} // RTC
} // Shell
} // Boussinesq
} // Model
} // QuICC

#endif // QUICC_MODEL_BOUSSINESQ_SHELL_RTC_IRTCBACKEND_HPP
