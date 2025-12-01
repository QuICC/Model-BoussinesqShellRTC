/**
 * @file IRTCState.hpp
 * @brief Implementation of the Boussinesq rotating thermal convection in a
 * spherical shell (Toroidal/Poloidal formulation)
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SHELL_RTC_IRTCSTATE_HPP
#define QUICC_MODEL_BOUSSINESQ_SHELL_RTC_IRTCSTATE_HPP

// System includes
//
#include <string>

// Project includes
//
#include "QuICC/Generator/StateGenerator.hpp"
#include "QuICC/Model/IStateGeneratorBuilder.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Shell {

namespace RTC {

/**
 * @brief Implementation of the Boussinesq rotating thermal convection spherical
 * shell model (Toroidal/Poloidal formulation)
 */
class IRTCState : public IStateGeneratorBuilder<StateGenerator>
{
public:
   /**
    * @brief Constructor
    */
   IRTCState() = default;

   /**
    * @brief Destructor
    */
   virtual ~IRTCState() = default;

   /// Formulation used for vector fields
   virtual VectorFormulation::Id SchemeFormulation() override;

   /**
    * @brief Version string
    */
   std::string version() const final;

   /**
    * @brief Add the initial state generation equations
    *
    * @param spGen   Shared generator object
    */
   virtual void addStates(SharedStateGenerator spGen) override;

protected:
private:
};

} // namespace RTC
} // namespace Shell
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC

#endif // QUICC_MODEL_BOUSSINESQ_SHELL_RTC_IRTCSTATE_HPP
