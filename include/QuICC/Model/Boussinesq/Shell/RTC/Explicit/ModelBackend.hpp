/**
 * @file ModelBackend.hpp
 * @brief Model backend with explicit treatment of the Coriolis term
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SHELL_RTC_EXPLICIT_MODELBACKEND_HPP
#define QUICC_MODEL_BOUSSINESQ_SHELL_RTC_EXPLICIT_MODELBACKEND_HPP

// System includes
//
#include <string>
#include <vector>
#include <map>
#include <memory>

// Project includes
//
#include "QuICC/Model/Boussinesq/Shell/RTC/IRTCBackend.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Shell {

namespace RTC {

namespace Explicit {

   /**
    * @brief Interface for model backend
    */
   class ModelBackend: public IRTCBackend
   {
      public:
         /**
          * @brief Constructor
          */
         ModelBackend();

         /**
          * @brief Destructor
          */
         virtual ~ModelBackend() = default;

         /**
          * @brief Get equation information
          */
         virtual void equationInfo(EquationInfo& info, const SpectralFieldId& fId, const Resolution& res) const override;

         /**
          * @brief Get operator information
          */
         virtual void operatorInfo(OperatorInfo& info, const SpectralFieldId& fId, const Resolution& res, const Equations::Tools::ICoupling& coupling, const BcMap& bcs) const override;

         /**
          * @brief Build model matrix
          */
         virtual void modelMatrix(DecoupledZSparse& rModelMatrix, const std::size_t opId, const Equations::CouplingInformation::FieldId_range imRange, const int matIdx, const std::size_t bcType, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const override;

         /**
          * @brief Build galerkin stencil
          */
         virtual void galerkinStencil(SparseMatrix& mat, const SpectralFieldId& fId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const bool makeSquare, const BcMap& bcs, const NonDimensional::NdMap& nds) const override;

         /**
          * @brief Build explicit block
          */
         virtual void explicitBlock(DecoupledZSparse& mat, const SpectralFieldId& fId, const std::size_t opId,  const SpectralFieldId fieldId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const override;

      protected:
         /**
          * @brief Set field coupling in implicit model matrix
          */
         SpectralFieldIds implicitFields(const SpectralFieldId& fId) const;

         /**
          * @brief Set field coupling in explicit linear terms
          */
         SpectralFieldIds explicitLinearFields(const SpectralFieldId& fId) const;

         /**
          * @brief Set field coupling in explicit nonlinear terms
          */
         SpectralFieldIds explicitNonlinearFields(const SpectralFieldId& fId) const;

         /**
          * @brief Get operator information
          */
         void blockSize(int& tN, int& gN, ArrayI& shift, int& rhs, const SpectralFieldId& fId, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs) const;

         /**
          * @brief Build implicit matrix block
          */
         void implicitBlock(DecoupledZSparse& decMat, const SpectralFieldId& rowId, const SpectralFieldId& colId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const NonDimensional::NdMap& nds, const bool isSplitOperator) const;

         /**
          * @brief Build time matrix block
          */
         void timeBlock(DecoupledZSparse& decMat, const SpectralFieldId& fieldId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const NonDimensional::NdMap& nds) const;

         /**
          * @brief Build inhomogeneous boundary value for split equation
          */
         void splitBoundaryValueBlock(DecoupledZSparse& decMat, const SpectralFieldId& fieldId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const NonDimensional::NdMap& nds) const;

      private:
         /**
          * @brief Truncate quasi-inverse operators?
          */
         const bool mcTruncateQI;

   };

} // Explicit
} // RTC
} // Shell
} // Boussinesq
} // Model
} // QuICC

#endif // QUICC_MODEL_BOUSSINESQ_SHELL_RTC_EXPLICIT_MODELBACKEND_HPP
