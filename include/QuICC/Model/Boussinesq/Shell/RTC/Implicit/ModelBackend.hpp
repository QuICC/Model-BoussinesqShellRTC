/**
 * @file ModelBackend.hpp
 * @brief Model backend for the implicit treatment for Coriolis term
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SHELL_RTC_IMPLICIT_MODELBACKEND_HPP
#define QUICC_MODEL_BOUSSINESQ_SHELL_RTC_IMPLICIT_MODELBACKEND_HPP

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

namespace Implicit {

   namespace internal {
      struct SystemInfo
      {
         int systemSize;
         int blockRows;
         int blockCols;
         int startRow;
         int startCol;

         SystemInfo(const int size, const int rows, const int cols, const int row, const int col)
            : systemSize(size), blockRows(rows), blockCols(cols), startRow(row), startCol(col)
         {
         };
      };
   }

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
         virtual void equationInfo(bool& isComplex, SpectralFieldIds& im, SpectralFieldIds& exL, SpectralFieldIds& exNL, SpectralFieldIds& exNS, int& indexMode, const SpectralFieldId& fId, const Resolution& res) const override;

         /**
          * @brief Get operator information
          */
         virtual void operatorInfo(ArrayI& tauN, ArrayI& galN, MatrixI& galShift, ArrayI& rhsCols, ArrayI& sysN, const SpectralFieldId& fId, const Resolution& res, const Equations::Tools::ICoupling& coupling, const BcMap& bcs) const override;

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
          * @brief Get operator block information
          */
         void blockInfo(int& tN, int& gN, ArrayI& shift, int& rhs, const SpectralFieldId& fId, const Resolution& res, const MHDFloat l, const BcMap& bcs) const;

         /**
          * @brief Build implicit matrix block
          */
         void implicitBlock(DecoupledZSparse& decMat, const SpectralFieldId& rowId, const SpectralFieldId& colId, const int matIdx, const std::size_t bcType, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const;

         /**
          * @brief Build time matrix block
          */
         void timeBlock(DecoupledZSparse& decMat, const SpectralFieldId& fieldId, const int matIdx, const std::size_t bcType, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const;

         /**
          * @brief Build boundary matrix block
          */
         void boundaryBlock(DecoupledZSparse& decMat, const SpectralFieldId& rowId, const SpectralFieldId& colId, const int matIdx, const std::size_t bcType, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const;

         /**
          * @brief Apply boundary condition
          */
         void applyBoundary(DecoupledZSparse& decMat, const SpectralFieldId& rowId, const SpectralFieldId& colId, const int matIdx, const std::size_t bcType, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const;

         /**
          * @brief Apply galerkin stencil for boundary condition
          */
         void applyGalerkinStencil(SparseMatrix& decMat, const SpectralFieldId& rowId, const SpectralFieldId& colId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const;

         /**
          * @brief Apply tau line for boundary condition
          */
         void applyTau(SparseMatrix& mat, const SpectralFieldId& rowId, const SpectralFieldId& colId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const;

         /**
          * @brief Boundary condition stencil
          */
         virtual void stencil(SparseMatrix& mat, const SpectralFieldId& fId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const bool makeSquare, const BcMap& bcs, const NonDimensional::NdMap& nds) const;

      private:
         /**
          * @brief Add block matrix to full system matrix
          */
         void addBlock(SparseMatrix& mat, const SparseMatrix& block, const int rowShift, const int colShift, const MHDFloat coeff = 1.0) const;

         /**
          * @brief Get operator block size
          */
         int blockSize(const SpectralFieldId& fId, const int m, const Resolution& res, const BcMap& bcs, const bool isGalerkin) const;

         /**
          * @brief Get operator block shape
          */
         std::pair<int,int> blockShape(const SpectralFieldId& rowId, const SpectralFieldId& colId, const int m, const Resolution& res, const BcMap& bcs, const bool isGalerkin, const bool dropRows) const;

         /**
          * @brief Compute size information of full system
          *
          * @param dropRows Number of rows to drop
          */
         internal::SystemInfo systemInfo(const SpectralFieldId& colId, const SpectralFieldId& rowId, const int m, const Resolution& res, const BcMap& bcs, const bool isGalerkin, const bool dropRows) const;
   };

} // Implicit
} // RTC
} // Shell
} // Boussinesq
} // Model
} // QuICC

#endif // QUICC_MODEL_BOUSSINESQ_SHELL_RTC_IMPLICIT_MODELBACKEND_HPP
