/**
 * @file ModelBackend.cpp
 * @brief Source of the interface for model backend
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "Model/Boussinesq/Shell/RTC/Explicit/ModelBackend.hpp"
#include "QuICC/ModelOperator/Time.hpp"
#include "QuICC/ModelOperator/ImplicitLinear.hpp"
#include "QuICC/ModelOperator/ExplicitLinear.hpp"
#include "QuICC/ModelOperator/ExplicitNonlinear.hpp"
#include "QuICC/ModelOperator/ExplicitNextstep.hpp"
#include "QuICC/ModelOperator/Stencil.hpp"
#include "QuICC/ModelOperator/Boundary.hpp"
#include "QuICC/ModelOperator/SplitImplicitLinear.hpp"
#include "QuICC/ModelOperator/SplitBoundary.hpp"
#include "QuICC/ModelOperator/SplitBoundaryValue.hpp"
#include "QuICC/ModelOperatorBoundary/FieldToRhs.hpp"
#include "QuICC/ModelOperatorBoundary/SolverHasBc.hpp"
#include "QuICC/ModelOperatorBoundary/SolverNoTau.hpp"
#include "QuICC/ModelOperatorBoundary/Stencil.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Bc/Name/FixedTemperature.hpp"
#include "QuICC/Bc/Name/FixedFlux.hpp"
#include "QuICC/Bc/Name/StressFree.hpp"
#include "QuICC/Bc/Name/NoSlip.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/NonDimensional/Prandtl.hpp"
#include "QuICC/NonDimensional/Rayleigh.hpp"
#include "QuICC/NonDimensional/Ekman.hpp"
#include "QuICC/NonDimensional/Heating.hpp"
#include "QuICC/NonDimensional/RRatio.hpp"
#include "QuICC/NonDimensional/Lower1d.hpp"
#include "QuICC/NonDimensional/Upper1d.hpp"
#include "QuICC/NonDimensional/CflInertial.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/ICondition.hpp"
#include "QuICC/Tools/IdToHuman.hpp"
#include "QuICC/Resolutions/Tools/IndexCounter.hpp"
#include "QuICC/Equations/CouplingIndexType.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Id.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2Y2.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2Y3.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2Y2SphLapl.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2Y3SphLapl.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I4Y4.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I4Y4SphLapl.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I4Y4SphLapl2.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/Operator.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/Value.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/D1.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/D2.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/R1D1DivR1.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/InsulatingShell.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Stencil/Value.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Stencil/D1.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Stencil/R1D1DivR1.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Stencil/ValueD1.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Stencil/ValueD2.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Shell {

namespace RTC {

namespace Explicit {

namespace implDetails {

/**
 * @brief Specific options for current model
 */
struct BlockOptionsImpl : public details::BlockOptions
{
    /**
     * @brief default ctor
     */
    BlockOptionsImpl() = default;

    /**
     * @brief default dtor
     */
    virtual ~BlockOptionsImpl() = default;

    /// Harmonic degree l
    int l;
    /// Use truncated quasi-inverse?
    bool truncateQI;
    /// Boundary condition
    std::size_t bcId;
    /// Split operator for influence matrix?
    bool isSplitOperator;
    /// Use split equation for influence matrix?
    bool useSplitEquation;
};
} // namespace implDetails

   ModelBackend::ModelBackend()
      : IRTCBackend(),
#ifdef QUICC_TRANSFORM_CHEBYSHEV_TRUNCATE_QI
      mcTruncateQI(true)
#else
      mcTruncateQI(false)
#endif // QUICC_TRANSFORM_CHEBYSHEV_TRUNCATE_QI
   {
   }

   bool ModelBackend::isComplex(const SpectralFieldId& fId) const
   {
       return false;
   }

   ModelBackend::SpectralFieldIds ModelBackend::implicitFields(const SpectralFieldId& fId) const
   {
      SpectralFieldIds fields = {fId};

      return fields;
   }

   ModelBackend::SpectralFieldIds ModelBackend::explicitLinearFields(const SpectralFieldId& fId) const
   {
      SpectralFieldIds fields;
      if(fId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::POL))
      {
         fields.push_back(std::make_pair(PhysicalNames::Temperature::id(),FieldComponents::Spectral::SCALAR));
      }
      else if(fId == std::make_pair(PhysicalNames::Temperature::id(),FieldComponents::Spectral::SCALAR))
      {
         fields.push_back(std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::POL));
      }

      return fields;
   }

   ModelBackend::SpectralFieldIds ModelBackend::explicitNonlinearFields(const SpectralFieldId& fId) const
   {
      SpectralFieldIds fields;
      if(fId == std::make_pair(PhysicalNames::Temperature::id(),FieldComponents::Spectral::SCALAR))
      {
         fields.push_back(std::make_pair(PhysicalNames::Temperature::id(),FieldComponents::Spectral::SCALAR));
      }

      return fields;
   }

   void ModelBackend::equationInfo(EquationInfo& info, const SpectralFieldId& fId, const Resolution& res) const
   {
      // Operators are real
      info.isComplex = this->isComplex(fId);

      // Splitting 4th poloidal equation into two systems
      if(fId == std::make_pair(PhysicalNames::Velocity::id(), FieldComponents::Spectral::POL))
      {
         info.isSplitEquation = this->useSplitEquation();
      }
      else
      {
         info.isSplitEquation = false;
      }

      // Implicit coupled fields
      info.im = this->implicitFields(fId);

      // Explicit linear terms
      info.exL = this->explicitLinearFields(fId);

      // Explicit nonlinear terms
      info.exNL = this->explicitNonlinearFields(fId);

      // Explicit nextstep terms
      info.exNS.clear();

      // Index mode
      info.indexMode = static_cast<int>(Equations::CouplingIndexType::SLOWEST_MULTI_RHS);
   }

   void ModelBackend::operatorInfo(OperatorInfo& info, const SpectralFieldId& fId, const Resolution& res, const Equations::Tools::ICoupling& coupling, const BcMap& bcs) const
   {
      // Loop overall matrices/eigs
      for(int idx = 0; idx < info.tauN.size(); ++idx)
      {
         auto eigs = coupling.getIndexes(res, idx);

         int tN, gN, rhs;
         ArrayI shift(3);

         this->blockInfo(tN, gN, shift, rhs, fId, res, eigs.at(0), bcs);

         info.tauN(idx) = tN;
         info.galN(idx) = gN;
         info.galShift.row(idx) = shift;
         info.rhsCols(idx) = rhs;

         // Compute system size
         int sN = 0;
         for(auto f: this->implicitFields(fId))
         {
            this->blockInfo(tN, gN, shift, rhs, f, res, eigs.at(0), bcs);
            sN += gN;
         }

         if(sN == 0)
         {
            sN = info.galN(idx);
         }

         info.sysN(idx) = sN;
      }
   }

   std::vector<details::BlockDescription> ModelBackend::implicitBlockBuilder(
         const SpectralFieldId& rowId, const SpectralFieldId& colId,
         const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs,
         const NonDimensional::NdMap& nds, const bool isSplitOperator) const
   {
      std::vector<details::BlockDescription> descr;

      // Create description with common options
      auto getDescription = [&]() -> details::BlockDescription&
      {
         descr.push_back({});
         auto& d = descr.back();
         auto opts = std::make_shared<implDetails::BlockOptionsImpl>();
         opts->l = eigs.at(0);
         opts->bcId = bcs.find(colId.first)->second;
         opts->truncateQI = this->mcTruncateQI;
         opts->isSplitOperator = isSplitOperator;
         opts->useSplitEquation = this->useSplitEquation();
         d.opts = opts;

         return d;
      };

      if (rowId == std::make_pair(PhysicalNames::Velocity::id(),
               FieldComponents::Spectral::TOR) &&
            rowId == colId)
      {
         // Real part of operator
         auto realOp = [](const int nNr, const int nNc, const int l,
               std::shared_ptr<details::BlockOptions> opts,
               const NonDimensional::NdMap& nds)
         {
            SparseMatrix bMat(nNr, nNc);

            if (l > 0)
            {
               auto& o =
                  *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(
                        opts);
               const auto ri = nds.find(NonDimensional::Lower1d::id())->second->value();
               const auto ro = nds.find(NonDimensional::Upper1d::id())->second->value();
               SparseSM::Chebyshev::LinearMap::I2Y2SphLapl spasm(nNr, nNc, ri, ro, l);
               bMat = spasm.mat();
            }

            return bMat;
         };

         // Create block diagonal operator
         auto& d = getDescription();
         d.nRowShift = 0;
         d.nColShift = 0;
         d.realOp = realOp;
         d.imagOp = nullptr;
      }
      else if (rowId == std::make_pair(PhysicalNames::Velocity::id(),
               FieldComponents::Spectral::POL) &&
            rowId == colId)
      {
         // Real part of block
         auto realOp = [](const int nNr, const int nNc, const int l,
               std::shared_ptr<details::BlockOptions> opts,
               const NonDimensional::NdMap& nds)
         {
            SparseMatrix bMat(nNr, nNc);

            if (l > 0)
            {
               auto& o =
                  *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(
                        opts);
               const auto ri = nds.find(NonDimensional::Lower1d::id())->second->value();
               const auto ro = nds.find(NonDimensional::Upper1d::id())->second->value();
               if (o.useSplitEquation)
               {
                  if (o.isSplitOperator)
                  {
                     SparseSM::Chebyshev::LinearMap::I2Y2SphLapl spasm(nNr, nNc, ri, ro, l);
                     bMat = spasm.mat();
                  }
                  else
                  {
                     SparseSM::Chebyshev::LinearMap::I2Y2SphLapl spasm(nNr, nNc, ri, ro, l);
                     bMat = spasm.mat();
                  }
               }
               else
               {
                  SparseSM::Chebyshev::LinearMap::I4Y4SphLapl2 spasm(nNr, nNc, ri, ro, l);
                  bMat = spasm.mat();
               }
            }

            return bMat;
         };

         // Create diagonal block
         auto& d = getDescription();
         d.nRowShift = 0;
         d.nColShift = 0;
         d.realOp = realOp;
         d.imagOp = nullptr;
      }
      else if (rowId == std::make_pair(PhysicalNames::Temperature::id(),
               FieldComponents::Spectral::SCALAR) &&
            rowId == colId)
      {
         // Creat real part of block
         auto realOp = [](const int nNr, const int nNc, const int l,
               std::shared_ptr<details::BlockOptions> opts,
               const NonDimensional::NdMap& nds)
         {
            SparseMatrix bMat(nNr, nNc);

            auto& o =
               *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);

            const auto ri = nds.find(NonDimensional::Lower1d::id())->second->value();
            const auto ro = nds.find(NonDimensional::Upper1d::id())->second->value();
            const auto heatingMode = nds.find(NonDimensional::Heating::id())->second->value();
            const auto Pr =
               nds.find(NonDimensional::Prandtl::id())->second->value();

            if(heatingMode == 0)
            {
               SparseSM::Chebyshev::LinearMap::I2Y2SphLapl spasm(nNr, nNc, ri, ro, l);
               bMat = (1.0/Pr)*spasm.mat();
            }
            else
            {
               SparseSM::Chebyshev::LinearMap::I2Y3SphLapl spasm(nNr, nNc, ri, ro, l);
               bMat = (1.0/Pr)*spasm.mat();
            }

            return bMat;
         };

         // Create diagonal block
         auto& d = getDescription();
         d.nRowShift = 0;
         d.nColShift = 0;
         d.realOp = realOp;
         d.imagOp = nullptr;
      }
      else
      {
         // throw std::logic_error("Equations are not setup properly");
      }

      return descr;
   }

   std::vector<details::BlockDescription> ModelBackend::timeBlockBuilder(
         const SpectralFieldId& rowId, const SpectralFieldId& colId,
         const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs,
         const NonDimensional::NdMap& nds) const
   {
      assert(rowId == colId);
      auto fieldId = rowId;

      std::vector<details::BlockDescription> descr;

      // Create description with common options
      auto getDescription = [&]() -> details::BlockDescription&
      {
         descr.push_back({});
         auto& d = descr.back();
         auto opts = std::make_shared<implDetails::BlockOptionsImpl>();
         opts->l = eigs.at(0);
         opts->bcId = bcs.find(colId.first)->second;
         opts->truncateQI = this->mcTruncateQI;
         opts->isSplitOperator = false;
         opts->useSplitEquation = this->useSplitEquation();
         d.opts = opts;

         return d;
      };

      if (fieldId == std::make_pair(PhysicalNames::Velocity::id(),
               FieldComponents::Spectral::TOR))
      {
         // Real part of operator
         auto realOp = [](const int nNr, const int nNc, const int l,
               std::shared_ptr<details::BlockOptions> opts,
               const NonDimensional::NdMap& nds)
         {
            assert(nNr == nNc);

            SparseMatrix bMat(nNr, nNc);
            auto& o =
               *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);

            const auto ri = nds.find(NonDimensional::Lower1d::id())->second->value();
            const auto ro = nds.find(NonDimensional::Upper1d::id())->second->value();

            if (l > 0)
            {
               SparseSM::Chebyshev::LinearMap::I2Y2 spasm(nNr, nNc, ri, ro);
               bMat = spasm.mat();
            }
            else
            {
               SparseSM::Chebyshev::LinearMap::Id qid(nNr, nNc, ri, ro);
               bMat = qid.mat();
            }

            return bMat;
         };

         // Create block diagonal operator
         auto& d = getDescription();
         d.nRowShift = 0;
         d.nColShift = 0;
         d.realOp = realOp;
         d.imagOp = nullptr;
      }
      else if (fieldId == std::make_pair(PhysicalNames::Velocity::id(),
               FieldComponents::Spectral::POL))
      {
         // Real part of operator
         auto realOp = [](const int nNr, const int nNc, const int l,
               std::shared_ptr<details::BlockOptions> opts,
               const NonDimensional::NdMap& nds)
         {
            assert(nNr == nNc);

            SparseMatrix bMat(nNr, nNc);
            auto& o =
               *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);

            const auto ri = nds.find(NonDimensional::Lower1d::id())->second->value();
            const auto ro = nds.find(NonDimensional::Upper1d::id())->second->value();

            if (l > 0)
            {
               if (o.useSplitEquation)
               {
                  SparseSM::Chebyshev::LinearMap::I2Y2 spasm(nNr, nNc, ri, ro);
                  bMat = spasm.mat();
               }
               else
               {

                  SparseSM::Chebyshev::LinearMap::I4Y4SphLapl spasm(nNr, nNc, ri, ro, l);

                  // Correction is not yet tested
#if 0
                  // Correct Laplacian for 4th order system according to:
                  // McFadden,Murray,Boisvert,
                  // Elimination of Spurious Eigenvalues in the
                  // Chebyshev Tau Spectral Method,
                  // JCP 91, 228-239 (1990)
                  // We simply drop the last two column
                  if (o.bcId == Bc::Name::NoSlip::id())
                  {
                     SparseSM::Chebyshev::LinearMap::Id qid(nNr, nNc, ri, ro, -2);
                     bMat = spasm.mat() * qid.mat();
                  }
                  else
                  {
                     bMat = spasm.mat();
                  }
#else
                  bMat = spasm.mat();
#endif
               }
            }
            else
            {
               SparseSM::Chebyshev::LinearMap::Id qid(nNr, nNc, ri, ro);
               bMat = qid.mat();
            }

            return bMat;
         };

         // Create block diagonal operator
         auto& d = getDescription();
         d.nRowShift = 0;
         d.nColShift = 0;
         d.realOp = realOp;
         d.imagOp = nullptr;
      }
      else if (fieldId == std::make_pair(PhysicalNames::Temperature::id(),
               FieldComponents::Spectral::SCALAR))
      {
         // Real part of operator
         auto realOp = [](const int nNr, const int nNc, const int l,
               std::shared_ptr<details::BlockOptions> opts,
               const NonDimensional::NdMap& nds)
         {
            auto& o =
               *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);

            const auto ri = nds.find(NonDimensional::Lower1d::id())->second->value();
            const auto ro = nds.find(NonDimensional::Upper1d::id())->second->value();
            const auto heatingMode = nds.find(NonDimensional::Heating::id())->second->value();

            SparseMatrix bMat(nNr, nNc);
            if(heatingMode == 0)
            {
               SparseSM::Chebyshev::LinearMap::I2Y2 spasm(nNr, nNc, ri, ro);
               bMat = spasm.mat();
            }
            else
            {
               SparseSM::Chebyshev::LinearMap::I2Y3 spasm(nNr, nNc, ri, ro);
               bMat = spasm.mat();
            }

            return bMat;
         };

         // Create block diagonal operator
         auto& d = getDescription();
         d.nRowShift = 0;
         d.nColShift = 0;
         d.realOp = realOp;
         d.imagOp = nullptr;
      }

      return descr;
   }

   std::vector<details::BlockDescription> ModelBackend::boundaryBlockBuilder(
         const SpectralFieldId& rowId, const SpectralFieldId& colId,
         const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs,
         const NonDimensional::NdMap& nds, const bool isSplit) const
   {
      std::vector<details::BlockDescription> descr;

      // Create description with common options
      auto getDescription = [&]() -> details::BlockDescription&
      {
         descr.push_back({});
         auto& d = descr.back();
         auto opts = std::make_shared<implDetails::BlockOptionsImpl>();
         opts->l = eigs.at(0);
         opts->bcId = bcs.find(colId.first)->second;
         opts->truncateQI = this->mcTruncateQI;
         opts->isSplitOperator = isSplit;
         opts->useSplitEquation = this->useSplitEquation();
         d.opts = opts;

         return d;
      };

      if (rowId == colId)
      {
         // Real part of operator
         auto realOp = [](const int nNr, const int nNc, const int l,
               std::shared_ptr<details::BlockOptions> opts,
               const NonDimensional::NdMap& nds)
         {
            SparseMatrix bMat(nNr, nNc);

            return bMat;
         };

         // Create block diagonal operator
         auto& d = getDescription();
         d.nRowShift = 0;
         d.nColShift = 0;
         d.realOp = realOp;
         d.imagOp = nullptr;
      }

      return descr;
   }

   std::vector<details::BlockDescription>
   ModelBackend::splitBoundaryValueBlockBuilder(const SpectralFieldId& rowId,
         const SpectralFieldId& colId, const Resolution& res,
         const std::vector<MHDFloat>& eigs, const BcMap& bcs,
         const NonDimensional::NdMap& nds) const
   {
      assert(rowId == colId);
      auto fieldId = rowId;

      std::vector<details::BlockDescription> descr;

      // Create description with common options
      auto getDescription = [&]() -> details::BlockDescription&
      {
         descr.push_back({});
         auto& d = descr.back();
         auto opts = std::make_shared<implDetails::BlockOptionsImpl>();
         opts->l = eigs.at(0);
         opts->bcId = bcs.find(colId.first)->second;
         opts->truncateQI = this->mcTruncateQI;
         opts->isSplitOperator = false;
         opts->useSplitEquation = this->useSplitEquation();
         d.opts = opts;

         return d;
      };

      if (fieldId == std::make_pair(PhysicalNames::Velocity::id(),
               FieldComponents::Spectral::POL))
      {
         // Boundary value operator
         auto bcValOp = [](const int nNr, const int nNc, const int l,
               std::shared_ptr<details::BlockOptions> opts,
               const NonDimensional::NdMap& nds)
         {
            assert(nNr == nNc);

            SparseMatrix bMat(nNr, 2);

            Eigen::Triplet<MHDFloat> valTop = {0, 0, 1.0};
            Eigen::Triplet<MHDFloat> valBottom = {1, 1, 1.0};
            std::vector<Eigen::Triplet<MHDFloat>> triplets = {valTop, valBottom};
            bMat.setFromTriplets(triplets.begin(), triplets.end());

            return bMat;
         };

         // Create block diagonal operator
         auto& d = getDescription();
         d.nRowShift = 0;
         d.nColShift = 0;
         d.realOp = bcValOp;
         d.imagOp = bcValOp;
      }

      return descr;
   }

   void ModelBackend::modelMatrix(DecoupledZSparse& rModelMatrix, const std::size_t opId, const Equations::CouplingInformation::FieldId_range imRange, const int matIdx, const std::size_t bcType, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const
   {
      assert(eigs.size() == 1);
      int l = eigs.at(0);

      // Time operator
      if(opId == ModelOperator::Time::id())
      {
         for (auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
         {
            auto rowId = *pRowId;
            auto colId = rowId;
            auto descr = timeBlockBuilder(rowId, colId, res, eigs, bcs, nds);
            buildBlock(rModelMatrix, descr, rowId, colId, matIdx, bcType, res,
                  l, l, bcs, nds, false);
         }
      }
      // Linear operator
      else if(opId == ModelOperator::ImplicitLinear::id() || opId == ModelOperator::SplitImplicitLinear::id())
      {
         bool isSplit = (opId == ModelOperator::SplitImplicitLinear::id());

         for (auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
         {
            for (auto pColId = imRange.first; pColId != imRange.second;
                  pColId++)
            {
               auto rowId = *pRowId;
               auto colId = *pColId;
               auto descr = implicitBlockBuilder(rowId, colId, res, eigs, bcs,
                     nds, isSplit);
               buildBlock(rModelMatrix, descr, rowId, colId, matIdx, bcType,
                     res, l, l, bcs, nds, isSplit);
            }
         }
      }
      // Boundary operator
      else if(opId == ModelOperator::Boundary::id() || opId == ModelOperator::SplitBoundary::id())
      {
         bool isSplit = (opId == ModelOperator::SplitBoundary::id());

         for (auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
         {
            for (auto pColId = imRange.first; pColId != imRange.second;
                  pColId++)
            {
               auto rowId = *pRowId;
               auto colId = *pColId;
               auto descr = boundaryBlockBuilder(rowId, colId, res, eigs, bcs,
                     nds, isSplit);
               buildBlock(rModelMatrix, descr, rowId, colId, matIdx, bcType,
                     res, l, l, bcs, nds, isSplit);
            }
         }
      }
      // Split equation boundary value
      else if(opId == ModelOperator::SplitBoundaryValue::id())
      {
         for (auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
         {
            for (auto pColId = imRange.first; pColId != imRange.second;
                  pColId++)
            {
               auto rowId = *pRowId;
               auto colId = *pColId;
               auto descr = splitBoundaryValueBlockBuilder(rowId, colId, res,
                     eigs, bcs, nds);
               buildBlock(rModelMatrix, descr, rowId, colId, matIdx, bcType,
                     res, l, l, bcs, nds, false);
            }
         }
      }
      else
      {
         throw std::logic_error("Requested operator type is not implemented");
      }
   }

   void ModelBackend::galerkinStencil(SparseMatrix& mat, const SpectralFieldId& fieldId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const bool makeSquare, const BcMap& bcs, const NonDimensional::NdMap& nds) const
   {
      assert(eigs.size() == 1);
      int l = eigs.at(0);
      this->stencil(mat, fieldId, l, res, makeSquare, bcs, nds);
   }

   void ModelBackend::explicitBlock(DecoupledZSparse& decMat, const SpectralFieldId& rowId, const std::size_t opId,  const SpectralFieldId colId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const
   {
      assert(eigs.size() == 1);
      int l = eigs.at(0);

      auto bcType = ModelOperatorBoundary::SolverNoTau::id();

      // Explicit linear operator
      if(opId == ModelOperator::ExplicitLinear::id())
      {
         auto descr = explicitLinearBlockBuilder(rowId, colId, res, eigs, bcs, nds);
         buildBlock(decMat, descr, rowId, colId, matIdx, bcType, res,
               l, l, bcs, nds, false, true);
      }
      // Explicit nonlinear operator
      else if(opId == ModelOperator::ExplicitNonlinear::id())
      {
         auto descr = explicitNonlinearBlockBuilder(rowId, colId, res, eigs, bcs, nds);
         buildBlock(decMat, descr, rowId, colId, matIdx, bcType, res,
               l, l, bcs, nds, false, true);
      }
      // Explicit nextstep operator
      else if(opId == ModelOperator::ExplicitNextstep::id())
      {
         throw std::logic_error("There are no explicit nextstep operators");
      }
   }

   std::vector<details::BlockDescription> ModelBackend::explicitLinearBlockBuilder(const SpectralFieldId& rowId,  const SpectralFieldId& colId, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const
   {
      assert(eigs.size() == 1);

      std::vector<details::BlockDescription> descr;

      // Create description with common options
      auto getDescription = [&]() -> details::BlockDescription&
      {
         descr.push_back({});
         auto& d = descr.back();
         auto opts = std::make_shared<implDetails::BlockOptionsImpl>();
         opts->l = eigs.at(0);
         opts->bcId = bcs.find(colId.first)->second;
         opts->truncateQI = this->mcTruncateQI;
         opts->isSplitOperator = false;
         opts->useSplitEquation = this->useSplitEquation();
         d.opts = opts;

         return d;
      };

      if(rowId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::POL) &&
            colId == std::make_pair(PhysicalNames::Temperature::id(),FieldComponents::Spectral::SCALAR))
      {
         // Real part of operator
         auto realOp = [](const int nNr, const int nNc, const int l,
               std::shared_ptr<details::BlockOptions> opts,
               const NonDimensional::NdMap& nds)
         {
            SparseMatrix bMat(nNr, nNc);

            auto& o =
               *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);

            const auto ri = nds.find(NonDimensional::Lower1d::id())->second->value();
            const auto ro = nds.find(NonDimensional::Upper1d::id())->second->value();
            auto Ra = nds.find(NonDimensional::Rayleigh::id())->second->value();
            // Scaled on gap width
            if(ro != 1.0)
            {
               const auto T = 1.0/nds.find(NonDimensional::Ekman::id())->second->value();
               Ra *= T/ro;
            }

            if(o.useSplitEquation)
            {
               SparseSM::Chebyshev::LinearMap::I2Y2 spasm(nNr, nNc, ri, ro);
               bMat = Ra*spasm.mat();
            }
            else
            {
               SparseSM::Chebyshev::LinearMap::I4Y4 spasm(nNr, nNc, ri, ro);
               bMat = Ra*spasm.mat();
            }

            return bMat;
         };

         // Create block diagonal operator
         auto& d = getDescription();
         d.nRowShift = 0;
         d.nColShift = 0;
         d.realOp = realOp;
         d.imagOp = nullptr;
      }
      else if(rowId == std::make_pair(PhysicalNames::Temperature::id(), FieldComponents::Spectral::SCALAR) && 
            colId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::POL))
      {
         // Real part of operator
         auto realOp = [](const int nNr, const int nNc, const int l,
               std::shared_ptr<details::BlockOptions> opts,
               const NonDimensional::NdMap& nds)
         {
            SparseMatrix bMat(nNr, nNc);
            auto& o =
               *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);

            const auto ri = nds.find(NonDimensional::Lower1d::id())->second->value();
            const auto ro = nds.find(NonDimensional::Upper1d::id())->second->value();
            const auto heatingMode = nds.find(NonDimensional::Heating::id())->second->value();
            const auto rratio = nds.find(NonDimensional::RRatio::id())->second->value();
            MHDFloat bg = 1.0;
            if(ro == 1.0)
            {
               // Nothing
            }
            // gap width and internal heating
            else if(heatingMode == 0)
            {
               bg = 2.0/(ro*(1.0 + rratio));
            }
            // gap width and differential heating
            else if(heatingMode == 1)
            {
               bg = ro*ro*rratio;
            }
            const auto dl = static_cast<MHDFloat>(l);
            const auto ll1 = dl*(dl + 1.0);

            if(heatingMode == 0)
            {
               SparseSM::Chebyshev::LinearMap::I2Y2 spasm(nNr, nNc, ri, ro);
               bMat = -bg*ll1*spasm.mat();
            }
            else
            {
               SparseSM::Chebyshev::LinearMap::I2 spasm(nNr, nNc, ri, ro);
               bMat = -bg*ll1*spasm.mat();
            }

            return bMat;
         };

         // Create block diagonal operator
         auto& d = getDescription();
         d.nRowShift = 0;
         d.nColShift = 0;
         d.realOp = realOp;
         d.imagOp = nullptr;
      }
      else
      {
         // Nothing to be done
         throw std::logic_error("There are no explicit linear operators");
      }

      return descr;
   }

   std::vector<details::BlockDescription> ModelBackend::explicitNonlinearBlockBuilder(const SpectralFieldId& rowId,  const SpectralFieldId& colId, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const
   {
      assert(eigs.size() == 1);

      std::vector<details::BlockDescription> descr;

      // Create description with common options
      auto getDescription = [&]() -> details::BlockDescription&
      {
         descr.push_back({});
         auto& d = descr.back();
         auto opts = std::make_shared<implDetails::BlockOptionsImpl>();
         opts->l = eigs.at(0);
         opts->bcId = bcs.find(colId.first)->second;
         opts->truncateQI = this->mcTruncateQI;
         opts->isSplitOperator = false;
         opts->useSplitEquation = this->useSplitEquation();
         d.opts = opts;

         return d;
      };

      if(rowId == std::make_pair(PhysicalNames::Temperature::id(), FieldComponents::Spectral::SCALAR) && rowId == colId)
      {
         // Real part of operator
         auto realOp = [](const int nNr, const int nNc, const int l,
               std::shared_ptr<details::BlockOptions> opts,
               const NonDimensional::NdMap& nds)
         {
            SparseMatrix bMat(nNr, nNc);
            auto& o =
               *std::dynamic_pointer_cast<implDetails::BlockOptionsImpl>(opts);

            const auto ri = nds.find(NonDimensional::Lower1d::id())->second->value();
            const auto ro = nds.find(NonDimensional::Upper1d::id())->second->value();
            const auto heatingMode = nds.find(NonDimensional::Heating::id())->second->value();

            if(heatingMode == 0)
            {
               SparseSM::Chebyshev::LinearMap::I2Y2 spasm(nNr, nNc, ri, ro);
               bMat = spasm.mat();
            }
            else
            {
               SparseSM::Chebyshev::LinearMap::I2Y3 spasm(nNr, nNc, ri, ro);
               bMat = spasm.mat();
            }

            return bMat;
         };

         // Create block diagonal operator
         auto& d = getDescription();
         d.nRowShift = 0;
         d.nColShift = 0;
         d.realOp = realOp;
         d.imagOp = nullptr;
      }
      else
      {
         throw std::logic_error("There are no explicit nonlinear operators");
      }

      return descr;
   }

} // namespace Explicit
} // namespace RTC
} // namespace Shell
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC
