/**
 * @file ModelBackend.cpp
 * @brief Source of the interface for model backend
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "QuICC/Model/Boussinesq/Shell/RTC/Explicit/ModelBackend.hpp"
#include "QuICC/ModelOperator/Time.hpp"
#include "QuICC/ModelOperator/ImplicitLinear.hpp"
#include "QuICC/ModelOperator/ExplicitLinear.hpp"
#include "QuICC/ModelOperator/ExplicitNonlinear.hpp"
#include "QuICC/ModelOperator/ExplicitNextstep.hpp"
#include "QuICC/ModelOperator/Stencil.hpp"
#include "QuICC/ModelOperator/Boundary.hpp"
#include "QuICC/ModelOperator/SplitImplicitLinear.hpp"
#include "QuICC/ModelOperator/SplitBoundary.hpp"
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

   ModelBackend::ModelBackend()
      : IRTCBackend(), mcTruncateQI(false)
   {
      this->enableSplitEquation(false);
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
      info.isComplex = false;

      // Use split operators
      info.isSplitEquation = this->useSplitEquation();

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

   void ModelBackend::blockSize(int& tN, int& gN, ArrayI& shift, int& rhs, const SpectralFieldId& fId, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs) const
   {
      assert(eigs.size() == 1);
      int l = eigs.at(0);
      this->blockInfo(tN, gN, shift, rhs, fId, res, l, bcs);
   }

   void ModelBackend::operatorInfo(OperatorInfo& info, const SpectralFieldId& fId, const Resolution& res, const Equations::Tools::ICoupling& coupling, const BcMap& bcs) const
   {
      // Loop overall matrices/eigs
      for(int idx = 0; idx < info.tauN.size(); ++idx)
      {
         auto eigs = coupling.getIndexes(res, idx);

         int tN, gN, rhs;
         ArrayI shift(3);

         this->blockSize(tN, gN, shift, rhs, fId, res, eigs, bcs);

         info.tauN(idx) = tN;
         info.galN(idx) = gN;
         info.galShift.row(idx) = shift;
         info.rhsCols(idx) = rhs;

         // Compute system size
         int sN = 0;
         for(auto f: this->implicitFields(fId))
         {
            this->blockSize(tN, gN, shift, rhs, f, res, eigs, bcs);
            sN += gN;
         }

         if(sN == 0)
         {
            sN = info.galN(idx);
         }

         info.sysN(idx) = sN;
      }
   }

   void ModelBackend::implicitBlock(DecoupledZSparse& decMat, const SpectralFieldId& rowId, const SpectralFieldId& colId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const NonDimensional::NdMap& nds, const bool isSplitOperator) const
   {
      assert(eigs.size() == 1);
      int l = eigs.at(0);

      auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, eigs.at(0))(0);

      auto ri = nds.find(NonDimensional::Lower1d::id())->second->value();
      auto ro = nds.find(NonDimensional::Upper1d::id())->second->value();

      if(rowId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::TOR) && rowId == colId)
      {
         SparseSM::Chebyshev::LinearMap::I2Y2SphLapl spasm(nN, nN, ri, ro, l);
         decMat.real() = spasm.mat();
      }
      else if(rowId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::POL) && rowId == colId)
      {
         if(this->useSplitEquation())
         {
            if(isSplitOperator)
            {
               SparseSM::Chebyshev::LinearMap::I2Y2SphLapl spasm(nN, nN, ri, ro, l);
               decMat.real() = spasm.mat();
            }
            else
            {
               SparseSM::Chebyshev::LinearMap::I2Y2SphLapl spasm(nN, nN, ri, ro, l);
               decMat.real() = spasm.mat();
            }
         }
         else
         {
            SparseSM::Chebyshev::LinearMap::I4Y4SphLapl2 spasm(nN, nN, ri, ro, l);
            decMat.real() = spasm.mat();
         }
      }
      else if(rowId == std::make_pair(PhysicalNames::Temperature::id(), FieldComponents::Spectral::SCALAR) && rowId == colId)
      {
         auto Pr = nds.find(NonDimensional::Prandtl::id())->second->value();
         auto heatingMode = nds.find(NonDimensional::Heating::id())->second->value();

         if(heatingMode == 0)
         {
            SparseSM::Chebyshev::LinearMap::I2Y2SphLapl spasm(nN, nN, ri, ro, l);
            decMat.real() = (1.0/Pr)*spasm.mat();
         }
         else
         {
            SparseSM::Chebyshev::LinearMap::I2Y3SphLapl spasm(nN, nN, ri, ro, l);
            decMat.real() = (1.0/Pr)*spasm.mat();
         }
      }
      else
      {
         throw std::logic_error("Equations are not setup properly");
      }
   }

   void ModelBackend::timeBlock(DecoupledZSparse& decMat, const SpectralFieldId& fieldId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const NonDimensional::NdMap& nds) const
   {
      assert(eigs.size() == 1);
      int l = eigs.at(0);

      auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, eigs.at(0))(0);

      auto ri = nds.find(NonDimensional::Lower1d::id())->second->value();
      auto ro = nds.find(NonDimensional::Upper1d::id())->second->value();

      if(fieldId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::TOR))
      {
         SparseSM::Chebyshev::LinearMap::I2Y2 spasm(nN, nN, ri, ro);
         decMat.real() = spasm.mat();
      }
      else if(fieldId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::POL))
      {
         if(this->useSplitEquation())
         {
            SparseSM::Chebyshev::LinearMap::I2Y2 spasm(nN, nN, ri, ro);
            decMat.real() = spasm.mat();
         }
         else
         {
            SparseSM::Chebyshev::LinearMap::I4Y4SphLapl spasm(nN, nN, ri, ro, l);
            decMat.real() = spasm.mat();
         }
      }
      else if(fieldId == std::make_pair(PhysicalNames::Temperature::id(), FieldComponents::Spectral::SCALAR))
      {
         auto heatingMode = nds.find(NonDimensional::Heating::id())->second->value();

         if(heatingMode == 0)
         {
            SparseSM::Chebyshev::LinearMap::I2Y2 spasm(nN, nN, ri, ro);
            decMat.real() = spasm.mat();
         }
         else
         {
            SparseSM::Chebyshev::LinearMap::I2Y3 spasm(nN, nN, ri, ro);
            decMat.real() = spasm.mat();
         }
      }
   }

   void ModelBackend::modelMatrix(DecoupledZSparse& rModelMatrix, const std::size_t opId, const Equations::CouplingInformation::FieldId_range imRange, const int matIdx, const std::size_t bcType, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const
   {
      assert(eigs.size() == 1);
      int l = eigs.at(0);

      // Time operator
      if(opId == ModelOperator::Time::id())
      {
         bool needStencil = (this->useGalerkin() && bcType == ModelOperatorBoundary::SolverNoTau::id());
         bool needTau = bcType == ModelOperatorBoundary::SolverHasBc::id();

         for(auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
         {
            this->timeBlock(rModelMatrix, *pRowId, matIdx, res, eigs, nds);

            // Apply boundary condition
            if(needStencil)
            {
               this->applyGalerkinStencil(rModelMatrix.real(), *pRowId, *pRowId, l, res, bcs, nds);
            }
            else if(needTau)
            {
               this->applyTau(rModelMatrix.real(), *pRowId, *pRowId, l, res, bcs, nds, false);
            }
         }
      }
      // Linear operator
      else if(opId == ModelOperator::ImplicitLinear::id() || opId == ModelOperator::SplitImplicitLinear::id())
      {
         bool isSplit = (opId == ModelOperator::SplitImplicitLinear::id());
         bool needStencil = (this->useGalerkin() && bcType == ModelOperatorBoundary::SolverNoTau::id());
         bool needTau = bcType == ModelOperatorBoundary::SolverHasBc::id();

         for(auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
         {
            for(auto pColId = imRange.first; pColId != imRange.second; pColId++)
            {
               this->implicitBlock(rModelMatrix, *pRowId, *pColId, matIdx, res, eigs, nds, isSplit);

               // Apply boundary condition
               if(needStencil)
               {
                  this->applyGalerkinStencil(rModelMatrix.real(), *pRowId, *pColId, l, res, bcs, nds);
               }
               else if(needTau)
               {
                  this->applyTau(rModelMatrix.real(), *pRowId, *pColId, l, res, bcs, nds, isSplit);
               }
            }
         }
      }
      // Boundary operator
      else if(opId == ModelOperator::Boundary::id() || opId == ModelOperator::SplitBoundary::id())
      {
         bool isSplit = (opId == ModelOperator::SplitBoundary::id());
         bool needStencil = this->useGalerkin();
         bool needTau = (bcType == ModelOperatorBoundary::SolverHasBc::id());

         auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, eigs.at(0))(0);

         for(auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
         {
            for(auto pColId = imRange.first; pColId != imRange.second; pColId++)
            {
               rModelMatrix.real().resize(nN, nN);

               // Apply boundary condition
               if(needStencil)
               {
                  this->applyGalerkinStencil(rModelMatrix.real(), *pRowId, *pColId, l, res, bcs, nds);
               }
               else if(needTau)
               {
                  this->applyTau(rModelMatrix.real(), *pRowId, *pColId, l, res, bcs, nds, isSplit);
               }
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

      auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, eigs.at(0))(0);

      auto ri = nds.find(NonDimensional::Lower1d::id())->second->value();
      auto ro = nds.find(NonDimensional::Upper1d::id())->second->value();

      int s = this->nBc(rowId);
      // Explicit linear operator
      if(opId == ModelOperator::ExplicitLinear::id())
      {
         if(rowId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::POL) &&
               colId == std::make_pair(PhysicalNames::Temperature::id(),FieldComponents::Spectral::SCALAR))
         {
            auto Ra = effectiveRa(nds);

            SparseSM::Chebyshev::LinearMap::I4Y4 spasm(nN, nN, ri, ro);
            decMat.real() = Ra*spasm.mat();
         }
         else if(rowId == std::make_pair(PhysicalNames::Temperature::id(), FieldComponents::Spectral::SCALAR) && 
               colId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::POL))
         {
            auto heatingMode = nds.find(NonDimensional::Heating::id())->second->value();
            auto bg = effectiveBg(nds);
            auto dl = static_cast<MHDFloat>(l);
            auto ll1 = dl*(dl + 1.0);

            if(heatingMode == 0)
            {
               SparseSM::Chebyshev::LinearMap::I2Y2 spasm(nN, nN, ri, ro);
               decMat.real() = -bg*ll1*spasm.mat();
            }
            else
            {
               SparseSM::Chebyshev::LinearMap::I2 spasm(nN, nN, ri, ro);
               decMat.real() = -bg*ll1*spasm.mat();
            }
         }
         else
         {
            // Nothing to be done
            throw std::logic_error("There are no explicit linear operators");
         }
      }
      // Explicit nonlinear operator
      else if(opId == ModelOperator::ExplicitNonlinear::id())
      {
         if(rowId == std::make_pair(PhysicalNames::Temperature::id(), FieldComponents::Spectral::SCALAR) && rowId == colId)
         {
            auto heatingMode = nds.find(NonDimensional::Heating::id())->second->value();

            if(heatingMode == 0)
            {
               SparseSM::Chebyshev::LinearMap::I2Y2 spasm(nN, nN, ri, ro);
               decMat.real() = spasm.mat();
            }
            else
            {
               SparseSM::Chebyshev::LinearMap::I2Y3 spasm(nN, nN, ri, ro);
               decMat.real() = spasm.mat();
            }
         }
         else
         {
            throw std::logic_error("There are no explicit nonlinear operators");
         }
      }
      // Explicit nextstep operator
      else if(opId == ModelOperator::ExplicitNextstep::id())
      {
         throw std::logic_error("There are no explicit nextstep operators");
      }

      if(this->useGalerkin())
      {
         SparseSM::Chebyshev::LinearMap::Id qId(nN-s, nN, ri, ro, 0, s);
         decMat.real() = qId.mat()*decMat.real();
      }
   }

} // Explicit
} // RTC
} // Shell
} // Boussinesq
} // Model
} // QuICC
