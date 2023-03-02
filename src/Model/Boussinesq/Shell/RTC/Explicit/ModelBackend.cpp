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

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Shell {

namespace RTC {

namespace Explicit {

   ModelBackend::ModelBackend()
      : IModelBackend(), mUseGalerkin(false)
   {
   }

   std::vector<std::string> ModelBackend::fieldNames() const
   {
      std::vector<std::string> names = {
         PhysicalNames::Velocity().tag(),
         PhysicalNames::Temperature().tag()
      };

      return names;
   }

   std::vector<std::string> ModelBackend::paramNames() const
   {
      std::vector<std::string> names = {
         NonDimensional::Prandtl().tag(),
         NonDimensional::Rayleigh().tag(),
         NonDimensional::Ekman().tag(),
         NonDimensional::Heating().tag(),
         NonDimensional::RRatio().tag()
      };

      return names;
   }

   std::vector<bool> ModelBackend::isPeriodicBox() const
   {
      std::vector<bool> periodic = {false, false, false};

      return periodic;
   }

   void ModelBackend::enableGalerkin(const bool flag)
   {
      this->mUseGalerkin = flag;
   }

   std::map<std::string,MHDFloat> ModelBackend::automaticParameters(const std::map<std::string,MHDFloat>& cfg) const
   {
      auto E = cfg.find(NonDimensional::Ekman().tag())->second;
      auto rratio = cfg.find(NonDimensional::RRatio().tag())->second;

      std::map<std::string,MHDFloat> params = {
         {NonDimensional::CflInertial().tag(), 0.1*E}
      };

      bool useGapWidth = true;
      if(useGapWidth)
      {
         params.emplace(NonDimensional::Lower1d().tag(), rratio/(1.0 - rratio));
         params.emplace(NonDimensional::Upper1d().tag(), 1.0/(1.0 - rratio));
      }
      else
      {
         params.emplace(NonDimensional::Lower1d().tag(), rratio);
         params.emplace(NonDimensional::Upper1d().tag(), 1.0);
      }

      return params;
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

   void ModelBackend::equationInfo(bool& isComplex, SpectralFieldIds& im, SpectralFieldIds& exL, SpectralFieldIds& exNL, SpectralFieldIds& exNS, int& indexMode, const SpectralFieldId& fId, const Resolution& res) const
   {
      // Operators are real
      isComplex = false;

      // Implicit coupled fields
      im = this->implicitFields(fId);

      // Explicit linear terms
      exL = this->explicitLinearFields(fId);

      // Explicit nonlinear terms
      exNL = this->explicitNonlinearFields(fId);

      // Explicit nextstep terms
      exNS.clear();

      // Index mode
      indexMode = static_cast<int>(Equations::CouplingIndexType::SLOWEST_MULTI_RHS);
   }

   void ModelBackend::blockSize(int& tN, int& gN, ArrayI& shift, int& rhs, const SpectralFieldId& fId, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs) const
   {
      auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, eigs.at(0))(0);
      tN = nN;

      int shiftR;
      if(this->mUseGalerkin)
      {
         if(fId == std::make_pair(PhysicalNames::Velocity::id(), FieldComponents::Spectral::TOR) ||
               fId == std::make_pair(PhysicalNames::Temperature::id(), FieldComponents::Spectral::SCALAR))
         {
            shiftR = 2;
         }
         else if(fId == std::make_pair(PhysicalNames::Velocity::id(), FieldComponents::Spectral::POL))
         {
            shiftR = 4;
            //shiftR = 2;
         }
         else
         {
            shiftR = 0;
         }

         gN = (nN - shiftR);
      }
      else
      {
         shiftR = 0;
         gN = nN;
      }

      // Set galerkin shifts
      shift(0) = shiftR;
      shift(1) = 0;
      shift(2) = 0;

      rhs = 1;
   }

   void ModelBackend::operatorInfo(ArrayI& tauN, ArrayI& galN, MatrixI& galShift, ArrayI& rhsCols, ArrayI& sysN, const SpectralFieldId& fId, const Resolution& res, const Equations::Tools::ICoupling& coupling, const BcMap& bcs) const
   {
      // Loop overall matrices/eigs
      for(int idx = 0; idx < tauN.size(); ++idx)
      {
         auto eigs = coupling.getIndexes(res, idx);

         int tN, gN, rhs;
         ArrayI shift(3);

         this->blockSize(tN, gN, shift, rhs, fId, res, eigs, bcs);

         tauN(idx) = tN;
         galN(idx) = gN;
         galShift.row(idx) = shift;
         rhsCols(idx) = rhs;

         // Compute system size
         int sN = 0;
         for(auto f: this->implicitFields(fId))
         {
            this->blockSize(tN, gN, shift, rhs, f, res, eigs, bcs);
            sN += gN;
         }

         if(sN == 0)
         {
            sN = galN(idx);
         }

         sysN(idx) = sN;
      }
   }

   void ModelBackend::implicitBlock(DecoupledZSparse& decMat, const SpectralFieldId& rowId, const SpectralFieldId& colId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const NonDimensional::NdMap& nds) const
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
         SparseSM::Chebyshev::LinearMap::I4Y4SphLapl2 spasm(nN, nN, ri, ro, l);
         //SparseSM::Chebyshev::LinearMap::I2Y2SphLapl spasm(nN, nN, ri, ro, l);
         decMat.real() = spasm.mat();
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
         SparseSM::Chebyshev::LinearMap::I4Y4SphLapl spasm(nN, nN, ri, ro, l);
         //SparseSM::Chebyshev::LinearMap::I2Y2 spasm(nN, nN, ri, ro);
         decMat.real() = spasm.mat();
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

   void ModelBackend::applyGalerkinStencil(DecoupledZSparse& decMat, const SpectralFieldId& rowId, const SpectralFieldId& colId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const
   {
      assert(eigs.size() == 1);

      auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, eigs.at(0))(0);

      auto ri = nds.find(NonDimensional::Lower1d::id())->second->value();
      auto ro = nds.find(NonDimensional::Upper1d::id())->second->value();

      auto stencil = decMat.real();
      this->galerkinStencil(stencil, colId, matIdx, res, eigs, false, bcs, nds);

      auto s = stencil.rows() - stencil.cols();
      SparseSM::Chebyshev::LinearMap::Id qId(nN-s, nN, ri, ro, 0, s);
      decMat.real() = qId.mat()*(decMat.real()*stencil);
   }

   void ModelBackend::applyTau(DecoupledZSparse& decMat, const SpectralFieldId& rowId, const SpectralFieldId& colId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const
   {
      assert(eigs.size() == 1);

      auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, eigs.at(0))(0);

      auto bcId = bcs.find(rowId.first)->second;

      auto ri = nds.find(NonDimensional::Lower1d::id())->second->value();
      auto ro = nds.find(NonDimensional::Upper1d::id())->second->value();

      typedef SparseSM::Chebyshev::LinearMap::Boundary::ICondition::Position Position;

      SparseSM::Chebyshev::LinearMap::Boundary::Operator bcOp(nN, nN, ri, ro);

      if(rowId == std::make_pair(PhysicalNames::Velocity::id(), FieldComponents::Spectral::TOR) && rowId == colId)
      {
         if(bcId == Bc::Name::NoSlip::id())
         {
            bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::Value>(Position::TOP);
            bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::Value>(Position::BOTTOM);
         }
         else if(bcId == Bc::Name::StressFree::id())
         {
            bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::R1D1DivR1>(Position::TOP);
            bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::R1D1DivR1>(Position::BOTTOM);
         }
         else
         {
            throw std::logic_error("Boundary conditions for Velocity Toroidal component not implemented");
         }
      }
      else if(rowId == std::make_pair(PhysicalNames::Velocity::id(), FieldComponents::Spectral::POL) && rowId == colId)
      {
         if(bcId == Bc::Name::NoSlip::id())
         {
            bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::Value>(Position::TOP);
            bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::Value>(Position::BOTTOM);
            bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::D1>(Position::TOP);
            bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::D1>(Position::BOTTOM);
         }
         else if(bcId == Bc::Name::StressFree::id())
         {
            bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::Value>(Position::TOP);
            bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::Value>(Position::BOTTOM);
            bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::D2>(Position::TOP);
            bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::D2>(Position::BOTTOM);
         }
         else
         {
            throw std::logic_error("Boundary conditions for Velocity Poloidal component not implemented");
         }
      }
      else if(rowId == std::make_pair(PhysicalNames::Temperature::id(), FieldComponents::Spectral::SCALAR) && rowId == colId)
      {
         if(bcId == Bc::Name::FixedTemperature::id())
         {
            bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::Value>(Position::TOP);
            bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::Value>(Position::BOTTOM);
         }
         else if(bcId == Bc::Name::FixedFlux::id())
         {
            bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::D1>(Position::TOP);
            bcOp.addRow<SparseSM::Chebyshev::LinearMap::Boundary::D1>(Position::BOTTOM);
         }
         else
         {
            throw std::logic_error("Boundary conditions for Temperature not implemented (" + std::to_string(bcId) + ")");
         }
      }

      decMat.real() += bcOp.mat();
   }

   void ModelBackend::modelMatrix(DecoupledZSparse& rModelMatrix, const std::size_t opId, const Equations::CouplingInformation::FieldId_range imRange, const int matIdx, const std::size_t bcType, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const
   {

      // Time operator
      if(opId == ModelOperator::Time::id())
      {
         bool needStencil = (this->mUseGalerkin && bcType == ModelOperatorBoundary::SolverNoTau::id());
         bool needTau = bcType == ModelOperatorBoundary::SolverHasBc::id();

         for(auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
         {
            this->timeBlock(rModelMatrix, *pRowId, matIdx, res, eigs, nds);

            // Apply boundary condition
            if(needStencil)
            {
               this->applyGalerkinStencil(rModelMatrix, *pRowId, *pRowId, matIdx, res, eigs, bcs, nds);
            }
            else if(needTau)
            {
               this->applyTau(rModelMatrix, *pRowId, *pRowId, matIdx, res, eigs, bcs, nds);
            }
         }
      }
      // Linear operator
      else if(opId == ModelOperator::ImplicitLinear::id())
      {
         bool needStencil = (this->mUseGalerkin && bcType == ModelOperatorBoundary::SolverNoTau::id());
         bool needTau = bcType == ModelOperatorBoundary::SolverHasBc::id();

         for(auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
         {
            for(auto pColId = imRange.first; pColId != imRange.second; pColId++)
            {
               this->implicitBlock(rModelMatrix, *pRowId, *pColId, matIdx, res, eigs, nds);

               // Apply boundary condition
               if(needStencil)
               {
                  this->applyGalerkinStencil(rModelMatrix, *pRowId, *pColId, matIdx, res, eigs, bcs, nds);
               }
               else if(needTau)
               {
                  this->applyTau(rModelMatrix, *pRowId, *pColId, matIdx, res, eigs, bcs, nds);
               }
            }
         }
      }
      // Boundary operator
      else if(opId == ModelOperator::Boundary::id())
      {
         bool needStencil = this->mUseGalerkin;
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
                  this->applyGalerkinStencil(rModelMatrix, *pRowId, *pColId, matIdx, res, eigs, bcs, nds);
               }
               else if(needTau)
               {
                  this->applyTau(rModelMatrix, *pRowId, *pColId, matIdx, res, eigs, bcs, nds);
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
      throw std::logic_error("Not yet implemented!");
      
#if 0
      assert(eigs.size() == 1);
      int l = eigs.at(0);

      auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, eigs.at(0))(0);

      auto bcId = bcs.find(fieldId.first)->second;

      int s = 0;
      if(fieldId == std::make_pair(PhysicalNames::Velocity::id(), FieldComponents::Spectral::TOR))
      {
         s = 1;
         if(bcId == Bc::Name::NoSlip::id())
         {
            SparseSM::Chebyshev::LinearMap::Stencil::Value bc(nN, nN-1, a, b, l);
            mat = bc.mat();
         }
         else if(bcId == Bc::Name::StressFree::id())
         {
            SparseSM::Chebyshev::LinearMap::Stencil::R1D1DivR1 bc(nN, nN-1, a, b, l);
            mat = bc.mat();
         }
         else
         {
            throw std::logic_error("Galerkin boundary conditions for Velocity Toroidal component not implemented");
         }
      }
      else if(fieldId == std::make_pair(PhysicalNames::Velocity::id(), FieldComponents::Spectral::POL))
      {
         s = 2;
         if(bcId == Bc::Name::NoSlip::id())
         {
            SparseSM::Chebyshev::LinearMap::Stencil::ValueD1 bc(nN, nN-2, a, b, l);
            mat = bc.mat();
         }
         else if(bcId == Bc::Name::StressFree::id())
         {
            SparseSM::Chebyshev::LinearMap::Stencil::ValueD2 bc(nN, nN-2, a, b, l);
            mat = bc.mat();
         }
         else
         {
            throw std::logic_error("Galerin boundary conditions for Velocity Poloidal component not implemented");
         }
      }
      else if(fieldId == std::make_pair(PhysicalNames::Temperature::id(), FieldComponents::Spectral::SCALAR))
      {
         s = 1;
         if(bcId == Bc::Name::FixedTemperature::id())
         {
            SparseSM::Chebyshev::LinearMap::Stencil::Value bc(nN, nN-1, a, b, l);
            mat = bc.mat();
         }
         else if(bcId == Bc::Name::FixedFlux::id())
         {
            SparseSM::Chebyshev::LinearMap::Stencil::D1 bc(nN, nN-1, a, b, l);
            mat = bc.mat();
         }
         else
         {
            throw std::logic_error("Galerkin boundary conditions for Temperature not implemented");
         }
      }

      if(makeSquare)
      {
         SparseSM::Chebyshev::LinearMap::Id qId(nN-s, nN, a, b, l);
         mat = qId.mat()*mat;
      }
#endif
   }

   void ModelBackend::explicitBlock(DecoupledZSparse& decMat, const SpectralFieldId& rowId, const std::size_t opId,  const SpectralFieldId colId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const
   {
      assert(eigs.size() == 1);
      int l = eigs.at(0);

      auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, eigs.at(0))(0);

      auto ri = nds.find(NonDimensional::Lower1d::id())->second->value();
      auto ro = nds.find(NonDimensional::Upper1d::id())->second->value();

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
   }

   MHDFloat ModelBackend::effectiveRa(const NonDimensional::NdMap& nds) const
   {
      auto effRa = nds.find(NonDimensional::Rayleigh::id())->second->value();
      auto ro = nds.find(NonDimensional::Upper1d::id())->second->value();

      // Scaled on gap width
      if(ro != 1.0)
      {
         auto T = 1.0/nds.find(NonDimensional::Ekman::id())->second->value();
         effRa *= T/ro;
      }

      return effRa;
   }

   MHDFloat ModelBackend::effectiveBg(const NonDimensional::NdMap& nds) const
   {
      MHDFloat effBg = 1.0;
      auto ro = nds.find(NonDimensional::Upper1d::id())->second->value();
      auto rratio = nds.find(NonDimensional::RRatio::id())->second->value();
      auto heatingMode = nds.find(NonDimensional::Heating::id())->second->value();

      if(ro == 1.0)
      {
         // Nothing
         //
      }
      // gap width and internal heating
      else if(heatingMode == 0)
      {
         effBg = 2.0/(ro*(1.0 + rratio));
      }
      // gap width and differential heating
      else if(heatingMode == 1)
      {
         effBg = ro*ro*rratio;
      }

      return effBg;
   }


} // Explicit
} // RTC
} // Shell
} // Boussinesq
} // Model
} // QuICC
