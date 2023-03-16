/**
 * @file ModelBackend.cpp
 * @brief Source of the interface for model backend
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "QuICC/Model/Boussinesq/Shell/RTC/Implicit/ModelBackend.hpp"
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
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2Y1.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2Y2.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2Y2D1.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2Y3.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2Y2SphLapl.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2Y3SphLapl.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I4Y3.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I4Y4D1.hpp"
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
#include "QuICC/Math/Constants.hpp"
#include "QuICC/Precision.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Shell {

namespace RTC {

namespace Implicit {

   ModelBackend::ModelBackend()
      : IRTCBackend(),
#ifdef QUICC_TRANSFORM_CHEBYSHEV_TRUNCATE_QI
      mcTruncateQI(true)
#else
      mcTruncateQI(false)
#endif // QUICC_TRANSFORM_CHEBYSHEV_TRUNCATE_QI
   {
      this->enableSplitEquation(false);
   }

   ModelBackend::SpectralFieldIds ModelBackend::implicitFields(const SpectralFieldId& fId) const
   {
      SpectralFieldId velTor = std::make_pair(PhysicalNames::Velocity::id(), FieldComponents::Spectral::TOR);
      SpectralFieldId velPol = std::make_pair(PhysicalNames::Velocity::id(), FieldComponents::Spectral::POL);
      SpectralFieldId temp = std::make_pair(PhysicalNames::Temperature::id(), FieldComponents::Spectral::SCALAR);
      SpectralFieldIds fields = {velTor, velPol, temp};

      return fields;
   }

   ModelBackend::SpectralFieldIds ModelBackend::explicitLinearFields(const SpectralFieldId& fId) const
   {
      SpectralFieldIds fields;
      if(fId == std::make_pair(PhysicalNames::Temperature::id(),FieldComponents::Spectral::SCALAR))
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
      // Operators are complex
      info.isComplex = true;

      // Use split operators
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
      info.indexMode = static_cast<int>(Equations::CouplingIndexType::SLOWEST_SINGLE_RHS);
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

   void ModelBackend::implicitBlock(DecoupledZSparse& decMat, const SpectralFieldId& rowId, const SpectralFieldId& colId, const int matIdx, const std::size_t bcType, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds, const bool isSplitEquation) const
   {
      assert(eigs.size() == 1);
      int m = eigs.at(0);

      // Compute system size
      const auto sysInfo = systemInfo(rowId, colId, m, res, bcs, this->useGalerkin(), false);
      const auto& sysN = sysInfo.systemSize;
      const auto& baseRowShift = sysInfo.startRow;
      const auto& baseColShift = sysInfo.startCol;
      auto nL = res.counter().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL, m);

      auto ri = nds.find(NonDimensional::Lower1d::id())->second->value();
      auto ro = nds.find(NonDimensional::Upper1d::id())->second->value();

      if(decMat.real().size() == 0)
      {
         decMat.real().resize(sysN,sysN);
         decMat.imag().resize(sysN,sysN);
      }
      assert(decMat.real().rows() == sysN);
      assert(decMat.real().cols() == sysN);
      assert(decMat.imag().rows() == sysN);
      assert(decMat.imag().cols() == sysN);

      int tN, gN, rhs;
      ArrayI shift(3);

      int rowShift = baseRowShift;
      int colShift = baseColShift;
      if(rowId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::TOR))
      {
         if(rowId == colId)
         {
            const auto Ek = nds.find(NonDimensional::Ekman::id())->second->value();
            const auto T = 1.0/Ek;

            for(int l = m; l < nL; l++)
            {
               auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);
               if(l > 0)
               {
                  const auto dl = static_cast<MHDFloat>(l);
                  const auto invlapl = 1.0/(dl*(dl + 1.0));
                  SparseSM::Chebyshev::LinearMap::I2Y2SphLapl i2r2lapl(nN, nN, ri, ro, l);
                  SparseMatrix bMat = i2r2lapl.mat(); 
                  if(this->useGalerkin())
                  {
                     this->applyGalerkinStencil(bMat, rowId, colId, l, res, bcs, nds);
                  }
                  this->addBlock(decMat.real(), bMat, rowShift, colShift);

                  SparseSM::Chebyshev::LinearMap::I2Y2 i2r2(nN, nN, ri, ro);
                  bMat = m*T*invlapl*i2r2.mat();
                  if(this->useGalerkin())
                  {
                     this->applyGalerkinStencil(bMat, rowId, colId, l, res, bcs, nds);
                  }
                  this->addBlock(decMat.imag(), bMat, rowShift, colShift);
               }
               this->blockInfo(tN, gN, shift, rhs, rowId, res, l, bcs);
               rowShift += gN;

               this->blockInfo(tN, gN, shift, rhs, colId, res, l, bcs);
               colShift += gN;
            }
         }
         else if(colId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::POL))
         {
            auto coriolis = [](const int l, const int m){
               return (l - MHD_MP(1.0))*(l + MHD_MP(1.0))*precision::sqrt(((l - m)*(l + m))/((MHD_MP(2.0)*l - MHD_MP(1.0))*(MHD_MP(2.0)*l + MHD_MP(1.0))));
            };

            const auto Ek = nds.find(NonDimensional::Ekman::id())->second->value();
            const auto T = 1.0/Ek;

            this->blockInfo(tN, gN, shift, rhs, rowId, res, m, bcs);

            rowShift = baseRowShift + gN;
            colShift = baseColShift;
            for(int l = m+1; l < nL; l++)
            {
               auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);
               if(l > 0)
               {
                  const auto dl = static_cast<QuICC::internal::MHDFloat>(l);
                  const auto invlapl = 1.0/(dl*(dl + 1.0));
                  SparseSM::Chebyshev::LinearMap::I2Y1 cor_r(nN, nN, ri, ro);
                  auto norm = (dl - MHD_MP(1.0))*coriolis(l, m);
                  SparseMatrix bMat = -static_cast<MHDFloat>(norm*T*invlapl)*cor_r.mat();

                  SparseSM::Chebyshev::LinearMap::I2Y2D1 cordr(nN, nN, ri, ro);
                  norm = -coriolis(l, m);
                  bMat += -static_cast<MHDFloat>(norm*T*invlapl)*cordr.mat();
                  if(this->useGalerkin())
                  {
                     this->applyGalerkinStencil(bMat, rowId, colId, l, res, bcs, nds);
                  }
                  this->addBlock(decMat.real(), bMat, rowShift, colShift);

               }

               this->blockInfo(tN, gN, shift, rhs, rowId, res, l, bcs);
               rowShift += gN;

               this->blockInfo(tN, gN, shift, rhs, colId, res, l, bcs);
               colShift += gN;
            }

            this->blockInfo(tN, gN, shift, rhs, colId, res, m, bcs);
            rowShift = baseRowShift;
            colShift = baseColShift + gN;
            for(int l = m; l < nL - 1; l++)
            {
               auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);
               if(l > 0)
               {
                  const auto dl = static_cast<QuICC::internal::MHDFloat>(l);
                  const auto invlapl = MHD_MP(1.0)/(dl*(dl + MHD_MP(1.0)));
                  SparseSM::Chebyshev::LinearMap::I2Y1 cor_r(nN, nN, ri, ro);
                  auto norm = -(dl + MHD_MP(2.0))*coriolis(l+1, m);
                  SparseMatrix bMat = -static_cast<MHDFloat>(norm*T*invlapl)*cor_r.mat();

                  SparseSM::Chebyshev::LinearMap::I2Y2D1 cordr(nN, nN, ri, ro);
                  norm = -coriolis(l+1, m);
                  bMat += -static_cast<MHDFloat>(norm*T*invlapl)*cordr.mat();
                  if(this->useGalerkin())
                  {
                     this->applyGalerkinStencil(bMat, rowId, colId, l, res, bcs, nds);
                  }
                  this->addBlock(decMat.real(), bMat, rowShift, colShift);

               }

               this->blockInfo(tN, gN, shift, rhs, rowId, res, l, bcs);
               rowShift += gN;

               this->blockInfo(tN, gN, shift, rhs, colId, res, l, bcs);
               colShift += gN;
            }
         }
      }
      else if(rowId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::POL))
      {
         if(rowId == colId)
         {
            const auto Ek = nds.find(NonDimensional::Ekman::id())->second->value();
            const auto T = 1.0/Ek;

            for(int l = m; l < nL; l++)
            {
               auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);
               if(l > 0)
               {
                  const auto dl = static_cast<QuICC::internal::MHDFloat>(l);
                  const auto invlapl = MHD_MP(1.0)/(dl*(dl + MHD_MP(1.0)));
                  SparseSM::Chebyshev::LinearMap::I4Y4SphLapl2 diffusion(nN, nN, ri, ro, l);
                  SparseMatrix bMat = diffusion.mat();
                  if(this->useGalerkin())
                  {
                     this->applyGalerkinStencil(bMat, rowId, colId, l, res, bcs, nds);
                  }

                  this->addBlock(decMat.real(), bMat, rowShift, colShift);
                  SparseSM::Chebyshev::LinearMap::I4Y4SphLapl coriolis(nN, nN, ri, ro, l);
                  bMat = static_cast<MHDFloat>(m*T*invlapl)*coriolis.mat();
                  if(this->useGalerkin())
                  {
                     this->applyGalerkinStencil(bMat, rowId, colId, l, res, bcs, nds);
                  }
                  this->addBlock(decMat.imag(), bMat, rowShift, colShift);
               }
               this->blockInfo(tN, gN, shift, rhs, rowId, res, l, bcs);
               rowShift += gN;

               this->blockInfo(tN, gN, shift, rhs, colId, res, l, bcs);
               colShift += gN;
            }
         }
         else if(colId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::TOR))
         {
            auto coriolis = [](const int l, const int m){
               return (l - MHD_MP(1.0))*(l + MHD_MP(1.0))*precision::sqrt(((l - m)*(l + m))/((MHD_MP(2.0)*l - MHD_MP(1.0))*(MHD_MP(2.0)*l + MHD_MP(1.0))));
            };

            const auto Ek = nds.find(NonDimensional::Ekman::id())->second->value();
            const auto T = 1.0/Ek;

            this->blockInfo(tN, gN, shift, rhs, rowId, res, m, bcs);
            rowShift = baseRowShift + gN;
            colShift = baseColShift;
            for(int l = m+1; l < nL; l++)
            {
               auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);
               if(l > 0)
               {
                  const auto dl = static_cast<MHDFloat>(l);
                  const auto invlapl = 1.0/(dl*(dl + 1.0));
                  SparseSM::Chebyshev::LinearMap::I4Y3 cor_r(nN, nN, ri, ro);
                  auto norm = (dl - 1.0)*coriolis(l, m);
                  SparseMatrix bMat = static_cast<MHDFloat>(norm*T*invlapl)*cor_r.mat();
                  SparseSM::Chebyshev::LinearMap::I4Y4D1 cordr(nN, nN, ri, ro);
                  norm = -coriolis(l, m);
                  bMat += static_cast<MHDFloat>(norm*T*invlapl)*cordr.mat();
                  if(this->useGalerkin())
                  {
                     this->applyGalerkinStencil(bMat, rowId, colId, l, res, bcs, nds);
                  }
                  this->addBlock(decMat.real(), bMat, rowShift, colShift);
               }
               this->blockInfo(tN, gN, shift, rhs, rowId, res, l, bcs);
               rowShift += gN;

               this->blockInfo(tN, gN, shift, rhs, colId, res, l, bcs);
               colShift += gN;
            }

            this->blockInfo(tN, gN, shift, rhs, colId, res, m, bcs);
            rowShift = baseRowShift;
            colShift = baseColShift + gN;
            for(int l = m; l < nL - 1; l++)
            {
               auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);
               if(l > 0)
               {
                  const auto dl = static_cast<MHDFloat>(l);
                  const auto invlapl = 1.0/(dl*(dl + 1.0));
                  SparseSM::Chebyshev::LinearMap::I4Y3 cor_r(nN, nN, ri, ro);
                  auto norm = -(dl + 2.0)*coriolis(l+1, m);
                  SparseMatrix bMat = static_cast<MHDFloat>(norm*T*invlapl)*cor_r.mat();
                  SparseSM::Chebyshev::LinearMap::I4Y4D1 cordr(nN, nN, ri, ro);
                  norm = -coriolis(l+1, m);
                  bMat += static_cast<MHDFloat>(norm*T*invlapl)*cordr.mat();
                  if(this->useGalerkin())
                  {
                     this->applyGalerkinStencil(bMat, rowId, colId, l, res, bcs, nds);
                  }
                  this->addBlock(decMat.real(), bMat, rowShift, colShift);
               }
               this->blockInfo(tN, gN, shift, rhs, rowId, res, l, bcs);
               rowShift += gN;

               this->blockInfo(tN, gN, shift, rhs, colId, res, l, bcs);
               colShift += gN;
            }
         }
         else if(colId == std::make_pair(PhysicalNames::Temperature::id(),FieldComponents::Spectral::SCALAR))
         {
            const auto Ra = this->effectiveRa(nds);
            for(int l = m; l < nL; l++)
            {
               auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);
               if(l > 0)
               {
                  SparseSM::Chebyshev::LinearMap::I4Y4 i4r4(nN, nN, ri, ro);
                  SparseMatrix bMat = -Ra*i4r4.mat();
                  if(this->useGalerkin())
                  {
                     this->applyGalerkinStencil(bMat, rowId, colId, l, res, bcs, nds);
                  }
                  this->addBlock(decMat.real(), bMat, rowShift, colShift);
               }
               this->blockInfo(tN, gN, shift, rhs, rowId, res, l, bcs);
               rowShift += gN;

               this->blockInfo(tN, gN, shift, rhs, colId, res, l, bcs);
               colShift += gN;
            }
         }
      }
      else if(rowId == std::make_pair(PhysicalNames::Temperature::id(), FieldComponents::Spectral::SCALAR))
      {
         if(rowId == colId)
         {
            auto Pr = nds.find(NonDimensional::Prandtl::id())->second->value();
            auto heatingMode = nds.find(NonDimensional::Heating::id())->second->value();

            for(int l = m; l < nL; l++)
            {
               auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);
               SparseMatrix bMat;
               if(heatingMode == 0)
               {
                  SparseSM::Chebyshev::LinearMap::I2Y2SphLapl spasm(nN, nN, ri, ro, l);
                  bMat = (1.0/Pr)*spasm.mat();
               }
               else
               {
                  SparseSM::Chebyshev::LinearMap::I2Y3SphLapl spasm(nN, nN, ri, ro, l);
                  bMat = (1.0/Pr)*spasm.mat();
               }
               if(this->useGalerkin())
               {
                  this->applyGalerkinStencil(bMat, rowId, colId, l, res, bcs, nds);
               }
               this->addBlock(decMat.real(), bMat, rowShift, colShift);

               this->blockInfo(tN, gN, shift, rhs, rowId, res, l, bcs);
               rowShift += gN;

               this->blockInfo(tN, gN, shift, rhs, colId, res, l, bcs);
               colShift += gN;
            }
         }
      }
      else
      {
         //throw std::logic_error("Equations are not setup properly");
      }
   }

   void ModelBackend::timeBlock(DecoupledZSparse& decMat, const SpectralFieldId& fieldId, const int matIdx, const std::size_t bcType, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const
   {
      assert(eigs.size() == 1);
      int m = eigs.at(0);

      // Compute system size
      const auto sysInfo = systemInfo(fieldId, fieldId, m, res, bcs, this->useGalerkin(), false);
      const auto& sysN = sysInfo.systemSize;
      const auto& baseShift = sysInfo.startRow;
      auto nL = res.counter().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL, m);

      auto ri = nds.find(NonDimensional::Lower1d::id())->second->value();
      auto ro = nds.find(NonDimensional::Upper1d::id())->second->value();

      // Resize matrices the first time
      if(decMat.real().size() == 0)
      {
         decMat.real().resize(sysN,sysN);
         decMat.imag().resize(sysN,sysN);
      }
      assert(decMat.real().rows() == sysN);
      assert(decMat.imag().rows() == sysN);

      int tN, gN, rhs;
      ArrayI shift(3);

      if(fieldId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::TOR))
      {
         int rowShift = baseShift;
         for(int l = m; l < nL; l++)
         {
            auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);
            this->blockInfo(tN, gN, shift, rhs, fieldId, res, l, bcs);
            SparseMatrix bMat;
            if(l > 0)
            {
               SparseSM::Chebyshev::LinearMap::I2Y2 spasm(nN, nN, ri, ro);
               bMat = spasm.mat();
               if(this->useGalerkin())
               {
                  this->applyGalerkinStencil(bMat, fieldId, fieldId, l, res, bcs, nds);
               }
            }
            else
            {
               SparseSM::Chebyshev::LinearMap::Id qid(gN, gN, ri, ro);
               bMat = qid.mat();
            }
            this->addBlock(decMat.real(), bMat, rowShift, rowShift);
            rowShift += gN;
         }
      }
      else if(fieldId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::POL))
      {
         int rowShift = baseShift;
         for(int l = m; l < nL; l++)
         {
            auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);
            this->blockInfo(tN, gN, shift, rhs, fieldId, res, l, bcs);
            SparseMatrix bMat;
            if(l > 0)
            {
               SparseSM::Chebyshev::LinearMap::I4Y4SphLapl spasm(nN, nN, ri, ro, l);
               bMat = spasm.mat();
               if(this->useGalerkin())
               {
                  this->applyGalerkinStencil(bMat, fieldId, fieldId, l, res, bcs, nds);
               }
            }
            else
            {
               SparseSM::Chebyshev::LinearMap::Id qid(gN, gN, ri, ro);
               bMat = qid.mat();
            }
            this->addBlock(decMat.real(), bMat, rowShift, rowShift);
            rowShift += gN;
         }
      }
      else if(fieldId == std::make_pair(PhysicalNames::Temperature::id(), FieldComponents::Spectral::SCALAR))
      {
         auto heatingMode = nds.find(NonDimensional::Heating::id())->second->value();

         int rowShift = baseShift;
         for(int l = m; l < nL; l++)
         {
            auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);
            SparseMatrix bMat;
            if(heatingMode == 0)
            {
               SparseSM::Chebyshev::LinearMap::I2Y2 spasm(nN, nN, ri, ro);
               bMat = spasm.mat();
            }
            else
            {
               SparseSM::Chebyshev::LinearMap::I2Y3 spasm(nN, nN, ri, ro);
               bMat = spasm.mat();
            }
            if(this->useGalerkin())
            {
               this->applyGalerkinStencil(bMat, fieldId, fieldId, l, res, bcs, nds);
            }
            this->addBlock(decMat.real(), bMat, rowShift, rowShift);
            this->blockInfo(tN, gN, shift, rhs, fieldId, res, l, bcs);
            rowShift += gN;
         }
      }
   }

   void ModelBackend::boundaryBlock(DecoupledZSparse& decMat, const SpectralFieldId& rowId, const SpectralFieldId& colId, const int matIdx, const std::size_t bcType, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds, const bool isSplit) const
   {
      bool needTau = (bcType == ModelOperatorBoundary::SolverHasBc::id());

      assert(eigs.size() == 1);
      int m = eigs.at(0);

      // Compute system size
      const auto sysInfo = systemInfo(rowId, colId, m, res, bcs, this->useGalerkin(), false);
      const auto sysN = sysInfo.systemSize;
      const auto baseRowShift = sysInfo.startRow;
      const auto baseColShift = sysInfo.startCol;
      auto nL = res.counter().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL, m);

      if(decMat.real().size() == 0)
      {
         decMat.real().resize(sysN,sysN);
         decMat.imag().resize(sysN,sysN);
      }
      assert(decMat.real().rows() == sysN);
      assert(decMat.real().cols() == sysN);
      assert(decMat.imag().rows() == sysN);
      assert(decMat.imag().cols() == sysN);

      int tN, gN, rhs;
      ArrayI shift(3);

      int rowShift = baseRowShift;
      int colShift = baseColShift;
      // Apply boundary condition
      if(this->useGalerkin())
      {
         for(int l = m; l < nL; l++)
         {
            auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);
            SparseMatrix mat(nN, nN);
            this->applyGalerkinStencil(mat, rowId, colId, l, res, bcs, nds);
            this->addBlock(decMat.real(), mat, rowShift, colShift);

            this->blockInfo(tN, gN, shift, rhs, rowId, res, l, bcs);
            rowShift += gN;

            this->blockInfo(tN, gN, shift, rhs, colId, res, l, bcs);
            colShift += gN;
         }
      }
      else if(needTau)
      {
         int minL = m;
         if(m == 0 && (rowId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::TOR) || rowId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::POL)))
         {
            minL = m + 1;
            auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, m)(0);
            rowShift += nN;
            colShift += nN;
         }
         for(int l = minL; l < nL; l++)
         {
            auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);
            SparseMatrix mat(nN, nN);
            this->applyTau(mat, rowId, colId, l, res, bcs, nds, isSplit);
            this->addBlock(decMat.real(), mat, rowShift, colShift);

            rowShift += nN;
            colShift += nN;
         }
      }
   }

   void ModelBackend::modelMatrix(DecoupledZSparse& rModelMatrix, const std::size_t opId, const Equations::CouplingInformation::FieldId_range imRange, const int matIdx, const std::size_t bcType, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const
   {
      // Time operator
      if(opId == ModelOperator::Time::id())
      {
         for(auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
         {
            this->timeBlock(rModelMatrix, *pRowId, matIdx, bcType, res, eigs, bcs, nds);
         }
      }
      // Linear operator
      else if(opId == ModelOperator::ImplicitLinear::id() || opId == ModelOperator::SplitImplicitLinear::id())
      {
         bool isSplit = (opId == ModelOperator::SplitImplicitLinear::id());

         for(auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
         {
            for(auto pColId = imRange.first; pColId != imRange.second; pColId++)
            {
               this->implicitBlock(rModelMatrix, *pRowId, *pColId, matIdx, bcType, res, eigs, bcs, nds, isSplit);
            }
         }
      }
      // Boundary operator
      else if(opId == ModelOperator::Boundary::id() || opId == ModelOperator::SplitBoundary::id())
      {
         bool isSplit = (opId == ModelOperator::SplitBoundary::id());

         for(auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
         {
            for(auto pColId = imRange.first; pColId != imRange.second; pColId++)
            {
               this->boundaryBlock(rModelMatrix, *pRowId, *pColId, matIdx, bcType, res, eigs, bcs, nds, isSplit);
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
      assert(this->useGalerkin());
      assert(eigs.size() == 1);
      int m = eigs.at(0);

      // Compute system size
      const auto sysRows = systemInfo(fieldId, fieldId, m, res, bcs, makeSquare, false).blockRows;
      const auto sysCols = systemInfo(fieldId, fieldId, m, res, bcs, true, false).blockCols;

      auto nL = res.counter().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL, m);

      if(mat.size() == 0)
      {
         mat.resize(sysRows,sysCols);
      }
      assert(mat.rows() == sysRows);
      assert(mat.cols() == sysCols);

      assert(eigs.size() == 1);

      int rowShift = 0;
      int colShift = 0;
      for(int l = m; l < nL; l++)
      {
         SparseMatrix S;
         this->stencil(S, fieldId, l, res, makeSquare, bcs, nds);
         this->addBlock(mat, S, rowShift, colShift);

         rowShift += S.rows();
         colShift += S.cols();
      }
   }

   void ModelBackend::explicitBlock(DecoupledZSparse& decMat, const SpectralFieldId& rowId, const std::size_t opId,  const SpectralFieldId colId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const
   {
      assert(eigs.size() == 1);
      int m = eigs.at(0);

      int tN, gN, rhs;
      ArrayI shift(3);
      this->blockInfo(tN, gN, shift, rhs, colId, res, m, bcs);

      // Compute system size
      const auto sysInfo = systemInfo(rowId, colId, m, res, bcs, false, this->useGalerkin());
      const auto sysRows = sysInfo.blockRows;
      const auto sysCols = sysInfo.blockCols;
      const auto baseShift = 0;
      auto nL = res.counter().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL, m);

      auto ri = nds.find(NonDimensional::Lower1d::id())->second->value();
      auto ro = nds.find(NonDimensional::Upper1d::id())->second->value();

      if(decMat.real().size() == 0)
      {
         decMat.real().resize(sysRows,sysCols);
         decMat.imag().resize(sysRows,sysCols);
      }
      assert(decMat.real().rows() == sysRows);
      assert(decMat.real().cols() == sysCols);
      assert(decMat.imag().rows() == sysRows);
      assert(decMat.imag().cols() == sysCols);

      // Explicit linear operator
      if(opId == ModelOperator::ExplicitLinear::id())
      {
         if(rowId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::POL) &&
               colId == std::make_pair(PhysicalNames::Temperature::id(),FieldComponents::Spectral::SCALAR))
         {
            auto Ra = effectiveRa(nds);

            int rowShift = baseShift;
            int colShift = baseShift;
            for(int l = m; l < nL; l++)
            {
               auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);
               SparseSM::Chebyshev::LinearMap::I4Y4 spasm(nN, nN, ri, ro);
               SparseMatrix bMat = Ra*spasm.mat();
               if(this->useGalerkin())
               {
                  this->blockInfo(tN, gN, shift, rhs, rowId, res, eigs.at(0), bcs);
                  SparseSM::Chebyshev::LinearMap::Id qId(nN-shift(0), nN, ri, ro, 0, shift(0));
                  bMat = qId.mat()*bMat;
               }
               this->addBlock(decMat.real(), bMat, rowShift, colShift);
               //this->addBlock(decMat.real(), spasm.mat(), rowShift, colShift, Ra);

               this->blockInfo(tN, gN, shift, rhs, rowId, res, eigs.at(0), bcs);
               rowShift += gN;

               this->blockInfo(tN, gN, shift, rhs, colId, res, eigs.at(0), bcs);
               colShift += tN;
            }
         }
         else if(rowId == std::make_pair(PhysicalNames::Temperature::id(), FieldComponents::Spectral::SCALAR) && 
               colId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::POL))
         {
            auto heatingMode = nds.find(NonDimensional::Heating::id())->second->value();
            auto bg = effectiveBg(nds);

            int rowShift = baseShift;
            int colShift = baseShift;
            for(int l = m; l < nL; l++)
            {
               auto dl = static_cast<MHDFloat>(l);
               auto ll1 = dl*(dl + 1.0);
               auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);
               SparseMatrix bMat;
               if(heatingMode == 0)
               {
                  SparseSM::Chebyshev::LinearMap::I2Y2 spasm(nN, nN, ri, ro);
                  bMat = -bg*ll1*spasm.mat();
                  if(this->useGalerkin())
                  {
                     this->blockInfo(tN, gN, shift, rhs, rowId, res, eigs.at(0), bcs);
                     SparseSM::Chebyshev::LinearMap::Id qId(nN-shift(0), nN, ri, ro, 0, shift(0));
                     bMat = qId.mat()*bMat;
                  }
                  //this->addBlock(decMat.real(), spasm.mat(), rowShift, colShift, -bg*ll1);
               }
               else
               {
                  SparseSM::Chebyshev::LinearMap::I2 spasm(nN, nN, ri, ro);
                  bMat = -bg*ll1*spasm.mat();
                  if(this->useGalerkin())
                  {
                     this->blockInfo(tN, gN, shift, rhs, rowId, res, eigs.at(0), bcs);
                     SparseSM::Chebyshev::LinearMap::Id qId(nN-shift(0), nN, ri, ro, 0, shift(0));
                     bMat = qId.mat()*bMat;
                  }
                  //this->addBlock(decMat.real(), spasm.mat(), rowShift, colShift, -bg*ll1);
               }
               this->addBlock(decMat.real(), bMat, rowShift, colShift);

               this->blockInfo(tN, gN, shift, rhs, rowId, res, eigs.at(0), bcs);
               rowShift += gN;

               this->blockInfo(tN, gN, shift, rhs, colId, res, eigs.at(0), bcs);
               colShift += tN;
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

            int rowShift = baseShift;
            int colShift = baseShift;
            for(int l = m; l < nL; l++)
            {
               auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);
               SparseMatrix bMat;
               if(heatingMode == 0)
               {
                  SparseSM::Chebyshev::LinearMap::I2Y2 spasm(nN, nN, ri, ro);
                  bMat = spasm.mat();
                  if(this->useGalerkin())
                  {
                     this->blockInfo(tN, gN, shift, rhs, rowId, res, eigs.at(0), bcs);
                     SparseSM::Chebyshev::LinearMap::Id qId(nN-shift(0), nN, ri, ro, 0, shift(0));
                     bMat = qId.mat()*bMat;
                  }
                  //this->addBlock(decMat.real(), spasm.mat(), rowShift, colShift);
               }
               else
               {
                  SparseSM::Chebyshev::LinearMap::I2Y3 spasm(nN, nN, ri, ro);
                  bMat = spasm.mat();
                  if(this->useGalerkin())
                  {
                     this->blockInfo(tN, gN, shift, rhs, rowId, res, eigs.at(0), bcs);
                     SparseSM::Chebyshev::LinearMap::Id qId(nN-shift(0), nN, ri, ro, 0, shift(0));
                     bMat = qId.mat()*bMat;
                  }
                  //this->addBlock(decMat.real(), spasm.mat(), rowShift, colShift);
               }
               this->addBlock(decMat.real(), bMat, rowShift, colShift);

               this->blockInfo(tN, gN, shift, rhs, rowId, res, eigs.at(0), bcs);
               rowShift += gN;

               this->blockInfo(tN, gN, shift, rhs, colId, res, eigs.at(0), bcs);
               colShift += tN;
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

   void ModelBackend::addBlock(SparseMatrix& mat, const SparseMatrix& block, const int rowShift, const int colShift, const MHDFloat coeff) const
   {
      std::vector<Eigen::Triplet<MHDFloat> > triplets;
      triplets.reserve(block.nonZeros());
      for(int k = 0; k < block.outerSize(); ++k)
      {
         for(SparseMatrix::InnerIterator it(block, k); it; ++it)
         {
            triplets.emplace_back(Eigen::Triplet<MHDFloat>(it.row() + rowShift, it.col() + colShift, coeff*it.value()));
         }
      }
      SparseMatrix full(mat.rows(), mat.cols());
      full.setFromTriplets(triplets.begin(), triplets.end());
      mat += full;
   }

   int ModelBackend::blockSize(const SpectralFieldId& fId, const int m, const Resolution& res, const BcMap& bcs, const bool isGalerkin) const
   {
      auto nL = res.counter().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL, m);

      // Compute size
      auto s = 0;
      for(int l = m; l < nL; l++)
      {
         int tN, gN, rhs;
         ArrayI shift(3);
         this->blockInfo(tN, gN, shift, rhs, fId, res, l, bcs);
         if(isGalerkin)
         {
            s += gN;
         }
         else
         {
            s += tN;
         }
      }

      return s;
   }

   std::pair<int, int> ModelBackend::blockShape(const SpectralFieldId& rowId, const SpectralFieldId& colId, const int m, const Resolution& res, const BcMap& bcs, const bool isGalerkin, const bool dropRows) const
   {
      // Compute number of rows
      auto rows = this->blockSize(rowId, m, res, bcs, isGalerkin || dropRows);

      // Compute number of cols
      int cols = this->blockSize(colId, m, res, bcs, isGalerkin);

      return std::make_pair(rows, cols);
   }

   internal::SystemInfo ModelBackend::systemInfo(const SpectralFieldId& rowId, const SpectralFieldId& colId, const int m, const Resolution& res, const BcMap& bcs, const bool isGalerkin, const bool dropRows) const
   {
      auto shape = this->blockShape(rowId, colId, m, res, bcs, isGalerkin, dropRows);

      int sysN = 0;
      bool rowCount = true;
      bool colCount = true;
      int rowIdx = 0;
      int colIdx = 0;
      const auto& fields = this->implicitFields(rowId);
      for(auto it = fields.begin(); it != fields.end(); ++it)
      {
         int s = this->blockSize(*it, m, res, bcs, isGalerkin);
         sysN += s;

         // Get block index of rowId
         if(rowCount && rowId != *it)
         {
            rowIdx += s;
         }
         else if(rowId == *it)
         {
            rowCount = false;
         }

         // Get block index of colId
         if(colCount && colId != *it)
         {
            colIdx += s;
         }
         else if(colId == *it)
         {
            colCount = false;
         }
      }

      internal::SystemInfo info(sysN, shape.first, shape.second, rowIdx, colIdx);
      return info;
   }


} // Implicit
} // RTC
} // Shell
} // Boussinesq
} // Model
} // QuICC
