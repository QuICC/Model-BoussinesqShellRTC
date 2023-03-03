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
#include "QuICC/Math/Constants.hpp"

#include <iostream>
namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Shell {

namespace RTC {

namespace Implicit {

   ModelBackend::ModelBackend()
      : IRTCBackend()
   {
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

   void ModelBackend::equationInfo(bool& isComplex, SpectralFieldIds& im, SpectralFieldIds& exL, SpectralFieldIds& exNL, SpectralFieldIds& exNS, int& indexMode, const SpectralFieldId& fId, const Resolution& res) const
   {
      // Operators are complex
      isComplex = true;

      // Implicit coupled fields
      im = this->implicitFields(fId);

      // Explicit linear terms
      exL = this->explicitLinearFields(fId);

      // Explicit nonlinear terms
      exNL = this->explicitNonlinearFields(fId);

      // Explicit nextstep terms
      exNS.clear();

      // Index mode
      indexMode = static_cast<int>(Equations::CouplingIndexType::SLOWEST_SINGLE_RHS);
   }

   void ModelBackend::blockSize(int& tN, int& gN, ArrayI& shift, int& rhs, const SpectralFieldId& fId, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs) const
   {
      auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, eigs.at(0))(0);
      tN = nN;

      int shiftR;
      if(this->useGalerkin())
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

   void ModelBackend::implicitBlock(DecoupledZSparse& decMat, const SpectralFieldId& rowId, const SpectralFieldId& colId, const int matIdx, const std::size_t bcType, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const
   {
      bool needStencil = (this->useGalerkin() && bcType == ModelOperatorBoundary::SolverNoTau::id());
      bool needTau = (bcType == ModelOperatorBoundary::SolverHasBc::id());

      assert(eigs.size() == 1);
      int m = eigs.at(0);

      // Compute system size
      const auto sysInfo = systemSize(rowId, colId, m, res);
      const auto sysN = std::get<0>(sysInfo)*std::get<1>(sysInfo);
      const auto baseRowShift = std::get<0>(sysInfo)*std::get<2>(sysInfo);
      const auto baseColShift = std::get<0>(sysInfo)*std::get<3>(sysInfo);
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
                  this->addBlock(decMat.real(), i2r2lapl.mat(), rowShift, colShift);
                  SparseSM::Chebyshev::LinearMap::I2Y2 i2r2(nN, nN, ri, ro);
                  this->addBlock(decMat.imag(), i2r2.mat(), rowShift, colShift, m*T*invlapl);
               }
               rowShift += nN;
               colShift += nN;
            }
         }
         else if(colId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::POL))
         {
            auto coriolis = [](const int l, const int m){
               return (l - 1.0)*(l + 1.0)*precision::sqrt(((l - m)*(l + m))/((2.0*l - 1.0)*(2.0*l + 1.0)));
            };

            const auto Ek = nds.find(NonDimensional::Ekman::id())->second->value();
            const auto T = 1.0/Ek;
            auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, m)(0);
            rowShift = baseRowShift + nN;
            colShift = baseColShift;
            for(int l = m+1; l < nL; l++)
            {
               auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);
               if(l > 0)
               {
                  const auto dl = static_cast<MHDFloat>(l);
                  const auto invlapl = 1.0/(dl*(dl + 1.0));
                  SparseSM::Chebyshev::LinearMap::I2Y1 cor_r(nN, nN, ri, ro);
                  auto norm = (dl - 1.0)*coriolis(l, m);
                  this->addBlock(decMat.real(), cor_r.mat(), rowShift, colShift, -norm*T*invlapl);

                  SparseSM::Chebyshev::LinearMap::I2Y2D1 cordr(nN, nN, ri, ro);
                  norm = -coriolis(l, m);
                  this->addBlock(decMat.real(), cordr.mat(), rowShift, colShift, -norm*T*invlapl);
               }
               rowShift += nN;
               colShift += nN;
            }
            nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, m)(0);
            rowShift = baseRowShift;
            colShift = baseColShift + nN;
            for(int l = m; l < nL - 1; l++)
            {
               auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);
               if(l > 0)
               {
                  const auto dl = static_cast<MHDFloat>(l);
                  const auto invlapl = 1.0/(dl*(dl + 1.0));
                  SparseSM::Chebyshev::LinearMap::I2Y1 cor_r(nN, nN, ri, ro);
                  auto norm = -(dl + 2.0)*coriolis(l+1, m);
                  this->addBlock(decMat.real(), cor_r.mat(), rowShift, colShift, -norm*T*invlapl);

                  SparseSM::Chebyshev::LinearMap::I2Y2D1 cordr(nN, nN, ri, ro);
                  norm = -coriolis(l+1, m);
                  this->addBlock(decMat.real(), cordr.mat(), rowShift, colShift, -norm*T*invlapl);
               }
               rowShift += nN;
               colShift += nN;
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
                  const auto dl = static_cast<MHDFloat>(l);
                  const auto invlapl = 1.0/(dl*(dl + 1.0));
                  SparseSM::Chebyshev::LinearMap::I4Y4SphLapl2 i4r4lapl2(nN, nN, ri, ro, l);
                  //SparseSM::Chebyshev::LinearMap::I2Y2SphLapl i2r2lapl(nN, nN, ri, ro, l);
                  this->addBlock(decMat.real(), i4r4lapl2.mat(), rowShift, colShift);
                  SparseSM::Chebyshev::LinearMap::I4Y4SphLapl i4r4lapl(nN, nN, ri, ro, l);
                  this->addBlock(decMat.imag(), i4r4lapl.mat(), rowShift, colShift, m*T*invlapl);
               }
               rowShift += nN;
               colShift += nN;
            }
         }
         else if(colId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::TOR))
         {
            auto coriolis = [](const int l, const int m){
               return (l - 1.0)*(l + 1.0)*precision::sqrt(((l - m)*(l + m))/((2.0*l - 1.0)*(2.0*l + 1.0)));
            };

            const auto Ek = nds.find(NonDimensional::Ekman::id())->second->value();
            const auto T = 1.0/Ek;
            auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, m)(0);
            rowShift = baseRowShift + nN;
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
                  this->addBlock(decMat.real(), cor_r.mat(), rowShift, colShift, norm*T*invlapl);
                  SparseSM::Chebyshev::LinearMap::I4Y4D1 cordr(nN, nN, ri, ro);
                  norm = -coriolis(l, m);
                  this->addBlock(decMat.real(), cordr.mat(), rowShift, colShift, norm*T*invlapl);
               }
               rowShift += nN;
               colShift += nN;
            }
            nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, m)(0);
            rowShift = baseRowShift;
            colShift = baseColShift + nN;
            for(int l = m; l < nL - 1; l++)
            {
               auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);
               if(l > 0)
               {
                  const auto dl = static_cast<MHDFloat>(l);
                  const auto invlapl = 1.0/(dl*(dl + 1.0));
                  SparseSM::Chebyshev::LinearMap::I4Y3 cor_r(nN, nN, ri, ro);
                  auto norm = -(dl + 2.0)*coriolis(l+1, m);
                  this->addBlock(decMat.real(), cor_r.mat(), rowShift, colShift, norm*T*invlapl);
                  SparseSM::Chebyshev::LinearMap::I4Y4D1 cordr(nN, nN, ri, ro);
                  norm = -coriolis(l+1, m);
                  this->addBlock(decMat.real(), cordr.mat(), rowShift, colShift, norm*T*invlapl);
               }
               rowShift += nN;
               colShift += nN;
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
                  this->addBlock(decMat.real(), i4r4.mat(), rowShift, colShift, -Ra);
               }
               rowShift += nN;
               colShift += nN;
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
               if(heatingMode == 0)
               {
                  SparseSM::Chebyshev::LinearMap::I2Y2SphLapl spasm(nN, nN, ri, ro, l);
                  this->addBlock(decMat.real(), spasm.mat(), rowShift, colShift, (1.0/Pr));
               }
               else
               {
                  SparseSM::Chebyshev::LinearMap::I2Y3SphLapl spasm(nN, nN, ri, ro, l);
                  this->addBlock(decMat.real(), spasm.mat(), rowShift, colShift, (1.0/Pr));
               }
               rowShift += nN;
               colShift += nN;
            }
         }
      }
      else
      {
         throw std::logic_error("Equations are not setup properly");
      }
   }

   void ModelBackend::timeBlock(DecoupledZSparse& decMat, const SpectralFieldId& fieldId, const int matIdx, const std::size_t bcType, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const
   {
      bool needStencil = (this->useGalerkin() && bcType == ModelOperatorBoundary::SolverNoTau::id());
      bool needTau = bcType == ModelOperatorBoundary::SolverHasBc::id();

      assert(eigs.size() == 1);
      int m = eigs.at(0);

      // Compute system size
      const auto sysInfo = systemSize(fieldId, fieldId, m, res);
      const auto sysN = std::get<0>(sysInfo)*std::get<1>(sysInfo);
      const auto baseShift = std::get<0>(sysInfo)*std::get<2>(sysInfo);
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

      if(fieldId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::TOR))
      {
         int shift = baseShift;
         for(int l = m; l < nL; l++)
         {
            auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);
            if(l > 0)
            {
               SparseSM::Chebyshev::LinearMap::I2Y2 spasm(nN, nN, ri, ro);
               this->addBlock(decMat.real(), spasm.mat(), shift, shift);
            }
            else
            {
               SparseSM::Chebyshev::LinearMap::Id qid(nN, nN, ri, ro);
               this->addBlock(decMat.real(), qid.mat(), shift, shift);
            }
            shift += nN;
         }
      }
      else if(fieldId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::POL))
      {
         int shift = baseShift;
         for(int l = m; l < nL; l++)
         {
            auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);
            if(l > 0)
            {
               SparseSM::Chebyshev::LinearMap::I4Y4SphLapl spasm(nN, nN, ri, ro, l);
               //SparseSM::Chebyshev::LinearMap::I2Y2 spasm(nN, nN, ri, ro);
               this->addBlock(decMat.real(), spasm.mat(), shift, shift);
            }
            else
            {
               SparseSM::Chebyshev::LinearMap::Id qid(nN, nN, ri, ro);
               this->addBlock(decMat.real(), qid.mat(), shift, shift);
            }
            shift += nN;
         }
      }
      else if(fieldId == std::make_pair(PhysicalNames::Temperature::id(), FieldComponents::Spectral::SCALAR))
      {
         auto heatingMode = nds.find(NonDimensional::Heating::id())->second->value();

         int shift = baseShift;
         for(int l = m; l < nL; l++)
         {
            auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);
            if(heatingMode == 0)
            {
               SparseSM::Chebyshev::LinearMap::I2Y2 spasm(nN, nN, ri, ro);
               this->addBlock(decMat.real(), spasm.mat(), shift, shift);
            }
            else
            {
               SparseSM::Chebyshev::LinearMap::I2Y3 spasm(nN, nN, ri, ro);
               this->addBlock(decMat.real(), spasm.mat(), shift, shift);
            }
            shift += nN;
         }
      }
   }

   void ModelBackend::boundaryBlock(DecoupledZSparse& decMat, const SpectralFieldId& rowId, const SpectralFieldId& colId, const int matIdx, const std::size_t bcType, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const
   {
      bool needStencil = this->useGalerkin();
      bool needTau = (bcType == ModelOperatorBoundary::SolverHasBc::id());

      assert(eigs.size() == 1);
      int m = eigs.at(0);

      // Compute system size
      const auto sysInfo = systemSize(rowId, colId, m, res);
      const auto sysN = std::get<0>(sysInfo)*std::get<1>(sysInfo);
      const auto baseRowShift = std::get<0>(sysInfo)*std::get<2>(sysInfo);
      const auto baseColShift = std::get<0>(sysInfo)*std::get<3>(sysInfo);
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

      int rowShift = baseRowShift;
      int colShift = baseColShift;
      // Apply boundary condition
      if(needStencil)
      {
         for(int l = m; l < nL; l++)
         {
            auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);
            SparseMatrix mat(nN, nN);
            this->applyGalerkinStencil(mat, rowId, colId, matIdx, res, eigs, bcs, nds);
            this->addBlock(decMat.real(), mat, rowShift, colShift);
            rowShift += nN;
            colShift += nN;
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
            this->applyTau(mat, rowId, colId, matIdx, res, eigs, bcs, nds);
            this->addBlock(decMat.real(), mat, rowShift, colShift);
            rowShift += nN;
            colShift += nN;
         }
      }
   }

   void ModelBackend::applyGalerkinStencil(SparseMatrix& mat, const SpectralFieldId& rowId, const SpectralFieldId& colId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const
   {
      assert(eigs.size() == 1);

      auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, eigs.at(0))(0);

      auto ri = nds.find(NonDimensional::Lower1d::id())->second->value();
      auto ro = nds.find(NonDimensional::Upper1d::id())->second->value();

      auto stencil = mat;
      this->galerkinStencil(stencil, colId, matIdx, res, eigs, false, bcs, nds);

      auto s = stencil.rows() - stencil.cols();
      SparseSM::Chebyshev::LinearMap::Id qId(nN-s, nN, ri, ro, 0, s);
      mat = qId.mat()*(mat*stencil);
   }

   void ModelBackend::applyTau(SparseMatrix& mat, const SpectralFieldId& rowId, const SpectralFieldId& colId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const
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

      mat += bcOp.mat();
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
      else if(opId == ModelOperator::ImplicitLinear::id())
      {
         for(auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
         {
            for(auto pColId = imRange.first; pColId != imRange.second; pColId++)
            {
               this->implicitBlock(rModelMatrix, *pRowId, *pColId, matIdx, bcType, res, eigs, bcs, nds);
            }
         }
      }
      // Boundary operator
      else if(opId == ModelOperator::Boundary::id())
      {
         for(auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
         {
            for(auto pColId = imRange.first; pColId != imRange.second; pColId++)
            {
               this->boundaryBlock(rModelMatrix, *pRowId, *pColId, matIdx, bcType, res, eigs, bcs, nds);
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
      int m = eigs.at(0);

      // Compute system size
      const auto sysInfo = systemSize(rowId, colId, m, res);
      const auto sysN = std::get<0>(sysInfo);
      const auto baseShift = 0;
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

      // Explicit linear operator
      if(opId == ModelOperator::ExplicitLinear::id())
      {
         if(rowId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::POL) &&
               colId == std::make_pair(PhysicalNames::Temperature::id(),FieldComponents::Spectral::SCALAR))
         {
            auto Ra = effectiveRa(nds);

            int shift = baseShift;
            for(int l = m; l < nL; l++)
            {
               auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);
               SparseSM::Chebyshev::LinearMap::I4Y4 spasm(nN, nN, ri, ro);
               this->addBlock(decMat.real(), spasm.mat(), shift, shift, Ra);
               shift += nN;
            }
         }
         else if(rowId == std::make_pair(PhysicalNames::Temperature::id(), FieldComponents::Spectral::SCALAR) && 
               colId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::POL))
         {
            auto heatingMode = nds.find(NonDimensional::Heating::id())->second->value();
            auto bg = effectiveBg(nds);

            int shift = baseShift;
            for(int l = m; l < nL; l++)
            {
               auto dl = static_cast<MHDFloat>(l);
               auto ll1 = dl*(dl + 1.0);
               auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);
               if(heatingMode == 0)
               {
                  SparseSM::Chebyshev::LinearMap::I2Y2 spasm(nN, nN, ri, ro);
                  this->addBlock(decMat.real(), spasm.mat(), shift, shift, -bg*ll1);
               }
               else
               {
                  SparseSM::Chebyshev::LinearMap::I2 spasm(nN, nN, ri, ro);
                  this->addBlock(decMat.real(), spasm.mat(), shift, shift, -bg*ll1);
               }
               shift += nN;
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

            int shift = baseShift;
            for(int l = m; l < nL; l++)
            {
               auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);
               if(heatingMode == 0)
               {
                  SparseSM::Chebyshev::LinearMap::I2Y2 spasm(nN, nN, ri, ro);
                  this->addBlock(decMat.real(), spasm.mat(), shift, shift);
               }
               else
               {
                  SparseSM::Chebyshev::LinearMap::I2Y3 spasm(nN, nN, ri, ro);
                  this->addBlock(decMat.real(), spasm.mat(), shift, shift);
               }
               shift += nN;
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

   std::tuple<int, int, int, int> ModelBackend::systemSize(const SpectralFieldId& rowId, const SpectralFieldId& colId, const int m, const Resolution& res) const
   {
      auto sysN = res.counter().dimensions(Dimensions::Space::SPECTRAL, m)(0);
      auto nL = res.counter().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL, m);

      for(int l = m+1; l < nL; l++)
      {
         auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);
         sysN += nN;
      }

      int rowIdx = 0;
      int colIdx = 0;
      int idx = 0;
      const auto& fields = this->implicitFields(rowId);
      for(auto it = fields.begin(); it != fields.end(); ++it)
      {
         // Get block index of rowId
         if(rowId == *it)
         {
            rowIdx = idx;
         }

         // Get block index of colId
         if(colId == *it)
         {
            colIdx = idx;
         }
         idx++;
      }
      int nFields = fields.size();

      return std::make_tuple(sysN, nFields, rowIdx, colIdx);
   }


} // Implicit
} // RTC
} // Shell
} // Boussinesq
} // Model
} // QuICC
