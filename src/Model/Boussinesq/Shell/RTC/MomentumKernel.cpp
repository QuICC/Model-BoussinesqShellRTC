/**
 * @file MomentumKernel.cpp
 * @brief Source of physical space kernel for the Momentum equation
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Model/Boussinesq/Shell/RTC/MomentumKernel.hpp"

// Project includes
//
#include "QuICC/SpatialScheme/3D/SLFl.hpp"
#include "QuICC/PhysicalOperators/Cross.hpp"
#include "QuICC/PhysicalOperators/SphericalCoriolis.hpp"

namespace QuICC {

namespace Physical {

namespace Kernel {

   MomentumKernel::MomentumKernel()
      : IPhysicalKernel()
   {
   }

   MomentumKernel::~MomentumKernel()
   {
   }

   std::size_t MomentumKernel::name() const
   {
      return this->mName;
   }

   void MomentumKernel::setVelocity(std::size_t name, Framework::Selector::VariantSharedVectorVariable spField)
   {
      // Safety assertion
      assert(this->mScalars.count(name) + this->mVectors.count(name) == 0);

      this->mName = name;

      this->setField(name, spField);
   }

   void MomentumKernel::init(const MHDFloat inertia, const MHDFloat coriolis)
   {
      // Set scaling constants
      this->mInertia = inertia;
      this->mCoriolis = coriolis;
   }

   void MomentumKernel::setMesh(std::shared_ptr<std::vector<Array> > spMesh)
   {
      IPhysicalKernel::setMesh(spMesh);

      if(std::visit([&](auto&& v)->bool{return (v->dom(0).res().sim().ss().has(SpatialScheme::Feature::SpectralOrdering132));},this->vector(this->name())))
      {
            this->mCosTheta = this->mspMesh->at(1).array().cos();
            this->mSinTheta = this->mspMesh->at(1).array().sin();
      }
   }

   void MomentumKernel::compute(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id id) const
   {
      ///
      /// Compute \f$\left(\nabla\wedge\vec u\right)\wedge\vec u\f$
      ///
      switch(id)
      {
         case(FieldComponents::Physical::R):
            std::visit([&](auto&& v){Physical::Cross<FieldComponents::Physical::THETA,FieldComponents::Physical::PHI>::set(rNLComp, v->dom(0).curl(), v->dom(0).phys(), this->mInertia);}, this->vector(this->name()));
            break;
         case(FieldComponents::Physical::THETA):
            std::visit([&](auto&& v){Physical::Cross<FieldComponents::Physical::PHI,FieldComponents::Physical::R>::set(rNLComp, v->dom(0).curl(), v->dom(0).phys(), this->mInertia);}, this->vector(this->name()));
            break;
         case(FieldComponents::Physical::PHI):
            std::visit([&](auto&& v){Physical::Cross<FieldComponents::Physical::R,FieldComponents::Physical::THETA>::set(rNLComp, v->dom(0).curl(), v->dom(0).phys(), this->mInertia);}, this->vector(this->name()));
            break;
         default:
            assert(false);
            break;
      }

      if(std::visit([&](auto&& v)->bool{return (v->dom(0).res().sim().ss().has(SpatialScheme::Feature::SpectralOrdering132));},this->vector(this->name())))
      {
         ///
         /// Compute Coriolis term
         ///
         std::visit([&](auto&& v){Physical::SphericalCoriolis::add(rNLComp, id, v->dom(0).res(), this->mCosTheta, this->mSinTheta, v->dom(0).phys(), this->mCoriolis);}, this->vector(this->name()));
      }
   }

}
}
}
