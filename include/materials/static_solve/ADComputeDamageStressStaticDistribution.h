//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADComputeStressBase.h"

/**
 * ADComputeDamageStressStaticDistribution computes the stress following linear elasticity theory (small strains)
 */
class ADComputeDamageStressStaticDistribution : public ADComputeStressBase
{
public:
  static InputParameters validParams();

  ADComputeDamageStressStaticDistribution(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  virtual void computeQpStress() override;

  /// @brief Compute gamma_r
  /// @return gamma_r
  ADReal computegammar();

  /// Material property initial damage profile

  /// initial lambda value 
  Real _lambda_o;

  /// initial shear modulus value
  Real _shear_modulus_o;

  /// xi_o value
  Real _xi_o;

  /// @brief initial damage value
  const ADMaterialProperty<Real> & _initial_damage_val;

};
