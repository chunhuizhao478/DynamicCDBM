//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"

/**
 *  Created by Chunhui Zhao, Jun 9th, 2025
 *  Initial damage profile near the fault plane 
 */
class NearFaultInitialDamage : public Material
{
public:
  static InputParameters validParams();

  NearFaultInitialDamage(const InputParameters & parameters);

  virtual void computeQpProperties() override;

  virtual void initQpStatefulProperties() override;

protected:

  /// Material property initial damage profile
  MaterialProperty<Real> & _initial_damage;

  /// Material property old initial damage profile
  const MaterialProperty<Real> & _initial_damage_old;

  /// Distance from the fault plane along the normal direction
  const Real _distance_from_fault;

  /// Initial damage value within the fault zone
  const Real _initial_damage_value;

};