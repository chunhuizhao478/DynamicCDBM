//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NearFaultInitialDamage.h"

/**
 *  Created by Chunhui Zhao, Nov 26th, 2024
 *  Material used in Create Time Dependent Damage/Shear Stress Perturbation in the Dynamic Solve
 */
registerMooseObject("DynamicCDBMApp", NearFaultInitialDamage);

InputParameters
NearFaultInitialDamage::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Material used in Create Time Dependent Damage Perturbation in the Dynamic Solve");
  params.addRequiredParam<Real>("distance_from_fault","distance from the fault plane along the normal direction");
  params.addRequiredParam<Real>("initial_damage_value","initial damage value within the fault zone");
  return params;
}

NearFaultInitialDamage::NearFaultInitialDamage(const InputParameters & parameters)
  : Material(parameters),
  _initial_damage(declareProperty<Real>("initial_damage")),
  _initial_damage_old(getMaterialPropertyOldByName<Real>("initial_damage")),
  _distance_from_fault(getParam<Real>("distance_from_fault")),
  _initial_damage_value(getParam<Real>("initial_damage_value"))
{
}

void
NearFaultInitialDamage::initQpStatefulProperties()
{
  _initial_damage[_qp] = 0.0;
}

void
NearFaultInitialDamage::computeQpProperties()
{

  // Get the current point coordinates in the mesh
  // const Real xcoord = _q_point[_qp](0); // strike direction
  // const Real ycoord = _q_point[_qp](1); // dip direction
  const Real zcoord = _q_point[_qp](2); // normal direction

  // Set the initial damage within the fault zone
  if (zcoord <= _distance_from_fault*0.5 && zcoord >= -_distance_from_fault*0.5)
    _initial_damage[_qp] = _initial_damage_value; // Set to 1.0 within the fault zone
  else
    _initial_damage[_qp] = 0.0; // Set to 0.0 outside the fault zone

}