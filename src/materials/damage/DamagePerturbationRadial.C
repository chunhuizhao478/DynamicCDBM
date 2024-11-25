//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DamagePerturbationRadial.h"

/**
 *  Created by Chunhui Zhao, Aug 27th, 2024
 *  Material used in Create Time Dependent Damage Perturbation in the Dynamic Solve
 */
registerMooseObject("DynamicCDBMApp", DamagePerturbationRadial);

InputParameters
DamagePerturbationRadial::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Material used in Create Time Dependent Damage Perturbation in the Dynamic Solve");
  params.addRequiredParam<std::vector<Real>>("nucl_center", "nucleation center (x,y,z)");
  params.addRequiredParam<Real>("e_damage","the peak damage value apply region normal to the fault (exponential decay)");
  params.addRequiredParam<Real>("thickness","the standard deviation used in apply damage value normal to the fault (exponential decay)");
  params.addRequiredParam<Real>("length","the standard deviation used in apply damage value normal to the fault (exponential decay)");
  params.addRequiredParam<Real>("duration","duration to reach peak damage");
  return params;
}

DamagePerturbationRadial::DamagePerturbationRadial(const InputParameters & parameters)
  : Material(parameters),
  _damage_perturbation(declareProperty<Real>("damage_perturbation")),
  _damage_perturbation_old(getMaterialPropertyOldByName<Real>("damage_perturbation")),
  _nucl_center(getParam<std::vector<Real>>("nucl_center")),
  _peak_damage(getParam<Real>("e_damage")),
  _thickness(getParam<Real>("thickness")),
  _length(getParam<Real>("length")),
  _duration(getParam<Real>("duration"))
{
  //in case I'm stupid
  if (_nucl_center.size() != 3){
    mooseError("fault plane parameter must accepts 3 numbers!");
  }
}

void
DamagePerturbationRadial::initQpStatefulProperties()
{
  _damage_perturbation[_qp] = 0.0;
}

void
DamagePerturbationRadial::computeQpProperties()
{

  //Get coordinates
  //here no rotation is applied yet
  Real xcoord = _q_point[_qp](0); //strike
  Real ycoord = _q_point[_qp](1); //dip
  Real zcoord = _q_point[_qp](2); //normal

  Real sigma_x = _length / 2.0;
  Real sigma_y = _length / 2.0;
  Real gaussian_factor = _peak_damage; // Maximum perturbation

  // Compute distance from nucleation center (XY plane only)
  Real dx = xcoord - _nucl_center[0];
  Real dy = ycoord - _nucl_center[1];

  // 2D Gaussian distribution in XY plane
  Real gaussian_value = gaussian_factor * std::exp(
    - (dx * dx) / (2.0 * sigma_x * sigma_x)
    - (dy * dy) / (2.0 * sigma_y * sigma_y)
  );

  // Scale Gaussian value over time
  Real scaled_gaussian_value = 0.0;
  if (_t <= _duration)
  {
    scaled_gaussian_value = gaussian_value * (_t / _duration);
  }
  else
  {
    scaled_gaussian_value = gaussian_value;
  }

  // Check Z direction constraint and apply perturbation
  Real dalpha = 0.0;
  if (zcoord >= _nucl_center[2] - _thickness / 2.0 && zcoord <= _nucl_center[2] + _thickness / 2.0)
  {
    dalpha = _damage_perturbation_old[_qp] + scaled_gaussian_value;
  }
  else
  {
    dalpha = _damage_perturbation_old[_qp];
  }

  _damage_perturbation[_qp] = dalpha;
}