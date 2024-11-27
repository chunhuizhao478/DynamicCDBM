//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PerturbationRadial.h"

/**
 *  Created by Chunhui Zhao, Nov 26th, 2024
 *  Material used in Create Time Dependent Damage/Shear Stress Perturbation in the Dynamic Solve
 */
registerMooseObject("DynamicCDBMApp", PerturbationRadial);

InputParameters
PerturbationRadial::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Material used in Create Time Dependent Damage Perturbation in the Dynamic Solve");
  params.addRequiredParam<std::vector<Real>>("nucl_center", "nucleation center (x,y,z)");
  params.addRequiredParam<Real>("peak_value","the peak damage value apply region normal to the fault (exponential decay)");
  params.addRequiredParam<Real>("thickness","the standard deviation used in apply damage value normal to the fault (exponential decay)");
  params.addRequiredParam<Real>("length","the standard deviation used in apply damage value normal to the fault (exponential decay)");
  params.addRequiredParam<Real>("duration","duration to reach peak damage");
  params.addRequiredParam<std::string>("perturbation_type", "Type of perturbation: 'damage' or 'shear_stress'"); // New parameter
  params.addParam<Real>("sigma_divisor", 2.0, "sigma value = (length / sigma_divisor)");
  return params;
}

PerturbationRadial::PerturbationRadial(const InputParameters & parameters)
  : Material(parameters),
  _damage_perturbation(declareProperty<Real>("damage_perturbation")),
  _damage_perturbation_old(getMaterialPropertyOldByName<Real>("damage_perturbation")),
  _shear_stress_perturbation(declareProperty<Real>("shear_stress_perturbation")),
  _shear_stress_perturbation_old(getMaterialPropertyOldByName<Real>("shear_stress_perturbation")),
  _nucl_center(getParam<std::vector<Real>>("nucl_center")),
  _peak_value(getParam<Real>("peak_value")),
  _thickness(getParam<Real>("thickness")),
  _length(getParam<Real>("length")),
  _duration(getParam<Real>("duration")),
  _perturbation_type(getParam<std::string>("perturbation_type")), // Initialize new parameter
  _sigma_divisor(getParam<Real>("sigma_divisor"))
{
  //in case I'm stupid
  if (_nucl_center.size() != 3){
    mooseError("fault plane parameter must accepts 3 numbers!");
  }
}

void
PerturbationRadial::initQpStatefulProperties()
{
  _damage_perturbation[_qp] = 0.0;
  _shear_stress_perturbation[_qp] = 0.0;
}

void
PerturbationRadial::computeQpProperties()
{

  //Get coordinates
  //here no rotation is applied yet
  Real xcoord = _q_point[_qp](0); //strike
  Real ycoord = _q_point[_qp](1); //dip
  Real zcoord = _q_point[_qp](2); //normal

  Real sigma_x = _length / _sigma_divisor;
  Real sigma_y = _length / _sigma_divisor;
  Real gaussian_factor = _peak_value; // Maximum perturbation

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
    scaled_gaussian_value = gaussian_value / (_duration / _dt);
  }
  else
  {
    scaled_gaussian_value = 0;
  }

  // Check Z direction constraint and apply perturbation
  Real dalpha_damage = 0.0;
  Real dalpha_stress = 0.0;
  if (zcoord >= _nucl_center[2] - _thickness / 2.0 && zcoord <= _nucl_center[2] + _thickness / 2.0)
  {
    dalpha_damage = _damage_perturbation_old[_qp] + scaled_gaussian_value;
    dalpha_stress = _shear_stress_perturbation_old[_qp] + scaled_gaussian_value;
  }
  else
  {
    dalpha_damage = _damage_perturbation_old[_qp];
    dalpha_stress = _shear_stress_perturbation_old[_qp];
  }

  if (_perturbation_type == "damage")
  {
    _damage_perturbation[_qp] = dalpha_damage;
    _shear_stress_perturbation[_qp] = 0.0;
  }
  else if (_perturbation_type == "shear_stress")
  {
    _damage_perturbation[_qp] = 0.0;
    _shear_stress_perturbation[_qp] = dalpha_stress;
  }
  else
  {
    mooseError("Invalid perturbation type: " + _perturbation_type);
  }

}