//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADComputeDamageStressStaticDistribution.h"

registerMooseObject("DynamicCDBMApp", ADComputeDamageStressStaticDistribution);

InputParameters
ADComputeDamageStressStaticDistribution::validParams()
{
  InputParameters params = ADComputeStressBase::validParams();
  params.addClassDescription("Compute stress using elasticity for small strains");
  params.addRequiredParam<Real>("lambda_o","initial lambda value");
  params.addRequiredParam<Real>("shear_modulus_o","initial shear modulus value");
  params.addRequiredParam<Real>("xi_o","xi_o value");
  return params;
}

ADComputeDamageStressStaticDistribution::ADComputeDamageStressStaticDistribution(const InputParameters & parameters)
  : ADComputeStressBase(parameters),
  _lambda_o(getParam<Real>("lambda_o")),
  _shear_modulus_o(getParam<Real>("shear_modulus_o")),
  _xi_o(getParam<Real>("xi_o")),
  _initial_damage_val(getADMaterialPropertyByName<Real>("initial_damage"))   
{
}

void
ADComputeDamageStressStaticDistribution::initialSetup()
{
  // _base_name + "unstabilized_deformation_gradient" is only declared if we're
  // using the Lagrangian kernels.  It's okay to invoke this small strain
  // material if you are using that kernel system and the
  // ComputeLagrangianWrappedStress wrapper
  if (hasBlockMaterialProperty<RankTwoTensor>(_base_name + "strain_increment") &&
      !hasBlockMaterialProperty<RankTwoTensor>(_base_name + "unstabilized_deformation_gradient"))
    mooseError("This linear elastic stress calculation only works for small strains; use "
               "ComputeFiniteStrainElasticStress for simulations using incremental and finite "
               "strains.");
}

void
ADComputeDamageStressStaticDistribution::computeQpStress()
{
  
  // Compute gammar
  ADReal gamma_r = computegammar();

  // Evaluate shear modulus
  ADReal shear_modulus = _shear_modulus_o + _xi_o * _initial_damage_val[_qp] * gamma_r;
  ADReal gamma_damaged_out = _initial_damage_val[_qp] * gamma_r;

  //
  const Real epsilon = 1e-12;
  ADReal I1 = epsilon + _mechanical_strain[_qp](0,0) + _mechanical_strain[_qp](1,1) + _mechanical_strain[_qp](2,2);
  ADReal I2 = epsilon + _mechanical_strain[_qp](0,0) * _mechanical_strain[_qp](0,0) + _mechanical_strain[_qp](1,1) * _mechanical_strain[_qp](1,1) + _mechanical_strain[_qp](2,2) * _mechanical_strain[_qp](2,2) + 2 * _mechanical_strain[_qp](1,2) * _mechanical_strain[_qp](1,2) + 2 * _mechanical_strain[_qp](0,1) * _mechanical_strain[_qp](0,1) + 2 * _mechanical_strain[_qp](0,2) * _mechanical_strain[_qp](0,2);
  ADReal xi = I1 / std::sqrt(I2);

  // stress = C * e
  _stress[_qp](0,0) = ( _lambda_o - gamma_damaged_out / xi ) * I1 + ( 2 * shear_modulus - gamma_damaged_out * xi ) * _mechanical_strain[_qp](0,0);
  _stress[_qp](1,1) = ( _lambda_o - gamma_damaged_out / xi ) * I1 + ( 2 * shear_modulus - gamma_damaged_out * xi ) * _mechanical_strain[_qp](1,1);
  _stress[_qp](2,2) = ( _lambda_o - gamma_damaged_out / xi ) * I1 + ( 2 * shear_modulus - gamma_damaged_out * xi ) * _mechanical_strain[_qp](2,2);
  _stress[_qp](0,1) = ( 2 * shear_modulus - gamma_damaged_out * xi ) * _mechanical_strain[_qp](0,1); 
  _stress[_qp](1,0) = ( 2 * shear_modulus - gamma_damaged_out * xi ) * _mechanical_strain[_qp](1,0);
  _stress[_qp](0,2) = ( 2 * shear_modulus - gamma_damaged_out * xi ) * _mechanical_strain[_qp](0,2); 
  _stress[_qp](2,0) = ( 2 * shear_modulus - gamma_damaged_out * xi ) * _mechanical_strain[_qp](2,0);
  _stress[_qp](1,2) = ( 2 * shear_modulus - gamma_damaged_out * xi ) * _mechanical_strain[_qp](1,2); 
  _stress[_qp](2,1) = ( 2 * shear_modulus - gamma_damaged_out * xi ) * _mechanical_strain[_qp](2,1);

  // Assign value for elastic strain, which is equal to the mechanical strain
  _elastic_strain[_qp] = _mechanical_strain[_qp];

}

ADReal
ADComputeDamageStressStaticDistribution::computegammar()
{
  // Calculate each part of the expression
  ADReal term1 = -_xi_o * (-_lambda_o * pow(_xi_o, 2) + 6 * _lambda_o + 2 * _shear_modulus_o);
  ADReal term2_sqrt = sqrt((_lambda_o * pow(_xi_o, 2) + 2 * _shear_modulus_o) * 
                            (_lambda_o * pow(_xi_o, 4) - 12 * _lambda_o * pow(_xi_o, 2) + 36 * _lambda_o 
                            - 6 * _shear_modulus_o * pow(_xi_o, 2) + 24 * _shear_modulus_o));
  ADReal denominator = 2 * (pow(_xi_o, 2) - 3);
  
  // Calculate gamma_r
  ADReal gamma_r = (term1 - term2_sqrt) / denominator;
  
  return gamma_r;
}