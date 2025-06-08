/*
Define Function for Initial Shear Stress for benchmark
*/

#include "InitialShearStressCDBM.h"

#include <string.h>

registerMooseObject("DynamicCDBMApp", InitialShearStressCDBM);

InputParameters
InitialShearStressCDBM::validParams()
{
  InputParameters params = Function::validParams();
  params.addRequiredParam<Real>("peak_value", "Peak value of the initial shear stress");
  params.addRequiredParam<Real>("domain_value", "Domain value of the initial shear stress");
  params.addRequiredParam<Real>("nucl_center_x", "X coordinate of the nucleation center");
  params.addRequiredParam<Real>("nucl_center_y", "Y coordinate of the nucleation center");
  params.addRequiredParam<Real>("nucl_size", "Size of the nucleation zone");
  params.addRequiredParam<Real>("elem_size", "Element size for the simulation");
  return params;
}

InitialShearStressCDBM::InitialShearStressCDBM(const InputParameters & parameters)
  : Function(parameters),
  _peak_value(getParam<Real>("peak_value")),
  _domain_value(getParam<Real>("domain_value")),
  _nucl_center_x(getParam<Real>("nucl_center_x")),
  _nucl_center_y(getParam<Real>("nucl_center_y")),
  _nucl_size(getParam<Real>("nucl_size")),
  _elem_size(getParam<Real>("elem_size"))
{
}

Real
InitialShearStressCDBM::value(Real /*t*/, const Point & p) const
{

  Real x_coord = p(0); //along the strike direction
  Real y_coord = p(1); //along the dip direction
  Real z_coord = p(2); //along the normal direction

  Real T1_o = 0.0;
    
  //tpv205
  if ((x_coord<=(_nucl_center_x+_nucl_size*0.5))&&(x_coord>=(_nucl_center_x-_nucl_size*0.5))&& (y_coord<=(_nucl_center_y+_nucl_size*0.5))&&(y_coord>=(_nucl_center_y-_nucl_size*0.5))&&(z_coord>=-1*_elem_size)&&(z_coord<=_elem_size))
  {
      T1_o = _peak_value;
  }
  else
  {
      T1_o = _domain_value;
  }

  return T1_o;

}