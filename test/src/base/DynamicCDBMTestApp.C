//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "DynamicCDBMTestApp.h"
#include "DynamicCDBMApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
DynamicCDBMTestApp::validParams()
{
  InputParameters params = DynamicCDBMApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

DynamicCDBMTestApp::DynamicCDBMTestApp(InputParameters parameters) : MooseApp(parameters)
{
  DynamicCDBMTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

DynamicCDBMTestApp::~DynamicCDBMTestApp() {}

void
DynamicCDBMTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  DynamicCDBMApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"DynamicCDBMTestApp"});
    Registry::registerActionsTo(af, {"DynamicCDBMTestApp"});
  }
}

void
DynamicCDBMTestApp::registerApps()
{
  registerApp(DynamicCDBMApp);
  registerApp(DynamicCDBMTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
DynamicCDBMTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  DynamicCDBMTestApp::registerAll(f, af, s);
}
extern "C" void
DynamicCDBMTestApp__registerApps()
{
  DynamicCDBMTestApp::registerApps();
}
