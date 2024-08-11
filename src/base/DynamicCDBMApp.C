#include "DynamicCDBMApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
DynamicCDBMApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

DynamicCDBMApp::DynamicCDBMApp(InputParameters parameters) : MooseApp(parameters)
{
  DynamicCDBMApp::registerAll(_factory, _action_factory, _syntax);
}

DynamicCDBMApp::~DynamicCDBMApp() {}

void
DynamicCDBMApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAllObjects<DynamicCDBMApp>(f, af, s);
  Registry::registerObjectsTo(f, {"DynamicCDBMApp"});
  Registry::registerActionsTo(af, {"DynamicCDBMApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
DynamicCDBMApp::registerApps()
{
  registerApp(DynamicCDBMApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
DynamicCDBMApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  DynamicCDBMApp::registerAll(f, af, s);
}
extern "C" void
DynamicCDBMApp__registerApps()
{
  DynamicCDBMApp::registerApps();
}
