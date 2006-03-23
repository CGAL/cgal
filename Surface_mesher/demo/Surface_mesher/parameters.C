// Copyright (c) 2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://lrineau@scm.gforge.inria.fr/svn/cgal/trunk/Surface_mesher/examples/Surface_mesher/implicit_functions.C $
// $Id: implicit_functions.C 29116 2006-03-07 10:01:08Z lrineau $
//
//
// Author(s)     : Laurent Rineau

#include <assert.h>

#include "parameters.h"
#include "implicit_functions.h"

String_options string_options;
Double_options double_options;
Implicit_function_map functions;

void init_string_options()
{
  string_options["surface_mesh"] = "";
  string_options["surface_off"] = "";
  string_options["read_initial_points"] = "";
  string_options["inrimage"] = "";
}

void init_double_options()
{
  double_options["distance_bound"] = 1000000000000000.;
  
  // bound on radius of surface Delaunay balls
  double_options["radius_bound"] = 10000000000000.;

  double_options["angle_bound"] = 0.; // in degrees
  double_options["enclosing_sphere_radius"] = 6.;
  double_options["precision"] = 1e-3;
  double_options["center_x"] = 0.0;
  double_options["center_y"] = 0.0;
  double_options["center_z"] = 0.0;
  double_options["number_of_initial_points"] = 20;
  double_options["iso_value"] = 0.;
}

void init_functions()
{
  functions["generic_inrimage"] = &generic_inrimage_function;
  functions["chair"] = &chair_function;
  functions["ellipsoid"] = &ellipsoid_function;
  functions["false_knot"] = &false_knot_function;
  functions["heart"] = &heart_function;
  functions["klein"] = &klein_function;
  functions["knot1"] = &knot1_function;
  functions["knot2"] = &knot2_function;
  functions["knot3"] = &knot3_function;
  functions["octic"] = &octic_function;
  functions["ring"] = &ring_function;
  functions["sphere"] = &sphere_function;
  functions["tanglecube"] = &tanglecube_function;
  functions["torus"] = &torus_function;
}

void init_parameters()
{
  init_string_options();
  init_double_options();
  init_functions();
}

Gray_image* isosurface = 0;

double generic_inrimage_function(double x, double y, double z)
{
  assert(isosurface != 0);

  return (*isosurface)(x, y, z);
}
