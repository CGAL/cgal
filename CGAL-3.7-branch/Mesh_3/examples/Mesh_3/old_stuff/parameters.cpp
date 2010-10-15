// Copyright (c) 2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Laurent Rineau

#include <assert.h>
#include <iostream>

#include "parameters.h"
#include "implicit_functions.h"

String_options string_options;
Double_options double_options;
Implicit_function_map functions;
Check_strings check_strings;

Usage* usage_ptr = 0;

/** debugging function, auxiliary function for init_check_strings() */
template <class Map>
void init_check_strings_aux(const Map& map)
{
  for(typename Map::const_iterator it = map.begin();
      it != map.end();
      ++it)
  {
    Check_strings::const_iterator chk_it = check_strings.find(it->first);
    if( chk_it != check_strings.end() )
      std::cerr << "Warning: duplicate option --"  << it->first << "\n";
    else
      check_strings.insert(it->first); // TODO: can be improved with lower_bound()
  }
}

/** Debugging function.
    Check that no option is used twice.
*/
void init_check_strings()
{
  init_check_strings_aux(string_options);
  init_check_strings_aux(double_options);
  init_check_strings_aux(functions);
}

void check_all_options_have_been_used()
{
  for(Check_strings::const_iterator it = check_strings.begin();
      it!= check_strings.end();
      ++it)
    if(functions.find(*it) == functions.end())
      std::cerr << "Warning: option --" << *it << " was not used.\n";
}

void init_string_options()
{
  string_options["dump_of_initial_points"] = "";
  string_options["read_initial_points"] = "";
  string_options["inrimage"] = "";

  string_options["initial_surface_off"] = "";

  string_options["cgal_mesh_before_exudation"] = "";

  string_options["mesh_after_filename"] = "";
  string_options["cgal_mesh_after_filename"] = "";

  string_options["slivers_off"] = "";

  string_options["surface_mesh"] = "";
  string_options["surface_off"] = "";
}

void init_double_options()
{
  double_options["r"] = 0;
  double_options["special"] = 0.;
  double_options["distance_bound"] = 1000000000000000.;

  // bound on radius of surface Delaunay balls
  double_options["radius_bound"] = 10000000000000.;

  // bound on the radius of circumspheres
  double_options["tets_radius_bound"] = 100000000000000.;
  double_options["angle_bound"] = 30.; // in degrees
  double_options["tets_aspect_ratio_bound"] = 0.;
  double_options["enclosing_sphere_radius"] = 6.;
  double_options["precision"] = 1e-6;
  double_options["center_x"] = 0.0;
  double_options["center_y"] = 0.0;
  double_options["center_z"] = 0.0;
  double_options["number_of_initial_points"] = 20;
  double_options["pumping_bound"] = 0;
  double_options["sliver_test"] = 0.15;
  double_options["iso_value"] = 0.;

  double_options["edges_radius_bound"] = 0.;
  double_options["edges_distance_bound"] = 0.;
  double_options["sharp_edge_cosine_bound"] = 0.5;
}

void init_functions()
{
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
  functions["cube"] = &cube_function;
}

void init_parameters()
{
  init_string_options();
  init_double_options();
  init_functions();

  init_check_strings(); // for debugging
}

template <class Map>
typename Map::mapped_type
get_option(const Map& map, const std::string& s)
{
  typename Map::const_iterator it = map.find(s);
  if( it == map.end() )
  {
    std::cerr << "Error: invalid option \"--" << s << "\"\n";
    if( usage_ptr ) usage_ptr("");
    return typename Map::mapped_type();
  }
  else
  {
    check_strings.erase(s);
    return it->second;
  }
}

std::string get_string_option(const std::string s)
{
  return get_option(string_options, s);
}

double get_double_option(std::string s)
{
  return get_option(double_options, s);
}

Implicit_function* get_function(std::string s)
{
  return get_option(functions, s);
}
