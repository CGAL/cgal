// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Stéphane Tayeb
//
//******************************************************************************
// File Description :
// Outputs to out.mesh a mesh of implicit domains. These domains are defined
// by a vector of functions. Each n-uplet of sign of function values defines a
// subdomain.
//******************************************************************************



#include <debug.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Mesh_3/implicit_to_labeled_function_wrapper.h>
#include <CGAL/Mesh_3/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include "../examples/Mesh_3/implicit_functions.h"

// IO
#include <CGAL/IO/File_medit.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;


// Domain
struct K: public CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef FT_to_point_function_wrapper<K::FT, K::Point_3> Function;
typedef CGAL::Mesh_3::Implicit_vector_to_labeled_function_wrapper<Function, K>
                                                        Function_wrapper;
typedef Function_wrapper::Function_vector Function_vector;
typedef CGAL::Mesh_3::Labeled_mesh_domain_3<Function_wrapper, K> Mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Mesh Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
typedef Mesh_criteria::Facet_criteria    Facet_criteria;
typedef Mesh_criteria::Cell_criteria     Cell_criteria;


template <typename T>
double set_arg(const std::string& param_name,
               const std::string& param_string,
               const po::variables_map& vm)
{
  T param_value(0);
  
  if ( vm.count(param_name) )
  {
    param_value = vm[param_name].as<T>();
    std::cout << param_string << ": " << param_value << "\n";
  }
  else
  {
    std::cout << param_string << " ignored.\n";
  }
  
  return param_value;
}


void set_function(Function_vector& v,
                  Function& f,
                  const std::string& function_name,
                  const po::variables_map& vm)
{
  if ( vm.count(function_name) )
  {
    v.push_back(&f);
    std::cout << function_name << " ";
  }
}


int main(int argc, char* argv[])
{
  po::options_description generic("Generic options");
  generic.add_options() ("help", "Produce help message");
  
  po::options_description mesh("Mesh generation parameters");
  mesh.add_options()("facet_angle", po::value<double>(), "Set facet angle bound")
  ("facet_size", po::value<double>(), "Set facet size bound")
  ("facet_error", po::value<double>(), "Set facet approximation error bound")
  ("tet_shape", po::value<double>(), "Set tet radius-edge bound")
  ("tet_size", po::value<double>(), "Set tet size bound");
  
  po::options_description functions("Implicit functions");
  functions.add_options()("torus", "Mesh torus function")
  ("sphere", "Mesh sphere function")
  ("chair", "Mesh chair function")
  ("tanglecube", "Mesh tanglecube function")
  ("cube", "Mesh cube function")
  ("ellipsoid", "Mesh ellipsoid function")
  ("heart", "Mesh heart function")
  ("octic", "Mesh octic function");
  
  po::options_description desc("Options");
  desc.add_options()
  ("exude", "Exude mesh after refinement")
  ("no_label_rebind", "Don't rebind cell labels in medit output")
  ("show_patches", "Show surface patches in medit output");

  
  po::options_description cmdline_options("Usage",1);
  cmdline_options.add(generic).add(mesh).add(functions).add(desc);

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
  po::notify(vm);

  if (vm.count("help") || argc < 2)
  {
    std::cout << cmdline_options << std::endl;
//    std::cout << "Generic options:\n"
//    <<"  --help                Produce help message\n\n"
//    <<"Mesh generation parameters:\n"
//    <<"  --facet_angle arg     Set facet angle bound\n"
//    <<"  --facet_size arg      Set facet size bound\n"
//    <<"  --facet_error arg     Set facet approximation error bound\n"
//    <<"  --tet_shape arg       Set tet radius-edge bound\n"
//    <<"  --tet_size arg        Set tet size bound\n\n"
//    <<"Implicit functions:\n"
//    <<"  --torus               Mesh torus function\n"
//    <<"  --sphere              Mesh sphere function\n"
//    <<"  --chair               Mesh chair function\n"
//    <<"  --tanglecube          Mesh tanglecube function\n"
//    <<"  --cube                Mesh cube function\n"
//    <<"  --ellipsoid           Mesh ellipsoid function\n"
//    <<"  --heart               Mesh heart function\n"
//    <<"  --octic               Mesh octic function\n\n"
//    <<"Options:\n"
//    <<"  --exude               Exude mesh after refinement\n"
//    <<"  --no_label_rebind     Don't rebind cell labels in medit output\n"
//    <<"  --show_patches        Show surface patches in medit output\n\n";
     
    return 1;
  }

  std::cout << "=========== Params ==========="<< std::endl;
  
  double facet_angle = set_arg<double>("facet_angle","Facet angle",vm);
  double facet_size = set_arg<double>("facet_size","Facet size",vm);
  double facet_error = set_arg<double>("facet_error","Facet approximation error",vm);

  double tet_shape = set_arg<double>("tet_shape","Tet shape (radius-edge)",vm);
  double tet_size = set_arg<double>("tet_size","Tet size",vm);
  

  
  // Define functions
  Function f1(&torus_function);
  Function f2(&sphere_function<3>);
  Function f3(&chair_function);
  Function f4(&tanglecube_function);
  Function f5(&cube_function);
  Function f6(&ellipsoid_function);
  Function f7(&heart_function);
  Function f8(&octic_function);

  Function_vector v;
  
  std::cout << "\nFunction(s): ";
  
  set_function(v,f1,"torus",vm);
  set_function(v,f2,"sphere",vm);
  set_function(v,f3,"chair",vm);
  set_function(v,f4,"tanglecube",vm);
  set_function(v,f5,"cube",vm);
  set_function(v,f6,"ellipsoid",vm);
  set_function(v,f7,"heart",vm);
  set_function(v,f8,"octic",vm);
    
  std::cout << "\n=============================="<< std::endl;
  std::cout << std::endl;
  
  if ( v.empty() )
  {
    std::cout << "No function set. Exit.\n";
    return 0;
  }
  
  // Domain (Warning: Sphere_3 constructor uses square radius !)
  Mesh_domain domain(v, K::Sphere_3(CGAL::ORIGIN, 7.*7.), 1e-8);

  // Set mesh criteria
  Facet_criteria facet_criteria(facet_angle, facet_size, facet_error); // angle, size, approximation
  Cell_criteria cell_criteria(tet_shape, tet_size); // radius-edge ratio, size
  Mesh_criteria criteria(facet_criteria, cell_criteria);

  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, vm.count("exude"));

  // Output
  std::ofstream medit_file("out.mesh");
  CGAL::output_to_medit(medit_file, c3t3, !vm.count("no_label_rebind"), vm.count("show_patches"));

  return 0;
}
