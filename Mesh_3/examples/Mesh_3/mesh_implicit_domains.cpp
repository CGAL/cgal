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
// $URL: $
// $Id: $
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
#include "implicit_functions.h"

// IO
#include <CGAL/IO/File_medit.h>

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


int main()
{
  // Define functions
  Function f1(&torus_function);
  Function f2(&sphere_function<3>);

  Function_vector v;
  v.push_back(&f1);
  v.push_back(&f2);

  // Domain (Warning: Sphere_3 constructor uses square radius !)
  Mesh_domain domain(v, K::Sphere_3(CGAL::ORIGIN, 5.*5.), 1e-6);

  // Set mesh criteria
  Facet_criteria facet_criteria(30, 0.15, 0.01); // angle, size, approximation
  Cell_criteria cell_criteria(3., 0.5); // radius-edge ratio, size
  Mesh_criteria criteria(facet_criteria, cell_criteria);

  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

  // Output
  std::ofstream medit_file("out.mesh");
  CGAL::output_to_medit(medit_file, c3t3);

  return 0;
}
