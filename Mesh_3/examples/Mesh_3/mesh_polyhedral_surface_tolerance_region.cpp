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
//
//******************************************************************************


#include <CGAL/AABB_intersections.h>

#include <debug.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Mesh_3/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_3/polyhedral_to_labeled_function_wrapper.h>
#include <CGAL/make_mesh_3.h>

// IO
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/File_medit.h>


// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::Mesh_3::Polyhedral_tolerance_to_labeled_function_wrapper<Polyhedron, K> Polyhedral_wrapper;
typedef CGAL::Mesh_3::Labeled_mesh_domain_3<Polyhedral_wrapper, K> Mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Mesh Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
typedef Mesh_criteria::Facet_criteria Facet_criteria;
typedef Mesh_criteria::Cell_criteria Cell_criteria;


int main()
{
  // Create and fill polyhedron
  Polyhedron polyhedron;
  std::ifstream input("data/elephant.off");
  input >> polyhedron;
  input.close();

  // Implicit function built by AABB_tree projection queries
  Polyhedral_wrapper polyhedral_wrapper(polyhedron, 0.03);
  Mesh_domain domain(polyhedral_wrapper, polyhedral_wrapper.bounding_sphere(), 1e-4);

  // Set mesh criteria
  Facet_criteria facet_criteria(20, 5, 0.01); // angle, size, approximation
  Cell_criteria cell_criteria(4, 0.05); // radius-edge ratio, size
  Mesh_criteria criteria(facet_criteria, cell_criteria);

  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

  // Output
  std::ofstream medit_file("out.mesh");
  CGAL::output_to_medit(medit_file, c3t3);

  return 0;
}
