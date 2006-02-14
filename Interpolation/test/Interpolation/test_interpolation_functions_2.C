// Copyright (c) 2003   INRIA Sophia-Antipolis (France).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Julia Floetotto

#include <CGAL/basic.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/_test_interpolation_functions_2.C>

struct K : CGAL::Exact_predicates_exact_constructions_kernel {};
typedef CGAL::Delaunay_triangulation_2<K>            Dt;

struct K2 : CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef CGAL::Delaunay_triangulation_2<K2>            Dt2;

int main()
{
  std::cout << "Testing interpolation functions with 2D NN neighbors "
	    << std::endl;
  std::cout << " using Exact_predicates_exact_constructions_kernel: "
	    << std::endl ;
  _test_interpolation_functions_2_delaunay( Dt(), K::FT(0));

  std::cout << "Testing interpolation functions with 2D NN neighbors "
	    << std::endl;
  std::cout << " using Exact_predicates_inexact_constructions_kernel: "
	    << std::endl ;
  _test_interpolation_functions_2_delaunay( Dt2(), K2::FT(1e-10));

  std::cout << "test_interpolation_functions_2 is finished" << std::endl;
  return 0;
}
