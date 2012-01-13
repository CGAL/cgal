// Copyright (c) 2003   INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Julia Floetotto

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/_test_interpolation_functions_2.cpp>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K>            Dt;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K2;
typedef CGAL::Delaunay_triangulation_2<K2>            Dt2;

int main()
{
  std::cout << "Testing interpolation functions with 2D NN neighbors "
	    << std::endl;
  std::cout << " using Exact_predicates_exact_constructions_kernel: "
	    << std::endl ;
  _test_interpolation_functions_2_delaunay( Dt(), K::FT(1e-10));
  
  std::cout << "Testing interpolation functions with 2D NN neighbors "
	    << std::endl;
  std::cout << " using Exact_predicates_inexact_constructions_kernel: "
	    << std::endl ;
  _test_interpolation_functions_2_delaunay( Dt2(), K2::FT(1e-10));

  std::cout << "test_interpolation_functions_2 is finished" << std::endl;
  return 0;
}
