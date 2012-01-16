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

#include <CGAL/basic.h>
#include <CGAL/double.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_euclidean_traits_2.h>
#include <CGAL/Regular_triangulation_2.h>

#include <CGAL/_test_regular_neighbors_2.cpp>


struct K : CGAL::Exact_predicates_exact_constructions_kernel {};
typedef double W;

typedef CGAL::Regular_triangulation_euclidean_traits_2<K,W>  Gt1;
typedef CGAL::Regular_triangulation_2<Gt1>                   Rt1;

int main()
{
  std::cout << "Testing NN_neighbors_2 " << std::endl;
  std::cout << " with Exact_predicates_exact_constructions_kernel: " << std::endl ;
  _test_regular_neighbors_2( Rt1() );

  return 0;
}
