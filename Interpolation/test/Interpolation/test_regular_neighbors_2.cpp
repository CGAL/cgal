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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Julia Floetotto

#include <CGAL/_test_regular_neighbors_2.cpp>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/double.h>
#include <CGAL/Regular_triangulation_2.h>

#include <iostream>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Regular_triangulation_2<K>                  Rt1;

int main()
{
  std::cout << "Testing NN_neighbors_2 " << std::endl;
  std::cout << " using Exact_predicates_exact_constructions_kernel: " << std::endl;
  _test_regular_neighbors_2( Rt1() );

  return EXIT_SUCCESS;
}
