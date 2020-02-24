// Copyright (c) 2003   INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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
