// Copyright (c) 2003   INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Julia Floetotto

#include <CGAL/_test_interpolation_functions_2.cpp>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/array.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Origin.h>

#include <iostream>

typedef CGAL::Exact_predicates_exact_constructions_kernel      EPECK;
typedef CGAL::Exact_predicates_inexact_constructions_kernel    EPICK;

template <typename V, typename G>
struct Value_and_gradient
{
  Value_and_gradient() : value(), gradient(CGAL::NULL_VECTOR) {}

  V value;
  G gradient;
};

template<typename Kernel>
void test_interpolation_functions()
{
  // For the Delaunay triangulation, values and gradients (three different data sets)
  // are stored directly in the vertices
  typedef typename Kernel::FT                                      Coord_type;
  typedef typename Kernel::Vector_2                                Vector;

  typedef CGAL::Triangulation_vertex_base_with_info_2<
                  std::array<
                    Value_and_gradient<Coord_type, Vector>, 3>,
                  Kernel>                                          Vb;
  typedef CGAL::Triangulation_data_structure_2<Vb>                 Tds;
  typedef CGAL::Delaunay_triangulation_2<Kernel, Tds>              Delaunay_triangulation;

  typedef CGAL::Regular_triangulation_2<Kernel>                    Regular_triangulation;

  std::cout << "Testing interpolation functions with 2D NN neighbors " << std::endl;
  _test_interpolation_functions_2_Delaunay(Delaunay_triangulation(), Coord_type(1e-10));

  std::cout << "Testing interpolation functions with 2D RN neighbors " << std::endl;
  _test_interpolation_functions_2_regular(Regular_triangulation(), Coord_type(1e-10));
}

int main()
{
  std::cout << "--------------------------------------------" << std::endl;
  std::cout << "Testing with EPECK" << std::endl;
  test_interpolation_functions<EPECK>();

  std::cout << "--------------------------------------------" << std::endl;
  std::cout << "Testing with EPICK" << std::endl;
  test_interpolation_functions<EPICK>();

  std::cout << "test_interpolation_functions_2 is finished" << std::endl;

  return EXIT_SUCCESS;
}
