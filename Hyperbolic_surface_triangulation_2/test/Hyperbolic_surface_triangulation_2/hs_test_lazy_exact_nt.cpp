// Copyright (c) 2024
// INRIA Nancy (France), and Université Gustave Eiffel Marne-la-Vallee (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Vincent Despré, Loïc Dubois, Monique Teillaud

#include <iostream>
#include <sstream>

#include <CGAL/Gmpq.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/Hyperbolic_surfaces_traits_2.h>
#include <CGAL/Hyperbolic_fundamental_domain_2.h>
#include <CGAL/Hyperbolic_fundamental_domain_factory_2.h>
#include <CGAL/Hyperbolic_surface_triangulation_2.h>

using namespace CGAL;

typedef Cartesian<Lazy_exact_nt<Gmpq>>                                  Kernel;
typedef Hyperbolic_Delaunay_triangulation_traits_2<Kernel>              ParentTraits;
typedef Hyperbolic_surfaces_traits_2<ParentTraits>                      Traits;
typedef Hyperbolic_fundamental_domain_2<Traits>                         Domain;
typedef Hyperbolic_fundamental_domain_factory_2<Traits>                 Factory;
typedef Hyperbolic_surface_triangulation_2<Traits>                      Triangulation;

typedef typename Traits::Hyperbolic_point_2                             Point;

int main() {
  Factory factory (3459);
  Domain domain = factory.generate_domain_g2();
  Triangulation triangulation0 = Triangulation(domain);

  assert( triangulation0.is_valid() );

  Triangulation triangulation1 = Triangulation(triangulation0.get_combinatorial_map_ref());
  assert( ! triangulation1.has_anchor() );

  Triangulation triangulation (triangulation0);
  assert( triangulation.has_anchor() );

  std::stringstream buffer;
  buffer << triangulation;
  buffer >> triangulation;

  std::vector<std::tuple<typename Triangulation::Combinatorial_map_with_cross_ratios::Dart_const_handle,Point,Point,Point>> input_not_centered;
  std::vector<std::tuple<typename Triangulation::Combinatorial_map_with_cross_ratios::Dart_const_handle,Point,Point,Point>> input_centered;

  input_not_centered = triangulation.lift(false);
  input_centered = triangulation.lift();

  triangulation.make_delaunay();

  assert( triangulation.is_delaunay() );

  std::vector<std::tuple<typename Triangulation::Combinatorial_map_with_cross_ratios::Dart_const_handle,Point,Point,Point>> output_not_centered;
  std::vector<std::tuple<typename Triangulation::Combinatorial_map_with_cross_ratios::Dart_const_handle,Point,Point,Point>> output_centered;

  output_not_centered = triangulation.lift(false);
  output_centered = triangulation.lift();

  Triangulation::Combinatorial_map_with_cross_ratios& cmap = triangulation.get_combinatorial_map_ref();
  Triangulation::Anchor& anchor = triangulation.get_anchor_ref();
  assert( cmap.is_dart_used(anchor.dart) );

  std::cout << "printing triangulation for test purposes : " << std::endl << triangulation;

  return 0;
}
