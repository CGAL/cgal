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

#include <CGAL/Complex_without_sqrt.h>
#include <CGAL/Hyperbolic_fundamental_domain_2.h>
#include <CGAL/Hyperbolic_surface_triangulation_2.h>
#include <CGAL/Hyperbolic_surfaces_traits_2.h>

#include <iostream>
#include <sstream>

#include <CGAL/Gmpq.h>

using namespace CGAL;

typedef Gmpq FT;

typedef Hyperbolic_surfaces_traits_2<FT>                            Traits;
typedef Hyperbolic_fundamental_domain_2<Traits>                     Domain;
typedef Hyperbolic_surface_triangulation_2<Traits>                  Triangulation;

typedef typename Traits::Point_2                                    Point;
typedef Complex_without_sqrt<FT>                                    Complex;

Domain build_domain(){
  //
  std::vector<Point> vertices;
  Complex z0 = Complex (FT(809,10000),FT(0));
  Complex z1 = Complex (FT(7359,10000),FT(1877,10000));
  Complex z2 = Complex (FT(-999,2500),FT(881,1000));
  Complex z3 = Complex (FT("-22088524601252853411192791001942853611410938513/24711029456888649611435724068315791591836010000"),FT("9482675065452890527617859332378101016513362487/24711029456888649611435724068315791591836010000"));
  vertices.push_back( Point(z0) );
  vertices.push_back( Point(z1) );
  vertices.push_back( Point(z2) );
  vertices.push_back( Point(z3) );
  vertices.push_back( Point(-z0) );
  vertices.push_back( Point(-z1) );
  vertices.push_back( Point(-z2) );
  vertices.push_back( Point(-z3) );

  std::vector<int> pairings;
  for (int k=0; k<8; k++){
    pairings.push_back((k+4)%8);
  }

  return Domain(vertices, pairings);
}

int main() {
  Domain domain = build_domain();
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
