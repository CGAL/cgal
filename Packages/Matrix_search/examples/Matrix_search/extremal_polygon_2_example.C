// Copyright (c) 1998-2003  ETH Zurich (Switzerland).
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
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>

#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>
#include <CGAL/extremal_polygon_2.h>
#include <iostream>
#include <vector>

typedef double                                    FT;

struct Kernel : public CGAL::Cartesian<FT> {};

typedef Kernel::Point_2                           Point;
typedef std::vector<int>                          Index_cont;
typedef CGAL::Polygon_2<Kernel>                   Polygon;
typedef CGAL::Random_points_in_square_2<Point>    Generator;

int main() {

  int n = 10;
  int k = 5;

  // generate random convex polygon:
  Polygon p;
  CGAL::random_convex_set_2(n, back_inserter(p), Generator(1));
  std::cout << "Generated Polygon:\n" << p << std::endl;

  // compute maximum area incribed k-gon of p:
  Polygon k_gon;
  CGAL::maximum_area_inscribed_k_gon_2(
    p.vertices_begin(), p.vertices_end(), k, std::back_inserter(k_gon));
  std::cout << "Maximum area " << k << "-gon:\n"
            << k_gon << std::endl;

  return 0;
} 
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

