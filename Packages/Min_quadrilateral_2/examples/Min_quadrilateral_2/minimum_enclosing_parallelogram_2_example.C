// Copyright (c) 1999-2003  ETH Zurich (Switzerland).
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
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch> and
//                 Emo Welzl <emo@inf.ethz.ch>

#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>
#include <CGAL/min_quadrilateral_2.h>
#include <iostream>

struct Kernel : public CGAL::Cartesian<double> {};

typedef Kernel::Point_2                           Point_2;
typedef Kernel::Line_2                            Line_2;
typedef CGAL::Polygon_2<Kernel>                   Polygon_2;
typedef CGAL::Random_points_in_square_2<Point_2>  Generator;

int main()
{
  // build a random convex 20-gon p
  Polygon_2 p;
  CGAL::random_convex_set_2(20, std::back_inserter(p), Generator(1.0));
  std::cout << p << std::endl;

  // compute the minimal enclosing parallelogram p_m of p
  Polygon_2 p_m;
  CGAL::min_parallelogram_2(
    p.vertices_begin(), p.vertices_end(), std::back_inserter(p_m));
  std::cout << p_m << std::endl;

  return 0;
} 
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

