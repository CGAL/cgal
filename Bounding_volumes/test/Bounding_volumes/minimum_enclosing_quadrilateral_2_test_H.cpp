// Copyright (c) 1999-2003  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch> and
//                 Emo Welzl <emo@inf.ethz.ch>

#include <CGAL/Homogeneous.h>
#include <CGAL/Point_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>
#include <CGAL/min_quadrilateral_2.h>
#include <iostream>


using CGAL::random_convex_set_2;
using CGAL::min_rectangle_2;
using CGAL::min_parallelogram_2;
using CGAL::min_strip_2;
using std::back_inserter;
using std::cout;
using std::endl;
typedef CGAL::Homogeneous< double >                       Kernel;
typedef Kernel::Point_2                           Point_2;
typedef Kernel::Line_2                            Line_2;
typedef CGAL::Polygon_2<Kernel>                   Polygon_2;
typedef CGAL::Random_points_in_square_2<Point_2>  Point_generator;

int main()
{
  CGAL::IO::set_pretty_mode(cout);

  // build a random convex 20-gon p
  Polygon_2 p;
  p.push_back(Point_2(1,2,3));
  p.push_back(Point_2(3,5,2));
  p.push_back(Point_2(-1,10,4));
  cout << "Input:\n" << p << endl;

  // compute the minimal enclosing rectangle p_m of p
  Polygon_2 p_r;
  min_rectangle_2(
    p.vertices_begin(), p.vertices_end(), back_inserter(p_r));
  cout << "Min_rectangle:\n" << p_r << endl;

  // compute the minimal enclosing parallelogram p_p of p
  Polygon_2 p_p;
  min_parallelogram_2(
    p.vertices_begin(), p.vertices_end(), back_inserter(p_p));
  cout << "Min_parallelogram:\n" << p_p << endl;

  // compute the minimal enclosing strip p_s of p
  Line_2 lines[2];
  min_strip_2(p.vertices_begin(), p.vertices_end(), lines);
  cout << "Min_strip:\n" << lines[0] << "\n" << lines[1] << endl;

  return 0;
}
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

