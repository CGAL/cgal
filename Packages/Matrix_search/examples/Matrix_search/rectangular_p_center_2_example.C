// ============================================================================
//
// Copyright (c) 1998, 1999, 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision $
// release_date  : $CGAL_Date $
//
// file          : rectangular_p_center_2_example.C
// chapter       : $CGAL_Chapter: Geometric Optimisation $
// package       : $CGAL_Package: Matrix_search $
// source        : pcenter.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// maintainer    : Michael Hoffmann <hoffmann@inf.ethz.ch>
// coordinator   : ETH
//
// 2-4-Centering Axis-Parallel 2D-Rectangles - example program
// ============================================================================

#include <CGAL/Cartesian.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/rectangular_p_center_2.h>
#include <CGAL/IO/Ostream_iterator.h>
#include <CGAL/algorithm.h>
#include <iostream>
#include <algorithm>
#include <vector>

using namespace std;
using CGAL::set_pretty_mode;
using CGAL::rectangular_p_center_2;

typedef double                                             FT;
typedef CGAL::Cartesian< FT >                              K;
typedef K::Point_2                                         Point;
typedef std::vector< Point >                               Cont;
typedef CGAL::Creator_uniform_2< FT, Point >               Creator;
typedef CGAL::Random_points_in_square_2< Point, Creator >  Point_generator;
typedef CGAL::Ostream_iterator< Point, ostream >    Ostream_iterator_point;

int main() {

  int number_of_points(10);
  int p(2);
  Ostream_iterator_point cout_ip(cout);
  set_pretty_mode(cout);

  Cont points;
  CGAL::copy_n(Point_generator(1),
               number_of_points,
               back_inserter(points));
  cout << "Generated Point Set:" << endl;
  copy(points.begin(), points.end(), cout_ip);

  Cont centers;
  FT p_radius;
  rectangular_p_center_2(
    points.begin(),
    points.end(),
    back_inserter(centers),
    p_radius,
    3);

  cout << "\n\n" << p << "-centers:" << endl;
  copy(centers.begin(), centers.end(), cout_ip);
  cout << "\n\n" << p << "-radius = " << p_radius << endl;

  return 0;
} 
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

