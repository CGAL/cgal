// ============================================================================
//
// Copyright (c) 1999, 2000 The CGAL Consortium
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
// file          : minimum_enclosing_strip_2_example.C
// chapter       : $CGAL_Chapter: Geometric Optimisation $
// package       : $CGAL_Package: Min_quadrilaterals $
// source        : oops.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch> and
//                 Emo Welzl <emo@inf.ethz.ch>
//
// maintainer    : Michael Hoffmann <hoffmann@inf.ethz.ch>
// coordinator   : ETH
//
// Example Program: Computing a minimum enclosing strip
// ============================================================================

#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>
#include <CGAL/min_quadrilateral_2.h>
#include <vector>
#include <iostream>

using CGAL::Random_points_in_square_2;
using CGAL::random_convex_set_2;
using CGAL::min_strip_2;
using std::back_inserter;
using std::cout;
using std::endl;

typedef CGAL::Cartesian< double >                      R;
typedef R::Point_2                                     Point_2;
typedef R::Line_2                                      Line_2;
typedef CGAL::Polygon_traits_2< R >                    P_traits;
typedef std::vector< Point_2 >                         Cont;
typedef CGAL::Polygon_2< P_traits, Cont >              Polygon_2;
typedef CGAL::Creator_uniform_2< double, Point_2 >     Creator;
typedef Random_points_in_square_2< Point_2, Creator >  Point_generator;

int main()
{
  // build a random convex 20-gon p
  Polygon_2 p;
  random_convex_set_2(20, back_inserter(p), Point_generator(1.0));
  cout << p << endl;

  // compute the minimal enclosing strip p_m of p
  Line_2 p_m[2];
  min_strip_2(p.vertices_begin(), p.vertices_end(), p_m);
  cout << p_m[0] << "\n" << p_m[1] << endl;

  return 0;
} 
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

