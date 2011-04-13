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
// file          : minimum_enclosing_quadrilateral_2_test_H.C
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
// Test Program: Computing minimum enclosing quadrilaterals
// ============================================================================

#include <CGAL/Homogeneous.h>
#include <CGAL/Point_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>
#include <CGAL/min_quadrilateral_2.h>
#include <vector>
#include <iostream>


using CGAL::Polygon_traits_2;
using CGAL::Creator_uniform_2;
using CGAL::Random_points_in_square_2;
using CGAL::random_convex_set_2;
using CGAL::min_rectangle_2;
using CGAL::min_parallelogram_2;
using CGAL::min_strip_2;
using std::back_inserter;
using std::vector;
using std::cout;
using std::endl;
typedef CGAL::Homogeneous< double >                       R;
typedef R::Point_2                                     Point_2;
typedef R::Line_2                                      Line_2;
typedef Polygon_traits_2< R >                          P_traits;
typedef vector< Point_2 >                              Cont;
typedef CGAL::Polygon_2< P_traits, Cont >              Polygon_2;
typedef Creator_uniform_2< double, Point_2 >           Creator;
typedef Random_points_in_square_2< Point_2, Creator >  Point_generator;

int main()
{
  CGAL::set_pretty_mode(cout);

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

