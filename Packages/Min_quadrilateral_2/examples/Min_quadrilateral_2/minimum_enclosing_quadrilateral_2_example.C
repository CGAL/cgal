#line 1634 "oops.aw"
#line 18 "code_formatting.awi"
// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// file          : minimum_enclosing_quadrilateral_2_example.C
// chapter       : $CGAL_Chapter: Geometric Optimisation $
// package       : $CGAL_Package: Min_quadrilaterals $
// source        : oops.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch> and
//                 Emo Welzl <emo@inf.ethz.ch>
//
// coordinator   : ETH Zurich (Bernd Gaertner <gaertner@inf.ethz.ch>)
//
// Example Program: Computing a minimum enclosing four-gon
// ============================================================================

#line 1638 "oops.aw"
#line 1226 "oops.aw"
#line 1546 "oops.aw"
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>
#include <CGAL/minimum_enclosing_quadrilateral_2.h>
#include <vector>
#include <iostream>


#line 1581 "oops.aw"
using CGAL::Polygon_traits_2;
using CGAL::Creator_uniform_2;
using CGAL::Random_points_in_square_2;
using CGAL::random_convex_set_2;
using CGAL::minimum_enclosing_rectangle_2;
using CGAL::minimum_enclosing_parallelogramm_2;
using CGAL::minimum_enclosing_strip_2;
using std::back_inserter;
using std::vector;
using std::cout;
using std::endl;
#line 1228 "oops.aw"
typedef CGAL::Cartesian< double >                            R;
#line 1595 "oops.aw"
typedef R::Point_2                                     Point_2;
typedef R::Line_2                                      Line_2;
typedef Polygon_traits_2< R >                          P_traits;
typedef vector< Point_2 >                              Cont;
typedef CGAL::Polygon_2< P_traits, Cont >              Polygon_2;
typedef Creator_uniform_2< double, Point_2 >           Creator;
typedef Random_points_in_square_2< Point_2, Creator >  Point_generator;
#line 1230 "oops.aw"

int main()
{
  // build a random convex 20-gon p
  Polygon_2 p;
  random_convex_set_2(20, back_inserter(p), Point_generator(1.0));
  cout << p << endl;

  // compute the minimal enclosing rectangle p_m of p
  Polygon_2 p_m;
  minimum_enclosing_rectangle_2(
    p.vertices_begin(), p.vertices_end(), back_inserter(p_m));
  cout << p_m << endl;

  return 0;
} 
#line 1639 "oops.aw"
#line 12 "code_formatting.awi"
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

