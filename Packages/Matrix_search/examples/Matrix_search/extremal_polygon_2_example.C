// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
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
// file          : extremal_polygon_2_example.C
// chapter       : $CGAL_Chapter: Geometric Optimisation $
// package       : $CGAL_Package: Matrix_search $
// source        : mon_search.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// coordinator   : ETH Zurich (Bernd Gaertner <gaertner@inf.ethz.ch>)
//
// Example program: Compute extremal polygons of a convex polygon
// ============================================================================

#ifndef CGAL_CARTESIAN_H
#include <CGAL/Cartesian.h>
#endif // CGAL_CARTESIAN_H
#ifndef CGAL_POINT_2_H
#include <CGAL/Point_2.h>
#endif // CGAL_POINT_2_H
#ifndef CGAL_POLYGON_2_H
#include <CGAL/Polygon_2.h>
#endif // CGAL_POLYGON_2_H
#ifndef CGAL_POINT_GENERATORS_2_H
#include <CGAL/point_generators_2.h>
#endif // CGAL_POINT_GENERATORS_2_H
#ifndef CGAL_RANDOM_CONVEX_SET_2_H
#include <CGAL/random_convex_set_2.h>
#endif // CGAL_RANDOM_CONVEX_SET_2_H
#ifndef CGAL_EXTREMAL_POLYGON_2_H
#include <CGAL/extremal_polygon_2.h>
#endif // CGAL_EXTREMAL_POLYGON_2_H
#include <iostream>
#include <vector>

using namespace std;
using CGAL::random_convex_set_2;
using CGAL::maximum_area_inscribed_k_gon;

typedef double                                FT;
typedef CGAL::Cartesian< FT >                 R;
typedef CGAL::Point_2< R >                    Point;
typedef CGAL::Polygon_traits_2< R >           P_traits;
typedef vector< Point >                       Cont;
typedef CGAL::Polygon_2< P_traits, Cont >     Polygon;
typedef CGAL::Creator_uniform_2< FT, Point >  Creator;
typedef CGAL::Random_points_in_square_2< Point, Creator >
  Point_generator;

int main() {

  Polygon p;
  int number_of_points( 10);
  int k( 5);

  random_convex_set_2( number_of_points,
                       back_inserter( p),
                       Point_generator( 1));
  cout << "Generated Polygon:\n" << p << endl;

  Polygon k_gon;
  maximum_area_inscribed_k_gon(
    p.vertices_begin(),
    p.vertices_end(),
    k,
    back_inserter( k_gon));
  cout << "Maximum area " << k << "-gon:\n" << k_gon << endl;

} 
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

