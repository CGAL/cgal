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

#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>
#include <CGAL/extremal_polygon_2.h>
#include <iostream.h>
#include <vector.h>

int main() {

  typedef double                            FT;
  typedef CGAL_Cartesian< FT >              R;
  typedef CGAL_Point_2< R >                 Point_2;
  typedef CGAL_Polygon_traits_2< R >        P_traits;
  typedef vector< Point_2 >                 Cont;
  typedef CGAL_Polygon_2< P_traits, Cont >  Polygon_2;
  typedef CGAL_Random_points_in_square_2<
    Point_2,
    CGAL_Creator_uniform_2< FT, Point_2 > >
  Point_generator;

  Polygon_2 p;
  int number_of_points( 10);
  int k( 5);

  CGAL_random_convex_set_2( number_of_points,
                            back_inserter( p),
                            Point_generator( 1));
  cout << "Generated Polygon:\n" << p << endl;

  Polygon_2 k_gon;
  CGAL_maximum_area_inscribed_k_gon(
    p.vertices_begin(),
    p.vertices_end(),
    k,
    back_inserter( k_gon));
  cout << "Maximum area " << k << "-gon:\n" << k_gon << endl;

} 
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

