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
// file          : all_furthest_neighbors_2_demo.C
// chapter       : $CGAL_Chapter: Geometric Optimisation $
// package       : $CGAL_Package: Matrix_search $
// source        : mon_search.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// coordinator   : ETH Zurich (Bernd Gaertner <gaertner@inf.ethz.ch>)
//
// Demo program: All Furthest Neighbors for a Convex Polygon
// ============================================================================

#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>
#include <CGAL/all_furthest_neighbors_2.h>
#include <vector.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/IO/Window_stream.h>
#include <LEDA/point_set.h>
#include <iostream.h>

typedef double                                 FT;
typedef CGAL_Cartesian< FT >                   R;
typedef CGAL_Point_2< R >                      Point_2;
typedef CGAL_Polygon_traits_2< R >             P_traits;
typedef vector< Point_2 >                      Point_cont;
typedef vector< int >                          Index_cont;
typedef CGAL_Polygon_2< P_traits, Point_cont > Polygon_2;
typedef CGAL_Random_points_in_square_2<
  Point_2,
  CGAL_Creator_uniform_2< FT, Point_2 > >
Point_generator;
#include <LEDA/REDEFINE_NAMES.h>
typedef point                              LEDA_Point;
typedef point_set< Point_2 >               LEDA_Point_set_Point;
#include <LEDA/UNDEFINE_NAMES.h>
void
wait_for_button_release( leda_window& W)
{
  // wait until mouse button is released
  double x, y;
  int v;
  do {}
  while ( W.read_event( v, x, y) != button_release_event);
}

int
main()
{
  leda_window W;
  W.init( -1.25, 1.25, -1.25);
  W.display();

  // get points, last point with middle button:
  cout << "\nCGAL ALL FURTHEST NEIGHBORS DEMO\n"
       << "================================\n\n"
       << " compute all points furthest neighbors\n"
       << " for the vertices of a convex polygon\n\n"
       << " < click left mouse button to get the furthest"
       << " neighbor of a vertex >\n"
       << " < click any other mouse button to quit >\n"
       << endl;
  int number_of_points( 30);
  // generate random convex polygon:
  Polygon_2 p;
  CGAL_random_convex_set_2( number_of_points,
                            back_inserter( p),
                            Point_generator( 1));
  W << CGAL_RED << p;
  // compute all furthest neighbors:
  Index_cont neighbors;
  CGAL_all_furthest_neighbors(
    p.vertices_begin(),
    p.vertices_end(),
    back_inserter( neighbors));
  // output solution:
  LEDA_Point_set_Point query_points;
  for ( int k( 0); k < p.size(); ++k)
    query_points.insert(
      LEDA_Point( p[k].x(), p[k].y()),
      p[neighbors[k]]);
  
  W.set_mode( leda_xor_mode);
  W.set_fg_color( leda_blue);
  ps_item qp;
  // first click, no need to clear old query from screen:
  int last_button;
  {
    double x, y;
    last_button = W.read_mouse( x, y);
    qp = query_points.nearest_neighbor( LEDA_Point( x, y));
  }
  for (;;) {
    double x, y;
    if ( last_button != MOUSE_BUTTON( 1))
      break;
    W.draw_disc( query_points.key( qp), .02);
    W.draw_disc( query_points.inf( qp).x(),
                 query_points.inf( qp).y(),
                 .03);
    last_button = W.read_mouse( x, y);
    W.draw_disc( query_points.key( qp), .02);
    W.draw_disc( query_points.inf( qp).x(),
                 query_points.inf( qp).y(),
                 .03);
    qp = query_points.nearest_neighbor( LEDA_Point( x, y));
  }
  
  wait_for_button_release( W);

} // int main()
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

