#line 1475 "mon_search.aw"
#line 18 "code_formatting.awi"
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

#line 1479 "mon_search.aw"
#line 509 "afn.awi"
#line 542 "afn.awi"
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>
#include <CGAL/distance_predicates_2.h>
#include <CGAL/all_furthest_neighbors_2.h>
#include <vector>
#line 510 "afn.awi"
#line 553 "afn.awi"
#include <CGAL/squared_distance_2.h>
#include <CGAL/IO/Window_stream.h>
#include <iostream>
#line 511 "afn.awi"

using std::cout;
using std::endl;
#line 559 "afn.awi"
using std::vector;
using std::back_inserter;
using CGAL::Cartesian;
using CGAL::Polygon_traits_2;
using CGAL::Creator_uniform_2;
using CGAL::Random_points_in_square_2;
using CGAL::random_convex_set_2;
using CGAL::has_smaller_dist_to_point;
using CGAL::all_furthest_neighbors;
#line 515 "afn.awi"
using CGAL::cgalize;
using CGAL::RED;

typedef double                                 FT;
#line 571 "afn.awi"
typedef Cartesian< FT >                              R;
typedef CGAL::Point_2< R >                           Point;
typedef Polygon_traits_2< R >                        P_traits;
typedef vector< Point >                              Point_cont;
typedef vector< int >                                Index_cont;
typedef CGAL::Polygon_2< P_traits, Point_cont >      Polygon;
typedef Creator_uniform_2< FT, Point >               Creator;
typedef Random_points_in_square_2< Point, Creator >  Point_generator;
#line 520 "afn.awi"
#line 659 "afn.awi"
void
wait_for_button_release( leda_window& W)
{
  // wait until mouse button is released
  double x, y;
  int v;
  do {}
  while ( W.read_event( v, x, y) != button_release_event);
}
#line 521 "afn.awi"

int
main()
{
  leda_window W;
  cgalize(W);
  W.init(-1.25, 1.25, -1.25);
  W.display();

  #line 582 "afn.awi"
  // get points, last point with middle button:
  cout << "\nCGAL ALL FURTHEST NEIGHBORS DEMO\n"
       << "================================\n\n"
       << " compute all points furthest neighbors\n"
       << " for the vertices of a convex polygon\n\n"
       << " < click left mouse button to get the furthest"
       << " neighbor of a vertex >\n"
       << " < click any other mouse button to quit >\n"
       << endl;
#line 531 "afn.awi"
  int number_of_points(30);
  #line 594 "afn.awi"
  // generate random convex polygon:
  Polygon p;
  random_convex_set_2( number_of_points,
                       back_inserter( p),
                       Point_generator( 1));
#line 533 "afn.awi"
  W << RED << p;
  #line 602 "afn.awi"
  // compute all furthest neighbors:
  Index_cont neighbors;
  all_furthest_neighbors(
    p.vertices_begin(),
    p.vertices_end(),
    back_inserter( neighbors));
#line 535 "afn.awi"
  #line 617 "afn.awi"
  // output solution:
  W.set_mode(leda_xor_mode);
  W.set_fg_color(leda_blue);
  // first click, no need to clear old query from screen:
  int last_button;
  int nearest = 0;
  {
    double x, y;
    last_button = W.read_mouse(x, y);
    Point sp(x, y);
    for (int k(1); k < p.size(); ++k)
      if (has_smaller_dist_to_point(sp, p[k], p[nearest]))
        nearest = k;
  }
  for (;;) {
    double x, y;
    if (last_button != MOUSE_BUTTON(1))
      break;
    W.draw_disc(p[nearest].x(), p[nearest].y(), .02);
    W.draw_disc(p[neighbors[nearest]].x(),
                p[neighbors[nearest]].y(),
                .03);
    last_button = W.read_mouse(x, y);
    W.draw_disc(p[nearest].x(), p[nearest].y(), .02);
    W.draw_disc(p[neighbors[nearest]].x(),
                p[neighbors[nearest]].y(),
                .03);
    nearest = 0;
    Point sp(x, y);
    for (int k(1); k < p.size(); ++k)
      if (has_smaller_dist_to_point(sp, p[k], p[nearest]))
        nearest = k;
  }
  
  wait_for_button_release(W);
#line 536 "afn.awi"

  return 0;
} // int main()
#line 1480 "mon_search.aw"
#line 12 "code_formatting.awi"
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

