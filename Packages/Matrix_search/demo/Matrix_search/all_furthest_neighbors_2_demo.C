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
// file          : all_furthest_neighbors_2_demo.C
// chapter       : $CGAL_Chapter: Geometric Optimisation $
// package       : $CGAL_Package: Matrix_search $
// source        : mon_search.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// maintainer    : Michael Hoffmann <hoffmann@inf.ethz.ch>
// coordinator   : ETH
//
// Demo program: All Furthest Neighbors for a Convex Polygon
// ============================================================================


#include <CGAL/basic.h>

#if (defined(CGAL_USE_LEDA) || defined(CGAL_USE_CGAL_WINDOW))

#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>
#include <CGAL/distance_predicates_2.h>
#include <CGAL/all_furthest_neighbors_2.h>
#include <vector>
#include <CGAL/squared_distance_2.h>
#include <CGAL/IO/Window_stream.h>
#include <iostream>

using std::cout;
using std::endl;
using std::vector;
using std::back_inserter;
using CGAL::Cartesian;
using CGAL::Polygon_traits_2;
using CGAL::Creator_uniform_2;
using CGAL::Random_points_in_square_2;
using CGAL::random_convex_set_2;
using CGAL::has_smaller_dist_to_point;
using CGAL::all_furthest_neighbors_2;
using CGAL::cgalize;
using CGAL::RED;


typedef double                                 FT;
typedef CGAL::Window_stream                    Window;
typedef Cartesian< FT >                              R;
typedef CGAL::Point_2< R >                           Point;
typedef Polygon_traits_2< R >                        P_traits;
typedef vector< Point >                              Point_cont;
typedef vector< int >                                Index_cont;
typedef CGAL::Polygon_2< P_traits, Point_cont >      Polygon;
typedef Creator_uniform_2< FT, Point >               Creator;
typedef Random_points_in_square_2< Point, Creator >  Point_generator;
void
wait_for_button_release( Window& W)
{
  // wait until mouse button is released
  double x, y;
  int v;
  do {}
#ifdef CGAL_USE_LEDA
  while ( W.read_event( v, x, y) != button_release_event);
#else
  while ( W.read_event( v, x, y) != CGAL::button_release_event);
#endif // CGAL_USE_LEDA
}

int
main()
{
  Window W;
  cgalize(W);
  W.init(-1.25, 1.25, -1.25);
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
  int number_of_points(30);
  // generate random convex polygon:
  Polygon p;
  random_convex_set_2( number_of_points,
                       back_inserter( p),
                       Point_generator( 1));
  W << RED << p;
  // compute all furthest neighbors:
  Index_cont neighbors;
  all_furthest_neighbors_2(
    p.vertices_begin(),
    p.vertices_end(),
    back_inserter( neighbors));
  // output solution:
  #ifdef CGAL_USE_LEDA
  W.set_mode(leda_xor_mode);
  #else
  W.set_mode(CGAL::xor_mode);
  #endif // CGAL_USE_LEDA
  W << CGAL::BLUE;
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

  return 0;
} // int main()

#else

#include <iostream>

int main()
{
  std::cerr << "This demo requires LEDA." << std::endl;
  return 0;
}

#endif

// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

