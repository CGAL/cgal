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
// file          : demo/Generator/rcs_demo.C
// package       : $CGAL_Package: Generator 2.12 (28 Jul 1999) $
// source        : src/rcs/rcs.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// coordinator   : ETH Zurich (Bernd Gaertner <gaertner@inf.ethz.ch>)
//
// Random Convex Point Sets: Demo Program
// ============================================================================

#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <CGAL/IO/Window_stream.h>

#if defined(CGAL_USE_CGAL_WINDOW)
#define leda_window CGAL::window
#endif

using CGAL::cgalize;

using std::cout;
using std::cerr;
using std::endl;
using std::flush;
using std::vector;
using std::back_inserter;
using CGAL::Cartesian;
using CGAL::Polygon_traits_2;
using CGAL::Creator_uniform_2;
using CGAL::Random_points_in_square_2;
using CGAL::set_pretty_mode;
using CGAL::random_convex_set_2;

typedef Cartesian< double >                          R;
typedef CGAL::Point_2< R >                           Point;
typedef Polygon_traits_2< R >                        P_traits;
typedef vector< Point >                              Cont;
typedef CGAL::Polygon_2< P_traits, Cont >            Polygon_2;
typedef Creator_uniform_2< double, Point >           Creator;
typedef Random_points_in_square_2< Point, Creator >  Point_generator;

int
main(int argc, char* argv[])
{
  // this is not initialized on MIPSPRO:
  set_pretty_mode( cout);
  set_pretty_mode( cerr);

  // take #points from command line:
  int n;
  if ( argc < 2 || (n = atoi( argv[1])) < 3) {
    cerr << "usage: " << argv[0] << " \"#points\" (>= 3)" << endl;
    return 1;
  }

  cout << "Test random_convex_set_2:\n" << endl;

  // build random n-gon:
  cout << "constructing random " << n << "-gon ..." << flush;
  Polygon_2 p;
  random_convex_set_2( n, back_inserter( p), Point_generator( 1));
  cout << " done." << endl;

  // output polygon:
  cout << "\nHere is the result:" << endl;

  leda_window W;
  cgalize( W);
  W.display();
  W.init( -1.05, 1.05, -1.05);
  W << p;

  // wait for mouse-click:
  W.read_mouse();

  // check convexity:
  if ( ! p.is_convex()) {
    cerr << "ERROR: polygon is not convex." << endl;
    return 1;
  }

  cout << "done." << endl;
  return 0;
} // int main( argc, argv)

// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

