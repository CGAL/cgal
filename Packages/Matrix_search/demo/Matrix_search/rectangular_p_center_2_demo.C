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
// file          : rectangular_p_center_2_demo.C
// chapter       : $CGAL_Chapter: Geometric Optimisation $
// package       : $CGAL_Package: Matrix_search $
// source        : pcenter.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// coordinator   : ETH Zurich (Bernd Gaertner <gaertner@inf.ethz.ch>)
//
// Demo: 2-4-Centering Axis-Parallel 2D-Rectangles
// ============================================================================

#ifndef CGAL_CARTESIAN_H
#include <CGAL/Cartesian.h>
#endif // CGAL_CARTESIAN_H
#ifndef CGAL_ISO_RECTANGLE_2_H
#include <CGAL/Iso_rectangle_2.h>
#endif // CGAL_ISO_RECTANGLE_2_H
#ifndef CGAL_POINT_2_H
#include <CGAL/Point_2.h>
#endif // CGAL_POINT_2_H
#ifndef CGAL_IO_LEDA_WINDOW_H
#include <CGAL/IO/leda_window.h>
#endif // CGAL_IO_LEDA_WINDOW_H
#ifndef CGAL_IO_OSTREAM_ITERATOR_H
#include <CGAL/IO/Ostream_iterator.h>
#endif // CGAL_IO_OSTREAM_ITERATOR_H
#ifndef CGAL_IO_ISTREAM_ITERATOR_H
#include <CGAL/IO/Istream_iterator.h>
#endif // CGAL_IO_ISTREAM_ITERATOR_H
#ifndef CGAL_RECTANGULAR_P_CENTER_2_H
#include <CGAL/rectangular_p_center_2.h>
#endif // CGAL_RECTANGULAR_P_CENTER_2_H
#ifndef CGAL_POINT_GENERATORS_2_H
#include <CGAL/point_generators_2.h>
#endif // CGAL_POINT_GENERATORS_2_H
#ifndef CGAL_LEDA_REAL_H
#include <CGAL/leda_real.h>
#endif // CGAL_LEDA_REAL_H
#ifndef CGAL_COPY_N_H
#include <CGAL/copy_n.h>
#endif // CGAL_COPY_N_H
#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <cstdlib>

using std::vector;
using std::copy;
using std::back_inserter;
using std::ostream_iterator;
using std::transform;
using std::bind2nd;
using CGAL::Cartesian;
using CGAL::Random;
using CGAL::rectangular_p_center_2;
using CGAL::default_random;
using CGAL::Creator_uniform_2;
using CGAL::Random_points_in_square_2;
using CGAL::Istream_iterator;
using CGAL::Ostream_iterator;
using CGAL::set_pretty_mode;
using CGAL::cgalize;
using CGAL::BLUE;
using CGAL::RED;
using CGAL::ORANGE;

typedef double                            FT;
// typedef leda_real                      FT;
typedef Cartesian< FT >                   R;
typedef CGAL::Point_2< R >                Point;
typedef CGAL::Iso_rectangle_2< R >        Square_2;
typedef vector< Point >                   Point_cont;
typedef Creator_uniform_2< FT, Point >    Creator;
typedef Random_points_in_square_2< Point, Creator >
  Point_generator;
typedef Ostream_iterator< Point, leda_window >
  Window_stream_iterator_point;
typedef Ostream_iterator< Square_2, leda_window >
  Window_stream_iterator_square;
typedef ostream_iterator< Point >       Ostream_iterator_point;
typedef ostream_iterator< Square_2 >    Ostream_iterator_square;
typedef Istream_iterator< Point, leda_window>
  Istream_iterator_point;

#include <ctime>
static time_t Measure;
static long long int measure;
#define MEASURE(comm) \
Measure = clock(); \
  comm; \
  measure = (long long int)((float)(clock() - Measure) \
  * 1000 / CLOCKS_PER_SEC); \
    cout << "[time: " << measure << " msec]\n";
#define MEASURE_NO_OUTPUT(comm) \
Measure = clock(); \
  comm; \
  measure = (long long int)((float)(clock() - Measure) \
  * 1000 / CLOCKS_PER_SEC);

// function class to construct a box
// around a point p with radius r
template < class Point, class FT, class Box >
struct Build_box
: public CGAL_STD::binary_function< Point, FT, Box >
{
  Box
  operator()( const Point& p, const FT& r) const
  {
    return Box( Point( p.x() - r, p.y() - r),
                Point( p.x() + r, p.y() + r));
  }
};

int
main( int argc, char* argv[])
{
  int number_of_points;
  Point_cont input_points;
  Point_cont output_points;
  typedef Build_box< Point, FT, Square_2 >  Build_square;

  // init CGAL stuff:
  leda_window W;
  cgalize( W);
  W.init( -1.5, 1.5, -1.5);
  W.display();
  set_pretty_mode( cout);
  set_pretty_mode( cerr);
  Window_stream_iterator_point wout_p( W);
  Window_stream_iterator_square wout_s( W);
  Ostream_iterator_point cout_p( cout, "\n");
  Ostream_iterator_square cout_s( cout, "\n");

  if ( argc < 2 || (number_of_points = atoi(argv[1])) <= 0) {
    cout << "-- reading input point set\n"
         << "-- press middle mouse button for last point"
         << endl;
    W << BLUE;
    copy( Istream_iterator_point( W),
          Istream_iterator_point(),
          back_inserter( input_points));
  }
  else {
    int random_seed( default_random.get_int( 0, (1 << 31)));
    if ( argc >= 3)
      // get seed from command line
      random_seed = atoi(argv[2]);

    Random my_rnd( random_seed);
    cout << "***********************************************\n"
         << "PCENTER - test with " << number_of_points
         << " points\n  random seed is " << random_seed
         << "\n***********************************************"
         << endl;

    // generate point set:
    Point_generator gen( 1.0, my_rnd);
    CGAL::copy_n( gen,
                  number_of_points,
                  back_inserter( input_points));
  } // else

  // show point set:
  W.clear();
  W << BLUE;
  copy( input_points.begin(), input_points.end(), wout_p);

  FT result;
  if ( !input_points.empty()) {
    MEASURE(
      rectangular_p_center_2(
        input_points.begin(),
        input_points.end(),
        back_inserter( output_points),
        result,
        4);
      )
  } // if ( !input_points.empty())

#if defined(CGAL_PCENTER_TRACE)
  int number_of_piercing_points( output_points.size());
  cerr << "finished with diameter " << result
       << " and " << number_of_piercing_points
       << " points." << endl;
#endif

  // show center points
  W << RED;
  copy( output_points.begin(), output_points.end(), wout_p);
  copy( output_points.begin(), output_points.end(), cout_p);
  cout << endl;

  // ... and the corresponding squares:
  W << ORANGE;
  transform( output_points.begin(),
             output_points.end(),
             wout_s,
             bind2nd( Build_square(), result / FT( 2)));
  transform( output_points.begin(),
             output_points.end(),
             cout_s,
             bind2nd( Build_square(), result / FT( 2)));
  cout << endl;

  double x, y;
  while ( W.read_mouse( x, y) != MOUSE_BUTTON( 1)) {}
}
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

