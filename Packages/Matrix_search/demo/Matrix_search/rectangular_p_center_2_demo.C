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
// 2-4-Centering Axis-Parallel 2D-Rectangles - demo program
// ============================================================================

#include <CGAL/Cartesian.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/IO/leda_window.h>
#include <CGAL/IO/Ostream_iterator.h>
#include <CGAL/IO/Istream_iterator.h>
#include <CGAL/rectangular_p_center_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/copy_n.h>
#include <CGAL/leda_real.h>
#include <iostream.h>
#include <vector.h>
#include <algo.h>
#include <function.h>
#include <stdlib.h>

typedef double                     FT;
// typedef leda_real                  FT;
typedef CGAL_Cartesian< FT >       R;
typedef CGAL_Point_2< R >          Point_2;
typedef CGAL_Iso_rectangle_2< R >  Square_2;
typedef vector< Point_2 >          Point_cont;
typedef CGAL_Random_points_in_square_2<
  Point_2,
  CGAL_Creator_uniform_2< FT, Point_2 > >
Point_generator;
typedef CGAL_Ostream_iterator< Point_2, leda_window >
  Window_stream_iterator_point;
typedef CGAL_Ostream_iterator< Square_2, leda_window >
  Window_stream_iterator_square;
typedef ostream_iterator< Point_2 >    Ostream_iterator_point;
typedef ostream_iterator< Square_2 >   Ostream_iterator_square;
typedef CGAL_Istream_iterator< Point_2, leda_window>
  Istream_iterator_point;

#include <time.h>
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
: public binary_function< Point, FT, Box >
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
  typedef Build_box< Point_2, FT, Square_2 >  Build_square;

  // init CGAL stuff:
  leda_window W;
  CGAL_cgalize( W);
  W.init( -1.5, 1.5, -1.5);
  W.display();
  CGAL_set_pretty_mode( cout);
  CGAL_set_pretty_mode( cerr);
  Window_stream_iterator_point wout_p( W);
  Window_stream_iterator_square wout_s( W);
  Ostream_iterator_point cout_p( cout, "\n");
  Ostream_iterator_square cout_s( cout, "\n");

  if ( argc < 2 || (number_of_points = atoi(argv[1])) <= 0) {
    cout << "-- reading input point set\n"
         << "-- press middle mouse button for last point"
         << endl;
    W << CGAL_BLUE;
    copy( Istream_iterator_point( W),
          Istream_iterator_point(),
          back_inserter( input_points));
  }
  else {
    int random_seed( CGAL_random.get_int( 0, (1 << 31)));
    if ( argc >= 3)
      // get seed from command line
      random_seed = atoi(argv[2]);

    CGAL_Random my_rnd( random_seed);
    cout << "***********************************************\n"
         << "PCENTER - test with " << number_of_points
         << " points\n  random seed is " << random_seed
         << "\n***********************************************"
         << endl;

    // generate point set:
    Point_generator gen( 1.0, my_rnd);
    CGAL_copy_n( gen,
                 number_of_points,
                 back_inserter( input_points));
  } // else

  // show point set:
  W.clear();
  W << CGAL_BLUE;
  copy( input_points.begin(), input_points.end(), wout_p);

  FT result;
  if ( !input_points.empty()) {
    MEASURE(
      CGAL_rectangular_p_center_2(
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
  W << CGAL_RED;
  copy( output_points.begin(), output_points.end(), wout_p);
  copy( output_points.begin(), output_points.end(), cout_p);
  cout << endl;

  // ... and the corresponding squares:
  W << CGAL_ORANGE;
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

