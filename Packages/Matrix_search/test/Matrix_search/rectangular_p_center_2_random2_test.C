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
// file          : rectangular_p_center_2_random2_test.C
// chapter       : $CGAL_Chapter: Geometric Optimisation $
// package       : $CGAL_Package: Matrix_search $
// source        : pcenter.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// coordinator   : ETH Zurich (Bernd Gaertner <gaertner@inf.ethz.ch>)
//
// 2-4-Centering Axis-Parallel 2D-Rectangles - test program
// ============================================================================

#ifndef CGAL_CARTESIAN_H
#include <CGAL/Cartesian.h>
#endif // CGAL_CARTESIAN_H
#ifndef CGAL_POINT_2_H
#include <CGAL/Point_2.h>
#endif // CGAL_POINT_2_H
#ifndef CGAL_VECTOR_2_H
#include <CGAL/Vector_2.h>
#endif // CGAL_VECTOR_2_H
#ifndef CGAL_ISO_RECTANGLE_2_H
#include <CGAL/Iso_rectangle_2.h>
#endif // CGAL_ISO_RECTANGLE_2_H
#ifndef CGAL_POINT_GENERATORS_2_H
#include <CGAL/point_generators_2.h>
#endif // CGAL_POINT_GENERATORS_2_H
#ifndef CGAL_RECTANGULAR_P_CENTER_2_H
#include <CGAL/rectangular_p_center_2.h>
#endif // CGAL_RECTANGULAR_P_CENTER_2_H
#ifndef CGAL_RANDOM_H
#include <CGAL/Random.h>
#endif // CGAL_RANDOM_H
#ifndef CGAL_COPY_N_H
#include <CGAL/copy_n.h>
#endif // CGAL_COPY_N_H
#include <vector>
#include <functional>
#include <algorithm>
#ifdef OUTPUT
#include <iostream>
#endif
#include <cstdlib>

using std::vector;
using std::back_inserter;
using CGAL::Cartesian;
using CGAL::Creator_uniform_2;
using CGAL::Random_points_in_square_2;
using CGAL::Random;
using CGAL::default_random;
using CGAL::rectangular_p_center_2;
using CGAL::ORIGIN;

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

typedef double                             FT;
typedef Cartesian< FT >                    R;
typedef CGAL::Point_2< R >                 Point;
typedef CGAL::Vector_2< R >                Vector;
typedef CGAL::Iso_rectangle_2< R >         Square_2;
typedef Build_box< Point, FT, Square_2 >   Build_square;
typedef vector< Point >                    PCont;
typedef PCont::iterator                    Piter;
typedef vector< Square_2 >                 SCont;
typedef SCont::iterator                    Siter;
typedef Creator_uniform_2< FT, Point >     Creator;
typedef Random_points_in_square_2< Point, Creator >
                                             Point_generator;

// translate a range of points by v:
template < class ForwardIterator, class OutputIterator >
void
translate_it( ForwardIterator f,
              ForwardIterator l,
              OutputIterator o,
              const Vector& v)
{
  for ( ForwardIterator i( f); i != l; ++i)
    *o++ = *i + v;
}

int
main( int argc, char* argv[])
{
#ifdef OUTPUT
  CGAL::set_pretty_mode( cerr);
#endif

  int number_of_points;
  int random_seed;

  // handle command line arguments:
  if ( argc < 2 || (number_of_points = atoi(argv[1])) <= 0) {
    cerr << "usage: " << argv[0]
         << " num [random_seed]" << endl;
    exit(1);
  }
  if ( argc < 3) {

#ifdef OUTPUT
  cerr << "warning: no random seed specified\n"
         << "generating random seed" << endl;
#endif

    // generate random seed
    random_seed = default_random.get_int( 0, (1 << 30));
  }
  else
    random_seed = atoi(argv[2]);

  // define random source:
  Random rnd( random_seed);

#ifdef OUTPUT
  cerr << "random seed is " << random_seed << endl;
#endif

  PCont points;
  Vector t;
  Point_generator ptgen( 1, rnd);
  FT p_radius;

  // generate a random cluster of size number_of_points:
  CGAL::copy_n( ptgen, number_of_points, back_inserter( points));

  // and add the most extreme points:
  points.push_back( Point( 1, 1));
  points.push_back( Point( -1, -1));

#ifdef OUTPUT
  cerr << "** check two center **" << endl;
#endif // OUTPUT
  {
    // vectors to translate the clusters:
    Vector v1( 0,2);
    Vector v2( 0,-2);
    PCont pts, centers;
    translate_it( points.begin(), points.end(), back_inserter( pts), v1);
    translate_it( points.begin(), points.end(), back_inserter( pts), v2);
  
    rectangular_p_center_2( pts.begin(),
                            pts.end(),
                            back_inserter( centers),
                            p_radius,
                            2);
  
    assert( p_radius == 2);
    assert( centers[0] == ORIGIN + v1 ||
            centers[0] == ORIGIN + v2);
    assert( centers[1] == ORIGIN + v1 ||
            centers[1] == ORIGIN + v2);
  
  #ifdef OUTPUT
      cerr << "." << flush;
  #endif // OUTPUT
    }
  {
    // vectors to translate the clusters:
    Vector v1( 2,0);
    Vector v2( -2,0);
    PCont pts, centers;
    translate_it( points.begin(), points.end(), back_inserter( pts), v1);
    translate_it( points.begin(), points.end(), back_inserter( pts), v2);
  
    rectangular_p_center_2( pts.begin(),
                            pts.end(),
                            back_inserter( centers),
                            p_radius,
                            2);
  
    assert( p_radius == 2);
    assert( centers[0] == ORIGIN + v1 ||
            centers[0] == ORIGIN + v2);
    assert( centers[1] == ORIGIN + v1 ||
            centers[1] == ORIGIN + v2);
  
  #ifdef OUTPUT
      cerr << "." << flush;
  #endif // OUTPUT
    }
  {
    // vectors to translate the clusters:
    Vector v1( 2,2);
    Vector v2( -2,-2);
    PCont pts, centers;
    translate_it( points.begin(), points.end(), back_inserter( pts), v1);
    translate_it( points.begin(), points.end(), back_inserter( pts), v2);
  
    rectangular_p_center_2( pts.begin(),
                            pts.end(),
                            back_inserter( centers),
                            p_radius,
                            2);
  
    assert( p_radius == 2);
    assert( centers[0] == ORIGIN + v1 ||
            centers[0] == ORIGIN + v2);
    assert( centers[1] == ORIGIN + v1 ||
            centers[1] == ORIGIN + v2);
  
  #ifdef OUTPUT
      cerr << "." << flush;
  #endif // OUTPUT
    }
  {
    // vectors to translate the clusters:
    Vector v1( -2,2);
    Vector v2( 2,-2);
    PCont pts, centers;
    translate_it( points.begin(), points.end(), back_inserter( pts), v1);
    translate_it( points.begin(), points.end(), back_inserter( pts), v2);
  
    rectangular_p_center_2( pts.begin(),
                            pts.end(),
                            back_inserter( centers),
                            p_radius,
                            2);
  
    assert( p_radius == 2);
    assert( centers[0] == ORIGIN + v1 ||
            centers[0] == ORIGIN + v2);
    assert( centers[1] == ORIGIN + v1 ||
            centers[1] == ORIGIN + v2);
  
  #ifdef OUTPUT
      cerr << "." << flush;
  #endif // OUTPUT
    }

#ifdef OUTPUT
  cerr << "\n** check three center **" << endl;
#endif // OUTPUT
  {
    // vectors to translate the clusters:
    Vector v1( 0,4);
    Vector v2( 0,0);
    Vector v3( 0,-4);
    PCont pts, centers;
    translate_it( points.begin(), points.end(), back_inserter( pts), v1);
    translate_it( points.begin(), points.end(), back_inserter( pts), v2);
    translate_it( points.begin(), points.end(), back_inserter( pts), v3);
  
    rectangular_p_center_2( pts.begin(),
                            pts.end(),
                            back_inserter( centers),
                            p_radius,
                            3);
  
    assert( p_radius == 2);
    assert( centers[0] == ORIGIN + v1 ||
            centers[0] == ORIGIN + v2 ||
            centers[0] == ORIGIN + v3 );
    assert( centers[1] == ORIGIN + v1 ||
            centers[1] == ORIGIN + v2 ||
            centers[1] == ORIGIN + v3 );
  
  #ifdef OUTPUT
      cerr << "." << flush;
  #endif // OUTPUT
    }
  {
    // vectors to translate the clusters:
    Vector v1( -4,0);
    Vector v2( 0,0);
    Vector v3( 4,0);
    PCont pts, centers;
    translate_it( points.begin(), points.end(), back_inserter( pts), v1);
    translate_it( points.begin(), points.end(), back_inserter( pts), v2);
    translate_it( points.begin(), points.end(), back_inserter( pts), v3);
  
    rectangular_p_center_2( pts.begin(),
                            pts.end(),
                            back_inserter( centers),
                            p_radius,
                            3);
  
    assert( p_radius == 2);
    assert( centers[0] == ORIGIN + v1 ||
            centers[0] == ORIGIN + v2 ||
            centers[0] == ORIGIN + v3 );
    assert( centers[1] == ORIGIN + v1 ||
            centers[1] == ORIGIN + v2 ||
            centers[1] == ORIGIN + v3 );
  
  #ifdef OUTPUT
      cerr << "." << flush;
  #endif // OUTPUT
    }
  {
    // vectors to translate the clusters:
    Vector v1( -4,0);
    Vector v2( 0,2);
    Vector v3( 4,0);
    PCont pts, centers;
    translate_it( points.begin(), points.end(), back_inserter( pts), v1);
    translate_it( points.begin(), points.end(), back_inserter( pts), v2);
    translate_it( points.begin(), points.end(), back_inserter( pts), v3);
  
    rectangular_p_center_2( pts.begin(),
                            pts.end(),
                            back_inserter( centers),
                            p_radius,
                            3);
  
    assert( p_radius == 2);
    assert( centers[0] == ORIGIN + v1 ||
            centers[0] == ORIGIN + v2 ||
            centers[0] == ORIGIN + v3 );
    assert( centers[1] == ORIGIN + v1 ||
            centers[1] == ORIGIN + v2 ||
            centers[1] == ORIGIN + v3 );
  
  #ifdef OUTPUT
      cerr << "." << flush;
  #endif // OUTPUT
    }
  {
    // vectors to translate the clusters:
    Vector v1( -4,-4);
    Vector v2( 0,2);
    Vector v3( 4,0);
    PCont pts, centers;
    translate_it( points.begin(), points.end(), back_inserter( pts), v1);
    translate_it( points.begin(), points.end(), back_inserter( pts), v2);
    translate_it( points.begin(), points.end(), back_inserter( pts), v3);
  
    rectangular_p_center_2( pts.begin(),
                            pts.end(),
                            back_inserter( centers),
                            p_radius,
                            3);
  
    assert( p_radius == 2);
    assert( centers[0] == ORIGIN + v1 ||
            centers[0] == ORIGIN + v2 ||
            centers[0] == ORIGIN + v3 );
    assert( centers[1] == ORIGIN + v1 ||
            centers[1] == ORIGIN + v2 ||
            centers[1] == ORIGIN + v3 );
  
  #ifdef OUTPUT
      cerr << "." << flush;
  #endif // OUTPUT
    }
  {
    // vectors to translate the clusters:
    Vector v1( -4,-4);
    Vector v2( 2,0);
    Vector v3( 4,-4);
    PCont pts, centers;
    translate_it( points.begin(), points.end(), back_inserter( pts), v1);
    translate_it( points.begin(), points.end(), back_inserter( pts), v2);
    translate_it( points.begin(), points.end(), back_inserter( pts), v3);
  
    rectangular_p_center_2( pts.begin(),
                            pts.end(),
                            back_inserter( centers),
                            p_radius,
                            3);
  
    assert( p_radius == 2);
    assert( centers[0] == ORIGIN + v1 ||
            centers[0] == ORIGIN + v2 ||
            centers[0] == ORIGIN + v3 );
    assert( centers[1] == ORIGIN + v1 ||
            centers[1] == ORIGIN + v2 ||
            centers[1] == ORIGIN + v3 );
  
  #ifdef OUTPUT
      cerr << "." << flush;
  #endif // OUTPUT
    }

#ifdef OUTPUT
  cerr << "\n** check four center **" << endl;
#endif // OUTPUT
  {
    // vectors to translate the clusters:
    Vector v1( 0,-4);
    Vector v2( 0,-2);
    Vector v3( 0,2);
    Vector v4( 0,4);
    PCont pts, centers;
    translate_it( points.begin(), points.end(), back_inserter( pts), v1);
    translate_it( points.begin(), points.end(), back_inserter( pts), v2);
    translate_it( points.begin(), points.end(), back_inserter( pts), v3);
    translate_it( points.begin(), points.end(), back_inserter( pts), v4);
  
    rectangular_p_center_2( pts.begin(),
                            pts.end(),
                            back_inserter( centers),
                            p_radius,
                            4);
  
    assert( p_radius == 2);
    assert( centers[0] == ORIGIN + v1 ||
            centers[0] == ORIGIN + v2 ||
            centers[0] == ORIGIN + v3 ||
            centers[0] == ORIGIN + v4);
    assert( centers[1] == ORIGIN + v1 ||
            centers[1] == ORIGIN + v2 ||
            centers[1] == ORIGIN + v3 ||
            centers[1] == ORIGIN + v4);
  
  #ifdef OUTPUT
      cerr << "." << flush;
  #endif // OUTPUT
    }
  {
    // vectors to translate the clusters:
    Vector v1( -4,0);
    Vector v2( -2,0);
    Vector v3( 2,0);
    Vector v4( 4,0);
    PCont pts, centers;
    translate_it( points.begin(), points.end(), back_inserter( pts), v1);
    translate_it( points.begin(), points.end(), back_inserter( pts), v2);
    translate_it( points.begin(), points.end(), back_inserter( pts), v3);
    translate_it( points.begin(), points.end(), back_inserter( pts), v4);
  
    rectangular_p_center_2( pts.begin(),
                            pts.end(),
                            back_inserter( centers),
                            p_radius,
                            4);
  
    assert( p_radius == 2);
    assert( centers[0] == ORIGIN + v1 ||
            centers[0] == ORIGIN + v2 ||
            centers[0] == ORIGIN + v3 ||
            centers[0] == ORIGIN + v4);
    assert( centers[1] == ORIGIN + v1 ||
            centers[1] == ORIGIN + v2 ||
            centers[1] == ORIGIN + v3 ||
            centers[1] == ORIGIN + v4);
  
  #ifdef OUTPUT
      cerr << "." << flush;
  #endif // OUTPUT
    }
  {
    // vectors to translate the clusters:
    Vector v1( 0,-2);
    Vector v2( 0,2);
    Vector v3( -20,0);
    Vector v4( 20,1);
    PCont pts, centers;
    translate_it( points.begin(), points.end(), back_inserter( pts), v1);
    translate_it( points.begin(), points.end(), back_inserter( pts), v2);
    translate_it( points.begin(), points.end(), back_inserter( pts), v3);
    translate_it( points.begin(), points.end(), back_inserter( pts), v4);
  
    rectangular_p_center_2( pts.begin(),
                            pts.end(),
                            back_inserter( centers),
                            p_radius,
                            4);
  
    assert( p_radius == 2);
    assert( centers[0] == ORIGIN + v1 ||
            centers[0] == ORIGIN + v2 ||
            centers[0] == ORIGIN + v3 ||
            centers[0] == ORIGIN + v4);
    assert( centers[1] == ORIGIN + v1 ||
            centers[1] == ORIGIN + v2 ||
            centers[1] == ORIGIN + v3 ||
            centers[1] == ORIGIN + v4);
  
  #ifdef OUTPUT
      cerr << "." << flush;
  #endif // OUTPUT
    }
  {
    // vectors to translate the clusters:
    Vector v1( 0,-2);
    Vector v2( 0,2);
    Vector v3( -200,0);
    Vector v4( 200,1);
    PCont pts, centers;
    translate_it( points.begin(), points.end(), back_inserter( pts), v1);
    translate_it( points.begin(), points.end(), back_inserter( pts), v2);
    translate_it( points.begin(), points.end(), back_inserter( pts), v3);
    translate_it( points.begin(), points.end(), back_inserter( pts), v4);
  
    rectangular_p_center_2( pts.begin(),
                            pts.end(),
                            back_inserter( centers),
                            p_radius,
                            4);
  
    assert( p_radius == 2);
    assert( centers[0] == ORIGIN + v1 ||
            centers[0] == ORIGIN + v2 ||
            centers[0] == ORIGIN + v3 ||
            centers[0] == ORIGIN + v4);
    assert( centers[1] == ORIGIN + v1 ||
            centers[1] == ORIGIN + v2 ||
            centers[1] == ORIGIN + v3 ||
            centers[1] == ORIGIN + v4);
  
  #ifdef OUTPUT
      cerr << "." << flush;
  #endif // OUTPUT
    }
  {
    // vectors to translate the clusters:
    Vector v1( -2,0);
    Vector v2( 2,0);
    Vector v3( 0,-200);
    Vector v4( 1,200);
    PCont pts, centers;
    translate_it( points.begin(), points.end(), back_inserter( pts), v1);
    translate_it( points.begin(), points.end(), back_inserter( pts), v2);
    translate_it( points.begin(), points.end(), back_inserter( pts), v3);
    translate_it( points.begin(), points.end(), back_inserter( pts), v4);
  
    rectangular_p_center_2( pts.begin(),
                            pts.end(),
                            back_inserter( centers),
                            p_radius,
                            4);
  
    assert( p_radius == 2);
    assert( centers[0] == ORIGIN + v1 ||
            centers[0] == ORIGIN + v2 ||
            centers[0] == ORIGIN + v3 ||
            centers[0] == ORIGIN + v4);
    assert( centers[1] == ORIGIN + v1 ||
            centers[1] == ORIGIN + v2 ||
            centers[1] == ORIGIN + v3 ||
            centers[1] == ORIGIN + v4);
  
  #ifdef OUTPUT
      cerr << "." << flush;
  #endif // OUTPUT
    }

#ifdef OUTPUT
  cerr << "\n*** done ***" << endl;
#endif // OUTPUT

  return 0;
} 
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

