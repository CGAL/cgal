#line 1524 "mon_search.aw"
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
// file          : extremal_polygon_2_test.C
// chapter       : $CGAL_Chapter: Geometric Optimisation $
// package       : $CGAL_Package: Matrix_search $
// source        : mon_search.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// coordinator   : ETH Zurich (Bernd Gaertner <gaertner@inf.ethz.ch>)
//
// Test program: Compute extremal polygons of a convex polygon
// ============================================================================

#line 1528 "mon_search.aw"
#line 516 "testprog.awi"
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>
#include <CGAL/extremal_polygon_2.h>
#include <functional>
#include <vector>

/*
#include <CGAL/leda_real.h>
#include <algorithm>
*/

using std::vector;
using std::back_inserter;
using CGAL::Cartesian;
using CGAL::Polygon_traits_2;
using CGAL::Creator_uniform_2;
using CGAL::Random_points_in_square_2;
using CGAL::random_convex_set_2;
using CGAL::maximum_area_inscribed_k_gon;
using CGAL::maximum_perimeter_inscribed_k_gon;

// typedefs:
typedef double                             FT;
//typedef leda_real                        FT;
typedef Cartesian< FT >                    R;
typedef CGAL::Point_2< R >                 Point;
typedef Polygon_traits_2< R >              P_traits;
typedef vector< Point >                  Cont;
typedef CGAL::Polygon_2< P_traits, Cont >  Polygon;

// do random convex set generation always with double
// (coordinates get too long with exact computation)
typedef CGAL::Point_2< Cartesian< double > >        Point_double;
typedef vector< Point_double >                    Cont_double;
typedef Creator_uniform_2< double, Point_double > Creator;
typedef Random_points_in_square_2< Point_double, Creator >
  Point_generator;
#line 668 "testprog.awi"
template < class RandomAccessIC,
           class OutputIterator >
OutputIterator
brute_force_area_3( RandomAccessIC b,
                    RandomAccessIC e,
                    OutputIterator o)
{
  RandomAccessIC i1( b);
  RandomAccessIC m1, m2, m3;
  FT a_m( -1);
  do {
    RandomAccessIC i2( i1);
    do {
      RandomAccessIC i3( i2);
      do {
        FT a(
          abs(
            (*i1).x() * ( (*i2).y() - (*i3).y()) +
            (*i2).x() * ( (*i3).y() - (*i1).y()) +
            (*i3).x() * ( (*i1).y() - (*i2).y())));
        if ( a > a_m) {
          a_m = a;
          m1 = i1;
          m2 = i2;
          m3 = i3;
        }
      } while ( ++i3 != e);
    } while ( ++i2 != e);
  } while ( ++i1 != e);

  *o++ = *m1;
  *o++ = *m2;
  *o++ = *m3;
  return o;

} // brute_force_3( b, e, o)
#line 707 "testprog.awi"
template < class RandomAccessIC,
           class OutputIterator >
OutputIterator
brute_force_area_4( RandomAccessIC b,
                    RandomAccessIC e,
                    OutputIterator o)
{
  RandomAccessIC i1( b);
  RandomAccessIC m1, m2, m3, m4;
  FT a_m( -1);
  do {
    RandomAccessIC i2( i1);
    do {
      RandomAccessIC i3( i2);
      do {
        RandomAccessIC i4( i3);
        do {
          FT a(
            abs(
              (*i1).x() * ( (*i4).y() - (*i3).y()) +
              (*i4).x() * ( (*i3).y() - (*i1).y()) +
              (*i3).x() * ( (*i1).y() - (*i4).y())) +
            abs(
              (*i1).x() * ( (*i2).y() - (*i3).y()) +
              (*i2).x() * ( (*i3).y() - (*i1).y()) +
              (*i3).x() * ( (*i1).y() - (*i2).y())));
          if ( a > a_m) {
            a_m = a;
            m1 = i1;
            m2 = i2;
            m3 = i3;
            m4 = i4;
          }
        } while ( ++i4 != e);
      } while ( ++i3 != e);
    } while ( ++i2 != e);
  } while ( ++i1 != e);

  *o++ = *m1;
  *o++ = *m2;
  *o++ = *m3;
  *o++ = *m4;
  return o;

} // brute_force_4( b, e, o)

#line 556 "testprog.awi"

/*
struct D2R : public unary_function< Point_double, Point >
{
  Point
  operator()( const Point_double& p)
  { return Point( p.x(), p.y()); }
};
*/

int main() {
  // set_pretty_mode( cout);

  int number_of_points [] = { 20, 51, 102, 500 };
  int k [] = { 3, 7, 12, 27 };
  int j;

  for ( int n( 0); n < 4; ++n) {
    /*
    cout << " : " << number_of_points[n] << endl;

    Cont_double p_d;
    random_convex_set_2( number_of_points[n],
                         back_inserter( p_d),
                         Point_generator( 1));

    // build polygon:
    Polygon p;
    transform( p_d.begin(),
               p_d.end(),
               back_inserter( p),
               D2R());
    */

    Polygon p;
    random_convex_set_2( number_of_points[n],
                         back_inserter( p),
                         Point_generator( 1));
    CGAL_assertion( p.is_convex());

    for ( j = 0; j < 4; ++j) {
      // maximum area:
      Cont k_gon;
      k_gon.reserve( k[j]);

      maximum_area_inscribed_k_gon(
        p.vertices_begin(),
        p.vertices_end(),
        k[j],
        back_inserter( k_gon));

      /* TAKES TOO LONG:
      // check it:
      Cont k_gon2;
      brute_force_area_3(
        p.vertices_begin(),
        p.vertices_end(),
        back_inserter( k_gon2));

      cout << "k_gon:\n";
      Polygon pp( k_gon.begin(), k_gon.end());
      cout << pp << endl;

      FT area_ms(
        abs(
          (*(k_gon.begin())).x() *
          ( (*(k_gon.begin()+1)).y() - (*(k_gon.begin()+2)).y()) +
          (*(k_gon.begin()+1)).x() *
          ( (*(k_gon.begin()+2)).y() - (*(k_gon.begin())).y()) +
          (*(k_gon.begin()+2)).x() *
          ( (*(k_gon.begin())).y() - (*(k_gon.begin()+1)).y()))
        );
      cout << "area1 = " << area_ms << endl;

      cout << "k_gon2:\n";
      Polygon ppp( k_gon2.begin(), k_gon2.end());
      cout << ppp << endl;

      FT area_bf(
        abs(
          (*(k_gon2.begin())).x() *
          ( (*(k_gon2.begin()+1)).y() - (*(k_gon2.begin()+2)).y()) +
          (*(k_gon2.begin()+1)).x() *
          ( (*(k_gon2.begin()+2)).y() - (*(k_gon2.begin())).y()) +
          (*(k_gon2.begin()+2)).x() *
          ( (*(k_gon2.begin())).y() - (*(k_gon2.begin()+1)).y())));

      cout << "area2 = " << area_bf << endl;

      CGAL_assertion( area_bf == area_ms);
      */

    } // for ( j = 0; j < 4; ++j)

    for ( j = 0; j < 4; ++j) {
      // maximum perimeter:
      Cont k_gon;
      k_gon.reserve( k[j]);
      maximum_perimeter_inscribed_k_gon(
        p.vertices_begin(),
        p.vertices_end(),
        k[j],
        back_inserter( k_gon));
    } // for ( j = 0; j < 4; ++j)

  } // for ( int n( 0); n < 4; ++n)

  return 0;
} // int main()
#line 1529 "mon_search.aw"
#line 12 "code_formatting.awi"
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

