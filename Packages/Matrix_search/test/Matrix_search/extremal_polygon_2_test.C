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

#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>
#include <CGAL/extremal_polygon_2.h>
#include <function.h>
#include <vector.h>

#include <deque.h>

/*
#include <CGAL/leda_real.h>
#include <algo.h>
*/

// typedefs:
typedef double                            FT;
//typedef leda_real                         FT;
typedef CGAL_Cartesian< FT >              R;
typedef CGAL_Point_2< R >                 Point_2;
typedef CGAL_Polygon_traits_2< R >        P_traits;
typedef vector< Point_2 >                 Cont;
typedef CGAL_Polygon_2< P_traits, Cont >  Polygon_2;

// do random convex set generation always with double
// (coordinates get too long with exact computation)
typedef CGAL_Point_2< CGAL_Cartesian< double > >
  Point_2_double;
typedef vector< Point_2_double >          Cont_double;
typedef CGAL_Random_points_in_square_2<
  Point_2_double,
  CGAL_Creator_uniform_2< double, Point_2_double > >
Point_generator;
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
          CGAL_abs(
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
            CGAL_abs(
              (*i1).x() * ( (*i4).y() - (*i3).y()) +
              (*i4).x() * ( (*i3).y() - (*i1).y()) +
              (*i3).x() * ( (*i1).y() - (*i4).y())) +
            CGAL_abs(
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


/*
struct D2R : public unary_function< Point_2_double, Point_2 >
{
  Point_2
  operator()( const Point_2_double& p)
  { return Point_2( p.x(), p.y()); }
};
*/

int main() {
  // CGAL_set_pretty_mode( cout);

  int number_of_points [] = { 20, 51, 102, 500 };
  int k [] = { 3, 7, 12, 27 };
  int j;

  for ( int n( 0); n < 4; ++n) {
    /*
    cout << " : " << number_of_points[n] << endl;

    Cont_double p_d;
    CGAL_random_convex_set_2( number_of_points[n],
                              back_inserter( p_d),
                              Point_generator( 1));

    // build polygon:
    Polygon_2 p;
    transform( p_d.begin(),
               p_d.end(),
               back_inserter( p),
               D2R());
    */

    Polygon_2 p;
    CGAL_random_convex_set_2( number_of_points[n],
                              back_inserter( p),
                              Point_generator( 1));
    CGAL_assertion( p.is_convex());

    for ( j = 0; j < 4; ++j) {
      // maximum area:
      Cont k_gon;
      k_gon.reserve( k[j]);

      CGAL_maximum_area_inscribed_k_gon(
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
      Polygon_2 pp( k_gon.begin(), k_gon.end());
      cout << pp << endl;

      FT area_ms(
        CGAL_abs(
          (*(k_gon.begin())).x() *
          ( (*(k_gon.begin()+1)).y() - (*(k_gon.begin()+2)).y()) +
          (*(k_gon.begin()+1)).x() *
          ( (*(k_gon.begin()+2)).y() - (*(k_gon.begin())).y()) +
          (*(k_gon.begin()+2)).x() *
          ( (*(k_gon.begin())).y() - (*(k_gon.begin()+1)).y()))
        );
      cout << "area1 = " << area_ms << endl;

      cout << "k_gon2:\n";
      Polygon_2 ppp( k_gon2.begin(), k_gon2.end());
      cout << ppp << endl;

      FT area_bf(
        CGAL_abs(
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
      CGAL_maximum_perimeter_inscribed_k_gon(
        p.vertices_begin(),
        p.vertices_end(),
        k[j],
        back_inserter( k_gon));
    } // for ( j = 0; j < 4; ++j)

  } // for ( int n( 0); n < 4; ++n)

  return 0;
} // int main()
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

