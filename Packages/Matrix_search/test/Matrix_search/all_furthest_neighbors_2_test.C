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
// file          : all_furthest_neighbors_2_test.C
// chapter       : $CGAL_Chapter: Geometric Optimisation $
// package       : $CGAL_Package: Matrix_search $
// source        : mon_search.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// coordinator   : ETH Zurich (Bernd Gaertner <gaertner@inf.ethz.ch>)
//
// Test program: All Furthest Neighbors for a Convex Polygon
// ============================================================================

#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>
#include <CGAL/all_furthest_neighbors_2.h>
#include <vector.h>

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

#include <CGAL/squared_distance_2.h>
#include <CGAL/circulator.h>
#include <algo.h>

template < class RandomAccessIC,
           class OutputIterator >
OutputIterator
afn_brute_force( RandomAccessIC b,
                 RandomAccessIC e,
                 OutputIterator o)
{
  RandomAccessIC i1( b);
  do {
    RandomAccessIC i2( b);
    RandomAccessIC i( b);
    do {
      if ( CGAL_squared_distance( *i1, *i2) >
           CGAL_squared_distance( *i1, *i))
        i = i2;
    } while ( ++i2 != e);
    *o++ = CGAL_iterator_distance( b, i);
  } while ( ++i1 != e);
  return o;
} // afn_brute_force( b, e, o)

int
main()
{
  int size [] = { 3, 5, 20, 101, 534 };
  for ( int i( 0); i < 5; ++i) {
    int number_of_points( size[i]);
    // generate random convex polygon:
    Polygon_2 p;
    CGAL_random_convex_set_2( number_of_points,
                              back_inserter( p),
                              Point_generator( 1));
    // compute all furthest neighbors:
    Index_cont neighbors;
    CGAL_all_furthest_neighbors(
      p.vertices_begin(),
      p.vertices_end(),
      back_inserter( neighbors));
    // compute again brute force:
    Index_cont neighbors2;
    afn_brute_force(
      p.vertices_begin(),
      p.vertices_end(),
      back_inserter( neighbors2));
    
    // and compare both results:
    CGAL_optimisation_assertion( equal( neighbors.begin(),
                           neighbors.end(),
                           neighbors2.begin()));
  } // for ( int i( 0); i < 5; ++i)

#if !defined(CGAL_CFG_NO_ITERATOR_TRAITS) && \
  !defined(CGAL_CFG_NO_MEMBER_TEMPLATES)

  // try also once with a random-acccess iterator:
  int number_of_points( 222);
  // generate random convex polygon:
  Polygon_2 p;
  CGAL_random_convex_set_2( number_of_points,
                            back_inserter( p),
                            Point_generator( 1));
  Index_cont neighbors( number_of_points);
  CGAL_all_furthest_neighbors(
    p.vertices_begin(),
    p.vertices_end(),
    neighbors.begin());
#endif

  return 0;
} // int main()
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

