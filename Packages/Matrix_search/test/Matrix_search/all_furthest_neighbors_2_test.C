#line 1513 "mon_search.aw"
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

#line 1517 "mon_search.aw"
#line 692 "afn.awi"
#line 549 "afn.awi"
#ifndef CGAL_CARTESIAN_H
#include <CGAL/Cartesian.h>
#endif // CGAL_CARTESIAN_H
#ifndef CGAL_POINT_2_H
#include <CGAL/Point_2.h>
#endif // CGAL_POINT_2_H
#ifndef CGAL_POLYGON_2_H
#include <CGAL/Polygon_2.h>
#endif // CGAL_POLYGON_2_H
#ifndef CGAL_POINT_GENERATORS_2_H
#include <CGAL/point_generators_2.h>
#endif // CGAL_POINT_GENERATORS_2_H
#ifndef CGAL_RANDOM_CONVEX_SET_2_H
#include <CGAL/random_convex_set_2.h>
#endif // CGAL_RANDOM_CONVEX_SET_2_H
#ifndef CGAL_ALL_FURTHEST_NEIGHBORS_2_H
#include <CGAL/all_furthest_neighbors_2.h>
#endif // CGAL_ALL_FURTHEST_NEIGHBORS_2_H
#include <vector>
#line 693 "afn.awi"

using std::equal;
#line 566 "afn.awi"
using std::vector;
using std::back_inserter;
using CGAL::Cartesian;
using CGAL::Polygon_traits_2;
using CGAL::Creator_uniform_2;
using CGAL::Random_points_in_square_2;
using CGAL::random_convex_set_2;
using CGAL::all_furthest_neighbors;
#line 696 "afn.awi"
using CGAL::squared_distance;
using CGAL::iterator_distance;

typedef double                                 FT;
#line 577 "afn.awi"
typedef Cartesian< FT >                              R;
typedef CGAL::Point_2< R >                           Point;
typedef Polygon_traits_2< R >                        P_traits;
typedef vector< Point >                              Point_cont;
typedef vector< int >                                Index_cont;
typedef CGAL::Polygon_2< P_traits, Point_cont >      Polygon;
typedef Creator_uniform_2< FT, Point >               Creator;
typedef Random_points_in_square_2< Point, Creator >  Point_generator;
#line 701 "afn.awi"

#line 737 "afn.awi"
#ifndef CGAL_SQUARED_DISTANCE_2_H
#include <CGAL/squared_distance_2.h>
#endif // CGAL_SQUARED_DISTANCE_2_H
#ifndef CGAL_CIRCULATOR_H
#include <CGAL/circulator.h>
#endif // CGAL_CIRCULATOR_H
#include <algorithm>

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
      if ( squared_distance( *i1, *i2) >
           squared_distance( *i1, *i))
        i = i2;
    } while ( ++i2 != e);
    *o++ = iterator_distance( b, i);
  } while ( ++i1 != e);
  return o;
} // afn_brute_force( b, e, o)
#line 703 "afn.awi"

int
main()
{
  int size [] = { 3, 5, 20, 101, 534 };
  for ( int i( 0); i < 5; ++i) {
    int number_of_points( size[i]);
    #line 607 "afn.awi"
    // generate random convex polygon:
    Polygon p;
    random_convex_set_2( number_of_points,
                         back_inserter( p),
                         Point_generator( 1));
#line 711 "afn.awi"
    #line 615 "afn.awi"
    // compute all furthest neighbors:
    Index_cont neighbors;
    all_furthest_neighbors(
      p.vertices_begin(),
      p.vertices_end(),
      back_inserter( neighbors));
#line 712 "afn.awi"
    #line 764 "afn.awi"
    // compute again brute force:
    Index_cont neighbors2;
    afn_brute_force(
      p.vertices_begin(),
      p.vertices_end(),
      back_inserter( neighbors2));
    
    // and compare both results:
    assert( equal( neighbors.begin(),
                   neighbors.end(),
                   neighbors2.begin()));
#line 713 "afn.awi"
  } // for ( int i( 0); i < 5; ++i)

#if !defined(CGAL_CFG_NO_ITERATOR_TRAITS) && \
  !defined(CGAL_CFG_NO_MEMBER_TEMPLATES)

  // try also once with a random-acccess iterator:
  int number_of_points( 222);
  #line 607 "afn.awi"
  // generate random convex polygon:
  Polygon p;
  random_convex_set_2( number_of_points,
                       back_inserter( p),
                       Point_generator( 1));
#line 721 "afn.awi"
  Index_cont neighbors( number_of_points);
  all_furthest_neighbors(
    p.vertices_begin(),
    p.vertices_end(),
    neighbors.begin());
#endif

  return 0;
} // int main()
#line 1518 "mon_search.aw"
#line 12 "code_formatting.awi"
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

