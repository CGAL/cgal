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
// file          : include/CGAL/random_convex_set_2.h
// package       : $CGAL_Package: Generator 2.12 (28 Jul 1999) $
// source        : src/rcs/rcs.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// coordinator   : ETH Zurich (Bernd Gaertner <gaertner@inf.ethz.ch>)
//
// Random Convex Point Sets
// ============================================================================

#if ! (CGAL_RANDOM_CONVEX_SET_2_H)
#define CGAL_RANDOM_CONVEX_SET_2_H 1

#include <CGAL/basic.h>
#include <CGAL/copy_n.h>
#include <vector>
#include <algorithm>
#include <numeric>
#ifdef CGAL_REP_CLASS_DEFINED
#include <CGAL/Random_convex_set_traits_2.h>
#endif

CGAL_BEGIN_NAMESPACE
template < class OutputIterator, class Point_generator, class Traits >
OutputIterator
random_convex_set_2( int n,
                     OutputIterator o,
                     const Point_generator& pg,
                     Traits t)
{
  CGAL_precondition( n >= 3);

#ifndef CGAL_CFG_NO_NAMESPACE
  using std::vector;
  using std::back_inserter;
  using std::accumulate;
  using std::transform;
  using std::bind2nd;
  using std::sort;
  using std::partial_sum;
  using std::less;
  using std::max_element;
  using CGAL::copy_n;
#endif

  typedef typename Traits::Point_2         Point_2;
  typedef typename Traits::FT              FT;
  typedef vector< Point_2 >                Container;
  typedef typename Traits::Sum             Sum;
  typedef typename Traits::Scale           Scale;
  typedef typename Traits::Angle_less      Angle_less;
  typedef typename Traits::Max_coordinate  Max_coordinate;

  // GCC 2.8 and egcs-1.0.1 require these:
  // (does not accept s.l. Scale()( p, 1))
  Scale scale;
  Max_coordinate max_coordinate;
  Sum sum;

  // build random point set:
  Container points;
  points.reserve( n);
  CGAL::copy_n( pg, n, back_inserter( points));

  // compute centroid of points:
  Point_2 centroid(
    scale(
      accumulate( points.begin(), points.end(), t.origin(), Sum()),
      FT( 1) / FT( n)));

  // translate s.t. centroid == origin:
  transform(
    points.begin(),
    points.end(),
    points.begin(),
    bind2nd( Sum(), scale( centroid, FT( -1))));

  // sort them according to their direction's angle
  // w.r.t. the positive x-axis:
  sort( points.begin(), points.end(), Angle_less());

  // construct polygon:
  partial_sum(
    points.begin(), points.end(), points.begin(), Sum());

  // and compute its centroid:
  Point_2 new_centroid(
    scale(
      accumulate( points.begin(), points.end(), t.origin(), Sum()),
      FT( 1) / FT( n)));

  // translate s.t. centroids match:
  transform(
    points.begin(),
    points.end(),
    points.begin(),
    bind2nd( Sum(), sum( centroid,
                         scale( new_centroid, FT( -1)))));

  // compute maximal coordinate:
  FT maxcoord( max_coordinate(
    *max_element( points.begin(),
                  points.end(),
                  compose2_2( less< FT >(),
                                   Max_coordinate(),
                                   Max_coordinate()))));

  // and finally scale to fit into original grid:
  return transform(
    points.begin(),
    points.end(),
    o,
    bind2nd( Scale(), FT( pg.range()) / maxcoord));

} // random_convex_set_2( n, o, pg, t)
CGAL_END_NAMESPACE

#endif // ! (CGAL_RANDOM_CONVEX_SET_2_H)

// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

