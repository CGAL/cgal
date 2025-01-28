// Copyright (c) 1998
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>

#ifndef CGAL_RANDOM_CONVEX_SET_2_H
#define CGAL_RANDOM_CONVEX_SET_2_H 1

#include <CGAL/basic.h>
#include <CGAL/algorithm.h>
#include <vector>
#include <algorithm>
#include <numeric>
#include <CGAL/Random_convex_set_traits_2.h>

namespace CGAL {

template < class OutputIterator, class Point_generator, class Traits >
OutputIterator
random_convex_set_2( std::size_t n,
                     OutputIterator o,
                     const Point_generator& pg,
                     Traits t)
{
  CGAL_precondition( n >= 3);

  using std::vector;
  using std::back_inserter;
  using std::accumulate;
  using std::transform;
  using std::sort;
  using std::partial_sum;
  using std::less;
  using std::max_element;
  using std::copy_n;

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
  std::copy_n( pg, n, back_inserter( points));

  // compute centroid of points:
  // Point_2 centroid = CGAL::centroid( points.begin(), points.end(), t );

  Point_2 centroid = t.origin();

  for(const Point_2& p : points){
    centroid = sum(centroid, p);
  }
  centroid = scale(centroid, FT(1)/FT(n));

  // translate s.t. centroid == origin:
  transform(
    points.begin(),
    points.end(),
    points.begin(),
    [&centroid, &sum, &scale](const Point_2& p) { return sum(p, scale(centroid, FT( -1))); });

  // sort them according to their direction's angle
  // w.r.t. the positive x-axis:
  sort( points.begin(), points.end(), Angle_less());

  // construct polygon:
  partial_sum(
    points.begin(), points.end(), points.begin(), Sum());

  // and compute its centroid:
  Point_2 new_centroid = t.origin();

  for(const Point_2& p : points){
    new_centroid = sum(new_centroid, p);
  }
  new_centroid = scale(new_centroid, FT(1)/FT(n));
  // translate s.t. centroids match:
  transform(
    points.begin(),
    points.end(),
    points.begin(),
    [&centroid, &new_centroid, &sum, &scale](const Point_2& p)
    {return sum(p, sum( centroid, scale(new_centroid, FT( -1)))); }
  );

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
    [&pg, &maxcoord, &scale](const Point_2& p){ return scale(p, FT( pg.range()) / maxcoord); });

} // random_convex_set_2( n, o, pg, t)

} //namespace CGAL

#endif // ! (CGAL_RANDOM_CONVEX_SET_2_H)
