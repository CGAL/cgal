// Copyright (c) 1998-2003  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>

#ifndef CGAL_ALL_FURTHEST_NEIGHBORS_2_H
#define CGAL_ALL_FURTHEST_NEIGHBORS_2_H 1

#include <CGAL/license/Polytope_distance_d.h>


#include <CGAL/Optimisation/assertions.h>
#include <CGAL/Cartesian_matrix.h>
#include <CGAL/Dynamic_matrix.h>
#include <CGAL/monotone_matrix_search.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <algorithm>
#include <boost/bind.hpp>

namespace CGAL {
template < class Operation, class RandomAccessIC >
class All_furthest_neighbor_matrix
: public Cartesian_matrix< Operation, RandomAccessIC, RandomAccessIC >
// represents the matrix used for computing
// all furthest neighbors of a convex polygon
{
public:
  typedef
    Cartesian_matrix< Operation, RandomAccessIC, RandomAccessIC > Base;

  typedef typename Base::Value Value;

  using Base::number_of_rows;

  All_furthest_neighbor_matrix(RandomAccessIC f, RandomAccessIC l)
  : Base(f, l, f, l)
  {}

  All_furthest_neighbor_matrix(RandomAccessIC f,
                               RandomAccessIC l,
                               const Operation& o)
  : Base(f, l, f, l, o)
  {}

  int number_of_columns() const { return 2 * number_of_rows() - 1; }

  Value
  operator()( int r, int c) const
  {
    CGAL_optimisation_precondition(r >= 0 && r < number_of_rows());
    CGAL_optimisation_precondition(c >= 0 && c < number_of_columns());
    if (c <= r)
      return Value(c - r);
    else if (c >= r + number_of_rows())
      return Value(0);
    else if (c >= number_of_rows())
      return Base::operator()(r, c - number_of_rows());
    else
      return Base::operator()(r, c);
  }
};


} //namespace CGAL
#include <iterator>
namespace CGAL {

template < class RandomAccessIC,
           class OutputIterator,
           class Traits,
           class IteratorCategory >
OutputIterator
all_furthest_neighbors_2( RandomAccessIC points_begin,
                          RandomAccessIC points_end,
                          OutputIterator o,
                          const Traits& t,
                          IteratorCategory)
{
  using std::vector;
  using std::transform;
  using std::modulus;

  typedef All_furthest_neighbor_matrix<
    typename Traits::Compute_squared_distance_2, RandomAccessIC >
  Afn_matrix;

 // check preconditions:
  int number_of_points(
                       static_cast<int>(iterator_distance( points_begin, points_end)));
  CGAL_optimisation_precondition( number_of_points > 0);
  CGAL_optimisation_expensive_precondition(
    is_convex_2( points_begin, points_end, t));

  // prepare random access container:
  vector< int > v;
  v.reserve( number_of_points);
  for (int i = 0; i < number_of_points; ++i)
    v.push_back( 0);

  // compute maxima:
  monotone_matrix_search(
    dynamic_matrix(
      Afn_matrix(points_begin,
                 points_end,
                 t.compute_squared_distance_2_object())),
    v.begin());

  // output result:
  return transform(v.begin(),
		   v.end(),
		   o,
		   boost::bind(modulus<int>(), _1, number_of_points));
} // all_furthest_neighbors_2( ... )


template < class RandomAccessIC, class OutputIterator, class Traits >
OutputIterator
all_furthest_neighbors_2( RandomAccessIC points_begin,
                          RandomAccessIC points_end,
                          OutputIterator o,
                          const Traits&
                          CGAL_optimisation_expensive_precondition_code(t),
                          std::random_access_iterator_tag)
{
  typedef All_furthest_neighbor_matrix<
    typename Traits::Compute_squared_distance_2, RandomAccessIC >
  Afn_matrix;

  // check preconditions:
  int number_of_points(
                       static_cast<int>(iterator_distance( points_begin, points_end)));
  CGAL_optimisation_precondition( number_of_points > 0);
  CGAL_optimisation_expensive_precondition(
    is_convex_2( points_begin, points_end, t));

  // compute maxima:
  monotone_matrix_search(
    dynamic_matrix(
      Afn_matrix( points_begin, points_end)),
    o);

  return o + number_of_points;
} // all_furthest_neighbors_2( ... )

template < class RandomAccessIC, class OutputIterator, class Traits >
inline
OutputIterator
all_furthest_neighbors_2( RandomAccessIC points_begin,
                          RandomAccessIC points_end,
                          OutputIterator o,
                          const Traits& t)
{
  typedef typename
    std::iterator_traits< OutputIterator >::iterator_category
  iterator_category;

  return all_furthest_neighbors_2(
    points_begin, points_end, o, t, iterator_category());
} // all_furthest_neighbors_2( ... )


template < class RandomAccessIC, class OutputIterator >
inline
OutputIterator
all_furthest_neighbors_2( RandomAccessIC points_begin,
                          RandomAccessIC points_end,
                          OutputIterator o)
{
  typedef typename std::iterator_traits< RandomAccessIC >::value_type::R R;
  return all_furthest_neighbors_2( points_begin, points_end, o, R());
} // all_furthest_neighbors_2( ... )

// backwards compatibility
template < class RandomAccessIC, class OutputIterator >
inline
OutputIterator
all_furthest_neighbors( RandomAccessIC points_begin,
                        RandomAccessIC points_end,
                        OutputIterator o)
{ return all_furthest_neighbors_2( points_begin, points_end, o); }

} //namespace CGAL

#endif // ! (CGAL_ALL_FURTHEST_NEIGHBORS_2_H)
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------
