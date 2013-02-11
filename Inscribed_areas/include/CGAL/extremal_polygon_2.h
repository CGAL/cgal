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
// 
//
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>

#ifndef CGAL_EXTREMAL_POLYGON_2_H
#define CGAL_EXTREMAL_POLYGON_2_H 1

#include <CGAL/Optimisation/assertions.h>
#include <CGAL/monotone_matrix_search.h>
#include <CGAL/Dynamic_matrix.h>
#include <CGAL/Transform_iterator.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <vector>
#include <functional>
#include <algorithm>
#include <CGAL/Extremal_polygon_traits_2.h>

namespace CGAL {

//!!! This will eventually be integrated into function_objects.h
template < class Array, class Index, class Element >
struct Index_operator
: public std::binary_function< Array, Index, Element >
{

  Element&
  operator()( Array& a, const Index& i) const
  { return a[i]; }

  const Element&
  operator()( const Array& a, const Index& i) const
  { return a[i]; }
};

template < class RandomAccessIC_object_,
           class RandomAccessIC_value_,
           class Operation_ >
// This class describes the kind of matrices used for the
// computation of extremal polygons.
//
// RandomAccessIC_object is a random access iterator or circulator
//   with value type Object
// RandomAccessIC_value is a random access iterator or circulator
//   with value type Value
// Operation is an adatable binary function:
//   Object x Object -> Value
//
// objects can be constructed using the helper function
// extremal_polygon_matrix.
//
class Extremal_polygon_matrix {
public:
  typedef RandomAccessIC_object_ RandomAccessIC_object;
  typedef RandomAccessIC_value_  RandomAccessIC_value;
  typedef Operation_             Operation;

  typedef typename
    std::iterator_traits< RandomAccessIC_object >::value_type
  Object;
  typedef typename
    std::iterator_traits< RandomAccessIC_value >::value_type
  Value;

  Extremal_polygon_matrix(
    RandomAccessIC_object begin_row,
    RandomAccessIC_object end_row,
    RandomAccessIC_object begin_col,
    RandomAccessIC_object end_col,
    RandomAccessIC_value  begin_value,
    RandomAccessIC_value  CGAL_optimisation_precondition_code(end_value),
    const Operation&      o)
  // initialization with two ranges [begin_row, end_row) and
  // [begin_col, end_col) of Objects, a range [begin_value, end_value)
  // of Values and an Operation o.
  //
  // an entry (r, c) of this matrix is then defined as:
  //   begin_value[c] + op( begin_row[r], begin_col[c]).
  //
  : op( o),
    begin_row_( begin_row),
    begin_col_( begin_col),
    begin_value_( begin_value),
    n_rows( static_cast<int>(iterator_distance( begin_row, end_row))),
    n_cols( static_cast<int>(iterator_distance( begin_col, end_col)))
  {
    CGAL_optimisation_precondition(
      iterator_distance( begin_value, end_value) == n_cols);
    CGAL_optimisation_assertion( n_rows > 0 && n_cols > 0);
  }

  int
  number_of_rows() const
  { return n_rows; }

  int
  number_of_columns() const
  { return n_cols; }

  Value
  operator()( int r, int c) const
  {
    CGAL_optimisation_precondition( r >= 0 && r < n_rows);
    CGAL_optimisation_precondition( c >= 0 && c < n_cols);
    return begin_value_[c] + op( begin_row_[r], begin_col_[c]);
  }

private:
  Operation              op;
  RandomAccessIC_object  begin_row_;
  RandomAccessIC_object  begin_col_;
  RandomAccessIC_value   begin_value_;
  int                    n_rows;
  int                    n_cols;
};

template < class RandomAccessIC_object,
           class RandomAccessIC_value,
           class Operation >
inline
Extremal_polygon_matrix< RandomAccessIC_object,
                         RandomAccessIC_value,
                         Operation >
extremal_polygon_matrix(
  RandomAccessIC_object begin_row,
  RandomAccessIC_object end_row,
  RandomAccessIC_object begin_col,
  RandomAccessIC_object end_col,
  RandomAccessIC_value  begin_value,
  RandomAccessIC_value  end_value,
  const Operation&      o)
{
  return Extremal_polygon_matrix< RandomAccessIC_object,
                                  RandomAccessIC_value,
                                  Operation >
  ( begin_row, end_row,
    begin_col, end_col,
    begin_value, end_value,
    o);
}


template < class RandomAccessIC, class Outputiterator, class Traits >
Outputiterator
CGAL_maximum_inscribed_rooted_k_gon_2(
  RandomAccessIC points_begin,
  RandomAccessIC points_end,
  int k,
  typename Traits::FT& max_area,
  Outputiterator o,
  const Traits& t)
//
// preconditions:
// --------------
//  * Traits fulfills the requirements for an extremal polygon
//    traits class
//  * the range [points_begin, points_end) of size n > 0
//    describes the vertices of a convex polygon $P$
//    enumerated clock- or counterclockwise
// n > k,
//  * k >= t.min_k()
//  * value_type of RandomAccessIC is Traits::Point_2
//  * OutputIterator accepts Traits::Point_2 as value_type 
//
// functionality:
// --------------
// computes maximum (as specified by t) inscribed k-gon $P_k$
// of the polygon $P$,
// that is rooted at points_begin[0],
// sets max_area to its associated value (as specified by t)
// writes the indices (relative to points_begin)
// of $P_k$'s vertices to o and
// returns the past-the-end iterator of that sequence.
{
  // check preconditions:
  CGAL_optimisation_precondition( k >= t.min_k());
  int number_of_points(
                       static_cast<int>(iterator_distance( points_begin, 
                                                           points_end)));
  CGAL_optimisation_precondition( number_of_points > k);

  typedef std::vector< int > Index_cont;

  if ( k == t.min_k())
    // compute min_k gon:
    return t.compute_min_k_gon(
      points_begin, points_end, max_area, o);

  // current i-gon (i = 2/3...k)
  Index_cont gon( k + 1);

  // compute initial min_k-gon:
  int i( t.min_k());
  t.compute_min_k_gon(
    points_begin, points_end, max_area, gon.rbegin() + k + 1 - i);
  
  for (;;) {
    CGAL_optimisation_assertion( gon[0] == 0);
    gon[i] = number_of_points - 1;
    if ( ++i >= k)
      break;
    CGAL_maximum_inscribed_rooted_k_gon_2(
      points_begin,
      points_end,
      0,
      gon.begin(),
      gon.begin() + i - 1,
      gon.begin() + 1,
      gon.begin() + i,
      max_area,
      gon.rbegin() + k + 1 - i,
      t);
  } // for (;;)
  
  return CGAL_maximum_inscribed_rooted_k_gon_2(
    points_begin,
    points_end,
    0,
    gon.begin(),
    gon.begin() + k - 1,
    gon.begin() + 1,
    gon.begin() + k,
    max_area,
    o,
    t);

} // CGAL_maximum_inscribed_rooted_k_gon_2( ... )
template < class RandomAccessIC_point,
           class RandomAccessIC_int,
           class OutputIterator,
           class Traits >
OutputIterator
CGAL_maximum_inscribed_rooted_k_gon_2(
  RandomAccessIC_point points_begin,
  RandomAccessIC_point points_end,
  int root,
  RandomAccessIC_int left_c_begin,
  RandomAccessIC_int CGAL_optimisation_precondition_code(left_c_end),
  RandomAccessIC_int right_c_begin,
  RandomAccessIC_int right_c_end,
  typename Traits::FT& max_area,
  OutputIterator o,
  const Traits& t)
//
// preconditions:
// --------------
//  * Traits fulfills the requirements for an extremal polygon
//    traits class
//  * the range [points_begin, points_end) of size n > 0
//    describes the vertices of a convex polygon $P$
//    enumerated clock- or counterclockwise
//  * value_type of RandomAccessIC_point is Traits::Point
//  * value_type of RandomAccessIC_int is int
//  * OutputIterator accepts int as value type
//  * length := right_c_end - right_c_begin == left_c_end - left_c_begin
//    >= t.min_k() - 1 (the root is already fixed)
//  * [left_c_begin, left_c_end) resp. [right_c_begin, right_c_end)
//    describe two subpolygons of $P$ by giving the indices of its
//    vertices relative to points_begin and for any 0 <= i < length:
//    left_c_begin[i] <= right_c_begin[i]
//  * for any 0 <= i < length: o + i must not be contained in
//    the range [right_c_begin, right_c_begin + length - i - 2].
//    (NOT checked!)
//
// functionality:
// --------------
// computes maximum (as specified by t) inscribed k-gon $P_k$
// of the polygon $P$,
// that is rooted at points_begin[left_c_begin[0]]
// such that for any 0 <= i < length:
//    left_c_begin[i] <= vertex i of $P_k$ <= right_c_begin[i],
// sets max_area to its associated value (as specified by t),
// writes the indices (relative to points_begin)
// of $P_k$'s vertices to o and
// returns the past-the-end iterator of that sequence.
{
  using std::max_element;

  // counter :)
  int i;

  // compute size of ranges:
  int number_of_points = static_cast<int>(iterator_distance( points_begin, 
                                                             points_end));
  int size_of_gon = static_cast<int>(iterator_distance( right_c_begin,
                                                        right_c_end));

  // check preconditions:
  CGAL_optimisation_precondition( number_of_points > t.min_k());
  CGAL_optimisation_precondition( size_of_gon >= t.min_k() - 1);
  CGAL_optimisation_precondition(
    iterator_distance( left_c_begin, left_c_end) ==
    iterator_distance( right_c_begin, right_c_end));
  CGAL_optimisation_precondition( left_c_begin[0] >= 0);
  CGAL_optimisation_precondition( right_c_begin[0] >= 0);
  CGAL_optimisation_precondition(
    left_c_begin[size_of_gon-1] < number_of_points);
  CGAL_optimisation_precondition(
    right_c_begin[size_of_gon-1] < number_of_points);
  CGAL_optimisation_expensive_precondition_code(
    for ( i = 0; i < size_of_gon; ++i) {
      CGAL_optimisation_expensive_precondition( left_c_begin[i] >= 0);
      CGAL_optimisation_expensive_precondition( right_c_begin[i] >= 0);
      CGAL_optimisation_expensive_precondition(
        left_c_begin[i] < number_of_points);
      CGAL_optimisation_expensive_precondition(
        right_c_begin[i] < number_of_points);
      CGAL_optimisation_expensive_precondition(
        left_c_begin[i] <= right_c_begin[i]);
    })

  typedef typename Traits::FT               FT;
  typedef std::vector< FT >                 FT_cont;
  typedef std::vector< int >                Index_cont;
  typedef typename Traits::Operation        Operation;
  //!!! static ???
  // area container:
  FT_cont area( number_of_points);
  
  // last vertex container:
  Index_cont last_vertex( number_of_points);
  
  // matrix operation:
  Operation op( t.operation( points_begin[root]));
  // initialize area and last vertex containers:
  for ( i = left_c_begin[0]; i <= right_c_begin[0]; ++i) {
    area[i] = t.init( points_begin[i], points_begin[root]);
    last_vertex[i] = root;
  }
  
  
  for ( i = 1; i < size_of_gon; ++i) {
  
    monotone_matrix_search(
      dynamic_matrix(
        extremal_polygon_matrix(
          points_begin + left_c_begin[i],
          points_begin + right_c_begin[i] + 1,
          points_begin + left_c_begin[i-1],
          points_begin + right_c_begin[i-1] + 1,
          area.begin() + left_c_begin[i-1],
          area.begin() + right_c_begin[i-1] + 1,
          op)),
          last_vertex.begin() + left_c_begin[i]);
  
    // compute new area values and adjust last_vertex values
    // (they are relative to left_c_begin[i-1] now)
    int j;
    for ( j = right_c_begin[i]; j >= left_c_begin[i]; --j) {
      last_vertex[j] += left_c_begin[i-1];
      area[j] = area[last_vertex[j]] +
        op( points_begin[j], points_begin[last_vertex[j]]);
    }
  
  } // for ( i = 1; i < size_of_gon; ++i)
  
  // find maximum in last range:
  int maxi =
    static_cast<int>(iterator_distance(
      area.begin(),
      max_element( area.begin() + left_c_begin[size_of_gon - 1],
                   area.begin() + right_c_begin[size_of_gon - 1] + 1)));
  // set max_area:
  max_area = area[maxi];
  
  // construct gon:
  *o++ = maxi;
  maxi = last_vertex[maxi];
  for ( i = size_of_gon - 1; i > 0; --i) {
    // We must not place the "*o++ = maxi" here,
    // since o might be the same as left_c_begin + i ...
    if ( maxi != right_c_begin[i-1]) {
      *o++ = maxi;
      maxi = last_vertex[maxi];
    }
    else {
      *o++ = maxi;
      maxi = right_c_begin[i-2];
    }
  } // for ( i = size_of_gon - 1; i > 0; --i)
  
  *o++ = root;
  return o;
  

} // CGAL_maximum_inscribed_rooted_k_gon_2( p, k, result)


template < class RandomAccessIC,
           class OutputIterator,
           class Traits >
inline
OutputIterator
extremal_polygon_2(
  RandomAccessIC points_begin,
  RandomAccessIC points_end,
  int k,
  OutputIterator o,
  const Traits& t)
//
// preconditions:
// --------------
//  * Traits fulfills the requirements for an extremal polygon
//    traits class
//  * the range [points_begin, points_end) of size n > 0
//    describes the vertices of a convex polygon $P$
//    enumerated clock- or counterclockwise
//  * k >= t.min_k()
//  * value_type of RandomAccessIC is Traits::Point_2
//  * OutputIterator accepts Traits::Point_2 as value_type 
//
// functionality:
// --------------
// computes maximum (as specified by t) inscribed k-gon $P_k$
// of the polygon $P$,
// writes the indices (relative to points_begin)
// of $P_k$'s vertices to o and
// returns the past-the-end iterator of that sequence.
{
  using std::bind1st;

  // check preconditions:
  CGAL_optimisation_precondition_code(
    int number_of_points(
                         static_cast<int>(iterator_distance( points_begin, 
                                                             points_end)));)
  CGAL_optimisation_precondition( number_of_points >= t.min_k());
  CGAL_optimisation_expensive_precondition(
    is_convex_2( points_begin, points_end, t));

  typedef typename Traits::Point_2 Point_2;
  return CGAL_maximum_inscribed_k_gon_2(
    points_begin,
    points_end,
    k,
    transform_iterator(
      o,
      bind1st(
        Index_operator< RandomAccessIC, int, Point_2 >(),
        points_begin)),
    t);
}

// backwards compatibility
template < class RandomAccessIC,
           class OutputIterator,
           class Traits >
inline
OutputIterator
extremal_polygon(
  RandomAccessIC points_begin,
  RandomAccessIC points_end,
  int k,
  OutputIterator o,
  const Traits& t)
{ return extremal_polygon_2(points_begin, points_end, k, o, t); }
template < class RandomAccessIC,
           class OutputIterator,
           class Traits >
OutputIterator
CGAL_maximum_inscribed_k_gon_2(
  RandomAccessIC points_begin,
  RandomAccessIC points_end,
  int k,
  OutputIterator o,
  const Traits& t)
//
// preconditions:
// --------------
//  * Traits fulfills the requirements for an extremal polygon
//    traits class
//  * the range [points_begin, points_end) of size n > 0
//    describes the vertices of a convex polygon $P$
//    enumerated clock- or counterclockwise
//  * k >= t.min_k()
//  * value_type of RandomAccessIC is Traits::Point_2
//  * OutputIterator accepts Traits::Point_2 as value_type 
//
// functionality:
// --------------
// computes maximum (as specified by t) inscribed k-gon $P_k$
// of the polygon $P$,
// writes the indices (relative to points_begin)
// of $P_k$'s vertices to o and
// returns the past-the-end iterator of that sequence.
{
  // check preconditions:
  CGAL_optimisation_precondition( k >= t.min_k());
  int number_of_points(
                       static_cast<int>(iterator_distance( points_begin, 
                                                           points_end)));
  CGAL_optimisation_precondition( number_of_points > 0);

  using std::copy;

  typedef typename Traits::FT   FT;
  typedef std::vector< int >    Index_cont;

  if ( number_of_points <= k) {
    for ( int j( k - 1); j >= 0; --j)
      *o++ = (std::min)( j, number_of_points - 1);
    return o;
  }
  // compute k-gon rooted at points_begin[0]
  Index_cont P_0( k + 1);
  FT area_0;
  CGAL_maximum_inscribed_rooted_k_gon_2(
    points_begin,
    points_end,
    k,
    area_0,
    P_0.rbegin() + 1,
    t);
  P_0[k] = number_of_points - 1;
  CGAL_optimisation_assertion( P_0[0] == 0);
  // compute k-gon rooted at points_begin[P_0[1]]
  Index_cont P_1( k);
  FT area_1;
  
  CGAL_maximum_inscribed_rooted_k_gon_2(
    points_begin,
    points_end,
    P_0[1],
    P_0.begin() + 1,
    P_0.begin() + k,
    P_0.begin() + 2,
    P_0.begin() + k + 1,
    area_1,
    P_1.rbegin(),
    t);
  
  CGAL_optimisation_assertion( P_1[0] == P_0[1]);
  
  
  // start recursive computation:
  FT area_r( 0);
  Index_cont P_r( k);
  if ( P_0[1] - P_0[0] > 1) {
    CGAL_maximum_inscribed_k_gon_2(
      points_begin,
      points_end,
      P_0[0] + 1,
      P_0[1] - 1,
      P_0.begin() + 1,
      P_0.begin() + k,
      P_0.begin() + 2,
      P_0.begin() + k + 1,
      k,
      area_r,
      P_r.rbegin(),
      t);
  }
  
  if ( area_r > area_0)
    if ( area_r > area_1)
      // recursive is maximum
      copy( P_r.begin(), P_r.end(), o);
    else
      // P_1 is maximum
      copy( P_1.begin(), P_1.end(), o);
  else if ( area_0 > area_1)
    // P_0 is maximum
    copy( P_0.begin(), P_0.begin() + k, o);
  else
    // P_1 is maximum
    copy( P_1.begin(), P_1.end(), o);

  return o;
} // CGAL_maximum_inscribed_k_gon_2( ... )
template < class RandomAccessIC_point,
           class RandomAccessIC_int,
           class OutputIterator,
           class Traits >
OutputIterator
CGAL_maximum_inscribed_k_gon_2(
  RandomAccessIC_point points_begin,
  RandomAccessIC_point points_end,
  int left_index,
  int right_index,
  RandomAccessIC_int left_c_begin,
  RandomAccessIC_int left_c_end,
  RandomAccessIC_int right_c_begin,
  RandomAccessIC_int right_c_end,
  int k,
  typename Traits::FT& max_area,
  OutputIterator o,
  const Traits& t)
//
// preconditions:
// --------------
//  * Traits fulfills the requirements for an extremal polygon
//    traits class
//  * the range [points_begin, points_end) of size n > 0
//    describes the vertices of a convex polygon $P$
//    enumerated clock- or counterclockwise
//  * value_type of RandomAccessIC_point is Traits::Point
//  * value_type of RandomAccessIC_int is int
//  * OutputIterator accepts int as value type
//  * 0 <= left_index <= right_index < |points_end - points_begin|
//  * |left_c_end - left_c_begin| == |right_c_end - right_c_begin| == k - 1
//  * [left_c_begin, left_c_end) resp. [right_c_begin, right_c_end)
//    describe two subpolygons $P_l$ resp $P_r$ of $P$ by giving
//    the indices of its vertices relative to points_begin and
//    for any 0 <= i < k - 1:
//      left_c_begin[i] <= right_c_begin[i]
//  * k >= t.min_k()
//
// functionality:
// --------------
// computes maximum (as specified by t) inscribed k-gon $P_k$
// of the polygon $P$,
//  * that is rooted at one of the vertices [points_begin[left_index],
//    points_begin[right_index]] and
//  * interleaves with both $P_l$ and $P_r$,
// sets max_area to its associated value (as specified by t),
// writes the indices (relative to points_begin)
// of $P_k$'s vertices to o and
// returns the past-the-end iterator of that sequence.
{
  // typedefs
  typedef typename Traits::FT               FT;
  typedef std::vector< int >        Index_cont;

  using std::copy;

  // check preconditions:
  CGAL_optimisation_precondition( k >= t.min_k());
  CGAL_optimisation_precondition( left_index <= right_index);
  CGAL_optimisation_precondition( left_index >= 0);
  CGAL_optimisation_precondition( right_index >= 0);
  CGAL_optimisation_precondition_code(
    int number_of_points(
                         static_cast<int>(iterator_distance( points_begin, 
                                                             points_end)));)
  CGAL_optimisation_precondition( left_index < number_of_points);
  CGAL_optimisation_precondition( right_index < number_of_points);
  CGAL_optimisation_precondition(
    iterator_distance( left_c_begin, left_c_end) == k - 1);
  CGAL_optimisation_precondition(
    iterator_distance( right_c_begin, right_c_end) == k - 1);
  CGAL_optimisation_expensive_precondition_code(
    for ( int i( 0); i < k - 1; ++i) {
      CGAL_optimisation_expensive_precondition( left_c_begin[i] >= 0);
      CGAL_optimisation_expensive_precondition( right_c_begin[i] >= 0);
      CGAL_optimisation_expensive_precondition(
        left_c_begin[i] < number_of_points);
      CGAL_optimisation_expensive_precondition(
        right_c_begin[i] < number_of_points);
      CGAL_optimisation_expensive_precondition(
        left_c_begin[i] <= right_c_begin[i]);
    })


  int middle_index( (left_index + right_index) >> 1);
  Index_cont P_m( k);
  FT area_middle;
  CGAL_maximum_inscribed_rooted_k_gon_2(
    points_begin,
    points_end,
    middle_index,
    left_c_begin,
    left_c_end,
    right_c_begin,
    right_c_end,
    area_middle,
    P_m.rbegin(),
    t);
  CGAL_optimisation_assertion( P_m[0] == middle_index);
  // left recursive branch:
  FT area_left( 0);
  Index_cont P_l( k);
  if ( left_index < middle_index) {
    CGAL_maximum_inscribed_k_gon_2(
      points_begin,
      points_end,
      left_index,
      middle_index - 1,
      left_c_begin,
      left_c_end,
      P_m.begin() + 1,
      P_m.end(),
      k,
      area_left,
      P_l.rbegin(),
      t);
  } // if ( left_index < middle_index)
  
  
  // right recursive branch:
  FT area_right( 0);
  Index_cont P_r( k);
  if ( right_index > middle_index) {
    CGAL_maximum_inscribed_k_gon_2(
      points_begin,
      points_end,
      middle_index + 1,
      right_index,
      P_m.begin() + 1,
      P_m.end(),
      right_c_begin,
      right_c_end,
      k,
      area_right,
      P_r.rbegin(),
      t);
  } // if ( right_index > middle_index)
  
  

  if ( area_left > area_right)
    if ( area_left > area_middle) {
      // left is maximum
      max_area = area_left;
      copy( P_l.begin(), P_l.end(), o);
    }
    else {
      // middle is maximum
      max_area = area_middle;
      copy( P_m.begin(), P_m.end(), o);
    }
  else if ( area_right > area_middle) {
    // right is maximum
    max_area = area_right;
    copy( P_r.begin(), P_r.end(), o);
  }
  else {
    // middle is maximum
    max_area = area_middle;
    copy( P_m.begin(), P_m.end(), o);
  }

  return o;
} // CGAL_maximum_inscribed_k_gon_2( ... )

} //namespace CGAL

#endif // ! (CGAL_EXTREMAL_POLYGON_2_H)
