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
// file          : all_furthest_neighbors_2.h
// chapter       : $CGAL_Chapter: Geometric Optimisation $
// package       : $CGAL_Package: Matrix_search $
// source        : mon_search.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// coordinator   : ETH Zurich (Bernd Gaertner <gaertner@inf.ethz.ch>)
//
// Compute all furthest neighbors for the vertices of a convex polygon
// ============================================================================

#if ! (CGAL_ALL_FURTHEST_NEIGHBORS_2_H)
#define CGAL_ALL_FURTHEST_NEIGHBORS_2_H 1

#ifndef CGAL_OPTIMISATION_ASSERTIONS_H
#include <CGAL/optimisation_assertions.h>
#endif // CGAL_OPTIMISATION_ASSERTIONS_H

#ifdef CGAL_REP_CLASS_DEFINED
#ifndef CGAL_ALL_FURTHEST_NEIGHBORS_TRAITS_2_H
#include <CGAL/All_furthest_neighbors_traits_2.h>
#endif // CGAL_ALL_FURTHEST_NEIGHBORS_TRAITS_2_H
#endif // CGAL_REP_CLASS_DEFINED

#ifndef CGAL_CARTESIAN_MATRIX_H
#include <CGAL/Cartesian_matrix.h>
#endif // CGAL_CARTESIAN_MATRIX_H
#ifndef CGAL_DYNAMIC_MATRIX_H
#include <CGAL/Dynamic_matrix.h>
#endif // CGAL_DYNAMIC_MATRIX_H
#ifndef CGAL_MONOTONE_MATRIX_SEARCH_H
#include <CGAL/monotone_matrix_search.h>
#endif // CGAL_MONOTONE_MATRIX_SEARCH_H
#ifndef CGAL_PROTECT_FUNCTION_H
#include <function.h>
#define CGAL_PROTECT_FUNCTION_H
#endif // CGAL_PROTECT_FUNCTION_H
#ifndef CGAL_PROTECT_ALGO_H
#include <algo.h>
#define CGAL_PROTECT_ALGO_H
#endif // CGAL_PROTECT_ALGO_H

template < class Operation, class RandomAccessIC >
class CGAL__All_furthest_neighbor_matrix
: public CGAL_Cartesian_matrix< Operation,
                                RandomAccessIC,
                                RandomAccessIC >
// represents the matrix used for computing
// all furthest neighbors of a convex polygon
{
public:
  typedef
    CGAL_Cartesian_matrix< Operation, RandomAccessIC, RandomAccessIC >
  BaseClass;

  CGAL__All_furthest_neighbor_matrix( RandomAccessIC f,
                                      RandomAccessIC l)
  : BaseClass( f, l, f, l)
  {}

  int
  number_of_columns() const
  { return 2 * number_of_rows() - 1; }

  Value
  operator()( int r, int c) const
  {
    CGAL_optimisation_precondition( r >= 0 && r < number_of_rows());
    CGAL_optimisation_precondition( c >= 0 && c < number_of_columns());
    if ( c <= r)
      return Value( c - r);
    else if ( c >= r + number_of_rows())
      return Value( 0);
    else if ( c >= number_of_rows())
      return BaseClass::operator()( r, c - number_of_rows());
    else
      return BaseClass::operator()( r, c);
  }
};

#if !defined(CGAL_CFG_NO_ITERATOR_TRAITS) && \
!defined(CGAL_CFG_MATCHING_BUG_2)

#ifndef CGAL_PROTECT_ITERATOR_H
#include <iterator.h>
#define CGAL_PROTECT_ITERATOR_H
#endif // CGAL_PROTECT_ITERATOR_H

template < class RandomAccessIC, class OutputIterator, class Traits >
inline
OutputIterator
CGAL_all_furthest_neighbors( RandomAccessIC points_begin,
                             RandomAccessIC points_end,
                             OutputIterator o,
                             const Traits& t)
{
  return
  CGAL_all_furthest_neighbors(
    points_begin,
    points_end,
    o,
    t,
    iterator_traits< OutputIterator >::iterator_category());
} // CGAL_all_furthest_neighbors( ... )

template < class RandomAccessIC,
           class OutputIterator,
           class Traits >
OutputIterator
CGAL_all_furthest_neighbors( RandomAccessIC points_begin,
                             RandomAccessIC points_end,
                             OutputIterator o,
                             const Traits& t,
                             random_access_iterator_tag)
{
  typedef CGAL__All_furthest_neighbor_matrix<
    typename Traits::Distance,
    RandomAccessIC
    >
  Afn_matrix;

  // check preconditions:
  int number_of_points(
    CGAL_iterator_distance( points_begin, points_end));
  CGAL_optimisation_precondition( number_of_points > 0);
  CGAL_optimisation_expensive_precondition(
    t.is_convex( points_begin, points_end));

  // compute maxima:
  CGAL_monotone_matrix_search(
    CGAL_dynamic_matrix(
      Afn_matrix( points_begin, points_end)),
    o);

  return o + number_of_points;
} // CGAL_all_furthest_neighbors( ... )

template < class RandomAccessIC,
           class OutputIterator,
           class Traits,
           class IteratorCategory >
OutputIterator
CGAL_all_furthest_neighbors( RandomAccessIC points_begin,
                             RandomAccessIC points_end,
                             OutputIterator o,
                             const Traits& t,
                             IteratorCategory)
#else
template < class RandomAccessIC,
           class OutputIterator,
           class Traits >
OutputIterator
CGAL_all_furthest_neighbors( RandomAccessIC points_begin,
                             RandomAccessIC points_end,
                             OutputIterator o,
                             const Traits& t)
#endif
{
  typedef CGAL__All_furthest_neighbor_matrix<
    typename Traits::Distance,
    RandomAccessIC
    >
  Afn_matrix;

 // check preconditions:
  int number_of_points(
    CGAL_iterator_distance( points_begin, points_end));
  CGAL_optimisation_precondition( number_of_points > 0);
  CGAL_optimisation_expensive_precondition(
    t.is_convex( points_begin, points_end));

  // prepare random access container:
  vector< int > v( number_of_points);

  // compute maxima:
  CGAL_monotone_matrix_search(
    CGAL_dynamic_matrix(
      Afn_matrix( points_begin, points_end)),
    v.begin());

  // output result:
  return transform( v.begin(),
                    v.end(),
                    o,
                    bind2nd( modulus< int >(), number_of_points));
} // CGAL_all_furthest_neighbors( ... )

#endif // ! (CGAL_ALL_FURTHEST_NEIGHBORS_2_H)

// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

