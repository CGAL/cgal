// ============================================================================
//
// Copyright (c) 1998, 1999, 2000 The CGAL Consortium
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
// file          : rectangular_p_center_2.h
// chapter       : $CGAL_Chapter: Geometric Optimisation $
// package       : $CGAL_Package: Matrix_search $
// source        : pcenter.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// maintainer    : Michael Hoffmann <hoffmann@inf.ethz.ch>
// coordinator   : ETH
//
// 2-4-Center Computation for Axis-Parallel 2D-Rectangles
// ============================================================================

#if ! (CGAL_RECTANGULAR_P_CENTER_2_H)
#define CGAL_RECTANGULAR_P_CENTER_2_H 1

#include <CGAL/pierce_rectangles_2.h>
#include <CGAL/functional.h>
#include <CGAL/sorted_matrix_search.h>
#include <CGAL/rectangular_3_center_2.h>
#include <algorithm>
#include <CGAL/number_utils_classes.h>
#ifdef CGAL_REP_CLASS_DEFINED
#include <CGAL/Rectangular_p_center_traits_2.h>
#endif // CGAL_REP_CLASS_DEFINED

CGAL_BEGIN_NAMESPACE
CGAL_END_NAMESPACE

#include <CGAL/Cartesian_matrix.h>

CGAL_BEGIN_NAMESPACE

#ifdef _MSC_VER
// that compiler cannot even distinguish between global
// and class scope, so ...
#define Base B_B_Base
#endif // _MSC_VER

template < class Operation,
           class RandomAccessIC_row,
           class RandomAccessIC_column >
class Cartesian_matrix_horizontally_flipped
: public Cartesian_matrix< Operation,
                                RandomAccessIC_row,
                                RandomAccessIC_column >
{
public:
  typedef
    Cartesian_matrix< Operation,
                           RandomAccessIC_row,
                           RandomAccessIC_column >
  Base;
  typedef typename Base::Value      Value;


  /*
  Cartesian_matrix_horizontally_flipped(
    Operation o = Operation())
  : Base( o)
  {}
  */

  Cartesian_matrix_horizontally_flipped(
    RandomAccessIC_row r_f,
    RandomAccessIC_row r_l,
    RandomAccessIC_column c_f,
    RandomAccessIC_column c_l)
  : Base( r_f, r_l, c_f, c_l)
  {}

  Cartesian_matrix_horizontally_flipped(
    RandomAccessIC_row r_f,
    RandomAccessIC_row r_l,
    RandomAccessIC_column c_f,
    RandomAccessIC_column c_l,
    const Operation& o)
  : Base( r_f, r_l, c_f, c_l, o)
  {}

  Value
  operator()( int r, int c) const
  {
    CGAL_optimisation_precondition( r >= 0 && r < number_of_rows());
    CGAL_optimisation_precondition( c >= 0 && c < number_of_columns());
    return Base::operator()( r, number_of_columns() - 1 - c);
  }
};

#ifdef _MSC_VER
#undef Base
#endif // _MSC_VER

template < class Operation,
           class RandomAccessIC_row,
           class RandomAccessIC_column >
inline
Cartesian_matrix_horizontally_flipped<
  Operation,
  RandomAccessIC_row,
  RandomAccessIC_column >
cartesian_matrix_horizontally_flipped(
  RandomAccessIC_row r_f,
  RandomAccessIC_row r_l,
  RandomAccessIC_column c_f,
  RandomAccessIC_column c_l,
  const Operation& o)
{
  return
  Cartesian_matrix_horizontally_flipped<
    Operation,
    RandomAccessIC_row,
    RandomAccessIC_column >
  ( r_f, r_l, c_f, c_l, o);
}
/*
template < class ForwardIterator,
           class OutputIterator,
           class FT,
           class PiercingFunction >
inline
OutputIterator
rectangular_p_center_2_binary_search(
  ForwardIterator f,
  ForwardIterator l,
  OutputIterator o,
  FT& r,
  const PiercingFunction& pf)
{
  return rectangular_p_center_2_binary_search(
    f,
    l,
    o,
    r,
    pf,
    CGAL_reinterpret_cast(
      Rectangular_p_center_matrix_search_traits_2< PiercingFunction >, 0));
} // rectangular_p_center_2_binary_search( ... )

template < class ForwardIterator,
           class OutputIterator,
#ifdef CGAL_CFG_MATCHING_BUG_1
           class PiercingFunction,
           class FT,
#endif
           class Traits >
OutputIterator
rectangular_p_center_2_binary_search(
  ForwardIterator f,
  ForwardIterator l,
  OutputIterator o,
#ifndef CGAL_CFG_MATCHING_BUG_1
  typename Traits::FT& r,
  const typename Traits::PiercingFunction& pf,
#else
  FT& r,
  PiercingFunction& pf,
#endif
  Traits*)
//
// preconditions:
// --------------
//  * Traits fulfills the requirements for PCenter traits classes
//  * value type of ForwardIterator is Traits::Point_2
//  * OutputIterator accepts Traits::Point_2 as value type
//  * the range [f,l) is not empty
//
// functionality:
// --------------
//
{
  CGAL_optimisation_precondition( f != l);

  // typedefs:
#ifndef CGAL_CFG_MATCHING_BUG_1
  typedef typename Traits::FT    FT;
#endif
  typedef typename Traits::X     X;
  typedef typename Traits::Y     Y;

  // create Traits object:
  Traits pierce_it( f, l, pf);

  // check, if input data is trivial
  bool ok;
  OutputIterator oi = pierce_it(FT(0), o, ok);
  if (ok) {
    r = 0;
    return oi;
  }
  // create vector with absolute coordinate differences:
  std::vector< FT > c_diffs;
  c_diffs.reserve( pierce_it.number_of_points() *
                   (pierce_it.number_of_points() - 1));
  for ( ForwardIterator i( f); i != l; ++i)
    for ( ForwardIterator j( i + 1); j != l; ++j) {
      c_diffs.push_back( CGAL_NTS abs( i->x() - j->x()));
      c_diffs.push_back( CGAL_NTS abs( i->y() - j->y()));
    }
  CGAL_optimisation_assertion(
    c_diffs.size() == pierce_it.number_of_points() *
    (pierce_it.number_of_points() - 1));
  
  // sort it:
  sort( c_diffs.begin(), c_diffs.end());
  // search it:
  int b( 0);
  int e( c_diffs.size());
  
  // invariant of the following loop:
  // forall 0 <= a < b: c_diffs[a] is infeasible  AND
  // forall e <= a < c_diffs.size(): c_diffs[a] is feasible
  while ( e > b) {
    int c = ( e + b) >> 1;
    if ( pierce_it( c_diffs[c])) {
      // c_diffs[c] is feasible
      e = c;
    }
    else {
      // c_diffs[c] is infeasible
      b = c + 1;
    }
  } // while ( e > b)
  CGAL_optimisation_assertion( e == b);
  
  // return the result:
  r = c_diffs[e];
  OutputIterator o_return( pierce_it( r, o, ok));
  CGAL_optimisation_assertion( ok);
  return o_return;

} // rectangular_p_center_2_binary_search( ... )
*/
template < class RandomAccessIC,
           class OutputIterator,
#ifdef CGAL_CFG_MATCHING_BUG_1
           class FT,
#endif
           class PiercingFunction,
           class Traits,
           class MatrixOperator >
OutputIterator
rectangular_p_center_2_matrix_search(
  RandomAccessIC f,
  RandomAccessIC l,
  OutputIterator o,
#ifndef CGAL_CFG_MATCHING_BUG_1
  typename Traits::FT& r,
#else
  FT& r,
#endif
  PiercingFunction pf,
  const Traits& t,
  const MatrixOperator& mop)
{
  int number_of_points( iterator_distance( f, l));
  CGAL_optimisation_precondition( number_of_points > 0);

#ifndef CGAL_CFG_NO_NAMESPACE
  using std::minus;
  using std::sort;
#endif

  // typedefs:
#ifndef CGAL_CFG_MATCHING_BUG_1
  typedef typename Traits::FT        FT;
#endif
  typedef std::vector< FT >          FT_cont;
  typedef typename FT_cont::iterator FT_iterator;
  typedef Cartesian_matrix_horizontally_flipped<
    MatrixOperator,
    FT_iterator,
    FT_iterator >
  Matrix;
  typedef std::vector< Matrix >      MatrixContainer;
  typedef
    Rectangular_p_center_matrix_search_traits_2< Traits, PiercingFunction >
  MSTraits;
  typedef Sorted_matrix_search_traits_adaptor< MSTraits&, Matrix >
    Matrix_search_traits;

  // create Traits object:
  MSTraits pierce_it(f, l, t, pf);

  // check, if input data is trivial
  bool ok;
  OutputIterator oi = pierce_it(FT(0), o, ok);
  if (ok) {
    r = 0;
    return oi;
  }

  // create matrix search traits:
  Matrix_search_traits search_it(pierce_it);

  // copy x and y coordinates from [f,l):
  std::vector< FT > x_coords;
  std::vector< FT > y_coords;
  x_coords.reserve( number_of_points);
  y_coords.reserve( number_of_points);
  for ( RandomAccessIC p( f); p != l; ++p) {
    x_coords.push_back(p->x());
    y_coords.push_back(p->y());
  }
  
  // sort coordinates:
  sort( x_coords.begin(), x_coords.end());
  sort( y_coords.begin(), y_coords.end());
  
  // create matrices:
  MatrixContainer matrices;
  
  // create matrix of x-differences:
  matrices.push_back(
    Matrix( x_coords.begin(),
            x_coords.end(),
            x_coords.begin(),
            x_coords.end(),
            mop));
  
  // create matrix of y-differences:
  matrices.push_back(
    Matrix( y_coords.begin(),
            y_coords.end(),
            y_coords.begin(),
            y_coords.end(),
            mop));

  // do the actual search:
  r = sorted_matrix_search(matrices.begin(),
                           matrices.end(),
                           search_it);

  // return result:
  OutputIterator o_return(pierce_it(r, o, ok));
  CGAL_optimisation_assertion(ok);
  return o_return;

} // P_center_matrix_search

template < class RandomAccessIC,
           class OutputIterator,
#ifdef CGAL_CFG_MATCHING_BUG_1
           class FT,
#endif
           class PiercingFunction,
           class Traits >
inline
OutputIterator
rectangular_p_center_2_matrix_search(
  RandomAccessIC f,
  RandomAccessIC l,
  OutputIterator o,
#ifndef CGAL_CFG_MATCHING_BUG_1
  typename Traits::FT& r,
#else
  FT& r,
#endif
  const PiercingFunction& pf,
  const Traits& t)
{
#ifndef CGAL_CFG_MATCHING_BUG_1
  typedef typename Traits::FT FT;
#endif
#ifndef CGAL_CFG_NO_NAMESPACE
  using std::minus;
#endif

  return rectangular_p_center_2_matrix_search(
    f,
    l,
    o,
    r,
    pf,
    t,
    compose(bind_1(Max< FT >(), 0), minus< FT >()));

} // Pcenter_matrix_search( ... )



template < class ForwardIterator, class OutputIterator, class FT >
inline OutputIterator
rectangular_p_center_matrix_search_2(
  ForwardIterator f,
  ForwardIterator l,
  OutputIterator o,
  FT& r,
  int p)
{
  CGAL_optimisation_precondition(p >= 2 && p < 5);
  if (p == 2)
    return rectangular_p_center_2_matrix_search(
      f, l, o, r, Two_covering_algorithm(), t);
  else if (p == 3)
    return rectangular_p_center_2_matrix_search(
      f, l, o, r, Three_covering_algorithm(), t);
  return rectangular_p_center_2_matrix_search(
    f, l, o, r, Four_covering_algorithm(), t);
} // rectangular_p_center_matrix_search_2(f, l, o, r, p)

template < class ForwardIterator,
           class OutputIterator,
#ifdef CGAL_CFG_MATCHING_BUG_1
           class FT,
#endif
           class Traits >
inline OutputIterator
rectangular_p_center_2(ForwardIterator f,
                       ForwardIterator l,
                       OutputIterator o,
#ifndef CGAL_CFG_MATCHING_BUG_1
                       typename Traits::FT& r,
#else
                       FT& r,
#endif
                       int p,
                       Traits& t)
{
  CGAL_optimisation_precondition(p >= 2 && p < 5);

  if (p == 2)
    return rectangular_2_center_2(f, l, o, r, t);
  else if (p == 3)
    return rectangular_3_center_2(f, l, o, r, t);
  return rectangular_p_center_2_matrix_search(
    f, l, o, r, Four_covering_algorithm(), t);

} // rectangular_p_center_2( ... )

template < class ForwardIterator, class OutputIterator, class FT >
inline OutputIterator
rectangular_p_center_2(ForwardIterator f,
                       ForwardIterator l,
                       OutputIterator o,
                       FT& r,
                       int p)
{
  CGAL_optimisation_precondition(p >= 2 && p < 5);
  typedef typename
    std::iterator_traits< ForwardIterator >::value_type::R R;
  Rectangular_p_center_default_traits_2< R > t;

  return rectangular_p_center_2(f, l, o, r, p, t);
} // rectangular_p_center_2( ... )


CGAL_END_NAMESPACE

#endif // ! (CGAL_RECTANGULAR_P_CENTER_2_H)

// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

