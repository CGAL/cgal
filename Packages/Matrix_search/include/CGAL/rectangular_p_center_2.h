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
// file          : rectangular_p_center_2.h
// chapter       : $CGAL_Chapter: Geometric Optimisation $
// package       : $CGAL_Package: Matrix_search $
// source        : pcenter.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// coordinator   : ETH Zurich (Bernd Gaertner <gaertner@inf.ethz.ch>)
//
// 2-4-Center Computation for Axis-Parallel 2D-Rectangles
// ============================================================================

#if ! (CGAL_RECTANGULAR_P_CENTER_2_H)
#define CGAL_RECTANGULAR_P_CENTER_2_H 1

#include <CGAL/pierce_rectangles_2.h>
#include <CGAL/function_objects.h>
#include <CGAL/sorted_matrix_search.h>
#include <algo.h>
#ifdef CGAL_REP_CLASS_DEFINED
#include <CGAL/rectangular_p_center_2_traits.h>
#endif // CGAL_REP_CLASS_DEFINED

#include <CGAL/Cartesian_matrix.h>

template < class Operation,
           class RandomAccessIC_row,
           class RandomAccessIC_column >
class CGAL_Cartesian_matrix_horizontally_flipped
: public CGAL_Cartesian_matrix< Operation,
                                RandomAccessIC_row,
                                RandomAccessIC_column >
{
public:
  typedef
    CGAL_Cartesian_matrix< Operation,
                           RandomAccessIC_row,
                           RandomAccessIC_column >
  Base;

  /*
  CGAL_Cartesian_matrix_horizontally_flipped(
    Operation o = Operation())
  : Base( o)
  {}
  */

  CGAL_Cartesian_matrix_horizontally_flipped(
    RandomAccessIC_row r_f,
    RandomAccessIC_row r_l,
    RandomAccessIC_column c_f,
    RandomAccessIC_column c_l)
  : Base( r_f, r_l, c_f, c_l)
  {}

  CGAL_Cartesian_matrix_horizontally_flipped(
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

template < class Operation,
           class RandomAccessIC_row,
           class RandomAccessIC_column >
inline
CGAL_Cartesian_matrix_horizontally_flipped<
  Operation,
  RandomAccessIC_row,
  RandomAccessIC_column >
CGAL_cartesian_matrix_horizontally_flipped(
  RandomAccessIC_row r_f,
  RandomAccessIC_row r_l,
  RandomAccessIC_column c_f,
  RandomAccessIC_column c_l,
  const Operation& o)
{
  return
  CGAL_Cartesian_matrix_horizontally_flipped<
    Operation,
    RandomAccessIC_row,
    RandomAccessIC_column >
  ( r_f, r_l, c_f, c_l, o);
}
template < class ForwardIterator,
           class OutputIterator,
           class FT,
           class PiercingFunction >
inline
OutputIterator
CGAL_rectangular_p_center_2_binary_search(
  ForwardIterator f,
  ForwardIterator l,
  OutputIterator o,
  FT& r,
  const PiercingFunction& pf)
{
  return CGAL_rectangular_p_center_2_binary_search(
    f,
    l,
    o,
    r,
    pf,
    CGAL_reinterpret_cast(
      CGAL_Pcenter_default_traits< PiercingFunction >, 0));
} // CGAL_rectangular_p_center_2_binary_search( ... )

template < class ForwardIterator,
           class OutputIterator,
#ifdef CGAL_CFG_MATCHING_BUG_1
           class PiercingFunction,
           class FT,
#endif
           class Traits >
OutputIterator
CGAL_rectangular_p_center_2_binary_search(
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
  #ifdef CGAL_PCENTER_TRACE
  cerr << "Pcenter binary search" << endl;
  #endif

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
  OutputIterator oi( pierce_it( FT( 0), o, ok));
  if ( ok)
    return oi;
  // create vector with absolute coordinate differences:
  vector< FT > c_diffs;
  X x;
  Y y;
  c_diffs.reserve( pierce_it.number_of_points() *
                   (pierce_it.number_of_points() - 1));
  for ( ForwardIterator i( f); i != l; ++i)
    for ( ForwardIterator j( i + 1); j != l; ++j) {
      c_diffs.push_back( CGAL_abs( x(*i) - x(*j)));
      c_diffs.push_back( CGAL_abs( y(*i) - y(*j)));
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
    #ifdef CGAL_PCENTER_TRACE
    cerr << "* (" << b << "-" << e << ") try diameter c_diffs[" <<
      c << "] = " << c_diffs[c] << endl;
    #endif
    if ( pierce_it( c_diffs[c])) {
      // c_diffs[c] is feasible
      e = c;
      #ifdef CGAL_PCENTER_TRACE
      cerr << "feasible" << endl;
      #endif
    }
    else {
      // c_diffs[c] is infeasible
      b = c + 1;
      #ifdef CGAL_PCENTER_TRACE
      cerr << "infeasible" << endl;
      #endif
    }
  } // while ( e > b)
  CGAL_optimisation_assertion( e == b);
  
  // return the result:
  r = c_diffs[e];
  OutputIterator o_return( pierce_it( r, o, ok));
  CGAL_optimisation_assertion( ok);
  return o_return;

} // CGAL_rectangular_p_center_2_binary_search( ... )
template < class RandomAccessIC,
           class OutputIterator,
#ifdef CGAL_CFG_MATCHING_BUG_1
           class FT,
#endif
           class PiercingFunction >
inline
OutputIterator
CGAL_rectangular_p_center_2_matrix_search(
  RandomAccessIC f,
  RandomAccessIC l,
  OutputIterator o,
#ifndef CGAL_CFG_MATCHING_BUG_1
  typename PiercingFunction::FT& r,
#else
  FT& r,
#endif
  const PiercingFunction& pf)
{
#ifndef CGAL_CFG_MATCHING_BUG_1
  typedef typename PiercingFunction::FT FT;
#endif
  return CGAL_rectangular_p_center_2_matrix_search(
    f,
    l,
    o,
    r,
    pf,
    CGAL_reinterpret_cast(
      CGAL_Pcenter_default_traits< PiercingFunction >*, 0),
    CGAL_compose1_2(
      bind1st( CGAL_Max< FT >(), 0),
      minus< FT >()));

} // CGAL_Pcenter_matrix_search( ... )

template < class RandomAccessIC,
           class OutputIterator,
#ifdef CGAL_CFG_MATCHING_BUG_1
           class FT,
           class PiercingFunction,
#endif
           class Traits,
           class MatrixOperator >
OutputIterator
CGAL_rectangular_p_center_2_matrix_search(
  RandomAccessIC f,
  RandomAccessIC l,
  OutputIterator o,
#ifndef CGAL_CFG_MATCHING_BUG_1
  typename Traits::FT& r,
  const typename Traits::PiercingFunction& pf,
#else
  FT& r,
  PiercingFunction& pf,
#endif
  Traits*,
  const MatrixOperator&)
{
  int number_of_points( CGAL_iterator_distance( f, l));
  CGAL_optimisation_precondition( number_of_points > 0);
  #ifdef CGAL_PCENTER_TRACE
  cerr << "Pcenter matrix search" << endl;
  #endif

  // typedefs:
#ifndef CGAL_CFG_MATCHING_BUG_1
  typedef typename Traits::FT       FT;
  typedef typename Traits::PiercingFunction
    PiercingFunction;
#endif
  typedef typename Traits::X        X;
  typedef typename Traits::Y        Y;
  typedef vector< FT >              FT_cont;
  typedef FT_cont::iterator         FT_iterator;
  typedef CGAL_Cartesian_matrix_horizontally_flipped<
    MatrixOperator,
    FT_iterator,
    FT_iterator >
  Matrix;
  typedef vector< Matrix >          MatrixContainer;
  typedef CGAL_Sorted_matrix_search_traits_adaptor<
    Traits, Matrix >
  Matrix_search_traits;

  // create Traits object:
  Traits pierce_it( f, l, pf);

  // check, if input data is trivial
  bool ok;
  OutputIterator oi( pierce_it( FT( 0), o, ok));
  if ( ok)
    return oi;

  // create matrix search traits:
  Matrix_search_traits search_it( pierce_it);

  // copy x and y coordinates from [f,l):
  vector< FT > x_coords;
  vector< FT > y_coords;
  X x;
  Y y;
  x_coords.reserve( number_of_points);
  y_coords.reserve( number_of_points);
  for ( RandomAccessIC p( f); p != l; ++p) {
    x_coords.push_back( x( *p));
    y_coords.push_back( y( *p));
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
            CGAL_compose1_2(
              bind1st( CGAL_Max< FT >(), 0),
              minus< FT >())));
  // create matrix of y-differences:
  matrices.push_back(
    Matrix( y_coords.begin(),
            y_coords.end(),
            y_coords.begin(),
            y_coords.end(),
            CGAL_compose1_2(
              bind1st( CGAL_Max< FT >(), 0),
              minus< FT >())));

  // do the actual search:
  r = CGAL_sorted_matrix_search( matrices.begin(),
                                 matrices.end(),
                                 search_it);

  // return result:
  OutputIterator o_return( pierce_it( r, o, ok));
  CGAL_optimisation_assertion( ok);
  return o_return;

} // CGAL_P_center_matrix_search


template < class ForwardIterator,
           class OutputIterator,
#ifdef CGAL_CFG_MATCHING_BUG_1
           class FT,
#endif
           class Traits >
inline
OutputIterator
CGAL_rectangular_p_center_2(
  ForwardIterator f,
  ForwardIterator l,
  OutputIterator o,
#ifndef CGAL_CFG_MATCHING_BUG_1
  typename Traits::FT& r,
#else
  FT& r,
#endif
  int p,
  const Traits&)
{
  CGAL_optimisation_precondition( p >= 2 && p <= 4);
  if ( p == 2)
    return CGAL_rectangular_p_center_2_matrix_search(
      f,
      l,
      o,
      r,
      CGAL_Two_piercing_algorithm< Traits >());
  else if ( p == 3)
    return CGAL_rectangular_p_center_2_matrix_search(
      f,
      l,
      o,
      r,
      CGAL_Three_piercing_algorithm< Traits >());
  else
    return CGAL_rectangular_p_center_2_matrix_search(
      f,
      l,
      o,
      r,
      CGAL_Four_piercing_algorithm< Traits >());
} // CGAL_rectangular_p_center_2( ... )

template < class ForwardIterator,
           class OutputIterator,
           class FT >
inline
OutputIterator
CGAL_rectangular_p_center_2( ForwardIterator f,
                             ForwardIterator l,
                             OutputIterator o,
                             FT& r,
                             int p)
{
  return _CGAL_rectangular_p_center_2(
    f, l, o, r, p, value_type( f));
} // CGAL_rectangular_p_center_2( ... )
 

#endif // ! (CGAL_RECTANGULAR_P_CENTER_2_H)

// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

