// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// release       : 
// release_date  : 2000, August 03
//
// file          : ch_graham_andrew.C
// package       : Convex_hull (3.3)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// source        : convex_hull_2.lw
// revision      : 3.3
// revision_date : 03 Aug 2000
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================


#ifndef CGAL_CH_GRAHAM_ANDREW_C
#define CGAL_CH_GRAHAM_ANDREW_C

#include <CGAL/stl_extensions.h>
#ifndef CGAL_CH_GRAHAM_ANDREW_H
#include <CGAL/ch_graham_andrew.h>
#endif // CGAL_CH_GRAHAM_ANDREW_H

CGAL_BEGIN_NAMESPACE

template <class BidirectionalIterator, class OutputIterator, class Traits>
OutputIterator
ch_graham_andrew_scan( BidirectionalIterator first,
                       BidirectionalIterator last,
                       OutputIterator        result,
                       const Traits&         ch_traits)
{
  typedef  typename Traits::Less_xy_2   Less_xy;
  typedef  typename Traits::Point_2     Point_2;
  typedef  typename Traits::Leftturn_2  Leftturn;

  std::vector< BidirectionalIterator >    S;
  BidirectionalIterator              alpha;
  BidirectionalIterator              beta;
  BidirectionalIterator              iter;
  CGAL_ch_precondition( first != last );
  CGAL_ch_precondition( successor(first) != last );

  --last;
  CGAL_ch_precondition( *first != *last );
  S.push_back( last  );
  S.push_back( first );
  Leftturn    leftturn = ch_traits.leftturn_2_object();


  iter = first;
  do
  {
      ++iter;
  }
  while (( iter != last ) && !leftturn(*last, *first, *iter) );

  if ( iter != last )
  {
      S.push_back( iter );
      typedef typename std::vector<BidirectionalIterator>::reverse_iterator  
              rev_iterator;
      rev_iterator  stack_rev_iter = S.rbegin(); 
      alpha = iter;
      beta  = *++stack_rev_iter;

      for ( ++iter ; iter != last; ++iter )
      {
          if ( leftturn(*alpha, *iter, *last) )
          {
              while ( !leftturn(*beta, *alpha, *iter) )
              {
                  S.pop_back();
                  alpha = beta;
                  stack_rev_iter = S.rbegin();
                  beta  = *++stack_rev_iter;
                  CGAL_ch_assertion(S.size() >= 2);
              }
              S.push_back( iter  );
              beta = alpha;
              alpha = iter;
          }
      }

  }

  typedef typename std::vector< BidirectionalIterator >::iterator std_iterator;
  std_iterator  stack_iter = S.begin();
  #if defined(CGAL_CH_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
    || defined(NDEBUG)
  OutputIterator  res(result);
  #else
  Tee_for_output_iterator<OutputIterator,Point_2> res(result);
  #endif // no postconditions ...
  for ( ++stack_iter;  stack_iter != S.end(); ++stack_iter)
  { *res =  **stack_iter;  ++res; }
  CGAL_ch_postcondition( \
      is_ccw_strongly_convex_2( res.output_so_far_begin(), \
                                     res.output_so_far_end(), \
                                     ch_traits));
  CGAL_ch_expensive_postcondition( \
      ch_brute_force_chain_check_2( \
          first, last, \
          res.output_so_far_begin(), res.output_so_far_end(), \
          ch_traits));
  #if defined(CGAL_CH_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
    || defined(NDEBUG)
  return res;
  #else
  return res.to_output_iterator();
  #endif // no postconditions ...
}

template <class BidirectionalIterator, class OutputIterator, class Traits>
OutputIterator
ch__ref_graham_andrew_scan( BidirectionalIterator first,
                                 BidirectionalIterator last,
                                 OutputIterator&       result,
                                 const Traits&         ch_traits)
{
  typedef  typename Traits::Less_xy_2   Less_xy;
  typedef  typename Traits::Point_2     Point_2;
  typedef  typename Traits::Leftturn_2  Leftturn;

  std::vector< BidirectionalIterator >    S;
  BidirectionalIterator              alpha;
  BidirectionalIterator              beta;
  BidirectionalIterator              iter;
  CGAL_ch_precondition( first != last );
  CGAL_ch_precondition( successor(first) != last );

  --last;
  CGAL_ch_precondition( *first != *last );
  S.push_back( last  );
  S.push_back( first );
  Leftturn    leftturn = ch_traits.leftturn_2_object();


  iter = first;
  do
  {
      ++iter;
  }
  while (( iter != last ) && !leftturn(*last, *first, *iter) );

  if ( iter != last )
  {
      S.push_back( iter );
      typedef typename std::vector<BidirectionalIterator>::reverse_iterator  
              rev_iterator;
      rev_iterator  stack_rev_iter = S.rbegin(); 
      alpha = iter;
      beta  = *++stack_rev_iter;

      for ( ++iter ; iter != last; ++iter )
      {
          if ( leftturn(*alpha, *iter, *last) )
          {
              while ( !leftturn(*beta, *alpha, *iter) )
              {
                  S.pop_back();
                  alpha = beta;
                  stack_rev_iter = S.rbegin();
                  beta  = *++stack_rev_iter;
                  CGAL_ch_assertion(S.size() >= 2);
              }
              S.push_back( iter  );
              beta = alpha;
              alpha = iter;
          }
      }

  }

  typedef typename std::vector< BidirectionalIterator >::iterator std_iterator;
  std_iterator  stack_iter = S.begin();
  for ( ++stack_iter;  stack_iter != S.end(); ++stack_iter)
  { *result =  **stack_iter;  ++result; }
  return result;
}

template <class InputIterator, class OutputIterator, class Traits>
OutputIterator
ch_graham_andrew( InputIterator  first,
                       InputIterator  last,
                       OutputIterator result,
                       const Traits&  ch_traits)
{
  typedef  typename Traits::Less_xy_2   Less_xy;
  typedef  typename Traits::Point_2     Point_2;
  typedef  typename Traits::Leftturn_2  Leftturn;

  if (first == last) return result;
  std::vector< Point_2 >  V;
  std::copy( first, last, std::back_inserter(V) );
  std::sort( V.begin(), V.end(), ch_traits.less_xy_2_object() );
  if ( *(V.begin()) == *(V.rbegin()) )
  {
      *result = *(V.begin());  ++result;
      return result;
  }

  #if defined(CGAL_CH_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
    || defined(NDEBUG)
  OutputIterator  res(result);
  #else
  Tee_for_output_iterator<OutputIterator,Point_2> res(result);
  #endif // no postconditions ...
  ch__ref_graham_andrew_scan( V.begin(), V.end(),  res, ch_traits);
  ch__ref_graham_andrew_scan( V.rbegin(), V.rend(), res, ch_traits);
  CGAL_ch_postcondition( \
      is_ccw_strongly_convex_2( res.output_so_far_begin(), \
                                     res.output_so_far_end(), \
                                     ch_traits));
  CGAL_ch_expensive_postcondition( \
      ch_brute_force_check_2( \
          V.begin(), V.end(), \
          res.output_so_far_begin(), res.output_so_far_end(), \
          ch_traits));
  #if defined(CGAL_CH_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
    || defined(NDEBUG)
  return res;
  #else
  return res.to_output_iterator();
  #endif // no postconditions ...


}

template <class InputIterator, class OutputIterator, class Traits>
OutputIterator
ch_lower_hull_scan( InputIterator  first,
                         InputIterator  last,
                         OutputIterator result,
                         const Traits&  ch_traits)
{
  typedef  typename Traits::Less_xy_2   Less_xy;
  typedef  typename Traits::Point_2     Point_2;
  typedef  typename Traits::Leftturn_2  Leftturn;

  if (first == last) return result;
  std::vector< Point_2 >  V;
  std::copy( first, last, std::back_inserter(V) );
  std::sort( V.begin(), V.end(), ch_traits.less_xy_2_object() );
  if ( *(V.begin()) == *(V.rbegin()) )
  {
      *result = *(V.begin());  ++result;
      return result;
  }

  #if defined(CGAL_CH_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
    || defined(NDEBUG)
  OutputIterator  res(result);
  #else
  Tee_for_output_iterator<OutputIterator,Point_2> res(result);
  #endif // no postconditions ...
  ch_graham_andrew_scan( V.begin(), V.end(), res, ch_traits);
  #if defined(CGAL_CH_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
    || defined(NDEBUG)
  return res;
  #else
  return res.to_output_iterator();
  #endif // no postconditions ...
}
template <class InputIterator, class OutputIterator, class Traits>
OutputIterator
ch_upper_hull_scan( InputIterator  first,
                         InputIterator  last,
                         OutputIterator result,
                         const Traits&  ch_traits)
{
  typedef  typename Traits::Less_xy_2   Less_xy;
  typedef  typename Traits::Point_2     Point_2;
  typedef  typename Traits::Leftturn_2  Leftturn;

  if (first == last) return result;
  std::vector< Point_2 >  V;
  std::copy( first, last, std::back_inserter(V) );
  std::sort( V.begin(), V.end(), ch_traits.less_xy_2_object() );
  if ( *(V.begin()) == *(V.rbegin()) )
  { return result; }
  #if defined(CGAL_CH_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
    || defined(NDEBUG)
  OutputIterator  res(result);
  #else
  Tee_for_output_iterator<OutputIterator,Point_2> res(result);
  #endif // no postconditions ...
  ch_graham_andrew_scan( V.rbegin(), V.rend(), res, ch_traits);
  #if defined(CGAL_CH_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
    || defined(NDEBUG)
  return res;
  #else
  return res.to_output_iterator();
  #endif // no postconditions ...
}
CGAL_END_NAMESPACE

#endif // CGAL_CH_GRAHAM_ANDREW_C
