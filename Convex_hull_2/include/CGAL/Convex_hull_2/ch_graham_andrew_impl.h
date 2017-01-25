// Copyright (c) 1999  Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Stefan Schirra


#ifndef CGAL_CH_GRAHAM_ANDREW_C
#define CGAL_CH_GRAHAM_ANDREW_C

#include <CGAL/license/Convex_hull_2.h>


#ifndef CGAL_CH_NO_POSTCONDITIONS
#include <CGAL/convexity_check_2.h>
#endif // CGAL_CH_NO_POSTCONDITIONS

#include <CGAL/Convex_hull_2/ch_assertions.h>
#include <CGAL/algorithm.h>
#include <CGAL/IO/Tee_for_output_iterator.h>
#include <vector>
#include <algorithm>

namespace CGAL {

template <class BidirectionalIterator, class OutputIterator, class Traits>
OutputIterator
ch_graham_andrew_scan( BidirectionalIterator first,
                       BidirectionalIterator last,
                       OutputIterator        result,
                       const Traits&         ch_traits)
{
  typedef  typename Traits::Left_turn_2  Left_turn;

  std::vector< BidirectionalIterator >    S;
  BidirectionalIterator              alpha;
  BidirectionalIterator              beta;
  BidirectionalIterator              iter;
  CGAL_ch_precondition( first != last );
  CGAL_ch_precondition( cpp11::next(first) != last );

  --last;
  CGAL_ch_precondition( *first != *last );
  S.push_back( last  );
  S.push_back( first );
  Left_turn    left_turn = ch_traits.left_turn_2_object();


  iter = first;
  do
  {
      ++iter;
  }
  while (( iter != last ) && !left_turn(*last, *first, *iter) );

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
          if ( left_turn(*alpha, *iter, *last) )
          {
              while ( !left_turn(*beta, *alpha, *iter) )
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
  typedef  typename Traits::Point_2     Point_2;
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
  typedef  typename Traits::Left_turn_2  Left_turn;

  CGAL_ch_precondition_code(
  typedef  typename Traits::Equal_2      Equal_2;
  Equal_2      equal_points = ch_traits.equal_2_object();
  )

  Left_turn    left_turn    = ch_traits.left_turn_2_object();

  std::vector< BidirectionalIterator >    S;
  BidirectionalIterator              alpha;
  BidirectionalIterator              beta;
  BidirectionalIterator              iter;
  CGAL_ch_precondition( first != last );
  CGAL_ch_precondition( cpp11::next(first) != last );

  --last;
  CGAL_ch_precondition(! equal_points(*first,*last) );
  S.push_back( last  );
  S.push_back( first );

  iter = first;
  do
  {
      ++iter;
  }
  while (( iter != last ) && !left_turn(*last, *first, *iter) );

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
          if ( left_turn(*alpha, *iter, *last) )
          {
              while ( !left_turn(*beta, *alpha, *iter) )
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
  typedef  typename Traits::Point_2     Point_2;
  typedef  typename Traits::Equal_2      Equal_2;  
  
  Equal_2      equal_points = ch_traits.equal_2_object();  

  if (first == last) return result;
  std::vector< Point_2 >  V (first, last);
  std::sort( V.begin(), V.end(), ch_traits.less_xy_2_object() );
  if (equal_points( *(V.begin()), *(V.rbegin())) )
  {
      *result++ = *(V.begin());
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
  typedef  typename Traits::Point_2      Point_2;
  typedef  typename Traits::Equal_2      Equal_2;  
  
  Equal_2      equal_points = ch_traits.equal_2_object();    

  if (first == last) return result;
  std::vector< Point_2 >  V (first, last);
  std::sort( V.begin(), V.end(), ch_traits.less_xy_2_object() );
  if (equal_points( *(V.begin()), *(V.rbegin())) )
  {
      *result++ = *(V.begin());
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
  typedef  typename Traits::Point_2      Point_2;
  typedef  typename Traits::Equal_2      Equal_2;  
  
  Equal_2      equal_points = ch_traits.equal_2_object();     

  if (first == last) return result;
  std::vector< Point_2 >  V (first, last);
  std::sort( V.begin(), V.end(), ch_traits.less_xy_2_object() );
  if (equal_points( *(V.begin()), *(V.rbegin())) )
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
} //namespace CGAL

#endif // CGAL_CH_GRAHAM_ANDREW_C
