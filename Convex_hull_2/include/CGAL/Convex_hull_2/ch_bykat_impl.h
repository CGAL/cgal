// Copyright (c) 1999  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stefan Schirra

#ifndef CGAL_CH_BYKAT_C
#define CGAL_CH_BYKAT_C

#include <CGAL/license/Convex_hull_2.h>


#ifndef CGAL_CH_NO_POSTCONDITIONS
#include <CGAL/convexity_check_2.h>
#endif // CGAL_CH_NO_POSTCONDITIONS

#include <CGAL/Convex_hull_2/ch_assertions.h>
#include <CGAL/ch_selected_extreme_points_2.h>
#include <CGAL/ch_graham_andrew.h>
#include <CGAL/algorithm.h>
#include <list>
#include <algorithm>

namespace CGAL {
template <class InputIterator, class OutputIterator, class Traits>
OutputIterator
ch_bykat(InputIterator first, InputIterator last,
              OutputIterator  result,
              const Traits& ch_traits)
{
  using namespace boost;

  typedef typename Traits::Point_2                         Point_2;
  typedef typename Traits::Left_turn_2                     Left_turn_2;
  typedef typename Traits::Less_signed_distance_to_line_2  Less_dist;
  typedef typename Traits::Equal_2                         Equal_2;

  Left_turn_2 left_turn    = ch_traits.left_turn_2_object();
  Equal_2     equal_points = ch_traits.equal_2_object();

  if (first == last) return result;

  std::vector< Point_2 >       P (first, last);      // Points in subsets
  std::vector< Point_2 >       H;      // right endpoints of subproblems
  H.reserve(16);

  typedef typename std::vector< Point_2 >::iterator   PointIterator;
  std::vector< PointIterator > L;      // start of subset range
  std::vector< PointIterator > R;      // end of subset range
  L.reserve(16);
  R.reserve(16);
  PointIterator           l;
  PointIterator           r;

  ch_we_point(P.begin(), P.end(), l, r, ch_traits);
  Point_2 a = *l;
  Point_2 b = *r;
  if (equal_points(a,b))
  {
      *result = a;  ++result;
      return result;
  }
  #if defined(CGAL_CH_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS)
  OutputIterator  res(result);
  #else
  Tee_for_output_iterator<OutputIterator,Point_2> res(result);
  #endif // no postconditions ...
  H.push_back( a );
  L.push_back( P.begin() );
  R.push_back( l = std::partition(P.begin(), P.end(),
                                  [&left_turn, &a, &b](const Point_2& p){ return left_turn(a,b,p); }) );
  r = std::partition( l, P.end(), [&left_turn, &a, &b](const Point_2& p){ return left_turn(b,a,p); });

  for (;;)
  {
    //  This functor must be in the for loop so that the Convex_hull_constructive traits_2 works correctly
    Less_dist less_dist = ch_traits.less_signed_distance_to_line_2_object();
      if ( l != r)
      {
        Point_2 c = *std::min_element( l, r, [&less_dist,&a,&b](const Point_2&p1, const Point_2& p2)
                                             { return less_dist(a, b, p1, p2); });
          H.push_back( b );
          L.push_back( l );
          R.push_back( l = std::partition(l, r, [&left_turn,&c,&b](const Point_2&p)
                                                { return left_turn(b, c, p); }));
          r = std::partition(l, r, [&left_turn,&c,&a](const Point_2&p)
                                   { return left_turn(c, a, p); });
          b = c;
      }
      else
      {
          *res = a;  ++res;
          if ( L.empty() ) break;
          a = b;
          b = H.back(); H.pop_back();
          l = L.back(); L.pop_back();
          r = R.back(); R.pop_back();
      }
  }
  CGAL_ch_postcondition( \
      is_ccw_strongly_convex_2( res.output_so_far_begin(), \
                                     res.output_so_far_end(), \
                                     ch_traits));
  CGAL_ch_expensive_postcondition( \
      ch_brute_force_check_2( \
          P.begin(), P.end(), \
          res.output_so_far_begin(), res.output_so_far_end(), \
          ch_traits));
  #if defined(CGAL_CH_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS)
  return res;
  #else
  return res.to_output_iterator();
  #endif // no postconditions ...
}

#define CGAL_ch_THRESHOLD 10
template <class InputIterator, class OutputIterator, class Traits>
OutputIterator
ch_bykat_with_threshold(InputIterator   first, InputIterator last,
                             OutputIterator  result,
                             const Traits&   ch_traits)
{
  using namespace boost;

  typedef typename Traits::Point_2               Point_2;
  typedef typename Traits::Left_turn_2            Left_turn_2;
  typedef typename Traits::Less_signed_distance_to_line_2
                                                 Less_dist;
  typedef typename std::vector< Point_2 >::iterator
                                                 PointIterator;
  typedef typename Traits::Equal_2                         Equal_2;

  Equal_2     equal_points = ch_traits.equal_2_object();

  if (first == last) return result;

  std::vector< Point_2 >       P;      // points in subsets
  std::vector< Point_2 >       H;      // right endpoints of subproblems
  P.reserve(16);
  H.reserve(16);
  std::vector< PointIterator > L;      // start of subset range
  std::vector< PointIterator > R;      // end of subset range
  L.reserve(16);
  R.reserve(16);
  PointIterator           l;
  PointIterator           r;
  PointIterator           Pbegin, Pend;

  P.push_back(Point_2() );
  std::copy(first,last,std::back_inserter(P));
  P.push_back(Point_2() );
  Pbegin = std::next(P.begin());
  Pend   = std::prev(P.end());
  ch_we_point(Pbegin, Pend, l, r, ch_traits);
  Point_2 a = *l;
  Point_2 b = *r;
  if (equal_points(a,b))
  {
      *result = a;  ++result;
      return result;
  }
  #if defined(CGAL_CH_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS)
  OutputIterator  res(result);
  #else
  Tee_for_output_iterator<OutputIterator,Point_2> res(result);
  #endif // no postconditions ...
  H.push_back( a );
  L.push_back( Pbegin );
  Left_turn_2 left_turn = ch_traits.left_turn_2_object();
  R.push_back( l = std::partition( Pbegin, Pend,  [&left_turn,&a,&b](const Point_2&p)
                                                  { return left_turn(a, b, p); }));
  r = std::partition( l, Pend, [&left_turn,&a,&b](const Point_2&p)
                               { return left_turn(b, a, p); });

  Less_dist less_dist = ch_traits.less_signed_distance_to_line_2_object();
  for (;;)
  {
      if ( l != r)
      {
          if ( r-l > CGAL_ch_THRESHOLD )
          {
              Point_2 c = *std::min_element( l, r, [&less_dist,&a,&b](const Point_2&p1, const Point_2& p2)
                                                   { return less_dist(a, b, p1, p2); });
              H.push_back( b );
              L.push_back( l );
              R.push_back( l = std::partition(l, r, [&left_turn,&c,&b](const Point_2&p)
                                                    { return left_turn(b, c, p); }));
              r = std::partition(l, r, [&left_turn,&a,&c](const Point_2&p)
                                       { return left_turn(c, a, p); });
              b = c;
          }
          else
          {
              std::swap( a, *--l);
              std::swap( b, *++r);
              if ( ch_traits.less_xy_2_object()(*l,*r) )
              {
                std::sort(std::next(l), r,
                          ch_traits.less_xy_2_object() );
              }
              else
              {
                std::sort(std::next(l), r, [&ch_traits](const Point_2&p1, const Point_2& p2)
                                           { return ch_traits.less_xy_2_object()(p2, p1); });
              }
              ch__ref_graham_andrew_scan(l, std::next(r), res, ch_traits);
              std::swap( a, *l);
              std::swap( b, *r);
              if ( L.empty() ) break;
              a = b;
              b = H.back(); H.pop_back();
              l = L.back(); L.pop_back();
              r = R.back(); R.pop_back();
          }

      }
      else
      {
          *res = a;  ++res;
          if ( L.empty() ) break;
          a = b;
          b = H.back(); H.pop_back();
          l = L.back(); L.pop_back();
          r = R.back(); R.pop_back();
      }
  }
  CGAL_ch_postcondition( \
      is_ccw_strongly_convex_2( res.output_so_far_begin(), \
                                     res.output_so_far_end(), \
                                     ch_traits));
  CGAL_ch_expensive_postcondition( \
      ch_brute_force_check_2( \
          Pbegin, Pend, \
          res.output_so_far_begin(), res.output_so_far_end(), \
          ch_traits));
  #if defined(CGAL_CH_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS)
  return res;
  #else
  return res.to_output_iterator();
  #endif // no postconditions ...
}

} //namespace CGAL

#endif // CGAL_CH_BYKAT_C
