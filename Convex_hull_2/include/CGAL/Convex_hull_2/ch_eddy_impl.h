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

#ifndef CGAL_CH_EDDY_C
#define CGAL_CH_EDDY_C

#include <CGAL/license/Convex_hull_2.h>


#ifndef CGAL_CH_NO_POSTCONDITIONS
#include <CGAL/convexity_check_2.h>
#endif // CGAL_CH_NO_POSTCONDITIONS

#include <CGAL/assertions.h>
#include <CGAL/ch_selected_extreme_points_2.h>
#include <CGAL/algorithm.h>
#include <list>
#include <algorithm>

namespace CGAL {

template <class List, class ListIterator, class Traits>
void
ch__recursive_eddy(List& L,
                        ListIterator  a_it, ListIterator  b_it,
                        const Traits& ch_traits)
{
  typedef typename Traits::Point_2                            Point_2;

  typedef typename Traits::Compare_signed_distance_to_line_2  Compare_dist_2;
  typedef typename Traits::Less_xy_2                          Less_xy_2;
  typedef typename Traits::Left_turn_2                        Left_turn_2;

  Compare_dist_2 cmp_dist = ch_traits.compare_signed_distance_to_line_2_object();
  Left_turn_2 left_turn = ch_traits.left_turn_2_object();
  Less_xy_2 less_xy = ch_traits.less_xy_2_object();

  CGAL_precondition( \
    std::find_if(a_it, b_it, \
                 [&left_turn, a_it, b_it](const Point_2& p)
                 { return left_turn(*b_it, *a_it, p); }) \
    != b_it );

  const Point_2& a = *a_it;
  const Point_2& b = *b_it;

  ListIterator f_it = std::next(a_it);

  // We need the farthest point, but since we are on the right side of the line,
  // signed distances are negative. Hence std::min_element.
  auto less_dist = [&a, &b, &cmp_dist, &less_xy](const Point_2&p1, const Point_2& p2) -> bool
  {
    CGAL::Comparison_result res = cmp_dist(a, b, p1, p2);
    if(res == CGAL::EQUAL)
      return less_xy(p1, p2);

    return (res == CGAL::SMALLER);
  };

  ListIterator c_it = std::min_element( f_it, b_it, less_dist);
  Point_2 c = *c_it;

  c_it = std::partition(f_it, b_it, [&left_turn, &c, a_it](const Point_2& p)
                                    {return left_turn(c, *a_it, p);});
  f_it = std::partition(c_it, b_it, [&left_turn, &c, b_it](const Point_2& p)
                                    {return left_turn(*b_it, c, p);});
  c_it = L.insert(c_it, c);
  L.erase( f_it, b_it );

  if ( std::next(a_it) != c_it )
  {
      ch__recursive_eddy( L, a_it, c_it, ch_traits);
  }
  if ( std::next(c_it) != b_it )
  {
      ch__recursive_eddy( L, c_it, b_it, ch_traits);
  }
}

template <class InputIterator, class OutputIterator, class Traits>
OutputIterator
ch_eddy(InputIterator first, InputIterator last,
             OutputIterator  result,
             const Traits& ch_traits)
{
  typedef  typename Traits::Point_2                         Point_2;
  typedef  typename Traits::Left_turn_2                     Left_turn_2;
  typedef  typename Traits::Equal_2                         Equal_2;

  Left_turn_2 left_turn    = ch_traits.left_turn_2_object();
  Equal_2     equal_points = ch_traits.equal_2_object();

  if (first == last) return result;
  std::list< Point_2 >   L (first, last);

  typedef typename std::list< Point_2 >::iterator  list_iterator;
  list_iterator   w, e;
  ch_we_point(L.begin(), L.end(), w, e, ch_traits);
  Point_2 wp = *w;
  Point_2 ep = *e;
  if (equal_points(wp,ep) )
  {
      *result = wp;  ++result;
      return result;
  }

  L.erase(w);
  L.erase(e);

  e = std::partition(L.begin(), L.end(),
                     [&left_turn, &wp, &ep](const Point_2& p)
                     {return left_turn(ep, wp, p);} );
  L.push_front(wp);
  e = L.insert(e, ep);

  if ( std::next(L.begin()) != e )
  {
      ch__recursive_eddy( L, L.begin(), e, ch_traits);
  }
  w = std::find_if( e, L.end(), [&left_turn, &wp, &ep](const Point_2& p)
                                { return left_turn(wp, ep, p); });
  if ( w == L.end() )
  {
      L.erase( ++e, L.end() );
      return std::copy( L.begin(), L.end(), result );
  }
  w = L.insert(L.end(), wp);
  ch__recursive_eddy( L, e, w, ch_traits);


  CGAL_postcondition( \
      is_ccw_strongly_convex_2( L.begin(), w, ch_traits) );
  CGAL_expensive_postcondition( \
      ch_brute_force_check_2( first, last, \
                                   L.begin(), w, ch_traits ) );

  return std::copy( L.begin(), w, result );
}

} //namespace CGAL

#endif // CGAL_CH_EDDY_C
