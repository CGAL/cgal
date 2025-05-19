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

#ifndef CGAL_CH_JARVIS_C
#define CGAL_CH_JARVIS_C

#include <CGAL/license/Convex_hull_2.h>


#ifndef CGAL_CH_NO_POSTCONDITIONS
#include <CGAL/convexity_check_2.h>
#endif // CGAL_CH_NO_POSTCONDITIONS

#include <CGAL/IO/Tee_for_output_iterator.h>
#include <CGAL/assertions.h>
#include <CGAL/ch_selected_extreme_points_2.h>
#include <algorithm>

namespace CGAL {

template <class ForwardIterator, class OutputIterator,
          class Point, class Traits>
OutputIterator
ch_jarvis_march(ForwardIterator first, ForwardIterator last,
                const Point& start_p,
                const Point& stop_p,
                OutputIterator  result,
                const Traits& ch_traits)
{
  if (first == last) return result;
  typedef   typename Traits::Less_rotate_ccw_2     Less_rotate_ccw;
  typedef   typename Traits::Equal_2               Equal_2;

  Equal_2     equal_points = ch_traits.equal_2_object();

  #if defined(CGAL_CH_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS)
  OutputIterator  res(result);
  #else
  typedef   typename Traits::Point_2               Point_2;
  Tee_for_output_iterator<OutputIterator,Point_2> res(result);
  #endif // no postconditions ...
  CGAL_assertion_code( \
      int count_points = 0; )
  CGAL_assertion_code( \
      for (ForwardIterator fit = first; fit!= last; ++fit) ++count_points; )
  Less_rotate_ccw
      rotation_predicate = ch_traits.less_rotate_ccw_2_object( );
  *res = start_p;  ++res;
  CGAL_assertion_code( \
      int constructed_points = 1; )
  CGAL_exactness_assertion_code( \
      Point previous_point = start_p; )

  ForwardIterator it = std::min_element( first, last,
                                         [&start_p, &rotation_predicate](const Point& p1, const Point& p2)
                                         {return rotation_predicate(start_p, p1, p2);} );
  while (! equal_points(*it, stop_p) )
  {
      CGAL_exactness_assertion( \
          *it != previous_point );
      CGAL_exactness_assertion_code( \
          previous_point = *it; )

      *res = *it;  ++res;
      CGAL_assertion_code( \
          ++constructed_points;)
      CGAL_assertion( \
          constructed_points <= count_points + 1 );

      it = std::min_element( first, last,
                             [it, &rotation_predicate](const Point& p1, const Point& p2)
                             {return rotation_predicate(*it, p1, p2);} );
  }
  CGAL_postcondition( \
      is_ccw_strongly_convex_2( res.output_so_far_begin(), \
                                     res.output_so_far_end(), \
                                     ch_traits));
  CGAL_expensive_postcondition( \
      ch_brute_force_check_2(
          first, last, \
          res.output_so_far_begin(), res.output_so_far_end(), \
          ch_traits));
  #if defined(CGAL_CH_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS)
  return res;
  #else
  return res.to_output_iterator();
  #endif // no postconditions ...
}

template <class ForwardIterator, class OutputIterator, class Traits>
OutputIterator
ch_jarvis(ForwardIterator first, ForwardIterator last,
               OutputIterator  result,
               const Traits& ch_traits)
{
  if (first == last) return result;
  ForwardIterator start;
  ch_w_point(first, last, start, ch_traits);
  return ch_jarvis_march( first, last, *start, *start, result, ch_traits);
}

} //namespace CGAL

#endif // CGAL_CH_JARVIS_C
