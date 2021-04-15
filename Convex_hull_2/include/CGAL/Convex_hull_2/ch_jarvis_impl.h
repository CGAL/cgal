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
#include <CGAL/Convex_hull_2/ch_assertions.h>
#include <CGAL/ch_selected_extreme_points_2.h>
#include <algorithm>
#include <boost/bind.hpp>

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
  using namespace boost;

  if (first == last) return result;
  typedef   typename Traits::Less_rotate_ccw_2     Less_rotate_ccw;
  typedef   typename Traits::Equal_2               Equal_2;

  Equal_2     equal_points = ch_traits.equal_2_object();

  #if defined(CGAL_CH_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS) \
    || defined(NDEBUG)
  OutputIterator  res(result);
  #else
  typedef   typename Traits::Point_2               Point_2;
  Tee_for_output_iterator<OutputIterator,Point_2> res(result);
  #endif // no postconditions ...
  CGAL_ch_assertion_code( \
      int count_points = 0; )
  CGAL_ch_assertion_code( \
      for (ForwardIterator fit = first; fit!= last; ++fit) ++count_points; )
  Less_rotate_ccw
      rotation_predicate = ch_traits.less_rotate_ccw_2_object( );
  *res = start_p;  ++res;
  CGAL_ch_assertion_code( \
      int constructed_points = 1; )
  CGAL_ch_exactness_assertion_code( \
      Point previous_point = start_p; )

  ForwardIterator it = std::min_element( first, last,
                                         boost::bind(rotation_predicate, boost::cref(start_p), _1, _2) );
  while (! equal_points(*it, stop_p) )
  {
      CGAL_ch_exactness_assertion( \
          *it != previous_point );
      CGAL_ch_exactness_assertion_code( \
          previous_point = *it; )

      *res = *it;  ++res;
      CGAL_ch_assertion_code( \
          ++constructed_points;)
      CGAL_ch_assertion( \
          constructed_points <= count_points + 1 );

      it = std::min_element( first, last,
                             boost::bind(rotation_predicate, *it, _1, _2) );
  }
  CGAL_ch_postcondition( \
      is_ccw_strongly_convex_2( res.output_so_far_begin(), \
                                     res.output_so_far_end(), \
                                     ch_traits));
  CGAL_ch_expensive_postcondition( \
      ch_brute_force_check_2(
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
