// Copyright (c) 1998  Max-Planck-Institute Saarbruecken (Germany).
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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Stefan Schirra 

#ifndef CGAL_SEGMENT_INTERSECTION_POINTS_2_H
#define CGAL_SEGMENT_INTERSECTION_POINTS_2_H

#include <CGAL/license/Convex_hull_2.h>


#include <CGAL/basic.h>
#include <CGAL/algorithm.h>
#include <iterator>

namespace CGAL {

/*
#include <CGAL/Segment_2_Segment_2_intersection.h>
template <class ForwardIterator, class OutputIterator, class R>
OutputIterator
si_brute_force(ForwardIterator first, ForwardIterator last,
                    OutputIterator result,
                    const R& )
{
  typedef Point_2<R>         Point;
  ForwardIterator inner, outer;
  Object i_obj;
  Point p;
  for ( outer = first; outer != last; ++outer)
      for ( inner = successor(outer); inner != last; ++inner)
      {
          i_obj = intersection( *outer, *inner);
          if ( assign( p, i_obj) )
              result++ = p;
      }
  return result;
}
*/

template <class ForwardIterator, class OutputIterator, class Traits>
OutputIterator
si_brute_force_II(ForwardIterator first, ForwardIterator last,
                  OutputIterator result,
                  const Traits& traits)
{
  typedef typename Traits::Point_2         Point;
  typedef typename Traits::Line_2          Line;
  typedef typename Traits::Orientation_2   Orientation;
  Orientation orientation = traits.orientation_2_object();

  for ( ForwardIterator outer = first; outer != last; ++outer)
    for ( ForwardIterator inner = cpp11::next(outer); inner != last; ++inner)
      {
          Point s1 = (*outer).source();
          Point e1 = (*outer).target();
          Point s2 = (*inner).source();
          Point e2 = (*inner).target();
          if ( (orientation( s1, e1, s2) != orientation( s1, e1, e2))
             &&(orientation( s2, e2, s1) != orientation( s2, e2, e1)))
          {
              Line l1( s1, e1);
              Line l2( s2, e2);
              result++ =  Point( l1.b()*l2.c() - l2.b()*l1.c(),
                                 l2.a()*l1.c() - l1.a()*l2.c(),
                                 l1.a()*l2.b() - l2.a()*l1.b());
          }
      }
  return result;
}

template <class ForwardIterator, class OutputIterator, class Traits>
OutputIterator
segment_intersection_points_2(ForwardIterator first, ForwardIterator last,
                              OutputIterator result,
                              const Traits& traits)
{ return si_brute_force_II( first, last, result, traits); }

template <class ForwardIterator, class OutputIterator, class Traits>
OutputIterator
segment_intersection_points_2(ForwardIterator first, ForwardIterator last,
                              OutputIterator result)
{
    typedef std::iterator_traits<ForwardIterator> ITraits;
    typedef typename ITraits::value_type          value_type;
    typedef CGAL::Kernel_traits<value_type>       KTraits;
    typedef typename KTraits::Kernel              Kernel;
    return segment_intersection_points_2( first, last, result, Kernel()); 
}

} // namespace CGAL

#endif // CGAL_SEGMENT_INTERSECTION_POINTS_2_H
