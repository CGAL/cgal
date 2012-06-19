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

#ifndef CGAL_CONVEX_HULL_2_H
#define CGAL_CONVEX_HULL_2_H

#include <CGAL/basic.h>
#include <CGAL/convex_hull_traits_2.h>
#include <CGAL/ch_akl_toussaint.h>
#include <CGAL/ch_bykat.h>
#include <iterator> 

namespace CGAL {

template <class InputIterator, class OutputIterator, class Traits>
inline
OutputIterator
CGAL_convex_hull_points_2(InputIterator first, InputIterator last,
                          OutputIterator  result,
                          const Traits& ch_traits,
                          std::input_iterator_tag )
{ return ch_bykat(first, last, result, ch_traits); }

template <class InputIterator, class OutputIterator, class Traits>
inline
OutputIterator
CGAL_convex_hull_points_2(InputIterator first, InputIterator last,
                          OutputIterator  result,
                          const Traits& ch_traits,
                          std::forward_iterator_tag )
{ return ch_akl_toussaint(first, last, result, ch_traits); }

template <class InputIterator, class OutputIterator, class Traits>
inline
OutputIterator
CGAL_convex_hull_points_2(InputIterator first, InputIterator last,
                          OutputIterator  result,
                          const Traits& ch_traits,
                          std::bidirectional_iterator_tag )
{ return ch_akl_toussaint(first, last, result, ch_traits); }

template <class InputIterator, class OutputIterator, class Traits>
inline
OutputIterator
CGAL_convex_hull_points_2(InputIterator first, InputIterator last,
                          OutputIterator  result,
                          const Traits& ch_traits,
                          std::random_access_iterator_tag )
{ return ch_akl_toussaint(first, last, result, ch_traits); }


template <class InputIterator, class OutputIterator, class Traits>
inline
OutputIterator
convex_hull_points_2(InputIterator first, InputIterator last,
                     OutputIterator  result,
                     const Traits& ch_traits)
{
    typedef std::iterator_traits<InputIterator>   ITraits;
    typedef typename ITraits::iterator_category   Category;
    return CGAL_convex_hull_points_2(first, last, result, ch_traits,
                                     Category());
}


template <class ForwardIterator, class OutputIterator>
inline
OutputIterator 
convex_hull_points_2(ForwardIterator first, ForwardIterator last, 
                     OutputIterator  result )
{ 
    typedef std::iterator_traits<ForwardIterator> ITraits;
    typedef typename ITraits::value_type          value_type;
    typedef typename ITraits::iterator_category   Category;
    typedef CGAL::Kernel_traits<value_type>       KTraits;
    typedef typename KTraits::Kernel              Kernel;
    return CGAL_convex_hull_points_2(first, last, result, Kernel(), 
                                     Category());
}


// generates the counterclockwise sequence of extreme points
// of the points in the range [|first|,|last|). The resulting sequence
// is placed starting at position |result|, and the past-the-end iterator
// for the resulting sequence is returned. It is not specified, at which
// point the cyclic sequence of extreme points is cut into a linear
// sequence.
// {\it Preconditions:}
// [|first|,|last|) does not contain |result|.
// {\sc traits}: operates on |Traits::Point_2| using |Traits::Less_xy_2|, 
// |Traits::Equal_2|, |Traits::Less_yx_2|, and |Traits::Left_turn_2|.
template <class InputIterator, class OutputIterator, class Traits>
inline
OutputIterator
convex_hull_2(InputIterator first, InputIterator last,
              OutputIterator  result, const Traits& ch_traits)
{
    return convex_hull_points_2(first, last, result, ch_traits);
}

template <class ForwardIterator, class OutputIterator>
inline
OutputIterator 
convex_hull_2(ForwardIterator first, ForwardIterator last, 
              OutputIterator  result )
{
    return convex_hull_points_2(first, last, result);
}


// generates the counterclockwise sequence of extreme points
// on the lower hull of the points in the range [|first|,|last|). 
// The resulting sequence is placed starting at position |result|, 
// and the past-the-end iterator for the resulting sequence is returned. 
// The sequence starts with the leftmost point, the rightmost point is
// not included.
// {\it Preconditions:}
// [|first|,|last|) does not contain |result|.
// {\sc traits}: operates on |Traits::Point_2| using |Traits::Less_xy_2|
// |Traits::Equal_2| and |Traits::Left_turn_2|.
template <class InputIterator, class OutputIterator, class Traits>
inline
OutputIterator
lower_hull_points_2(InputIterator first, InputIterator last,
                    OutputIterator  result,
                    const Traits& ch_traits)
{ return ch_lower_hull_scan(first, last, result, ch_traits); }

template <class ForwardIterator, class OutputIterator>
inline
OutputIterator 
lower_hull_points_2(ForwardIterator first, ForwardIterator last, 
                    OutputIterator  result )
{ 
    typedef std::iterator_traits<ForwardIterator> ITraits;
    typedef typename ITraits::value_type          value_type;
    typedef CGAL::Kernel_traits<value_type>       KTraits;
    typedef typename KTraits::Kernel              Kernel;
    return lower_hull_points_2(first, last, result, Kernel());
}


// generates the counterclockwise sequence of extreme points
// on the upper hull of the points in the range [|first|,|last|). 
// The resulting sequence is placed starting at position |result|, 
// and the past-the-end iterator for the resulting sequence is returned. 
// The sequence starts with the rightmost point, the leftmost point is
// not included.
// {\it Preconditions:}
// [|first|,|last|) does not contain |result|.
// {\sc traits}: operates on |Traits::Point_2| using |Traits::Less_xy_2|,
// |Traits::Equal_2| and |Traits::Left_turn_2|.
template <class InputIterator, class OutputIterator, class Traits>
inline
OutputIterator
upper_hull_points_2(InputIterator first, InputIterator last,
                    OutputIterator  result,
                    const Traits& ch_traits)
{ return ch_upper_hull_scan(first, last, result, ch_traits); }


template <class ForwardIterator, class OutputIterator>
inline
OutputIterator 
upper_hull_points_2(ForwardIterator first, ForwardIterator last, 
                    OutputIterator  result )
{ 
    typedef std::iterator_traits<ForwardIterator> ITraits;
    typedef typename ITraits::value_type          value_type;
    typedef CGAL::Kernel_traits<value_type>       KTraits;
    typedef typename KTraits::Kernel              Kernel;
    return upper_hull_points_2(first, last, result, Kernel());
}

} //namespace CGAL

#endif // CGAL_CONVEX_HULL_2_H
