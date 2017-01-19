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


#ifndef CGAL_CH_JARVIS_H
#define CGAL_CH_JARVIS_H

#include <CGAL/license/Convex_hull_2.h>


#include <CGAL/basic.h>
#include <iterator>


namespace CGAL {

// generates the counterclockwise ordered subsequence of
// extreme points between |start_p| and |stop_p| of the points in the
// range [|first|,|last|), starting at position result with point |start_p|.
// The last point generated is the point preceding |stop_p| in the
// counterclockwise order of extreme points.
// {\it Precondition:} |start_p| and |stop_p| are extreme points with respect
// to the points in the range [|first|,|last|) and |stop_p| is an element of
// range [|first|,|last|).
// {\sc traits}: uses |Traits::Point_2| $\equiv$ |Point|, |Traits::Equal_2|
// and |Traits::Less_rotate_ccw_2|.
template <class ForwardIterator, class OutputIterator, 
          class Point, class Traits>
OutputIterator
ch_jarvis_march(ForwardIterator first, ForwardIterator last,
                const Point& start_p, const Point& stop_p,
                OutputIterator  result,
                const Traits& ch_traits);

template <class ForwardIterator, class OutputIterator, class Point>
inline
OutputIterator
ch_jarvis_march(ForwardIterator first, ForwardIterator last,
                const Point& start_p,
                const Point& stop_p,
                OutputIterator  result )
{
    typedef CGAL::Kernel_traits<Point>  KTraits;
    typedef typename KTraits::Kernel    Kernel;
    return ch_jarvis_march( first, last, start_p, stop_p, result, Kernel());
}


// same as |convex_hull_2(first,last,result)|.
// {\sc traits}: uses |Traits::Point_2|, |Traits::Less_rotate_ccw_2|,
// |Traits::Equal_2| and |Traits::Less_xy_2|.
template <class ForwardIterator, class OutputIterator, class Traits>
OutputIterator
ch_jarvis(ForwardIterator first, ForwardIterator last, 
               OutputIterator  result,
               const Traits& ch_traits);
template <class ForwardIterator, class OutputIterator>
inline
OutputIterator
ch_jarvis(ForwardIterator first, ForwardIterator last, OutputIterator  result)
{ 
    typedef std::iterator_traits<ForwardIterator> ITraits;
    typedef typename ITraits::value_type          value_type;
    typedef CGAL::Kernel_traits<value_type>       KTraits;
    typedef typename KTraits::Kernel              Kernel;
    return ch_jarvis( first, last, result, Kernel()); 
}



} //namespace CGAL

#include <CGAL/Convex_hull_2/ch_jarvis_impl.h>

#endif // CGAL_CH_JARVIS_H
