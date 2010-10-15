// Copyright (c) 1997  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Wieger Wesselink <wieger@cs.ruu.nl>

#ifndef CGAL_POLYGON_2_ALGORITHMS_H
#define CGAL_POLYGON_2_ALGORITHMS_H

#include <CGAL/basic.h>
#include <CGAL/enum.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Polygon_2/polygon_assertions.h>

namespace CGAL {

//-----------------------------------------------------------------------//
//                  algorithms for sequences of 2D points
//-----------------------------------------------------------------------//

template <class ForwardIterator, class Traits>
ForwardIterator left_vertex_2(ForwardIterator first,
			      ForwardIterator last,
			      const Traits& traits);

template <class ForwardIterator, class Traits>
ForwardIterator right_vertex_2(ForwardIterator first,
			       ForwardIterator last,
			       const Traits& traits);

template <class ForwardIterator, class Traits>
ForwardIterator top_vertex_2(ForwardIterator first,
			     ForwardIterator last,
			     const Traits& traits);

template <class ForwardIterator, class Traits>
ForwardIterator bottom_vertex_2(ForwardIterator first,
				ForwardIterator last,
				const Traits& traits);

template <class InputIterator, class Traits>
Bbox_2 bbox_2(InputIterator first, 
	      InputIterator last,
	      const Traits& traits);


template <class ForwardIterator, class Traits>
void 
area_2( ForwardIterator first, ForwardIterator last,
   	typename Traits::FT &result,
        const Traits& traits)
{
   typedef typename Traits::FT FT;
   result = FT(0);
   // check if the polygon is empty
   if (first == last) return;
   ForwardIterator second = first; ++second;
   // check if the polygon has only one point
   if (second == last) return;
   typename Traits::Compute_area_2 compute_area_2 =
            traits.compute_area_2_object();
   ForwardIterator third = second;
   while (++third != last) {
	result = result + compute_area_2(*first, *second, *third);
	second = third;
   }
}

template <class ForwardIterator, class Traits>
typename Traits::FT 
polygon_area_2( ForwardIterator first, ForwardIterator last,
		const Traits& traits)
{
   typedef typename Traits::FT FT;
   FT result = FT(0);
   // check if the polygon is empty
   if (first == last) return result;
   ForwardIterator second = first; ++second;
   // check if the polygon has only one point
   if (second == last) return result;
   typename Traits::Compute_area_2 compute_area_2 =
            traits.compute_area_2_object();
   ForwardIterator third = second;
   while (++third != last) {
	result = result + compute_area_2(*first, *second, *third);
	second = third;
   }
   return result;
}

template <class ForwardIterator, class Traits>
bool is_convex_2(ForwardIterator first,
		 ForwardIterator last,
		 const Traits& traits);

template <class ForwardIterator, class Traits>
bool is_simple_2(ForwardIterator first,
		 ForwardIterator last,
		 const Traits& traits);

// In the following two functions we would like to use Traits::Point_2 instead
// of Point, but this is not allowed by g++ 2.7.2.

template <class ForwardIterator, class Point, class Traits>
Oriented_side oriented_side_2(ForwardIterator first,
			      ForwardIterator last,
			      const Point& point,
			      const Traits& traits);

template <class ForwardIterator, class Point, class Traits>
Bounded_side bounded_side_2(ForwardIterator first,
			    ForwardIterator last,
			    const Point& point,
			    const Traits& traits);

template <class ForwardIterator, class Traits>
Orientation orientation_2(ForwardIterator first,
			  ForwardIterator last,
			  const Traits& traits);

//-----------------------------------------------------------------------//
//                         implementation
//-----------------------------------------------------------------------//

template <class ForwardIterator>
inline
ForwardIterator left_vertex_2(ForwardIterator first,
			      ForwardIterator last)
{  
  typedef typename Kernel_traits<
    typename std::iterator_traits<ForwardIterator>::value_type>::Kernel K; 
  return left_vertex_2(first, last, K());
}


template <class ForwardIterator>
inline
ForwardIterator right_vertex_2(ForwardIterator first,
			       ForwardIterator last)
{ 
  typedef typename Kernel_traits<
    typename std::iterator_traits<ForwardIterator>::value_type>::Kernel K; 
  return right_vertex_2(first, last, K());
}


template <class ForwardIterator>
inline
ForwardIterator top_vertex_2(ForwardIterator first,
			     ForwardIterator last)
{
  typedef typename Kernel_traits<
    typename std::iterator_traits<ForwardIterator>::value_type>::Kernel K; 
  return top_vertex_2(first, last, K());
}


template <class ForwardIterator>
inline
ForwardIterator bottom_vertex_2(ForwardIterator first,
				ForwardIterator last)
{
  typedef typename Kernel_traits<
    typename std::iterator_traits<ForwardIterator>::value_type>::Kernel K; 
  return bottom_vertex_2(first, last, K());
}

template <class InputIterator>
inline
Bbox_2 bbox_2(InputIterator first,
	      InputIterator last)
{
  typedef typename Kernel_traits<
    typename std::iterator_traits<InputIterator>::value_type>::Kernel K; 
  return bbox_2(first, last, K());
}


template <class ForwardIterator, class Numbertype>
inline
void area_2(ForwardIterator first,
	    ForwardIterator last,
	    Numbertype& result)
{
  typedef typename Kernel_traits<
    typename std::iterator_traits<ForwardIterator>::value_type>::Kernel K; 
  area_2(first, last, result, K());
}



template <class ForwardIterator>
inline
bool is_convex_2(ForwardIterator first,
		 ForwardIterator last)
{
  typedef typename Kernel_traits<
    typename std::iterator_traits<ForwardIterator>::value_type>::Kernel K; 
  return is_convex_2(first, last, K());
}



template <class ForwardIterator>
inline
bool is_simple_2(ForwardIterator first,
		 ForwardIterator last)
{
  typedef typename Kernel_traits<
    typename std::iterator_traits<ForwardIterator>::value_type>::Kernel K; 
  return is_simple_2(first, last, K());
}

template <class ForwardIterator>
inline
Oriented_side oriented_side_2(
  ForwardIterator first,
  ForwardIterator last,
  const typename std::iterator_traits<ForwardIterator>::value_type& point)
{
  typedef typename Kernel_traits<
    typename std::iterator_traits<ForwardIterator>::value_type>::Kernel K; 
  return oriented_side_2(first, last, point, K());
}


template <class ForwardIterator>
inline
Bounded_side bounded_side_2(
  ForwardIterator first,
  ForwardIterator last,
  const typename std::iterator_traits<ForwardIterator>::value_type& point)
{
  typedef typename Kernel_traits<
    typename std::iterator_traits<ForwardIterator>::value_type>::Kernel K; 
  return bounded_side_2(first, last, point, K());
}



template <class ForwardIterator>
inline
Orientation orientation_2(ForwardIterator first,
			  ForwardIterator last)
{
  typedef typename Kernel_traits<
    typename std::iterator_traits<ForwardIterator>::value_type>::Kernel K; 
  return orientation_2(first, last, K());
}

} //namespace CGAL

#include <CGAL/Polygon_2/Polygon_2_algorithms_impl.h>

#endif // CGAL_POLYGON_2_ALGORITHMS_H
