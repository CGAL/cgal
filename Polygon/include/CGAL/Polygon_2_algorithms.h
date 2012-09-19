// Copyright (c) 1997  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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

/*!
  \file Polygon_2_algorithms.h
 */

#ifndef CGAL_POLYGON_2_ALGORITHMS_H
#define CGAL_POLYGON_2_ALGORITHMS_H

#include <CGAL/basic.h>
#include <CGAL/enum.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Polygon_2/polygon_assertions.h>

///
namespace CGAL {

//-----------------------------------------------------------------------//
//                  algorithms for sequences of 2D points
//-----------------------------------------------------------------------//


/// \addtogroup PkgPolygon2Functions
/// @{

/// Returns an iterator to the leftmost point from the range
/// `[first,last)`. In case of a tie, the point
/// with the smallest `y`-coordinate is taken.
///
/// \requires `Traits` is a model of the concept `PolygonTraits_2`.
/// 	  In fact, only the members `Less_xy_2` and
/// 	  `less_xy_2_object()` are used.
/// \requires The value type of `ForwardIterator` must be `Traits::Point_2`,
///
/// \sa `CGAL::right_vertex_2()`
/// \sa `CGAL::top_vertex_2()`
/// \sa `CGAL::bottom_vertex_2()`
/// \sa `CGAL::Polygon_2`
template <class ForwardIterator, class Traits>
ForwardIterator left_vertex_2(ForwardIterator first,
			      ForwardIterator last,
			      const Traits& traits);

/// Returns an iterator to the rightmost point from the range
/// `[first,last)`. In case of a tie, the point
/// with the largest `y`-coordinate is taken.
/// 
/// \requires `Traits` is a model of the concept 
/// 	  `PolygonTraits_2`
/// 	  In fact, only the members `Less_xy_2` and
/// 	  `less_xy_2_object()` are used.
/// \requires The value type of `ForwardIterator` must be `Traits::Point_2`,
/// 
/// 
/// \sa `CGAL::left_vertex_2()`
/// \sa `CGAL::top_vertex_2`
/// \sa `CGAL::bottom_vertex_2`
/// \sa `CGAL::Polygon_2`
template <class ForwardIterator, class Traits>
ForwardIterator right_vertex_2(ForwardIterator first,
			       ForwardIterator last,
			       const Traits& traits);

/// Returns an iterator to the topmost point from the range
/// `[first,last)`. In case of a tie, the point
/// with the largest `x`-coordinate is taken.
///
/// \requires `Traits` is a model of the concept 
/// 	  `PolygonTraits_2`
/// 	  In fact, only the members `Less_yx_2` and
/// 	  `less_yx_2_object()` are used.
/// \requires The value type of `ForwardIterator` must be `Traits::Point_2`,
/// 
/// \sa `CGAL::left_vertex_2()`
/// \sa `CGAL::right_vertex_2()`
/// \sa `CGAL::bottom_vertex_2()`
/// \sa `CGAL::Polygon_2`
template <class ForwardIterator, class Traits>
ForwardIterator top_vertex_2(ForwardIterator first,
			     ForwardIterator last,
			     const Traits& traits);

/// Returns an iterator to the bottommost point from the range
/// `[first,last)`. In case of a tie, the point
/// with the smallest `x`-coordinate is taken.
///
/// \requires `Traits` is a model of the concept 
/// 	  `PolygonTraits_2`
/// 	  In fact, only the members `Less_yx_2` and
/// 	  `less_yx_2_object()` are used.
/// \requires The value type of `ForwardIterator` must be `Traits::Point_2`,
/// 
/// \ref CGAL::left_vertex_2
/// \sa `CGAL::left_vertex_2()`
/// \sa `CGAL::right_vertex_2()`
/// \sa `CGAL::top_vertex_2()`
/// \sa `CGAL::Polygon_2`
/// \sa `PolygonTraits_2`
template <class ForwardIterator, class Traits>
ForwardIterator bottom_vertex_2(ForwardIterator first,
				ForwardIterator last,
				const Traits& traits);

/// Returns the bounding box of the range `[first,last)`. 
/// 
/// \requires `Traits` is a model of the concept 
/// 	  `PolygonTraits_2`
/// 	  In fact, only the members `Construct_bbox_2` and
/// 	  `construct_bbox_2_object()` are used.
/// \requires The value type of `InputIterator` must be `Traits::Point_2`.
/// 
/// \sa `CGAL::Polygon_2`
template <class InputIterator, class Traits>
Bbox_2 bbox_2(InputIterator first, 
	      InputIterator last,
	      const Traits& traits);


/// Computes the signed area of the polygon defined by the range of points
/// `[first,last)`. The area is returned in the parameter
/// `result`. The sign is positive for counterclockwise polygons, negative for
/// clockwise polygons. If the polygon is not simple, the area is not well defined.
/// The functionality is also available by the `polygon_area_2` function, which
/// returns the area instead of taking it as a parameter.
/// 
/// \requires `Traits` is a model of the concept 
/// 	  `PolygonTraits_2
/// 	  Only the following members of this traits class are used:
///   - `Compute_area_2` : Computes the signed area of the
/// 	    oriented triangle defined by 3 `Point_2` passed as arguments.
///   - `FT`
///   - `compute_area_2_object()`
/// \requires The value type of `ForwardIterator` must be `Traits::Point_2`,
/// 
/// \sa `CGAL::polygon_area_2`
/// \sa `PolygonTraits_2`
/// \sa `CGAL::orientation_2`
/// \sa `CGAL::Polygon_2`
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

/// Computes the signed area of the polygon defined by the range of points
/// `[first,last)`.
/// The sign is positive for counterclockwise polygons, negative for
/// clockwise polygons. If the polygon is not simple, the area is not well defined.
/// 
/// \requires `Traits` is a model of the concept `PolygonTraits_2`. Only the following members of this traits class are used:
///   - `Compute_area_2` : Computes the signed area of the
/// 	    oriented triangle defined by 3 `Point_2` passed as arguments.
///   - `FT`
///   - `compute_area_2_object`
/// \requires `ForwardIterator::value_type` should be `Traits::Point_2`,
/// 
/// 
/// \sa `PolygonTraits_2 `
/// \sa `CGAL::orientation_2 `
/// \sa `CGAL::Polygon_2 `
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

/// The function is_convex_2 computes if a polygon is convex.
/// 
/// \requires `Traits` is a model of the concept 
/// 	  `PolygonTraits_2`
/// 	  Only the following members of this traits class are used:
///   - `Less_xy_2`
///   - `Orientation_2`
///   - `less_xy_2_object`
///   - `orientation_2_object`
/// \requires `ForwardIterator::value_type` should be `Traits::Point_2`,
///
/// \sa `PolygonTraits_2 `
/// \sa `CGAL::Polygon_2 `
template <class ForwardIterator, class Traits>
bool is_convex_2(ForwardIterator first,
		 ForwardIterator last,
		 const Traits& traits);

/// The function is_simple_2 computes if a polygon is simple, that is, if the edges 
/// do not intersect, except consecutive edges in their common vertex.
/// 
/// Checks if the polygon defined by the
/// iterator range `[first,last)` is simple.
/// 
/// \requires `Traits` is a model of the concept 
/// 	  `PolygonTraits_2`
/// 	  Only the following members of this traits class are used:
///   - `Point_2`
///   - `Less_xy_2`
///   - `Orientation_2`
///   - `less_xy_2_object()`
///   - `orientation_2_object()`
/// \requires The value type of `ForwardIterator` must be `Traits::Point_2`,
/// 
/// Implementation
/// --------------
/// 
/// The simplicity test is implemented by means of a plane sweep algorithm.
/// The algorithm is quite robust when used with inexact number types.
/// The running time is O(n log n), where n is the number of vertices of the
/// polygon.
/// 
/// \sa `PolygonTraits_2 `
/// \sa `CGAL::Polygon_2 `
template <class ForwardIterator, class Traits>
bool is_simple_2(ForwardIterator first,
		 ForwardIterator last,
		 const Traits& traits);

// In the following two functions we would like to use Traits::Point_2
// instead of Point, but this is not allowed by g++ 2.7.2.
///
/// The function oriented_side_2 computes on which side of a polygon a point lies.
/// \requires `Traits` is a model of the concept 
/// 	  `PolygonTraits_2`
/// 	  Only the following members of this traits class are used:
///   - `Less_xy_2`
///   - `Compare_x_2`
///   - `Compare_y_2`
///   - `Orientation_2`
///   - `less_xy_2_object()`
///   - `compare_x_2_object()`
///   - `compare_y_2_object()`
///   - `orientation_2_object()`
/// \requires The value type of `ForwardIterator` must be `Traits::Point_2`,
///
/// \sa `PolygonTraits_2 `
/// \sa `CGAL::bounded_side_2 `
/// \sa `CGAL::is_simple_2 `
/// \sa `CGAL::Polygon_2 `
/// \sa `Oriented_side`
template <class ForwardIterator, class Point, class Traits>
Oriented_side oriented_side_2(ForwardIterator first,
			      ForwardIterator last,
			      const Point& point,
			      const Traits& traits);

/// The function  bounded_side_2 computes if a point lies inside a polygon.
/// The function bounded_side_2 computes if a point lies inside a polygon. 
/// The polygon is defined by the sequence of points `[first,last)`.
/// Being inside is defined by the odd-even rule. If we take a ray starting at the
/// point and extending to infinity (in any direction), we count the number of
/// intersections. If this number is odd, the point is inside, otherwise it is
/// outside. If the point is on a polygon edge, a special value is returned.  A
/// simple polygon divides the plane in an unbounded and a bounded region.
/// According to the definition points in the bounded region are inside the polygon.
/// 
/// 
/// \requires `Traits` is a model of the concept 
/// 	  `PolygonTraits_2`
/// 	  Only the following members of this traits class are used:
///   - `Compare_x_2`
///   - `Compare_y_2`
///   - `Orientation_2`
///   - `compare_x_2_object()`
///   - `compare_y_2_object()`
///   - `orientation_2_object()`
/// \requires The value type of  `ForwardIterator` must be `Traits::Point_2`,
/// 
/// Implementation
/// --------------
/// 
/// The running time is linear in the number of vertices of the polygon.
/// A horizontal ray is taken to count the number of intersections.
/// Special care is taken that the result is correct even if there are degeneracies
/// (if the ray passes through a vertex).
/// 
/// \sa `PolygonTraits_2 `
/// \sa `CGAL::oriented_side_2 `
/// \sa `CGAL::Polygon_2 `
/// \sa `CGAL::Bounded_side`
template <class ForwardIterator, class Point, class Traits>
Bounded_side bounded_side_2(ForwardIterator first,
			    ForwardIterator last,
			    const Point& point,
			    const Traits& traits);

/// The function orientation_2 computes if a polygon is clockwise or counterclockwise
/// oriented.
/// \pre `is_simple_2(first, last, traits);`
/// 
/// \requires `Traits` is a model of the concept 
/// 	  `PolygonTraits_2`
/// 	  Only the following members of this traits class are used:
///   - `Less_xy_2`
///   - `less_xy_2_object()`
///   - `orientation_2_object()`
/// \requires The value type of `ForwardIterator` must be `Traits::Point_2`,
/// 
/// 
/// 
/// \sa `PolygonTraits_2`
/// \sa `CGAL::is_simple_2`
/// \sa `CGAL::Polygon_2`
/// \sa `CGAL::Orientation`
template <class ForwardIterator, class Traits>
Orientation orientation_2(ForwardIterator first,
			  ForwardIterator last,
			  const Traits& traits);

/// @}

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
