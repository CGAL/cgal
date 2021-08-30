// Copyright (c) 1997
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Wieger Wesselink <wieger@cs.ruu.nl>

/*!
  \file Polygon_2_algorithms.h
*/

#ifndef CGAL_POLYGON_2_ALGORITHMS_H
#define CGAL_POLYGON_2_ALGORITHMS_H

#include <CGAL/config.h>
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
/// \tparam Traits is a model of the concept `PolygonTraits_2`.
///           Only the members `Less_xy_2` and
///           `less_xy_2_object()` are used.
/// \tparam ForwardIterator must have `Traits::Point_2` as value type.
///
///
/// \sa `CGAL::right_vertex_2()`
/// \sa `CGAL::top_vertex_2()`
/// \sa `CGAL::bottom_vertex_2()`
/// \sa `CGAL::Polygon_2`
template <class ForwardIterator, class PolygonTraits>
ForwardIterator left_vertex_2(ForwardIterator first,
                              ForwardIterator last,
                              const PolygonTraits& traits);

/// Returns an iterator to the rightmost point from the range
/// `[first,last)`. In case of a tie, the point
/// with the largest `y`-coordinate is taken.
///
/// \tparam Traits is a model of the concept
///           `PolygonTraits_2`.
///           In fact, only the members `Less_xy_2` and
///           `less_xy_2_object()` are used.
/// \tparam ForwardIterator must have`Traits::Point_2` as value type.
///
///
/// \sa `CGAL::left_vertex_2()`
/// \sa `CGAL::top_vertex_2()`
/// \sa `CGAL::bottom_vertex_2()`
/// \sa `CGAL::Polygon_2`
template <class ForwardIterator, class PolygonTraits>
ForwardIterator right_vertex_2(ForwardIterator first,
                               ForwardIterator last,
                               const PolygonTraits& traits);

/// Returns an iterator to the topmost point from the range
/// `[first,last)`. In case of a tie, the point
/// with the largest `x`-coordinate is taken.
///
/// \tparam Traits is a model of the concept
///           `PolygonTraits_2`.
///           Only the members `Less_yx_2` and
///           `less_yx_2_object()` are used.
/// \tparam ForwardIterator must have `Traits::Point_2` as value type.
///
/// \sa `CGAL::left_vertex_2()`
/// \sa `CGAL::right_vertex_2()`
/// \sa `CGAL::bottom_vertex_2()`
/// \sa `CGAL::Polygon_2`
template <class ForwardIterator, class PolygonTraits>
ForwardIterator top_vertex_2(ForwardIterator first,
                             ForwardIterator last,
                             const PolygonTraits& traits);

/// Returns an iterator to the bottommost point from the range
/// `[first,last)`. In case of a tie, the point
/// with the smallest `x`-coordinate is taken.
///
/// \tparam Traits is a model of the concept
///           `PolygonTraits_2`.
///           Only the members `Less_yx_2` and
///           `less_yx_2_object()` are used.
/// \tparam ForwardIterator must have `Traits::Point_2` as value type.
///
/// \sa `CGAL::left_vertex_2()`
/// \sa `CGAL::right_vertex_2()`
/// \sa `CGAL::top_vertex_2()`
/// \sa `CGAL::Polygon_2`
/// \sa `PolygonTraits_2`
template <class ForwardIterator, class PolygonTraits>
ForwardIterator bottom_vertex_2(ForwardIterator first,
                                ForwardIterator last,
                                const PolygonTraits& traits);

/// Computes the signed area of the polygon defined by the range of points
/// `[first,last)`. The area is returned in the parameter
/// `result`. The sign is positive for counterclockwise polygons, negative for
/// clockwise polygons. If the polygon is not simple, the area is not well defined.
/// The functionality is also available by the `polygon_area_2()` function, which
/// returns the area instead of taking it as a parameter.
///
/// \tparam Traits is a model of the concept
///           `PolygonTraits_2`.
///           Only the following members of this traits class are used:
///   - `Compute_area_2` : Computes the signed area of the
///             oriented triangle defined by 3 `Point_2` passed as arguments.
///   - `FT`
///   - `compute_area_2_object()`
/// \tparam ForwardIterator must have `Traits::Point_2` as value type.
///
/// \sa `CGAL::polygon_area_2()`
/// \sa `PolygonTraits_2`
/// \sa `CGAL::orientation_2()`
/// \sa `CGAL::Polygon_2`
template <class ForwardIterator, class PolygonTraits>
void
area_2( ForwardIterator first, ForwardIterator last,
           typename PolygonTraits::FT &result,
        const PolygonTraits& traits)
{
  typedef typename PolygonTraits::FT FT;
   result = FT(0);
   // check if the polygon is empty
   if (first == last) return;
   ForwardIterator second = first; ++second;
   // check if the polygon has only one point
   if (second == last) return;
   typename PolygonTraits::Compute_area_2 compute_area_2 =
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
/// \tparam Traits is a model of the concept `PolygonTraits_2`. Only the following members of this traits class are used:
///   - `Compute_area_2` : Computes the signed area of the
///             oriented triangle defined by 3 `Point_2` passed as arguments.
///   - `FT`
///   - `compute_area_2_object`
/// \tparam ForwardIterator must have `Traits::Point_2` as value type.
///
///
/// \sa `PolygonTraits_2 `
/// \sa `CGAL::orientation_2()`
/// \sa `CGAL::Polygon_2 `
template <class ForwardIterator, class PolygonTraits>
typename PolygonTraits::FT
polygon_area_2( ForwardIterator first, ForwardIterator last,
                const PolygonTraits& traits)
{
   typedef typename PolygonTraits::FT FT;
   FT result = FT(0);
   // check if the polygon is empty
   if (first == last) return result;
   ForwardIterator second = first; ++second;
   // check if the polygon has only one point
   if (second == last) return result;
   typename PolygonTraits::Compute_area_2 compute_area_2 =
            traits.compute_area_2_object();
   ForwardIterator third = second;
   while (++third != last) {
        result = result + compute_area_2(*first, *second, *third);
        second = third;
   }
   return result;
}

/// Checks if the polygon is convex.
///
/// \tparam Traits is a model of the concept
///           `PolygonTraits_2`.
///           Only the following members of this traits class are used:
///   - `Less_xy_2`
///   - `Orientation_2`
///   - `less_xy_2_object`
///   - `orientation_2_object`
/// \tparam ForwardIterator must have `PolygonTraits::Point_2` as value type.
///
/// \sa `PolygonTraits_2 `
/// \sa `CGAL::Polygon_2 `
template <class ForwardIterator, class PolygonTraits>
bool is_convex_2(ForwardIterator first,
                 ForwardIterator last,
                 const PolygonTraits& traits);

/// Checks if the polygon defined by the
/// iterator range `[first,last)` is simple, that is, if the edges
/// do not intersect, except consecutive edges in their common vertex.
///
/// \tparam Traits is a model of the concept
///           `PolygonTraits_2`.
///           Only the following members of this traits class are used:
///   - `Point_2`
///   - `Less_xy_2`
///   - `Orientation_2`
///   - `less_xy_2_object()`
///   - `orientation_2_object()`
/// \tparam ForwardIterator must have `PolygonTraits::Point_2` as value type.
///
/// \cgalHeading{Implementation}
///
/// The simplicity test is implemented by means of a plane sweep algorithm.
/// The algorithm is quite robust when used with inexact number types.
/// The running time is `O(n log n)`, where n is the number of vertices of the
/// polygon.
///
/// \sa `PolygonTraits_2`
/// \sa `CGAL::Polygon_2`
template <class ForwardIterator, class PolygonTraits>
bool is_simple_2(ForwardIterator first,
                 ForwardIterator last,
                 const PolygonTraits& traits);

// In the following two functions we would like to use Traits::Point_2
// instead of Point, but this is not allowed by g++ 2.7.2.
///
/// Computes on which side of a polygon a point lies.
/// \tparam Traits is a model of the concept
///           `PolygonTraits_2`.
///           Only the following members of this traits class are used:
///   - `Less_xy_2`
///   - `Compare_x_2`
///   - `Compare_y_2`
///   - `Orientation_2`
///   - `less_xy_2_object()`
///   - `compare_x_2_object()`
///   - `compare_y_2_object()`
///   - `orientation_2_object()`
/// \tparam ForwardIterator must have `PolygonTraits::Point_2` as value type.
///
/// \sa `PolygonTraits_2`
/// \sa `CGAL::bounded_side_2()`
/// \sa `CGAL::is_simple_2()`
/// \sa `CGAL::Polygon_2`
/// \sa `Oriented_side`
template <class ForwardIterator, class Point, class Traits>
Oriented_side oriented_side_2(ForwardIterator first,
                              ForwardIterator last,
                              const Point& point,
                              const Traits& traits);

/// Computes if a point lies inside a polygon.
/// The polygon is defined by the sequence of points `[first,last)`.
/// Being inside is defined by the odd-even rule. If we take a ray starting at the
/// point and extending to infinity (in any direction), we count the number of
/// intersections. If this number is odd, the point is inside, otherwise it is
/// outside. If the point is on a polygon edge, a special value is returned.  A
/// simple polygon divides the plane in an unbounded and a bounded region.
/// According to the definition points in the bounded region are inside the polygon.
///
///
/// \tparam Traits is a model of the concept
///           `PolygonTraits_2`.
///           Only the following members of this traits class are used:
///   - `Compare_x_2`
///   - `Compare_y_2`
///   - `Orientation_2`
///   - `compare_x_2_object()`
///   - `compare_y_2_object()`
///   - `orientation_2_object()`
/// \tparam ForwardIterator must have `Traits::Point_2` as value type.
///
/// \cgalHeading{Implementation}
///
/// The running time is linear in the number of vertices of the polygon.
/// A horizontal ray is taken to count the number of intersections.
/// Special care is taken that the result is correct even if there are degeneracies
/// (if the ray passes through a vertex).
///
///
/// \sa `PolygonTraits_2`
/// \sa `CGAL::oriented_side_2()`
/// \sa `CGAL::Polygon_2 `
/// \sa `CGAL::Bounded_side`
template <class ForwardIterator, class Point, class PolygonTraits>
Bounded_side bounded_side_2(ForwardIterator first,
                            ForwardIterator last,
                            const Point& point,
                            const PolygonTraits& traits);

/// Computes if a polygon is clockwise or counterclockwise oriented.
/// \pre `is_simple_2(first, last, traits);`
///
/// \tparam Traits is a model of the concept
///           `PolygonTraits_2`.
///           Only the following members of this traits class are used:
///   - `Less_xy_2`
///   - `less_xy_2_object()`
///   - `orientation_2_object()`
/// \tparam ForwardIterator must have`Traits::Point_2` as value type.
///
///
///
/// \sa `PolygonTraits_2`
/// \sa `CGAL::is_simple_2()`
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
