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

#include <CGAL/algorithm.h>
#include <CGAL/assertions.h>
#include <CGAL/determinant.h>
#include <CGAL/number_utils.h>

#include <CGAL/Polygon_2/Polygon_2_simplicity.h>

#include <algorithm>
#include <cstdlib>
#include <iterator>
#include <limits>
#include <set>
#include <vector>

/// \cond SKIP_IN_MANUAL

namespace CGAL {

namespace internal {
namespace Polygon_2 {

// Filter a range of points to simplify sequences of collinear (or almost) points.
// A point is removed if the two segments, formed using its previous and next points
// in the range are collinear segments, up to a given tolerance.
//
// \tparam K must be a model of `Kernel`
// \tparam InputForwardIterator must be a model of `ForwardIterator`
//                              with value type `K::Point_2`
// \tparam OutputForwardIterator must be a model of `OutputIterator`
//                               with value type `K::Point_2`
//
// \param first, beyond the range
// \param out points that are not removed are output in `out`
// \param tolerance a tolerance on the collinearity of the two segments formed
//                  by three consecutive points of the range (more specifically,
//                  on the value of the determinant).
//
// \pre The range `(first, beyond)` is composed of at least three points.
// \pre Not all points in the range `(first, beyond)` are (almost) collinear.
template<typename K, typename InputForwardIterator, typename OutputForwardIterator>
OutputForwardIterator filter_collinear_points(InputForwardIterator first,
                                              InputForwardIterator beyond,
                                              OutputForwardIterator out,
                                              const typename K::FT tolerance =
                                                std::numeric_limits<typename K::FT>::epsilon())
{
  CGAL_precondition(std::distance(first, beyond) >= 3);

  typedef typename K::FT                              FT;
  typedef typename K::Point_2                         Point;

  InputForwardIterator last = std::prev(beyond);

  InputForwardIterator vit = first, vit_next = vit, vit_next_2 = vit, vend = vit;
  ++vit_next;
  ++(++vit_next_2);

  bool stop = false;

  do
  {
    CGAL_assertion(vit != vit_next);
    CGAL_assertion(vit_next != vit_next_2);
    CGAL_assertion(vit != vit_next_2);

    const Point& o = *vit;
    const Point& p = *vit_next;
    const Point& q = *vit_next_2;

    // Stop when 'p' is the starting point. It does not matter whether we are
    // in a collinear case or not.
    stop = (vit_next == vend);

    const FT det = CGAL::determinant(o.x() - q.x(), o.y() - q.y(),
                                     p.x() - q.x(), p.y() - q.y());

    if(CGAL::abs(det) <= tolerance)
    {
      // Only move 'p' and 'q' to ignore consecutive collinear points
      vit_next = (vit_next == last) ? first : ++vit_next;
      vit_next_2 = (vit_next_2 == last) ? first : ++vit_next_2;
    }
    else
    {
      // 'vit = vit_next' and not '++vit' because we don't necessarily have *(next(vit) == p)
      // and collinear points between 'o' and 'p' are ignored
      vit = vit_next;
      vit_next = (vit_next == last) ? first : ++vit_next;
      vit_next_2 = (vit_next_2 == last) ? first : ++vit_next_2;

      *out++ = p;
    }
  }
  while(!stop);

  return out;
}

} // namespace Polygon_2
} // namespace internal


//-----------------------------------------------------------------------//
//                          is_simple_2
//-----------------------------------------------------------------------//
// uses PolygonTraits::Less_xy_2
//      PolygonTraits::less_xy_2
//      PolygonTraits::Orientation_2
//      PolygonTraits::orientation_2
//      PolygonTraits::Point_2


template <class ForwardIterator, class PolygonTraits>
bool is_simple_2(ForwardIterator first,
                      ForwardIterator last,
                      const PolygonTraits& traits)
{
    if (first == last) return true;

    return is_simple_polygon(first, last, traits);
}

namespace internal { namespace Polygon_2 {

template <typename Traits>
class Compare_vertices {
    typedef typename Traits::Less_xy_2 Less_xy_2;
    typedef typename Traits::Point_2 Point_2;
    Less_xy_2 less;
public:
    Compare_vertices(Less_xy_2 less) : less(less) {}

    // `Point_like` derives from `Point_2`
    template <typename Point_like>
    bool operator()(const Point_like& p1, const Point_like& p2) {
        return less(Point_2(p1), Point_2(p2));
    }
}; // end Compare_vertices

} // end namespace Polygon_2
} // end namespace internal

//-----------------------------------------------------------------------//
//                          left_vertex_2
//-----------------------------------------------------------------------//
// uses PolygonTraits::Less_xy_2 and less_xy_2_object()

template <class ForwardIterator, class PolygonTraits>
ForwardIterator left_vertex_2(ForwardIterator first,
                                   ForwardIterator last,
                                   const PolygonTraits&traits)
{
    CGAL_polygon_precondition(first != last);
    internal::Polygon_2::Compare_vertices<PolygonTraits>
        less(traits.less_xy_2_object());
    return std::min_element(first, last, less);
}

//-----------------------------------------------------------------------//
//                          right_vertex_2
//-----------------------------------------------------------------------//
// uses PolygonTraits::Less_xy_2 and less_xy_2_object()

template <class ForwardIterator, class PolygonTraits>
ForwardIterator right_vertex_2(ForwardIterator first,
                                    ForwardIterator last,
                                    const PolygonTraits &traits)
{
    CGAL_polygon_precondition(first != last);
    internal::Polygon_2::Compare_vertices<PolygonTraits>
        less(traits.less_xy_2_object());
    return std::max_element(first, last, less);
}

//-----------------------------------------------------------------------//
//                          top_vertex_2
//-----------------------------------------------------------------------//
// uses PolygonTraits::Less_yx_2 and less_yx_2_object()

template <class ForwardIterator, class PolygonTraits>
ForwardIterator top_vertex_2(ForwardIterator first,
                                  ForwardIterator last,
                                  const PolygonTraits&traits)
{
    CGAL_polygon_precondition(first != last);
    return std::max_element(first, last, traits.less_yx_2_object());
}

//-----------------------------------------------------------------------//
//                          bottom_vertex_2
//-----------------------------------------------------------------------//
// uses PolygonTraits::Less_yx_2 and less_yx_2_object()

template <class ForwardIterator, class PolygonTraits>
ForwardIterator bottom_vertex_2(ForwardIterator first,
                                     ForwardIterator last,
                                     const PolygonTraits&traits)
{
    CGAL_polygon_precondition(first != last);
    return std::min_element(first, last, traits.less_yx_2_object());
}

//-----------------------------------------------------------------------//
//                          area_2
//-----------------------------------------------------------------------//
// uses Traits::
//  implemented in header file


//-----------------------------------------------------------------------//
//                          is_convex_2
//-----------------------------------------------------------------------//
// uses Traits::Less_xy_2 and less_xy_2_object()
//      Traits::Orientation_2 and orientation_2_object()
//      Traits::Equal_2 for filtering repeated points

template <class ForwardIterator, class Traits>
bool is_convex_2(ForwardIterator first,
                      ForwardIterator last,
                      const Traits& traits)
{
  ForwardIterator previous = first;
  if (previous == last) return true;

  ForwardIterator current = previous; ++current;
  if (current == last) return true;

  ForwardIterator next = current; ++next;
  if (next == last) return true;

  typename Traits::Equal_2 equal = traits.equal_2_object();

  while(equal(*previous, *current)) {
    current = next;
    ++next;
    if (next == last) return true;
  }

  typename Traits::Less_xy_2 less_xy_2 = traits.less_xy_2_object();
  typename Traits::Orientation_2 orientation = traits.orientation_2_object();
  // initialization
  bool HasClockwiseTriples = false;
  bool HasCounterClockwiseTriples = false;
  bool Order = less_xy_2(*previous, *current);
  int NumOrderChanges = 0;

  do {
  switch_orient:
    switch (orientation(*previous, *current, *next)) {
      case CLOCKWISE:
        HasClockwiseTriples = true;
        break;
      case COUNTERCLOCKWISE:
        HasCounterClockwiseTriples = true;
        break;
      case ZERO:
        if(equal(*current, *next)) {
          if(next == first) {
            first = current;
          }
          ++next;
          if (next == last)
            next = first;
          goto switch_orient;
        }
        break;
    }

    bool NewOrder = less_xy_2(*current, *next);
    if (Order != NewOrder) NumOrderChanges++;

    if (NumOrderChanges > 2) {
#ifdef CGAL_POLYGON_DEBUG
std::cout << "too many order changes: not convex!" << std::endl;
#endif
      return false;
    }

    if (HasClockwiseTriples && HasCounterClockwiseTriples) {
#ifdef CGAL_POLYGON_DEBUG
std::cout << "polygon not locally convex!" << std::endl;
#endif
      return false;
    }

    previous = current;
    current = next;
    ++next;
    if (next == last) next = first;
    Order = NewOrder;
  }
  while (previous != first);

  return true;
}

//-----------------------------------------------------------------------//
//                          oriented_side_2
//-----------------------------------------------------------------------//
// uses Traits::Less_xy_2
//      Traits::Compare_x_2 compare_x_2_object()
//      Traits::Compare_y_2 compare_y_2_object()
//      Traits::Orientation_2 and orientation_2_object()

template <class ForwardIterator, class Point, class Traits>
Oriented_side oriented_side_2(ForwardIterator first,
                                        ForwardIterator last,
                                        const Point& point,
                                        const Traits& traits)
{
  Orientation o = orientation_2(first, last, traits);
  CGAL_polygon_assertion(o != COLLINEAR);

  Bounded_side b = bounded_side_2(first, last, point, traits);
  switch (b) {
    case ON_BOUNDARY:
      return ON_ORIENTED_BOUNDARY;

    case ON_BOUNDED_SIDE:
      return (o == CLOCKWISE) ?  ON_NEGATIVE_SIDE : ON_POSITIVE_SIDE;

    default:
    //case ON_UNBOUNDED_SIDE:
      return (o == CLOCKWISE) ?  ON_POSITIVE_SIDE : ON_NEGATIVE_SIDE;
  }
}

//-----------------------------------------------------------------------//
//                          bounded_side_2
//-----------------------------------------------------------------------//
// uses Traits::Compare_x_2 compare_x_2_object()
//      Traits::Compare_y_2 compare_y_2_object()
//      Traits::Orientation_2 and orientation_2_object()
//
// returns ON_BOUNDED_SIDE, ON_BOUNDARY or ON_UNBOUNDED_SIDE

/*
   Implementation: we shoot a horizontal ray from the point to the right
   and count the number of intersections with polygon segments.
   If the number of intersections is odd, the point is inside.
   We don't count intersections with horizontal segments.
   With non-horizontal segments, the top vertex is considered to be part of
   the segment, but the bottom vertex is not. (Segments are half-closed).
*/

namespace i_polygon {

template <class Point, class Orientation_2, class CompareX_2>
int which_side_in_slab(Point const &point, Point const &low, Point const &high,
    Orientation_2 &orientation_2, CompareX_2 &compare_x_2)
// returns -1 if point is left of segment <low, high>, 0 if its on the segment
// and 1 if it is to the right
// precondition: low.y < point.y < high.y
{
    // first we try to decide on x coordinate values alone
    // This is an optimisation (whether this is really faster for
    // a homogeneous kernel is not clear, as comparisons can be expensive.
    Comparison_result low_x_comp_res = compare_x_2(point, low);
    Comparison_result high_x_comp_res = compare_x_2(point, high);
    if (low_x_comp_res == SMALLER) {
        if (high_x_comp_res == SMALLER)
            return -1;
    } else {
        switch (high_x_comp_res) {
          case LARGER: return 1;
          case SMALLER: break;
          case EQUAL: return (low_x_comp_res == EQUAL) ? 0 : 1;
        }
    }
    switch (orientation_2(low, point, high)) {
      case LEFT_TURN: return 1;
      case RIGHT_TURN: return -1;
      default: return 0;
    }
}

}  // end namespace i_polygon

template <class ForwardIterator, class Point, class PolygonTraits>
Bounded_side bounded_side_2(ForwardIterator first,
                                      ForwardIterator last,
                                      const Point& point,
                                      const PolygonTraits& traits)
{

  ForwardIterator current = first;
  if (current == last) return ON_UNBOUNDED_SIDE;

  ForwardIterator next = current; ++next;
  if (next == last) return ON_UNBOUNDED_SIDE;

  typename PolygonTraits::Compare_x_2 compare_x_2 = traits.compare_x_2_object();
  typename PolygonTraits::Compare_y_2 compare_y_2 = traits.compare_y_2_object();
  typename PolygonTraits::Orientation_2 orientation_2 = traits.orientation_2_object();
  bool IsInside = false;
  Comparison_result cur_y_comp_res = compare_y_2(*current, point);

  do // check if the segment (current,next) intersects
     // the ray { (t,point.y()) | t >= point.x() }
  {
    Comparison_result next_y_comp_res = compare_y_2(*next, point);

    switch (cur_y_comp_res) {
      case SMALLER:
        switch (next_y_comp_res) {
          case SMALLER:
            break;
          case EQUAL:
            switch (compare_x_2(point, *next)) {
              case SMALLER: IsInside = !IsInside; break;
              case EQUAL:   return ON_BOUNDARY;
              case LARGER:  break;
            }
            break;
          case LARGER:
            switch (i_polygon::which_side_in_slab(point, *current, *next,
                        orientation_2, compare_x_2)) {
              case -1: IsInside = !IsInside; break;
              case  0: return ON_BOUNDARY;
            }
            break;
        }
        break;
      case EQUAL:
        switch (next_y_comp_res) {
          case SMALLER:
            switch (compare_x_2(point, *current)) {
              case SMALLER: IsInside = !IsInside; break;
              case EQUAL:   return ON_BOUNDARY;
              case LARGER:  break;
            }
            break;
          case EQUAL:
            switch (compare_x_2(point, *current)) {
              case SMALLER:
                if (compare_x_2(point, *next) != SMALLER)
                    return ON_BOUNDARY;
                break;
              case EQUAL: return ON_BOUNDARY;
              case LARGER:
                if (compare_x_2(point, *next) != LARGER)
                    return ON_BOUNDARY;
                break;
            }
            break;
          case LARGER:
            if (compare_x_2(point, *current) == EQUAL) {
              return ON_BOUNDARY;
            }
            break;
        }
        break;
      case LARGER:
        switch (next_y_comp_res) {
          case SMALLER:
            switch (i_polygon::which_side_in_slab(point, *next, *current,
                        orientation_2, compare_x_2)) {
              case -1: IsInside = !IsInside; break;
              case  0: return ON_BOUNDARY;
            }
            break;
          case EQUAL:
            if (compare_x_2(point, *next) == EQUAL) {
              return ON_BOUNDARY;
            }
            break;
          case LARGER:
            break;
        }
        break;
    }

    current = next;
    cur_y_comp_res = next_y_comp_res;
    ++next;
    if (next == last) next = first;
  }
  while (current != first);

  return IsInside ? ON_BOUNDED_SIDE : ON_UNBOUNDED_SIDE;
}

//-----------------------------------------------------------------------//
//                          orientation_2
//-----------------------------------------------------------------------//
// uses Traits::Less_xy_2 (used by left_vertex_2)
//      Traits::orientation_2_object()

template <class ForwardIterator, class Traits>
Orientation orientation_2(ForwardIterator first,
                                    ForwardIterator last,
                                    const Traits& traits)
{
  CGAL_polygon_precondition(is_simple_2(first, last, traits));

  ForwardIterator i = left_vertex_2(first, last, traits);

  ForwardIterator prev = (i == first) ? last : i;
  --prev;

  ForwardIterator next = i;
  ++next;
  if (next == last)
    next = first;

  // if the range [first,last) contains less than three points, then some
  // of the points (prev,i,next) will coincide

  // return the orientation of the triple (prev,i,next)
  typedef typename Traits::Point_2 Point;
  return traits.orientation_2_object()(Point(*prev), Point(*i), Point(*next));
}

} //namespace CGAL

/// \endcond

// Local Variables:
// c-basic-offset: 4
// End:
