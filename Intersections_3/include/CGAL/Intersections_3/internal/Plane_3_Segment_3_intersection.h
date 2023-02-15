// Copyright (c) 2019 GeometryFactory(France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Geert-Jan Giezeman

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_PLANE_3_SEGMENT_3_INTERSECTION_H
#define CGAL_INTERNAL_INTERSECTIONS_3_PLANE_3_SEGMENT_3_INTERSECTION_H

#include <CGAL/Intersection_traits_3.h>
#include <CGAL/Intersections_3/internal/Line_3_Plane_3_intersection.h>

#include <CGAL/kernel_assertions.h>
#include <CGAL/enum.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
typename Intersection_traits<K, typename K::Plane_3, typename K::Segment_3>::result_type
intersection(const typename K::Plane_3& plane,
             const typename K::Segment_3& seg,
             const K& k)
{
  typedef typename K::Point_3 Point_3;

  const Point_3& source = seg.source();
  const Point_3& target = seg.target();

  CGAL::Oriented_side source_side = plane.oriented_side(source);
  CGAL::Oriented_side target_side = plane.oriented_side(target);

  switch(source_side)
  {
    case ON_ORIENTED_BOUNDARY:
      if (target_side == ON_ORIENTED_BOUNDARY)
        return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Segment_3>(seg);
      else
        return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Segment_3>(source);
    case ON_POSITIVE_SIDE:
      switch(target_side)
      {
        case ON_ORIENTED_BOUNDARY:
          return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Segment_3>(target);
        case ON_POSITIVE_SIDE:
          return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Segment_3>();
        default: // ON_NEGATIVE_SIDE:
        {
          // intersection object should be a point, but rounding errors
          // could lead to:
          // - a line: in such case, return seg,
          // - the empty set: return the empty set.
          typename Intersection_traits<K, typename K::Plane_3, typename K::Line_3>::result_type
              v = internal::intersection(plane, seg.supporting_line(), k);
          if(v)
          {
            if(const typename K::Point_3* p = intersect_get<typename K::Point_3>(v))
              return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Segment_3>(*p);
            else
              return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Segment_3>(seg);
          }
          else
          {
            return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Segment_3>();
          }
        }
      }
    case ON_NEGATIVE_SIDE:
      switch(target_side)
      {
        case ON_ORIENTED_BOUNDARY:
          return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Segment_3>(target);
        case ON_POSITIVE_SIDE:
        {
          // intersection object should be a point, but rounding errors
          // could lead to:
          // - a line: in such case, return seg,
          // - the empty set: return the empty set.
          typename Intersection_traits<K, typename K::Plane_3, typename K::Line_3>::result_type
              v = internal::intersection(plane, seg.supporting_line(), k);
          if(v) {
            if(const typename K::Point_3* p = intersect_get<typename K::Point_3>(v))
              return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Segment_3>(*p);
            else
              return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Segment_3>(seg);
          }
          else
            return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Segment_3>();
        }
        case ON_NEGATIVE_SIDE:
          return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Segment_3>();
      }
  }

  CGAL_kernel_assertion_msg(false, "Supposedly unreachable code.");
  return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Segment_3>();
}

template <class K>
inline
typename Intersection_traits<K, typename K::Segment_3, typename K::Plane_3>::result_type
intersection(const typename K::Segment_3& seg,
             const typename K::Plane_3& plane,
             const K& k)
{
  return intersection(plane, seg, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_PLANE_3_SEGMENT_3_INTERSECTION_H
