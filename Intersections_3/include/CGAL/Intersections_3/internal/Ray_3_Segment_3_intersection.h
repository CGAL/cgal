// Copyright (c) 2003  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sébastien Loriot

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_RAY_3_SEGMENT_3_INTERSECTION_H
#define CGAL_INTERNAL_INTERSECTIONS_3_RAY_3_SEGMENT_3_INTERSECTION_H

#include <CGAL/Intersection_traits_3.h>
#include <CGAL/Intersections_3/internal/Point_3_Ray_3_do_intersect.h>
#include <CGAL/Intersections_3/internal/Line_3_Segment_3_intersection.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
typename Intersection_traits<K, typename K::Segment_3, typename K::Ray_3>::result_type
intersection(const typename K::Segment_3& s,
             const typename K::Ray_3& r,
             const K& k)
{
  CGAL_precondition(!s.is_degenerate() && !r.is_degenerate());

  typename Intersection_traits<K, typename K::Line_3, typename K::Segment_3>::result_type
      v = internal::intersection(r.supporting_line(), s, k);

  if(v)
  {
    if(const typename K::Point_3* p = intersect_get<typename K::Point_3>(v))
    {
      if(Ray_3_has_on_collinear_Point_3(r, *p, k))
        return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Ray_3>(*p);
    }
    else if(const typename K::Segment_3* s2 = intersect_get<typename K::Segment_3>(v))
    {
      bool has_source = Ray_3_has_on_collinear_Point_3(r, s.source(), k);
      bool has_target = Ray_3_has_on_collinear_Point_3(r, s.target(), k);
      if(has_source)
      {
        if(has_target)
        {
          return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Ray_3>(*s2);
        }
        else
        {
          if(k.equal_3_object()(r.source(), s.source()))
            return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Ray_3>(r.source());
          else
            return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Ray_3>(
                     k.construct_segment_3_object()(r.source(), s.source()));
        }
      }
      else
      {
        if(has_target)
        {
          if(k.equal_3_object()(r.source(), s.target()))
            return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Ray_3>(r.source());
          else
            return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Ray_3>(
                     k.construct_segment_3_object()(r.source(), s.target()));
        }
      }
    }
  }

  return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Ray_3>();
}

template <class K>
typename Intersection_traits<K, typename K::Ray_3, typename K::Segment_3>::result_type
intersection(const typename K::Ray_3& r,
             const typename K::Segment_3& s,
             const K& k)
{
  return intersection(s, r, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_RAY_3_SEGMENT_3_INTERSECTION_H
