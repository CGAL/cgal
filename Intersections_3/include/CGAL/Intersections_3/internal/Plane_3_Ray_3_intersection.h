// Copyright (c) 1997-2010
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
// Author(s)     : Geert-Jan Giezeman <geert@cs.uu.nl>

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_PLANE_3_RAY_3_INTERSECTION_H
#define CGAL_INTERNAL_INTERSECTIONS_3_PLANE_3_RAY_3_INTERSECTION_H

#include <CGAL/Intersection_traits_3.h>
#include <CGAL/Intersections_3/internal/Line_3_Plane_3_intersection.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
typename Intersection_traits<K, typename K::Plane_3, typename K::Ray_3>::result_type
intersection(const typename K::Plane_3& plane,
             const typename K::Ray_3& ray,
             const K& k)
{
  typedef typename K::Point_3 Point_3;

  typename Intersection_traits<K, typename K::Plane_3, typename K::Line_3>::result_type
      v = internal::intersection(plane, ray.supporting_line(), k);

  if(v)
  {
    if(const Point_3* p = intersect_get<Point_3>(v))
    {
      if(ray.collinear_has_on(*p))
        return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Ray_3>(*p);
      else
        return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Ray_3>();
    }
  }
  else
  {
    return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Ray_3>();
  }

  return intersection_return<typename K::Intersect_3, typename K::Plane_3, typename K::Ray_3>(ray);
}

template <class K>
inline
typename Intersection_traits<K, typename K::Ray_3, typename K::Plane_3>::result_type
intersection(const typename K::Ray_3& ray,
             const typename K::Plane_3& plane,
             const K& k)
{
  return intersection(plane, ray, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_PLANE_3_RAY_3_INTERSECTION_H
