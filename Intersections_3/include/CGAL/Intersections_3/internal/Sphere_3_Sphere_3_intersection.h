// Copyright (c) 2001 GeometryFactory(France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Pedro Machado Manhaes

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_SPHERE_3_SPHERE_3_INTERSECTION_H
#define CGAL_INTERNAL_INTERSECTIONS_3_SPHERE_3_SPHERE_3_INTERSECTION_H

#include <CGAL/Intersection_traits_3.h>
#include <CGAL/Intersections_3/internal/Plane_3_Sphere_3_intersection.h>

#include <CGAL/number_utils.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
inline
typename Intersection_traits<K, typename K::Sphere_3, typename K::Sphere_3>::result_type
intersection(const typename K::Sphere_3& s1,
             const typename K::Sphere_3& s2,
             const K& k)
{
  typedef typename K::Plane_3 Plane_3;

  typename K::Construct_center_3 center = k.construct_center_3_object();
  typename K::Compute_squared_radius_3 sqr = k.compute_squared_radius_3_object();

  if(k.equal_3_object()(center(s1), center(s2)))
  {
    if(sqr(s1) == sqr(s2))
    {
      if(is_zero(sqr(s1)))
        return intersection_return<typename K::Intersect_3, typename K::Sphere_3, typename K::Sphere_3>(s1.center());
      else
        return intersection_return<typename K::Intersect_3, typename K::Sphere_3, typename K::Sphere_3>(s1);
    }
    else
    {
      // cocentrics
      return intersection_return<typename K::Intersect_3, typename K::Sphere_3, typename K::Sphere_3>();
    }
  }

  Plane_3 pl = k.construct_radical_plane_3_object()(s1,s2);
  typename Intersection_traits<K, typename K::Sphere_3, typename K::Plane_3>::result_type v = intersection(pl, s1, k);

  if(v)
  {
    if(const typename K::Point_3* p = intersect_get<typename K::Point_3>(v))
      return intersection_return<typename K::Intersect_3, typename K::Sphere_3, typename K::Sphere_3>(*p);
    else if(const typename K::Circle_3* c = intersect_get<typename K::Circle_3>(v))
      return intersection_return<typename K::Intersect_3, typename K::Sphere_3, typename K::Sphere_3>(*c);
  }

  return intersection_return<typename K::Intersect_3, typename K::Sphere_3, typename K::Sphere_3>();
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_SPHERE_3_SPHERE_3_INTERSECTION_H
