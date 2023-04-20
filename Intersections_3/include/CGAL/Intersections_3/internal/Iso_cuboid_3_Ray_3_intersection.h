// Copyright (c) 2001 ETH Zurich (Switzerland)
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

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_RAY_3_INTERSECTION_H
#define CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_RAY_3_INTERSECTION_H

#include <CGAL/Intersection_traits_3.h>

#include <CGAL/kernel_assertions.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
typename Intersection_traits<K, typename K::Ray_3, typename K::Iso_cuboid_3>::result_type
intersection(const typename K::Ray_3& ray,
             const typename K::Iso_cuboid_3& box,
             const K&)
{
  typedef typename K::Point_3 Point_3;
  typedef typename K::Vector_3 Vector_3;
  typedef typename K::Segment_3 Segment_3;
  typedef typename K::FT FT;

  bool all_values = true;
  FT _min = 0, _max = 0; // initialization to prevent compiler warning
  const Point_3& _ref_point = ray.source();
  const Vector_3& _dir = ray.direction().vector();
  const Point_3& _iso_min = (box.min)();
  const Point_3& _iso_max = (box.max)();

  for (int i=0; i< _ref_point.dimension(); ++i) {
    if (_dir.homogeneous(i) == 0) {
      if (_ref_point.cartesian(i) < _iso_min.cartesian(i)) {
        return intersection_return<typename K::Intersect_3, typename K::Ray_3, typename K::Iso_cuboid_3>();
      }
      if (_ref_point.cartesian(i) > _iso_max.cartesian(i)) {
        return intersection_return<typename K::Intersect_3, typename K::Ray_3, typename K::Iso_cuboid_3>();
      }
    } else {
      FT newmin, newmax;
      if (_dir.homogeneous(i) > 0) {
        newmin = (_iso_min.cartesian(i) - _ref_point.cartesian(i)) / _dir.cartesian(i);
        newmax = (_iso_max.cartesian(i) - _ref_point.cartesian(i)) / _dir.cartesian(i);
      } else {
        newmin = (_iso_max.cartesian(i) - _ref_point.cartesian(i)) / _dir.cartesian(i);
        newmax = (_iso_min.cartesian(i) - _ref_point.cartesian(i)) / _dir.cartesian(i);
      }
      if (all_values) {
        _max = newmax;
      } else {
        if (newmax < _max)
          _max = newmax;
      }
      if (newmin > _min)
        _min = newmin;
      if (_max < _min)
        return intersection_return<typename K::Intersect_3, typename K::Ray_3, typename K::Iso_cuboid_3>();
      all_values = false;
    }
  }
  CGAL_kernel_assertion(!all_values);
  if (_max == _min) {
    return intersection_return<typename K::Intersect_3, typename K::Ray_3, typename K::Iso_cuboid_3>(Point_3(_ref_point + _dir * _min ));
  }
  return intersection_return<typename K::Intersect_3, typename K::Ray_3, typename K::Iso_cuboid_3>(
        Segment_3(_ref_point + _dir*_min, _ref_point + _dir*_max));
}

template <class K>
inline
typename Intersection_traits<K, typename K::Iso_cuboid_3, typename K::Ray_3>::result_type
intersection(const typename K::Iso_cuboid_3& box,
             const typename K::Ray_3& ray,
             const K& k)
{
  return intersection(ray, box, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_RAY_3_INTERSECTION_H
