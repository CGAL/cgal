// Copyright (c) 2008  INRIA Sophia-Antipolis (France), ETH Zurich (Switzerland).
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

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_LINE_3_INTERSECTION_H
#define CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_LINE_3_INTERSECTION_H

#include <CGAL/Intersection_traits_3.h>

#include <CGAL/kernel_assertions.h>
#include <CGAL/number_utils.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
typename Intersection_traits<K, typename K::Line_3, typename K::Iso_cuboid_3>::result_type
intersection(const typename K::Line_3 &line,
             const typename K::Iso_cuboid_3 &box,
             const K& k)
{
  typedef typename K::Point_3 Point_3;
  typedef typename K::Vector_3 Vector_3;
  typedef typename K::FT FT;

  bool all_values = true;
  FT _min = 0, _max = 0; // initialization to stop compiler warning
  const Point_3& _ref_point = line.point();
  const Vector_3& _dir = line.direction().vector();
  const Point_3& _iso_min = (box.min)();
  const Point_3& _iso_max = (box.max)();

  for(int i=0; i< _ref_point.dimension(); ++i) {
    if(_dir.homogeneous(i) == 0) {
      if(_ref_point.cartesian(i) < _iso_min.cartesian(i))
        return intersection_return<typename K::Intersect_3, typename K::Line_3, typename K::Iso_cuboid_3>();

      if(_ref_point.cartesian(i) > _iso_max.cartesian(i))
        return intersection_return<typename K::Intersect_3, typename K::Line_3, typename K::Iso_cuboid_3>();
    } else {
      FT newmin, newmax;
      if(_dir.homogeneous(i) > 0) {
        newmin = (_iso_min.cartesian(i) - _ref_point.cartesian(i)) / _dir.cartesian(i);
        newmax = (_iso_max.cartesian(i) - _ref_point.cartesian(i)) / _dir.cartesian(i);
      } else {
        newmin = (_iso_max.cartesian(i) - _ref_point.cartesian(i)) / _dir.cartesian(i);
        newmax = (_iso_min.cartesian(i) - _ref_point.cartesian(i)) / _dir.cartesian(i);
      }

      if(all_values) {
        _min = newmin;
        _max = newmax;
      } else {
        if(newmin > _min)
          _min = newmin;
        if(newmax < _max)
          _max = newmax;
        if(_max < _min)
          return intersection_return<typename K::Intersect_3, typename K::Line_3, typename K::Iso_cuboid_3>();
      }
      all_values = false;
    }
  }

  CGAL_kernel_assertion(!all_values);
  if(_max == _min)
    return intersection_return<typename K::Intersect_3, typename K::Line_3, typename K::Iso_cuboid_3>(
             k.construct_translated_point_3_object()(_ref_point,
               k.construct_scaled_vector_3_object()(_dir, _min)));

  return intersection_return<typename K::Intersect_3, typename K::Line_3, typename K::Iso_cuboid_3>(
           k.construct_segment_3_object()(
             k.construct_translated_point_3_object()(_ref_point,
               k.construct_scaled_vector_3_object()(_dir, _min)),
             k.construct_translated_point_3_object()(_ref_point,
               k.construct_scaled_vector_3_object()(_dir, _max))));
}

template <class K>
inline
typename Intersection_traits<K, typename K::Iso_cuboid_3, typename K::Line_3>::result_type
intersection(const typename K::Iso_cuboid_3& box,
             const typename K::Line_3& line,
             const K& k)
{
  return intersection(line, box, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_LINE_3_INTERSECTION_H
