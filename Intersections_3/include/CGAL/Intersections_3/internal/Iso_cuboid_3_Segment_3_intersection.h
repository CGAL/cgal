// Copyright (c) 2008 ETH Zurich (Switzerland)
// Copyright (c) 2008-2009 INRIA Sophia-Antipolis (France)
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

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_SEGMENT_3_INTERSECTION_H
#define CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_SEGMENT_3_INTERSECTION_H

#include <CGAL/Intersection_traits_3.h>

#include <CGAL/kernel_assertions.h>
#include <CGAL/number_utils.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
typename Intersection_traits<K, typename K::Segment_3, typename K::Iso_cuboid_3>::result_type
intersection(const typename K::Segment_3& seg,
             const typename K::Iso_cuboid_3& box,
             const K&)
{
  typedef typename K::Point_3 Point_3;
  typedef typename K::Vector_3 Vector_3;
  typedef typename K::Segment_3 Segment_3;
  typedef typename K::FT FT;

  FT _min = 0, _max;
  const Point_3& _ref_point = seg.source();
  const Vector_3& _dir = seg.direction().vector();
  const Point_3& _iso_min = (box.min)();
  const Point_3& _iso_max = (box.max)();
  int main_dir = (CGAL_NTS abs(_dir.x()) > CGAL_NTS abs(_dir.y()) )
                   ? (CGAL_NTS abs(_dir.x()) > CGAL_NTS abs(_dir.z()) ? 0 : 2)
                   : (CGAL_NTS abs(_dir.y()) > CGAL_NTS abs(_dir.z()) ? 1 : 2);
  _max = (seg.target().cartesian(main_dir)-_ref_point.cartesian(main_dir)) / _dir.cartesian(main_dir);

  for(int i=0; i< _ref_point.dimension(); ++i) {
    if(_dir.homogeneous(i) == 0) {
      if(_ref_point.cartesian(i) < _iso_min.cartesian(i)) {
        return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Iso_cuboid_3>();
      }
      if(_ref_point.cartesian(i) > _iso_max.cartesian(i)) {
        return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Iso_cuboid_3>();
      }
    } else {
      FT newmin, newmax;
      if(_dir.homogeneous(i) > 0) {
        newmin = (_iso_min.cartesian(i) - _ref_point.cartesian(i)) / _dir.cartesian(i);
        newmax = (_iso_max.cartesian(i) - _ref_point.cartesian(i)) / _dir.cartesian(i);
      } else {
        newmin = (_iso_max.cartesian(i) - _ref_point.cartesian(i)) / _dir.cartesian(i);
        newmax = (_iso_min.cartesian(i) - _ref_point.cartesian(i)) / _dir.cartesian(i);
      }

      if(newmax < _max)
        _max = newmax;
      if(newmin > _min)
        _min = newmin;
      if(_max < _min)
        return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Iso_cuboid_3>();
    }
  }

  if(_max == _min)
    return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Iso_cuboid_3>(
             Point_3{_ref_point + _dir * _min });

  return intersection_return<typename K::Intersect_3, typename K::Segment_3, typename K::Iso_cuboid_3>(
           Segment_3{_ref_point + _dir*_min, _ref_point + _dir*_max});
}

template <class K>
inline
typename Intersection_traits<K, typename K::Iso_cuboid_3, typename K::Segment_3>::result_type
intersection(const typename K::Iso_cuboid_3& box,
             const typename K::Segment_3& seg,
             const K& k)
{
  return intersection(seg, box, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_SEGMENT_3_INTERSECTION_H
