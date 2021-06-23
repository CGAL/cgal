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

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_LINE_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_LINE_3_DO_INTERSECT_H

#include <CGAL/kernel_assertions.h>
#include <CGAL/number_utils.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
inline bool
do_intersect(const typename K::Line_3& line,
             const typename K::Iso_cuboid_3& box,
             const K&)
{
  typedef typename K::Point_3 Point_3;
  typedef typename K::Vector_3 Vector_3;
  typedef typename K::FT FT;

  bool all_values = true;
  FT _min = 0, _max = 0; // initialization to stop compiler warning
  FT _denum;
  const Point_3& _ref_point = line.point();
  const Vector_3& _dir = line.direction().vector();
  const Point_3& _iso_min = (box.min)();
  const Point_3& _iso_max = (box.max)();

  for(int i=0; i< _ref_point.dimension(); ++i) {
    if(_dir.homogeneous(i) == 0) {
      if(_ref_point.cartesian(i) < _iso_min.cartesian(i))
        return false;

      if(_ref_point.cartesian(i) > _iso_max.cartesian(i))
        return false;
    } else {
      FT newmin, newmax;
      FT newdenum = _dir.cartesian(i);
      if(_dir.homogeneous(i) > 0) {
        newmin = (_iso_min.cartesian(i) - _ref_point.cartesian(i));
        newmax = (_iso_max.cartesian(i) - _ref_point.cartesian(i));
      } else {
        newmin = (_iso_max.cartesian(i) - _ref_point.cartesian(i));
        newmax = (_iso_min.cartesian(i) - _ref_point.cartesian(i));
      }

      if(all_values) {
        _min = newmin;
        _max = newmax;
        _denum = newdenum;
      } else {
        if(compare_quotients(newmin, newdenum, _min, _denum) == LARGER)
          _min = newmin;
        if(compare_quotients(newmax, newdenum, _max, _denum) == LARGER)
          _max = newmax;
        if(compare_quotients(_max, _denum, _min, _denum) == SMALLER)
          return false;

        _denum = newdenum;

      }
      all_values = false;
    }
  }

  CGAL_kernel_assertion(!all_values);
  return true;
}

template <class K>
inline bool
do_intersect(const typename K::Iso_cuboid_3& ic,
             const typename K::Line_3& l,
             const K& k)
{
  return do_intersect(l, ic, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_LINE_3_DO_INTERSECT_H
