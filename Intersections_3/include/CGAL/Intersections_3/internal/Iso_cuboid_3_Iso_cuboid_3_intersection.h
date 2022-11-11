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
// Author(s)     : Andreas Fabri

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_ISO_CUBOID_3_INTERSECTION_H
#define CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_ISO_CUBOID_3_INTERSECTION_H

#include <CGAL/Intersection_traits_3.h>

#include <CGAL/enum.h>
#include <CGAL/wmult.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
typename Intersection_traits<K, typename K::Iso_cuboid_3, typename K::Iso_cuboid_3>::result_type
intersection(const typename K::Iso_cuboid_3& icub1,
             const typename K::Iso_cuboid_3& icub2,
             const K& k)
{
  typedef typename K::Point_3 Point_3;

  Point_3 min_points[2];
  Point_3 max_points[2];
  min_points[0] = (icub1.min)();
  min_points[1] = (icub2.min)();
  max_points[0] = (icub1.max)();
  max_points[1] = (icub2.max)();
  const int DIM = 3;
  int min_idx[DIM];
  int max_idx[DIM];

  Point_3 newmin;
  Point_3 newmax;
  for(int dim = 0; dim < DIM; ++dim)
  {
    min_idx[dim] = min_points[0].cartesian(dim) >= min_points[1].cartesian(dim) ? 0 : 1;
    max_idx[dim] = max_points[0].cartesian(dim) <= max_points[1].cartesian(dim) ? 0 : 1;
    if(min_idx[dim] != max_idx[dim] &&
       max_points[max_idx[dim]].cartesian(dim) < min_points[min_idx[dim]].cartesian(dim))
      return intersection_return<typename K::Intersect_3, typename K::Iso_cuboid_3, typename K::Iso_cuboid_3>();
  }

  if(min_idx[0] == min_idx[1] && min_idx[0] == min_idx[2])
  {
    newmin = min_points[min_idx[0]];
  }
  else
  {
    newmin = Point_3(min_idx[0] == 0 ? wmult_hw((K*)0, min_points[0].hx(), min_points[1])
                                     : wmult_hw((K*)0, min_points[1].hx(), min_points[0]),
                     min_idx[1] == 0 ? wmult_hw((K*)0, min_points[0].hy(), min_points[1])
                                     : wmult_hw((K*)0, min_points[1].hy(), min_points[0]),
                     min_idx[2] == 0 ? wmult_hw((K*)0, min_points[0].hz(), min_points[1])
                                     : wmult_hw((K*)0, min_points[1].hz(), min_points[0]),
                     wmult_hw((K*)0, min_points[0].hw(), min_points[1]));
  }

  if(max_idx[0] == max_idx[1] && max_idx[0] == max_idx[2])
  {
    newmax = max_points[max_idx[0]];
  }
  else
  {
    newmax = Point_3(max_idx[0] == 0 ? wmult_hw((K*)0, max_points[0].hx(), max_points[1])
                                     : wmult_hw((K*)0, max_points[1].hx(), max_points[0]),
                     max_idx[1] == 0 ? wmult_hw((K*)0, max_points[0].hy(), max_points[1])
                                     : wmult_hw((K*)0, max_points[1].hy(), max_points[0]),
                     max_idx[2] == 0 ? wmult_hw((K*)0, max_points[0].hz(), max_points[1])
                                     : wmult_hw((K*)0, max_points[1].hz(), max_points[0]),
                     wmult_hw((K*)0, max_points[0].hw(), max_points[1]));
  }

  return intersection_return<typename K::Intersect_3, typename K::Iso_cuboid_3, typename K::Iso_cuboid_3>(
           k.construct_iso_cuboid_3_object()(newmin, newmax));
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_ISO_CUBOID_3_INTERSECTION_H
