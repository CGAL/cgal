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

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_ISO_CUBOID_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_ISO_CUBOID_3_DO_INTERSECT_H

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
inline bool
do_intersect(const typename K::Iso_cuboid_3& icub1,
             const typename K::Iso_cuboid_3& icub2,
             const K&)
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

  for(int dim = 0; dim < DIM; ++dim)
  {
    min_idx[dim] = min_points[0].cartesian(dim) >= min_points[1].cartesian(dim) ? 0 : 1;
    max_idx[dim] = max_points[0].cartesian(dim) <= max_points[1].cartesian(dim) ? 0 : 1;
    if(min_idx[dim] != max_idx[dim] &&
       max_points[max_idx[dim]].cartesian(dim) < min_points[min_idx[dim]].cartesian(dim))
      return false;
  }

  return true;
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_ISO_CUBOID_3_DO_INTERSECT_H
