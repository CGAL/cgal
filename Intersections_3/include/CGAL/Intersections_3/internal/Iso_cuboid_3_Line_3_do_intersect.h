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

#include <CGAL/Intersections_3/internal/Bbox_3_Line_3_do_intersect.h>

#include <CGAL/kernel_assertions.h>
#include <CGAL/number_utils.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
inline bool
do_intersect(const typename K::Line_3& line,
             const typename K::Iso_cuboid_3& ic,
             const K&)
{
  typedef typename K::Point_3 Point_3;
  typedef typename K::Vector_3 Vector_3;

  const Point_3& point = line.point();
  const Vector_3& v = line.to_vector();

  return bbox_line_do_intersect_aux(point.x(), point.y(), point.z(),
                                    v.x(), v.y(), v.z(),
                                    (ic.min)().x(), (ic.min)().y(), (ic.min)().z(),
                                    (ic.max)().x(), (ic.max)().y(), (ic.max)().z());
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
