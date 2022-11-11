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
// Author(s)     : Geert-Jan Giezeman

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_PLANE_3_SEGMENT_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_3_PLANE_3_SEGMENT_3_DO_INTERSECT_H

#include <CGAL/enum.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
bool
do_intersect(const typename K::Plane_3& plane,
             const typename K::Segment_3& seg,
             const K&)
{
  typedef typename K::Point_3 Point_3;

  const Point_3& source = seg.source();
  const Point_3& target = seg.target();

  CGAL::Oriented_side source_side = plane.oriented_side(source);
  CGAL::Oriented_side target_side = plane.oriented_side(target);

  if(source_side == target_side && target_side != ON_ORIENTED_BOUNDARY)
    return false;

  return true;
}

template <class K>
inline
bool
do_intersect(const typename K::Segment_3& seg,
             const typename K::Plane_3& plane,
             const K& k)
{
  return do_intersect(plane, seg, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_PLANE_3_SEGMENT_3_DO_INTERSECT_H
