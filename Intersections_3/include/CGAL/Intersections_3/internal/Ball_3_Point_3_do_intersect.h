// Copyright (c) 2026 GeometryFactory Sarl
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sébastien Loriot
//

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_BALL_3_POINT_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_3_BALL_3_POINT_3_DO_INTERSECT_H


namespace CGAL {
namespace Intersections {
namespace internal {

// A point intersects a ball iff it is on its bounded side or on its boundary.
template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Ball_3& b,
             const typename K::Point_3& pt,
             const K& k)
{
  typename K::Comparison_result res = k.compare_squared_distance_3_object()(b.center(), pt, b.squared_radius());
  return res != LARGER;
}

template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Point_3& pt,
             const typename K::Ball_3& b,
             const K& k)
{
  return do_intersect(b, pt, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_BALL_3_POINT_3_DO_INTERSECT_H
