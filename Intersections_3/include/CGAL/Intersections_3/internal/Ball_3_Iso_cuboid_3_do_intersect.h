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

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_BALL_3_ISO_CUBOID_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_3_BALL_3_ISO_CUBOID_3_DO_INTERSECT_H

#include <CGAL/Intersections_3/internal/Iso_cuboid_3_Sphere_3_do_intersect.h>

namespace CGAL {
namespace Intersections {
namespace internal {

// The ball and the iso-cuboid intersect iff the cuboid is not entirely
// on the unbounded side of the ball.
template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Ball_3& ball,
             const typename K::Iso_cuboid_3& ic,
             const K& k)
{
  return do_intersect_sphere_box_3(ball.center(), ball.squared_radius(),
                                   (ic.min)().x(), (ic.min)().y(), (ic.min)().z(),
                                   (ic.max)().x(), (ic.max)().y(), (ic.max)().z(),
                                   k,
                                   false);
}

template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Iso_cuboid_3& ic,
             const typename K::Ball_3& ball,
             const K& k)
{
  return do_intersect(ball, ic, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_BALL_3_ISO_CUBOID_3_DO_INTERSECT_H
