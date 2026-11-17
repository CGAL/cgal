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

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_BALL_3_BBOX_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_3_BALL_3_BBOX_3_DO_INTERSECT_H

#include <CGAL/Bbox_3.h>
#include <CGAL/Intersections_3/internal/Iso_cuboid_3_Sphere_3_do_intersect.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
typename K::Boolean
do_intersect(const typename K::Ball_3& ball,
             const CGAL::Bbox_3& bbox,
             const K& k)
{
  return do_intersect_sphere_box_3(ball.center(), ball.squared_radius(),
                                   bbox.xmin(), bbox.ymin(), bbox.zmin(),
                                   bbox.xmax(), bbox.ymax(), bbox.zmax(),
                                   k, false);
}

template <class K>
typename K::Boolean
do_intersect(const CGAL::Bbox_3& bbox,
             const typename K::Ball_3& ball,
             const K& k)
{
  return do_intersect(ball, bbox, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_BALL_3_BBOX_3_DO_INTERSECT_H
