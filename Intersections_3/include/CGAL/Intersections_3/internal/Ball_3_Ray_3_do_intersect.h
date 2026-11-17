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

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_BALL_3_RAY_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_3_BALL_3_RAY_3_DO_INTERSECT_H

#include <CGAL/Rational_traits.h>
#include <CGAL/Intersections_3/internal/Ray_3_Sphere_3_do_intersect.h>

namespace CGAL {
namespace Intersections {
namespace internal {


template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Ray_3& ray,
             const typename K::Ball_3& ball,
             const K& k)
{
  return do_intersect_impl(ball.center(), ball.squared_radius(), ray, k);
}

template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Ball_3& ball,
             const typename K::Ray_3& ray,
             const K& k)
{
  return do_intersect(ray, ball, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_BALL_3_RAY_3_DO_INTERSECT_H
