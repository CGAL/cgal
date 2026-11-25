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

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_BALL_3_TRIANGLE_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_3_BALL_3_TRIANGLE_3_DO_INTERSECT_H

#include <CGAL/Rational_traits.h>
#include <CGAL/Distance_3/Point_3_Triangle_3.h>

namespace CGAL {
namespace Intersections {
namespace internal {

// The ball and the triangle intersect iff the distance between the center
// of the ball and the triangle is not larger than the radius.
template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Ball_3& b,
             const typename K::Triangle_3& tr,
             const K& k)
{
  typedef typename K::RT RT;
  RT num, den;

  CGAL::internal::squared_distance_RT(b.center(), tr, num, den, k);
  return !(compare_quotients<RT>(num, den,
                                 Rational_traits<typename K::FT>().numerator(b.squared_radius()),
                                 Rational_traits<typename K::FT>().denominator(b.squared_radius())) == LARGER);
}

template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Triangle_3& tr,
             const typename K::Ball_3& b,
             const K& k)
{
  return do_intersect(b, tr, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_BALL_3_TRIANGLE_3_DO_INTERSECT_H
