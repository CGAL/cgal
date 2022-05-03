// Copyright (c) 2018 GeometryFactory Sarl
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

#ifndef CGAL_INTERNAL_INTERSECTIONS_SPHERE_3_TRIANGLE_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_SPHERE_3_TRIANGLE_3_DO_INTERSECT_H

#include <CGAL/Rational_traits.h>
#include <CGAL/Distance_3/Point_3_Triangle_3.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Sphere_3& sp,
             const typename K::Triangle_3& tr,
             const K& k)
{
  typedef typename K::RT RT;
  typedef typename K::Bounded_side Bounded_side;

  typename K::Bounded_side_3 bounded_side = k.bounded_side_3_object();
  typename K::Compute_squared_radius_3 sq_radius = k.compute_squared_radius_3_object();
  typename K::Construct_center_3 center = k.construct_center_3_object();
  typename K::Construct_vertex_3 vertex = k.construct_vertex_3_object();

  const Bounded_side v0_side = bounded_side(sp, vertex(tr, 0));
  const Bounded_side v1_side = bounded_side(sp, vertex(tr, 1));
  const Bounded_side v2_side = bounded_side(sp, vertex(tr, 2));

  if(v0_side != v1_side || v0_side != v2_side || v1_side != v2_side)
    return true;
  else if(v0_side == ON_BOUNDED_SIDE) // 'else if' ==> all vertices are on the same side
    return false;
  else if(v0_side == ON_BOUNDARY)
    return true;

  // All vertices are out, but the triangle can still be intersecting the sphere
  RT num, den;
  CGAL::internal::squared_distance_RT(center(sp), tr, num, den, k);

  return !(compare_quotients<RT>(num, den,
                                 Rational_traits<typename K::FT>().numerator(sq_radius(sp)),
                                 Rational_traits<typename K::FT>().denominator(sq_radius(sp))) == LARGER);
}

template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Triangle_3& tr,
             const typename K::Sphere_3& sp,
             const K& k)
{
  return do_intersect(sp, tr, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_SPHERE_3_TRIANGLE_3_DO_INTERSECT_H
