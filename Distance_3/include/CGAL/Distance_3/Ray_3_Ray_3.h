// Copyright (c) 1998-2021
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Geert-Jan Giezeman, Andreas Fabri

#ifndef CGAL_DISTANCE_3_RAY_3_RAY_3_H
#define CGAL_DISTANCE_3_RAY_3_RAY_3_H

#include <CGAL/Distance_3/internal/squared_distance_utils_3.h>

#include <CGAL/Plane_3.h>
#include <CGAL/Vector_3.h>

namespace CGAL {
namespace internal {

template <class K>
typename K::FT
ray_ray_squared_distance_parallel(const typename K::Vector_3& ray1dir,
                                  const typename K::Vector_3& ray2dir,
                                  const typename K::Vector_3& s1_min_s2,
                                  const K& k)
{
  if (!is_acute_angle(ray2dir, s1_min_s2, k))
    if (!same_direction(ray1dir, ray2dir, k))
      return typename K::FT(s1_min_s2*s1_min_s2);

  return squared_distance_to_line(ray1dir, s1_min_s2, k);
}

template <class K>
typename K::FT
squared_distance(const typename K::Ray_3& ray1,
                 const typename K::Ray_3& ray2,
                 const K& k)
{
  typedef typename K::Vector_3 Vector_3;
  typedef typename K::Point_3 Point_3;
  typedef typename K::RT RT;
  typedef typename K::FT FT;

  typename K::Construct_vector_3 construct_vector;

  const Point_3& s1 = ray1.source();
  const Point_3& s2 = ray2.source();
  Vector_3 dir1, dir2, normal;
  dir1 = ray1.direction().vector();
  dir2 = ray2.direction().vector();
  normal = wcross(dir1, dir2, k);
  Vector_3 s1_min_s2 = construct_vector(s2, s1);
  if (is_null(normal, k))
    return ray_ray_squared_distance_parallel(dir1, dir2, s1_min_s2, k);

  bool crossing1, crossing2;
  RT sdm_s1_2, sdm_s2_1;
  Vector_3 perpend1, perpend2;
  perpend1 = wcross(dir1, normal, k);
  perpend2 = wcross(dir2, normal, k);

  sdm_s1_2 = wdot(perpend2, s1_min_s2, k);
  if (sdm_s1_2 < RT(0)) {
    crossing1 = (RT(wdot(perpend2, dir1, k)) >= RT(0));
  } else {
    if (RT(wdot(perpend2, dir1, k)) <= RT(0)) {
      crossing1 = true;
    } else {
      crossing1 = (sdm_s1_2 == RT(0));
    }
  }

  sdm_s2_1 = -RT(wdot(perpend1, s1_min_s2, k));
  if (sdm_s2_1 < RT(0)) {
    crossing2 = (RT(wdot(perpend1, dir2, k)) >= RT(0));
  } else {
    if (RT(wdot(perpend1, dir2, k)) <= RT(0)) {
      crossing2 = true;
    } else {
      crossing2 = (sdm_s2_1 == RT(0));
    }
  }

  if (crossing1) {
    if (crossing2)
      return squared_distance_to_plane(normal, s1_min_s2, k);
    return squared_distance(ray2.source(), ray1, k);
  } else {
    if (crossing2) {
      return squared_distance(ray1.source(), ray2, k);
    } else {
      FT min1, min2;
      min1 = squared_distance(ray1.source(), ray2, k);
      min2 = squared_distance(ray2.source(), ray1, k);
      return (min1 < min2) ? min1 : min2;
    }
  }
}

} // namespace internal

template <class K>
inline
typename K::FT
ray_ray_squared_distance_parallel(const Vector_3<K>& ray1dir,
                                  const Vector_3<K>& ray2dir,
                                  const Vector_3<K>& s1_min_s2)
{
  return internal::ray_ray_squared_distance_parallel(ray1dir, ray2dir, s1_min_s2, K());
}

template <class K>
inline
typename K::FT
squared_distance(const Ray_3<K>& ray1,
                 const Ray_3<K>& ray2)
{
  return internal::squared_distance(ray1, ray2, K());
}

} // namespace CGAL

#endif // CGAL_DISTANCE_3_RAY_3_RAY_3_H
