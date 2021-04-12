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

#ifndef CGAL_DISTANCE_3_RAY_3_LINE_3_H
#define CGAL_DISTANCE_3_RAY_3_LINE_3_H

#include <CGAL/Distance_3/internal/squared_distance_utils_3.h>

#include <CGAL/Line_3.h>
#include <CGAL/Ray_3.h>

namespace CGAL {
namespace internal {

template <class K>
typename K::FT
squared_distance(const typename K::Line_3& line,
                 const typename K::Ray_3& ray,
                 const K& k)
{
  typedef typename K::Vector_3 Vector_3;
  typedef typename K::Point_3 Point_3;
  typedef typename K::RT RT;

  typename K::Construct_vector_3 construct_vector;

  const Point_3& rs = ray.source();

  Vector_3 raydir, linedir, normal;
  linedir = line.direction().vector();
  raydir = ray.direction().vector();
  normal = wcross(raydir, linedir, k);

  Vector_3 rs_min_lp = construct_vector(line.point(), rs);
  if (is_null(normal, k))
    return squared_distance_to_line(linedir, rs_min_lp, k);

  bool crossing;
  RT sdm_sr_l;
  Vector_3 perpend2l;
  perpend2l = wcross(linedir, normal, k);

  sdm_sr_l = wdot(perpend2l, rs_min_lp, k);
  if (sdm_sr_l < RT(0))
  {
    crossing = (RT(wdot(perpend2l, raydir, k)) >= RT(0));
  }
  else
  {
    if (RT(wdot(perpend2l, raydir, k)) <= RT(0))
      crossing = true;
    else
      crossing = (sdm_sr_l == RT(0));
  }

  if (crossing)
    return squared_distance_to_plane(normal, rs_min_lp, k);
  else
    return squared_distance_to_line(linedir, rs_min_lp, k);
}

} // namespace internal

template <class K>
inline
typename K::FT
squared_distance(const Line_3<K> &line,
                 const Ray_3<K> &ray)
{
  return internal::squared_distance(line, ray, K());
}

template <class K>
inline
typename K::FT
squared_distance(const Ray_3<K>& ray,
                 const Line_3<K>& line)
{
  return internal::squared_distance(line, ray, K());
}

} // namespace CGAL

#endif // CGAL_DISTANCE_3_RAY_3_LINE_3_H
