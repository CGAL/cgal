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

#ifndef CGAL_DISTANCE_3_SEGMENT_3_RAY_3_H
#define CGAL_DISTANCE_3_SEGMENT_3_RAY_3_H

#include <CGAL/Distance_3/internal/squared_distance_utils_3.h>
#include <CGAL/Distance_3/Point_3_Line_3.h>
#include <CGAL/Distance_3/Point_3_Point_3.h>
#include <CGAL/Distance_3/Point_3_Ray_3.h>
#include <CGAL/Distance_3/Point_3_Segment_3.h>

#include <CGAL/Segment_3.h>
#include <CGAL/Ray_3.h>

namespace CGAL {
namespace internal {

template <class K>
typename K::FT
squared_distance_parallel(const typename K::Segment_3& seg,
                          const typename K::Ray_3& ray,
                          const K& k)
{
  typedef typename K::Vector_3 Vector_3;

  bool same_direction;
  const Vector_3 &dir1 = seg.direction().vector();
  const Vector_3 &dir2 = ray.direction().vector();

  if (CGAL_NTS abs(dir1.hx()) > CGAL_NTS abs(dir1.hy())) {
    same_direction = (CGAL_NTS sign(dir1.hx()) == CGAL_NTS sign(dir2.hx()));
  } else {
    same_direction = (CGAL_NTS sign(dir1.hy()) == CGAL_NTS sign(dir2.hy()));
  }

  if (same_direction) {
    if (!is_acute_angle(seg.source(), seg.target(), ray.source(), k))
      return squared_distance(seg.target(), ray.source(), k);
  } else {
    if (!is_acute_angle(seg.target(), seg.source(), ray.source(), k))
      return squared_distance(seg.source(), ray.source(), k);
  }

  return squared_distance(ray.source(), seg.supporting_line(), k);
}

template <class K>
typename K::FT
squared_distance(const typename K::Segment_3& seg,
                 const typename K::Ray_3& ray,
                 const K& k)
{
  typedef typename K::Point_3 Point_3;
  typedef typename K::Vector_3 Vector_3;
  typedef typename K::RT RT;
  typedef typename K::FT FT;

  typename K::Construct_vector_3 construct_vector;

  const Point_3& ss = seg.source();
  const Point_3& se = seg.target();

  if (ss == se)
    return squared_distance(ss, ray, k);

  Vector_3 raydir, segdir, normal;
  raydir = ray.direction().vector();
  segdir = seg.direction().vector();
  normal = wcross(segdir, raydir, k);

  if (is_null(normal, k))
    return squared_distance_parallel(seg, ray, k);

  bool crossing1, crossing2;
  RT sdm_ss2r, sdm_se2r, sdm_rs2s, sdm_re2s;
  Vector_3 perpend2seg, perpend2ray, ss_min_rs, se_min_rs;
  perpend2seg = wcross(segdir, normal, k);
  perpend2ray = wcross(raydir, normal, k);
  ss_min_rs = construct_vector(ray.source(), ss);
  se_min_rs = construct_vector(ray.source(), se);
  sdm_ss2r = wdot(perpend2ray, ss_min_rs, k);
  sdm_se2r = wdot(perpend2ray, se_min_rs, k);

  if (sdm_ss2r < RT(0)) {
    crossing1 = (sdm_se2r >= RT(0));
  } else {
    if (sdm_se2r <= RT(0)) {
      crossing1 = true;
    } else {
      crossing1 = (sdm_ss2r == RT(0));
    }
  }

  sdm_rs2s = -RT(wdot(perpend2seg, ss_min_rs, k));
  sdm_re2s = wdot(perpend2seg, raydir, k);
  if (sdm_rs2s < RT(0)) {
    crossing2 = (sdm_re2s >= RT(0));
  } else {
    if (sdm_re2s <= RT(0)) {
      crossing2 = true;
    } else {
      crossing2 = (sdm_rs2s == RT(0));
    }
  }

  if (crossing1) {
    if (crossing2) {
      return squared_distance_to_plane(normal, ss_min_rs, k);
    }
    return squared_distance(ray.source(), seg, k);
  } else {
    if (crossing2) {
      RT dm;
      dm = distance_measure_sub(sdm_ss2r, sdm_se2r, ss_min_rs, se_min_rs, k);
      if (dm < RT(0)) {
        return squared_distance(ss, ray, k);
      } else {
        if (dm > RT(0)) {
          return squared_distance(se, ray, k);
        } else {
          // parallel, should not happen (no crossing)
          return squared_distance_parallel(seg, ray, k);
        }
      }
    } else {
      FT min1, min2;
      RT dm;
      dm = distance_measure_sub(sdm_ss2r, sdm_se2r, ss_min_rs, se_min_rs, k);
      if (dm == RT(0))
        return squared_distance_parallel(seg, ray, k);
      min1 = (dm < RT(0))
             ? squared_distance(ss, ray, k)
             : squared_distance(se, ray, k);
      min2 = squared_distance(ray.source(), seg, k);
      return (min1 < min2) ? min1 : min2;
    }
  }
}

template <class K>
typename K::FT
squared_distance(const typename K::Ray_3& ray,
                 const typename K::Segment_3& seg,
                 const K& k)
{
  return squared_distance(seg, ray, k);
}

} // namespace internal

template <class K>
inline
typename K::FT
squared_distance_parallel(const Segment_3<K>& seg,
                          const Ray_3<K>& ray)
{
  return internal::squared_distance_parallel(ray,seg, K());
}

template <class K>
inline
typename K::FT
squared_distance(const Segment_3<K>& seg,
                 const Ray_3<K>& ray)
{
  return internal::squared_distance(seg, ray, K());
}

template <class K>
inline
typename K::FT
squared_distance(const Ray_3<K>& ray,
                 const Segment_3<K>& seg)
{
  return internal::squared_distance(seg, ray, K());
}

} // namespace CGAL

#endif // CGAL_DISTANCE_3_SEGMENT_3_RAY_3_H
