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

  const Vector_3 dir1 = seg.direction().vector();
  const Vector_3 dir2 = ray.direction().vector();

  bool same_direction;
  if(CGAL_NTS abs(dir1.hx()) > CGAL_NTS abs(dir1.hy()))
    same_direction = (CGAL_NTS sign(dir1.hx()) == CGAL_NTS sign(dir2.hx()));
  else
    same_direction = (CGAL_NTS sign(dir1.hy()) == CGAL_NTS sign(dir2.hy()));

  if(same_direction)
  {
    if(!is_acute_angle(seg.source(), seg.target(), ray.source(), k))
      return squared_distance(seg.target(), ray.source(), k);
  }
  else
  {
    if(!is_acute_angle(seg.target(), seg.source(), ray.source(), k))
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
  typedef typename K::RT RT;
  typedef typename K::FT FT;
  typedef typename K::Point_3 Point_3;
  typedef typename K::Vector_3 Vector_3;

  typename K::Construct_vector_3 vector = k.construct_vector_3_object();
  typename K::Compute_squared_distance_3 sq_dist = k.compute_squared_distance_3_object();

  const Point_3& ss = seg.source();
  const Point_3& se = seg.target();

  if(ss == se)
    return sq_dist(ss, ray);

  const Vector_3 raydir = ray.direction().vector();
  const Vector_3 segdir = seg.direction().vector();
  const Vector_3 normal = wcross(segdir, raydir, k);

  if(is_null(normal, k))
    return squared_distance_parallel(seg, ray, k);

  bool crossing1, crossing2;

  const Vector_3 perpend2seg = wcross(segdir, normal, k);
  const Vector_3 perpend2ray = wcross(raydir, normal, k);
  const Vector_3 ss_min_rs = vector(ray.source(), ss);
  const Vector_3 se_min_rs = vector(ray.source(), se);
  const RT sdm_ss2r = wdot(perpend2ray, ss_min_rs, k);
  const RT sdm_se2r = wdot(perpend2ray, se_min_rs, k);

  if(sdm_ss2r < RT(0))
  {
    crossing1 = (sdm_se2r >= RT(0));
  }
  else
  {
    if(sdm_se2r <= RT(0))
      crossing1 = true;
    else
      crossing1 = (sdm_ss2r == RT(0));
  }

  const RT sdm_rs2s = - wdot(perpend2seg, ss_min_rs, k);
  const RT sdm_re2s = wdot(perpend2seg, raydir, k);
  if(sdm_rs2s < RT(0))
  {
    crossing2 = (sdm_re2s >= RT(0));
  } else
  {
    if(sdm_re2s <= RT(0))
      crossing2 = true;
    else
      crossing2 = (sdm_rs2s == RT(0));
  }

  if(crossing1)
  {
    if(crossing2)
      return squared_distance_to_plane(normal, ss_min_rs, k);

    return sq_dist(ray.source(), seg);
  }
  else
  {
    if(crossing2)
    {
      const RT dm = distance_measure_sub(sdm_ss2r, sdm_se2r, ss_min_rs, se_min_rs, k);
      if(dm < RT(0))
      {
        return sq_dist(ss, ray);
      }
      else
      {
        if(dm > RT(0))
          return sq_dist(se, ray);
        else
          // parallel, should not happen (no crossing)
          return squared_distance_parallel(seg, ray, k);
      }
    }
    else
    {
      const RT dm = distance_measure_sub(sdm_ss2r, sdm_se2r, ss_min_rs, se_min_rs, k);
      if(dm == RT(0))
        return squared_distance_parallel(seg, ray, k);

      const FT min1 = (dm < RT(0)) ? sq_dist(ss, ray)
                                   : sq_dist(se, ray);
      const FT min2 = sq_dist(ray.source(), seg);

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
  return internal::squared_distance_parallel(ray, seg, K());
}

template <class K>
inline
typename K::FT
squared_distance(const Segment_3<K>& seg,
                 const Ray_3<K>& ray)
{
  return K().compute_squared_distance_3_object()(seg, ray);
}

template <class K>
inline
typename K::FT
squared_distance(const Ray_3<K>& ray,
                 const Segment_3<K>& seg)
{
  return K().compute_squared_distance_3_object()(ray, seg);
}

} // namespace CGAL

#endif // CGAL_DISTANCE_3_SEGMENT_3_RAY_3_H
