// Copyright (c) 1998-2004
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
// Author(s)     : Geert-Jan Giezeman
//                 Michel Hoffmann <hoffmann@inf.ethz.ch>
//                 Andreas Fabri <Andreas.Fabri@geometryfactory.com>

#ifndef CGAL_DISTANCE_2_SEGMENT_2_RAY_2_H
#define CGAL_DISTANCE_2_SEGMENT_2_RAY_2_H

#include <CGAL/Distance_2/internal/squared_distance_utils_2.h>

#include <CGAL/Distance_2/Point_2_Line_2.h>
#include <CGAL/Distance_2/Point_2_Point_2.h>
#include <CGAL/Distance_2/Point_2_Ray_2.h>

#include <CGAL/Segment_2.h>
#include <CGAL/Ray_2.h>

namespace CGAL {
namespace internal {

template <class K>
inline typename K::RT
_distance_measure_sub(const typename K::RT& startwcross,
                      const typename K::RT& endwcross,
                      const typename K::Vector_2& start,
                      const typename K::Vector_2& end)
{
  return CGAL_NTS abs(wmult((K*)0, startwcross, end.hw())) -
           CGAL_NTS abs(wmult((K*)0, endwcross, start.hw()));
}

template <class K>
typename K::FT
squared_distance_parallel(const typename K::Segment_2& seg,
                          const typename K::Ray_2& ray,
                          const K& k)
{
  typedef typename K::Vector_2 Vector_2;

  typename K::Compute_squared_distance_2 sq_dist = k.compute_squared_distance_2_object();

  const Vector_2 dir1 = seg.direction().vector();
  const Vector_2 dir2 = ray.direction().vector();

  if(same_direction(dir1, dir2, k))
  {
    if(!is_acute_angle(seg.source(), seg.target(), ray.source(), k))
      return sq_dist(seg.target(), ray.source());
  }
  else
  {
    if(!is_acute_angle(seg.target(), seg.source(), ray.source(), k))
      return sq_dist(seg.source(), ray.source());
  }

  return sq_dist(ray.source(), seg.supporting_line());
}

template <class K>
typename K::FT
squared_distance(const typename K::Segment_2& seg,
                 const typename K::Ray_2& ray,
                 const K& k)
{
  typedef typename K::RT RT;
  typedef typename K::FT FT;
  typedef typename K::Vector_2 Vector_2;

  typename K::Construct_vector_2 vector = k.construct_vector_2_object();
  typename K::Orientation_2 orientation = k.orientation_2_object();
  typename K::Compute_squared_distance_2 sq_dist = k.compute_squared_distance_2_object();

  const Vector_2& raydir = ray.direction().vector();
  const Vector_2 startvec = vector(ray.source(), seg.source());
  const Vector_2 endvec = vector(ray.source(), seg.target());

  if(seg.source() == seg.target())
    return internal::squared_distance(seg.source(), ray, k);

  bool crossing1, crossing2;
  const RT c1s = wcross(raydir, startvec, k);
  const RT c1e = wcross(raydir, endvec, k);
  if(c1s < RT(0))
  {
    crossing1 = (c1e >= RT(0));
  }
  else
  {
    if(c1e <= RT(0))
    {
      if(c1s == RT(0) && c1e == RT(0))
        return internal::squared_distance_parallel(seg, ray, k);

      crossing1 = true;
    }
    else
    {
      crossing1 = (c1s == RT(0));
    }
  }

  switch (orientation(seg.source(), seg.target(), ray.source()))
  {
    case LEFT_TURN:
      crossing2 = right_turn(vector(seg.source(), seg.target()), raydir, k);
      break;
    case RIGHT_TURN:
      crossing2 = left_turn(vector(seg.source(), seg.target()), raydir, k);
      break;
    default:
      crossing2 = true;
      break;
  }

  if(crossing1)
  {
    if(crossing2)
      return FT(0);

    return sq_dist(ray.source(), seg);
  }
  else
  {
    if(crossing2)
    {
      const RT dm = _distance_measure_sub<K>(c1s, c1e, startvec, endvec);
      if(dm < RT(0))
      {
        return sq_dist(seg.source(), ray);
      }
      else
      {
        if(dm > RT(0))
          return sq_dist(seg.target(), ray);
        else // parallel, should not happen (no crossing)
          return internal::squared_distance_parallel(seg, ray, k);
      }
    }
    else
    {
      const RT dm = _distance_measure_sub<K>(c1s, c1e, startvec, endvec);
      if(dm == RT(0))
        return internal::squared_distance_parallel(seg, ray, k);

      const FT min1 = (dm < RT(0)) ? sq_dist(seg.source(), ray)
                                   : sq_dist(seg.target(), ray);
      const FT min2 = sq_dist(ray.source(), seg);

      return (min1 < min2) ? min1 : min2;
    }
  }
}

template <class K>
inline typename K::FT
squared_distance(const typename K::Ray_2& ray,
                 const typename K::Segment_2& seg,
                 const K& k)
{
  return internal::squared_distance(seg, ray, k);
}

} // namespace internal

template <class K>
inline typename K::FT
squared_distance(const Segment_2<K>& seg,
                 const Ray_2<K>& ray)
{
  return K().compute_squared_distance_2_object()(seg, ray);
}

template <class K>
inline typename K::FT
squared_distance(const Ray_2<K>& ray,
                 const Segment_2<K>& seg)
{
  return K().compute_squared_distance_2_object()(ray, seg);
}

} // namespace CGAL

#endif // CGAL_DISTANCE_2_SEGMENT_2_RAY_2_H
