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
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_DISTANCE_3_SEGMENT_3_SEGMENT_3_H
#define CGAL_DISTANCE_3_SEGMENT_3_SEGMENT_3_H

#include <CGAL/internal/squared_distance_utils_3.h>
#include <CGAL/squared_distance_Point_3_Point_3.h>

#include <CGAL/Segment_3.h>

#include <boost/algorithm/clamp.hpp>

namespace CGAL {
namespace internal {

template <class K>
typename K::FT
squared_distance_parallel(const typename K::Segment_3& seg1,
                          const typename K::Segment_3& seg2,
                          const K& k)
{
  typedef typename K::Vector_3 Vector_3;

  const Vector_3& dir1 = seg1.direction().vector();
  const Vector_3& dir2 = seg2.direction().vector();

  if(same_direction(dir1, dir2, k))
  {
    if(!is_acute_angle(seg1.source(), seg1.target(), seg2.source(), k))
      return squared_distance(seg1.target(), seg2.source(), k);
    if(!is_acute_angle(seg1.target(), seg1.source(), seg2.target(), k))
      return squared_distance(seg1.source(), seg2.target(), k);
  }
  else
  {
    if(!is_acute_angle(seg1.source(), seg1.target(), seg2.target(), k))
      return squared_distance(seg1.target(), seg2.target(), k);
    if(!is_acute_angle(seg1.target(), seg1.source(), seg2.source(), k))
      return squared_distance(seg1.source(), seg2.source(), k);
  }
  return squared_distance(seg2.source(), seg1.supporting_line(), k);
}

template <class K>
typename K::FT
squared_distance(const typename K::Segment_3& seg1,
                 const typename K::Segment_3& seg2,
                 const K& k)
{
  typedef typename K::Vector_3 Vector_3;
  typedef typename K::Point_3 Point_3;
  typedef typename K::RT RT;
  typedef typename K::FT FT;

  typename K::Construct_vector_3 construct_vector;

  const Point_3& start1 = seg1.source();
  const Point_3& start2 = seg2.source();
  const Point_3& end1 = seg1.target();
  const Point_3& end2 = seg2.target();

  if(start1 == end1)
    return squared_distance(start1, seg2, k);
  if(start2 == end2)
    return squared_distance(start2, seg1, k);

  Vector_3 dir1, dir2, normal;
  dir1 = seg1.direction().vector();
  dir2 = seg2.direction().vector();
  normal = wcross(dir1, dir2, k);
  if(is_null(normal, k))
    return squared_distance_parallel(seg1, seg2, k);

  bool crossing1, crossing2;
  RT sdm_s1to2, sdm_e1to2, sdm_s2to1, sdm_e2to1;
  Vector_3 perpend1, perpend2, s2mins1, e2mins1, e1mins2;
  perpend1 = wcross(dir1, normal, k);
  perpend2 = wcross(dir2, normal, k);
  s2mins1 = construct_vector(start1, start2);
  e2mins1 = construct_vector(start1, end2);
  e1mins2 = construct_vector(start2, end1);
  sdm_s1to2 = -RT(wdot(perpend2, s2mins1, k));
  sdm_e1to2 = wdot(perpend2, e1mins2, k);
  sdm_s2to1 = wdot(perpend1, s2mins1, k);
  sdm_e2to1 = wdot(perpend1, e2mins1, k);

  if(sdm_s1to2 < RT(0)) {
    crossing1 = (sdm_e1to2 >= RT(0));
  } else {
    if(sdm_e1to2 <= RT(0)) {
      crossing1 = true;
    } else {
      crossing1 = (sdm_s1to2 == RT(0));
    }
  }
  if(sdm_s2to1 < RT(0)) {
    crossing2 = (sdm_e2to1 >= RT(0));
  } else {
    if(sdm_e2to1 <= RT(0)) {
      crossing2 = true;
    } else {
      crossing2 = (sdm_s2to1 == RT(0));
    }
  }

  if(crossing1) {
    if(crossing2) {
      return squared_distance_to_plane(normal, s2mins1, k);
    }

    RT dm;
    dm = distance_measure_sub(sdm_s2to1, sdm_e2to1, s2mins1, e2mins1, k);
    if(dm < RT(0)) {
      return squared_distance(start2, seg1, k);
    } else {
      if(dm > RT(0)) {
        return squared_distance(end2, seg1, k);
      } else {
        // should not happen with exact arithmetic.
        return squared_distance_parallel(seg1, seg2, k);
      }
    }
  } else {
    if(crossing2) {
      RT dm;
      dm =distance_measure_sub(sdm_s1to2, sdm_e1to2, s2mins1, e1mins2, k);
      if(dm < RT(0)) {
        return squared_distance(start1, seg2, k);
      } else {
        if(dm > RT(0)) {
          return squared_distance(end1, seg2, k);
        } else {
          // should not happen with exact arithmetic.
          return squared_distance_parallel(seg1, seg2, k);
        }
      }
    } else {
      FT min1, min2;
      RT dm;
      dm = distance_measure_sub(sdm_s1to2, sdm_e1to2, s2mins1, e1mins2, k);
      if(dm == RT(0)) // should not happen with exact arithmetic.
        return squared_distance_parallel(seg1, seg2, k);
      min1 = (dm < RT(0)) ?
               squared_distance(seg1.source(), seg2, k):
               squared_distance(end1, seg2, k);
      dm = distance_measure_sub(sdm_s2to1, sdm_e2to1, s2mins1, e2mins1, k);
      if(dm == RT(0)) // should not happen with exact arithmetic.
        return squared_distance_parallel(seg1, seg2, k);
      min2 = (dm < RT(0)) ?
               squared_distance(start2, seg1, k):
               squared_distance(end2, seg1, k);
      return (min1 < min2) ? min1 : min2;
    }
  }
}

} // namespace internal

template <class K>
inline
typename K::FT
squared_distance_parallel(const Segment_3<K>& seg1,
                          const Segment_3<K>& seg2)
{
  return internal::squared_distance_parallel(seg1, seg2, K());
}

template <class K>
inline
typename K::FT
squared_distance(const Segment_3<K>& seg1,
                 const Segment_3<K>& seg2)
{
  return internal::squared_distance(seg1, seg2, K());
}

} // namespace CGAL

#endif // CGAL_DISTANCE_3_SEGMENT_3_SEGMENT_3_H
