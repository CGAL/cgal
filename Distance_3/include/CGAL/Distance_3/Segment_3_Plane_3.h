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

#ifndef CGAL_DISTANCE_3_SEGMENT_3_PLANE_3_H
#define CGAL_DISTANCE_3_SEGMENT_3_PLANE_3_H

#include <CGAL/Distance_3/internal/squared_distance_utils_3.h>

#include <CGAL/Segment_3.h>
#include <CGAL/Plane_3.h>

namespace CGAL {
namespace internal {

template <class K>
typename K::FT
squared_distance(const typename K::Segment_3 &seg,
                 const typename K::Plane_3 &plane,
                 const K& k)
{
  typedef typename K::RT RT;
  typedef typename K::FT FT;
  typedef typename K::Vector_3 Vector_3;
  typedef typename K::Point_3 Point_3;

  typename K::Construct_vector_3 vector = k.construct_vector_3_object();

  const Point_3& start = seg.start();
  const Point_3& end = seg.end();

  if (start == end)
    return squared_distance(start, plane, k);

  const Point_3& planepoint = plane.point();
  const Vector_3 start_min_pp = vector(planepoint, start);
  const Vector_3 end_min_pp = vector(planepoint, end);
  const Vector_3& normal = plane.orthogonal_vector();

  const RT sdm_ss2pp = wdot(normal, start_min_pp, k);
  const RT sdm_se2pp = wdot(normal, end_min_pp, k);

  switch (CGAL_NTS sign(sdm_ss2pp))
  {
    case -1:
      if (sdm_se2pp >= RT(0))
        return FT(0);
      if (sdm_ss2pp * end_min_pp.hw() >= sdm_se2pp * start_min_pp.hw())
        return squared_distance_to_plane(normal, start_min_pp, k);
      else
        return squared_distance_to_plane(normal, end_min_pp, k);
    case 0:
    default:
      return FT(0);
    case 1:
      if (sdm_se2pp <= RT(0))
        return FT(0);
      if (sdm_ss2pp * end_min_pp.hw() <= sdm_se2pp * start_min_pp.hw())
        return squared_distance_to_plane(normal, start_min_pp, k);
      else
        return squared_distance_to_plane(normal, end_min_pp, k);
  }
}

template <class K>
inline typename K::FT
squared_distance(const typename K::Plane_3& plane,
                 const typename K::Segment_3& seg,
                 const K& k)
{
  return squared_distance(seg, plane, k);
}

} // namespace internal

template <class K>
inline
typename K::FT
squared_distance(const Segment_3<K>& seg,
                 const Plane_3<K>& plane)
{
  return K().compute_squared_distance_3_object()(seg, plane);
}

template <class K>
inline
typename K::FT
squared_distance(const Plane_3<K>& plane,
                 const Segment_3<K>& seg)
{
  return K().compute_squared_distance_3_object()(plane, seg);
}

} // namespace CGAL

#endif // CGAL_DISTANCE_3_SEGMENT_3_PLANE_3_H
