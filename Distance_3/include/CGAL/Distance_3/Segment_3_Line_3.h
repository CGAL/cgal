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

#ifndef CGAL_DISTANCE_3_SEGMENT_3_LINE_3_H
#define CGAL_DISTANCE_3_SEGMENT_3_LINE_3_H

#include <CGAL/Distance_3/internal/squared_distance_utils_3.h>

#include <CGAL/Line_3.h>
#include <CGAL/Segment_3.h>

namespace CGAL {
namespace internal {

template <class K>
typename K::FT
squared_distance(const typename K::Segment_3& seg,
                 const typename K::Line_3& line,
                 const K& k)
{
  typedef typename K::Vector_3 Vector_3;
  typedef typename K::Point_3 Point_3;
  typedef typename K::RT RT;

  typename K::Construct_vector_3 vector = k.construct_vector_3_object();

  const Point_3& linepoint = line.point();
  const Point_3& start = seg.source();
  const Point_3& end = seg.target();

  if (start == end)
    return squared_distance(start, line, k);

  Vector_3 linedir = line.direction().vector();
  Vector_3 segdir = seg.direction().vector();
  Vector_3 normal = wcross(segdir, linedir, k);

  if (is_null(normal, k))
    return squared_distance_to_line(linedir, vector(linepoint,start), k);

  bool crossing;
  RT sdm_ss2l, sdm_se2l;
  Vector_3 perpend2line, start_min_lp, end_min_lp;
  perpend2line = wcross(linedir, normal, k);
  start_min_lp = vector(linepoint, start);
  end_min_lp = vector(linepoint, end);
  sdm_ss2l = wdot(perpend2line, start_min_lp, k);
  sdm_se2l = wdot(perpend2line, end_min_lp, k);

  if (sdm_ss2l < RT(0)) {
    crossing = (sdm_se2l >= RT(0));
  } else {
    if (sdm_se2l <= RT(0)) {
      crossing = true;
    } else {
      crossing = (sdm_ss2l == RT(0));
    }
  }

  if (crossing) {
    return squared_distance_to_plane(normal, start_min_lp, k);
  } else {
    RT dm;
    dm = distance_measure_sub(sdm_ss2l, sdm_se2l, start_min_lp, end_min_lp, k);
    if (dm <= RT(0)) {
      return squared_distance_to_line(linedir, start_min_lp, k);
    } else {
      return squared_distance_to_line(linedir, end_min_lp, k);
    }
  }
}

template <class K>
typename K::FT
squared_distance(const typename K::Line_3& line,
                 const typename K::Segment_3& seg,
                 const K& k)
{
  return squared_distance(seg, line, k);
}

} // namespace internal

template <class K>
inline
typename K::FT
squared_distance(const Segment_3<K>& seg,
                 const Line_3<K>& line)
{
  return internal::squared_distance(seg, line, K());
}

template <class K>
inline
typename K::FT
squared_distance(const Line_3<K>& line,
                 const Segment_3<K>& seg)
{
  return internal::squared_distance(seg, line, K());
}

} // namespace CGAL

#endif // CGAL_DISTANCE_3_SEGMENT_3_LINE_3_H
