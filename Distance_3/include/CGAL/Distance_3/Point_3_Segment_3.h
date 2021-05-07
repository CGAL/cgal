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

#ifndef CGAL_DISTANCE_3_POINT_3_SEGMENT_3_H
#define CGAL_DISTANCE_3_POINT_3_SEGMENT_3_H

#include <CGAL/Distance_3/internal/squared_distance_utils_3.h>
#include <CGAL/Distance_3/Point_3_Point_3.h>

#include <CGAL/Point_3.h>
#include <CGAL/Segment_3.h>

namespace CGAL {
namespace internal {

template <class K>
typename K::FT
squared_distance(const typename K::Point_3& pt,
                 const typename K::Segment_3& seg,
                 const K& k,
                 const Homogeneous_tag&)
{
  typedef typename K::Vector_3 Vector_3;
  typedef typename K::RT RT;
  typedef typename K::FT FT;

  typename K::Construct_vector_3 vector = k.construct_vector_3_object();
  typename K::Compute_squared_distance_3 sq_dist = k.compute_squared_distance_3_object();

  // assert that the segment is valid (non zero length).
  Vector_3 diff = vector(seg.source(), pt);
  Vector_3 segvec = vector(seg.source(), seg.target());

  RT d = wdot(diff,segvec, k);
  if(d <= RT(0))
    return diff*diff;
  RT e = wdot(segvec,segvec, k);
  if((d * segvec.hw()) > (e * diff.hw()))
    return sq_dist(pt, seg.target());

  Vector_3 wcr = wcross(segvec, diff, k);
  return FT(wcr*wcr) / FT(e * diff.hw() * diff.hw());
}

template <class K>
typename K::FT
squared_distance(const typename K::Point_3& pt,
                 const typename K::Segment_3& seg,
                 const K& k,
                 const Cartesian_tag&)
{
  typedef typename K::Vector_3 Vector_3;
  typedef typename K::RT RT;
  typedef typename K::FT FT;

  typename K::Construct_vector_3 vector = k.construct_vector_3_object();
  typename K::Compute_squared_distance_3 sqd = k.compute_squared_distance_3_object();

  // assert that the segment is valid (non zero length).
  Vector_3 diff = vector(seg.source(), pt);
  Vector_3 segvec = vector(seg.source(), seg.target());

  RT d = wdot(diff,segvec, k);
  if(d <= RT(0))
    return diff*diff;

  RT e = wdot(segvec,segvec, k);
  if(d > e)
    return sqd(pt, seg.target());

  Vector_3 wcr = wcross(segvec, diff, k);
  return FT(wcr*wcr)/e;
}

template <class K>
inline
typename K::FT
squared_distance(const typename K::Point_3& pt,
                 const typename K::Segment_3& seg,
                 const K& k)
{
  typedef typename K::Kernel_tag Tag;
  Tag tag;
  return squared_distance(pt, seg, k, tag);
}

template <class K>
inline
typename K::FT
squared_distance(const typename K::Segment_3& seg,
                 const typename K::Point_3& pt,
                 const K& k)
{
  return squared_distance(pt, seg, k);
}

} // namespace internal

template <class K>
inline
typename K::FT
squared_distance(const Point_3<K>& pt,
                 const Segment_3<K>& seg)
{
  return K().compute_squared_distance_3_object()(pt, seg);
}

template <class K>
inline
typename K::FT
squared_distance(const Segment_3<K>& seg,
                 const Point_3<K>& pt)
{
  return K().compute_squared_distance_3_object()(seg, pt);
}

} // namespace CGAL

#endif // CGAL_DISTANCE_3_POINT_3_SEGMENT_3_H
