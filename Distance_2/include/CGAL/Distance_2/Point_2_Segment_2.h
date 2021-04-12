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

#ifndef CGAL_DISTANCE_2_POINT_2_SEGMENT_2_H
#define CGAL_DISTANCE_2_POINT_2_SEGMENT_2_H

#include <CGAL/Distance_2/internal/squared_distance_utils_2.h>
#include <CGAL/Distance_2/Point_2_Point_2.h>
#include <CGAL/Distance_2/Point_2_Line_2.h>

#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>

namespace CGAL {
namespace internal {

template <class K>
void
distance_index(int &ind,
               const typename K::Point_2 &pt,
               const typename K::Segment_2 &seg,
               const K& k)
{
  if(!is_acute_angle(seg.target(),seg.source(),pt, k))
  {
    ind = 0;
    return;
  }

  if(!is_acute_angle(seg.source(),seg.target(),pt, k))
  {
    ind = 1;
    return;
  }
  ind = -1;
}

template <class K>
inline typename K::FT
squared_distance_indexed(const typename K::Point_2 &pt,
                         const typename K::Segment_2 &seg,
                         int ind,
                         const K& k)
{
  if(ind == 0)
    return internal::squared_distance(pt, seg.source(), k);
  if(ind == 1)
    return internal::squared_distance(pt, seg.target(), k);
  return internal::squared_distance(pt, seg.supporting_line(), k);
}

template <class K>
typename K::FT
squared_distance(const typename K::Point_2& pt,
                 const typename K::Segment_2& seg,
                 const K& k)
{
  typedef typename K::Vector_2 Vector_2;
  typedef typename K::RT RT;

  typename K::Construct_vector_2 construct_vector = k.construct_vector_2_object();

  // assert that the segment is valid (non zero length).
  Vector_2 diff = construct_vector(seg.source(), pt);
  Vector_2 segvec = construct_vector(seg.source(), seg.target());

  RT d = wdot(diff, segvec, k);
  if(d <= RT(0))
    return k.compute_squared_length_2_object()(diff);

  RT e = wdot(segvec,segvec, k);
  if(wmult((K*)0 ,d, segvec.hw()) > wmult((K*)0, e, diff.hw()))
    return internal::squared_distance(pt, seg.target(), k);

  return internal::squared_distance(pt, seg.supporting_line(), k);
}

template <class K>
inline typename K::FT
squared_distance(const typename K::Segment_2& seg,
                 const typename K::Point_2& pt,
                 const K& k)
{
  return internal::squared_distance(pt, seg, k);
}

} // namespace internal

template <class K>
inline typename K::FT
squared_distance(const Point_2<K>& pt,
                 const Segment_2<K>& seg)
{
  return internal::squared_distance(pt, seg, K());
}

template <class K>
inline typename K::FT
squared_distance(const Segment_2<K>& seg,
                 const Point_2<K>& pt)
{
  return internal::squared_distance(pt, seg, K());
}

} // namespace CGAL

#endif // CGAL_DISTANCE_2_POINT_2_SEGMENT_2_H
