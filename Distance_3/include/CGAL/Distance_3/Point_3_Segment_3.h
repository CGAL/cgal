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
void
squared_distance_RT(const typename K::Point_3& pt,
                    const typename K::Segment_3& seg,
                    typename K::RT& num,
                    typename K::RT& den,
                    const K& k)
{
  typedef typename K::RT RT;
  typedef typename K::Vector_3 Vector_3;

  typename K::Construct_vector_3 vector = k.construct_vector_3_object();

  // assert that the segment is valid (non zero length).
  const Vector_3 diff_s = vector(seg.source(), pt);
  const Vector_3 segvec = vector(seg.source(), seg.target());

  const RT d = wdot(diff_s, segvec, k);
  if(d <= RT(0))
  {
    // this is squared_distance(pt, seg.source())
    num = wdot(diff_s, diff_s, k);
    den = wmult((K*)0, RT(1), diff_s.hw(), diff_s.hw());
    return;
  }

  const RT e = wdot(segvec, segvec, k);
  if(wmult((K*)0, d, segvec.hw()) > wmult((K*)0, e, diff_s.hw()))
  {
    // this is squared_distance(pt, seg.target())
    const Vector_3 diff_t = vector(seg.target(), pt);
    num = wdot(diff_t, diff_t, k);
    den = wmult((K*)0, RT(1), diff_t.hw(), diff_t.hw());
    return;
  }

  // This is an expanded call to squared_distance_to_line_RT() to avoid recomputing 'e'
  const Vector_3 wcr = wcross(segvec, diff_s, k);
  num = wdot(wcr, wcr, k);
  den = wmult((K*)0, e, diff_s.hw(), diff_s.hw());
}

template <class K>
typename K::FT
squared_distance(const typename K::Point_3& pt,
                 const typename K::Segment_3& seg,
                 const K& k)
{
  typedef typename K::RT RT;
  typedef typename K::FT FT;
  typedef typename K::Vector_3 Vector_3;

  typename K::Construct_vector_3 vector = k.construct_vector_3_object();
  typename K::Compute_squared_distance_3 sq_dist = k.compute_squared_distance_3_object();

  // assert that the segment is valid (non zero length).
  const Vector_3 diff = vector(seg.source(), pt);
  const Vector_3 segvec = vector(seg.source(), seg.target());

  const RT d = wdot(diff, segvec, k);
  if(d <= RT(0))
    return (FT(diff*diff));

  const RT e = wdot(segvec, segvec, k);
  if(wmult((K*)0, d, segvec.hw()) > wmult((K*)0, e, diff.hw()))
    return sq_dist(pt, seg.target());

  // This is an expanded call to squared_distance_to_line() to avoid recomputing 'e'
  const Vector_3 wcr = wcross(segvec, diff, k);

  return FT(wcr*wcr) / wmult((K*)0, e, diff.hw(), diff.hw());
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

template <class K>
typename K::Comparison_result
compare_squared_distance(const typename K::Point_3& pt,
                         const typename K::Segment_3& seg,
                         const K& k,
                         const typename K::FT& d2)
{
  typedef typename K::RT RT;
  typedef typename K::FT FT;
  typedef typename K::Vector_3 Vector_3;

  typename K::Construct_vector_3 vector = k.construct_vector_3_object();
  typename K::Compare_squared_distance_3 csq_dist = k.compare_squared_distance_3_object();

  // assert that the segment is valid (non zero length).
  const Vector_3 diff = vector(seg.source(), pt);
  const Vector_3 segvec = vector(seg.source(), seg.target());

  //Compare first the distance to the line, if larger we can exit early
  const typename K::Comparison_result res_pl= compare_squared_distance_to_line(segvec, diff, k, d2);
  if(res_pl==LARGER)
    return LARGER;

  //If the distance is realized by the source
  const RT d = wdot(diff, segvec, k);
  if(d <= RT(0))
    return compare(FT(diff*diff), d2);

  //If the distance is realized by the target
  const RT e = wdot(segvec, segvec, k);
  if(wmult((K*)0, d, segvec.hw()) > wmult((K*)0, e, diff.hw()))
    return csq_dist(pt, seg.target(), d2);

  //If the distance is realized by the interior
  return res_pl;
}

template <class K>
inline
typename K::Comparison_result
compare_squared_distance(const typename K::Segment_3& seg,
                         const typename K::Point_3& pt,
                         const K& k,
                         const typename K::FT &d2)
{
  return compare_squared_distance(pt, seg, k, d2);
}

} // namespace internal

} // namespace CGAL

#endif // CGAL_DISTANCE_3_POINT_3_SEGMENT_3_H
