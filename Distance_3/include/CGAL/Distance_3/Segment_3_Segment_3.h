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

#include <CGAL/Distance_3/internal/squared_distance_utils_3.h>

#include <CGAL/Segment_3.h>

#include <boost/algorithm/clamp.hpp>

namespace CGAL {
namespace Distance_3 {
namespace internal {

template <typename K>
struct Segment_3_Segment_3_Result
{
  typename K::FT x, y;
  typename K::FT squared_distance;
};

// Using Lumelsky, 'On Fast Computation of Distance Between Line Segments' 1984
template <typename K>
Segment_3_Segment_3_Result<K>
squared_distance(const typename K::Segment_3& s1,
                 const typename K::Segment_3& s2,
                 const K& k)
{
  typedef typename K::FT                                                  FT;
  typedef typename K::Point_3                                             Point_3;
  typedef typename K::Vector_3                                            Vector_3;

  typename K::Construct_vertex_3 vertex = k.construct_vertex_3_object();
  typename K::Construct_vector_3 cv = k.construct_vector_3_object();
  typename K::Compute_scalar_product_3 sp = k.compute_scalar_product_3_object();
  typename K::Compute_squared_distance_3 sq_dist = k.compute_squared_distance_3_object();

  Segment_3_Segment_3_Result<K> res;

  const Point_3& p1 = vertex(s1, 0);
  const Point_3& q1 = vertex(s1, 1);
  const Point_3& p2 = vertex(s2, 0);
  const Point_3& q2 = vertex(s2, 1);
  const Vector_3 v1 = cv(p1, q1), v2 = cv(p2, q2);
  const Vector_3 p1p2 = cv(p1, p2);

  // @todo compute these only when needed?
  const FT a =   sp(v1, v1);
  const FT b = - sp(v1, v2);
  const FT c = - b;
  const FT d = - sp(v2, v2);
  const FT e =   sp(v1, p1p2);
  const FT f =   sp(v2, p1p2);

  if(p1 == q1)
  {
    if(p2 == q2)
    {
      res.x = 0;
      res.y = 0;
      res.squared_distance = sq_dist(p1, p2);
      return res;
    }

    CGAL_assertion(d < 0);

    res.x = 0;
    res.y = boost::algorithm::clamp<FT>(f/d, 0, 1); // (f - x*c) / d
    res.squared_distance = sq_dist(p1, p2 + res.y*v2);

    return res;
  }
  else if(p2 == q2)
  {
    CGAL_assertion(a > 0);

    res.y = 0;
    res.x = boost::algorithm::clamp<FT>(e/a, 0, 1); // (e + y*c) / a
    res.squared_distance = sq_dist(p1 + res.x*v1, p2);

    return res;
  }

  CGAL_assertion(a > 0 && d < 0);

  const FT det = a*d - b*c;
  if(det == 0)
    res.x = 0;
  else
    res.x = boost::algorithm::clamp<FT>((e*d - b*f) / det, 0, 1);

  FT xc = res.x*c;
  // res.y = (f - xc) / d, by definition, but building it up more efficiently
  if(f > xc) // y < 0 <=> f - xc / d < 0 <=> f - xc > 0 (since d is -||v2||)
  {
    res.y = 0;
    res.x = boost::algorithm::clamp<FT>(e/a, 0, 1); // (e + y*c) / a
  }
  else // y >= 0
  {
    FT n = f - xc; // delay the division by d
    if(n < d) // y > 1 <=> n / d > 1 <=> n < d (once again, important to note that d is negative!)
    {
      res.y = 1;
      res.x = boost::algorithm::clamp<FT>((e + c) / a, 0, 1); // (e + y*c) / a
    }
    else // 0 <= y <= 1
    {
      res.y = n / d;
    }
  }

  CGAL_postcondition(FT(0) <= res.x && res.x <= FT(1));
  CGAL_postcondition(FT(0) <= res.y && res.y <= FT(1));

  res.squared_distance = sq_dist(p1 + res.x*v1, p2 + res.y*v2);

  CGAL_postcondition(res.squared_distance >= FT(0));

  return res;
}

} // namespace internal
} // namespace Distance_3

namespace internal {

template <class K>
typename K::FT
squared_distance_parallel(const typename K::Segment_3& seg1,
                          const typename K::Segment_3& seg2,
                          const K& k)
{
  typedef typename K::Vector_3 Vector_3;

  typename K::Compute_squared_distance_3 sq_dist = k.compute_squared_distance_3_object();

  const Vector_3 dir1 = seg1.direction().vector();
  const Vector_3 dir2 = seg2.direction().vector();

  if(same_direction(dir1, dir2, k))
  {
    if(!is_acute_angle(seg1.source(), seg1.target(), seg2.source(), k))
      return sq_dist(seg1.target(), seg2.source(), k);
    if(!is_acute_angle(seg1.target(), seg1.source(), seg2.target(), k))
      return sq_dist(seg1.source(), seg2.target(), k);
  }
  else
  {
    if(!is_acute_angle(seg1.source(), seg1.target(), seg2.target(), k))
      return sq_dist(seg1.target(), seg2.target(), k);
    if(!is_acute_angle(seg1.target(), seg1.source(), seg2.source(), k))
      return sq_dist(seg1.source(), seg2.source(), k);
  }

  return sq_dist(seg2.source(), seg1.supporting_line(), k);
}

template <typename K>
typename K::FT
squared_distance(const typename K::Segment_3& seg1,
                 const typename K::Segment_3& seg2,
                 const K& k)
{
  Distance_3::internal::Segment_3_Segment_3_Result<K> res =
      Distance_3::internal::squared_distance(seg1, seg2, k);

  return res.squared_distance;
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
  return K().compute_squared_distance_3_object()(seg1, seg2);
}

} // namespace CGAL

#endif // CGAL_DISTANCE_3_SEGMENT_3_SEGMENT_3_H
