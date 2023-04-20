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
// Author(s)     : Geert-Jan Giezeman, Andreas Fabri, Mael Rouxel-Labbé

#ifndef CGAL_DISTANCE_3_POINT_3_RAY_3_H
#define CGAL_DISTANCE_3_POINT_3_RAY_3_H

#include <CGAL/Distance_3/internal/squared_distance_utils_3.h>

#include <CGAL/Point_3.h>
#include <CGAL/Ray_3.h>

namespace CGAL {
namespace internal {

template <class K>
void
squared_distance_RT(const typename K::Point_3 &pt,
                    const typename K::Ray_3 &ray,
                    typename K::RT& num,
                    typename K::RT& den,
                    const K& k)
{
  typedef typename K::Vector_3 Vector_3;
  typedef typename K::RT RT;

  typename K::Construct_vector_3 vector = k.construct_vector_3_object();

  const Vector_3 dir = ray.direction().vector();
  const Vector_3 diff = vector(ray.source(), pt);

  if(!is_acute_angle(dir, diff, k))
  {
    num = wdot(diff, diff, k);
    den = wmult((K*)0, RT(1), diff.hw(), diff.hw());
    return;
  }

  squared_distance_to_line_RT(dir, diff, num, den, k);
}

template <class K>
typename K::FT
squared_distance(const typename K::Point_3& pt,
                 const typename K::Ray_3& ray,
                 const K& k)
{
  // This duplicates code from the _RT functions, but it is a slowdown to do something like:
  //
  //   RT num, den;
  //   squared_distance_RT(pt, ray, num, den, k);
  //   return Rational_traits<FT>().make_rational(num, den);
  //
  // See https://github.com/CGAL/cgal/pull/5680

  typedef typename K::Vector_3 Vector_3;

  typename K::Construct_vector_3 vector = k.construct_vector_3_object();

  const Vector_3 dir = ray.direction().vector();
  const Vector_3 diff = vector(ray.source(), pt);

  if(!is_acute_angle(dir, diff, k))
    return diff*diff;

  return squared_distance_to_line(dir, diff, k);
}

template <class K>
typename K::FT
squared_distance(const typename K::Ray_3& ray,
                 const typename K::Point_3& pt,
                 const K& k)
{
  return squared_distance(pt, ray, k);
}

} // namespace internal

template <class K>
inline
typename K::FT
squared_distance(const Point_3<K>& pt,
                 const Ray_3<K>& ray)
{
  return K().compute_squared_distance_3_object()(pt, ray);
}

template <class K>
inline
typename K::FT
squared_distance(const Ray_3<K>& ray,
                 const Point_3<K>& pt)
{
  return K().compute_squared_distance_3_object()(ray, pt);
}

} // namespace CGAL

#endif // CGAL_DISTANCE_3_POINT_3_RAY_3_H
