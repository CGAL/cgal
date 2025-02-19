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

template <class K>
typename K::Comparison_result
compare_squared_distance(const typename K::Point_3& pt,
                         const typename K::Ray_3& ray,
                         const K& k,
                         const typename K::FT &d2)
{
  typedef typename K::Vector_3 Vector_3;

  typename K::Construct_vector_3 vector = k.construct_vector_3_object();

  const Vector_3 dir = ray.direction().vector();
  const Vector_3 diff = vector(ray.source(), pt);

  //Compare first the distance to the line, if larger we can exit early
  const typename K::Comparison_result res_pl = compare_squared_distance_to_line(dir, diff, k, d2);
  if(res_pl==LARGER)
    return LARGER;

  if(!is_acute_angle(dir, diff, k))
    return compare(diff*diff, d2);

  return res_pl;
}

template <class K>
typename K::Comparison_result
compare_squared_distance(const typename K::Ray_3& ray,
                         const typename K::Point_3& pt,
                         const K& k,
                         const typename K::FT &d2)
{
  return compare_squared_distance(pt, ray, k, d2);
}

} // namespace internal

} // namespace CGAL

#endif // CGAL_DISTANCE_3_POINT_3_RAY_3_H
