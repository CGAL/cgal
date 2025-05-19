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

#ifndef CGAL_DISTANCE_2_POINT_2_RAY_2_H
#define CGAL_DISTANCE_2_POINT_2_RAY_2_H

#include <CGAL/Distance_2/internal/squared_distance_utils_2.h>
#include <CGAL/Distance_2/Point_2_Line_2.h>

#include <CGAL/Point_2.h>
#include <CGAL/Ray_2.h>

namespace CGAL {
namespace internal {

template <class K>
void
distance_index(int &ind,
               const typename K::Point_2 &pt,
               const typename K::Ray_2 &ray,
               const K& k)
{
  typename K::Construct_vector_2 construct_vector = k.construct_vector_2_object();
  if(!is_acute_angle(ray.direction().vector(), construct_vector(ray.source(), pt), k))
  {
    ind = 0;
    return;
  }

  ind = -1;
}

template <class K>
inline typename K::FT
squared_distance_indexed(const typename K::Point_2 &pt,
                         const typename K::Ray_2 &ray,
                         int ind,
                         const K& k)
{
  if(ind == 0)
    return internal::squared_distance(pt, ray.source(), k);

  return internal::squared_distance(pt, ray.supporting_line(), k);
}

template <class K>
typename K::FT
squared_distance(const typename K::Point_2& pt,
                 const typename K::Ray_2& ray,
                 const K& k)
{
  typedef typename K::Vector_2 Vector_2;

  typename K::Construct_vector_2 construct_vector = k.construct_vector_2_object();

  Vector_2 diff = construct_vector(ray.source(), pt);
  const Vector_2& dir = ray.direction().vector();
  if (!is_acute_angle(dir, diff, k))
    return k.compute_squared_length_2_object()(diff);

  return internal::squared_distance(pt, ray.supporting_line(), k);
}

template <class K>
inline typename K::FT
squared_distance(const typename K::Ray_2& ray,
                 const typename K::Point_2& pt,
                 const K& k)
{
  return internal::squared_distance(pt, ray, k);
}

} // namespace internal

template <class K>
inline typename K::FT
squared_distance(const Point_2<K>& pt,
                 const Ray_2<K>& ray)
{
  return internal::squared_distance(pt, ray, K());
}

template <class K>
inline typename K::FT
squared_distance(const Ray_2<K>& ray,
                 const Point_2<K>& pt)
{
  return internal::squared_distance(pt, ray, K());
}

} // namespace CGAL

#endif // CGAL_DISTANCE_2_POINT_2_RAY_2_H
