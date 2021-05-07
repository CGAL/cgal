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

#ifndef CGAL_DISTANCE_2_RAY_2_LINE_2_H
#define CGAL_DISTANCE_2_RAY_2_LINE_2_H

#include <CGAL/Distance_2/internal/squared_distance_utils_2.h>

#include <CGAL/Line_2.h>
#include <CGAL/Ray_2.h>

namespace CGAL {
namespace internal {

template <class K>
typename K::FT
squared_distance(const typename K::Line_2& line,
                 const typename K::Ray_2& ray,
                 const K& k)
{
  typedef typename K::FT FT;
  typedef typename K::Vector_2 Vector_2;

  typename K::Construct_vector_2 construct_vector = k.construct_vector_2_object();

  Vector_2 normalvec(line.a(), line.b());
  Vector_2 diff = construct_vector(line.point(), ray.source());

  FT sign_dist = k.compute_scalar_product_2_object()(diff, normalvec);
  if(sign_dist < FT(0))
  {
    if(is_acute_angle(normalvec, ray.direction().vector(), k))
      return FT(0);
  }
  else
  {
    if(is_obtuse_angle(normalvec, ray.direction().vector(), k))
      return FT(0);
  }

  return (square(sign_dist) / k.compute_squared_length_2_object()(normalvec));
}

template <class K>
inline typename K::FT
squared_distance(const typename K::Ray_2& ray,
                 const typename K::Line_2& line,
                 const K& k)
{
  return internal::squared_distance(line, ray, k);
}

} // namespace internal

template <class K>
inline typename K::FT
squared_distance(const Line_2<K>& line,
                 const Ray_2<K>& ray)
{
  return K().compute_squared_distance_2_object()(line, ray);
}

template <class K>
inline typename K::FT
squared_distance(const Ray_2<K>& ray,
                 const Line_2<K>& line)
{
  return K().compute_squared_distance_2_object()(ray, line);
}

} // namespace CGAL

#endif // CGAL_DISTANCE_2_RAY_2_LINE_2_H
