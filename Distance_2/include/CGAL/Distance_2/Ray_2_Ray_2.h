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

#ifndef CGAL_DISTANCE_2_RAY_2_RAY_2_H
#define CGAL_DISTANCE_2_RAY_2_RAY_2_H

#include <CGAL/Distance_2/internal/squared_distance_utils_2.h>
#include <CGAL/Distance_2/Point_2_Ray_2.h>

#include <CGAL/Ray_2.h>

namespace CGAL {
namespace internal {

template <class K>
typename K::FT
ray_ray_squared_distance_parallel(const typename K::Vector_2& ray1dir,
                                  const typename K::Vector_2& ray2dir,
                                  const typename K::Vector_2& from1to2,
                                  const K& k)
{
  typedef typename K::RT RT;
  typedef typename K::FT FT;

  if(!is_acute_angle(ray1dir, from1to2, k))
  {
    if(!same_direction(ray1dir, ray2dir, k))
      return k.compute_squared_length_2_object()(from1to2);
  }

  RT wcr = wcross(ray1dir, from1to2, k);
  RT w = from1to2.hw();

  return (square(wcr) / FT(wmult((K*)0, wdot(ray1dir, ray1dir, k), w, w)));
}

template <class K>
typename K::FT
squared_distance(const typename K::Ray_2& ray1,
                 const typename K::Ray_2& ray2,
                 const K& k)
{
  typedef typename K::Vector_2 Vector_2;
  typedef typename K::FT FT;

  typename K::Construct_vector_2 vector = k.construct_vector_2_object();
  typename K::Compute_squared_distance_2 sq_dist = k.compute_squared_distance_2_object();

  const Vector_2 ray1dir = ray1.direction().vector();
  const Vector_2 ray2dir = ray2.direction().vector();
  const Vector_2 diffvec = vector(ray1.source(),ray2.source());

  bool crossing1, crossing2;
  switch(orientation(ray1dir, ray2dir, k))
  {
    case COUNTERCLOCKWISE:
      crossing1 = !clockwise(diffvec, ray2dir, k);
      crossing2 = !counterclockwise(ray1dir, diffvec, k);
      break;
    case CLOCKWISE:
      crossing1 = !counterclockwise(diffvec, ray2dir, k);
      crossing2 = !clockwise(ray1dir, diffvec, k);
      break;
    default:
      return ray_ray_squared_distance_parallel(ray1dir, ray2dir, diffvec,k);
  }

  if(crossing1)
  {
    if(crossing2)
      return FT(0);
    return sq_dist(ray2.source(), ray1);
  }
  else
  {
    if(crossing2)
    {
      return sq_dist(ray1.source(), ray2);
    }
    else
    {
      FT min1 = sq_dist(ray1.source(), ray2);
      FT min2 = sq_dist(ray2.source(), ray1);
      return (min1 < min2) ? min1 : min2;
    }
  }
}

} // namespace internal

template <class K>
inline typename K::FT
squared_distance(const Ray_2<K> &ray1, const Ray_2<K> &ray2)
{
  return K().compute_squared_distance_2_object()(ray1, ray2);
}

} // namespace CGAL

#endif // CGAL_DISTANCE_2_RAY_2_RAY_2_H
