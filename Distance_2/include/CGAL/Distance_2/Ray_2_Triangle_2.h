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

#ifndef CGAL_DISTANCE_2_RAY_2_TRIANGLE_2_H
#define CGAL_DISTANCE_2_RAY_2_TRIANGLE_2_H

#include <CGAL/Distance_2/internal/squared_distance_utils_2.h>
#include <CGAL/Distance_2/Point_2_Triangle_2.h>

#include <CGAL/Ray_2.h>
#include <CGAL/Triangle_2.h>

namespace CGAL {
namespace internal {

template <class K>
typename K::FT
squared_distance(const typename K::Ray_2& ray,
                 const typename K::Triangle_2& triangle,
                 const K& k)
{
  typedef typename K::FT       FT;
  typedef typename K::Point_2  Point_2;
  typedef typename K::Line_2   Line_2;

  int ind_tr1, ind_tr2, ind_ray = 0, ind1;

  distance_index<K>(ind_tr1, ind_tr2, ray.source(), triangle, k);
  FT mindist = squared_distance_indexed(ray.source(), triangle, ind_tr1, ind_tr2, k);

  for(int i=0; i<3; ++i)
  {
    const Point_2& pt = triangle.vertex(i);
    distance_index<K>(ind1, pt, ray, k);
    FT dist = squared_distance_indexed(pt, ray, ind1, k);
    if(dist < mindist)
    {
      ind_ray = ind1;
      ind_tr1 = i; ind_tr2 = -1;
      mindist = dist;
    }
  }

  // now check if all vertices are on the right side of the separating line.
  // In case of vertex-vertex smallest distance this is the case.
  if(ind_tr2 == -1 && ind_ray != -1)
    return mindist;

  if(ind_tr2 != -1)
  {
    // Check if all the segment vertices lie at the same side of
    // the triangle segment.
    const Point_2& vt1 = triangle.vertex(ind_tr1);
    const Point_2& vt2 = triangle.vertex(ind_tr2);
    if(clockwise(ray.direction().vector(), vt2-vt1, k))
      mindist = FT(0);
  }
  else
  {
    // Check if all the triangle vertices lie
    // at the same side of the segment.
    const Line_2& sl = ray.supporting_line();
    Oriented_side or_s = sl.oriented_side(triangle.vertex(0));
    for(int i=1; i<3; ++i)
    {
      if(sl.oriented_side(triangle.vertex(i)) != or_s)
      {
        mindist = FT(0);
        break;
      }
    }
  }

  return mindist;
}

template <class K>
inline typename K::FT
squared_distance(const typename K::Triangle_2& triangle,
                 const typename K::Ray_2& ray,
                 const K& k)
{
  return internal::squared_distance(ray, triangle, k);
}

} // namespace internal

template <class K>
inline typename K::FT
squared_distance(const Ray_2<K>& ray,
                 const Triangle_2<K>& triangle)
{
  return K().compute_squared_distance_2_object()(ray, triangle);
}

template <class K>
inline typename K::FT
squared_distance(const Triangle_2<K>& triangle,
                 const Ray_2<K>& ray)
{
  return K().compute_squared_distance_2_object()(triangle, ray);
}

} // namespace CGAL

#endif // CGAL_DISTANCE_2_RAY_2_TRIANGLE_2_H
