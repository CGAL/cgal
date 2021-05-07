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

#ifndef CGAL_DISTANCE_2_TRIANGLE_2_TRIANGLE_2_H
#define CGAL_DISTANCE_2_TRIANGLE_2_TRIANGLE_2_H

#include <CGAL/Distance_2/internal/squared_distance_utils_2.h>
#include <CGAL/Distance_2/Point_2_Point_2.h>
#include <CGAL/Distance_2/Point_2_Triangle_2.h>

#include <CGAL/Triangle_2.h>

namespace CGAL {
namespace internal {

template <class K>
typename K::FT
squared_distance(const typename K::Triangle_2& triangle1,
                 const typename K::Triangle_2& triangle2,
                 const K& k)
{
  typedef typename K::FT       FT;
  typedef typename K::Point_2  Point_2;

  typename K::Orientation_2 orientation = k.orientation_2_object();

  int ind1_1 = 0, ind1_2 = -1, ind2_1 = 0, ind2_2 = -1, ind1, ind2;
  FT dist;

  FT mindist = internal::squared_distance(triangle1.vertex(0), triangle2.vertex(0), k);
  for(int i=0; i<3; ++i)
  {
    const Point_2& pt = triangle1.vertex(i);
    distance_index<K>(ind1, ind2, pt, triangle2, k);
    dist = squared_distance_indexed(pt, triangle2, ind1, ind2, k);
    if(dist < mindist)
    {
      ind1_1 = i; ind1_2 = -1;
      ind2_1 = ind1; ind2_2 = ind2;
      mindist = dist;
    }
  }

  for(int i=0; i<3; ++i)
  {
    const Point_2& pt = triangle2.vertex(i);
    distance_index<K>(ind1, ind2, pt, triangle1, k);
    dist = squared_distance_indexed(pt, triangle1, ind1, ind2, k);
    if(dist < mindist)
    {
      ind1_1 = ind1; ind1_2 = ind2;
      ind2_1 = i; ind2_2 = -1;
      mindist = dist;
    }
  }

  // now check if all vertices are on the right side of the separating line.
  if(ind1_2 == -1 && ind2_2 == -1)
    return mindist;

  // In case of point-segment closest distance, there is still the
  // possibility of overlapping triangles.  Check if all the
  // vertices lie at the same side of the segment.
  if(ind1_2 != -1)
  {
    const Point_2& vt1 = triangle1.vertex(ind1_1);
    const Point_2& vt2 = triangle1.vertex(ind1_2);
    const Orientation or_s = orientation(vt1, vt2, triangle2.vertex(0));
    for(int i=1; i<3; ++i)
    {
      if(orientation(vt1, vt2, triangle2.vertex(i)) != or_s)
      {
        mindist = FT(0);
        break;
      }
    }
  }
  else
  {
    const Point_2& vt1 = triangle2.vertex(ind2_1);
    const Point_2& vt2 = triangle2.vertex(ind2_2);
    const Orientation or_s = orientation(vt1, vt2, triangle1.vertex(0));
    for(int i=1; i<3; ++i)
    {
      if(orientation(vt1, vt2, triangle1.vertex(i)) != or_s)
      {
        mindist = FT(0);
        break;
      }
    }
  }

  return mindist;
}

} // namespace internal

template <class K>
inline typename K::FT
squared_distance(const Triangle_2<K>& triangle1,
                 const Triangle_2<K>& triangle2)
{
  return K().compute_squared_distance_2_object()(triangle1, triangle2);
}

} // namespace CGAL

#endif // CGAL_DISTANCE_2_TRIANGLE_2_TRIANGLE_2_H
