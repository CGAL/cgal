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

#ifndef CGAL_DISTANCE_2_SEGMENT_2_TRIANGLE_2_H
#define CGAL_DISTANCE_2_SEGMENT_2_TRIANGLE_2_H

#include <CGAL/Distance_2/internal/squared_distance_utils_2.h>
#include <CGAL/Distance_2/Point_2_Triangle_2.h>

#include <CGAL/Segment_2.h>
#include <CGAL/Triangle_2.h>

namespace CGAL {
namespace internal {

template <class K>
typename K::FT
squared_distance(const typename K::Segment_2& seg,
                 const typename K::Triangle_2& triangle,
                 const K& k)
{
  typedef typename K::FT       FT;
  typedef typename K::Point_2  Point_2;

  typename K::Orientation_2 orientation = k.orientation_2_object();

  int i, ind_tr1 = 0, ind_tr2 = -1, ind_seg = 0, ind1, ind2;
  FT dist;

  FT mindist = internal::squared_distance(seg.source(), triangle.vertex(0), k);
  for(i=0; i<2; ++i)
  {
    const Point_2& pt = seg.vertex(i);
    distance_index<K>(ind1, ind2, pt, triangle, k);
    dist = internal::squared_distance_indexed(pt, triangle, ind1, ind2, k);
    if(dist < mindist)
    {
      ind_seg = i;
      ind_tr1 = ind1; ind_tr2 = ind2;
      mindist = dist;
    }
  }

  for(i=0; i<3; ++i)
  {
    const Point_2& pt = triangle.vertex(i);
    distance_index<K>(ind1, pt, seg, k);
    dist = internal::squared_distance_indexed(pt, seg, ind1, k);
    if(dist < mindist)
    {
      ind_seg = ind1;
      ind_tr1 = i; ind_tr2 = -1;
      mindist = dist;
    }
  }

  // now check if all vertices are on the right side of the separating line.
  // In case of vertex-vertex smallest distance this is the case.
  if(ind_tr2 == -1 && ind_seg != -1)
    return mindist;

  if(ind_tr2 != -1)
  {
    // Check if all the segment vertices lie at the same side of
    // the triangle segment.
    const Point_2 &vt1 = triangle.vertex(ind_tr1);
    const Point_2 &vt2 = triangle.vertex(ind_tr2);
    const Orientation or_s = orientation(vt1, vt2, seg.source());
    if(orientation(vt1, vt2, seg.target()) != or_s)
      mindist = FT(0);
  }
  else
  {
    // Check if all the triangle vertices lie
    // at the same side of the segment.
    const Point_2& vt1 = seg.source();
    const Point_2& vt2 = seg.target();
    const Orientation or_s = orientation(vt1, vt2, triangle.vertex(0));
    for(i=1; i<3; ++i)
    {
      if(orientation(vt1, vt2, triangle.vertex(i)) != or_s)
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
squared_distance(const typename K::Triangle_2 & triangle,
                 const typename K::Segment_2 & seg,
                 const K& k)
{
  return internal::squared_distance(seg, triangle, k);
}

} // namespace internal

template <class K>
inline typename K::FT
squared_distance(const Segment_2<K>& seg,
                 const Triangle_2<K>& triangle)
{
  return K().compute_squared_distance_2_object()(seg, triangle);
}

template <class K>
inline typename K::FT
squared_distance(const Triangle_2<K>& triangle,
                 const Segment_2<K>& seg)
{
  return K().compute_squared_distance_2_object()(triangle, seg);
}

} // namespace CGAL

#endif // CGAL_DISTANCE_2_SEGMENT_2_TRIANGLE_2_H
