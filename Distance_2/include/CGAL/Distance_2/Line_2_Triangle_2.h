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

#ifndef CGAL_DISTANCE_2_LINE_2_TRIANGLE_2_H
#define CGAL_DISTANCE_2_LINE_2_TRIANGLE_2_H

#include <CGAL/Distance_2/internal/squared_distance_utils_2.h>
#include <CGAL/Distance_2/Point_2_Line_2.h>

#include <CGAL/enum.h>

#include <CGAL/Line_2.h>
#include <CGAL/Triangle_2.h>

namespace CGAL {
namespace internal {

template <class K>
typename K::FT
squared_distance(const typename K::Line_2& line,
                 const typename K::Triangle_2& triangle,
                 const K& k)
{
  typedef typename K::FT FT;

  Oriented_side side0 = line.oriented_side(triangle.vertex(0));
  if(line.oriented_side(triangle.vertex(1)) != side0)
    return FT(0);

  if(line.oriented_side(triangle.vertex(2)) != side0)
    return FT(0);

  FT mindist = internal::squared_distance(triangle.vertex(0), line, k);
  for(int i=1; i<3; ++i)
  {
    FT dist = internal::squared_distance(triangle.vertex(i), line, k);
    if(dist < mindist)
      mindist = dist;
  }

  return mindist;
}

template <class K>
inline typename K::FT
squared_distance(const typename K::Triangle_2& triangle,
                 const typename K::Line_2& line,
                 const K& k)
{
  return internal::squared_distance(line, triangle, k);
}

} // namespace internal

template <class K>
inline typename K::FT
squared_distance(const Line_2<K>& line,
                 const Triangle_2<K>& triangle)
{
  return internal::squared_distance(line, triangle, K());
}

template <class K>
inline typename K::FT
squared_distance(const Triangle_2<K>& triangle,
                 const Line_2<K>& line)
{
  return internal::squared_distance(line, triangle, K());
}

} // namespace CGAL

#endif // CGAL_DISTANCE_2_LINE_2_TRIANGLE_2_H
