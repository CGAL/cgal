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
// Author(s)     : Andreas Fabri

#ifndef CGAL_DISTANCE_3_POINT_3_TETRAHEDRON_3_H
#define CGAL_DISTANCE_3_POINT_3_TETRAHEDRON_3_H

#include <CGAL/Distance_3/Point_3_Triangle_3.h>

#include <CGAL/Point_3.h>
#include <CGAL/Tetrahedron_3.h>

namespace CGAL {
namespace internal {

template <class K>
inline
typename K::FT
squared_distance(const typename K::Point_3& pt,
                 const typename K::Tetrahedron_3& tet,
                 const K& k)
{
  typedef typename K::Point_3 Point_3;

  typename K::Construct_vertex_3 vertex = k.construct_vertex_3_object();
  typename K::Orientation_3 orientation = k.orientation_3_object();

  bool on_bounded_side = true;
  const Point_3& t0 = vertex(tet, 0);
  const Point_3& t1 = vertex(tet, 1);
  const Point_3& t2 = vertex(tet, 2);
  const Point_3& t3 = vertex(tet, 3);

  bool dmin_initialized = false;
  typename K::FT dmin;
  bool inside = false;
  if(orientation(t0,t1,t2, pt) == NEGATIVE)
  {
    on_bounded_side = false;
    dmin = squared_distance_to_triangle(pt, t0, t1, t2, k, inside);
    dmin_initialized = true;
    if(inside)
      return dmin;
  }

  if(orientation(t0,t3,t1, pt) == NEGATIVE)
  {
    on_bounded_side = false;
    const typename K::FT d = squared_distance_to_triangle(pt, t0, t3, t1, k, inside);
    if(inside)
      return d;

    if(!dmin_initialized)
    {
      dmin = d;
      dmin_initialized = true;
    }
    else
    {
      dmin = (std::min)(d, dmin);
    }
  }

  if(orientation(t1,t3,t2, pt) == NEGATIVE)
  {
    on_bounded_side = false;
    const typename K::FT d = squared_distance_to_triangle(pt, t1, t3, t2, k, inside);
    if(inside)
      return d;

    if(!dmin_initialized)
    {
      dmin = d;
      dmin_initialized = true;
    }
    else
    {
      dmin = (std::min)(d, dmin);
    }
  }

  if(orientation(t2,t3,t0, pt) == NEGATIVE)
  {
    on_bounded_side = false;
    const typename K::FT d = squared_distance_to_triangle(pt, t2, t3, t0, k, inside);
    if(inside)
      return d;

    if(!dmin_initialized)
    {
      dmin = d;
      dmin_initialized = true;
    }
    else
    {
      dmin = (std::min)(d, dmin);
    }
  }

  if(on_bounded_side)
    return typename K::FT(0);

  return dmin;
}

template <class K>
inline
typename K::FT
squared_distance(const typename K::Tetrahedron_3& tet,
                 const typename K::Point_3& pt,
                 const K& k)
{
  return squared_distance(pt, tet, k);
}

} // namespace internal

template <class K>
typename K::FT
squared_distance(const Tetrahedron_3<K>& tet,
                 const Point_3<K>& pt)
{
  return K().compute_squared_distance_3_object()(tet, pt);
}

template <class K>
typename K::FT
squared_distance(const Point_3<K>& pt,
                 const Tetrahedron_3<K>& tet)
{
  return K().compute_squared_distance_3_object()(pt, tet);
}

} // namespace CGAL

#endif // CGAL_DISTANCE_3_POINT_3_TETRAHEDRON_3_H
