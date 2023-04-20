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
// Author(s)     : Geert-Jan Giezeman, Andreas Fabri

#ifndef CGAL_DISTANCE_3_LINE_3_PLANE_3_H
#define CGAL_DISTANCE_3_LINE_3_PLANE_3_H

#include <CGAL/Distance_3/internal/squared_distance_utils_3.h>

#include <CGAL/Line_3.h>
#include <CGAL/Plane_3.h>

namespace CGAL {
namespace internal {

template <class K>
bool
contains_vector(const typename K::Plane_3& pl,
                const typename K::Vector_3& vec,
                const K&)
{
  typedef typename K::RT RT;

  return pl.a()*vec.hx() + pl.b()*vec.hy() + pl.c() * vec.hz() == RT(0);
}

template <class K>
typename K::FT
squared_distance(const typename K::Line_3& l,
                 const typename K::Plane_3& pl,
                 const K& k)
{
  typedef typename K::FT FT;

  typename K::Compute_squared_distance_3 sq_dist = k.compute_squared_distance_3_object();

  if(contains_vector(pl, l.direction().vector(), k))
    return sq_dist(pl, l.point());

  return FT(0);
}

template <class K>
inline typename K::FT
squared_distance(const typename K::Plane_3& pl,
                 const typename K::Line_3& l,
                 const K& k)
{
  return squared_distance(l, pl, k);
}

} // namespace internal

template <class K>
inline
typename K::FT
squared_distance(const Line_3<K>& line,
                 const Plane_3<K>& plane)
{
  return K().compute_squared_distance_3_object()(line, plane);
}

template <class K>
inline
typename K::FT
squared_distance(const Plane_3<K>& plane,
                 const Line_3<K>& line)
{
  return K().compute_squared_distance_3_object()(plane, line);
}

} // namespace CGAL

#endif // CGAL_DISTANCE_3_LINE_3_PLANE_3_H
