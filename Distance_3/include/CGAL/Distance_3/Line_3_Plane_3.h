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

template <class K>
inline typename K::Comparison_result
compare_squared_distance(const typename K::Line_3& l,
                         const typename K::Plane_3& pl,
                         const K& k,
                         const typename K::FT& d2)
{
  return compare(squared_distance(l, pl, k), d2);
}

template <class K>
inline typename K::Comparison_result
compare_squared_distance(const typename K::Plane_3& pl,
                         const typename K::Line_3& l,
                         const K& k,
                         const typename K::FT& d2)
{
  return compare_squared_distance(l, pl, k, d2);
}

} // namespace internal

} // namespace CGAL

#endif // CGAL_DISTANCE_3_LINE_3_PLANE_3_H
