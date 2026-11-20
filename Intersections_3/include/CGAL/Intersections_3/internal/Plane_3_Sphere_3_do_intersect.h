// Copyright (c) 2003  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Pedro Machado Manhaes de Castro

#ifndef CGAL_INTERNAL_INTERSECTIONS_PLANE_3_SPHERE_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_PLANE_3_SPHERE_3_DO_INTERSECT_H

#include <CGAL/number_utils.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
inline
typename K::Boolean
do_intersect_impl(const typename K::Plane_3& p,
                  const typename K::Point_3& c,
                  const typename K::FT& r2,
                  const K&)
{
  typedef typename K::FT FT;

  const FT d2 = CGAL::square(p.a()*c.x() +
                             p.b()*c.y() +
                             p.c()*c.z() + p.d());

  return (d2 <= r2 * (square(p.a()) + square(p.b()) + square(p.c())));
}

template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Sphere_3& s,
             const typename K::Plane_3& p,
             const K& k)
{
  return do_intersect_impl(p,s.center(), s.squared_radius(), k);
}

template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Plane_3& p,
             const typename K::Sphere_3& s,
             const K& k)
{
  return do_intersect(s, p, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_PLANE_3_SPHERE_3_DO_INTERSECT_H
