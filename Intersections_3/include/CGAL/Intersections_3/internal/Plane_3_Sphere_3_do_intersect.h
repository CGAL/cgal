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
bool
do_intersect(const typename K::Plane_3& p,
             const typename K::Sphere_3& s,
             const K&)
{
  typedef typename K::FT FT;

  const FT d2 = CGAL::square(p.a()*s.center().x() +
                             p.b()*s.center().y() +
                             p.c()*s.center().z() + p.d());

  return (d2 <= s.squared_radius() * (square(p.a()) + square(p.b()) + square(p.c())));
}

template <class K>
inline
bool
do_intersect(const typename K::Sphere_3& s,
             const typename K::Plane_3& p,
             const K&)
{
  return do_intersect(p,s);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_PLANE_3_SPHERE_3_DO_INTERSECT_H
