// Copyright (c) 2001 GeometryFactory(France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Pedro Machado Manhaes

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_SPHERE_3_SPHERE_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_3_SPHERE_3_SPHERE_3_DO_INTERSECT_H

#include <CGAL/Intersections_3/internal/Plane_3_Sphere_3_do_intersect.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
inline
bool
do_intersect(const typename K::Sphere_3& s1,
             const typename K::Sphere_3& s2,
             const K& k)
{
  typedef typename K::Plane_3 Plane_3;

  typename K::Construct_center_3 center = k.construct_center_3_object();
  typename K::Compute_squared_radius_3 sqr = k.compute_squared_radius_3_object();

  //  Weaker than s1 == s2 since it allows opposite orientations
  if(k.equal_3_object()(center(s1), center(s2)))
    return (sqr(s1) == sqr(s2));

  Plane_3 pl = k.construct_radical_plane_3_object()(s1, s2);

  return do_intersect(pl, s1, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_SPHERE_3_SPHERE_3_DO_INTERSECT_H
