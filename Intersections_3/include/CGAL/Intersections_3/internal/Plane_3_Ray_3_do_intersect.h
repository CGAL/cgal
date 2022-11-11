// Copyright (c) 1997-2010
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
// Author(s)     : Geert-Jan Giezeman <geert@cs.uu.nl>

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_PLANE_3_RAY_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_3_PLANE_3_RAY_3_DO_INTERSECT_H

#include <CGAL/enum.h>
#include <CGAL/number_utils.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
bool
do_intersect(const typename K::Plane_3& plane,
             const typename K::Ray_3& ray,
             const K& k)
{
  typename K::Oriented_side_3 oriented_side_3 = k.oriented_side_3_object();
  CGAL::Oriented_side os = oriented_side_3(plane, ray.source());
  if(os == ON_ORIENTED_BOUNDARY)
    return true;

  return (sign(ray.to_vector()* plane.orthogonal_vector()) * os == -1);
}

template <class K>
inline
bool
do_intersect(const typename K::Ray_3& ray,
             const typename K::Plane_3& plane,
             const K& k)
{
  return do_intersect(plane, ray, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_PLANE_3_RAY_3_DO_INTERSECT_H
