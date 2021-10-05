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
// Author(s)     : Andreas Fabri

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_POINT_3_SPHERE_3_INTERSECTION_H
#define CGAL_INTERNAL_INTERSECTIONS_3_POINT_3_SPHERE_3_INTERSECTION_H

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
inline
bool
do_intersect(const typename K::Point_3& pt,
             const typename K::Sphere_3& sphere,
             const K& k)
{
  return k.has_on_boundary_3_object()(sphere, pt);
}

template <class K>
inline
bool
do_intersect(const typename K::Sphere_3& sphere,
             const typename K::Point_3& pt,
             const K& k)
{
  return k.has_on_boundary_3_object()(sphere, pt);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_POINT_3_SPHERE_3_INTERSECTION_H
