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
// Author(s)     : Sébastien Loriot

#ifndef CGAL_INTERNAL_INTERSECTIONS_POINT_3_RAY_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_POINT_3_RAY_3_DO_INTERSECT_H

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
inline
bool
do_intersect(const typename K::Point_3& pt,
             const typename K::Ray_3& ray,
             const K& k)
{
  return k.has_on_3_object()(ray, pt);
}

template <class K>
inline
bool
do_intersect(const typename K::Ray_3& ray,
             const typename K::Point_3& pt,
             const K& k)
{
  return k.has_on_3_object()(ray, pt);
}


template <class K>
bool
Ray_3_has_on_collinear_Point_3(const typename K::Ray_3& r,
                               const typename K::Point_3& p,
                               const K& k)
{
  return k.equal_3_object()(r.source(), p) ||
           k.equal_3_object()(
             k.construct_direction_3_object()(
               k.construct_vector_3_object()(r.source(), p)), r.direction());
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_POINT_3_RAY_3_DO_INTERSECT_H
