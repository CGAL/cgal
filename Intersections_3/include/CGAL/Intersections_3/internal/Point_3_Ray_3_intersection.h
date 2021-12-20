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

#ifndef CGAL_INTERNAL_INTERSECTIONS_POINT_3_RAY_3_INTERSECTION_H
#define CGAL_INTERNAL_INTERSECTIONS_POINT_3_RAY_3_INTERSECTION_H

#include <CGAL/Intersection_traits_3.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
typename CGAL::Intersection_traits<K, typename K::Point_3, typename K::Ray_3>::result_type
intersection(const typename K::Point_3& pt,
             const typename K::Ray_3& ray,
             const K& k)
{
  if(do_intersect(pt, ray, k))
    return intersection_return<typename K::Intersect_3, typename K::Point_3, typename K::Ray_3>(pt);

  return intersection_return<typename K::Intersect_3, typename K::Point_3, typename K::Ray_3>();
}

template <class K>
typename CGAL::Intersection_traits<K, typename K::Ray_3, typename K::Point_3>::result_type
intersection(const typename K::Ray_3& ray,
             const typename K::Point_3& pt,
             const K& k)
{
  return internal::intersection(pt, ray, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_POINT_3_RAY_3_INTERSECTION_H
