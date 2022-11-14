// Copyright (c) 2018 INRIA Sophia-Antipolis (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Maxime Gimeno

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_POINT_3_INTERSECTION_H
#define CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_POINT_3_INTERSECTION_H

#include <CGAL/Intersection_traits_3.h>
#include <CGAL/Intersections_3/internal/Iso_cuboid_3_Point_3_do_intersect.h>

#include <CGAL/kernel_assertions.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
typename CGAL::Intersection_traits<K, typename K::Point_3, typename K::Iso_cuboid_3>::result_type
intersection(const typename K::Point_3& pt,
             const typename K::Iso_cuboid_3 &iso,
             const K& k)
{
  if(internal::do_intersect(pt,iso,k))
    return intersection_return<typename K::Intersect_3, typename K::Point_3, typename K::Iso_cuboid_3>(pt);

  return intersection_return<typename K::Intersect_3, typename K::Point_3, typename K::Iso_cuboid_3>();
}

template <class K>
typename CGAL::Intersection_traits
<K, typename K::Point_3, typename K::Iso_cuboid_3>::result_type
intersection(const typename K::Iso_cuboid_3 &iso,
             const typename K::Point_3& pt,
             const K& k)
{
  return internal::intersection(pt, iso, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_ISO_CUBOID_3_POINT_3_INTERSECTION_H
