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

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_LINE_3_POINT_3_INTERSECTION_H
#define CGAL_INTERNAL_INTERSECTIONS_3_LINE_3_POINT_3_INTERSECTION_H

#include <CGAL/Intersection_traits_3.h>
#include <CGAL/Intersections_3/internal/Line_3_Point_3_do_intersect.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
typename CGAL::Intersection_traits<K, typename K::Point_3, typename K::Line_3>::result_type
intersection(const typename K::Point_3& pt,
             const typename K::Line_3& line,
             const K& k)
{
  if(do_intersect(pt,line,k))
    return intersection_return<typename K::Intersect_3, typename K::Point_3, typename K::Line_3>(pt);

  return intersection_return<typename K::Intersect_3, typename K::Point_3, typename K::Line_3>();
}

template <class K>
typename CGAL::Intersection_traits<K, typename K::Line_3, typename K::Point_3>::result_type
intersection(const typename K::Line_3& line,
             const typename K::Point_3& pt,
             const K& k)
{
  return intersection(pt, line, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_LINE_3_POINT_3_INTERSECTION_H
