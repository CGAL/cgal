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
// Author(s)     : Sébastien Loriot

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_LINE_3_RAY_3_INTERSECTION_H
#define CGAL_INTERNAL_INTERSECTIONS_3_LINE_3_RAY_3_INTERSECTION_H

#include <CGAL/Intersection_traits_3.h>
#include <CGAL/Intersections_3/internal/Line_3_Line_3_intersection.h>
#include <CGAL/Intersections_3/internal/Point_3_Ray_3_do_intersect.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
typename Intersection_traits<K, typename K::Line_3, typename K::Ray_3>::result_type
intersection(const typename K::Line_3& l,
             const typename K::Ray_3& r,
             const K& k)
{
  CGAL_precondition(!l.is_degenerate() && !r.is_degenerate());

  typename Intersection_traits<K, typename K::Line_3, typename K::Line_3>::result_type
      v = internal::intersection(l, r.supporting_line(), k);

  if(v)
  {
    if(const typename K::Point_3* p = intersect_get<typename K::Point_3>(v))
    {
      if(Ray_3_has_on_collinear_Point_3(r, *p, k))
        return intersection_return<typename K::Intersect_3, typename K::Line_3, typename K::Ray_3>(*p);
    }
    else if(intersect_get<typename K::Line_3>(v))
    {
      return intersection_return<typename K::Intersect_3, typename K::Line_3, typename K::Ray_3>(r);
    }
  }

  return intersection_return<typename K::Intersect_3, typename K::Line_3, typename K::Ray_3>();
}

template <class K>
typename Intersection_traits<K, typename K::Ray_3, typename K::Line_3>::result_type
intersection(const typename K::Ray_3& r,
             const typename K::Line_3& l,
             const K& k)
{
  return intersection(l, r, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_LINE_3_RAY_3_INTERSECTION_H
