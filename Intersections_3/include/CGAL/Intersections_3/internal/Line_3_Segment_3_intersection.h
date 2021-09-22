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

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_LINE_3_SEGMENT_3_INTERSECTION_H
#define CGAL_INTERNAL_INTERSECTIONS_3_LINE_3_SEGMENT_3_INTERSECTION_H

#include <CGAL/Intersection_traits_3.h>
#include <CGAL/Intersections_3/internal/Line_3_Line_3_intersection.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
typename Intersection_traits<K, typename K::Line_3, typename K::Segment_3>::result_type
intersection(const typename K::Line_3& l,
             const typename K::Segment_3& s,
             const K& k)
{
  CGAL_precondition(!l.is_degenerate() && !s.is_degenerate());

  typename Intersection_traits<K, typename K::Line_3, typename K::Line_3>::result_type
      v = internal::intersection(l, s.supporting_line(), k);

  if(v)
  {
    if(const typename K::Point_3* p = intersect_get<typename K::Point_3> (v))
    {
      typename K::Collinear_are_ordered_along_line_3 cln_order = k.collinear_are_ordered_along_line_3_object();
      if(cln_order(s[0], *p, s[1]))
        return intersection_return<typename K::Intersect_3, typename K::Line_3, typename K::Segment_3>(*p);
    } else {
      return intersection_return<typename K::Intersect_3, typename K::Line_3, typename K::Segment_3>(s);
    }
  }

  return intersection_return<typename K::Intersect_3, typename K::Line_3, typename K::Segment_3>();
}

template <class K>
typename Intersection_traits<K, typename K::Line_3, typename K::Segment_3>::result_type
intersection(const typename K::Segment_3& s,
             const typename K::Line_3& l,
             const K& k)
{
  return intersection(l, s, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_LINE_3_SEGMENT_3_INTERSECTION_H
