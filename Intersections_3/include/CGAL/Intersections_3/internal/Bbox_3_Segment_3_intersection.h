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
// Author(s)     : Michael.Hemmer@sophia.inria.fr

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_SEGMENT_3_INTERSECTION_H
#define CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_SEGMENT_3_INTERSECTION_H

#include <CGAL/Intersection_traits_3.h>
#include <CGAL/Intersections_3/internal/Iso_cuboid_3_Segment_3_intersection.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/Iso_cuboid_3.h>
#include <CGAL/number_utils.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
typename Intersection_traits<K, typename K::Segment_3, Bbox_3>::result_type
intersection(const typename K::Segment_3& seg,
             const Bbox_3& box,
             const K& k)
{
  // Delegate to the Iso_cuboid_3/Segment_3 intersection which uses exact
  // kernel arithmetic throughout, avoiding the to_double() precision loss
  // that the old intersection_bl() helper suffered from (issue #7124).
  typedef typename K::Iso_cuboid_3 Iso_cuboid_3;
  typedef typename Intersection_traits<K, typename K::Segment_3, Bbox_3>::result_type result_type;
  typedef typename Intersection_traits<K, typename K::Segment_3, Bbox_3>::variant_type variant_type;

  auto res = internal::intersection(seg, Iso_cuboid_3(box), k);
  if(!res) return result_type();
  return std::visit([](auto&& v) -> result_type {
    return result_type(variant_type(std::forward<decltype(v)>(v)));
  }, *res);
}

template <class K>
inline
typename Intersection_traits<K, Bbox_3, typename K::Segment_3>::result_type
intersection(const Bbox_3& box,
             const typename K::Segment_3& seg,
             const K& k)
{
  return intersection(seg, box, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_SEGMENT_3_INTERSECTION_H
