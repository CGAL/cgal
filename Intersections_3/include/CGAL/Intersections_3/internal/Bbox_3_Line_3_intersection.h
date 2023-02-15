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
//                 Sebastien Loriot <Sebastien.Loriot@geometryfactory.com>

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_LINE_3_INTERSECTION_H
#define CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_LINE_3_INTERSECTION_H

#include <CGAL/Intersection_traits_3.h>
#include <CGAL/Intersections_3/internal/Bbox_3_Segment_3_intersection.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/number_utils.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
typename Intersection_traits<K, typename K::Line_3, Bbox_3>::result_type
intersection(const typename K::Line_3& line,
             const Bbox_3& box,
             const K&)
{
  typedef typename K::Point_3 Point_3;
  typedef typename K::Direction_3 Direction_3;

  const Point_3& linepoint = line.point();
  const Direction_3& linedir = line.direction();

  return intersection_bl<K>(box,
                            CGAL::to_double(linepoint.x()),
                            CGAL::to_double(linepoint.y()),
                            CGAL::to_double(linepoint.z()),
                            CGAL::to_double(linedir.dx()),
                            CGAL::to_double(linedir.dy()),
                            CGAL::to_double(linedir.dz()),
                            true, true);
}

template <class K>
inline
typename Intersection_traits<K, Bbox_3, typename K::Line_3>::result_type
intersection(const Bbox_3& box,
             const typename K::Line_3& line,
             const K& k)
{
  return intersection(line, box, k);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_BBOX_3_LINE_3_INTERSECTION_H
