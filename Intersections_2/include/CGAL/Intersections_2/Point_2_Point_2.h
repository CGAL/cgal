// Copyright (c) 2000
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
// Author(s)     : Geert-Jan Giezeman


#ifndef CGAL_INTERSECTIONS_2_POINT_2_POINT_2_H
#define CGAL_INTERSECTIONS_2_POINT_2_POINT_2_H

#include <CGAL/Point_2.h>
#include <CGAL/Intersection_traits_2.h>

namespace CGAL {

namespace Intersections {

namespace internal {

template <class K>
inline bool
do_intersect(const typename K::Point_2 &pt1,
             const typename K::Point_2 &pt2,
             const K&)
{
    return pt1 == pt2;
}

template <class K>
typename CGAL::Intersection_traits
<K, typename K::Point_2, typename K::Point_2>::result_type
intersection(const typename K::Point_2 &pt1,
             const typename K::Point_2 &pt2,
             const K&)
{
    if (pt1 == pt2) {
      return intersection_return<typename K::Intersect_2, typename K::Point_2, typename K::Point_2>(pt1);
    }
    return intersection_return<typename K::Intersect_2, typename K::Point_2, typename K::Point_2>();
}

}// namespace internal
} // namespace Intersections

CGAL_INTERSECTION_FUNCTION_SELF(Point_2, 2)
CGAL_DO_INTERSECT_FUNCTION_SELF(Point_2, 2)


} //namespace CGAL

#endif
