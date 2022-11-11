// Copyright (c) 2008 INRIA Sophia-Antipolis (France), ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Pedro Machado Manhaes de Castro

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_LINE_3_LINE_3_DO_INTERSECT_H
#define CGAL_INTERNAL_INTERSECTIONS_3_LINE_3_LINE_3_DO_INTERSECT_H

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
bool
do_intersect(const typename K::Line_3& l1,
             const typename K::Line_3& l2,
             const K& k)
{
  typedef typename K::Point_3      Point_3;
  typedef typename K::Vector_3     Vector_3;

  if(k.has_on_3_object()(l1, l2.point()))
    return true;
  if(k.are_parallel_3_object()(l1,l2))
    return false;

  const Point_3& p1 = l1.point();
  const Point_3& p3 = l2.point();
  const Vector_3 &v1 = l1.to_vector();
  const Vector_3& v2 = l2.to_vector();
  const Point_3 p2 = p1 + v1;
  const Point_3 p4 = p2 + v2;
  return k.coplanar_3_object()(p1,p2,p3,p4);
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_LINE_3_LINE_3_DO_INTERSECT_H
