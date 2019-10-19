// Copyright (c) 2010 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot
//

#ifndef CGAL_INTERSECTIONS_3_BBOX_3_POINT_3_H
#define CGAL_INTERSECTIONS_3_BBOX_3_POINT_3_H

#include <CGAL/Bbox_3.h>
#include <CGAL/Point_3.h>

#include <CGAL/Intersections_3/internal/intersection_3_1_impl.h>
#include <CGAL/Intersections_3/Iso_cuboid_3_Point_3.h>

namespace CGAL {

template<typename K>
bool do_intersect(const CGAL::Bbox_3& a,
                  const Point_3<K>& b) {
  Point_3<K> bl(a.xmin(), a.ymin(),a.zmin()), tr(a.xmax(), a.ymax(),a.zmax());

  Iso_cuboid_3<K> ic(bl,tr);
  return K().do_intersect_3_object()(ic, b);
}


template<typename K>
bool do_intersect(const Point_3<K>& a,
                  const CGAL::Bbox_3& b) {
  return do_intersect(b,a);
}


template<typename K>
typename Intersection_traits<K, typename K::Point_3, CGAL::Bbox_3>::result_type
intersection(const Point_3<K>& a,
             const CGAL::Bbox_3& b)
{
  if (do_intersect(a,b))
    return Intersections::internal::intersection_return<typename K::Intersect_3, typename K::Point_3, CGAL::Bbox_3>(a);
  return Intersections::internal::intersection_return<typename K::Intersect_3, typename K::Point_3, CGAL::Bbox_3>();
}


template<typename K>
typename Intersection_traits<K, CGAL::Bbox_3, typename K::Point_3>::result_type
intersection(const CGAL::Bbox_3& b,
             const Point_3<K>& a)
{
  if (do_intersect(a,b))
    return Intersections::internal::intersection_return<typename K::Intersect_3, CGAL::Bbox_3, typename K::Point_3>(a);
  return Intersections::internal::intersection_return<typename K::Intersect_3, CGAL::Bbox_3, typename K::Point_3>();
}

} // namespace CGAL

#endif // CGAL_INTERSECTIONS_3_BBOX_3_POINT_3_H
