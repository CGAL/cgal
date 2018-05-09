// Copyright (c) 2010 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
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
             const CGAL::Bbox_3& b
             ) {
    if (do_intersect(a,b)) {
      return Intersections::internal::intersection_return<typename K::Intersect_3, typename K::Point_3, CGAL::Bbox_3>(a);
  }
    return Intersections::internal::intersection_return<typename K::Intersect_3, typename K::Point_3, CGAL::Bbox_3>();
}

  
template<typename K>
typename Intersection_traits<K, CGAL::Bbox_3, typename K::Point_3>::result_type
intersection(const CGAL::Bbox_3& b,
             const Point_3<K>& a) {
    if (do_intersect(a,b)) {
      return Intersections::internal::intersection_return<typename K::Intersect_3, CGAL::Bbox_3, typename K::Point_3>(a);
  }
    return Intersections::internal::intersection_return<typename K::Intersect_3, CGAL::Bbox_3, typename K::Point_3>();
}

} // namespace CGAL

#endif // CGAL_INTERSECTIONS_3_BBOX_3_POINT_3_H
