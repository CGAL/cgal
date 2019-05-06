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

#ifndef CGAL_INTERSECTIONS_2_BBOX_2_POINT_2_H
#define CGAL_INTERSECTIONS_2_BBOX_2_POINT_2_H

#include <CGAL/Bbox_2.h>
#include <CGAL/Point_2.h>

#include <CGAL/Intersections_2/Iso_rectangle_2_Point_2.h>

namespace CGAL {

template<typename K>
bool do_intersect(const CGAL::Bbox_2& a,
                  const Point_2<K>& b)
{
  Point_2<K> bl(a.xmin(), a.ymin()), tr(a.xmax(), a.ymax());

  Iso_rectangle_2<K> ic(bl,tr);
  return K().do_intersect_2_object()(ic, b);
}


template<typename K>
bool do_intersect(const Point_2<K>& a,
                  const CGAL::Bbox_2& b)
{
  return do_intersect(b,a);
}

namespace Intersections {

namespace internal {

template<typename K>
typename CGAL::Intersection_traits<K, typename K::Point_2, CGAL::Bbox_2>::result_type
intersection(const Point_2<K>& a,
             const CGAL::Bbox_2& b)
{
  if (do_intersect(a,b))
    return Intersections::internal::intersection_return<typename K::Intersect_2, typename K::Point_2, CGAL::Bbox_2>(a);
  return Intersections::internal::intersection_return<typename K::Intersect_2, typename K::Point_2, CGAL::Bbox_2>();
}


template<typename K>
typename CGAL::Intersection_traits<K, CGAL::Bbox_2, typename K::Point_2>::result_type
intersection(const CGAL::Bbox_2& b,
             const typename K::Point_2  & a,
             const K& /*k*/ )
{
  if (do_intersect(a,b))
    return Intersections::internal::intersection_return<typename K::Intersect_2, CGAL::Bbox_2, typename K::Point_2>(a);
  return Intersections::internal::intersection_return<typename K::Intersect_2, CGAL::Bbox_2, typename K::Point_2>();
}

} // namespace internal
} // namespace Intersections

template<typename K>
typename CGAL::Intersection_traits<K, Bbox_2, Point_2<K> >::result_type
intersection(const Bbox_2& b,
             const Point_2<K>  & a)
{
  return Intersections::internal::intersection(b,a,K());
}

template<typename K>
typename CGAL::Intersection_traits<K, Bbox_2, Point_2<K> >::result_type
  intersection(const Point_2<K>  & a,
               const Bbox_2& b)
{
  return Intersections::internal::intersection(b,a,K());
}
} // namespace CGAL

#endif // CGAL_INTERSECTIONS_2_BBOX_2_POINT_2_H
