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

#ifndef CGAL_INTERSECTIONS_2_BBOX_2_POINT_2_H
#define CGAL_INTERSECTIONS_2_BBOX_2_POINT_2_H

#include <CGAL/Bbox_2.h>
#include <CGAL/Point_2.h>

#include <CGAL/Intersections_2/Iso_rectangle_2_Point_2.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
inline bool do_intersect(const Bbox_2 &bbox,
                         const Point_2<K> &pt,
                         const K& k)
{
  Point_2<K> bl(bbox.xmin(), bbox.ymin()),
             tr(bbox.xmax(), bbox.ymax());
  Iso_rectangle_2<K> ic(bl,tr);

  return k.do_intersect_2_object()(ic, pt);
}

template <class K>
inline bool do_intersect(const Point_2<K> &pt,
                         const Bbox_2& bbox,
                         const K& k)
{
  return do_intersect(bbox, pt, k);
}

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
bool do_intersect(const CGAL::Bbox_2& a,
                  const Point_2<K>& b)
{
  return Intersections::internal::do_intersect(a,b,K());
}

template<typename K>
bool do_intersect(const Point_2<K>& a,
                  const CGAL::Bbox_2& b)
{
  return Intersections::internal::do_intersect(b,a,K());
}

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
