// Copyright (c) 2008  INRIA Sophia-Antipolis (France), ETH Zurich (Switzerland).
// Copyright (c) 2010, 2014  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Camille Wormser, Jane Tournois, Pierre Alliez

#ifndef CGAL_INTERNAL_INTERSECTIONS_3_LINE_3_LINE_3_INTERSECTION_H
#define CGAL_INTERNAL_INTERSECTIONS_3_LINE_3_LINE_3_INTERSECTION_H

#include <CGAL/Intersection_traits_3.h>

#include <CGAL/number_utils.h>
#include <CGAL/Uncertain.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
typename Intersection_traits<K, typename K::Line_3, typename K::Line_3>::result_type
intersection(const typename K::Line_3& l1,
             const typename K::Line_3& l2,
             const K& k)
{
  typedef typename K::FT           FT;
  typedef typename K::Point_3      Point_3;
  typedef typename K::Vector_3     Vector_3;

  if(k.has_on_3_object()(l1, l2.point()))
  {
    const Vector_3& v1 = l1.to_vector();
    const Vector_3& v2 = l2.to_vector();
    if((v1.x() * v2.y() == v1.y() * v2.x()) &&
       (v1.x() * v2.z() == v1.z() * v2.x()) &&
       (v1.y() * v2.z() == v1.z() * v2.y()))
      return intersection_return<typename K::Intersect_3, typename K::Line_3, typename K::Line_3>(l1);
  }

  if(k.are_parallel_3_object()(l1,l2))
    return intersection_return<typename K::Intersect_3, typename K::Line_3, typename K::Line_3>();

  const Point_3& p1 = l1.point();
  const Point_3& p3 = l2.point();
  const Vector_3& v1 = l1.to_vector();
  const Vector_3& v2 = l2.to_vector();
  const Point_3 p2 = p1 + v1;
  const Point_3 p4 = p2 + v2;
  if(!k.coplanar_3_object()(p1,p2,p3,p4))
    return intersection_return<typename K::Intersect_3, typename K::Line_3, typename K::Line_3>();

  const Vector_3 v3 = p3 - p1;
  const Vector_3 v3v2 = cross_product(v3,v2);
  const Vector_3 v1v2 = cross_product(v1,v2);
  const FT sl = v1v2.squared_length();
  if(certainly(is_zero(sl)))
    return intersection_return<typename K::Intersect_3, typename K::Line_3, typename K::Line_3>();

  const FT t = ((v3v2.x()*v1v2.x()) + (v3v2.y()*v1v2.y()) + (v3v2.z()*v1v2.z())) / sl;

  return intersection_return<typename K::Intersect_3, typename K::Line_3, typename K::Line_3>(p1 + (v1 * t));
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_INTERNAL_INTERSECTIONS_3_LINE_3_LINE_3_INTERSECTION_H
