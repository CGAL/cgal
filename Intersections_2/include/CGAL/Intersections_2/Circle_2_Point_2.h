// Copyright (c) 2018  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Maxime Gimeno

#ifndef CGAL_INTERSECTIONS_2_POINT_2_CIRCLE_2_H
#define CGAL_INTERSECTIONS_2_POINT_2_CIRCLE_2_H

#include <CGAL/Circle_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Intersection_traits_2.h>

namespace CGAL {
  
namespace Intersections {

namespace internal {

template <class K>
inline
bool
do_intersect(const typename K::Point_2 &pt,
             const typename K::Circle_2 &circle,
             const K&)
{
  return circle.has_on_boundary(pt);
}


template <class K>
inline
bool
do_intersect(const typename K::Circle_2 &circle,
             const typename K::Point_2 &pt,
             const K&)
{
  return circle.has_on_boundary(pt);
}


template <class K>
typename CGAL::Intersection_traits
<K, typename K::Point_2, typename K::Circle_2>::result_type
intersection(const typename K::Point_2 &pt,
             const typename K::Circle_2 &circle,
             const K& k)
{
  if (do_intersect(pt,circle, k)) {
    return intersection_return<typename K::Intersect_2, typename K::Point_2, typename K::Circle_2>(pt);
  }
  return intersection_return<typename K::Intersect_2, typename K::Point_2, typename K::Circle_2>();
}

template <class K>
typename CGAL::Intersection_traits
<K, typename K::Circle_2, typename K::Point_2>::result_type
intersection(const typename K::Circle_2 &circle,
             const typename K::Point_2 &pt,
             const K& k)
{
  return internal::intersection(pt, circle, k);
}

} // namespace internal
} // namespace Intersections
  
CGAL_INTERSECTION_FUNCTION(Point_2, Circle_2, 2)
CGAL_DO_INTERSECT_FUNCTION(Circle_2, Point_2, 2)

} //namespace CGAL
#endif // CGAL_INTERSECTIONS_2_POINT_2_CIRCLE_2_H
