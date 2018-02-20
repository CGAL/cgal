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

#ifndef CGAL_INTERSECTIONS_3_POINT_3_TETRAHEDRON_3_H
#define CGAL_INTERSECTIONS_3_POINT_3_TETRAHEDRON_3_H

#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/Intersection_traits_3.h>

namespace CGAL {
  
namespace Intersections {

namespace internal {

template <class K>
inline
bool
do_intersect(const typename K::Point_3 &pt,
             const typename K::Tetrahedron_3 &tetrahedron,
             const K&)
{
  return ! tetrahedron.has_on_unbounded_side(pt);
}


template <class K>
inline
bool
do_intersect(const typename K::Tetrahedron_3 &tetrahedron,
             const typename K::Point_3 &pt,
             const K&)
{
  return tetrahedron.has_on_unbounded_side(pt);
}


template <class K>
typename CGAL::Intersection_traits
<K, typename K::Point_3, typename K::Tetrahedron_3>::result_type
intersection(const typename K::Point_3 &pt,
             const typename K::Tetrahedron_3 &tetrahedron,
             const K& k)
{
  if (do_intersect(pt,tetrahedron, k)) {
    return intersection_return<typename K::Intersect_3, typename K::Point_3, typename K::Tetrahedron_3>(pt);
  }
  return intersection_return<typename K::Intersect_3, typename K::Point_3, typename K::Tetrahedron_3>();
}

template <class K>
typename CGAL::Intersection_traits
<K, typename K::Tetrahedron_3, typename K::Point_3>::result_type
intersection(const typename K::Tetrahedron_3 &tetrahedron,
             const typename K::Point_3 &pt,
             const K& k)
{
  return internal::intersection(pt, tetrahedron, k);
}

} // namespace internal
} // namespace Intersections
  
CGAL_INTERSECTION_FUNCTION(Point_3, Tetrahedron_3, 3)
CGAL_DO_INTERSECT_FUNCTION(Point_3, Tetrahedron_3, 3)


} //namespace CGAL
#endif // CGAL_INTERSECTIONS_3_POINT_3_TETRAHEDRON_3_H
