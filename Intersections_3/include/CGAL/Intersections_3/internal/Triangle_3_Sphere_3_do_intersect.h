// Copyright (c) 2018 GeometryFactory Sarl
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
// Author(s)     : Andreas Fabri

#ifndef CGAL_TRIANGLE_3_SPHERE_3_DO_INTERSECT_H
#define CGAL_TRIANGLE_3_SPHERE_3_DO_INTERSECT_H

#include <CGAL/squared_distance_3_2.h>
#include <CGAL/Intersection_traits_3.h>

namespace CGAL {

template <class K>
class Triangle_3;

template <class K>
class Sphere_3;

template <class K>
class Line_3;

namespace Intersections {

namespace internal {

template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Sphere_3 &sp,
             const typename K::Triangle_3 &tr,
	     const K & /* k */)
{
  return squared_distance(sp.center(), tr) <= sp.squared_radius();
}

template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Triangle_3 &tr,
             const typename K::Sphere_3 &sp,
	     const K & /* k */)
{
  return squared_distance(sp.center(), tr) <= sp.squared_radius();
}
template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Sphere_3 &sp,
             const typename K::Line_3 &lin,
	     const K & /* k */)
{
  return squared_distance(sp.center(), lin) <= sp.squared_radius();
}


template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Line_3 &lin,
	     const typename K::Sphere_3 &sp,
	     const K & /* k */)
{
  return squared_distance(sp.center(), lin) <= sp.squared_radius();
}



template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Sphere_3 &sp,
             const typename K::Ray_3 &lin,
	     const K & /* k */)
{
  return squared_distance(sp.center(), lin) <= sp.squared_radius();
}


template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Ray_3 &lin,
	     const typename K::Sphere_3 &sp,
	     const K & /* k */)
{
  return squared_distance(sp.center(), lin) <= sp.squared_radius();
}

template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Sphere_3 &sp,
             const typename K::Segment_3 &lin,
	     const K & /* k */)
{
  return squared_distance(sp.center(), lin) <= sp.squared_radius();
}


template <class K>
inline
typename K::Boolean
do_intersect(const typename K::Segment_3 &lin,
	     const typename K::Sphere_3 &sp,
	     const K & /* k */)
{
  return squared_distance(sp.center(), lin) <= sp.squared_radius();
}

} // namespace internal
} // namespace Intersections
} // namespace CGAL

#endif // CGAL_TRIANGLE_3_SPHERE_3_DO_INTERSECT_H
