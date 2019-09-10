// Copyright (c) 2000  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
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
// Author(s)     : Geert-Jan Giezeman


#ifndef CGAL_POINT_2_RAY_2_INTERSECTION_H
#define CGAL_POINT_2_RAY_2_INTERSECTION_H

#include <CGAL/Ray_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Intersection_traits_2.h>

namespace CGAL {

namespace internal {

template <class K>
inline 
bool
do_intersect(const typename K::Point_2 &pt, 
	     const typename K::Ray_2 &ray,
	     const K&)
{
  return ray.has_on(pt);
}


template <class K>
inline 
bool
do_intersect(const typename K::Ray_2 &ray,
	     const typename K::Point_2 &pt, 
	     const K&)
{
  return ray.has_on(pt);
}


template <class K>
typename CGAL::Intersection_traits
<K, typename K::Point_2, typename K::Ray_2>::result_type
intersection(const typename K::Point_2 &pt, 
	     const typename K::Ray_2 &ray,
	     const K& k)
{
  if (do_intersect(pt,ray, k)) {
    return intersection_return<typename K::Intersect_2, typename K::Point_2, typename K::Ray_2>(pt);
  }
  return intersection_return<typename K::Intersect_2, typename K::Point_2, typename K::Ray_2>();
}

template <class K>
typename CGAL::Intersection_traits
<K, typename K::Ray_2, typename K::Point_2>::result_type
intersection(const typename K::Ray_2 &ray,
	     const typename K::Point_2 &pt, 
	     const K& k)
{
  return internal::intersection(pt, ray, k);
}

} // namespace internal

CGAL_INTERSECTION_FUNCTION(Point_2, Ray_2, 2)
CGAL_DO_INTERSECT_FUNCTION(Point_2, Ray_2, 2)


} //namespace CGAL

#endif
