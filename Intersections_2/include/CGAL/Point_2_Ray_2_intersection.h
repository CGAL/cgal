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
// 
//
// Author(s)     : Geert-Jan Giezeman


#ifndef CGAL_POINT_2_RAY_2_INTERSECTION_H
#define CGAL_POINT_2_RAY_2_INTERSECTION_H

#include <CGAL/Ray_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Object.h>

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
Object
intersection(const typename K::Point_2 &pt, 
	     const typename K::Ray_2 &ray,
	     const K& k)
{
  if (do_intersect(pt,ray, k)) {
    return make_object(pt);
  }
  return Object();
}

template <class K>
Object
intersection(const typename K::Ray_2 &ray,
	     const typename K::Point_2 &pt, 
	     const K& k)
{
  if (do_intersect(pt,ray, k)) {
    return make_object(pt);
  }
  return Object();
}

} // namespace internal


template <class K>
inline
bool
do_intersect(const Ray_2<K> &ray, const Point_2<K> &pt)
{
  typedef typename K::Do_intersect_2 Do_intersect;
  return Do_intersect()(pt, ray);
}

template <class K>
inline
bool
do_intersect(const Point_2<K> &pt, const Ray_2<K> &ray)
{
  typedef typename K::Do_intersect_2 Do_intersect;
  return Do_intersect()(pt, ray);
}


template <class K>
inline Object
intersection(const Ray_2<K> &ray, const Point_2<K> &pt)
{
  typedef typename K::Intersect_2 Intersect;
  return Intersect()(pt, ray);
}

template <class K>
inline Object
intersection(const Point_2<K> &pt, const Ray_2<K> &ray)
{
  typedef typename K::Intersect_2 Intersect;
  return Intersect()(pt, ray);
}

} //namespace CGAL

#endif
