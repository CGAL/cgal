
// Copyright (c) 2000  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Geert-Jan Giezeman


#ifndef CGAL_POINT_2_RAY_2_INTERSECTION_H
#define CGAL_POINT_2_RAY_2_INTERSECTION_H

#include <CGAL/Ray_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Object.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

template <class K>
inline 
bool
do_intersect(const typename CGAL_WRAP(K)::Point_2 &pt, 
	     const typename CGAL_WRAP(K)::Ray_2 &ray,
	     const K&)
{
  return ray.has_on(pt);
}


template <class K>
inline 
bool
do_intersect(const typename CGAL_WRAP(K)::Ray_2 &ray,
	     const typename CGAL_WRAP(K)::Point_2 &pt, 
	     const K&)
{
  return ray.has_on(pt);
}


template <class K>
Object
intersection(const typename CGAL_WRAP(K)::Point_2 &pt, 
	     const typename CGAL_WRAP(K)::Ray_2 &ray,
	     const K& k)
{
  if (do_intersect(pt,ray, k)) {
    return make_object(pt);
  }
  return Object();
}

template <class K>
Object
intersection(const typename CGAL_WRAP(K)::Ray_2 &ray,
	     const typename CGAL_WRAP(K)::Point_2 &pt, 
	     const K& k)
{
  if (do_intersect(pt,ray, k)) {
    return make_object(pt);
  }
  return Object();
}

} // namespace CGALi


template <class K>
inline
bool
do_intersect(const Ray_2<K> &ray, const Point_2<K> &pt)
{
  return typename K::do_intersect_2()(pt, ray);
}

template <class K>
inline
bool
do_intersect(const Point_2<K> &pt, const Ray_2<K> &ray)
{
  return typename K::Do_intersect_2()(pt, ray);
}


template <class K>
inline Object
intersection(const Ray_2<K> &ray, const Point_2<K> &pt)
{
  return typename K::Intersect_2()(pt, ray);
}

template <class K>
inline Object
intersection(const Point_2<K> &pt, const Ray_2<K> &ray)
{
  return typename K::Intersect_2()(pt, ray);
}

CGAL_END_NAMESPACE

#endif
