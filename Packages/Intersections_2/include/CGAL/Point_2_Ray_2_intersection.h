
// ============================================================================
//
// Copyright (c) 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision:  $
// release_date  : $CGAL_Date:  $
//
// file          : include/CGAL/Point_2_Ray_2_intersection.h
// source        : intersection_2_1.fw
// author(s)     : Geert-Jan Giezeman
//
// coordinator   : Saarbruecken
//
// ============================================================================


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
  return CGALi::do_intersect(pt, ray, K());
}

template <class K>
inline
bool
do_intersect(const Point_2<K> &pt, const Ray_2<K> &ray)
{
  return CGALi::do_intersect(pt, ray, K());
}


template <class K>
inline Object
intersection(const Ray_2<K> &ray, const Point_2<K> &pt)
{
  return CGALi::intersection(pt, ray, K());
}

template <class K>
inline Object
intersection(const Point_2<K> &pt, const Ray_2<K> &ray)
{
  return CGALi::intersection(pt, ray, K());
}

CGAL_END_NAMESPACE

#endif
