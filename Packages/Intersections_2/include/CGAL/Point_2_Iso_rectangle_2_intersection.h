
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
// file          : include/CGAL/Point_2_Iso_rectangle_2_intersection.h
// source        : intersection_2_2.fw
// author(s)     : Geert-Jan Giezeman
//
// coordinator   : Saarbruecken
//
// ============================================================================


#ifndef CGAL_POINT_2_ISO_RECTANGLE_2_INTERSECTION_H
#define CGAL_POINT_2_ISO_RECTANGLE_2_INTERSECTION_H

#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Object.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

template <class K>
inline 
bool
do_intersect(const typename CGAL_WRAP(K)::Point_2 &pt,
	     const typename CGAL_WRAP(K)::Iso_rectangle_2 &iso,
	     const K&)
{
    return !iso.has_on_unbounded_side(pt);
}

template <class K>
inline 
bool
do_intersect(const typename CGAL_WRAP(K)::Iso_rectangle_2 &iso,
	     const typename CGAL_WRAP(K)::Point_2 &pt,
	     const K&)
{
    return !iso.has_on_unbounded_side(pt);
}

template <class K>
Object
intersection(const typename CGAL_WRAP(K)::Point_2 &pt,
	     const typename CGAL_WRAP(K)::Iso_rectangle_2 &iso,
	     const K& k)
{
  if (CGALi::do_intersect(pt,iso,k)) {
    return make_object(pt);
    }
    return Object();
}


template <class K>
Object
intersection(const typename CGAL_WRAP(K)::Iso_rectangle_2 &iso,
	     const typename CGAL_WRAP(K)::Point_2 &pt,
	     const K& k)
{
  if (CGALi::do_intersect(pt,iso,k)) {
    return make_object(pt);
    }
    return Object();
}

} // namespace CGALi


template <class K>
inline 
bool
do_intersect(const Iso_rectangle_2<K> &iso,
	     const Point_2<K> &pt)
{
  return CGALi::do_intersect(pt, iso, K());
}

template <class K>
inline 
bool
do_intersect(const Point_2<K> &pt,
	     const Iso_rectangle_2<K> &iso)
{
  return CGALi::do_intersect(pt, iso, K());
}

template <class K>
inline 
Object
intersection(const Iso_rectangle_2<K> &iso,
	     const Point_2<K> &pt)
{
  return CGALi::intersection(pt, iso, K());;
}
template <class K>
inline 
Object
intersection(const Point_2<K> &pt,
	     const Iso_rectangle_2<K> &iso)
{
  return CGALi::intersection(pt, iso, K());;
}
CGAL_END_NAMESPACE

#endif
