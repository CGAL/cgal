
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
// file          : include/CGAL/Point_2_Point_2_intersection.h
// source        : intersection_2_1.fw
// author(s)     : Geert-Jan Giezeman
//
// coordinator   : Saarbruecken
//
// ============================================================================


#ifndef CGAL_POINT_2_POINT_2_INTERSECTION_H
#define CGAL_POINT_2_POINT_2_INTERSECTION_H

#include <CGAL/Point_2.h>

#include <CGAL/Object.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

template <class K>
inline bool
do_intersect(const typename CGAL_WRAP(K)::Point_2 &pt1, 
	     const typename CGAL_WRAP(K)::Point_2 &pt2)
{
    return pt1 == pt2;
}
template <class K>
Object
intersection(const typename CGAL_WRAP(K)::Point_2 &pt1, 
	     const typename CGAL_WRAP(K)::Point_2 &pt2)
{
    if (pt1 == pt2) {
        return make_object(pt1);
    }
    return Object();
}

}// namespace CGALi


template <class K>
inline 
bool
do_intersect(const Point_2<K> &pt1, const Point_2<K> &pt2)
{
  return CGALi::do_intersect(pt1, pt2, K());
}


template <class K>
inline
Object
intersection(const Point_2<K> &pt1, const Point_2<K> &pt2)
{
  return CGALi::intersection(pt1, pt2, K());
}

CGAL_END_NAMESPACE

#endif
