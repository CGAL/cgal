
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
// file          : include/CGAL/Point_2_Line_2_intersection.h
// source        : intersection_2_1.fw
// author(s)     : Geert-Jan Giezeman
//
// coordinator   : Saarbruecken
//
// ============================================================================


#ifndef CGAL_POINT_2_LINE_2_INTERSECTION_H
#define CGAL_POINT_2_LINE_2_INTERSECTION_H

#include <CGAL/Line_2.h>
#include <CGAL/Point_2.h>

#include <CGAL/Object.h>
CGAL_BEGIN_NAMESPACE

namespace CGALi {

template <class K>
inline bool
do_intersect(const typename CGAL_WRAP(K)::Point_2 &pt, 
	     const typename CGAL_WRAP(K)::Line_2 &line,
	     const K&)
{
    return line.has_on(pt);
}

template <class K>
Object
intersection(const typename CGAL_WRAP(K)::Point_2 &pt, 
	     const typename CGAL_WRAP(K)::Line_2 &line,
	     const K& k)
{
    if (do_intersect(pt,line, k)) {
        return make_object(pt);
    }
    return Object();
}

} // namespace CGALi

template <class K>
inline 
bool
do_intersect(const Line_2<K> &line, 
	     const Point_2<K> &pt)
{
    return CGALi::do_intersect(pt, line, K());
}

template <class K>
inline 
bool
do_intersect(const Point_2<K> &pt,
	     const Line_2<K> &line)
{
    return CGALi::do_intersect(pt, line, K());
}

template <class K>
inline 
Object
intersection(const Line_2<K> &line, 
	     const Point_2<K> &pt)
{
  return CGALi::intersection(pt, line, K());
}

template <class K>
inline 
Object
intersection(const Point_2<K> &pt,
	     const Line_2<K> &line)
{
  return CGALi::intersection(pt, line, K());
}
CGAL_END_NAMESPACE

#endif
