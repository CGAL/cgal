
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
// file          : include/CGAL/Point_2_Segment_2_intersection.h
// source        : intersection_2_1.fw
// author(s)     : Geert-Jan Giezeman
//
// coordinator   : Saarbruecken
//
// ============================================================================


#ifndef CGAL_POINT_2_SEGMENT_2_INTERSECTION_H
#define CGAL_POINT_2_SEGMENT_2_INTERSECTION_H

#include <CGAL/Segment_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Object.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

template <class K>
inline 
bool
do_intersect(const typename CGAL_WRAP(K)::Point_2 &pt, 
	     const typename CGAL_WRAP(K)::Segment_2 &seg,
	     const K&)
{
    return seg.has_on(pt);
}

template <class K>
inline 
bool
do_intersect(const typename CGAL_WRAP(K)::Segment_2 &seg,
	     const typename CGAL_WRAP(K)::Point_2 &pt, 
	     const K&)
{
    return seg.has_on(pt);
}


template <class K>
inline
Object
intersection(const typename CGAL_WRAP(K)::Point_2 &pt, 
	     const typename CGAL_WRAP(K)::Segment_2 &seg, 
	     const K&)
{
    if (do_intersect(pt,seg)) {
        return make_object(pt);
    }
    return Object();
}

template <class K>
inline
Object
intersection( const typename CGAL_WRAP(K)::Segment_2 &seg, 
	      const typename CGAL_WRAP(K)::Point_2 &pt, 
	      const K&)
{
    if (do_intersect(pt,seg)) {
        return make_object(pt);
    }
    return Object();
}

} // namespace CGALi


template <class K>
inline bool
do_intersect(const Segment_2<K> &seg, const Point_2<K> &pt)
{
    return CGALi::do_intersect(pt, seg, K());
}

template <class K>
inline bool
do_intersect(const Point_2<K> &pt, const Segment_2<K> &seg)
{
    return CGALi::do_intersect(pt, seg, K());
}


template <class K>
inline Object
intersection(const Segment_2<K> &seg, const Point_2<K> &pt)
{
    return CGALi::intersection(pt, seg, K());
}

template <class K>
inline Object
intersection(const Point_2<K> &pt, const Segment_2<K> &seg)
{
    return CGALi::intersection(pt, seg, K());
}

CGAL_END_NAMESPACE

#endif
