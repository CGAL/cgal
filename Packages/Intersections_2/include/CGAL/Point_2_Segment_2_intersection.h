
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

CGAL_BEGIN_NAMESPACE

template <class R>
inline bool
do_intersect(const Point_2<R> &pt, const Segment_2<R> &seg)
{
    return seg.has_on(pt);
}

CGAL_END_NAMESPACE

#include <CGAL/Object.h>

CGAL_BEGIN_NAMESPACE

template <class R>
Object
intersection(const Point_2<R> &pt, const Segment_2<R> &seg)
{
    if (do_intersect(pt,seg)) {
        return Object(new Wrapper< Point_2<R> >(pt));
    }
    return Object();
}


template <class R>
inline bool
do_intersect(const Segment_2<R> &seg, const Point_2<R> &pt)
{
    return seg.has_on(pt);
}


template <class R>
inline Object
intersection(const Segment_2<R> &seg, const Point_2<R> &pt)
{
    if (do_intersect(pt,seg)) {
        return Object(new Wrapper< Point_2<R> >(pt));
    }
    return Object();
}

CGAL_END_NAMESPACE

#endif
