
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

CGAL_BEGIN_NAMESPACE

template <class R>
inline bool
do_intersect(const Point_2<R> &pt1, const Point_2<R> &pt2)
{
    return pt1 == pt2;
}

CGAL_END_NAMESPACE

#include <CGAL/Object.h>

CGAL_BEGIN_NAMESPACE

template <class R>
Object
intersection(const Point_2<R> &pt1, const Point_2<R> &pt2)
{
    if (pt1 == pt2) {
        return Object(new Wrapper< Point_2<R> >(pt1));
    }
    return Object();
}

CGAL_END_NAMESPACE

#endif
