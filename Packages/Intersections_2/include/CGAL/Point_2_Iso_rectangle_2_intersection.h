
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

CGAL_BEGIN_NAMESPACE

template <class R>
inline bool
do_intersect(
    const Point_2<R> &pt,
    const Iso_rectangle_2<R> &iso)
{
    return !iso.has_on_unbounded_side(pt);
}

CGAL_END_NAMESPACE

#include <CGAL/Object.h>

CGAL_BEGIN_NAMESPACE

template <class R>
Object
intersection(
    const Point_2<R> &pt,
    const Iso_rectangle_2<R> &iso)
{
    if (do_intersect(pt,iso)) {
        return Object(new Wrapper< Point_2<R> >(pt));
    }
    return Object();
}

template <class R>
inline bool
do_intersect(
    const Iso_rectangle_2<R> &iso,
    const Point_2<R> &pt)
{
    return !iso.has_on_unbounded_side(pt);
}


template <class R>
inline Object
intersection(
    const Iso_rectangle_2<R> &iso,
    const Point_2<R> &pt)
{
    if (do_intersect(pt, iso)) {
        return Object(new Wrapper< Point_2<R> >(pt));
    }
    return Object();
}

CGAL_END_NAMESPACE

#endif
