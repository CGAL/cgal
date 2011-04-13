
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
// file          : include/CGAL/Iso_rectangle_2_Iso_rectangle_2_intersection.h
// source        : intersection_2_2.fw
// author(s)     : Geert-Jan Giezeman
//
// coordinator   : Saarbruecken
//
// ============================================================================


#ifndef CGAL_ISO_RECTANGLE_2_ISO_RECTANGLE_2_INTERSECTION_H
#define CGAL_ISO_RECTANGLE_2_ISO_RECTANGLE_2_INTERSECTION_H

#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Object.h>

CGAL_BEGIN_NAMESPACE

template <class R>
Object
intersection(
    const Iso_rectangle_2<R> &irect1,
    const Iso_rectangle_2<R> &irect2)
{
    const Point_2<R> &min1 = irect1.min();
    const Point_2<R> &min2 = irect2.min();
    const Point_2<R> &max1 = irect1.max();
    const Point_2<R> &max2 = irect2.max();
    typename R::FT minx, miny, maxx, maxy;
    Point_2<R> newmin;
    Point_2<R> newmax;
    minx = (min1.x() >= min2.x()) ? min1.x() : min2.x();
    maxx = (max1.x() <= max2.x()) ? max1.x() : max2.x();
    if (maxx < minx)
        return Object();
    miny = (min1.y() >= min2.y()) ? min1.y() : min2.y();
    maxy = (max1.y() <= max2.y()) ? max1.y() : max2.y();
    if (maxy < miny)
        return Object();
    if (R::FT_denominator(minx) == R::FT_denominator(miny)) {
        newmin = Point_2<R>(R::FT_numerator(minx), R::FT_numerator(miny),
                    R::FT_denominator(minx));
    } else {
        newmin = Point_2<R>(R::FT_numerator(minx)*R::FT_denominator(miny),
                    R::FT_numerator(miny)*R::FT_denominator(minx),
                    R::FT_denominator(minx) * R::FT_denominator(miny));
    }
    if (R::FT_denominator(maxx) == R::FT_denominator(maxy)) {
        newmax = Point_2<R>(R::FT_numerator(maxx), R::FT_numerator(maxy),
                    R::FT_denominator(maxx));
    } else {
        newmax = Point_2<R>(R::FT_numerator(maxx)*R::FT_denominator(maxy),
                    R::FT_numerator(maxy)*R::FT_denominator(maxx),
                    R::FT_denominator(maxx) * R::FT_denominator(maxy));
    }
    return make_object(Iso_rectangle_2<R>(newmin, newmax));
}

template <class R>
inline bool
do_intersect(
    const Iso_rectangle_2<R> &irect1,
    const Iso_rectangle_2<R> &irect2)
{
    Object obj(intersection(irect1, irect2));
    Iso_rectangle_2<R> irect;
    return (assign(irect, obj));
}

CGAL_END_NAMESPACE

#endif
