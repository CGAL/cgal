// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
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
// file          : include/CGAL/intersection_3_1.h
// source        : web/intersection_3.fw
// author(s)     : Geert-Jan Giezeman <geert@cs.uu.nl>
//
// coordinator   : Saarbruecken
//
// ============================================================================


#ifndef CGAL_INTERSECTION_3_1_H
#define CGAL_INTERSECTION_3_1_H



#include <CGAL/Object.h>
#include <CGAL/bbox_intersection_3.h>


CGAL_BEGIN_NAMESPACE

template <class R>
Object
intersection(const Plane_3<R> &plane1, const Plane_3<R>&plane2);

template <class R>
inline bool
do_intersect(const Plane_3<R> &plane1, const Plane_3<R>&plane2)
{
    return ! intersection(plane1, plane2).is_empty();
}
CGAL_END_NAMESPACE




CGAL_BEGIN_NAMESPACE

template <class R>
Object
intersection(const Plane_3<R> &plane, const Line_3<R>&line);

template <class R>
inline
Object
intersection(const Line_3<R>&line, const Plane_3<R> &plane)
{
    return intersection(plane,line);
}

template <class R>
bool
do_intersect(const Plane_3<R> &p2, const Line_3<R> &p1);


template <class R>
inline bool
do_intersect(
    const Line_3<R> &p1,
    const Plane_3<R> &p2)
{
    return do_intersect(p2,p1);
}

CGAL_END_NAMESPACE




CGAL_BEGIN_NAMESPACE

template <class R>
Object
intersection(const Plane_3<R> &plane, const Ray_3<R>&ray);

template <class R>
inline
Object
intersection(const Ray_3<R>&ray, const Plane_3<R> &plane)
{
    return intersection(plane,ray);
}

template <class R>
bool
do_intersect(const Plane_3<R> &p1, const Ray_3<R> &p2);


template <class R>
inline bool
do_intersect(
    const Ray_3<R> &p1,
    const Plane_3<R> &p2)
{
    return do_intersect(p2,p1);
}
CGAL_END_NAMESPACE



CGAL_BEGIN_NAMESPACE

template <class R>
Object
intersection(const Plane_3<R> &plane, const Segment_3<R>&seg);

template <class R>
inline
Object
intersection(const Segment_3<R>&seg, const Plane_3<R> &plane)
{
    return intersection(plane,seg);
}

template <class R>
bool
do_intersect(const Plane_3<R> &p1, const Segment_3<R> &p2);


template <class R>
inline bool
do_intersect(
    const Segment_3<R> &p1,
    const Plane_3<R> &p2)
{
    return do_intersect(p2,p1);
}

CGAL_END_NAMESPACE




CGAL_BEGIN_NAMESPACE

template <class R>
Object
intersection(const Line_3<R> &line,
        const Bbox_3 &box) ;

template <class R>
inline Object
intersection(const Bbox_3 &box,
        const Line_3<R> &line)
{
    return intersection(line, box);
}

CGAL_END_NAMESPACE



CGAL_BEGIN_NAMESPACE

template <class R>
Object
intersection(const Ray_3<R> &ray,
        const Bbox_3 &box) ;

template <class R>
inline Object
intersection(const Bbox_3 &box,
        const Ray_3<R> &ray)
{
    return intersection(ray, box);
}

CGAL_END_NAMESPACE



CGAL_BEGIN_NAMESPACE

template <class R>
Object
intersection(const Segment_3<R> &seg,
        const Bbox_3 &box) ;

template <class R>
inline Object
intersection(const Bbox_3 &box,
        const Segment_3<R> &seg)
{
    return intersection(seg, box);
}
CGAL_END_NAMESPACE



CGAL_BEGIN_NAMESPACE

template <class R>
Object
intersection(const Line_3<R> &line,
        const Iso_cuboid_3<R> &box) ;

template <class R>
inline Object
intersection(const Iso_cuboid_3<R> &box,
        const Line_3<R> &line)
{
    return intersection(line, box);
}

CGAL_END_NAMESPACE



CGAL_BEGIN_NAMESPACE

template <class R>
Object
intersection(const Ray_3<R> &ray,
        const Iso_cuboid_3<R> &box) ;

template <class R>
inline Object
intersection(const Iso_cuboid_3<R> &box,
        const Ray_3<R> &ray)
{
    return intersection(ray, box);
}

CGAL_END_NAMESPACE



CGAL_BEGIN_NAMESPACE

template <class R>
Object
intersection(const Segment_3<R> &seg,
        const Iso_cuboid_3<R> &box) ;

template <class R>
inline Object
intersection(const Iso_cuboid_3<R> &box,
        const Segment_3<R> &seg)
{
    return intersection(seg, box);
}

CGAL_END_NAMESPACE



CGAL_BEGIN_NAMESPACE

template <class R>
Object
intersection(const Iso_cuboid_3<R> &box1,
        const Iso_cuboid_3<R> &box2) ;

CGAL_END_NAMESPACE




#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/intersection_3_1.C>
#endif

#endif // CGAL_INTERSECTION_3_1_H
