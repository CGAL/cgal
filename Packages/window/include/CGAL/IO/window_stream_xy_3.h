// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : IO/window_stream_xy_3.h
// package       : window
// chapter       : $CGAL_Chapter: Window Stream $
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner 
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
//
// Supports the display of 3D objects in the 2d CGAL_Window_stream.
// ======================================================================

// Note: This file could be included multiple times. Thus, the
// usual protection against multiple inclusion is not used.
// The following section protects itself against multiple inclusion
// with a non-standard macro name to hinder redundant protection
// guards in other files to exclude this file from inclusion.

// (This section is currently empty.)
#ifndef CGAL_IO_WINDOW_STREAM_XY_3_H_1
#define CGAL_IO_WINDOW_STREAM_XY_3_H_1 1
#endif // CGAL_IO_WINDOW_STREAM_XY_3_H_1 //


//  Each of the following operators is individually
//  protected against multiple inclusion.

// Check first, whether the related 2d files are included.
// -------------------------------------------------------
#if defined(CGAL_POINT_3_H) && ! defined(CGAL_POINT_2_H)
#include <CGAL/Point_2.h>
#endif
#if defined(CGAL_VECTOR_3_H) && ! defined(CGAL_VECTOR_2_H)
#include <CGAL/Vector_2.h>
#endif
#if defined(CGAL_DIRECTION_3_H) && ! defined(CGAL_DIRECTION_2_H)
#include <CGAL/Direction_2.h>
#endif
#if defined(CGAL_LINE_3_H) && ! defined(CGAL_LINE_2_H)
#include <CGAL/Line_2.h>
#endif
#if defined(CGAL_RAY_3_H) && ! defined(CGAL_RAY_2_H)
#include <CGAL/Ray_2.h>
#endif
#if defined(CGAL_SEGMENT_3_H) && ! defined(CGAL_SEGMENT_2_H)
#include <CGAL/Segment_2.h>
#endif
#if defined(CGAL_TRIANGLE_3_H) && ! defined(CGAL_TRIANGLE_2_H)
#include <CGAL/Triangle_2.h>
#endif
#if defined(CGAL_TETRAHEDRON_3_H) && ! defined(CGAL_SEGMENT_2_H)
#include <CGAL/Segment_2.h>
#endif
#if defined(CGAL_TETRAHEDRON_3_H) && ! defined(CGAL_TRIANGLE_2_H)
#include <CGAL/Triangle_2.h>
#endif
#if defined(CGAL_BBOX_3_H) && ! defined(CGAL_BBOX_2_H)
#include <CGAL/Bbox_2.h>
#endif

#if defined(CGAL_LINE_3_H) || defined(CGAL_RAY_3_H) \
    || defined(CGAL_SEGMENT_3_H) || defined(CGAL_TRIANGLE_3_H) \
    || defined(CGAL_TETRAHEDRON_3_H)
#ifndef CGAL_POINT_2_H
#include <CGAL/Point_2.h>
#endif
#ifndef CGAL_POINT_3_H
#include <CGAL/Point_3.h>
#endif
#endif

// Define necessary 2d stream operators.
// -------------------------------------
#ifndef CGAL_IO_WINDOW_STREAM_H
#include <CGAL/IO/Window_stream.h>
#endif // CGAL_IO_WINDOW_STREAM_H

CGAL_BEGIN_NAMESPACE

// Define the stream operators for the xy projected 3d objects.
// Note that data structures like polygons and triangulations
// work independant from the dimension of the stored geometry.
// ------------------------------------------------------------

#ifdef CGAL_POINT_3_H
#ifndef CGAL_WINDOW_STREAM_POINT_3
#define CGAL_WINDOW_STREAM_POINT_3
template< class R >
inline
Window_stream& operator<<( Window_stream &w, const Point_3<R> &p)
{
    return  w << Point_2<R>( p.hx(), p.hy(), p.hw());
}
template< class R >
Window_stream& operator>>( Window_stream &w, Point_3<R> &p)
{
    typedef typename R::RT   RT;
    Point_2<R> q;
    w >> q;
    p =  Point_3<R>( q.hx(), q.hy(), RT(0), q.hw());
    return w;
}
#endif // CGAL_WINDOW_STREAM_POINT_3
#endif // CGAL_POINT_3_H

#ifdef CGAL_VECTOR_3_H
#ifndef CGAL_WINDOW_STREAM_VECTOR_3
#define CGAL_WINDOW_STREAM_VECTOR_3
template< class R >
inline
Window_stream& operator<<( Window_stream &w, const Vector_3<R> &v)
{
    return  w << Vector_2<R>( v.hx(), v.hy(), v.hw());
}
template< class R >
Window_stream& operator>>( Window_stream &w, Vector_3<R> &v)
{
    typedef typename R::RT   RT;
    Vector_2<R> q;
    w >> q;
    v =  Vector_3<R>( q.hx(), q.hy(), RT(0), q.hw());
    return w;
}
#endif // CGAL_WINDOW_STREAM_VECTOR_3
#endif // CGAL_VECTOR_3_H

#ifdef CGAL_DIRECTION_3_H
#ifndef CGAL_WINDOW_STREAM_DIRECTION_3
#define CGAL_WINDOW_STREAM_DIRECTION_3
template< class R >
Window_stream& operator<<( Window_stream &w, const Direction_3<R> &d)
{
    return  w << Direction_2<R>( d.dx(), d.dy());
}
template< class R >
Window_stream& operator>>( Window_stream &w, Direction_3<R> &d)
{
    typedef typename R::RT   RT;
    Direction_2<R> q;
    w >> q;
    d =  Direction_3<R>( q.dx(), q.dy(), RT(0));
    return w;
}
#endif // CGAL_WINDOW_STREAM_DIRECTION_3
#endif // CGAL_DIRECTION_3_H

#ifdef CGAL_LINE_3_H
#ifndef CGAL_WINDOW_STREAM_LINE_3
#define CGAL_WINDOW_STREAM_LINE_3
template< class R >
Window_stream& operator<<( Window_stream &w, const Line_3<R> &l)
{
    return  w << Line_2<R>(
      Point_2<R>( l.point(0).hx(), l.point(0).hy(), l.point(0).hw()),
      Point_2<R>( l.point(1).hx(), l.point(1).hy(), l.point(1).hw()));
}
template< class R >
Window_stream& operator>>( Window_stream &w, Line_3<R> &l)
{
    Line_2<R> q;
    w >> q;
    l =  Line_3<R>(
      Point_3<R>( q.point(0).hx(),q.point(0).hy(),0,q.point(0).hw()),
      Point_3<R>( q.point(1).hx(),q.point(1).hy(),0,q.point(1).hw()));
    return w;
}
#endif // CGAL_WINDOW_STREAM_LINE_3
#endif // CGAL_LINE_3_H

#ifdef CGAL_RAY_3_H
#ifndef CGAL_WINDOW_STREAM_RAY_3
#define CGAL_WINDOW_STREAM_RAY_3
template< class R >
Window_stream& operator<<( Window_stream &w, const Ray_3<R> &r)
{
    return  w << Ray_2<R>(
      Point_2<R>( r.point(0).hx(), r.point(0).hy(), r.point(0).hw()),
      Point_2<R>( r.point(1).hx(), r.point(1).hy(), r.point(1).hw()));
}
template< class R >
Window_stream& operator>>( Window_stream &w, Ray_3<R> &r)
{
    Ray_2<R> q;
    w >> q;
    r =  Ray_3<R>(
      Point_3<R>( q.point(0).hx(),q.point(0).hy(),0,q.point(0).hw()),
      Point_3<R>( q.point(1).hx(),q.point(1).hy(),0,q.point(1).hw()));
    return w;
}
#endif // CGAL_WINDOW_STREAM_RAY_3
#endif // CGAL_RAY_3_H

#ifdef CGAL_SEGMENT_3_H
#ifndef CGAL_WINDOW_STREAM_SEGMENT_3
#define CGAL_WINDOW_STREAM_SEGMENT_3
template< class R >
Window_stream& operator<<( Window_stream &w, const Segment_3<R> &s)
{
    return  w << Segment_2<R>(
      Point_2<R>( s.source().hx(), s.source().hy(), s.source().hw()),
      Point_2<R>( s.target().hx(), s.target().hy(), s.target().hw()));
}
template< class R >
Window_stream& operator>>( Window_stream &w, Segment_3<R> &s)
{
    Segment_2<R> q;
    w >> q;
    s =  Segment_3<R>(
      Point_3<R>( q.source().hx(),q.source().hy(),0,q.source().hw()),
      Point_3<R>( q.target().hx(),q.target().hy(),0,q.target().hw()));
    return w;
}
#endif // CGAL_WINDOW_STREAM_SEGMENT_3
#endif // CGAL_SEGMENT_3_H

#ifdef CGAL_TRIANGLE_3_H
#ifndef CGAL_WINDOW_STREAM_TRIANGLE_3
#define CGAL_WINDOW_STREAM_TRIANGLE_3
template< class R >
Window_stream& operator<<( Window_stream &w, const Triangle_3<R> &t)
{
    return  w << Triangle_2<R>(
        Point_2<R>( t[0].hx(), t[0].hy(), t[0].hw()),
        Point_2<R>( t[1].hx(), t[1].hy(), t[1].hw()),
        Point_2<R>( t[2].hx(), t[2].hy(), t[2].hw()));
}
template< class R >
Window_stream& operator>>( Window_stream &w, Triangle_3<R> &t)
{
    Triangle_2<R> q;
    w >> q;
    t =  Triangle_3<R>(
        Point_3<R>( q[0].hx(), q[0].hy(), 0, q[0].hw()),
        Point_3<R>( q[1].hx(), q[1].hy(), 0, q[1].hw()),
        Point_3<R>( q[2].hx(), q[2].hy(), 0, q[2].hw()));
    return w;
}
#endif // CGAL_WINDOW_STREAM_TRIANGLE_3
#endif // CGAL_TRIANGLE_3_H

#ifdef CGAL_TETRAHEDRON_3_H
#ifndef CGAL_WINDOW_STREAM_TETRAHEDRON_3
#define CGAL_WINDOW_STREAM_TETRAHEDRON_3
template< class R >
Window_stream& operator<<( Window_stream &w, const Tetrahedron_3<R> &t)
{
    w << Segment_3<R>( t[0], t[1]);
    w << Segment_3<R>( t[1], t[2]);
    w << Segment_3<R>( t[2], t[0]);
    w << Segment_3<R>( t[0], t[3]);
    w << Segment_3<R>( t[1], t[3]);
    w << Segment_3<R>( t[2], t[3]);
    return  w;
}
template< class R >
Window_stream& operator>>( Window_stream &w, Tetrahedron_3<R> &t)
{
    typedef typename R::RT   RT;
    Triangle_3<R> q;
    w >> q;
    double x0 = to_double( q[0].x());
    double y0 = to_double( q[0].y());
    double x1, y1;
    w.read_mouse_seg(x0, y0, x1, y1);
    Point_3<R> p = Point_3<R>( RT(x1), RT(y1), RT(0));
    w << Segment_3<R>( q[0], p);
    w << Segment_3<R>( q[1], p);
    w << Segment_3<R>( q[2], p);
    t =  Tetrahedron_3<R>( q[0], q[1], q[2], p);
    return w;
}
#endif // CGAL_WINDOW_STREAM_TETRAHEDRON_3
#endif // CGAL_TETRAHEDRON_3_H

#ifdef CGAL_BBOX_3_H
#ifndef CGAL_WINDOW_STREAM_BBOX_3
#define CGAL_WINDOW_STREAM_BBOX_3
inline
Window_stream& operator<<( Window_stream &w, const Bbox_3 &b)
{
    return  w << Bbox_2( b.xmin(), b.ymin(), b.xmax(), b.ymax());
}
inline
Window_stream& operator>>( Window_stream &w, Bbox_3 &b)
{
    double x0, y0, x1, y1;
    w.read_mouse(x0,y0);
    w.read_mouse_rect(x0,y0, x1, y1);
    if ( x1 < x0) {
        double tmp = x0; x0 = x1; x1 = tmp;
    }
    if ( y1 < y0) {
        double tmp = y0; y0 = y1; y1 = tmp;
    }
    b = Bbox_3( x0, y0, 0, x1, y1, 0);
    w << b;
    return w;
}
#endif // CGAL_WINDOW_STREAM_BBOX_3
#endif // CGAL_BBOX_3_H

CGAL_END_NAMESPACE

// EOF //
