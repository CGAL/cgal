// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
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
// file          : include/CGAL/squared_distance_3_0.h
// source        : sqdistance_3.fw
// author(s)     : Geert-Jan Giezeman
//
// coordinator   : Saarbruecken
//
// ============================================================================


#ifndef CGAL_DISTANCE_3_0_H
#define CGAL_DISTANCE_3_0_H

#include <CGAL/Point_3.h>

#include <CGAL/utils.h>
#include <CGAL/enum.h>
#include <CGAL/Kernel/Wutils.h>

CGAL_BEGIN_NAMESPACE



template <class R>
bool is_null(const Vector_3<R> &v)
{
    typedef typename R::RT RT;
    return v.hx()==RT(0) && v.hy()==RT(0) && v.hz()==RT(0);
}


template <class R>
typename R::RT
wdot(const Vector_3<R> &u,
    const Vector_3<R> &v)
{
    return  (u.hx()*v.hx() + u.hy()*v.hy() + u.hz()*v.hz());
}


template <class R>
typename R::RT
wdot(const Point_3<R> &p,
     const Point_3<R> &q,
     const Point_3<R> &r)
{
    R* pR = 0;
    return  (wmult(pR, p.hx(),q.hw()) - wmult(pR, q.hx(),p.hw()))
          * (wmult(pR, r.hx(),q.hw()) - wmult(pR, q.hx(),r.hw()))
        +   (wmult(pR, p.hy(),q.hw()) - wmult(pR, q.hy(),p.hw()))
          * (wmult(pR, r.hy(),q.hw()) - wmult(pR, q.hy(),r.hw()))
        +   (wmult(pR, p.hz(),q.hw()) - wmult(pR, q.hz(),p.hw()))
          * (wmult(pR, r.hz(),q.hw()) - wmult(pR, q.hz(),r.hw()));
}




template <class R>
Vector_3<R> wcross(const Vector_3<R> &u,
    const Vector_3<R> &v)
{
    return Vector_3<R>(
        u.hy()*v.hz() - u.hz()*v.hy(),
        u.hz()*v.hx() - u.hx()*v.hz(),
        u.hx()*v.hy() - u.hy()*v.hx());
}


template <class R>
inline bool is_acute_angle(const Vector_3<R> &u,
    const Vector_3<R> &v)
{
    typedef typename R::RT RT;
    return RT(wdot(u, v)) > RT(0) ;
}

template <class R>
inline bool is_straight_angle(const Vector_3<R> &u,
    const Vector_3<R> &v)
{
    typedef typename R::RT RT;
    return RT(wdot(u, v)) == RT(0) ;
}

template <class R>
inline bool is_obtuse_angle(const Vector_3<R> &u,
    const Vector_3<R> &v)
{
    typedef typename R::RT RT;
    return RT(wdot(u, v)) < RT(0) ;
}

template <class R>
inline bool is_acute_angle(const Point_3<R> &p,
    const Point_3<R> &q, const Point_3<R> &r)
{
    typedef typename R::RT RT;
    return RT(wdot(p, q, r)) > RT(0) ;
}

template <class R>
inline bool is_straight_angle(const Point_3<R> &p,
    const Point_3<R> &q, const Point_3<R> &r)
{
    typedef typename R::RT RT;
    return RT(wdot(p, q, r)) == RT(0) ;
}

template <class R>
inline bool is_obtuse_angle(const Point_3<R> &p,
    const Point_3<R> &q, const Point_3<R> &r)
{
    typedef typename R::RT RT;
    return RT(wdot(p, q, r)) < RT(0) ;
}


template <class R>
inline typename R::FT
squared_distance(
    const Point_3<R> & pt1,
    const Point_3<R> & pt2)
{
    Vector_3<R> vec(pt1-pt2);
    return vec*vec;
}


template <class R>
typename R::FT
squared_distance_to_plane(
    const Vector_3<R> & normal,
    const Vector_3<R> & diff)
{
    typedef typename R::RT RT;
    typedef typename R::FT FT;
    RT dot, squared_length;
    dot = wdot(normal, diff);
    squared_length = wdot(normal, normal);
    return FT(dot*dot) /
        FT(wmult((R*)0, squared_length, diff.hw(), diff.hw()));
//    return R::make_FT((dot*dot),
//        wmult((R*)0, squared_length, diff.hw(), diff.hw()));
}


template <class R>
typename R::FT
squared_distance_to_line(
    const Vector_3<R> & dir,
    const Vector_3<R> & diff)
{
    typedef typename R::RT RT;
    typedef typename R::FT FT;
    Vector_3<R> wcr = wcross(dir, diff);
    return FT(wcr*wcr)/FT(wmult(
        (R*)0, RT(wdot(dir, dir)), diff.hw(), diff.hw()));
//    return R::make_FT((wcr*wcr),
//        wmult((R*)0, RT(wdot(dir, dir)), diff.hw(), diff.hw()));
}



CGAL_END_NAMESPACE


#endif
