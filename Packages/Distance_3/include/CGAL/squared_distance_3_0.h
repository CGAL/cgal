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
#include <CGAL/wmult.h>

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

#if defined CGAL_HOMOGENEOUS_H
template <class RT>
Vector_3< Homogeneous<RT> >
wcross(const Point_3< Homogeneous<RT> > &p,
    const Point_3< Homogeneous<RT> > &q,
    const Point_3< Homogeneous<RT> > &r)
{
    RT x,y,z;
    x =  p.hy() * (q.hz()*r.hw() - q.hw()*r.hz() )
       + p.hz() * (q.hw()*r.hy() - q.hy()*r.hw() )
       + p.hw() * (q.hy()*r.hz() - q.hz()*r.hy() );
    y =  p.hz() * (q.hx()*r.hw() - q.hw()*r.hx() )
       + p.hx() * (q.hw()*r.hz() - q.hz()*r.hw() )
       + p.hw() * (q.hz()*r.hx() - q.hx()*r.hz() );
    z =  p.hx() * (q.hy()*r.hw() - q.hw()*r.hy() )
       + p.hy() * (q.hw()*r.hx() - q.hx()*r.hw() )
       + p.hw() * (q.hx()*r.hy() - q.hy()*r.hx() );
    return Vector_3< Homogeneous<RT> >(x, y, z);
}
#endif // CGAL_HOMOGENEOUS_H

#if defined CGAL_SIMPLE_HOMOGENEOUS_H
template <class RT>
Vector_3< Simple_homogeneous<RT> >
wcross(const Point_3< Simple_homogeneous<RT> > &p,
    const Point_3< Simple_homogeneous<RT> > &q,
    const Point_3< Simple_homogeneous<RT> > &r)
{
    RT x,y,z;
    x =  p.hy() * (q.hz()*r.hw() - q.hw()*r.hz() )
       + p.hz() * (q.hw()*r.hy() - q.hy()*r.hw() )
       + p.hw() * (q.hy()*r.hz() - q.hz()*r.hy() );
    y =  p.hz() * (q.hx()*r.hw() - q.hw()*r.hx() )
       + p.hx() * (q.hw()*r.hz() - q.hz()*r.hw() )
       + p.hw() * (q.hz()*r.hx() - q.hx()*r.hz() );
    z =  p.hx() * (q.hy()*r.hw() - q.hw()*r.hy() )
       + p.hy() * (q.hw()*r.hx() - q.hx()*r.hw() )
       + p.hw() * (q.hx()*r.hy() - q.hy()*r.hx() );
    return Vector_3< Simple_homogeneous<RT> >(x, y, z);
}
#endif // CGAL_SIMPLE_HOMOGENEOUS_H

#if defined CGAL_CARTESIAN_H
template <class FT>
Vector_3< Cartesian<FT> >
wcross(const Point_3< Cartesian<FT> > &p,
    const Point_3< Cartesian<FT> > &q,
    const Point_3< Cartesian<FT> > &r)
{
    FT x,y,z;
    x = (q.y()-p.y())*(r.z()-q.z()) - (q.z()-p.z())*(r.y()-q.y());
    y = (q.z()-p.z())*(r.x()-q.x()) - (q.x()-p.x())*(r.z()-q.z());
    z = (q.x()-p.x())*(r.y()-q.y()) - (q.y()-p.y())*(r.x()-q.x());
    return Vector_3< Cartesian<FT> >(x, y, z);
}
#endif // CGAL_CARTESIAN_H

#if defined CGAL_SIMPLE_CARTESIAN_H
template <class FT>
Vector_3< Simple_cartesian<FT> >
wcross(const Point_3< Simple_cartesian<FT> > &p,
    const Point_3< Simple_cartesian<FT> > &q,
    const Point_3< Simple_cartesian<FT> > &r)
{
    FT x,y,z;
    x = (q.y()-p.y())*(r.z()-q.z()) - (q.z()-p.z())*(r.y()-q.y());
    y = (q.z()-p.z())*(r.x()-q.x()) - (q.x()-p.x())*(r.z()-q.z());
    z = (q.x()-p.x())*(r.y()-q.y()) - (q.y()-p.y())*(r.x()-q.x());
    return Vector_3< Simple_cartesian<FT> >(x, y, z);
}
#endif // CGAL_SIMPLE_CARTESIAN_H



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
//    return FT(dot*dot) /
//        FT(wmult((R*)0, squared_length, diff.hw(), diff.hw()));
    return R::make_FT((dot*dot),
        wmult((R*)0, squared_length, diff.hw(), diff.hw()));
}


template <class R>
typename R::FT
squared_distance_to_line(
    const Vector_3<R> & dir,
    const Vector_3<R> & diff)
{
    typedef typename R::RT RT;
    Vector_3<R> wcr = wcross(dir, diff);
    return R::make_FT((wcr*wcr),
        wmult((R*)0, RT(wdot(dir, dir)), diff.hw(), diff.hw()));
}



CGAL_END_NAMESPACE


#endif
