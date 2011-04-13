
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
// file          : include/CGAL/Ray_2_Triangle_2_intersection.h
// source        : intersection_2_2.fw
// author(s)     : Geert-Jan Giezeman
//
// coordinator   : Saarbruecken
//
// ============================================================================


#ifndef CGAL_RAY_2_TRIANGLE_2_INTERSECTION_H
#define CGAL_RAY_2_TRIANGLE_2_INTERSECTION_H

#include <CGAL/Segment_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Point_2.h>

CGAL_BEGIN_NAMESPACE

template <class R>
class Ray_2_Triangle_2_pair {
public:
    enum Intersection_results {NO, POINT, SEGMENT};
    Ray_2_Triangle_2_pair() ;
    Ray_2_Triangle_2_pair(Ray_2<R> const *ray,
                            Triangle_2<R> const *trian);
    ~Ray_2_Triangle_2_pair() {}
#ifdef CGAL_CFG_RETURN_TYPE_BUG_2
    Intersection_results intersection_type() const
    {
        if (_known)
            return _result;
    // The non const this pointer is used to cast away const.
        _known = true;
        Straight_2_<R> straight(*_ray);
    Line_2<R> l(_trian->vertex(0), _trian->vertex(1));
    if (l.oriented_side(_trian->vertex(2)) == ON_POSITIVE_SIDE) {
    //    if (_trian->is_counterclockwise()) {
            straight.cut_right_off(
                Line_2<R>(_trian->vertex(0), _trian->vertex(1)));
            straight.cut_right_off(
                Line_2<R>(_trian->vertex(1), _trian->vertex(2)));
            straight.cut_right_off(
                Line_2<R>(_trian->vertex(2), _trian->vertex(0)));
        } else {
            straight.cut_right_off(
                Line_2<R>(_trian->vertex(2), _trian->vertex(1)));
            straight.cut_right_off(
                Line_2<R>(_trian->vertex(1), _trian->vertex(0)));
            straight.cut_right_off(
                Line_2<R>(_trian->vertex(0), _trian->vertex(2)));
        }
        switch (straight.current_state()) {
        case Straight_2_<R>::EMPTY:
            _result = NO;
            return _result;
        case Straight_2_<R>::POINT: {
            straight.current(_intersection_point);
            _result = POINT;
            return _result;
            }
        case Straight_2_<R>::SEGMENT: {
            Segment_2<R> seg;
            straight.current(seg);
            _intersection_point = seg.start();
            _other_point = seg.end();
            _result = SEGMENT;
            return _result;
            }
        default:  // should not happen.
            CGAL_kernel_assertion_msg(false, "Internal CGAL error.");
            _result = NO;
            return _result;
        }
    }
    
#else
    Intersection_results intersection_type() const;
#endif // CGAL_CFG_RETURN_TYPE_BUG_2
    bool                intersection(Point_2<R> &result) const;
    bool                intersection(Segment_2<R> &result) const;
protected:
    Ray_2<R> const* _ray;
    Triangle_2<R> const *  _trian;
    mutable bool                    _known;
    mutable Intersection_results     _result;
    mutable Point_2<R>         _intersection_point;
    mutable Point_2<R>         _other_point;
};

template <class R>
inline bool do_intersect(
    const Ray_2<R> &p1,
    const Triangle_2<R> &p2)
{
    typedef Ray_2_Triangle_2_pair<R> pair_t;
    pair_t pair(&p1, &p2);
    return pair.intersection_type() != pair_t::NO;
}

CGAL_END_NAMESPACE



#include <CGAL/Line_2.h>
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/Straight_2.h>

CGAL_BEGIN_NAMESPACE

template <class R>
Ray_2_Triangle_2_pair<R>::
Ray_2_Triangle_2_pair()
{
    _known = false;
    _ray = 0;
    _trian = 0;
}

template <class R>
Ray_2_Triangle_2_pair<R>::
Ray_2_Triangle_2_pair(Ray_2<R> const *ray,
                            Triangle_2<R> const *trian)
{
    _known = false;
    _ray = ray;
    _trian = trian;
}

#ifndef CGAL_CFG_RETURN_TYPE_BUG_2
template <class R>
Ray_2_Triangle_2_pair<R>::Intersection_results
Ray_2_Triangle_2_pair<R>::intersection_type() const
{
    if (_known)
        return _result;
// The non const this pointer is used to cast away const.
    _known = true;
    Straight_2_<R> straight(*_ray);
Line_2<R> l(_trian->vertex(0), _trian->vertex(1));
if (l.oriented_side(_trian->vertex(2)) == ON_POSITIVE_SIDE) {
//    if (_trian->is_counterclockwise()) {
        straight.cut_right_off(
            Line_2<R>(_trian->vertex(0), _trian->vertex(1)));
        straight.cut_right_off(
            Line_2<R>(_trian->vertex(1), _trian->vertex(2)));
        straight.cut_right_off(
            Line_2<R>(_trian->vertex(2), _trian->vertex(0)));
    } else {
        straight.cut_right_off(
            Line_2<R>(_trian->vertex(2), _trian->vertex(1)));
        straight.cut_right_off(
            Line_2<R>(_trian->vertex(1), _trian->vertex(0)));
        straight.cut_right_off(
            Line_2<R>(_trian->vertex(0), _trian->vertex(2)));
    }
    switch (straight.current_state()) {
    case Straight_2_<R>::EMPTY:
        _result = NO;
        return _result;
    case Straight_2_<R>::POINT: {
        straight.current(_intersection_point);
        _result = POINT;
        return _result;
        }
    case Straight_2_<R>::SEGMENT: {
        Segment_2<R> seg;
        straight.current(seg);
        _intersection_point = seg.start();
        _other_point = seg.end();
        _result = SEGMENT;
        return _result;
        }
    default:  // should not happen.
        CGAL_kernel_assertion_msg(false, "Internal CGAL error.");
        _result = NO;
        return _result;
    }
}

#endif // CGAL_CFG_RETURN_TYPE_BUG_2

template <class R>
bool
Ray_2_Triangle_2_pair<R>::
intersection(Point_2<R> &result) const
{
    if (!_known)
        intersection_type();
    if (_result != POINT)
        return false;
    result = _intersection_point;
    return true;
}

template <class R>
bool
Ray_2_Triangle_2_pair<R>::
intersection(Segment_2<R> &result) const
{
    if (!_known)
        intersection_type();
    if (_result != SEGMENT)
        return false;
    result = Segment_2<R>(_intersection_point, _other_point);
    return true;
}

CGAL_END_NAMESPACE



#include <CGAL/Object.h>

CGAL_BEGIN_NAMESPACE

template <class R>
Object
intersection(const Ray_2<R> &ray, const Triangle_2<R>&tr)
{
    typedef Ray_2_Triangle_2_pair<R> is_t;
    is_t ispair(&ray, &tr);
    switch (ispair.intersection_type()) {
    case is_t::NO:
    default:
        return Object();
    case is_t::POINT: {
        Point_2<R> pt;
        ispair.intersection(pt);
        return Object(new Wrapper< Point_2<R> >(pt));
    }
    case is_t::SEGMENT: {
        Segment_2<R> iseg;
        ispair.intersection(iseg);
        return Object(new Wrapper< Segment_2<R> >(iseg));
    }
    }
}

template <class R>
class Triangle_2_Ray_2_pair
: public Ray_2_Triangle_2_pair<R> {
public:
    Triangle_2_Ray_2_pair(
            Triangle_2<R> const *trian,
            Ray_2<R> const *ray) :
                        Ray_2_Triangle_2_pair<R>(ray, trian) {}
};

template <class R>
inline bool do_intersect(
    const Triangle_2<R> &p1,
    const Ray_2<R> &p2)
{
    typedef Triangle_2_Ray_2_pair<R> pair_t;
    pair_t pair(&p1, &p2);
    return pair.intersection_type() != pair_t::NO;
}

template <class R>
inline Object
intersection(const Triangle_2<R> &tr, const Ray_2<R> &ray)
{
    return intersection(ray, tr);
}

CGAL_END_NAMESPACE

#endif
