
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
// file          : include/CGAL/Ray_2_Ray_2_intersection.h
// source        : intersection_2_1.fw
// author(s)     : Geert-Jan Giezeman
//
// coordinator   : Saarbruecken
//
// ============================================================================


#ifndef CGAL_RAY_2_RAY_2_INTERSECTION_H
#define CGAL_RAY_2_RAY_2_INTERSECTION_H

#include <CGAL/Ray_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>

CGAL_BEGIN_NAMESPACE

template <class R>
class Ray_2_Ray_2_pair {
public:
    enum Intersection_results {NO, POINT, SEGMENT, RAY};
    Ray_2_Ray_2_pair() ;
    Ray_2_Ray_2_pair(Ray_2<R> const *ray1,
                            Ray_2<R> const *ray2);
    ~Ray_2_Ray_2_pair() {}

#ifndef CGAL_CFG_RETURN_TYPE_BUG_2
    Intersection_results intersection_type() const;
#else
    Intersection_results intersection_type() const
{
    if (_known)
        return _result;
    // The non const this pointer is used to cast away const.
    _known = true;
//    if (!do_overlap(_ray1->bbox(), _ray2->bbox()))
//        return NO;
    const Line_2<R> &l1 = _ray1->supporting_line();
    const Line_2<R> &l2 = _ray2->supporting_line();
    Line_2_Line_2_pair<R> linepair(&l1, &l2);
    switch ( linepair.intersection_type()) {
    case Line_2_Line_2_pair<R>::NO:
        _result = NO;
        return _result;
    case Line_2_Line_2_pair<R>::POINT:
        linepair.intersection(_intersection_point);
        _result = (_ray1->collinear_has_on(_intersection_point)
                && _ray2->collinear_has_on(_intersection_point) )
            ? POINT :  NO;
        return _result;
    case Line_2_Line_2_pair<R>::LINE:
        {
        typedef typename R::RT RT;
        const Vector_2<R> &dir1 = _ray1->direction().to_vector();
        const Vector_2<R> &dir2 = _ray2->direction().to_vector();
        if (CGAL_NTS abs(dir1.x()) > CGAL_NTS abs(dir1.y())) {
            typedef typename R::FT FT;
            if (dir1.x() > FT(0)) {
                if (dir2.x() > FT(0)) {
                    _intersection_point =
                            (_ray1->start().x() < _ray2->start().x())
                            ? _ray2->start() : _ray1->start();
                    _result = RAY;
                    return _result;
                } else {
                    if (_ray1->start().x() > _ray2->start().x()) {
                        _result = NO;
                        return _result;
                    }
                    if (_ray1->start().x() == _ray2->start().x()) {
                        _intersection_point = _ray1->start();
                        _result = POINT;
                        return _result;
                    }
                    _result = SEGMENT;
                    return _result;
                }
            } else {
                if (dir2.x() < FT(0)) {
                    _intersection_point =
                            (_ray1->start().x() > _ray2->start().x())
                            ? _ray2->start() : _ray1->start();
                    _result = RAY;
                    return _result;
                } else {
                    if (_ray1->start().x() < _ray2->start().x()) {
                        _result = NO;
                        return _result;
                    }
                    if (_ray1->start().x() == _ray2->start().x()) {
                        _intersection_point = _ray1->start();
                        _result = POINT;
                        return _result;
                    }
                    _result = SEGMENT;
                    return _result;
                }
            }
            
        } else {
            typedef typename R::FT FT;
            if (dir1.y() > FT(0)) {
                if (dir2.y() > FT(0)) {
                    _intersection_point =
                            (_ray1->start().y() < _ray2->start().y())
                            ? _ray2->start() : _ray1->start();
                    _result = RAY;
                    return _result;
                } else {
                    if (_ray1->start().y() > _ray2->start().y()) {
                        _result = NO;
                        return _result;
                    }
                    if (_ray1->start().y() == _ray2->start().y()) {
                        _intersection_point = _ray1->start();
                        _result = POINT;
                        return _result;
                    }
                    _result = SEGMENT;
                    return _result;
                }
            } else {
                if (dir2.y() < FT(0)) {
                    _intersection_point =
                            (_ray1->start().y() > _ray2->start().y())
                            ? _ray2->start() : _ray1->start();
                    _result = RAY;
                    return _result;
                } else {
                    if (_ray1->start().y() < _ray2->start().y()) {
                        _result = NO;
                        return _result;
                    }
                    if (_ray1->start().y() == _ray2->start().y()) {
                        _intersection_point = _ray1->start();
                        _result = POINT;
                        return _result;
                    }
                    _result = SEGMENT;
                    return _result;
                }
            }
            
        }
        } 
    default:
        CGAL_kernel_assertion(false); // should not be reached:
        return _result;
    }
}

#endif // CGAL_CFG_RETURN_TYPE_BUG_2

    bool                intersection(Point_2<R> &result) const;
    bool                intersection(Segment_2<R> &result) const;
    bool                intersection(Ray_2<R> &result) const;
protected:
    Ray_2<R> const*    _ray1;
    Ray_2<R> const *   _ray2;
    mutable bool                    _known;
    mutable Intersection_results    _result;
    mutable Point_2<R>         _intersection_point, _other_point;
};

template <class R>
inline bool do_intersect(
    const Ray_2<R> &p1,
    const Ray_2<R> &p2)
{
    typedef Ray_2_Ray_2_pair<R> pair_t;
    pair_t pair(&p1, &p2);
    return pair.intersection_type() != pair_t::NO;
}

CGAL_END_NAMESPACE



#include <CGAL/Line_2.h>
#include <CGAL/Line_2_Line_2_intersection.h>

CGAL_BEGIN_NAMESPACE

template <class R>
Ray_2_Ray_2_pair<R>::Ray_2_Ray_2_pair()
{
    _ray1 = 0;
    _ray2 = 0;
    _known = false;
}

template <class R>
Ray_2_Ray_2_pair<R>::Ray_2_Ray_2_pair(
    Ray_2<R> const *ray1, Ray_2<R> const *ray2)
{
    _ray1 = ray1;
    _ray2 = ray2;
    _known = false;
}

#ifndef CGAL_CFG_RETURN_TYPE_BUG_2
template <class R>
Ray_2_Ray_2_pair<R>::Intersection_results
Ray_2_Ray_2_pair<R>::intersection_type() const
{
    if (_known)
        return _result;
    // The non const this pointer is used to cast away const.
    _known = true;
//    if (!do_overlap(_ray1->bbox(), _ray2->bbox()))
//        return NO;
    const Line_2<R> &l1 = _ray1->supporting_line();
    const Line_2<R> &l2 = _ray2->supporting_line();
    Line_2_Line_2_pair<R> linepair(&l1, &l2);
    switch ( linepair.intersection_type()) {
    case Line_2_Line_2_pair<R>::NO:
        _result = NO;
        return _result;
    case Line_2_Line_2_pair<R>::POINT:
        linepair.intersection(_intersection_point);
        _result = (_ray1->collinear_has_on(_intersection_point)
                && _ray2->collinear_has_on(_intersection_point) )
            ? POINT :  NO;
        return _result;
    case Line_2_Line_2_pair<R>::LINE:
        {
        typedef typename R::RT RT;
        const Vector_2<R> &dir1 = _ray1->direction().to_vector();
        const Vector_2<R> &dir2 = _ray2->direction().to_vector();
        if (CGAL_NTS abs(dir1.x()) > CGAL_NTS abs(dir1.y())) {
            typedef typename R::FT FT;
            if (dir1.x() > FT(0)) {
                if (dir2.x() > FT(0)) {
                    _intersection_point =
                            (_ray1->start().x() < _ray2->start().x())
                            ? _ray2->start() : _ray1->start();
                    _result = RAY;
                    return _result;
                } else {
                    if (_ray1->start().x() > _ray2->start().x()) {
                        _result = NO;
                        return _result;
                    }
                    if (_ray1->start().x() == _ray2->start().x()) {
                        _intersection_point = _ray1->start();
                        _result = POINT;
                        return _result;
                    }
                    _result = SEGMENT;
                    return _result;
                }
            } else {
                if (dir2.x() < FT(0)) {
                    _intersection_point =
                            (_ray1->start().x() > _ray2->start().x())
                            ? _ray2->start() : _ray1->start();
                    _result = RAY;
                    return _result;
                } else {
                    if (_ray1->start().x() < _ray2->start().x()) {
                        _result = NO;
                        return _result;
                    }
                    if (_ray1->start().x() == _ray2->start().x()) {
                        _intersection_point = _ray1->start();
                        _result = POINT;
                        return _result;
                    }
                    _result = SEGMENT;
                    return _result;
                }
            }
            
        } else {
            typedef typename R::FT FT;
            if (dir1.y() > FT(0)) {
                if (dir2.y() > FT(0)) {
                    _intersection_point =
                            (_ray1->start().y() < _ray2->start().y())
                            ? _ray2->start() : _ray1->start();
                    _result = RAY;
                    return _result;
                } else {
                    if (_ray1->start().y() > _ray2->start().y()) {
                        _result = NO;
                        return _result;
                    }
                    if (_ray1->start().y() == _ray2->start().y()) {
                        _intersection_point = _ray1->start();
                        _result = POINT;
                        return _result;
                    }
                    _result = SEGMENT;
                    return _result;
                }
            } else {
                if (dir2.y() < FT(0)) {
                    _intersection_point =
                            (_ray1->start().y() > _ray2->start().y())
                            ? _ray2->start() : _ray1->start();
                    _result = RAY;
                    return _result;
                } else {
                    if (_ray1->start().y() < _ray2->start().y()) {
                        _result = NO;
                        return _result;
                    }
                    if (_ray1->start().y() == _ray2->start().y()) {
                        _intersection_point = _ray1->start();
                        _result = POINT;
                        return _result;
                    }
                    _result = SEGMENT;
                    return _result;
                }
            }
            
        }
        } 
    default:
        CGAL_kernel_assertion(false); // should not be reached:
        return _result;
    }
}

#endif // CGAL_CFG_RETURN_TYPE_BUG_2

template <class R>
bool
Ray_2_Ray_2_pair<R>::intersection(Point_2<R> &result) const
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
Ray_2_Ray_2_pair<R>::intersection(Segment_2<R> &result) const
{
    if (!_known)
        intersection_type();
    if (_result != RAY)
        return false;
    result = Segment_2<R>(_ray1->start(), _ray2->start());
    return true;
}

template <class R>
bool
Ray_2_Ray_2_pair<R>::intersection(Ray_2<R> &result) const
{
    if (!_known)
        intersection_type();
    if (_result != RAY)
        return false;
    result = Ray_2<R>(_intersection_point, _ray1->direction());
    return true;
}

CGAL_END_NAMESPACE



#include <CGAL/Object.h>

CGAL_BEGIN_NAMESPACE

template <class R>
Object
intersection(const Ray_2<R> &ray1, const Ray_2<R>&ray2)
{
    typedef Ray_2_Ray_2_pair<R> is_t;
    is_t ispair(&ray1, &ray2);
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
    case is_t::RAY: {
        Segment_2<R> iray;
        ispair.intersection(iray);
        return Object(new Wrapper< Segment_2<R> >(iray));
    }
    }
}

CGAL_END_NAMESPACE

#endif
