
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
// file          : include/CGAL/Ray_2_Segment_2_intersection.h
// source        : intersection_2_1.fw
// author(s)     : Geert-Jan Giezeman
//
// coordinator   : Saarbruecken
//
// ============================================================================


#ifndef CGAL_RAY_2_SEGMENT_2_INTERSECTION_H
#define CGAL_RAY_2_SEGMENT_2_INTERSECTION_H

#include <CGAL/Ray_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>

CGAL_BEGIN_NAMESPACE

template <class R>
class Ray_2_Segment_2_pair {
public:
    enum Intersection_results {NO, POINT, SEGMENT};
    Ray_2_Segment_2_pair() ;
    Ray_2_Segment_2_pair(Ray_2<R> const *ray,
                            Segment_2<R> const *line);
    ~Ray_2_Segment_2_pair() {}

#ifndef CGAL_CFG_RETURN_TYPE_BUG_2
    Intersection_results intersection_type() const;
#else
    Intersection_results intersection_type() const
{
    if (_known)
        return _result;
    // The non const this pointer is used to cast away const.
    _known = true;
//    if (!do_overlap(_ray->bbox(), _seg->bbox()))
//        return NO;
    const Line_2<R> &l1 = _ray->supporting_line();
    const Line_2<R> &l2 = _seg->supporting_line();
    Line_2_Line_2_pair<R> linepair(&l1, &l2);
    switch ( linepair.intersection_type()) {
    case Line_2_Line_2_pair<R>::NO:
        _result = NO;
        return _result;
    case Line_2_Line_2_pair<R>::POINT:
        linepair.intersection(_intersection_point);
        _result = (_ray->collinear_has_on(_intersection_point)
                && _seg->collinear_has_on(_intersection_point) )
            ? POINT :  NO;
        return _result;
    case Line_2_Line_2_pair<R>::LINE: {
        typedef typename R::RT RT;
        const Point_2<R> &start1 = _seg->start();
        const Point_2<R> &end1 = _seg->end();
        const Point_2<R> &start2 = _ray->start();
        const Point_2<R> *minpt,  *maxpt;
        Vector_2<R> diff1 = end1-start1;
        if (CGAL_NTS abs(diff1.x()) > CGAL_NTS abs(diff1.y())) {
            typedef typename R::FT FT;
            if (start1.x() < end1.x()) {
                minpt = &start1;
                maxpt = &end1;
            } else {
                minpt = &end1;
                maxpt = &start1;
            }
            if (_ray->direction().to_vector().x() > FT(0)) {
                if (maxpt->x() < start2.x()) {
                    _result = NO;
                    return _result;
                }
                if (maxpt->x() == start2.x()) {
                    _intersection_point = *maxpt;
                    _result = POINT;
                    return _result;
                }
                if (minpt->x() < start2.x()) {
                    _intersection_point = start2;
                    _other_point = *maxpt;
                } else {
                    _intersection_point = _seg->start();
                    _other_point = _seg->end();
                }
                _result = SEGMENT;
                return _result;
            } else {
                if (minpt->x() > start2.x()) {
                    _result = NO;
                    return _result;
                }
                if (minpt->x() == start2.x()) {
                    _intersection_point = *minpt;
                    _result = POINT;
                    return _result;
                }
                if (maxpt->x() > start2.x()) {
                    _intersection_point = start2;
                    _other_point = *maxpt;
                } else {
                    _intersection_point = _seg->start();
                    _other_point = _seg->end();
                }
                _result = SEGMENT;
                return _result;
            }
        } else {
            typedef typename R::FT FT;
            if (start1.y() < end1.y()) {
                minpt = &start1;
                maxpt = &end1;
            } else {
                minpt = &end1;
                maxpt = &start1;
            }
            if (_ray->direction().to_vector().y() > FT(0)) {
                if (maxpt->y() < start2.y()) {
                    _result = NO;
                    return _result;
                }
                if (maxpt->y() == start2.y()) {
                    _intersection_point = *maxpt;
                    _result = POINT;
                    return _result;
                }
                if (minpt->y() < start2.y()) {
                    _intersection_point = start2;
                    _other_point = *maxpt;
                } else {
                    _intersection_point = _seg->start();
                    _other_point = _seg->end();
                }
                _result = SEGMENT;
                return _result;
            } else {
                if (minpt->y() > start2.y()) {
                    _result = NO;
                    return _result;
                }
                if (minpt->y() == start2.y()) {
                    _intersection_point = *minpt;
                    _result = POINT;
                    return _result;
                }
                if (maxpt->y() > start2.y()) {
                    _intersection_point = start2;
                    _other_point = *maxpt;
                } else {
                    _intersection_point = _seg->start();
                    _other_point = _seg->end();
                }
                _result = SEGMENT;
                return _result;
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
protected:
    Ray_2<R> const *   _ray;
    Segment_2<R> const *  _seg;
    mutable bool                    _known;
    mutable Intersection_results    _result;
    mutable Point_2<R>         _intersection_point, _other_point;
};

template <class R>
inline bool do_intersect(
    const Ray_2<R> &p1,
    const Segment_2<R> &p2)
{
    typedef Ray_2_Segment_2_pair<R> pair_t;
    pair_t pair(&p1, &p2);
    return pair.intersection_type() != pair_t::NO;
}

CGAL_END_NAMESPACE


#include <CGAL/Line_2.h>
#include <CGAL/Line_2_Line_2_intersection.h>

CGAL_BEGIN_NAMESPACE

template <class R>
Ray_2_Segment_2_pair<R>::Ray_2_Segment_2_pair()
{
    _ray = 0;
    _seg = 0;
    _known = false;
}

template <class R>
Ray_2_Segment_2_pair<R>::Ray_2_Segment_2_pair(
    Ray_2<R> const *ray, Segment_2<R> const *seg)
{
    _ray = ray;
    _seg = seg;
    _known = false;
}

#ifndef CGAL_CFG_RETURN_TYPE_BUG_2
template <class R>
Ray_2_Segment_2_pair<R>::Intersection_results
Ray_2_Segment_2_pair<R>::intersection_type() const
{
    if (_known)
        return _result;
    // The non const this pointer is used to cast away const.
    _known = true;
//    if (!do_overlap(_ray->bbox(), _seg->bbox()))
//        return NO;
    const Line_2<R> &l1 = _ray->supporting_line();
    const Line_2<R> &l2 = _seg->supporting_line();
    Line_2_Line_2_pair<R> linepair(&l1, &l2);
    switch ( linepair.intersection_type()) {
    case Line_2_Line_2_pair<R>::NO:
        _result = NO;
        return _result;
    case Line_2_Line_2_pair<R>::POINT:
        linepair.intersection(_intersection_point);
        _result = (_ray->collinear_has_on(_intersection_point)
                && _seg->collinear_has_on(_intersection_point) )
            ? POINT :  NO;
        return _result;
    case Line_2_Line_2_pair<R>::LINE: {
        typedef typename R::RT RT;
        const Point_2<R> &start1 = _seg->start();
        const Point_2<R> &end1 = _seg->end();
        const Point_2<R> &start2 = _ray->start();
        const Point_2<R> *minpt,  *maxpt;
        Vector_2<R> diff1 = end1-start1;
        if (CGAL_NTS abs(diff1.x()) > CGAL_NTS abs(diff1.y())) {
            typedef typename R::FT FT;
            if (start1.x() < end1.x()) {
                minpt = &start1;
                maxpt = &end1;
            } else {
                minpt = &end1;
                maxpt = &start1;
            }
            if (_ray->direction().to_vector().x() > FT(0)) {
                if (maxpt->x() < start2.x()) {
                    _result = NO;
                    return _result;
                }
                if (maxpt->x() == start2.x()) {
                    _intersection_point = *maxpt;
                    _result = POINT;
                    return _result;
                }
                if (minpt->x() < start2.x()) {
                    _intersection_point = start2;
                    _other_point = *maxpt;
                } else {
                    _intersection_point = _seg->start();
                    _other_point = _seg->end();
                }
                _result = SEGMENT;
                return _result;
            } else {
                if (minpt->x() > start2.x()) {
                    _result = NO;
                    return _result;
                }
                if (minpt->x() == start2.x()) {
                    _intersection_point = *minpt;
                    _result = POINT;
                    return _result;
                }
                if (maxpt->x() > start2.x()) {
                    _intersection_point = start2;
                    _other_point = *maxpt;
                } else {
                    _intersection_point = _seg->start();
                    _other_point = _seg->end();
                }
                _result = SEGMENT;
                return _result;
            }
        } else {
            typedef typename R::FT FT;
            if (start1.y() < end1.y()) {
                minpt = &start1;
                maxpt = &end1;
            } else {
                minpt = &end1;
                maxpt = &start1;
            }
            if (_ray->direction().to_vector().y() > FT(0)) {
                if (maxpt->y() < start2.y()) {
                    _result = NO;
                    return _result;
                }
                if (maxpt->y() == start2.y()) {
                    _intersection_point = *maxpt;
                    _result = POINT;
                    return _result;
                }
                if (minpt->y() < start2.y()) {
                    _intersection_point = start2;
                    _other_point = *maxpt;
                } else {
                    _intersection_point = _seg->start();
                    _other_point = _seg->end();
                }
                _result = SEGMENT;
                return _result;
            } else {
                if (minpt->y() > start2.y()) {
                    _result = NO;
                    return _result;
                }
                if (minpt->y() == start2.y()) {
                    _intersection_point = *minpt;
                    _result = POINT;
                    return _result;
                }
                if (maxpt->y() > start2.y()) {
                    _intersection_point = start2;
                    _other_point = *maxpt;
                } else {
                    _intersection_point = _seg->start();
                    _other_point = _seg->end();
                }
                _result = SEGMENT;
                return _result;
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
Ray_2_Segment_2_pair<R>::intersection(Point_2<R> &result) const
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
Ray_2_Segment_2_pair<R>::intersection(Segment_2<R> &result) const
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
intersection(const Ray_2<R> &ray, const Segment_2<R>&seg)
{
    typedef Ray_2_Segment_2_pair<R> is_t;
    is_t ispair(&ray, &seg);
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
class Segment_2_Ray_2_pair: public Ray_2_Segment_2_pair<R> {
public:
    Segment_2_Ray_2_pair(
            Segment_2<R> const *seg,
            Ray_2<R> const *ray) :
                                Ray_2_Segment_2_pair<R>(ray, seg) {}
};

template <class R>
inline bool do_intersect(
    const Segment_2<R> &p1,
    const Ray_2<R> &p2)
{
    typedef Segment_2_Ray_2_pair<R> pair_t;
    pair_t pair(&p1, &p2);
    return pair.intersection_type() != pair_t::NO;
}

template <class R>
inline Object
intersection(const Segment_2<R> &seg, const Ray_2<R> &ray)
{
    return intersection(ray, seg);
}

CGAL_END_NAMESPACE

#endif
