
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
// file          : include/CGAL/Point_2_Triangle_2_intersection.h
// source        : intersection_2_2.fw
// author(s)     : Geert-Jan Giezeman
//
// coordinator   : Saarbruecken
//
// ============================================================================


#ifndef CGAL_POINT_2_TRIANGLE_2_INTERSECTION_H
#define CGAL_POINT_2_TRIANGLE_2_INTERSECTION_H

#include <CGAL/Point_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Point_2.h>

CGAL_BEGIN_NAMESPACE

template <class R>
class Point_2_Triangle_2_pair {
public:
    enum Intersection_results {NO, POINT};
    Point_2_Triangle_2_pair() ;
    Point_2_Triangle_2_pair(Point_2<R> const *pt,
                            Triangle_2<R> const *trian);
    ~Point_2_Triangle_2_pair() {}
#ifdef CGAL_CFG_RETURN_TYPE_BUG_2
    Intersection_results intersection_type() const
    {
        typedef Line_2<R> line_t;
        if (_known)
            return _result;
    // The non const this pointer is used to cast away const.
        _known = true;
        if (_trian->has_on_unbounded_side(*_pt)) {
            _result = NO;
        } else {
            _result = POINT;
        }
        return _result;
    /*
        line_t l(_trian->vertex(0), _trian->vertex(1));
        if (l.has_on_positive_side(_trian->vertex(2))) {
            for (int i=0; i<3; i++) {
                if (line_t(_trian->vertex(i), _trian->vertex(i+1)).
                                    has_on_negative_side(*_pt)) {
                    _result = NO;
                    return _result;
                }
            }
        } else {
            for (int i=0; i<3; i++)
                if(line_t(_trian->vertex(i), _trian->vertex(i-1)).
                                    has_on_negative_side(*_pt)){
                    _result = NO;
                    return _result;
                }
        }
    */
    }
    
#else
    Intersection_results intersection_type() const;
#endif // CGAL_CFG_RETURN_TYPE_BUG_2
    bool                intersection(Point_2<R> &result) const;
protected:
    Point_2<R> const *    _pt;
    Triangle_2<R> const * _trian;
    mutable bool                       _known;
    mutable Intersection_results       _result;
    mutable Point_2<R>            _intersection_point;
    mutable Point_2<R>            _other_point;
};

template <class R>
inline bool do_intersect(
    const Point_2<R> &p1,
    const Triangle_2<R> &p2)
{
    typedef Point_2_Triangle_2_pair<R> pair_t;
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
Point_2_Triangle_2_pair<R>::
Point_2_Triangle_2_pair()
{
    _known = false;
    _pt = 0;
    _trian = 0;
}

template <class R>
Point_2_Triangle_2_pair<R>::
Point_2_Triangle_2_pair(Point_2<R> const *pt,
                            Triangle_2<R> const *trian)
{
    _known = false;
    _pt = pt;
    _trian = trian;
}

#ifndef CGAL_CFG_RETURN_TYPE_BUG_2
template <class R>
Point_2_Triangle_2_pair<R>::Intersection_results
Point_2_Triangle_2_pair<R>::intersection_type() const
{
    typedef Line_2<R> line_t;
    if (_known)
        return _result;
// The non const this pointer is used to cast away const.
    _known = true;
    if (_trian->has_on_unbounded_side(*_pt)) {
        _result = NO;
    } else {
        _result = POINT;
    }
    return _result;
/*
    line_t l(_trian->vertex(0), _trian->vertex(1));
    if (l.has_on_positive_side(_trian->vertex(2))) {
        for (int i=0; i<3; i++) {
            if (line_t(_trian->vertex(i), _trian->vertex(i+1)).
                                has_on_negative_side(*_pt)) {
                _result = NO;
                return _result;
            }
        }
    } else {
        for (int i=0; i<3; i++)
            if(line_t(_trian->vertex(i), _trian->vertex(i-1)).
                                has_on_negative_side(*_pt)){
                _result = NO;
                return _result;
            }
    }
*/
}

#endif // CGAL_CFG_RETURN_TYPE_BUG_2



template <class R>
bool
Point_2_Triangle_2_pair<R>::
intersection(Point_2<R> &result) const
{
    if (!_known)
        intersection_type();
    if (_result != POINT)
        return false;
    result = *_pt;
    return true;
}

CGAL_END_NAMESPACE



#include <CGAL/Object.h>

CGAL_BEGIN_NAMESPACE

template <class R>
Object
intersection(const Point_2<R> &pt, const Triangle_2<R>&tr)
{
    typedef Point_2_Triangle_2_pair<R> is_t;
    is_t ispair(&pt, &tr);
    switch (ispair.intersection_type()) {
    case is_t::NO:
    default:
        return Object();
    case is_t::POINT: {
        return Object(new Wrapper< Point_2<R> >(pt));
    }
    }
}

template <class R>
class Triangle_2_Point_2_pair
: public Point_2_Triangle_2_pair<R> {
public:
    Triangle_2_Point_2_pair(
            Triangle_2<R> const *trian,
            Point_2<R> const *pt) :
                        Point_2_Triangle_2_pair<R>(pt, trian) {}
};

template <class R>
inline bool do_intersect(
    const Triangle_2<R> &p1,
    const Point_2<R> &p2)
{
    typedef Triangle_2_Point_2_pair<R> pair_t;
    pair_t pair(&p1, &p2);
    return pair.intersection_type() != pair_t::NO;
}

template <class R>
inline Object
intersection(const Triangle_2<R> &tr, const Point_2<R> &pt)
{
    return intersection(pt, tr);
}

CGAL_END_NAMESPACE

#endif
