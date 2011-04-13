
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
// file          : include/CGAL/Line_2_Line_2_intersection.h
// source        : intersection_2_1.fw
// author(s)     : Geert-Jan Giezeman
//
// coordinator   : Saarbruecken
//
// ============================================================================


#ifndef CGAL_LINE_2_LINE_2_INTERSECTION_H
#define CGAL_LINE_2_LINE_2_INTERSECTION_H

#include <CGAL/Line_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>

CGAL_BEGIN_NAMESPACE

template <class R>
class Line_2_Line_2_pair {
public:
    enum Intersection_results {NO, POINT, LINE};
    Line_2_Line_2_pair() ;
    Line_2_Line_2_pair(Line_2<R> const *line1,
                            Line_2<R> const *line2);
    ~Line_2_Line_2_pair() {}

#ifndef CGAL_CFG_RETURN_TYPE_BUG_2
    Intersection_results intersection_type() const;

#else
    Intersection_results intersection_type() const
{
    typedef typename R::RT RT;
    if (_known)
        return _result;
    RT nom1, nom2, denom;
    // The non const this pointer is used to cast away const.
    _known = true;
    denom = _line1->a()*_line2->b() - _line2->a()*_line1->b();
    if (denom == RT(0)) {
        if (RT(0) == (_line1->a()*_line2->c() - _line2->a()*_line1->c()) &&
            RT(0) == (_line1->b()*_line2->c() - _line2->b()*_line1->c()))
            _result = LINE;
        else
            _result = NO;
        return _result;
    }
    nom1 = (_line1->b()*_line2->c() - _line2->b()*_line1->c());
    if (!::CGAL::is_finite(nom1)) {
        _result = NO;
        return _result;
    }
    nom2 = (_line2->a()*_line1->c() - _line1->a()*_line2->c());
    if (!::CGAL::is_finite(nom2)) {
        _result = NO;
        return _result;
    }
    R dummyR;
    if (!construct_if_finite(_intersection_point,
                            nom1, nom2, denom, dummyR)){
        _result = NO;
        return _result;
    }
    _result = POINT;
    return _result;
}

#endif // CGAL_CFG_RETURN_TYPE_BUG_2

    bool                intersection(Point_2<R> &result) const;
    bool                intersection(Line_2<R> &result) const;
protected:
    Line_2<R> const*   _line1;
    Line_2<R> const *  _line2;
    mutable bool                    _known;
    mutable Intersection_results    _result;
    mutable Point_2<R>         _intersection_point;
};

template <class R>
inline bool do_intersect(
    const Line_2<R> &p1,
    const Line_2<R> &p2)
{
    typedef Line_2_Line_2_pair<R> pair_t;
    pair_t pair(&p1, &p2);
    return pair.intersection_type() != pair_t::NO;
}

CGAL_END_NAMESPACE

#include <CGAL/Object.h>

CGAL_BEGIN_NAMESPACE

template <class R>
Object
intersection(const Line_2<R> &line1, const Line_2<R> &line2)
{
    typedef Line_2_Line_2_pair<R> is_t;
    is_t linepair(&line1, &line2);
    switch (linepair.intersection_type()) {
    case is_t::NO:
    default:
        return Object();
    case is_t::POINT: {
        Point_2<R> pt;
        linepair.intersection(pt);
        return Object(new Wrapper< Point_2<R> >(pt));
    }
    case is_t::LINE:
        return Object(new Wrapper< Line_2<R> >(line1));
    }
}

CGAL_END_NAMESPACE



CGAL_BEGIN_NAMESPACE

template <class R, class POINT, class RT>
bool construct_if_finite(POINT &pt, RT x, RT y, RT w, R &)
{
    typedef typename R::FT FT;
    CGAL_kernel_precondition(::CGAL::is_finite(x)
                             && ::CGAL::is_finite(y)
                             && w != RT(0));

    if (!::CGAL::is_finite(FT(x)/FT(w)) || !::CGAL::is_finite(FT(y)/FT(w)))
        return false;
    pt = POINT(x, y, w);
    return true;
}

CGAL_END_NAMESPACE


CGAL_BEGIN_NAMESPACE

template <class R>
Line_2_Line_2_pair<R>::Line_2_Line_2_pair()
{
    _line1 = 0;
    _line2 = 0;
    _known = false;
}

template <class R>
Line_2_Line_2_pair<R>::Line_2_Line_2_pair(
    Line_2<R> const *line1, Line_2<R> const *line2)
{
    _line1 = line1;
    _line2 = line2;
    _known = false;
}

#ifndef CGAL_CFG_RETURN_TYPE_BUG_2
template <class R>
Line_2_Line_2_pair<R>::Intersection_results
Line_2_Line_2_pair<R>::intersection_type() const
{
    typedef typename R::RT RT;
    if (_known)
        return _result;
    RT nom1, nom2, denom;
    // The non const this pointer is used to cast away const.
    _known = true;
    denom = _line1->a()*_line2->b() - _line2->a()*_line1->b();
    if (denom == RT(0)) {
        if (RT(0) == (_line1->a()*_line2->c() - _line2->a()*_line1->c()) &&
            RT(0) == (_line1->b()*_line2->c() - _line2->b()*_line1->c()))
            _result = LINE;
        else
            _result = NO;
        return _result;
    }
    nom1 = (_line1->b()*_line2->c() - _line2->b()*_line1->c());
    if (!::CGAL::is_finite(nom1)) {
        _result = NO;
        return _result;
    }
    nom2 = (_line2->a()*_line1->c() - _line1->a()*_line2->c());
    if (!::CGAL::is_finite(nom2)) {
        _result = NO;
        return _result;
    }
    R dummyR;
    if (!construct_if_finite(_intersection_point,
                            nom1, nom2, denom, dummyR)){
        _result = NO;
        return _result;
    }
    _result = POINT;
    return _result;
}

#endif // CGAL_CFG_RETURN_TYPE_BUG_2

template <class R>
bool
Line_2_Line_2_pair<R>::intersection(Point_2<R> &pt) const
{
    if (!_known)
        intersection_type();
    if (_result != POINT)
        return false;
    pt = _intersection_point;
    return true;
}

template <class R>
bool
Line_2_Line_2_pair<R>::intersection(Line_2<R> &l) const
{
    if (!_known)
        intersection_type();
    if (_result != LINE)
        return false;
    l = *_line1;
    return true;
}

CGAL_END_NAMESPACE



#endif
