
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
// file          : include/CGAL/Bbox_2_Line_2_intersection.h
// source        : intersection_2_1.fw
// author(s)     : Geert-Jan Giezeman
//
// coordinator   : MPI, Saarbruecken
//
// ============================================================================


#ifndef CGAL_BBOX_2_LINE_2_INTERSECTION_H
#define CGAL_BBOX_2_LINE_2_INTERSECTION_H

#include <CGAL/Cartesian.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>

CGAL_BEGIN_NAMESPACE

template <class R>
class Bbox_2_Line_2_pair {
public:
    typedef Cartesian<double> Rcart;
    enum Intersection_results {NO, POINT, SEGMENT};
    Bbox_2_Line_2_pair() ;
    Bbox_2_Line_2_pair(Bbox_2 const *bbox,
                            Line_2<R> const *line);
    ~Bbox_2_Line_2_pair() {}

#ifndef CGAL_CFG_RETURN_TYPE_BUG_2
  Intersection_results intersection_type() const;

#else

  Intersection_results intersection_type() const
{
    if (_known)
        return _result;
    // The non const this pointer is used to cast away const.
    _known = true;
    const Point_2< Cartesian<double> > &ref_point = _line.point();
    const Vector_2< Cartesian<double> > &dir =
                               _line.direction().to_vector();
    bool to_infinity = true;
// first on x value
    if (dir.x() == 0.0) {
        if (ref_point.x() < _bbox->xmin()) {
            _result = NO;
            return _result;
        }
        if (ref_point.x() > _bbox->xmax()) {
            _result = NO;
            return _result;
        }
    } else {
        double newmin, newmax;
        if (dir.x() > 0.0) {
            newmin = (_bbox->xmin()-ref_point.x())/dir.x();
            newmax = (_bbox->xmax()-ref_point.x())/dir.x();
        } else {
            newmin = (_bbox->xmax()-ref_point.x())/dir.x();
            newmax = (_bbox->xmin()-ref_point.x())/dir.x();
        }
        if (to_infinity) {
            _min = newmin;
            _max = newmax;
        } else {
            if (newmin > _min)
                _min = newmin;
            if (newmax < _max)
                _max = newmax;
            if (_max < _min) {
                _result = NO;
                return _result;
            }
        }
        to_infinity = false;
    }
// now on y value
    if (dir.y() == 0.0) {
        if (ref_point.y() < _bbox->ymin()) {
            _result = NO;
            return _result;
        }
        if (ref_point.y() > _bbox->ymax()) {
            _result = NO;
            return _result;
        }
    } else {
        double newmin, newmax;
        if (dir.y() > 0.0) {
            newmin = (_bbox->ymin()-ref_point.y())/dir.y();
            newmax = (_bbox->ymax()-ref_point.y())/dir.y();
        } else {
            newmin = (_bbox->ymax()-ref_point.y())/dir.y();
            newmax = (_bbox->ymin()-ref_point.y())/dir.y();
        }
        if (to_infinity) {
            _min = newmin;
            _max = newmax;
        } else {
            if (newmin > _min)
                _min = newmin;
            if (newmax < _max)
                _max = newmax;
            if (_max < _min) {
                _result = NO;
                return _result;
            }
        }
        to_infinity = false;
    }
    CGAL_kernel_assertion(!to_infinity);
    if (_max == _min) {
        _result = POINT;
        return _result;
    }
    _result = SEGMENT;
    return _result;
}

#endif // CGAL_CFG_RETURN_TYPE_BUG_2

    bool intersection(Point_2< Rcart > &result) const;
    bool intersection(Segment_2< Rcart > &result) const;
protected:
    Bbox_2 const *      _bbox;
    Line_2< Rcart > _line;
    mutable bool                     _known;
    mutable Intersection_results     _result;
    mutable double                   _min, _max;
};

template <class R>
inline bool do_intersect(
    const Bbox_2 &box,
    const Line_2<R> &line)
{
    typedef Bbox_2_Line_2_pair<R> pair_t;
    pair_t pair(&box, &line);
    return pair.intersection_type() != pair_t::NO;
}

CGAL_END_NAMESPACE



CGAL_BEGIN_NAMESPACE

template <class R>
Bbox_2_Line_2_pair<R>::Bbox_2_Line_2_pair()
{
    _bbox = 0;
    _known = false;
}

template <class R>
Bbox_2_Line_2_pair<R>::Bbox_2_Line_2_pair(
    Bbox_2 const *bbox, Line_2<R> const *line)
{
    _bbox = bbox;
    _line = Line_2< Rcart >(to_double(line->a()),
            to_double(line->b()), to_double(line->c()));
    _known = false;
}

#ifndef CGAL_CFG_RETURN_TYPE_BUG_2
template <class R>
Bbox_2_Line_2_pair<R>::Intersection_results
Bbox_2_Line_2_pair<R>::intersection_type() const
{
    if (_known)
        return _result;
    // The non const this pointer is used to cast away const.
    _known = true;
    const Point_2< Cartesian<double> > &ref_point = _line.point();
    const Vector_2< Cartesian<double> > &dir =
                               _line.direction().to_vector();
    bool to_infinity = true;
// first on x value
    if (dir.x() == 0.0) {
        if (ref_point.x() < _bbox->xmin()) {
            _result = NO;
            return _result;
        }
        if (ref_point.x() > _bbox->xmax()) {
            _result = NO;
            return _result;
        }
    } else {
        double newmin, newmax;
        if (dir.x() > 0.0) {
            newmin = (_bbox->xmin()-ref_point.x())/dir.x();
            newmax = (_bbox->xmax()-ref_point.x())/dir.x();
        } else {
            newmin = (_bbox->xmax()-ref_point.x())/dir.x();
            newmax = (_bbox->xmin()-ref_point.x())/dir.x();
        }
        if (to_infinity) {
            _min = newmin;
            _max = newmax;
        } else {
            if (newmin > _min)
                _min = newmin;
            if (newmax < _max)
                _max = newmax;
            if (_max < _min) {
                _result = NO;
                return _result;
            }
        }
        to_infinity = false;
    }
// now on y value
    if (dir.y() == 0.0) {
        if (ref_point.y() < _bbox->ymin()) {
            _result = NO;
            return _result;
        }
        if (ref_point.y() > _bbox->ymax()) {
            _result = NO;
            return _result;
        }
    } else {
        double newmin, newmax;
        if (dir.y() > 0.0) {
            newmin = (_bbox->ymin()-ref_point.y())/dir.y();
            newmax = (_bbox->ymax()-ref_point.y())/dir.y();
        } else {
            newmin = (_bbox->ymax()-ref_point.y())/dir.y();
            newmax = (_bbox->ymin()-ref_point.y())/dir.y();
        }
        if (to_infinity) {
            _min = newmin;
            _max = newmax;
        } else {
            if (newmin > _min)
                _min = newmin;
            if (newmax < _max)
                _max = newmax;
            if (_max < _min) {
                _result = NO;
                return _result;
            }
        }
        to_infinity = false;
    }
    CGAL_kernel_assertion(!to_infinity);
    if (_max == _min) {
        _result = POINT;
        return _result;
    }
    _result = SEGMENT;
    return _result;
}

#endif // CGAL_CFG_RETURN_TYPE_BUG_2

template <class R>
bool
Bbox_2_Line_2_pair<R>::intersection(
    Segment_2< Cartesian<double> > &seg) const
{
    if (!_known)
        intersection_type();
    if (_result != SEGMENT)
        return false;
    Point_2< Cartesian<double> > p1(_line.point()
                + _min*_line.direction().to_vector());
    Point_2< Cartesian<double> > p2(_line.point()
                + _max*_line.direction().to_vector());
    seg = Segment_2< Cartesian<double> >(p1, p2);
    return true;
}

template <class R>
bool
Bbox_2_Line_2_pair<R>::intersection(
    Point_2< Cartesian<double> > &pt) const
{
    if (!_known)
        intersection_type();
    if (_result != POINT)
        return false;
    pt = (_line.point() + _min*_line.direction().to_vector());
    return true;
}

CGAL_END_NAMESPACE



CGAL_BEGIN_NAMESPACE

template <class R>
class Line_2_Bbox_2_pair: public Bbox_2_Line_2_pair<R> {
public:
    Line_2_Bbox_2_pair(
            Line_2<R> const *line,
            Bbox_2 const *bbox) :
                                Bbox_2_Line_2_pair<R>(bbox, line) {}
};

template <class R>
inline bool do_intersect(
    const Line_2<R> &line,
    const Bbox_2 &box)
{
    typedef Line_2_Bbox_2_pair<R> pair_t;
    pair_t pair(&line, &box);
    return pair.intersection_type() != pair_t::NO;
}

CGAL_END_NAMESPACE

#endif
