
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
// file          : include/CGAL/Segment_2_Iso_rectangle_2_intersection.h
// source        : intersection_2_2.fw
// author(s)     : Geert-Jan Giezeman
//
// coordinator   : Saarbruecken
//
// ============================================================================


#ifndef CGAL_SEGMENT_2_ISO_RECTANGLE_2_INTERSECTION_H
#define CGAL_SEGMENT_2_ISO_RECTANGLE_2_INTERSECTION_H

#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>

CGAL_BEGIN_NAMESPACE

template <class R>
class Segment_2_Iso_rectangle_2_pair {
public:
    enum Intersection_results {NO, POINT, SEGMENT};
    Segment_2_Iso_rectangle_2_pair() ;
    Segment_2_Iso_rectangle_2_pair(Segment_2<R> const *seg,
                          Iso_rectangle_2<R> const *rect) ;

#ifndef CGAL_CFG_RETURN_TYPE_BUG_2
  Intersection_results intersection_type() const;

#else
  Intersection_results intersection_type() const
 {
     typedef typename R::RT RT;
     typedef typename R::FT FT;
     if (_known)
         return _result;
     _known = true;
     for (int i=0; i<_ref_point.dimension(); i++) {
         if (_dir.homogeneous(i) == RT(0)) {
             if (_ref_point.cartesian(i) < _isomin.cartesian(i)) {
                 _result = NO;
                 return _result;
             }
             if (_ref_point.cartesian(i) > _isomax.cartesian(i)) {
                 _result = NO;
                 return _result;
             }
         } else {
             FT newmin, newmax;
             if (_dir.homogeneous(i) > RT(0)) {
                 newmin = (_isomin.cartesian(i)-_ref_point.cartesian(i)) /
                     _dir.cartesian(i);
                 newmax = (_isomax.cartesian(i)-_ref_point.cartesian(i)) /
                     _dir.cartesian(i);
             } else {
                 newmin = (_isomax.cartesian(i)-_ref_point.cartesian(i)) /
                     _dir.cartesian(i);
                 newmax = (_isomin.cartesian(i)-_ref_point.cartesian(i)) /
                     _dir.cartesian(i);
             }
             if (newmin > _min)
                 _min = newmin;
             if (newmax < _max)
                 _max = newmax;
             if (_max < _min) {
                 _result = NO;
                 return _result;
             }
         }
     }
     if (_max == _min) {
         _result = POINT;
         return _result;
     }
     _result = SEGMENT;
     return _result;
 }
 
#endif // CGAL_CFG_RETURN_TYPE_BUG_2

    bool                       intersection(
                                    Point_2<R> &result) const;
    bool                       intersection(
                                    Segment_2<R> &result) const;
protected:
    mutable bool                       _known;
    mutable Intersection_results       _result;
    mutable Point_2<R>            _ref_point;
    mutable Vector_2<R>           _dir;
    mutable Point_2<R>            _isomin;
    mutable Point_2<R>            _isomax;
    mutable typename R::FT              _min,
                               _max;
};

template <class R>
inline bool do_intersect(
    const Segment_2<R> &p1,
    const Iso_rectangle_2<R> &p2)
{
    typedef Segment_2_Iso_rectangle_2_pair<R> pair_t;
    pair_t pair(&p1, &p2);
    return pair.intersection_type() != pair_t::NO;
}

CGAL_END_NAMESPACE

#include <CGAL/Object.h>

CGAL_BEGIN_NAMESPACE

template <class R>
Object
intersection(
    const Segment_2<R> &seg,
    const Iso_rectangle_2<R> &iso)
{
    typedef Segment_2_Iso_rectangle_2_pair<R> is_t;
    is_t ispair(&seg, &iso);
    switch (ispair.intersection_type()) {
    case is_t::NO:
    default:
        return Object();
    case is_t::POINT: {
        Point_2<R> ipt;
        ispair.intersection(ipt);
        return Object(new Wrapper< Point_2<R> >(ipt));
    }
    case is_t::SEGMENT: {
        Segment_2<R> iseg;
        ispair.intersection(iseg);
        return Object(new Wrapper< Segment_2<R> >(iseg));
    }
    }
}

CGAL_END_NAMESPACE




CGAL_BEGIN_NAMESPACE

template <class R>
Segment_2_Iso_rectangle_2_pair<R>::Segment_2_Iso_rectangle_2_pair()
{
    _known = false;
}

template <class R>
Segment_2_Iso_rectangle_2_pair<R>::
Segment_2_Iso_rectangle_2_pair(
        Segment_2<R> const *seg,
        Iso_rectangle_2<R> const *iso)
{
    _known = false;
    _isomin = iso->min();
    _isomax = iso->max();
    _ref_point = seg->source();
    _dir = seg->direction().to_vector();
    _min = (typename R::FT)(0);
    int main_dir = (CGAL_NTS abs(_dir.x()) > CGAL_NTS abs(_dir.y()) ) ? 0 : 1;
    _max = (seg->target().cartesian(main_dir)-_ref_point.cartesian(main_dir)) /
            _dir.cartesian(main_dir);
}

#ifndef CGAL_CFG_RETURN_TYPE_BUG_2
template <class R>
Segment_2_Iso_rectangle_2_pair<R>::Intersection_results
Segment_2_Iso_rectangle_2_pair<R>::intersection_type() const
{
    typedef typename R::RT RT;
    typedef typename R::FT FT;
    if (_known)
        return _result;
    _known = true;
    for (int i=0; i<_ref_point.dimension(); i++) {
        if (_dir.homogeneous(i) == RT(0)) {
            if (_ref_point.cartesian(i) < _isomin.cartesian(i)) {
                _result = NO;
                return _result;
            }
            if (_ref_point.cartesian(i) > _isomax.cartesian(i)) {
                _result = NO;
                return _result;
            }
        } else {
            FT newmin, newmax;
            if (_dir.homogeneous(i) > RT(0)) {
                newmin = (_isomin.cartesian(i)-_ref_point.cartesian(i)) /
                    _dir.cartesian(i);
                newmax = (_isomax.cartesian(i)-_ref_point.cartesian(i)) /
                    _dir.cartesian(i);
            } else {
                newmin = (_isomax.cartesian(i)-_ref_point.cartesian(i)) /
                    _dir.cartesian(i);
                newmax = (_isomin.cartesian(i)-_ref_point.cartesian(i)) /
                    _dir.cartesian(i);
            }
            if (newmin > _min)
                _min = newmin;
            if (newmax < _max)
                _max = newmax;
            if (_max < _min) {
                _result = NO;
                return _result;
            }
        }
    }
    if (_max == _min) {
        _result = POINT;
        return _result;
    }
    _result = SEGMENT;
    return _result;
}

#endif // CGAL_CFG_RETURN_TYPE_BUG_2

template <class R>
bool Segment_2_Iso_rectangle_2_pair<R>::
intersection(Segment_2<R> &seg) const
{
    if (!_known)
        intersection_type();
    if (_result != SEGMENT)
        return false;
    Point_2<R> p1(_ref_point + _dir*_min);
    Point_2<R> p2(_ref_point + _dir*_max);
    seg = Segment_2<R>(p1, p2);
    return true;
}

template <class R> bool Segment_2_Iso_rectangle_2_pair<R>::
intersection(Point_2<R> &pt) const
{
    if (!_known)
        intersection_type();
    if (_result != POINT)
        return false;
    pt = Point_2<R>(_ref_point + _dir*_min);
    return true;
}

CGAL_END_NAMESPACE



CGAL_BEGIN_NAMESPACE

template <class R>
class Iso_rectangle_2_Segment_2_pair:
          public Segment_2_Iso_rectangle_2_pair<R> {
public:
    Iso_rectangle_2_Segment_2_pair() {}
    Iso_rectangle_2_Segment_2_pair(Iso_rectangle_2<R> const *rect,
                               Segment_2<R> const *seg)
                :Segment_2_Iso_rectangle_2_pair<R> (seg, rect){}
};

template <class R>
inline bool do_intersect(
    const Iso_rectangle_2<R> &p1,
    const Segment_2<R> &p2)
{
    typedef Iso_rectangle_2_Segment_2_pair<R> pair_t;
    pair_t pair(&p1, &p2);
    return pair.intersection_type() != pair_t::NO;
}

template <class R>
inline Object
intersection(
    const Iso_rectangle_2<R>&iso,
    const Segment_2<R>&seg)
{
    return intersection(seg, iso);
}

CGAL_END_NAMESPACE

#endif // CGAL_SEGMENT_2_ISO_RECTANGLE_2_INTERSECTION_H
