
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
// file          : include/CGAL/Line_2_Iso_rectangle_2_intersection.h
// source        : intersection_2_2.fw
// author(s)     : Geert-Jan Giezeman
//
// coordinator   : Saarbruecken
//
// ============================================================================


#ifndef CGAL_LINE_2_ISO_RECTANGLE_2_INTERSECTION_H
#define CGAL_LINE_2_ISO_RECTANGLE_2_INTERSECTION_H

#include <CGAL/Line_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/Object.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

template <class K>
class Line_2_Iso_rectangle_2_pair {
public:
    enum Intersection_results {NO, POINT, SEGMENT};
    Line_2_Iso_rectangle_2_pair() ;
    Line_2_Iso_rectangle_2_pair(typename K::Line_2 const *pt,
                            typename K::Iso_rectangle_2 const *iso);
    ~Line_2_Iso_rectangle_2_pair() {}
#ifndef CGAL_CFG_RETURN_TYPE_BUG_2
  Intersection_results intersection_type() const;
#else
  Intersection_results intersection_type() const
{
    typedef K::Line_2 line_t;
    if (_known)
        return _result;
// The non const this pointer is used to cast away const.
    _known = true;
    typedef typename K::FT FT;
    typedef typename K::RT RT;
    bool all_values = true;
    int i;
    for (i=0; i< _ref_point.dimension(); i++) {
        if (_dir.homogeneous(i) == RT(0)) {
            if (_ref_point.cartesian(i) < _isomin.cartesian(i)) {
                _result = NO;
                return NO;
            }
            if (_ref_point.cartesian(i) > _isomax.cartesian(i)) {
                _result = NO;
                return NO;
            }
        } else {
            FT newmin, newmax;
            if (_dir.homogeneous(i) > RT(0)) {
                newmin = (_isomin.cartesian(i) - _ref_point.cartesian(i)) /
                    _dir.cartesian(i);
                newmax = (_isomax.cartesian(i) - _ref_point.cartesian(i)) /
                    _dir.cartesian(i);
            } else {
                newmin = (_isomax.cartesian(i) - _ref_point.cartesian(i)) /
                    _dir.cartesian(i);
                newmax = (_isomin.cartesian(i) - _ref_point.cartesian(i)) /
                    _dir.cartesian(i);
            }
            if (all_values) {
                _min = newmin;
                _max = newmax;
            } else {
                if (newmin > _min)
                    _min = newmin;
                if (newmax < _max)
                    _max = newmax;
                if (_max < _min) {
                    _result = NO;
                    return NO;
                }
            }
            all_values = false;
        }
    }
    CGAL_kernel_assertion(!all_values);
    if (_max == _min) {
        _result = POINT;
        return POINT;
    }
    _result = SEGMENT;
    return SEGMENT;
}

#endif // CGAL_CFG_RETURN_TYPE_BUG_2
    bool                intersection(typename K::Point_2 &result) const;
    bool                intersection(typename K::Segment_2 &result) const;
protected:
    typename K::Point_2             _ref_point;
    typename K::Vector_2            _dir;
    typename K::Point_2             _isomin;
    typename K::Point_2             _isomax;
    mutable bool                        _known;
    mutable Intersection_results        _result;
    mutable typename K::FT              _min, _max;
};

template <class K>
inline bool do_intersect(const typename CGAL_WRAP(K)::Line_2 &p1,
			 const typename CGAL_WRAP(K)::Iso_rectangle_2 &p2,
			 const K&)
{
    typedef Line_2_Iso_rectangle_2_pair<K> pair_t;
    pair_t pair(&p1, &p2);
    return pair.intersection_type() != pair_t::NO;
}

template <class K>
inline bool do_intersect(const typename CGAL_WRAP(K)::Iso_rectangle_2 &p2,
			 const typename CGAL_WRAP(K)::Line_2 &p1,
			 const K& k)
{
  return CGALi::do_intersect(p1, p2, k);
}


template <class K>
Line_2_Iso_rectangle_2_pair<K>::
Line_2_Iso_rectangle_2_pair()
{
    _known = false;
}

template <class K>
Line_2_Iso_rectangle_2_pair<K>::
Line_2_Iso_rectangle_2_pair(typename K::Line_2 const *line,
                            typename K::Iso_rectangle_2 const *iso)
  : _known(false),
    _ref_point(line->point()),
    _dir(line->direction().to_vector()),
    _isomin(iso->min()),
    _isomax(iso->max())
{}

#ifndef CGAL_CFG_RETURN_TYPE_BUG_2
template <class K>
typename Line_2_Iso_rectangle_2_pair<K>::Intersection_results
Line_2_Iso_rectangle_2_pair<K>::intersection_type() const
{
    typedef typename K::Line_2 line_t;
    if (_known)
        return _result;
// The non const this pointer is used to cast away const.
    _known = true;
    typedef typename K::FT FT;
    typedef typename K::RT RT;
    bool all_values = true;
    int i;
    for (i=0; i< _ref_point.dimension(); i++) {
        if (_dir.homogeneous(i) == RT(0)) {
            if (_ref_point.cartesian(i) < _isomin.cartesian(i)) {
                _result = NO;
                return NO;
            }
            if (_ref_point.cartesian(i) > _isomax.cartesian(i)) {
                _result = NO;
                return NO;
            }
        } else {
            FT newmin, newmax;
            if (_dir.homogeneous(i) > RT(0)) {
                newmin = (_isomin.cartesian(i) - _ref_point.cartesian(i)) /
                    _dir.cartesian(i);
                newmax = (_isomax.cartesian(i) - _ref_point.cartesian(i)) /
                    _dir.cartesian(i);
            } else {
                newmin = (_isomax.cartesian(i) - _ref_point.cartesian(i)) /
                    _dir.cartesian(i);
                newmax = (_isomin.cartesian(i) - _ref_point.cartesian(i)) /
                    _dir.cartesian(i);
            }
            if (all_values) {
                _min = newmin;
                _max = newmax;
            } else {
                if (newmin > _min)
                    _min = newmin;
                if (newmax < _max)
                    _max = newmax;
                if (_max < _min) {
                    _result = NO;
                    return NO;
                }
            }
            all_values = false;
        }
    }
    CGAL_kernel_assertion(!all_values);
    if (_max == _min) {
        _result = POINT;
        return POINT;
    }
    _result = SEGMENT;
    return SEGMENT;
}

#endif


template <class K>
bool
Line_2_Iso_rectangle_2_pair<K>::
intersection(typename K::Point_2 &result) const
{
    if (!_known)
        intersection_type();
    if (_result != POINT)
        return false;
    result = _ref_point + _dir * _min;
    return true;
}

template <class K>
bool
Line_2_Iso_rectangle_2_pair<K>::
intersection(typename K::Segment_2 &result) const
{
  typename K::Construct_segment_2 construct_segment_2;
    if (!_known)
        intersection_type();
    if (_result != SEGMENT)
        return false;
    result = construct_segment_2(_ref_point + _dir*_min, 
				 _ref_point + _dir*_max);
    return true;
}



template <class K>
Object
intersection(const typename CGAL_WRAP(K)::Line_2 &line, 
	     const typename CGAL_WRAP(K)::Iso_rectangle_2 &iso,
	     const K&)
{
    typename K::Construct_object_2 construct_object;
    typedef Line_2_Iso_rectangle_2_pair<K> is_t;
    is_t ispair(&line, &iso);
    switch (ispair.intersection_type()) {
    case is_t::NO:
    default:
        return Object();
    case is_t::POINT: {
        typename K::Point_2 ipt;
        ispair.intersection(ipt);
        return construct_object(ipt);
    }
    case is_t::SEGMENT: {
        typename K::Segment_2 iseg;
        ispair.intersection(iseg);
        return construct_object(iseg);
    }
    }
}

template <class K>
inline
Object
intersection(const typename CGAL_WRAP(K)::Iso_rectangle_2 &iso,
	     const typename CGAL_WRAP(K)::Line_2 &line, 
	     const K& k)
{
  return CGALi::intersection(line, iso, k);
}

template <class K>
class Iso_rectangle_2_Line_2_pair
: public Line_2_Iso_rectangle_2_pair<K> {
public:
    Iso_rectangle_2_Line_2_pair(
            typename K::Iso_rectangle_2 const *iso,
            typename K::Line_2 const *line) :
                        Line_2_Iso_rectangle_2_pair<K>(line, iso) {}
};

} // namespace CGALi




template <class K>
inline bool do_intersect(
    const Line_2<K> &p1,
    const Iso_rectangle_2<K> &p2)
{
  return CGALi::do_intersect(p1, p2, K());
}

template <class K>
inline bool do_intersect(
    const Iso_rectangle_2<K> &p1,
    const Line_2<K> &p2)
{
  return CGALi::do_intersect(p2, p1, K());
}

template <class K>
Object
intersection(const Line_2<K> &line, const Iso_rectangle_2<K> &iso)
{
  return CGALi::intersection(line, iso, K());
}

template <class K>
inline Object
intersection(const Iso_rectangle_2<K> &iso, const Line_2<K> &line)
{
    return CGALi::intersection(line, iso, K());
}

CGAL_END_NAMESPACE

#endif
