
// Copyright (c) 2000  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Geert-Jan Giezeman


#ifndef CGAL_LINE_2_TRIANGLE_2_INTERSECTION_H
#define CGAL_LINE_2_TRIANGLE_2_INTERSECTION_H

#include <CGAL/Line_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Object.h>
#include <CGAL/Straight_2.h>
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

template <class K>
class Line_2_Triangle_2_pair {
public:
    enum Intersection_results {NO, POINT, SEGMENT};
    Line_2_Triangle_2_pair() ;
    Line_2_Triangle_2_pair(typename K::Line_2 const *line,
			   typename K::Triangle_2 const *trian);
    ~Line_2_Triangle_2_pair() {}

    Intersection_results intersection_type() const;

    bool                intersection(typename K::Point_2 &result) const;
    bool                intersection(typename K::Segment_2 &result) const;
protected:
    typename K::Line_2 const*_line;
    typename K::Triangle_2 const *  _trian;
    mutable bool                    _known;
    mutable Intersection_results     _result;
    mutable typename K::Point_2         _intersection_point;
    mutable typename K::Point_2        _other_point;
};

template <class K>
inline 
bool 
do_intersect(const typename CGAL_WRAP(K)::Line_2 &p1,
	     const typename CGAL_WRAP(K)::Triangle_2 &p2,
	     const K&)
{
    typedef Line_2_Triangle_2_pair<K> pair_t;
    pair_t pair(&p1, &p2);
    return pair.intersection_type() != pair_t::NO;
}

template <class K>
inline 
bool 
do_intersect(const typename CGAL_WRAP(K)::Triangle_2 &p2,
	     const typename CGAL_WRAP(K)::Line_2 &p1,
	     const K& k)
{
  return CGALi::do_intersect(p1, p2, k);
}

template <class K>
Line_2_Triangle_2_pair<K>::
Line_2_Triangle_2_pair()
{
    _known = false;
    _line = 0;
    _trian = 0;
}

template <class K>
Line_2_Triangle_2_pair<K>::
Line_2_Triangle_2_pair(typename K::Line_2 const *line,
		       typename K::Triangle_2 const *trian)
{
    _known = false;
    _line = line;
    _trian = trian;
}

template <class K>
typename Line_2_Triangle_2_pair<K>::Intersection_results
Line_2_Triangle_2_pair<K>::intersection_type() const
{
  typedef typename K::Line_2 Line_2;
    if (_known)
        return _result;
// The non const this pointer is used to cast away const.
    _known = true;
    Straight_2_<K> straight(*_line);
    Line_2 l(_trian->vertex(0), _trian->vertex(1));
if (l.oriented_side(_trian->vertex(2)) == ON_POSITIVE_SIDE) {
//    if (_trian->is_counterclockwise()) {
        straight.cut_right_off(
            Line_2(_trian->vertex(0), _trian->vertex(1)));
        straight.cut_right_off(
            Line_2(_trian->vertex(1), _trian->vertex(2)));
        straight.cut_right_off(
            Line_2(_trian->vertex(2), _trian->vertex(0)));
    } else {
        straight.cut_right_off(
            Line_2(_trian->vertex(2), _trian->vertex(1)));
        straight.cut_right_off(
            Line_2(_trian->vertex(1), _trian->vertex(0)));
        straight.cut_right_off(
            Line_2(_trian->vertex(0), _trian->vertex(2)));
    }
    switch (straight.current_state()) {
    case Straight_2_<K>::EMPTY:
        _result = NO;
        return _result;
    case Straight_2_<K>::POINT: {
        straight.current(_intersection_point);
        _result = POINT;
        return _result;
        }
    case Straight_2_<K>::SEGMENT: {
        typename K::Segment_2 seg;
        straight.current(seg);
        _intersection_point = seg.source();
        _other_point = seg.target();
        _result = SEGMENT;
        return _result;
        }
    default:  // should not happen.
        CGAL_kernel_assertion_msg(false, "Internal CGAL error.");
        _result = NO;
        return _result;
    }
}


template <class K>
bool
Line_2_Triangle_2_pair<K>::
intersection(typename K::Point_2 &result) const
{
    if (!_known)
        intersection_type();
    if (_result != POINT)
        return false;
    result = _intersection_point;
    return true;
}

template <class K>
bool
Line_2_Triangle_2_pair<K>::
intersection(typename K::Segment_2 &result) const
{
  typedef typename K::Segment_2 Segment_2;
    if (!_known)
        intersection_type();
    if (_result != SEGMENT)
        return false;
    result = Segment_2(_intersection_point, _other_point);
    return true;
}




template <class K>
Object
intersection(const typename CGAL_WRAP(K)::Line_2 &line, 
	     const typename CGAL_WRAP(K)::Triangle_2 &tr,
	     const K&)
{
    typedef Line_2_Triangle_2_pair<K> is_t;
    is_t ispair(&line, &tr);
    switch (ispair.intersection_type()) {
    case is_t::NO:
    default:
        return Object();
    case is_t::POINT: {
        typename K::Point_2 pt;
        ispair.intersection(pt);
        return make_object(pt);
    }
    case is_t::SEGMENT: {
        typename K::Segment_2 iseg;
        ispair.intersection(iseg);
        return make_object(iseg);
    }
    }
}


template <class K>
inline
Object
intersection(const typename CGAL_WRAP(K)::Triangle_2 &tr,
	     const typename CGAL_WRAP(K)::Line_2 &line,
	     const K& k)
{
  return intersection(line, tr, k);
}


template <class K>
class Triangle_2_Line_2_pair
: public Line_2_Triangle_2_pair<K> {
public:
    Triangle_2_Line_2_pair(
            typename K::Triangle_2 const *trian,
            typename K::Line_2 const *line) :
                        Line_2_Triangle_2_pair<K>(line, trian) {}
};

} // namespace CGALi

template <class K>
inline 
bool 
do_intersect(const Line_2<K> &p1,
	     const Triangle_2<K> &p2)
{
  typedef typename K::Do_intersect_2 Do_intersect;
  return Do_intersect()(p1, p2);
}

template <class K>
inline bool do_intersect(
    const Triangle_2<K> &p1,
    const Line_2<K> &p2)
{
  typedef typename K::Do_intersect_2 Do_intersect;
    return Do_intersect()(p2, p1);
}

template <class K>
inline Object
intersection(const Line_2<K> &line, const Triangle_2<K> &tr)
{
  typedef typename K::Intersect_2 Intersect;
    return Intersect()(line, tr);
}

template <class K>
inline Object
intersection(const Triangle_2<K> &tr, const Line_2<K> &line)
{
  typedef typename K::Intersect_2 Intersect;
    return Intersect()(line, tr);
}

CGAL_END_NAMESPACE

#endif
