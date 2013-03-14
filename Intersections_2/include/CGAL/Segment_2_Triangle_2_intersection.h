// Copyright (c) 2000  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Geert-Jan Giezeman


#ifndef CGAL_SEGMENT_2_TRIANGLE_2_INTERSECTION_H
#define CGAL_SEGMENT_2_TRIANGLE_2_INTERSECTION_H

#include <CGAL/Segment_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/kernel_assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/Straight_2.h>
#include <CGAL/Intersection_traits_2.h>

namespace CGAL {

namespace internal {

template <class K>
class Segment_2_Triangle_2_pair {
public:
    enum Intersection_results {NO_INTERSECTION, POINT, SEGMENT};
    Segment_2_Triangle_2_pair(typename K::Segment_2 const *seg,
                              typename K::Triangle_2 const *trian)
      : _seg(seg), _trian(trian), _known(false) {}

    Intersection_results intersection_type() const;

    typename K::Point_2    intersection_point() const;
    typename K::Segment_2  intersection_segment() const;
protected:
    typename K::Segment_2 const *  _seg;
    typename K::Triangle_2 const * _trian;
    mutable bool                       _known;
    mutable Intersection_results       _result;
    mutable typename K::Point_2            _intersection_point;
    mutable typename K::Point_2            _other_point;
};

template <class K>
inline bool do_intersect(
    const typename K::Segment_2 &p1,
    const typename K::Triangle_2 &p2,
    const K&)
{
    typedef Segment_2_Triangle_2_pair<K> pair_t;
    pair_t pair(&p1, &p2);
    return pair.intersection_type() != pair_t::NO_INTERSECTION;
}





template <class K>
typename Segment_2_Triangle_2_pair<K>::Intersection_results
Segment_2_Triangle_2_pair<K>::intersection_type() const
{
    if (_known)
        return _result;
// The non const this pointer is used to cast away const.
    _known = true;
    Straight_2_<K> straight(*_seg);
    typedef typename K::Line_2 Line_2;

Line_2 l(_trian->vertex(0), _trian->vertex(1));
if (l.oriented_side(_trian->vertex(2)) == ON_POSITIVE_SIDE) {
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
        _result = NO_INTERSECTION;
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
        _result = NO_INTERSECTION;
        return _result;
    }
}


template <class K>
typename K::Point_2
Segment_2_Triangle_2_pair<K>::
intersection_point() const
{
    if (!_known)
        intersection_type();
    CGAL_kernel_assertion(_result == POINT);
    return _intersection_point;
}

template <class K>
typename K::Segment_2
Segment_2_Triangle_2_pair<K>::
intersection_segment() const
{
  typedef typename K::Segment_2 Segment_2; 
    if (!_known)
        intersection_type();
    CGAL_kernel_assertion(_result == SEGMENT);
    return Segment_2(_intersection_point, _other_point);
}




template <class K>
typename CGAL::Intersection_traits
<K, typename K::Segment_2, typename K::Triangle_2>::result_type
intersection(const typename K::Segment_2 &seg, 
	     const typename K::Triangle_2&tr,
	     const K&)
{
    typedef Segment_2_Triangle_2_pair<K> is_t;
    is_t ispair(&seg, &tr);
    switch (ispair.intersection_type()) {
    case is_t::NO_INTERSECTION:
    default:
        return intersection_return<typename K::Intersect_2, typename K::Segment_2, typename K::Triangle_2>();
    case is_t::POINT:
        return intersection_return<typename K::Intersect_2, typename K::Segment_2, typename K::Triangle_2>(ispair.intersection_point());
    case is_t::SEGMENT:
        return intersection_return<typename K::Intersect_2, typename K::Segment_2, typename K::Triangle_2>(ispair.intersection_segment());
    }
}


template <class K>
typename CGAL::Intersection_traits
<K, typename K::Segment_2, typename K::Triangle_2>::result_type
intersection(const typename K::Triangle_2&tr,
	     const typename K::Segment_2 &seg, 
	     const K& k)
{
  return internal::intersection(seg, tr, k);
}


template <class K>
inline bool do_intersect(
    const typename K::Triangle_2 &p1,
    const typename K::Segment_2 &p2,
    const K&)
{
    typedef Segment_2_Triangle_2_pair<K> pair_t;
    pair_t pair(&p2, &p1);
    return pair.intersection_type() != pair_t::NO_INTERSECTION;
}

} // namespace internal

CGAL_INTERSECTION_FUNCTION(Segment_2, Triangle_2, 2)
CGAL_DO_INTERSECT_FUNCTION(Segment_2, Triangle_2, 2)

} //namespace CGAL

#endif
