
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


#ifndef CGAL_RAY_2_SEGMENT_2_INTERSECTION_H
#define CGAL_RAY_2_SEGMENT_2_INTERSECTION_H

#include <CGAL/Ray_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/Line_2.h>
#include <CGAL/Line_2_Line_2_intersection.h>
#include <CGAL/Object.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

template <class K>
class Ray_2_Segment_2_pair {
public:
    enum Intersection_results {NO, POINT, SEGMENT};
    Ray_2_Segment_2_pair() ;
    Ray_2_Segment_2_pair(typename K::Ray_2 const *ray,
			 typename K::Segment_2 const *line);
    ~Ray_2_Segment_2_pair() {}

    Intersection_results intersection_type() const;

    bool                intersection(typename K::Point_2 &result) const;
    bool                intersection(typename K::Segment_2 &result) const;
protected:
    typename K::Ray_2 const *   _ray;
    typename K::Segment_2 const *  _seg;
    mutable bool                    _known;
    mutable Intersection_results    _result;
    mutable typename K::Point_2         _intersection_point, _other_point;
};

template <class K>
inline bool do_intersect(const typename CGAL_WRAP(K)::Ray_2 &p1,
			 const typename CGAL_WRAP(K)::Segment_2 &p2,
			 const K& k)
{
    typedef Ray_2_Segment_2_pair<K> pair_t;
    pair_t pair(&p1, &p2);
    return pair.intersection_type() != pair_t::NO;
}

template <class K>
inline bool do_intersect(const typename CGAL_WRAP(K)::Segment_2 &p2,
			 const typename CGAL_WRAP(K)::Ray_2 &p1,
			 const K& k)
{
  return CGALi::do_intersect(p1, p2, k);
}

template <class K>
Ray_2_Segment_2_pair<K>::Ray_2_Segment_2_pair()
{
    _ray = 0;
    _seg = 0;
    _known = false;
}

template <class K>
Ray_2_Segment_2_pair<K>::Ray_2_Segment_2_pair(
    typename K::Ray_2 const *ray, typename K::Segment_2 const *seg)
{
    _ray = ray;
    _seg = seg;
    _known = false;
}

template <class K>
typename Ray_2_Segment_2_pair<K>::Intersection_results
Ray_2_Segment_2_pair<K>::intersection_type() const
{
    if (_known)
        return _result;
    // The non const this pointer is used to cast away const.
    _known = true;
//    if (!do_overlap(_ray->bbox(), _seg->bbox()))
//        return NO;
    const typename K::Line_2 &l1 = _ray->supporting_line();
    const typename K::Line_2 &l2 = _seg->supporting_line();
    Line_2_Line_2_pair<K> linepair(&l1, &l2);
    switch ( linepair.intersection_type()) {
    case Line_2_Line_2_pair<K>::NO:
        _result = NO;
        return _result;
    case Line_2_Line_2_pair<K>::POINT:
        linepair.intersection(_intersection_point);
        _result = (_ray->collinear_has_on(_intersection_point)
                && _seg->collinear_has_on(_intersection_point) )
            ? POINT :  NO;
        return _result;
    case Line_2_Line_2_pair<K>::LINE: {
        typedef typename K::RT RT;
        const typename K::Point_2 &start1 = _seg->source();
        const typename K::Point_2 &end1 = _seg->target();
        const typename K::Point_2 &start2 = _ray->source();
        const typename K::Point_2 *minpt,  *maxpt;
        typename K::Vector_2 diff1 = end1-start1;
        if (CGAL_NTS abs(diff1.x()) > CGAL_NTS abs(diff1.y())) {
            typedef typename K::FT FT;
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
                    _intersection_point = _seg->source();
                    _other_point = _seg->target();
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
                    _intersection_point = _seg->source();
                    _other_point = _seg->target();
                }
                _result = SEGMENT;
                return _result;
            }
        } else {
            typedef typename K::FT FT;
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
                    _intersection_point = _seg->source();
                    _other_point = _seg->target();
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
                    _intersection_point = _seg->source();
                    _other_point = _seg->target();
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

template <class K>
bool
Ray_2_Segment_2_pair<K>::intersection(typename K::Point_2 &result) const
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
Ray_2_Segment_2_pair<K>::intersection(typename K::Segment_2 &result) const
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
intersection(const typename CGAL_WRAP(K)::Ray_2 &ray, 
	     const typename CGAL_WRAP(K)::Segment_2&seg,
	     const K&)
{
    typedef Ray_2_Segment_2_pair<K> is_t;
    is_t ispair(&ray, &seg);
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
Object
intersection(const typename CGAL_WRAP(K)::Segment_2 &seg,
	     const typename CGAL_WRAP(K)::Ray_2 &ray, 
	     const K& k)
{
  return CGALi::intersection(ray, seg, k);
}


template <class K>
class Segment_2_Ray_2_pair: public Ray_2_Segment_2_pair<K> {
public:
    Segment_2_Ray_2_pair(
            typename K::Segment_2 const *seg,
            typename K::Ray_2 const *ray) :
                                Ray_2_Segment_2_pair<K>(ray, seg) {}
};


} // namespace CGALi

template <class K>
inline bool do_intersect(
    const Segment_2<K> &p1,
    const Ray_2<K> &p2)
{
  return typename K::Do_intersect_2()(p1, p2);
}

template <class K>
inline bool do_intersect(
    const Ray_2<K> &p1,
    const Segment_2<K> &p2)
{
  return typename K::Do_intersect_2()(p2, p1);
}

template <class K>
inline Object
intersection(const Segment_2<K> &seg, const Ray_2<K> &ray)
{
    return typename K::Intersect_2()(ray, seg);
}

template <class K>
inline Object
intersection(const Ray_2<K> &ray, const Segment_2<K> &seg)
{
    return typename K::Intersect_2()(ray, seg);
}
CGAL_END_NAMESPACE

#endif
