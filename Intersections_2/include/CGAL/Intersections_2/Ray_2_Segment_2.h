// Copyright (c) 2000
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Geert-Jan Giezeman


#ifndef CGAL_INTERSECTIONS_2_RAY_2_SEGMENT_2_H
#define CGAL_INTERSECTIONS_2_RAY_2_SEGMENT_2_H

#include <CGAL/Ray_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/kernel_assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/Intersections_2/Line_2_Line_2.h>
#include <CGAL/Intersection_traits_2.h>


namespace CGAL {

namespace Intersections {

namespace internal {

template <class K>
class Ray_2_Segment_2_pair {
public:
    enum Intersection_results {NO_INTERSECTION, POINT, SEGMENT};
    Ray_2_Segment_2_pair(typename K::Ray_2 const *ray,
                         typename K::Segment_2 const *seg)
            : _ray(ray), _seg(seg), _known(false) {}

    Intersection_results intersection_type() const;

    typename K::Point_2    intersection_point() const;
    typename K::Segment_2  intersection_segment() const;
protected:
    typename K::Ray_2 const *   _ray;
    typename K::Segment_2 const *  _seg;
    mutable bool                    _known;
    mutable Intersection_results    _result;
    mutable typename K::Point_2         _intersection_point, _other_point;
};

template <class K>
inline bool do_intersect(const typename K::Ray_2 &p1,
                         const typename K::Segment_2 &p2,
                         const K&)
{
    typedef Ray_2_Segment_2_pair<K> pair_t;
    pair_t pair(&p1, &p2);
    return pair.intersection_type() != pair_t::NO_INTERSECTION;
}

template <class K>
inline bool do_intersect(const typename K::Segment_2 &p2,
                         const typename K::Ray_2 &p1,
                         const K& k)
{
  return internal::do_intersect(p1, p2, k);
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
//        return NO_INTERSECTION;
    const typename K::Line_2 &l1 = _ray->supporting_line();
    const typename K::Line_2 &l2 = _seg->supporting_line();
    Line_2_Line_2_pair<K> linepair(&l1, &l2);
    switch ( linepair.intersection_type()) {
    case Line_2_Line_2_pair<K>::NO_INTERSECTION:
        _result = NO_INTERSECTION;
        return _result;
    case Line_2_Line_2_pair<K>::POINT:
        _intersection_point = linepair.intersection_point();
        _result = (_ray->collinear_has_on(_intersection_point)
                && _seg->collinear_has_on(_intersection_point) )
            ? POINT :  NO_INTERSECTION;
        return _result;
    case Line_2_Line_2_pair<K>::LINE: {
        //typedef typename K::RT RT;
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
                    _result = NO_INTERSECTION;
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
                    _result = NO_INTERSECTION;
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
                    _result = NO_INTERSECTION;
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
                    _result = NO_INTERSECTION;
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
typename K::Point_2
Ray_2_Segment_2_pair<K>::intersection_point() const
{
    if (!_known)
        intersection_type();
    CGAL_kernel_assertion(_result == POINT);
    return _intersection_point;
}

template <class K>
typename K::Segment_2
Ray_2_Segment_2_pair<K>::intersection_segment() const
{
  typedef typename K::Segment_2 Segment_2;
    if (!_known)
        intersection_type();
    CGAL_kernel_assertion(_result == SEGMENT);
    return Segment_2(_intersection_point, _other_point);
}



template <class K>
typename CGAL::Intersection_traits
<K, typename K::Ray_2, typename K::Segment_2>::result_type
intersection(const typename K::Ray_2 &ray,
             const typename K::Segment_2&seg,
             const K&)
{
    typedef Ray_2_Segment_2_pair<K> is_t;
    is_t ispair(&ray, &seg);
    switch (ispair.intersection_type()) {
    case is_t::NO_INTERSECTION:
    default:
        return intersection_return<typename K::Intersect_2, typename K::Ray_2, typename K::Segment_2>();
    case is_t::POINT:
        return intersection_return<typename K::Intersect_2, typename K::Ray_2, typename K::Segment_2>(ispair.intersection_point());
    case is_t::SEGMENT:
        return intersection_return<typename K::Intersect_2, typename K::Ray_2, typename K::Segment_2>(ispair.intersection_segment());
    }
}


template <class K>
typename CGAL::Intersection_traits
<K, typename K::Ray_2, typename K::Segment_2>::result_type
intersection(const typename K::Segment_2 &seg,
             const typename K::Ray_2 &ray,
             const K& k)
{
  return internal::intersection(ray, seg, k);
}

} // namespace internal
} // namespace Intersections

CGAL_INTERSECTION_FUNCTION(Ray_2, Segment_2, 2)
CGAL_DO_INTERSECT_FUNCTION(Ray_2, Segment_2, 2)

} //namespace CGAL

#endif
