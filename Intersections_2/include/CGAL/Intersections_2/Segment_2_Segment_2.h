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


#ifndef CGAL_INTERSECTIONS_2_SEGMENT_2_SEGMENT_2_H
#define CGAL_INTERSECTIONS_2_SEGMENT_2_SEGMENT_2_H

#include <CGAL/Segment_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/kernel_assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/Intersections_2/Line_2_Line_2.h>
#include <CGAL/Uncertain.h>
#include <CGAL/Intersection_traits_2.h>

namespace CGAL {

namespace Intersections {

namespace internal {

template <class K>
inline bool
do_intersect(const typename K::Segment_2 &seg1, const typename K::Segment_2 &seg2);




template <class K>
bool seg_seg_do_intersect_crossing(
        const typename K::Point_2  &p1, const typename K::Point_2 &p2,
        const typename K::Point_2 &p3, const typename K::Point_2 &p4,
        const K& k)
{
    switch (make_certain(k.orientation_2_object()(p1,p2,p3))) {
    case LEFT_TURN:
      return ! (k.orientation_2_object()(p3,p4,p2) == RIGHT_TURN); //   right_turn(p3,p4,p2);
    case RIGHT_TURN:
        return ! (k.orientation_2_object()(p3,p4,p2) == LEFT_TURN); //left_turn(p3,p4,p2);
    case COLLINEAR:
        return true;
    }
    CGAL_kernel_assertion(false);
    return false;
}


template <class K>
bool seg_seg_do_intersect_contained(
        const typename K::Point_2  &p1, const typename K::Point_2 &p2,
        const typename K::Point_2 &p3, const typename K::Point_2 &p4,
        const K& k)
{
    switch (make_certain(k.orientation_2_object()(p1,p2,p3))) {
    case LEFT_TURN:
      return ! (k.orientation_2_object()(p1,p2,p4) == LEFT_TURN); // left_turn(p1,p2,p4);
    case RIGHT_TURN:
        return ! (k.orientation_2_object()(p1,p2,p4) == RIGHT_TURN); // right_turn(p1,p2,p4);
    case COLLINEAR:
        return true;
    }
    CGAL_kernel_assertion(false);
    return false;
}


template <class K>
bool
do_intersect(const typename K::Segment_2 &seg1,
             const typename K::Segment_2 &seg2,
             const K& k)
{
    typename K::Point_2 const & A1 = seg1.source();
    typename K::Point_2 const & A2 = seg1.target();
    typename K::Point_2 const & B1 = seg2.source();
    typename K::Point_2 const & B2 = seg2.target();
    typename K::Less_xy_2 less_xy;
    typename K::Compare_xy_2 compare_xy;

    if (less_xy(A1,A2)) {
        if (less_xy(B1,B2)) {
            if (less_xy(A2,B1)
             || less_xy(B2,A1))
                return false;
        } else {
            if (less_xy(A2,B2)
             || less_xy(B1,A1))
                return false;
        }
    } else {
        if (less_xy(B1,B2)) {
            if (less_xy(A1,B1)
             || less_xy(B2,A2))
                return false;
        } else {
            if (less_xy(A1,B2)
             || less_xy(B1,A2))
                return false;
        }
    }
    if (less_xy(A1,A2)) {
        if (less_xy(B1,B2)) {
            switch(make_certain(compare_xy(A1,B1))) {
            case SMALLER:
                switch(make_certain(compare_xy(A2,B1))) {
                case SMALLER:
                    return false;
                case EQUAL:
                    return true;
                default: // LARGER
                    switch(make_certain(compare_xy(A2,B2))) {
                    case SMALLER:
                        return seg_seg_do_intersect_crossing(A1,A2,B1,B2, k);
                    case EQUAL:
                        return true;
                    default: // LARGER
                        return seg_seg_do_intersect_contained(A1,A2,B1,B2, k);
                    }
                }
            case EQUAL:
                return true;
            default: // LARGER
                switch(make_certain(compare_xy(B2,A1))) {
                case SMALLER:
                    return false;
                case EQUAL:
                    return true;
                default: // LARGER
                    switch(make_certain(compare_xy(B2,A2))) {
                    case SMALLER:
                        return seg_seg_do_intersect_crossing(B1,B2,A1,A2, k);
                    case EQUAL:
                        return true;
                    default: // LARGER
                        return seg_seg_do_intersect_contained(B1,B2,A1,A2, k);
                    }
                }
            }
        } else {
            switch(make_certain(compare_xy(A1,B2))) {
            case SMALLER:
                switch(make_certain(compare_xy(A2,B2))) {
                case SMALLER:
                    return false;
                case EQUAL:
                    return true;
                default: // LARGER
                    switch(make_certain(compare_xy(A2,B1))) {
                    case SMALLER:
                        return seg_seg_do_intersect_crossing(A1,A2,B2,B1, k);
                    case EQUAL:
                        return true;
                    default: // LARGER
                        return seg_seg_do_intersect_contained(A1,A2,B2,B1, k);
                    }
                }
            case EQUAL:
                return true;
            default: // LARGER
                switch(make_certain(compare_xy(B1,A1))) {
                case SMALLER:
                    return false;
                case EQUAL:
                    return true;
                default: // LARGER
                    switch(make_certain(compare_xy(B1,A2))) {
                    case SMALLER:
                        return seg_seg_do_intersect_crossing(B2,B1,A1,A2, k);
                    case EQUAL:
                        return true;
                    default: // LARGER
                        return seg_seg_do_intersect_contained(B2,B1,A1,A2, k);
                    }
                }
            }
        }
    } else {
        if (less_xy(B1,B2)) {
            switch(make_certain(compare_xy(A2,B1))) {
            case SMALLER:
                switch(make_certain(compare_xy(A1,B1))) {
                case SMALLER:
                    return false;
                case EQUAL:
                    return true;
                default: // LARGER
                    switch(make_certain(compare_xy(A1,B2))) {
                    case SMALLER:
                        return seg_seg_do_intersect_crossing(A2,A1,B1,B2, k);
                    case EQUAL:
                        return true;
                    default: // LARGER
                        return seg_seg_do_intersect_contained(A2,A1,B1,B2, k);
                    }
                }
            case EQUAL:
                return true;
            default: // LARGER
                switch(make_certain(compare_xy(B2,A2))) {
                case SMALLER:
                    return false;
                case EQUAL:
                    return true;
                default: // LARGER
                    switch(make_certain(compare_xy(B2,A1))) {
                    case SMALLER:
                        return seg_seg_do_intersect_crossing(B1,B2,A2,A1, k);
                    case EQUAL:
                        return true;
                    default: // LARGER
                        return seg_seg_do_intersect_contained(B1,B2,A2,A1, k);
                    }
                }
            }
        } else {
            switch(make_certain(compare_xy(A2,B2))) {
            case SMALLER:
                switch(make_certain(compare_xy(A1,B2))) {
                case SMALLER:
                    return false;
                case EQUAL:
                    return true;
                default: // LARGER
                    switch(make_certain(compare_xy(A1,B1))) {
                    case SMALLER:
                        return seg_seg_do_intersect_crossing(A2,A1,B2,B1, k);
                    case EQUAL:
                        return true;
                    default: // LARGER
                        return seg_seg_do_intersect_contained(A2,A1,B2,B1, k);
                    }
                }
            case EQUAL:
                return true;
            default: // LARGER
                switch(make_certain(compare_xy(B1,A2))) {
                case SMALLER:
                    return false;
                case EQUAL:
                    return true;
                default: // LARGER
                    switch(make_certain(compare_xy(B1,A1))) {
                    case SMALLER:
                        return seg_seg_do_intersect_crossing(B2,B1,A2,A1, k);
                    case EQUAL:
                        return true;
                    default: // LARGER
                        return seg_seg_do_intersect_contained(B2,B1,A2,A1, k);
                    }
                }
            }
        }
    }
    CGAL_kernel_assertion(false);
    return false;
}


template <class K>
class Segment_2_Segment_2_pair {
public:
    enum Intersection_results {NO_INTERSECTION, POINT, SEGMENT};
    Segment_2_Segment_2_pair(typename K::Segment_2 const *seg1,
                            typename K::Segment_2 const *seg2)
            : _seg1(seg1), _seg2(seg2), _known(false) {}

    Intersection_results intersection_type() const;

    typename K::Point_2    intersection_point() const;
    typename K::Segment_2  intersection_segment() const;
protected:
    typename K::Segment_2 const*   _seg1;
    typename K::Segment_2 const *  _seg2;
    mutable bool                       _known;
    mutable Intersection_results       _result;
    mutable typename K::Point_2            _intersection_point, _other_point;
};

template <class K>
typename Segment_2_Segment_2_pair<K>::Intersection_results
Segment_2_Segment_2_pair<K>::intersection_type() const
{
  typename K::Construct_vector_2 construct_vector;
    if (_known)
        return _result;
    _known = true;
    if (!internal::do_intersect(*_seg1, *_seg2, K())) {
        _result = NO_INTERSECTION;
        return _result;
    }
    typename K::Line_2 const &l1 = _seg1->supporting_line();
    typename K::Line_2 const &l2 = _seg2->supporting_line();
    Line_2_Line_2_pair<K> linepair(&l1, &l2);
    switch ( linepair.intersection_type()) {
    case Line_2_Line_2_pair<K>::NO_INTERSECTION:
    default:
        _result = NO_INTERSECTION;
        break;
    case Line_2_Line_2_pair<K>::POINT:
        _intersection_point = linepair.intersection_point();
        _result = POINT;
        break;
    case Line_2_Line_2_pair<K>::LINE:
        {
        //typedef typename K::RT RT;
        typename K::Point_2 const &start1 = _seg1->source();
        typename K::Point_2 const &end1   = _seg1->target();
        typename K::Point_2 const &start2 = _seg2->source();
        typename K::Point_2 const &end2   = _seg2->target();
        typename K::Vector_2 diff1 = construct_vector(start1, end1);
        typename K::Point_2 const *minpt;
        typename K::Point_2 const *maxpt;
        if (CGAL_NTS abs(diff1.x()) > CGAL_NTS abs(diff1.y())) {
            if (start1.x() < end1.x()) {
                minpt = &start1;
                maxpt = &end1;
            } else {
                minpt = &end1;
                maxpt = &start1;
            }
            if (start2.x() < end2.x()) {
                if (start2.x() > minpt->x()) {
                    minpt = &start2;
                }
                if (end2.x() < maxpt->x()) {
                    maxpt = &end2;
                }
            } else {
                if (end2.x() > minpt->x()) {
                    minpt = &end2;
                }
                if (start2.x() < maxpt->x()) {
                    maxpt = &start2;
                }
            }
            if (maxpt->x() < minpt->x()) {
                _result = NO_INTERSECTION;
                return _result;
            }
            if (maxpt->x() == minpt->x()) {
                _intersection_point = *minpt;
                _result = POINT;
                return _result;
            }
            _intersection_point = *minpt;
            _other_point = *maxpt;
            _result = SEGMENT;
            return _result;
        } else {
            if (start1.y() < end1.y()) {
                minpt = &start1;
                maxpt = &end1;
            } else {
                minpt = &end1;
                maxpt = &start1;
            }
            if (start2.y() < end2.y()) {
                if (start2.y() > minpt->y()) {
                    minpt = &start2;
                }
                if (end2.y() < maxpt->y()) {
                    maxpt = &end2;
                }
            } else {
                if (end2.y() > minpt->y()) {
                    minpt = &end2;
                }
                if (start2.y() < maxpt->y()) {
                    maxpt = &start2;
                }
            }
            if (maxpt->y() < minpt->y()) {
                _result = NO_INTERSECTION;
                return _result;
            }
            if (maxpt->y() == minpt->y()) {
                _intersection_point = *minpt;
                _result = POINT;
                return _result;
            }
            _intersection_point = *minpt;
            _other_point = *maxpt;
            _result = SEGMENT;
            return _result;
        }
        }
    }
    return _result;
}


template <class K>
typename K::Point_2
Segment_2_Segment_2_pair<K>::intersection_point() const
{
    if (!_known)
        intersection_type();
    CGAL_kernel_assertion(_result == POINT);
    return _intersection_point;
}

template <class K>
typename K::Segment_2
Segment_2_Segment_2_pair<K>::intersection_segment() const
{
  typedef typename K::Segment_2 Segment_2;
    if (!_known)
        intersection_type();
    CGAL_kernel_assertion(_result == SEGMENT);
    return Segment_2(_intersection_point, _other_point);
}


template <class K>
typename CGAL::Intersection_traits
<K, typename K::Segment_2, typename K::Segment_2>::result_type
intersection(const typename K::Segment_2 &seg1,
             const typename K::Segment_2 &seg2,
             const K&)
{
    typedef Segment_2_Segment_2_pair<K> is_t;
    is_t ispair(&seg1, &seg2);
    switch (ispair.intersection_type()) {
    case is_t::NO_INTERSECTION:
    default:
        return intersection_return<typename K::Intersect_2, typename K::Segment_2, typename K::Segment_2>();
    case is_t::POINT:
        return intersection_return<typename K::Intersect_2, typename K::Segment_2, typename K::Segment_2>(ispair.intersection_point());
    case is_t::SEGMENT:
        return intersection_return<typename K::Intersect_2, typename K::Segment_2, typename K::Segment_2>(ispair.intersection_segment());
    }
}

} // namespace internal
} // namespace Intersections

CGAL_INTERSECTION_FUNCTION_SELF(Segment_2, 2)
CGAL_DO_INTERSECT_FUNCTION_SELF(Segment_2, 2)

} //namespace CGAL

#endif
