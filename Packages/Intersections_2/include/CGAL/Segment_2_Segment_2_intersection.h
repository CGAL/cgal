
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


#ifndef CGAL_SEGMENT_2_SEGMENT_2_INTERSECTION_H
#define CGAL_SEGMENT_2_SEGMENT_2_INTERSECTION_H

#include <CGAL/Segment_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>

#include <cassert>
#include <CGAL/predicates_on_points_2.h>



#include <CGAL/Line_2.h>
#include <CGAL/Line_2_Line_2_intersection.h>

#include <CGAL/Object.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

template <class K>
class Segment_2_Segment_2_pair {
public:
    enum Intersection_results {NO, POINT, SEGMENT};
    Segment_2_Segment_2_pair() ;
    Segment_2_Segment_2_pair(typename K::Segment_2 const *seg1,
                            typename K::Segment_2 const *seg2);
    ~Segment_2_Segment_2_pair() {}

    Intersection_results intersection_type() const;

    bool                intersection(typename K::Point_2 &result) const;
    bool                intersection(typename K::Segment_2 &result) const;
protected:
    typename K::Segment_2 const*   _seg1;
    typename K::Segment_2 const *  _seg2;
    mutable bool                       _known;
    mutable Intersection_results       _result;
    mutable typename K::Point_2            _intersection_point, _other_point;
};

template <class K>
inline bool
do_intersect(const typename K::Segment_2 &seg1, const typename K::Segment_2 &seg2);







template <class K>
bool seg_seg_do_intersect_crossing(
        const typename K::Point_2  &p1, const typename K::Point_2 &p2, 
	const typename K::Point_2 &p3, const typename K::Point_2 &p4,
	const K& k)
{
  typename K::Orientation_2 orientation;
    switch (orientation(p1,p2,p3)) {
    case LEFT_TURN:
      return ! (orientation(p3,p4,p2) == RIGHT_TURN); //   right_turn(p3,p4,p2);
    case RIGHT_TURN:
        return ! (orientation(p3,p4,p2) == LEFT_TURN); //left_turn(p3,p4,p2);
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
  typename K::Orientation_2 orientation;
    switch (orientation(p1,p2,p3)) {
    case LEFT_TURN:
      return ! (orientation(p1,p2,p4) == LEFT_TURN); // left_turn(p1,p2,p4);
    case RIGHT_TURN:
        return ! (orientation(p1,p2,p4) == RIGHT_TURN); // right_turn(p1,p2,p4);
    case COLLINEAR:
        return true;
    }
    CGAL_kernel_assertion(false);
    return false;
}


template <class K>
bool
do_intersect(const typename CGAL_WRAP(K)::Segment_2 &seg1, 
	     const typename CGAL_WRAP(K)::Segment_2 &seg2,
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
            switch(compare_xy(A1,B1)) {
            case SMALLER:
                switch(compare_xy(A2,B1)) {
                case SMALLER:
                    return false;
                case EQUAL:
                    return true;
                case LARGER:
                    switch(compare_xy(A2,B2)) {
                    case SMALLER:
                        return seg_seg_do_intersect_crossing(A1,A2,B1,B2, k);
                    case EQUAL:
                        return true;
                    case LARGER:
                        return seg_seg_do_intersect_contained(A1,A2,B1,B2, k);
                    }
                }
            case EQUAL:
                return true;
            case LARGER:
                switch(compare_xy(B2,A1)) {
                case SMALLER:
                    return false;
                case EQUAL:
                    return true;
                case LARGER:
                    switch(compare_xy(B2,A2)) {
                    case SMALLER:
                        return seg_seg_do_intersect_crossing(B1,B2,A1,A2, k);
                    case EQUAL:
                        return true;
                    case LARGER:
                        return seg_seg_do_intersect_contained(B1,B2,A1,A2, k);
                    }
                }
            }
        } else {
            switch(compare_xy(A1,B2)) {
            case SMALLER:
                switch(compare_xy(A2,B2)) {
                case SMALLER:
                    return false;
                case EQUAL:
                    return true;
                case LARGER:
                    switch(compare_xy(A2,B1)) {
                    case SMALLER:
                        return seg_seg_do_intersect_crossing(A1,A2,B2,B1, k);
                    case EQUAL:
                        return true;
                    case LARGER:
                        return seg_seg_do_intersect_contained(A1,A2,B2,B1, k);
                    }
                }
            case EQUAL:
                return true;
            case LARGER:
                switch(compare_xy(B1,A1)) {
                case SMALLER:
                    return false;
                case EQUAL:
                    return true;
                case LARGER:
                    switch(compare_xy(B1,A2)) {
                    case SMALLER:
                        return seg_seg_do_intersect_crossing(B2,B1,A1,A2, k);
                    case EQUAL:
                        return true;
                    case LARGER:
                        return seg_seg_do_intersect_contained(B2,B1,A1,A2, k);
                    }
                }
            }
        }
    } else {
        if (less_xy(B1,B2)) {
            switch(compare_xy(A2,B1)) {
            case SMALLER:
                switch(compare_xy(A1,B1)) {
                case SMALLER:
                    return false;
                case EQUAL:
                    return true;
                case LARGER:
                    switch(compare_xy(A1,B2)) {
                    case SMALLER:
                        return seg_seg_do_intersect_crossing(A2,A1,B1,B2, k);
                    case EQUAL:
                        return true;
                    case LARGER:
                        return seg_seg_do_intersect_contained(A2,A1,B1,B2, k);
                    }
                }
            case EQUAL:
                return true;
            case LARGER:
                switch(compare_xy(B2,A2)) {
                case SMALLER:
                    return false;
                case EQUAL:
                    return true;
                case LARGER:
                    switch(compare_xy(B2,A1)) {
                    case SMALLER:
                        return seg_seg_do_intersect_crossing(B1,B2,A2,A1, k);
                    case EQUAL:
                        return true;
                    case LARGER:
                        return seg_seg_do_intersect_contained(B1,B2,A2,A1, k);
                    }
                }
            }
        } else {
            switch(compare_xy(A2,B2)) {
            case SMALLER:
                switch(compare_xy(A1,B2)) {
                case SMALLER:
                    return false;
                case EQUAL:
                    return true;
                case LARGER:
                    switch(compare_xy(A1,B1)) {
                    case SMALLER:
                        return seg_seg_do_intersect_crossing(A2,A1,B2,B1, k);
                    case EQUAL:
                        return true;
                    case LARGER:
                        return seg_seg_do_intersect_contained(A2,A1,B2,B1, k);
                    }
                }
            case EQUAL:
                return true;
            case LARGER:
                switch(compare_xy(B1,A2)) {
                case SMALLER:
                    return false;
                case EQUAL:
                    return true;
                case LARGER:
                    switch(compare_xy(B1,A1)) {
                    case SMALLER:
                        return seg_seg_do_intersect_crossing(B2,B1,A2,A1, k);
                    case EQUAL:
                        return true;
                    case LARGER:
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
Segment_2_Segment_2_pair<K>::Segment_2_Segment_2_pair()
{
    _seg1 = 0;
    _seg2 = 0;
    _known = false;
}

template <class K>
Segment_2_Segment_2_pair<K>::Segment_2_Segment_2_pair(
    typename K::Segment_2 const *seg1, typename K::Segment_2 const *seg2)
{
    _seg1 = seg1;
    _seg2 = seg2;
    _known = false;
}

template <class K>
typename Segment_2_Segment_2_pair<K>::Intersection_results
Segment_2_Segment_2_pair<K>::intersection_type() const
{
  typename K::Construct_vector_2 construct_vector;
    if (_known)
        return _result;
    _known = true;
    if (!do_intersect(*_seg1, *_seg2, K())) {
        _result = NO;
        return _result;
    }
    typename K::Line_2 const &l1 = _seg1->supporting_line();
    typename K::Line_2 const &l2 = _seg2->supporting_line();
    Line_2_Line_2_pair<K> linepair(&l1, &l2);
    switch ( linepair.intersection_type()) {
    case Line_2_Line_2_pair<K>::NO:
        _result = NO;
        break;
    case Line_2_Line_2_pair<K>::POINT:
        linepair.intersection(_intersection_point);
        _result = POINT;
        break;
    case Line_2_Line_2_pair<K>::LINE:
        {
        typedef typename K::RT RT;
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
                _result = NO;
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
                _result = NO;
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
bool
Segment_2_Segment_2_pair<K>::intersection(typename K::Point_2 &result) const
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
Segment_2_Segment_2_pair<K>::intersection(typename K::Segment_2 &result) const
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
intersection(const typename CGAL_WRAP(K)::Segment_2 &seg1, 
	     const typename CGAL_WRAP(K)::Segment_2 &seg2,
	     const K&)
{
    typedef Segment_2_Segment_2_pair<K> is_t;
    is_t ispair(&seg1, &seg2);
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

} // namespace CGALi

template <class K>
inline
bool
do_intersect(const Segment_2<K> &seg1, 
	     const Segment_2<K> &seg2)
{
  typedef typename K::Do_intersect_2 Do_intersect;
  return Do_intersect()(seg1, seg2);
}


template <class K>
Object
intersection(const Segment_2<K> &seg1, 
	     const Segment_2<K> &seg2)
{
  typedef typename K::Intersect_2 Intersect;
  return Intersect()(seg1, seg2);
}

CGAL_END_NAMESPACE

#endif
