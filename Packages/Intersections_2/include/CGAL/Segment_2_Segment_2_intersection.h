
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
// file          : include/CGAL/Segment_2_Segment_2_intersection.h
// source        : intersection_2_1.fw
// author(s)     : Geert-Jan Giezeman
//
// coordinator   : Saarbruecken
//
// ============================================================================


#ifndef CGAL_SEGMENT_2_SEGMENT_2_INTERSECTION_H
#define CGAL_SEGMENT_2_SEGMENT_2_INTERSECTION_H

#include <CGAL/Segment_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>

CGAL_BEGIN_NAMESPACE

template <class R>
class Segment_2_Segment_2_pair {
public:
    enum Intersection_results {NO, POINT, SEGMENT};
    Segment_2_Segment_2_pair() ;
    Segment_2_Segment_2_pair(Segment_2<R> const *seg1,
                            Segment_2<R> const *seg2);
    ~Segment_2_Segment_2_pair() {}

#ifndef CGAL_CFG_RETURN_TYPE_BUG_2
    Intersection_results intersection_type() const;
#else
    Intersection_results intersection_type() const
{
    if (_known)
        return _result;
    _known = true;
    if (!do_intersect(*_seg1, *_seg2)) {
        _result = NO;
        return _result;
    }
    Line_2<R> const &l1 = _seg1->supporting_line();
    Line_2<R> const &l2 = _seg2->supporting_line();
    Line_2_Line_2_pair<R> linepair(&l1, &l2);
    switch ( linepair.intersection_type()) {
    case Line_2_Line_2_pair<R>::NO:
        _result = NO;
        break;
    case Line_2_Line_2_pair<R>::POINT:
        linepair.intersection(_intersection_point);
        _result = POINT;
        break;
    case Line_2_Line_2_pair<R>::LINE:
        {
        typedef typename R::RT RT;
        Point_2<R> const &start1 = _seg1->start();
        Point_2<R> const &end1   = _seg1->end();
        Point_2<R> const &start2 = _seg2->start();
        Point_2<R> const &end2   = _seg2->end();
        Vector_2<R> diff1 = end1-start1;
        Point_2<R> const *minpt;
        Point_2<R> const *maxpt;
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

#endif // CGAL_CFG_RETURN_TYPE_BUG_2
    bool                intersection(Point_2<R> &result) const;
    bool                intersection(Segment_2<R> &result) const;
protected:
    Segment_2<R> const*   _seg1;
    Segment_2<R> const *  _seg2;
    mutable bool                       _known;
    mutable Intersection_results       _result;
    mutable Point_2<R>            _intersection_point, _other_point;
};

template <class R>
inline bool
do_intersect(const Segment_2<R> &seg1, const Segment_2<R> &seg2);


CGAL_END_NAMESPACE


#include <cassert>
#include <CGAL/predicates_on_points_2.h>

namespace CGAL {

template <class PT>
bool seg_seg_do_intersect_crossing(
        PT const &p1, PT const &p2, PT const &p3, PT const &p4)
{
    switch (orientation(p1,p2,p3)) {
    case LEFTTURN:
        return !rightturn(p3,p4,p2);
    case RIGHTTURN:
        return !leftturn(p3,p4,p2);
    case COLLINEAR:
        return true;
    }
    assert(false);
    return false;
}


template <class PT>
bool seg_seg_do_intersect_contained(
        PT const &p1, PT const &p2, PT const &p3, PT const &p4)
{
    switch (orientation(p1,p2,p3)) {
    case LEFTTURN:
        return !leftturn(p1,p2,p4);
    case RIGHTTURN:
        return !rightturn(p1,p2,p4);
    case COLLINEAR:
        return true;
    }
    assert(false);
    return false;
}


template <class R>
bool
do_intersect(const Segment_2<R> &seg1, const Segment_2<R> &seg2)
{
    typename R::Point_2 const & A1 = seg1.source();
    typename R::Point_2 const & A2 = seg1.target();
    typename R::Point_2 const & B1 = seg2.source();
    typename R::Point_2 const & B2 = seg2.target();
    if (lexicographically_yx_smaller(A1,A2)) {
        if (lexicographically_yx_smaller(B1,B2)) {
            if (lexicographically_yx_smaller(A2,B1)
             || lexicographically_yx_smaller(B2,A1))
                return false;
        } else {
            if (lexicographically_yx_smaller(A2,B2)
             || lexicographically_yx_smaller(B1,A1))
                return false;
        }
    } else {
        if (lexicographically_yx_smaller(B1,B2)) {
            if (lexicographically_yx_smaller(A1,B1)
             || lexicographically_yx_smaller(B2,A2))
                return false;
        } else {
            if (lexicographically_yx_smaller(A1,B2)
             || lexicographically_yx_smaller(B1,A2))
                return false;
        }
    }
    if (lexicographically_xy_smaller(A1,A2)) {
        if (lexicographically_xy_smaller(B1,B2)) {
            switch(compare_lexicographically_xy(A1,B1)) {
            case SMALLER:
                switch(compare_lexicographically_xy(A2,B1)) {
                case SMALLER:
                    return false;
                case EQUAL:
                    return true;
                case LARGER:
                    switch(compare_lexicographically_xy(A2,B2)) {
                    case SMALLER:
                        return seg_seg_do_intersect_crossing(A1,A2,B1,B2);
                    case EQUAL:
                        return true;
                    case LARGER:
                        return seg_seg_do_intersect_contained(A1,A2,B1,B2);
                    }
                }
            case EQUAL:
                return true;
            case LARGER:
                switch(compare_lexicographically_xy(B2,A1)) {
                case SMALLER:
                    return false;
                case EQUAL:
                    return true;
                case LARGER:
                    switch(compare_lexicographically_xy(B2,A2)) {
                    case SMALLER:
                        return seg_seg_do_intersect_crossing(B1,B2,A1,A2);
                    case EQUAL:
                        return true;
                    case LARGER:
                        return seg_seg_do_intersect_contained(B1,B2,A1,A2);
                    }
                }
            }
        } else {
            switch(compare_lexicographically_xy(A1,B2)) {
            case SMALLER:
                switch(compare_lexicographically_xy(A2,B2)) {
                case SMALLER:
                    return false;
                case EQUAL:
                    return true;
                case LARGER:
                    switch(compare_lexicographically_xy(A2,B1)) {
                    case SMALLER:
                        return seg_seg_do_intersect_crossing(A1,A2,B2,B1);
                    case EQUAL:
                        return true;
                    case LARGER:
                        return seg_seg_do_intersect_contained(A1,A2,B2,B1);
                    }
                }
            case EQUAL:
                return true;
            case LARGER:
                switch(compare_lexicographically_xy(B1,A1)) {
                case SMALLER:
                    return false;
                case EQUAL:
                    return true;
                case LARGER:
                    switch(compare_lexicographically_xy(B1,A2)) {
                    case SMALLER:
                        return seg_seg_do_intersect_crossing(B2,B1,A1,A2);
                    case EQUAL:
                        return true;
                    case LARGER:
                        return seg_seg_do_intersect_contained(B2,B1,A1,A2);
                    }
                }
            }
        }
    } else {
        if (lexicographically_xy_smaller(B1,B2)) {
            switch(compare_lexicographically_xy(A2,B1)) {
            case SMALLER:
                switch(compare_lexicographically_xy(A1,B1)) {
                case SMALLER:
                    return false;
                case EQUAL:
                    return true;
                case LARGER:
                    switch(compare_lexicographically_xy(A1,B2)) {
                    case SMALLER:
                        return seg_seg_do_intersect_crossing(A2,A1,B1,B2);
                    case EQUAL:
                        return true;
                    case LARGER:
                        return seg_seg_do_intersect_contained(A2,A1,B1,B2);
                    }
                }
            case EQUAL:
                return true;
            case LARGER:
                switch(compare_lexicographically_xy(B2,A2)) {
                case SMALLER:
                    return false;
                case EQUAL:
                    return true;
                case LARGER:
                    switch(compare_lexicographically_xy(B2,A1)) {
                    case SMALLER:
                        return seg_seg_do_intersect_crossing(B1,B2,A2,A1);
                    case EQUAL:
                        return true;
                    case LARGER:
                        return seg_seg_do_intersect_contained(B1,B2,A2,A1);
                    }
                }
            }
        } else {
            switch(compare_lexicographically_xy(A2,B2)) {
            case SMALLER:
                switch(compare_lexicographically_xy(A1,B2)) {
                case SMALLER:
                    return false;
                case EQUAL:
                    return true;
                case LARGER:
                    switch(compare_lexicographically_xy(A1,B1)) {
                    case SMALLER:
                        return seg_seg_do_intersect_crossing(A2,A1,B2,B1);
                    case EQUAL:
                        return true;
                    case LARGER:
                        return seg_seg_do_intersect_contained(A2,A1,B2,B1);
                    }
                }
            case EQUAL:
                return true;
            case LARGER:
                switch(compare_lexicographically_xy(B1,A2)) {
                case SMALLER:
                    return false;
                case EQUAL:
                    return true;
                case LARGER:
                    switch(compare_lexicographically_xy(B1,A1)) {
                    case SMALLER:
                        return seg_seg_do_intersect_crossing(B2,B1,A2,A1);
                    case EQUAL:
                        return true;
                    case LARGER:
                        return seg_seg_do_intersect_contained(B2,B1,A2,A1);
                    }
                }
            }
        }
    }
    assert(false);
    return false;
}

} // end namespace CGAL




#include <CGAL/Line_2.h>
#include <CGAL/Line_2_Line_2_intersection.h>

CGAL_BEGIN_NAMESPACE

template <class R>
Segment_2_Segment_2_pair<R>::Segment_2_Segment_2_pair()
{
    _seg1 = 0;
    _seg2 = 0;
    _known = false;
}

template <class R>
Segment_2_Segment_2_pair<R>::Segment_2_Segment_2_pair(
    Segment_2<R> const *seg1, Segment_2<R> const *seg2)
{
    _seg1 = seg1;
    _seg2 = seg2;
    _known = false;
}

#ifndef CGAL_CFG_RETURN_TYPE_BUG_2
template <class R>
Segment_2_Segment_2_pair<R>::Intersection_results
Segment_2_Segment_2_pair<R>::intersection_type() const
{
    if (_known)
        return _result;
    _known = true;
    if (!do_intersect(*_seg1, *_seg2)) {
        _result = NO;
        return _result;
    }
    Line_2<R> const &l1 = _seg1->supporting_line();
    Line_2<R> const &l2 = _seg2->supporting_line();
    Line_2_Line_2_pair<R> linepair(&l1, &l2);
    switch ( linepair.intersection_type()) {
    case Line_2_Line_2_pair<R>::NO:
        _result = NO;
        break;
    case Line_2_Line_2_pair<R>::POINT:
        linepair.intersection(_intersection_point);
        _result = POINT;
        break;
    case Line_2_Line_2_pair<R>::LINE:
        {
        typedef typename R::RT RT;
        Point_2<R> const &start1 = _seg1->start();
        Point_2<R> const &end1   = _seg1->end();
        Point_2<R> const &start2 = _seg2->start();
        Point_2<R> const &end2   = _seg2->end();
        Vector_2<R> diff1 = end1-start1;
        Point_2<R> const *minpt;
        Point_2<R> const *maxpt;
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

#endif // CGAL_CFG_RETURN_TYPE_BUG_2

template <class R>
bool
Segment_2_Segment_2_pair<R>::intersection(Point_2<R> &result) const
{
    if (!_known)
        intersection_type();
    if (_result != POINT)
        return false;
    result = _intersection_point;
    return true;
}

template <class R>
bool
Segment_2_Segment_2_pair<R>::intersection(Segment_2<R> &result) const
{
    if (!_known)
        intersection_type();
    if (_result != SEGMENT)
        return false;
    result = Segment_2<R>(_intersection_point, _other_point);
    return true;
}

CGAL_END_NAMESPACE



#include <CGAL/Object.h>

CGAL_BEGIN_NAMESPACE

template <class R>
Object
intersection(const Segment_2<R> &seg1, const Segment_2<R>&seg2)
{
    typedef Segment_2_Segment_2_pair<R> is_t;
    is_t ispair(&seg1, &seg2);
    switch (ispair.intersection_type()) {
    case is_t::NO:
    default:
        return Object();
    case is_t::POINT: {
        Point_2<R> pt;
        ispair.intersection(pt);
        return Object(new Wrapper< Point_2<R> >(pt));
    }
    case is_t::SEGMENT: {
        Segment_2<R> iseg;
        ispair.intersection(iseg);
        return Object(new Wrapper< Segment_2<R> >(iseg));
    }
    }
}

CGAL_END_NAMESPACE

#endif
