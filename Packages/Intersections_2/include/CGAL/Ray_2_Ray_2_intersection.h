
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
// file          : include/CGAL/Ray_2_Ray_2_intersection.h
// source        : intersection_2_1.fw
// author(s)     : Geert-Jan Giezeman
//
// coordinator   : Saarbruecken
//
// ============================================================================


#ifndef CGAL_RAY_2_RAY_2_INTERSECTION_H
#define CGAL_RAY_2_RAY_2_INTERSECTION_H

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
class Ray_2_Ray_2_pair {
public:
    enum Intersection_results {NO, POINT, SEGMENT, RAY};
    Ray_2_Ray_2_pair() ;
    Ray_2_Ray_2_pair(typename K::Ray_2 const *ray1,
		     typename K::Ray_2 const *ray2);
    ~Ray_2_Ray_2_pair() {}

#ifndef CGAL_CFG_RETURN_TYPE_BUG_2
    Intersection_results intersection_type() const;
#else
    Intersection_results intersection_type() const
{
    if (_known)
        return _result;
    // The non const this pointer is used to cast away const.
    _known = true;
//    if (!do_overlap(_ray1->bbox(), _ray2->bbox()))
//        return NO;
    const typename K::Line_2 &l1 = _ray1->supporting_line();
    const typename K::Line_2 &l2 = _ray2->supporting_line();
    Line_2_Line_2_pair<K> linepair(&l1, &l2);
    switch ( linepair.intersection_type()) {
    case Line_2_Line_2_pair<K>::NO:
        _result = NO;
        return _result;
    case Line_2_Line_2_pair<K>::POINT:
        linepair.intersection(_intersection_point);
        _result = (_ray1->collinear_has_on(_intersection_point)
                && _ray2->collinear_has_on(_intersection_point) )
            ? POINT :  NO;
        return _result;
    case Line_2_Line_2_pair<K>::LINE:
        {
        typedef typename K::RT RT;
        const typename K::Vector_2 &dir1 = _ray1->direction().to_vector();
        const typename K::Vector_2 &dir2 = _ray2->direction().to_vector();
        if (CGAL_NTS abs(dir1.x()) > CGAL_NTS abs(dir1.y())) {
            typedef typename K::FT FT;
            if (dir1.x() > FT(0)) {
                if (dir2.x() > FT(0)) {
                    _intersection_point =
                            (_ray1->start().x() < _ray2->start().x())
                            ? _ray2->start() : _ray1->start();
                    _result = RAY;
                    return _result;
                } else {
                    if (_ray1->start().x() > _ray2->start().x()) {
                        _result = NO;
                        return _result;
                    }
                    if (_ray1->start().x() == _ray2->start().x()) {
                        _intersection_point = _ray1->start();
                        _result = POINT;
                        return _result;
                    }
                    _result = SEGMENT;
                    return _result;
                }
            } else {
                if (dir2.x() < FT(0)) {
                    _intersection_point =
                            (_ray1->start().x() > _ray2->start().x())
                            ? _ray2->start() : _ray1->start();
                    _result = RAY;
                    return _result;
                } else {
                    if (_ray1->start().x() < _ray2->start().x()) {
                        _result = NO;
                        return _result;
                    }
                    if (_ray1->start().x() == _ray2->start().x()) {
                        _intersection_point = _ray1->start();
                        _result = POINT;
                        return _result;
                    }
                    _result = SEGMENT;
                    return _result;
                }
            }
            
        } else {
            typedef typename K::FT FT;
            if (dir1.y() > FT(0)) {
                if (dir2.y() > FT(0)) {
                    _intersection_point =
                            (_ray1->start().y() < _ray2->start().y())
                            ? _ray2->start() : _ray1->start();
                    _result = RAY;
                    return _result;
                } else {
                    if (_ray1->start().y() > _ray2->start().y()) {
                        _result = NO;
                        return _result;
                    }
                    if (_ray1->start().y() == _ray2->start().y()) {
                        _intersection_point = _ray1->start();
                        _result = POINT;
                        return _result;
                    }
                    _result = SEGMENT;
                    return _result;
                }
            } else {
                if (dir2.y() < FT(0)) {
                    _intersection_point =
                            (_ray1->start().y() > _ray2->start().y())
                            ? _ray2->start() : _ray1->start();
                    _result = RAY;
                    return _result;
                } else {
                    if (_ray1->start().y() < _ray2->start().y()) {
                        _result = NO;
                        return _result;
                    }
                    if (_ray1->start().y() == _ray2->start().y()) {
                        _intersection_point = _ray1->start();
                        _result = POINT;
                        return _result;
                    }
                    _result = SEGMENT;
                    return _result;
                }
            }
            
        }
        } 
    default:
        CGAL_kernel_assertion(false); // should not be reached:
        return _result;
    }
}

#endif // CGAL_CFG_RETURN_TYPE_BUG_2

    bool                intersection(typename K::Point_2 &result) const;
    bool                intersection(typename K::Segment_2 &result) const;
    bool                intersection(typename K::Ray_2 &result) const;
protected:
    typename K::Ray_2 const*    _ray1;
    typename K::Ray_2 const *   _ray2;
    mutable bool                    _known;
    mutable Intersection_results    _result;
    mutable typename K::Point_2         _intersection_point, _other_point;
};

template <class K>
inline bool do_intersect(
    const typename CGAL_WRAP(K)::Ray_2 &p1,
    const typename CGAL_WRAP(K)::Ray_2 &p2,
    const K& k)
{
    typedef Ray_2_Ray_2_pair<K> pair_t;
    pair_t pair(&p1, &p2);
    return pair.intersection_type() != pair_t::NO;
}



template <class K>
Ray_2_Ray_2_pair<K>::Ray_2_Ray_2_pair()
{
    _ray1 = 0;
    _ray2 = 0;
    _known = false;
}

template <class K>
Ray_2_Ray_2_pair<K>::Ray_2_Ray_2_pair(
    typename K::Ray_2 const *ray1, typename K::Ray_2 const *ray2)
{
    _ray1 = ray1;
    _ray2 = ray2;
    _known = false;
}

#ifndef CGAL_CFG_RETURN_TYPE_BUG_2
template <class K>
typename Ray_2_Ray_2_pair<K>::Intersection_results
Ray_2_Ray_2_pair<K>::intersection_type() const
{
    if (_known)
        return _result;
    // The non const this pointer is used to cast away const.
    _known = true;
//    if (!do_overlap(_ray1->bbox(), _ray2->bbox()))
//        return NO;
    const typename K::Line_2 &l1 = _ray1->supporting_line();
    const typename K::Line_2 &l2 = _ray2->supporting_line();
    Line_2_Line_2_pair<K> linepair(&l1, &l2);
    switch ( linepair.intersection_type()) {
    case Line_2_Line_2_pair<K>::NO:
        _result = NO;
        return _result;
    case Line_2_Line_2_pair<K>::POINT:
        linepair.intersection(_intersection_point);
        _result = (_ray1->collinear_has_on(_intersection_point)
                && _ray2->collinear_has_on(_intersection_point) )
            ? POINT :  NO;
        return _result;
    case Line_2_Line_2_pair<K>::LINE:
        {
        typedef typename K::RT RT;
        const typename K::Vector_2 &dir1 = _ray1->direction().to_vector();
        const typename K::Vector_2 &dir2 = _ray2->direction().to_vector();
        if (CGAL_NTS abs(dir1.x()) > CGAL_NTS abs(dir1.y())) {
            typedef typename K::FT FT;
            if (dir1.x() > FT(0)) {
                if (dir2.x() > FT(0)) {
                    _intersection_point =
                            (_ray1->start().x() < _ray2->start().x())
                            ? _ray2->start() : _ray1->start();
                    _result = RAY;
                    return _result;
                } else {
                    if (_ray1->start().x() > _ray2->start().x()) {
                        _result = NO;
                        return _result;
                    }
                    if (_ray1->start().x() == _ray2->start().x()) {
                        _intersection_point = _ray1->start();
                        _result = POINT;
                        return _result;
                    }
                    _result = SEGMENT;
                    return _result;
                }
            } else {
                if (dir2.x() < FT(0)) {
                    _intersection_point =
                            (_ray1->start().x() > _ray2->start().x())
                            ? _ray2->start() : _ray1->start();
                    _result = RAY;
                    return _result;
                } else {
                    if (_ray1->start().x() < _ray2->start().x()) {
                        _result = NO;
                        return _result;
                    }
                    if (_ray1->start().x() == _ray2->start().x()) {
                        _intersection_point = _ray1->start();
                        _result = POINT;
                        return _result;
                    }
                    _result = SEGMENT;
                    return _result;
                }
            }
            
        } else {
            typedef typename K::FT FT;
            if (dir1.y() > FT(0)) {
                if (dir2.y() > FT(0)) {
                    _intersection_point =
                            (_ray1->start().y() < _ray2->start().y())
                            ? _ray2->start() : _ray1->start();
                    _result = RAY;
                    return _result;
                } else {
                    if (_ray1->start().y() > _ray2->start().y()) {
                        _result = NO;
                        return _result;
                    }
                    if (_ray1->start().y() == _ray2->start().y()) {
                        _intersection_point = _ray1->start();
                        _result = POINT;
                        return _result;
                    }
                    _result = SEGMENT;
                    return _result;
                }
            } else {
                if (dir2.y() < FT(0)) {
                    _intersection_point =
                            (_ray1->start().y() > _ray2->start().y())
                            ? _ray2->start() : _ray1->start();
                    _result = RAY;
                    return _result;
                } else {
                    if (_ray1->start().y() < _ray2->start().y()) {
                        _result = NO;
                        return _result;
                    }
                    if (_ray1->start().y() == _ray2->start().y()) {
                        _intersection_point = _ray1->start();
                        _result = POINT;
                        return _result;
                    }
                    _result = SEGMENT;
                    return _result;
                }
            }
            
        }
        } 
    default:
        CGAL_kernel_assertion(false); // should not be reached:
        return _result;
    }
}

#endif // CGAL_CFG_RETURN_TYPE_BUG_2

template <class K>
bool
Ray_2_Ray_2_pair<K>::intersection(typename K::Point_2 &result) const
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
Ray_2_Ray_2_pair<K>::intersection(typename K::Segment_2 &result) const
{
  typedef typename K::Segment_2 Segment_2;
    if (!_known)
        intersection_type();
    if (_result != SEGMENT)
        return false;
    result = Segment_2(_ray1->start(), _ray2->start());
    return true;
}

template <class K>
bool
Ray_2_Ray_2_pair<K>::intersection(typename K::Ray_2 &result) const
{
  typedef typename K::Ray_2 Ray_2;
    if (!_known)
        intersection_type();
    if (_result != RAY)
        return false;
    result = Ray_2(_intersection_point, _ray1->direction());
    return true;
}



template <class K>
Object
intersection(const typename CGAL_WRAP(K)::Ray_2 &ray1, 
	     const typename CGAL_WRAP(K)::Ray_2 &ray2,
	     const K&)
{
    typedef Ray_2_Ray_2_pair<K> is_t;
    is_t ispair(&ray1, &ray2);
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
    case is_t::RAY: {
        typename K::Ray_2 iray;
        ispair.intersection(iray);
        return make_object(iray);
    }
    }
}

} // namespace CGALi


template <class K>
inline
bool
do_intersect(const Ray_2<K> &ray1, 
	     const Ray_2<K> &ray2)
{
  return CGALi::do_intersect(ray1, ray2, K());
}

template <class K>
Object
intersection(const  Ray_2<K> &ray1, 
	     const Ray_2<K> &ray2)
{
  return CGALi::intersection(ray1, ray2, K());
}

CGAL_END_NAMESPACE

#endif
