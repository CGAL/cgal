
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

    Intersection_results intersection_type() const;

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
                            (_ray1->source().x() < _ray2->source().x())
                            ? _ray2->source() : _ray1->source();
                    _result = RAY;
                    return _result;
                } else {
                    if (_ray1->source().x() > _ray2->source().x()) {
                        _result = NO;
                        return _result;
                    }
                    if (_ray1->source().x() == _ray2->source().x()) {
                        _intersection_point = _ray1->source();
                        _result = POINT;
                        return _result;
                    }
                    _result = SEGMENT;
                    return _result;
                }
            } else {
                if (dir2.x() < FT(0)) {
                    _intersection_point =
                            (_ray1->source().x() > _ray2->source().x())
                            ? _ray2->source() : _ray1->source();
                    _result = RAY;
                    return _result;
                } else {
                    if (_ray1->source().x() < _ray2->source().x()) {
                        _result = NO;
                        return _result;
                    }
                    if (_ray1->source().x() == _ray2->source().x()) {
                        _intersection_point = _ray1->source();
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
                            (_ray1->source().y() < _ray2->source().y())
                            ? _ray2->source() : _ray1->source();
                    _result = RAY;
                    return _result;
                } else {
                    if (_ray1->source().y() > _ray2->source().y()) {
                        _result = NO;
                        return _result;
                    }
                    if (_ray1->source().y() == _ray2->source().y()) {
                        _intersection_point = _ray1->source();
                        _result = POINT;
                        return _result;
                    }
                    _result = SEGMENT;
                    return _result;
                }
            } else {
                if (dir2.y() < FT(0)) {
                    _intersection_point =
                            (_ray1->source().y() > _ray2->source().y())
                            ? _ray2->source() : _ray1->source();
                    _result = RAY;
                    return _result;
                } else {
                    if (_ray1->source().y() < _ray2->source().y()) {
                        _result = NO;
                        return _result;
                    }
                    if (_ray1->source().y() == _ray2->source().y()) {
                        _intersection_point = _ray1->source();
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
    result = Segment_2(_ray1->source(), _ray2->source());
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
  typedef typename K::Do_intersect_2 Do_intersect;
  return Do_intersect()(ray1, ray2);
}

template <class K>
Object
intersection(const  Ray_2<K> &ray1, 
	     const Ray_2<K> &ray2)
{
  typedef typename K::Intersect_2 Intersect;
  return Intersect()(ray1, ray2);
}

CGAL_END_NAMESPACE

#endif
