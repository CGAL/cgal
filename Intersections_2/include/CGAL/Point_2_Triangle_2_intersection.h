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


#ifndef CGAL_POINT_2_TRIANGLE_2_INTERSECTION_H
#define CGAL_POINT_2_TRIANGLE_2_INTERSECTION_H

#include <CGAL/Point_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/kernel_assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/Straight_2.h>
#include <CGAL/Object.h>

namespace CGAL {

namespace internal {

template <class K>
class Point_2_Triangle_2_pair {
public:
    enum Intersection_results {NO_INTERSECTION, POINT};
    Point_2_Triangle_2_pair(typename K::Point_2 const *pt,
                            typename K::Triangle_2 const *trian)
	    : _pt(pt), _trian(trian), _known(false) {}

    Intersection_results intersection_type() const;

    typename K::Point_2  intersection_point() const;
protected:
    typename K::Point_2 const *    _pt;
    typename K::Triangle_2 const * _trian;
    mutable bool                   _known;
    mutable Intersection_results   _result;
    mutable typename K::Point_2    _intersection_point;
    mutable typename K::Point_2    _other_point;
};

template <class K>
inline bool do_intersect(const typename K::Point_2 &p1,
			 const typename K::Triangle_2 &p2,
			 const K&)
{
  typedef Point_2_Triangle_2_pair<K> pair_t;
  pair_t pair(&p1, &p2);
  return pair.intersection_type() != pair_t::NO_INTERSECTION;
}

template <class K>
inline bool do_intersect(const typename K::Triangle_2 &p2,
			 const typename K::Point_2 &p1,
			 const K& k)
{
  return internal::do_intersect(p1, p2, k);
}


template <class K>
typename Point_2_Triangle_2_pair<K>::Intersection_results
Point_2_Triangle_2_pair<K>::intersection_type() const
{
    if (_known)
        return _result;
// The non const this pointer is used to cast away const.
    _known = true;
    if (_trian->has_on_unbounded_side(*_pt)) {
        _result = NO_INTERSECTION;
    } else {
        _result = POINT;
    }
    return _result;
/*
    typedef typename K::Line_2 line_t;
    line_t l(_trian->vertex(0), _trian->vertex(1));
    if (l.has_on_positive_side(_trian->vertex(2))) {
        for (int i=0; i<3; i++) {
            if (line_t(_trian->vertex(i), _trian->vertex(i+1)).
                                has_on_negative_side(*_pt)) {
                _result = NO_INTERSECTION;
                return _result;
            }
        }
    } else {
        for (int i=0; i<3; i++)
            if(line_t(_trian->vertex(i), _trian->vertex(i-1)).
                                has_on_negative_side(*_pt)){
                _result = NO_INTERSECTION;
                return _result;
            }
    }
*/
}



template <class K>
typename K::Point_2
Point_2_Triangle_2_pair<K>::
intersection_point() const
{
    if (!_known)
        intersection_type();
    CGAL_kernel_assertion(_result == POINT);
    return *_pt;
}



template <class K>
Object
intersection(const typename K::Point_2 &pt, 
	     const typename K::Triangle_2 &tr,
	     const K&)
{
    typedef Point_2_Triangle_2_pair<K> is_t;
    is_t ispair(&pt, &tr);
    switch (ispair.intersection_type()) {
    case is_t::NO_INTERSECTION:
    default:
        return Object();
    case is_t::POINT:
        return make_object(pt);
    }
}

template <class K>
inline
Object
intersection(const typename K::Triangle_2 &tr,
	     const typename K::Point_2 &pt, 
	     const K&k)
{
  return internal::intersection(pt, tr, k);
}

} // namespace internal


template <class K>
inline 
bool
do_intersect(const Triangle_2<K> &tr, const Point_2<K> &pt)
{
  typedef typename K::Do_intersect_2 Do_intersect;
  return Do_intersect()(pt, tr);
}

template <class K>
inline 
bool 
do_intersect(const Point_2<K> &pt, const Triangle_2<K> &tr)
{
  typedef typename K::Do_intersect_2 Do_intersect;
  return Do_intersect()(pt, tr);
}

template <class K>
inline Object
intersection(const Triangle_2<K> &tr, const Point_2<K> &pt)
{
  typedef typename K::Intersect_2 Intersect;
  return Intersect()(pt, tr);
}

template <class K>
inline Object
intersection(const Point_2<K> &pt, const Triangle_2<K> &tr)
{
  typedef typename K::Intersect_2 Intersect;
  return Intersect()(pt, tr);
}

} //namespace CGAL

#endif
