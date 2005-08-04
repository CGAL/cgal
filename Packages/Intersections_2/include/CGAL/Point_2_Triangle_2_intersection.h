
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


#ifndef CGAL_POINT_2_TRIANGLE_2_INTERSECTION_H
#define CGAL_POINT_2_TRIANGLE_2_INTERSECTION_H

#include <CGAL/Point_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/kernel_assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/Straight_2.h>
#include <CGAL/Object.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

template <class K>
class Point_2_Triangle_2_pair {
public:
    enum Intersection_results {NO, POINT};
    Point_2_Triangle_2_pair() ;
    Point_2_Triangle_2_pair(typename K::Point_2 const *pt,
                            typename K::Triangle_2 const *trian);
    ~Point_2_Triangle_2_pair() {}

    Intersection_results intersection_type() const;

    bool                intersection(typename K::Point_2 &result) const;
protected:
    typename K::Point_2 const *    _pt;
    typename K::Triangle_2 const * _trian;
    mutable bool                       _known;
    mutable Intersection_results       _result;
    mutable typename K::Point_2           _intersection_point;
    mutable typename K::Point_2            _other_point;
};

template <class K>
inline bool do_intersect(const typename CGAL_WRAP(K)::Point_2 &p1,
			 const typename CGAL_WRAP(K)::Triangle_2 &p2,
			 const K&)
{
  typedef Point_2_Triangle_2_pair<K> pair_t;
  pair_t pair(&p1, &p2);
  return pair.intersection_type() != pair_t::NO;
}
template <class K>
inline bool do_intersect(const typename CGAL_WRAP(K)::Triangle_2 &p2,
			 const typename CGAL_WRAP(K)::Point_2 &p1,
			 const K& k)
{
  return CGALi::do_intersect(p1, p2, k);
}


template <class K>
Point_2_Triangle_2_pair<K>::
Point_2_Triangle_2_pair()
{
    _known = false;
    _pt = 0;
    _trian = 0;
}

template <class K>
Point_2_Triangle_2_pair<K>::
Point_2_Triangle_2_pair(typename K::Point_2 const *pt,
			typename K::Triangle_2 const *trian)
{
    _known = false;
    _pt = pt;
    _trian = trian;
}

template <class K>
typename Point_2_Triangle_2_pair<K>::Intersection_results
Point_2_Triangle_2_pair<K>::intersection_type() const
{
    typedef typename K::Line_2 line_t;
    if (_known)
        return _result;
// The non const this pointer is used to cast away const.
    _known = true;
    if (_trian->has_on_unbounded_side(*_pt)) {
        _result = NO;
    } else {
        _result = POINT;
    }
    return _result;
/*
    line_t l(_trian->vertex(0), _trian->vertex(1));
    if (l.has_on_positive_side(_trian->vertex(2))) {
        for (int i=0; i<3; i++) {
            if (line_t(_trian->vertex(i), _trian->vertex(i+1)).
                                has_on_negative_side(*_pt)) {
                _result = NO;
                return _result;
            }
        }
    } else {
        for (int i=0; i<3; i++)
            if(line_t(_trian->vertex(i), _trian->vertex(i-1)).
                                has_on_negative_side(*_pt)){
                _result = NO;
                return _result;
            }
    }
*/
}



template <class K>
bool
Point_2_Triangle_2_pair<K>::
intersection(typename K::Point_2 &result) const
{
    if (!_known)
        intersection_type();
    if (_result != POINT)
        return false;
    result = *_pt;
    return true;
}



template <class K>
Object
intersection(const typename CGAL_WRAP(K)::Point_2 &pt, 
	     const typename CGAL_WRAP(K)::Triangle_2 &tr,
	     const K&)
{
    typedef Point_2_Triangle_2_pair<K> is_t;
    is_t ispair(&pt, &tr);
    switch (ispair.intersection_type()) {
    case is_t::NO:
    default:
        return Object();
    case is_t::POINT: {
        return make_object(pt);
    }
    }
}

template <class K>
inline
Object
intersection(const typename CGAL_WRAP(K)::Triangle_2 &tr,
	     const typename CGAL_WRAP(K)::Point_2 &pt, 
	     const K&k)
{
  return CGALi::intersection(pt, tr, k);
}


template <class K>
class Triangle_2_Point_2_pair
: public Point_2_Triangle_2_pair<K> {
public:
    Triangle_2_Point_2_pair(
            typename K::Triangle_2 const *trian,
            typename K::Point_2 const *pt) :
                        Point_2_Triangle_2_pair<K>(pt, trian) {}
};

} // namespace CGALi


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
CGAL_END_NAMESPACE

#endif
