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


#ifndef CGAL_RAY_2_LINE_2_INTERSECTION_H
#define CGAL_RAY_2_LINE_2_INTERSECTION_H

#include <CGAL/Line_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/kernel_assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/Object.h>
#include <CGAL/Line_2_Line_2_intersection.h>

namespace CGAL {

namespace internal {

template <class K>
class Ray_2_Line_2_pair {
public:
    enum Intersection_results {NO_INTERSECTION, POINT, RAY};
    Ray_2_Line_2_pair(typename K::Ray_2 const *ray,
		      typename K::Line_2 const *line)
      : _ray(ray), _line(line), _known(false) {}

    Intersection_results intersection_type() const;

    typename K::Point_2   intersection_point() const;
    typename K::Ray_2     intersection_ray() const;
protected:
    typename K::Ray_2 const *   _ray;
    typename K::Line_2 const *  _line;
    mutable bool                    _known;
    mutable Intersection_results    _result;
    mutable typename K::Point_2         _intersection_point;
};

template <class K>
inline bool do_intersect(
    const typename K::Ray_2 &p1,
    const typename K::Line_2 &p2)
{
    typedef Ray_2_Line_2_pair<K> pair_t;
    pair_t pair(&p1, &p2);
    return pair.intersection_type() != pair_t::NO_INTERSECTION;
}



template <class K>
Object
intersection(const typename K::Ray_2 &ray,
	     const typename K::Line_2 &line,
	     const K&)
{
    typedef Ray_2_Line_2_pair<K> is_t;
    is_t ispair(&ray, &line);
    switch (ispair.intersection_type()) {
    case is_t::NO_INTERSECTION:
    default:
        return Object();
    case is_t::POINT:
        return make_object(ispair.intersection_point());
    case is_t::RAY:
        return make_object(ray);
    }
}


template <class K>
inline
Object
intersection(const typename K::Line_2 &line,
	     const typename K::Ray_2 &ray,
	     const K& k)
{
  return internal::intersection(ray, line, k);
}


template <class K>
inline bool do_intersect(
    const typename K::Line_2 &p1,
    const typename K::Ray_2 &p2,
    const K&)
{
    typedef Ray_2_Line_2_pair<K> pair_t;
    pair_t pair(&p2, &p1);
    return pair.intersection_type() != pair_t::NO_INTERSECTION;
}



template <class K>
typename Ray_2_Line_2_pair<K>::Intersection_results
Ray_2_Line_2_pair<K>::intersection_type() const
{
    if (_known)
        return _result;
    // The non const this pointer is used to cast away const.
    _known = true;
    const typename K::Line_2 &l1 = _ray->supporting_line();
    Line_2_Line_2_pair<K> linepair(&l1, _line);
    switch ( linepair.intersection_type()) {
    case Line_2_Line_2_pair<K>::NO_INTERSECTION:
        _result = NO_INTERSECTION;
        break;
    case Line_2_Line_2_pair<K>::POINT:
        _intersection_point = linepair.intersection_point();
        _result = (_ray->collinear_has_on(_intersection_point) ) ?
                POINT : NO_INTERSECTION;
        break;
    case Line_2_Line_2_pair<K>::LINE:
        _result = RAY;
        break;
    }
    return _result;
}


template <class K>
typename K::Point_2
Ray_2_Line_2_pair<K>::intersection_point() const
{
    if (!_known)
        intersection_type();
    CGAL_kernel_assertion(_result == POINT);
    return _intersection_point;
}

template <class K>
typename K::Ray_2
Ray_2_Line_2_pair<K>::intersection_ray() const
{
    if (!_known)
        intersection_type();
    CGAL_kernel_assertion(_result == RAY);
    return *_ray;
}

} // namespace internal

template <class K>
inline bool do_intersect(const Line_2<K> &p1, const Ray_2<K> &p2)
{
  typedef typename K::Do_intersect_2 Do_intersect;
  return Do_intersect()(p1, p2);
}

template <class K>
inline bool do_intersect(const Ray_2<K> &p2, const Line_2<K> &p1)
{
  typedef typename K::Do_intersect_2 Do_intersect;
  return Do_intersect()(p1, p2);
}

template <class K>
inline Object
intersection(const Line_2<K> &line, const Ray_2<K> &ray)
{
  typedef typename K::Intersect_2 Intersect;
    return Intersect()(ray, line);
}

template <class K>
inline Object
intersection(const Ray_2<K> &ray, const Line_2<K> &line)
{
  typedef typename K::Intersect_2 Intersect;
    return Intersect()(ray, line);
}

} //namespace CGAL

#endif
