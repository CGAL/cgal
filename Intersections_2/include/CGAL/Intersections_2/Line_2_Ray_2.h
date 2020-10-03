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


#ifndef CGAL_INTERSECTIONS_2_LINE_2_RAY_2_H
#define CGAL_INTERSECTIONS_2_LINE_2_RAY_2_H

#include <CGAL/Line_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/kernel_assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/Intersections_2/Line_2_Line_2.h>
#include <CGAL/Intersection_traits_2.h>

namespace CGAL {

namespace Intersections {

namespace internal {

template <class K>
class Ray_2_Line_2_pair {
public:
    enum Intersection_results {NOT_COMPUTED_YET, NO_INTERSECTION, POINT, RAY};
    typedef typename K::FT FT;
    Ray_2_Line_2_pair(typename K::Ray_2 const *ray,
                      typename K::Line_2 const *line)
      : _ray(ray), _line(line), _result(NOT_COMPUTED_YET),
        _intersection_point(K().construct_point_2_object()(ORIGIN))
    {}

    Intersection_results intersection_type() const;

    typename K::Point_2   intersection_point() const;
    typename K::Ray_2     intersection_ray() const;
protected:
    typename K::Ray_2 const *   _ray;
    typename K::Line_2 const *  _line;
    mutable Intersection_results    _result;
    mutable typename K::Point_2         _intersection_point;
};

template <class K>
inline bool do_intersect(
    const typename K::Ray_2 &p1,
    const typename K::Line_2 &p2,
    const K&)
{
    typedef Ray_2_Line_2_pair<K> pair_t;
    pair_t pair(&p1, &p2);
    return pair.intersection_type() != pair_t::NO_INTERSECTION;
}



template <class K>
typename Intersection_traits
<K, typename K::Ray_2, typename K::Line_2>::result_type
intersection(const typename K::Ray_2 &ray,
             const typename K::Line_2 &line,
             const K&)
{
    typedef Ray_2_Line_2_pair<K> is_t;
    is_t ispair(&ray, &line);
    switch (ispair.intersection_type()) {
    case is_t::NO_INTERSECTION:
    default:
        return intersection_return<typename K::Intersect_2, typename K::Ray_2, typename K::Line_2>();
    case is_t::POINT:
        return intersection_return<typename K::Intersect_2, typename K::Ray_2, typename K::Line_2>(ispair.intersection_point());
    case is_t::RAY:
        return intersection_return<typename K::Intersect_2, typename K::Ray_2, typename K::Line_2>(ray);
    }
}


template <class K>
inline
typename Intersection_traits
<K, typename K::Line_2, typename K::Ray_2>::result_type
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
    if (_result != NOT_COMPUTED_YET)
        return _result;
    // The non const this pointer is used to cast away const.
    const typename K::Line_2 &l1 = _ray->supporting_line();
    Line_2_Line_2_pair<K> linepair(&l1, _line);
    switch ( linepair.intersection_type()) {
    case Line_2_Line_2_pair<K>::NO_INTERSECTION:
    default:
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
    if (_result == NOT_COMPUTED_YET)
        intersection_type();
    CGAL_kernel_assertion(_result == POINT);
    return _intersection_point;
}

template <class K>
typename K::Ray_2
Ray_2_Line_2_pair<K>::intersection_ray() const
{
    if (_result == NOT_COMPUTED_YET)
        intersection_type();
    CGAL_kernel_assertion(_result == RAY);
    return *_ray;
}

} // namespace internal
} // namespace Intersections

CGAL_INTERSECTION_FUNCTION(Ray_2, Line_2, 2)
CGAL_DO_INTERSECT_FUNCTION(Ray_2, Line_2, 2)

} //namespace CGAL

#endif
