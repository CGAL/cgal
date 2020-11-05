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


#ifndef CGAL_INTERSECTIONS_2_POINT_2_TRIANGLE_2_H
#define CGAL_INTERSECTIONS_2_POINT_2_TRIANGLE_2_H

#include <CGAL/Point_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/kernel_assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/Intersections_2/internal/Straight_2.h>

namespace CGAL {

namespace Intersections {

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
typename Intersection_traits<K, typename K::Point_2, typename K::Triangle_2>
::result_type
intersection(const typename K::Point_2 &pt,
             const typename K::Triangle_2 &tr,
             const K&)
{
    typedef Point_2_Triangle_2_pair<K> is_t;
    is_t ispair(&pt, &tr);
    switch (ispair.intersection_type()) {
    case is_t::NO_INTERSECTION:
    default:
        return intersection_return<typename K::Intersect_2, typename K::Point_2, typename K::Triangle_2>();
    case is_t::POINT:
        return intersection_return<typename K::Intersect_2, typename K::Point_2, typename K::Triangle_2>(pt);
    }
}

template <class K>
inline
typename Intersection_traits<K, typename K::Point_2, typename K::Triangle_2>
::result_type
intersection(const typename K::Triangle_2 &tr,
             const typename K::Point_2 &pt,
             const K&k)
{
  return internal::intersection(pt, tr, k);
}

} // namespace internal
} // namespace Intersections

CGAL_INTERSECTION_FUNCTION(Point_2, Triangle_2, 2)
CGAL_DO_INTERSECT_FUNCTION(Point_2, Triangle_2, 2)

} //namespace CGAL

#endif
