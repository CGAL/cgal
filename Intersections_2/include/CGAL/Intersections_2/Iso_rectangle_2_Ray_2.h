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


#ifndef CGAL_INTERSECTIONS_2_ISO_RECTANGLE_2_RAY_2_H
#define CGAL_INTERSECTIONS_2_ISO_RECTANGLE_2_RAY_2_H

#include <CGAL/disable_warnings.h>

#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/kernel_assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/Intersection_traits_2.h>

namespace CGAL {

namespace Intersections {

namespace internal {

template <class K>
class Ray_2_Iso_rectangle_2_pair {
public:
    enum Intersection_results {NO_INTERSECTION, POINT, SEGMENT};
    Ray_2_Iso_rectangle_2_pair(typename K::Ray_2 const *ray,
                          typename K::Iso_rectangle_2 const *iso)
      : _known(false), _ref_point(ray->source()), _dir(ray->direction().to_vector()),
        _isomin((iso->min)()), _isomax((iso->max)()), _min((typename K::FT)(0)) {}

    Intersection_results intersection_type() const;

    typename K::Point_2   intersection_point() const;
    typename K::Segment_2 intersection_segment() const;
protected:
    mutable bool                       _known;
    mutable Intersection_results       _result;
    mutable typename K::Point_2            _ref_point;
    mutable typename K::Vector_2           _dir;
    mutable typename K::Point_2            _isomin;
    mutable typename K::Point_2            _isomax;
    mutable typename K::FT                     _min, _max;
};

template <class K>
inline bool do_intersect(const typename K::Ray_2 &p1,
                         const typename K::Iso_rectangle_2 &p2,
                         const K&)
{
    typedef Ray_2_Iso_rectangle_2_pair<K> pair_t;
    pair_t pair(&p1, &p2);
    return pair.intersection_type() != pair_t::NO_INTERSECTION;
}

template <class K>
inline bool do_intersect(const typename K::Iso_rectangle_2 &p2,
                         const typename K::Ray_2 &p1,
                         const K& k)
{
  return do_intersect(p1, p2, k);
}

template <class K>
typename CGAL::Intersection_traits
<K, typename K::Ray_2, typename K::Iso_rectangle_2>::result_type
intersection(const typename K::Ray_2 &ray,
             const typename K::Iso_rectangle_2 &iso,
             const K& )
{
    typedef Ray_2_Iso_rectangle_2_pair<K> is_t;
    is_t ispair(&ray, &iso);
    switch (ispair.intersection_type()) {
    case is_t::NO_INTERSECTION:
    default:
        return intersection_return<typename K::Intersect_2, typename K::Ray_2, typename K::Iso_rectangle_2>();
    case is_t::POINT:
        return intersection_return<typename K::Intersect_2, typename K::Ray_2, typename K::Iso_rectangle_2>(ispair.intersection_point());
    case is_t::SEGMENT:
        return intersection_return<typename K::Intersect_2, typename K::Ray_2, typename K::Iso_rectangle_2>(ispair.intersection_segment());
    }
}

template <class K>
typename CGAL::Intersection_traits
<K, typename K::Ray_2, typename K::Iso_rectangle_2>::result_type
intersection(const typename K::Iso_rectangle_2 &iso,
             const typename K::Ray_2 &ray,
             const K& k)
{
  return intersection(ray, iso, k);
}

template <class K>
typename Ray_2_Iso_rectangle_2_pair<K>::Intersection_results
Ray_2_Iso_rectangle_2_pair<K>::intersection_type() const
{
    typedef typename K::RT RT;
    typedef typename K::FT FT;
    if (_known)
        return _result;
    _known = true;
    bool to_infinity = true;

    typename K::Construct_cartesian_const_iterator_2 construct_cccit;
    typename K::Cartesian_const_iterator_2 ref_point_it = construct_cccit(_ref_point);
    typename K::Cartesian_const_iterator_2 end = construct_cccit(_ref_point, 0);
    typename K::Cartesian_const_iterator_2 isomin_it = construct_cccit(_isomin);
    typename K::Cartesian_const_iterator_2 isomax_it = construct_cccit(_isomax);

    for (unsigned int i=0; ref_point_it != end; ++i, ++ref_point_it, ++isomin_it, ++isomax_it) {
        if (_dir.homogeneous(i) == RT(0)) {
            if ((*ref_point_it) < (*isomin_it)) {
                _result = NO_INTERSECTION;
                return _result;
            }
            if ((*ref_point_it) > (*isomax_it)) {
                _result = NO_INTERSECTION;
                return _result;
            }
        } else {
            FT newmin, newmax;
            if (_dir.homogeneous(i) > RT(0)) {
                newmin = (*isomin_it - *ref_point_it) /
                    _dir.cartesian(i);
                newmax = (*isomax_it - *ref_point_it) /
                    _dir.cartesian(i);
            } else {
                newmin = (*isomax_it - *ref_point_it) /
                    _dir.cartesian(i);
                newmax = (*isomin_it - *ref_point_it) /
                    _dir.cartesian(i);
            }
            if (newmin > _min)
                _min = newmin;
            if (to_infinity) {
                _max = newmax;
            } else {
                if (newmax < _max)
                    _max = newmax;
            }
            if (_max < _min) {
                _result = NO_INTERSECTION;
                return _result;
            }
            to_infinity = false;
        }
    }
    CGAL_kernel_assertion(!to_infinity);
    if (_max == _min) {
        _result = POINT;
        return _result;
    }
    _result = SEGMENT;
    return _result;
}


template <class K>
typename K::Segment_2
Ray_2_Iso_rectangle_2_pair<K>::intersection_segment() const
{
    typedef typename K::Segment_2 Segment_2;
    typename K::Construct_translated_point_2 translated_point;
    typename K::Construct_scaled_vector_2 construct_scaled_vector;
    if (!_known)
        intersection_type();
    CGAL_kernel_assertion(_result == SEGMENT);
    typename K::Point_2 p1(translated_point(_ref_point, construct_scaled_vector(_dir,_min)));
    typename K::Point_2 p2(translated_point(_ref_point, construct_scaled_vector(_dir,_max)));
    return Segment_2(p1, p2);
}

template <class K>
typename K::Point_2
Ray_2_Iso_rectangle_2_pair<K>::intersection_point() const
{
    typedef typename K::Point_2 Point_2;
    typename K::Construct_translated_point_2 translated_point;
    typename K::Construct_scaled_vector_2 construct_scaled_vector;
    if (!_known)
        intersection_type();
    CGAL_kernel_assertion(_result == POINT);
    return Point_2(translated_point(_ref_point, construct_scaled_vector(_dir, _min)));
}

} // namespace internal
} // namespace Intersections

CGAL_INTERSECTION_FUNCTION(Ray_2, Iso_rectangle_2, 2)
CGAL_DO_INTERSECT_FUNCTION(Ray_2, Iso_rectangle_2, 2)


} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_RAY_2_iSO_RECTANGLE_2_INTERSECTION_H
