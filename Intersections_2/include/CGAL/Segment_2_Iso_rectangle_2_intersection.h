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
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Geert-Jan Giezeman


#ifndef CGAL_SEGMENT_2_ISO_RECTANGLE_2_INTERSECTION_H
#define CGAL_SEGMENT_2_ISO_RECTANGLE_2_INTERSECTION_H

#include <CGAL/disable_warnings.h>

#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/kernel_assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/Intersection_traits_2.h>


namespace CGAL {
namespace internal {

template <class K>
class Segment_2_Iso_rectangle_2_pair {
public:
    enum Intersection_results {NO_INTERSECTION, POINT, SEGMENT};
    Segment_2_Iso_rectangle_2_pair(typename K::Segment_2 const *seg,
                          typename K::Iso_rectangle_2 const *rect) ;

    Intersection_results intersection_type() const;

    typename K::Point_2 intersection_point() const;
    typename K::Segment_2 intersection_segment() const;
protected:
    mutable bool                       _known;
    mutable Intersection_results       _result;
    mutable typename K::Point_2            _ref_point;
    mutable typename K::Vector_2           _dir;
    mutable typename K::Point_2            _isomin;
    mutable typename K::Point_2            _isomax;
    mutable typename K::FT              _min,
                               _max;
};

template <class K>
inline bool do_intersect(
    const typename K::Segment_2 &p1,
    const typename K::Iso_rectangle_2 &p2,
    const K&)
{
    typedef Segment_2_Iso_rectangle_2_pair<K> pair_t;
    pair_t pair(&p1, &p2);
    return pair.intersection_type() != pair_t::NO_INTERSECTION;
}





template <class K>
typename CGAL::Intersection_traits
<K, typename K::Segment_2, typename K::Iso_rectangle_2>::result_type
intersection(
    const typename K::Segment_2 &seg,
    const typename K::Iso_rectangle_2 &iso,
    const K&)
{
    typedef Segment_2_Iso_rectangle_2_pair<K> is_t;
    is_t ispair(&seg, &iso);
    switch (ispair.intersection_type()) {
    case is_t::NO_INTERSECTION:
    default:
        return intersection_return<typename K::Intersect_2, typename K::Segment_2, typename K::Iso_rectangle_2>();
    case is_t::POINT:
        return intersection_return<typename K::Intersect_2, typename K::Segment_2, typename K::Iso_rectangle_2>(ispair.intersection_point());
    case is_t::SEGMENT:
        return intersection_return<typename K::Intersect_2, typename K::Segment_2, typename K::Iso_rectangle_2>(ispair.intersection_segment());
    }
}

template <class K>
inline
typename CGAL::Intersection_traits
<K, typename K::Segment_2, typename K::Iso_rectangle_2>::result_type
intersection(const typename K::Iso_rectangle_2 &iso,
	     const typename K::Segment_2 &seg,
	     const K& k)
{
  return internal::intersection(seg, iso, k);
}

template <class K>
Segment_2_Iso_rectangle_2_pair<K>::
Segment_2_Iso_rectangle_2_pair(
        typename K::Segment_2 const *seg,
        typename K::Iso_rectangle_2 const *iso)
{
    _known = false;
    _isomin = (iso->min)();
    _isomax = (iso->max)();
    _ref_point = seg->source();
    _dir = seg->direction().to_vector();
    _min = (typename K::FT)(0);
    int main_dir = (CGAL_NTS abs(_dir.x()) > CGAL_NTS abs(_dir.y()) ) ? 0 : 1;

  typename K::Construct_cartesian_const_iterator_2 construct_cccit;
  typename K::Cartesian_const_iterator_2 seg_target_it = construct_cccit(seg->target()) + main_dir;
  typename K::Cartesian_const_iterator_2 ref_point_it = construct_cccit(_ref_point) + main_dir;

    _max = (*seg_target_it - *ref_point_it) /
            _dir.cartesian(main_dir);
}

template <class K>
typename Segment_2_Iso_rectangle_2_pair<K>::Intersection_results
Segment_2_Iso_rectangle_2_pair<K>::intersection_type() const
{
    typedef typename K::RT RT;
    typedef typename K::FT FT;
    if (_known)
        return _result;
    _known = true;

    typename K::Construct_cartesian_const_iterator_2 construct_cccit;
    typename K::Cartesian_const_iterator_2 ref_point_it = construct_cccit(_ref_point);
    typename K::Cartesian_const_iterator_2 end = construct_cccit(_ref_point, 0);
    typename K::Cartesian_const_iterator_2 isomin_it = construct_cccit(_isomin);
    typename K::Cartesian_const_iterator_2 isomax_it = construct_cccit(_isomax);

    for (unsigned int i=0; ref_point_it != end; ++i, ++ref_point_it, ++isomin_it, ++isomax_it) {
        if (_dir.homogeneous(i) == RT(0)) {
            if ( *(ref_point_it) < *(isomin_it) ) {
                _result = NO_INTERSECTION;
                return _result;
            }
            if ( *(ref_point_it) > *(isomax_it)) {
                _result = NO_INTERSECTION;
                return _result;
            }
        } else {
            FT newmin, newmax;
            if (_dir.homogeneous(i) > RT(0)) {
                newmin = ( *(isomin_it) - (*ref_point_it)) /
		  _dir.cartesian(i);
                newmax = ( *(isomax_it) - (*ref_point_it)) /
                    _dir.cartesian(i);
            } else {
                newmin = ( (*isomax_it) - (*ref_point_it)) /
                    _dir.cartesian(i);
                newmax = ( (*isomin_it) - (*ref_point_it)) /
                    _dir.cartesian(i);
            }
            if (newmin > _min)
                _min = newmin;
            if (newmax < _max)
                _max = newmax;
            if (_max < _min) {
                _result = NO_INTERSECTION;
                return _result;
            }
        }
    }
    if (_max == _min) {
        _result = POINT;
        return _result;
    }
    _result = SEGMENT;
    return _result;
}


template <class K>
typename K::Segment_2
Segment_2_Iso_rectangle_2_pair<K>::
intersection_segment() const
{
  typedef typename K::Segment_2 Segment_2;
  typename K::Construct_translated_point_2 translated_point;
  typename K::Construct_scaled_vector_2 construct_scaled_vector;
    if (!_known)
        intersection_type();
    CGAL_kernel_assertion(_result == SEGMENT);
    typename K::Point_2 p1(translated_point(_ref_point,  construct_scaled_vector(_dir,_min)));
    typename K::Point_2 p2(translated_point(_ref_point, construct_scaled_vector(_dir,_max)));
    return Segment_2(p1, p2);
}

template <class K>
typename K::Point_2
Segment_2_Iso_rectangle_2_pair<K>::
intersection_point() const
{
  typename K::Construct_translated_point_2 translated_point;
  typename K::Construct_scaled_vector_2 construct_scaled_vector;
    if (!_known)
        intersection_type();
    CGAL_kernel_assertion(_result == POINT);
    return translated_point(_ref_point, construct_scaled_vector(_dir,_min));
}



template <class K>
inline bool do_intersect(
    const typename K::Iso_rectangle_2 &p1,
    const typename K::Segment_2 &p2,
    const K&)
{
    typedef Segment_2_Iso_rectangle_2_pair<K> pair_t;
    pair_t pair(&p2, &p1);
    return pair.intersection_type() != pair_t::NO_INTERSECTION;
}

} // namespace internal

CGAL_INTERSECTION_FUNCTION(Segment_2, Iso_rectangle_2, 2)
CGAL_DO_INTERSECT_FUNCTION(Segment_2, Iso_rectangle_2, 2)

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_SEGMENT_2_ISO_RECTANGLE_2_INTERSECTION_H
