
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


#ifndef CGAL_LINE_2_ISO_RECTANGLE_2_INTERSECTION_H
#define CGAL_LINE_2_ISO_RECTANGLE_2_INTERSECTION_H

#include <CGAL/Line_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/Object.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

template <class K>
class Line_2_Iso_rectangle_2_pair {
public:
    enum Intersection_results {NO, POINT, SEGMENT};
    Line_2_Iso_rectangle_2_pair() ;
    Line_2_Iso_rectangle_2_pair(typename K::Line_2 const *pt,
                            typename K::Iso_rectangle_2 const *iso);
    ~Line_2_Iso_rectangle_2_pair() {}

  Intersection_results intersection_type() const;

    bool                intersection(typename K::Point_2 &result) const;
    bool                intersection(typename K::Segment_2 &result) const;
protected:
    mutable bool                        _known;
    mutable Intersection_results        _result;
    mutable typename K::FT              _min, _max;
    typename K::Point_2             _ref_point;
    typename K::Vector_2            _dir;
    typename K::Point_2             _isomin;
    typename K::Point_2             _isomax;
};

template <class K>
inline bool do_intersect(const typename CGAL_WRAP(K)::Line_2 &p1,
			 const typename CGAL_WRAP(K)::Iso_rectangle_2 &p2,
			 const K&)
{
    typedef Line_2_Iso_rectangle_2_pair<K> pair_t;
    pair_t pair(&p1, &p2);
    return pair.intersection_type() != pair_t::NO;
}

template <class K>
inline bool do_intersect(const typename CGAL_WRAP(K)::Iso_rectangle_2 &p2,
			 const typename CGAL_WRAP(K)::Line_2 &p1,
			 const K& k)
{
  return CGALi::do_intersect(p1, p2, k);
}


template <class K>
Line_2_Iso_rectangle_2_pair<K>::
Line_2_Iso_rectangle_2_pair()
{
    _known = false;
}

template <class K>
Line_2_Iso_rectangle_2_pair<K>::
Line_2_Iso_rectangle_2_pair(typename K::Line_2 const *line,
                            typename K::Iso_rectangle_2 const *iso)
  : _known(false),
    _ref_point(line->point()),
    _dir(line->direction().to_vector()),
    _isomin(iso->min()),
    _isomax(iso->max())
{}


template <class K>
typename Line_2_Iso_rectangle_2_pair<K>::Intersection_results
Line_2_Iso_rectangle_2_pair<K>::intersection_type() const
{
    typedef typename K::Line_2 line_t;
    if (_known)
        return _result;
// The non const this pointer is used to cast away const.
    _known = true;
    typedef typename K::FT FT;
    typedef typename K::RT RT;
    bool all_values = true;

    typename K::Construct_cartesian_const_iterator_2 construct_cccit;
    typename K::Cartesian_const_iterator_2 ref_point_it = construct_cccit(_ref_point);
    typename K::Cartesian_const_iterator_2 end = construct_cccit(_ref_point, 0);
    typename K::Cartesian_const_iterator_2 isomin_it = construct_cccit(_isomin);
    typename K::Cartesian_const_iterator_2 isomax_it = construct_cccit(_isomax);

    for (unsigned int i=0; ref_point_it != end; ++i, ++ref_point_it, ++isomin_it, ++isomax_it) {

        if (_dir.homogeneous(i) == RT(0)) {
            if (*ref_point_it < *isomin_it) {
                _result = NO;
                return NO;
            }
            if (*ref_point_it > *isomax_it) {
                _result = NO;
                return NO;
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
            if (all_values) {
                _min = newmin;
                _max = newmax;
            } else {
                if (newmin > _min)
                    _min = newmin;
                if (newmax < _max)
                    _max = newmax;
                if (_max < _min) {
                    _result = NO;
                    return NO;
                }
            }
            all_values = false;
        }
    }
    CGAL_kernel_assertion(!all_values);
    if (_max == _min) {
        _result = POINT;
        return POINT;
    }
    _result = SEGMENT;
    return SEGMENT;
}


template <class K>
bool
Line_2_Iso_rectangle_2_pair<K>::
intersection(typename K::Point_2 &result) const
{
  typename K::Construct_translated_point_2 translated_point;
  typename K::Construct_scaled_vector_2 construct_scaled_vector;

    if (!_known)
        intersection_type();
    if (_result != POINT)
        return false;
    result = translated_point(_ref_point, construct_scaled_vector(_dir, _min));
    return true;
}

template <class K>
bool
Line_2_Iso_rectangle_2_pair<K>::
intersection(typename K::Segment_2 &result) const
{
  typename K::Construct_segment_2 construct_segment_2;
  typename K::Construct_translated_point_2 translated_point;
  typename K::Construct_scaled_vector_2 construct_scaled_vector;
    if (!_known)
        intersection_type();
    if (_result != SEGMENT)
        return false;
    result = construct_segment_2(translated_point(_ref_point, construct_scaled_vector(_dir,_min)), 
				 translated_point(_ref_point, construct_scaled_vector(_dir,_max)));
    return true;
}



template <class K>
Object
intersection(const typename CGAL_WRAP(K)::Line_2 &line, 
	     const typename CGAL_WRAP(K)::Iso_rectangle_2 &iso,
	     const K&)
{
    typename K::Construct_object_2 construct_object;
    typedef Line_2_Iso_rectangle_2_pair<K> is_t;
    is_t ispair(&line, &iso);
    switch (ispair.intersection_type()) {
    case is_t::NO:
    default:
        return Object();
    case is_t::POINT: {
        typename K::Point_2 ipt;
        ispair.intersection(ipt);
        return construct_object(ipt);
    }
    case is_t::SEGMENT: {
        typename K::Segment_2 iseg;
        ispair.intersection(iseg);
        return construct_object(iseg);
    }
    }
}

template <class K>
inline
Object
intersection(const typename CGAL_WRAP(K)::Iso_rectangle_2 &iso,
	     const typename CGAL_WRAP(K)::Line_2 &line, 
	     const K& k)
{
  return CGALi::intersection(line, iso, k);
}

template <class K>
class Iso_rectangle_2_Line_2_pair
: public Line_2_Iso_rectangle_2_pair<K> {
public:
    Iso_rectangle_2_Line_2_pair(
            typename K::Iso_rectangle_2 const *iso,
            typename K::Line_2 const *line) :
                        Line_2_Iso_rectangle_2_pair<K>(line, iso) {}
};

} // namespace CGALi




template <class K>
inline bool do_intersect(
    const Line_2<K> &p1,
    const Iso_rectangle_2<K> &p2)
{
  typedef typename K::Do_intersect_2 Do_intersect;
  return Do_intersect()(p1, p2);
}

template <class K>
inline bool do_intersect(
    const Iso_rectangle_2<K> &p1,
    const Line_2<K> &p2)
{
  typedef typename K::Do_intersect_2 Do_intersect;
  return Do_intersect()(p2, p1);
}

template <class K>
Object
intersection(const Line_2<K> &line, const Iso_rectangle_2<K> &iso)
{
  typedef typename K::Intersect_2 Intersect;
  return Intersect()(line, iso);
}

template <class K>
inline Object
intersection(const Iso_rectangle_2<K> &iso, const Line_2<K> &line)
{
  typedef typename K::Intersect_2 Intersect;
    return Intersect()(line, iso);
}

CGAL_END_NAMESPACE

#endif
