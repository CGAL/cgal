
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


#ifndef CGAL_RAY_2_ISO_RECTANGLE_2_INTERSECTION_H
#define CGAL_RAY_2_ISO_RECTANGLE_2_INTERSECTION_H

#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/kernel_assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/Object.h>


CGAL_BEGIN_NAMESPACE

namespace CGALi {

template <class K>
class Ray_2_Iso_rectangle_2_pair {
public:
    enum Intersection_results {NO, POINT, SEGMENT};
    Ray_2_Iso_rectangle_2_pair() ;
    Ray_2_Iso_rectangle_2_pair(typename K::Ray_2 const *ray,
                          typename K::Iso_rectangle_2 const *rect) ;

  Intersection_results intersection_type() const;

    bool intersection( typename K::Point_2 &result) const;
    bool intersection( typename K::Segment_2 &result) const;
protected:
    mutable bool                       _known;
    mutable Intersection_results       _result;
    mutable typename K::Point_2            _ref_point;
    mutable typename K::Vector_2           _dir;
    mutable typename K::Point_2            _isomin;
    mutable typename K::Point_2            _isomax;
    mutable typename K::FT                     _min,
                               _max;
};

template <class K>
inline bool do_intersect(const typename CGAL_WRAP(K)::Ray_2 &p1,
			 const typename CGAL_WRAP(K)::Iso_rectangle_2 &p2,
			 const K& k)
{
    typedef Ray_2_Iso_rectangle_2_pair<K> pair_t;
    pair_t pair(&p1, &p2);
    return pair.intersection_type() != pair_t::NO;
}

template <class K>
inline bool do_intersect(const typename CGAL_WRAP(K)::Iso_rectangle_2 &p2,
			 const typename CGAL_WRAP(K)::Ray_2 &p1,
			 const K& k)
{
  return do_intersect(p1, p2, k);
}

template <class K>
Object
intersection(const typename CGAL_WRAP(K)::Ray_2 &ray,
	     const typename CGAL_WRAP(K)::Iso_rectangle_2 &iso,
	     const K& k)
{
    typedef Ray_2_Iso_rectangle_2_pair<K> is_t;
    is_t ispair(&ray, &iso);
    switch (ispair.intersection_type()) {
    case is_t::NO:
    default:
        return Object();
    case is_t::POINT: {
        typename K::Point_2 ipt;
        ispair.intersection(ipt);
        return make_object(ipt);
    }
    case is_t::SEGMENT: {
        typename K::Segment_2 iseg;
        ispair.intersection(iseg);
        return make_object(iseg);
    }
    }
}

template <class K>
Object
intersection(const typename CGAL_WRAP(K)::Iso_rectangle_2 &iso,
	     const typename CGAL_WRAP(K)::Ray_2 &ray,
	     const K& k)
{
  return intersection(ray, iso, k);
}


template <class K>
Ray_2_Iso_rectangle_2_pair<K>::Ray_2_Iso_rectangle_2_pair()
{
    _known = false;
}

template <class K>
Ray_2_Iso_rectangle_2_pair<K>::
Ray_2_Iso_rectangle_2_pair(
        typename K::Ray_2 const *ray,
        typename K::Iso_rectangle_2 const *iso)
{
    _known = false;
    _isomin = iso->min();
    _isomax = iso->max();
    _ref_point = ray->source();
    _dir = ray->direction().to_vector();
    _min = (typename K::FT)(0);
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
                _result = NO;
                return _result;
            }
            if ((*ref_point_it) > (*isomax_it)) {
                _result = NO;
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
                _result = NO;
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
bool Ray_2_Iso_rectangle_2_pair<K>::
intersection(typename K::Segment_2 &seg) const
{
  typedef typename K::Segment_2 Segment_2; 
  typename K::Construct_translated_point_2 translated_point;
  typename K::Construct_scaled_vector_2 construct_scaled_vector;
    if (!_known)
        intersection_type();
    if (_result != SEGMENT)
        return false;
    typename K::Point_2 p1(translated_point(_ref_point, construct_scaled_vector(_dir,_min)));
    typename K::Point_2 p2(translated_point(_ref_point, construct_scaled_vector(_dir,_max)));
    seg = Segment_2(p1, p2);
    return true;
}

template <class K> bool Ray_2_Iso_rectangle_2_pair<K>::
intersection(typename K::Point_2 &pt) const
{
  typedef typename K::Point_2 Point_2;
  typename K::Construct_translated_point_2 translated_point;
  typename K::Construct_scaled_vector_2 construct_scaled_vector;
    if (!_known)
        intersection_type();
    if (_result != POINT)
        return false;
    pt = Point_2(translated_point(_ref_point, construct_scaled_vector(_dir, _min)));
    return true;
}



template <class K>
class Iso_rectangle_2_Ray_2_pair:
          public Ray_2_Iso_rectangle_2_pair<K> {
public:
    Iso_rectangle_2_Ray_2_pair() {}
    Iso_rectangle_2_Ray_2_pair(typename K::Iso_rectangle_2 const *rect,
                               typename K::Ray_2 const *ray)
                :Ray_2_Iso_rectangle_2_pair<K> (ray, rect){}
};

} // namespace CGALi


template <class K>
inline bool do_intersect(const Iso_rectangle_2<K> &p1,
			 const Ray_2<K> &p2)
{
  typedef typename K::Do_intersect_2 Do_intersect;
  return Do_intersect()(p1, p2);
}

template <class K>
inline bool do_intersect(const Ray_2<K> &p2,
			 const Iso_rectangle_2<K> &p1)
{
  typedef typename K::Do_intersect_2 Do_intersect;
  return Do_intersect()(p1, p2);
}


template <class K>
inline Object
intersection(const Iso_rectangle_2<K>&iso, const Ray_2<K>&ray)
{
  typedef typename K::Intersect_2 Intersect;
  return Intersect()(ray, iso);
}

template <class K>
inline Object
intersection(const Ray_2<K>&ray, const Iso_rectangle_2<K>&iso)
{
  typedef typename K::Intersect_2 Intersect;
  return Intersect()(ray, iso);
}

CGAL_END_NAMESPACE

#endif // CGAL_RAY_2_iSO_RECTANGLE_2_INTERSECTION_H
