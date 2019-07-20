// Copyright (c) 2003-2004  
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
// Author(s)     : Sylvain Pion
 
#ifndef CGAL_KERNEL_GLOBAL_FUNCTIONS_2_H
#define CGAL_KERNEL_GLOBAL_FUNCTIONS_2_H

// Generic functions taking "user classes" as parameters, calling the
// internal functions (in *_internal*.h, namespace internal) taking a kernel as
// additional parameter, which themselves call the corresponding kernel
// functors.

#include <CGAL/user_classes.h>
#include <CGAL/Kernel/global_functions_internal_2.h>
#include <CGAL/Kernel/mpl.h>

namespace CGAL {

template < class K >
typename K::Boolean
operator==(const Point_2<K> &p, const Origin& o)
{
  return p == Point_2<K>(o);
}

template < class K >
typename K::Boolean
operator!=(const Point_2<K> &p, const Origin& o)
{
  return p != Point_2<K>(o);
}

template < class K >
inline
Angle
angle(const Vector_2<K> &u,
      const Vector_2<K> &v)
{
  return internal::angle(u, v, K());
}

template < class K >
inline
Angle
angle(const Point_2<K> &p,
      const Point_2<K> &q,
      const Point_2<K> &r)
{
  return internal::angle(p, q, r, K());
}

template < class K >
inline
Angle
angle(const Point_2<K> &p,
      const Point_2<K> &q,
      const Point_2<K> &r,
      const Point_2<K> &s)
{
  return internal::angle(p, q, r, s, K());
}

template < class K >
inline
typename K::Boolean
are_ordered_along_line(const Point_2<K> &p,
                       const Point_2<K> &q,
                       const Point_2<K> &r)
{
  return internal::are_ordered_along_line(p, q, r, K());
}

template < class K >
inline
typename K::Boolean
are_strictly_ordered_along_line(const Point_2<K> &p,
                                const Point_2<K> &q,
                                const Point_2<K> &r)
{
  return internal::are_strictly_ordered_along_line(p, q, r, K());
}

template < class K >
inline
typename K::FT
area(const Point_2<K> &p, const Point_2<K> &q, const Point_2<K> &r)
{
  return internal::area(p, q, r, K());
}

template < class K >
inline
typename K::Point_2
barycenter(const Point_2<K> &p1, const typename K::FT& w1,
           const Point_2<K> &p2)
{
  return internal::barycenter(p1, w1, p2, K());
}

template < class K >
inline
typename K::Point_2
barycenter(const Point_2<K> &p1, const typename K::FT& w1,
           const Point_2<K> &p2, const typename K::FT& w2)
{
  return internal::barycenter(p1, w1, p2, w2, K());
}

template < class K >
inline
typename K::Point_2
barycenter(const Point_2<K> &p1, const typename K::FT& w1,
           const Point_2<K> &p2, const typename K::FT& w2,
           const Point_2<K> &p3)
{
  return internal::barycenter(p1, w1, p2, w2, p3, K());
}

template < class K >
inline
typename K::Point_2
barycenter(const Point_2<K> &p1, const typename K::FT& w1,
           const Point_2<K> &p2, const typename K::FT& w2,
           const Point_2<K> &p3, const typename K::FT& w3)
{
  return internal::barycenter(p1, w1, p2, w2, p3, w3, K());
}

template < class K >
inline
typename K::Point_2
barycenter(const Point_2<K> &p1, const typename K::FT& w1,
           const Point_2<K> &p2, const typename K::FT& w2,
           const Point_2<K> &p3, const typename K::FT& w3,
           const Point_2<K> &p4)
{
  return internal::barycenter(p1, w1, p2, w2, p3, w3, p4, K());
}

template < class K >
inline
typename K::Point_2
barycenter(const Point_2<K> &p1, const typename K::FT& w1,
           const Point_2<K> &p2, const typename K::FT& w2,
           const Point_2<K> &p3, const typename K::FT& w3,
           const Point_2<K> &p4, const typename K::FT& w4)
{
  return internal::barycenter(p1, w1, p2, w2, p3, w3, p4, w4, K());
}

template <typename K>
inline
typename K::Line_2
bisector(const Point_2<K> &p, const Point_2<K> &q)
{
  return internal::bisector(p, q, K());
}

template <typename K>
inline
typename K::Line_2
bisector(const Line_2<K> &l1, const Line_2<K> &l2)
{
  return internal::bisector(l1, l2, K());
}

template < class K >
inline
typename K::Point_2
centroid(const Point_2<K> &p,
         const Point_2<K> &q,
         const Point_2<K> &r)
{
  return internal::centroid(p, q, r, K());
}

template < class K >
inline
typename K::Point_2
centroid(const Triangle_2<K> &t)
{
  return internal::centroid(t, K());
}

template < class K >
inline
typename K::Point_2
centroid(const Point_2<K> &p,
         const Point_2<K> &q,
         const Point_2<K> &r,
         const Point_2<K> &s)
{
  return internal::centroid(p, q, r, s, K());
}

template < class K >
inline
typename K::Point_2
circumcenter(const Point_2<K> &p,
             const Point_2<K> &q)
{
  return internal::circumcenter(p, q, K());
}

template < class K >
inline
typename K::Point_2
circumcenter(const Point_2<K> &p,
             const Point_2<K> &q,
             const Point_2<K> &r)
{
  return internal::circumcenter(p, q, r, K());
}

template < class K >
inline
typename K::Point_2
circumcenter(const Triangle_2<K> &t)
{
  return internal::circumcenter(t, K());
}

template < class K >
inline
typename K::Boolean
collinear(const Point_2<K> &p, const Point_2<K> &q, const Point_2<K> &r)
{
  return internal::collinear(p, q, r, K());
}

template < class K >
inline
typename K::Boolean
collinear_are_ordered_along_line(const Point_2<K> &p,
                                 const Point_2<K> &q,
                                 const Point_2<K> &r)
{
  return internal::collinear_are_ordered_along_line(p, q, r, K());
}

template < class K >
inline
typename K::Boolean
collinear_are_strictly_ordered_along_line(const Point_2<K> &p,
                                          const Point_2<K> &q,
                                          const Point_2<K> &r)
{
  return internal::collinear_are_strictly_ordered_along_line(p, q, r, K());
}

template < typename K >
inline
typename K::Comparison_result
compare_angle_with_x_axis(const Direction_2<K>& d1,
                          const Direction_2<K>& d2)
{
  return internal::compare_angle_with_x_axis(d1, d2, K());
}

template <class K >
inline
typename K::Comparison_result
compare_distance_to_point(const Point_2<K>& p,
                          const Point_2<K>& q,
                          const Point_2<K>& r)
{
  return internal::compare_distance_to_point(p, q, r, K());
}

template <class K >
inline
typename K::Comparison_result
compare_power_distance(const Point_2<K> &r,
                       const Weighted_point_2<K> &p,
                       const Weighted_point_2<K> &q)
{
  return internal::compare_power_distance(r, p, q, K());
}

template <class K >
inline
typename K::Comparison_result
compare_squared_distance(const Point_2<K>& p,
                         const Point_2<K>& q,
                         const typename K::FT& d2)
{
  return internal::compare_squared_distance(p, q, d2, K());
}

template <class K>
inline
typename K::Comparison_result
compare_signed_distance_to_line(const Point_2<K>& p,
				const Point_2<K>& q,
				const Point_2<K>& r,
				const Point_2<K>& s)
{
  return internal::compare_signed_distance_to_line(p, q, r, s, K());
}

template <class K>
inline
typename K::Comparison_result
compare_signed_distance_to_line(const Line_2<K>& l,
				const Point_2<K>& p,
				const Point_2<K>& q)
{
  return internal::compare_signed_distance_to_line(l, p, q, K());
}

/* FIXME : Undocumented, obsolete...
template < class K >
inline
typename K::Comparison_result
compare_lexicographically_xy(const Point_2<K> &p,
                             const Point_2<K> &q)
{
  return K().compare_xy_2_object()(p, q);
}
*/

template < class K >
inline
typename K::Comparison_result
compare_slope(const Line_2<K> &l1, const Line_2<K> &l2)
{
  return internal::compare_slope(l1, l2, K());
}

template < class K >
inline
typename K::Comparison_result
compare_slope(const Segment_2<K> &s1, const Segment_2<K> &s2)
{
  return internal::compare_slope(s1, s2, K());
}


#ifndef CGAL_NO_DEPRECATED_CODE
// kept for backward compatibility
template < class K >
CGAL_DEPRECATED_MSG("This function is deprecated. CGAL::compare_slope() should be used instead")
inline
typename K::Comparison_result
compare_slopes(const Line_2<K> &l1, const Line_2<K> &l2)
{
  return internal::compare_slope(l1, l2, K());
}

// kept for backward compatibility
template < class K >
CGAL_DEPRECATED_MSG("This function is deprecated. CGAL::compare_slope() should be used instead")
inline
typename K::Comparison_result
compare_slopes(const Segment_2<K> &s1, const Segment_2<K> &s2)
{
  return internal::compare_slope(s1, s2, K());
}
#endif

template < class K >
inline
typename K::Comparison_result
compare_x(const Point_2<K> &p, const Point_2<K> &q)
{
  return internal::compare_x(p, q, K());
}

template < class K >
inline
typename K::Comparison_result
compare_x(const Point_2<K>& p,
          const Line_2<K>& l1,
          const Line_2<K>& l2)
{
  return internal::compare_x(p, l1, l2, K());
}

template < class K >
inline
typename K::Comparison_result
compare_x(const Line_2<K> &l,
          const Line_2<K> &h1,
          const Line_2<K> &h2)
{
  return internal::compare_x(l, h1, h2, K());
}

template < class K >
inline
typename K::Comparison_result
compare_x(const Line_2<K> &l1,
          const Line_2<K> &h1,
          const Line_2<K> &l2,
          const Line_2<K> &h2)
{
  return internal::compare_x(l1, h1, l2, h2, K());
}

template < class K >
inline
typename K::Comparison_result
compare_x_at_y(const Point_2<K>& p, const Line_2<K>& h)
{
  return internal::compare_x_at_y(p, h, K());
}

/* Undocumented
template < class K >
inline
typename K::Comparison_result
compare_x_at_y(const Point_2<K>& p, const Segment_2<K>& s)
{
  return internal::compare_x_at_y(p, s, K());
}
*/

template < class K >
inline
typename K::Comparison_result
compare_x_at_y(const Point_2<K> &p,
               const Line_2<K> &h1,
               const Line_2<K> &h2)
{
  return internal::compare_x_at_y(p, h1, h2, K());
}

template < class K >
inline
typename K::Comparison_result
compare_x_at_y(const Line_2<K> &l1,
               const Line_2<K> &l2,
               const Line_2<K> &h)
{
  return internal::compare_x_at_y(l1, l2, h, K());
}

template < class K >
inline
typename K::Comparison_result
compare_x_at_y(const Line_2<K> &l1,
               const Line_2<K> &l2,
               const Line_2<K> &h1,
               const Line_2<K> &h2)
{
  return internal::compare_x_at_y(l1, l2, h1, h2, K());
}

template < class K >
inline
typename K::Comparison_result
compare_xy(const Point_2<K> &p, const Point_2<K> &q)
{
  return internal::compare_xy(p, q, K());
}

template < class K >
inline
typename K::Comparison_result
compare_lexicographically(const Point_2<K> &p, const Point_2<K> &q)
{
  return internal::compare_xy(p, q, K());
}

template < class K >
inline
typename K::Comparison_result
compare_y(const Point_2<K> &p, const Point_2<K> &q)
{
  return internal::compare_y(p, q, K());
}

template < class K >
inline
typename K::Comparison_result
compare_y(const Point_2<K> &p,
          const Line_2<K> &l1,
          const Line_2<K> &l2)
{
  return internal::compare_y(p, l1, l2, K());
}

template < class K >
inline
typename K::Comparison_result
compare_y(const Line_2<K> &l1,
          const Line_2<K> &l2,
          const Line_2<K> &h1,
          const Line_2<K> &h2)
{
  return internal::compare_y(l1, l2, h1, h2, K());
}

template < class K >
inline
typename K::Comparison_result
compare_y(const Line_2<K> &l,
          const Line_2<K> &h1,
          const Line_2<K> &h2)
{
  return internal::compare_y(l, h1, h2, K());
}

template < class K >
inline
typename K::Comparison_result
compare_y_at_x(const Point_2<K> &p, const Segment_2<K> &s)
{
  return internal::compare_y_at_x(p, s, K());
}

template < class K >
inline
typename K::Comparison_result
compare_y_at_x(const Point_2<K> &p,
               const Segment_2<K> &s1,
               const Segment_2<K> &s2)
{
  return internal::compare_y_at_x(p, s1, s2, K());
}

template < class K >
inline
typename K::Comparison_result
compare_y_at_x(const Point_2<K> &p, const Line_2<K> &h)
{
  return internal::compare_y_at_x(p, h, K());
}  

template < class K >
inline
typename K::Comparison_result
compare_y_at_x(const Point_2<K> &p,
               const Line_2<K> &h1,
               const Line_2<K> &h2)
{
  return internal::compare_y_at_x(p, h1, h2, K());
}

template < class K >
inline
typename K::Comparison_result
compare_y_at_x(const Line_2<K> &l1,
               const Line_2<K> &l2,
               const Line_2<K> &h)
{
  return internal::compare_y_at_x(l1, l2, h, K());
}

template < class K >
inline
typename K::Comparison_result
compare_y_at_x(const Line_2<K> &l1,
               const Line_2<K> &l2,
               const Line_2<K> &h1,
               const Line_2<K> &h2)
{
  return internal::compare_y_at_x(l1, l2, h1, h2, K());
}

template < class K >
inline
typename K::Comparison_result
compare_yx(const Point_2<K> &p, const Point_2<K> &q)
{
  return internal::compare_yx(p, q, K());
}

template < class K >
inline
typename K::FT
determinant(const Vector_2<K> &v0, const Vector_2<K> &v1)
{
  return internal::determinant(v0, v1, K());
}

template <class K>
inline
typename K::Boolean
has_larger_distance_to_point(const Point_2<K>& p,
			     const Point_2<K>& q,
			     const Point_2<K>& r)
{
  return internal::has_larger_distance_to_point(p, q, r, K());
}

template <class K>
inline
typename K::Boolean
has_smaller_distance_to_point(const Point_2<K>& p,
                              const Point_2<K>& q,
                              const Point_2<K>& r)
{
  return internal::has_smaller_distance_to_point(p, q, r, K());
}

template <class K>
inline
typename K::Boolean
has_smaller_signed_distance_to_line(const Line_2<K>& l,
                                    const Point_2<K>& p,
                                    const Point_2<K>& q)
{
  return internal::has_smaller_signed_distance_to_line(l, p, q, K());
}

template <class K>
inline
typename K::Boolean
has_larger_signed_distance_to_line(const Line_2<K>& l,
				   const Point_2<K>& p,
				   const Point_2<K>& q)
{
  return internal::has_larger_signed_distance_to_line(l, p, q, K());
}

template <class K>
inline
typename K::Boolean
has_larger_signed_distance_to_line(const Point_2<K>& p,
				   const Point_2<K>& q,
				   const Point_2<K>& r,
				   const Point_2<K>& s)
{
  return internal::has_larger_signed_distance_to_line(p, q, r, s, K());
}

template <class K>
inline
typename K::Boolean
has_smaller_signed_distance_to_line(const Point_2<K>& p,
                                    const Point_2<K>& q,
                                    const Point_2<K>& r,
                                    const Point_2<K>& s)
{
  return internal::has_smaller_signed_distance_to_line(p, q, r, s, K());
}

template < class K >
inline
typename K::Boolean
left_turn(const Point_2<K> &p, const Point_2<K> &q, const Point_2<K> &r)
{
  return internal::left_turn(p, q, r, K());
}

template < class K >
inline
typename K::Boolean
less_x(const Point_2<K> &p, const Point_2<K> &q)
{
  return internal::less_x(p, q, K());
}

template < class K >
inline
typename K::Boolean
less_y(const Point_2<K> &p, const Point_2<K> &q)
{
  return internal::less_y(p, q, K());
}

template < class K >
inline
typename K::Boolean
lexicographically_xy_larger(const Point_2<K> &p, const Point_2<K> &q)
{
  return internal::lexicographically_xy_larger(p, q, K());
}

template < class K >
inline
typename K::Boolean
lexicographically_xy_larger_or_equal(const Point_2<K> &p, const Point_2<K> &q)
{
  return internal::lexicographically_xy_larger_or_equal(p, q, K());
}

template < class K >
inline
typename K::Boolean
lexicographically_xy_smaller(const Point_2<K> &p, const Point_2<K> &q)
{
  return internal::lexicographically_xy_smaller(p, q, K());
}

template < class K >
inline
typename K::Boolean
lexicographically_xy_smaller_or_equal(const Point_2<K> &p,
                                      const Point_2<K> &q)
{
  return internal::lexicographically_xy_smaller_or_equal(p, q, K());
}

template < class K >
inline
typename K::Boolean
lexicographically_yx_smaller(const Point_2<K> &p, const Point_2<K> &q)
{
  return internal::lexicographically_yx_smaller(p, q, K());
}

template < class K >
inline
typename K::Boolean
lexicographically_yx_smaller_or_equal(const Point_2<K> &p, const Point_2<K> &q)
{
  return internal::lexicographically_yx_smaller_or_equal(p, q, K());
}

// FIXME : Undocumented
template < class K >
inline
typename K::Boolean
lexicographically_yx_larger(const Point_2<K> &p, const Point_2<K> &q)
{
  return internal::lexicographically_yx_larger(p, q, K());
}

// FIXME : Undocumented
template < class K >
inline
typename K::Boolean
lexicographically_yx_larger_or_equal(const Point_2<K> &p, const Point_2<K> &q)
{
  return internal::lexicographically_yx_larger_or_equal(p, q, K());
}

template < class K >
typename K::FT
l_infinity_distance(const Point_2<K> &p, const Point_2<K> &q)
{
  return internal::l_infinity_distance(p,q, K());
}

template < class K >
inline
typename K::Point_2
midpoint(const Point_2<K> &p, const Point_2<K> &q)
{
  return internal::midpoint(p, q, K());
}

template < class K >
inline
typename K::Point_2
max_vertex(const Iso_rectangle_2<K> &ir)
{
  return internal::max_vertex(ir, K());
}

template < class K >
inline
typename K::Point_2
min_vertex(const Iso_rectangle_2<K> &ir)
{
  return internal::min_vertex(ir, K());
}

// FIXME TODO : What do we do with the operators ?
// They have no counter part with the kernel argument...
template < class K >
inline
typename K::Boolean
operator<(const Direction_2<K>& d1, const Direction_2<K>& d2)
{ return compare_angle_with_x_axis(d1, d2) == SMALLER; }

template < class K >
inline
typename K::Boolean
operator>(const Direction_2<K>& d1, const Direction_2<K>& d2)
{ return compare_angle_with_x_axis(d1, d2) == LARGER; }

template < class K >
inline
typename K::Boolean
operator>=(const Direction_2<K>& d1, const Direction_2<K>& d2)
{ return compare_angle_with_x_axis(d1, d2) != SMALLER; }

template < class K >
inline
typename K::Boolean
operator<=(const Direction_2<K>& d1, const Direction_2<K>& d2)
{ return compare_angle_with_x_axis(d1, d2) != LARGER; }

template < class K >
inline
typename K::Boolean
operator==(const Point_2<K>& p, const Point_2<K>& q)
{ return K().equal_2_object()(p, q); }

template < class K >
inline
typename K::Boolean
operator!=(const Point_2<K>& p, const Point_2<K>& q)
{ return ! (p == q); }

template < class K >
inline
typename K::Boolean
operator<(const Point_2<K>& p, const Point_2<K>& q)
{ return K().less_xy_2_object()(p, q); }

template < class K >
inline
typename K::Boolean
operator>(const Point_2<K>& p, const Point_2<K>& q)
{ return K().less_xy_2_object()(q, p); }

template < class K >
inline
typename K::Boolean
operator<=(const Point_2<K>& p, const Point_2<K>& q)
{ return ! K().less_xy_2_object()(q, p); }

template < class K >
inline
typename K::Boolean
operator>=(const Point_2<K>& p, const Point_2<K>& q)
{ return ! K().less_xy_2_object()(p, q); }

template < class K >
inline
typename K::Boolean
operator==(const Vector_2<K>& v, const Vector_2<K>& w)
{ return K().equal_2_object()(v, w); }

template < class K >
inline
typename K::Boolean
operator!=(const Vector_2<K>& v, const Vector_2<K>& w)
{ return ! (v == w); }

template < class K >
inline
typename K::Vector_2
operator*(const typename K::FT &c, const Vector_2<K> &w)
{
  return K().construct_scaled_vector_2_object()(w, c);
}

template < class K >
inline
typename K::Vector_2
operator*(const Vector_2<K> &w, const typename K::FT &c)
{
  return K().construct_scaled_vector_2_object()(w, c);
}

template < class K >
inline
typename K::Vector_2
operator*(const typename First_if_different<typename K::RT,
                                            typename K::FT>::Type &c,
          const Vector_2<K> &w)
{
  return K().construct_scaled_vector_2_object()(w, c);
}

template < class K >
inline
typename K::Vector_2
operator*(const Vector_2<K> &w,
          const typename First_if_different<typename K::RT,
                                            typename K::FT>::Type &c)
{
  return K().construct_scaled_vector_2_object()(w, c);
}

template < class K >
inline
typename K::FT
operator*(const Vector_2<K> &v, const Vector_2<K> &w)
{
  return K().compute_scalar_product_2_object()(v, w);
}

template < class K >
inline
typename K::Point_2
operator+(const Point_2<K> &p, const Vector_2<K> &v)
{
  return K().construct_translated_point_2_object()(p, v);
}

template < class K >
inline
typename K::Point_2
operator+(const Origin &o, const Vector_2<K> &v)
{
  return K().construct_translated_point_2_object()(o, v);
}

template < class K >
inline
typename K::Point_2
operator-(const Point_2<K> &p, const Vector_2<K> &v)
{
  return K().construct_translated_point_2_object()
                (p, K().construct_opposite_vector_2_object()(v));
}

template < class K >
inline
typename K::Point_2
operator-(const Origin &o, const Vector_2<K> &v)
{
  return K().construct_translated_point_2_object()
               (o, K().construct_opposite_vector_2_object()(v));
}

template < class K >
inline
typename K::Vector_2
operator-(const Point_2<K> &p, const Point_2<K> &q)
{
  return K().construct_vector_2_object()(q, p);
}

template < class K >
inline
typename K::Vector_2
operator-(const Point_2<K> &p, const Origin &o)
{
  return K().construct_vector_2_object()(o, p);
}

template < class K >
inline
typename K::Vector_2
operator-(const Origin &o, const Point_2<K> &q)
{
  return K().construct_vector_2_object()(q, o);
}

template <typename K>
inline
typename K::Orientation
orientation(const Point_2<K> &p, const Point_2<K> &q, const Point_2<K> &r)
{
  return internal::orientation(p, q, r, K());
}

template <typename K>
inline
typename K::Orientation
orientation(const Vector_2<K> &u, const Vector_2<K> &v)
{
  return internal::orientation(u, v, K());
}

// parallel() functions are in global_functions.h

template <class K >
inline
typename K::FT
power_product(const Weighted_point_2<K> &p,
              const Weighted_point_2<K> &q)
{
  return internal::power_product(p, q, K());
}

template <class K >
inline
typename K::Bounded_side
power_side_of_bounded_power_circle(const Weighted_point_2<K> &p,
                                   const Weighted_point_2<K> &q)
{
  return internal::power_side_of_bounded_power_circle(p, q, K());
}

template <class K >
inline
typename K::Bounded_side
power_side_of_bounded_power_circle(const Weighted_point_2<K> &p,
                                   const Weighted_point_2<K> &q,
                                   const Weighted_point_2<K> &r)
{
  return internal::power_side_of_bounded_power_circle(p, q, r, K());
}

template <class K >
inline
typename K::Bounded_side
power_side_of_bounded_power_circle(const Weighted_point_2<K> &p,
                                   const Weighted_point_2<K> &q,
                                   const Weighted_point_2<K> &r,
                                   const Weighted_point_2<K> &s)
{
  return internal::power_side_of_bounded_power_circle(p, q, r, s, K());
}

template <typename K>
inline
typename K::Orientation
power_side_of_oriented_power_circle(const Weighted_point_2<K> &p,
                                    const Weighted_point_2<K> &q)
{
  return internal::power_side_of_oriented_power_circle(p, q, K());
}

template <typename K>
inline
typename K::Orientation
power_side_of_oriented_power_circle(const Weighted_point_2<K> &p,
                                    const Weighted_point_2<K> &q,
                                    const Weighted_point_2<K> &r)
{
  return internal::power_side_of_oriented_power_circle(p, q, r, K());
}

template <typename K>
inline
typename K::Orientation
power_side_of_oriented_power_circle(const Weighted_point_2<K> &p,
                                    const Weighted_point_2<K> &q,
                                    const Weighted_point_2<K> &r,
                                    const Weighted_point_2<K> &s)
{
  return internal::power_side_of_oriented_power_circle(p, q, r, s, K());
}

template <class K>
inline
typename K::Line_2
radical_axis(const Weighted_point_2<K> &p,
             const Weighted_point_2<K> &q)
{
  return internal::radical_axis(p, q, K());
}

template <class K>
inline
typename K::Line_2
radical_line(const Circle_2<K> &s1,
             const Circle_2<K> &s2)
{
  return K().construct_radical_line_2_object()(s1,s2);
}

template <typename K>
inline
typename K::Boolean
right_turn(const Point_2<K> &p, const Point_2<K> &q, const Point_2<K> &r)
{
  return internal::right_turn(p, q, r, K());
}

template < class K >
inline
typename K::FT
scalar_product(const Vector_2<K> &v, const Vector_2<K> &w)
{
  return K().compute_scalar_product_2_object()(v, w);
}

template <class K>
inline
typename K::Bounded_side
side_of_bounded_circle(const Point_2<K> &p,
                       const Point_2<K> &q,
                       const Point_2<K> &r,
                       const Point_2<K> &t)
{
  return internal::side_of_bounded_circle(p, q, r, t, K());
}

template <class K>
inline
typename K::Bounded_side
side_of_bounded_circle(const Point_2<K> &p,
                       const Point_2<K> &q,
                       const Point_2<K> &r)
{
  return internal::side_of_bounded_circle(p, q, r, K());
}

template <class K>
inline
typename K::Oriented_side
side_of_oriented_circle(const Point_2<K> &p,
                        const Point_2<K> &q,
                        const Point_2<K> &r,
                        const Point_2<K> &t)
{
  return internal::side_of_oriented_circle(p, q, r, t, K());
}

template < class K >
inline
typename K::FT
squared_radius(const Point_2<K> &p)
{
  return internal::squared_radius(p, K());
}

template < class K >
inline
typename K::FT
squared_radius(const Point_2<K> &p, const Point_2<K> &q)
{
  return internal::squared_radius(p, q, K());
}

template < class K >
CGAL_KERNEL_INLINE
typename K::FT
squared_radius(const Point_2<K>& p, const Point_2<K>& q, const Point_2<K>& r)
{
  return internal::squared_radius(p, q, r, K());
}

template < class K >
inline
typename K::FT
squared_radius_smallest_orthogonal_circle(const Weighted_point_2<K> &p)
{
  return internal::squared_radius_smallest_orthogonal_circle(p, K());
}

template < class K >
inline
typename K::FT
squared_radius_smallest_orthogonal_circle(const Weighted_point_2<K> &p,
                                          const Weighted_point_2<K> &q)
{
  return internal::squared_radius_smallest_orthogonal_circle(p, q, K());
}

template < class K >
inline
typename K::FT
squared_radius_smallest_orthogonal_circle(const Weighted_point_2<K> &p,
                                          const Weighted_point_2<K> &q,
                                          const Weighted_point_2<K> &r)
{
  return internal::squared_radius_smallest_orthogonal_circle(p, q, r, K());
}

template < class K >
inline
typename K::Point_2
weighted_circumcenter(const Weighted_point_2<K> &p,
                      const Weighted_point_2<K> &q,
                      const Weighted_point_2<K> &r)
{
  return internal::weighted_circumcenter(p, q, r, K());
}

template < class K >
inline
typename K::Boolean
x_equal(const Point_2<K> &p, const Point_2<K> &q)
{
  return internal::x_equal(p, q, K());
}

template < class K >
inline
typename K::Boolean
y_equal(const Point_2<K> &p, const Point_2<K> &q)
{
  return internal::y_equal(p, q, K());
}

} //namespace CGAL

#endif  // CGAL_KERNEL_GLOBAL_FUNCTIONS_2_H
