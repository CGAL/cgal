// Copyright (c) 2003-2004  Utrecht University (The Netherlands),
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
// Author(s)     : Sylvain Pion
 
#ifndef CGAL_KERNEL_GLOBAL_FUNCTIONS_2_H
#define CGAL_KERNEL_GLOBAL_FUNCTIONS_2_H

// Generic functions taking "user classes" as parameters, calling the
// internal functions (in *_internal*.h, namespace CGALi) taking a kernel as
// additional parameter, which themselves call the corresponding kernel
// functors.

#include <CGAL/user_classes.h>
#include <CGAL/Kernel/global_functions_internal_2.h>
#include <CGAL/Kernel/mpl.h>

CGAL_BEGIN_NAMESPACE

template < class K >
inline
Angle
angle(const Point_2<K> &p,
      const Point_2<K> &q,
      const Point_2<K> &r)
{
  return CGALi::angle(p, q, r, K());
}

template < class K >
inline
bool
are_ordered_along_line(const Point_2<K> &p,
                       const Point_2<K> &q,
                       const Point_2<K> &r)
{
  return CGALi::are_ordered_along_line(p, q, r, K());
}

template < class K >
inline
bool
are_strictly_ordered_along_line(const Point_2<K> &p,
                                const Point_2<K> &q,
                                const Point_2<K> &r)
{
  return CGALi::are_strictly_ordered_along_line(p, q, r, K());
}

template < class K >
inline
typename K::FT
area(const Point_2<K> &p, const Point_2<K> &q, const Point_2<K> &r)
{
  return CGALi::area(p, q, r, K());
}

template <typename K>
inline
typename K::Line_2
bisector(const Point_2<K> &p, const Point_2<K> &q)
{
  return CGALi::bisector(p, q, K());
}

template <typename K>
inline
typename K::Line_2
bisector(const Line_2<K> &l1, const Line_2<K> &l2)
{
  return CGALi::bisector(l1, l2, K());
}

template < class K >
inline
typename K::Point_2
centroid(const Point_2<K> &p,
         const Point_2<K> &q,
         const Point_2<K> &r)
{
  return CGALi::centroid(p, q, r, K());
}

template < class K >
inline
typename K::Point_2
centroid(const Point_2<K> &p,
         const Point_2<K> &q,
         const Point_2<K> &r,
         const Point_2<K> &s)
{
  return CGALi::centroid(p, q, r, s, K());
}

template < class K >
inline
typename K::Point_2
circumcenter(const Point_2<K> &p,
             const Point_2<K> &q,
             const Point_2<K> &r)
{
  return CGALi::circumcenter(p, q, r, K());
}

template < class K >
inline
typename K::Point_2
circumcenter(const Triangle_2<K> &t)
{
  return CGALi::circumcenter(t, K());
}

template < class K >
inline
bool
collinear(const Point_2<K> &p, const Point_2<K> &q, const Point_2<K> &r)
{
  return CGALi::collinear(p, q, r, K());
}

template < class K >
inline
bool
collinear_are_ordered_along_line(const Point_2<K> &p,
                                 const Point_2<K> &q,
                                 const Point_2<K> &r)
{
  return CGALi::collinear_are_ordered_along_line(p, q, r, K());
}

template < class K >
inline
bool
collinear_are_strictly_ordered_along_line(const Point_2<K> &p,
                                          const Point_2<K> &q,
                                          const Point_2<K> &r)
{
  return CGALi::collinear_are_strictly_ordered_along_line(p, q, r, K());
}

template < typename K >
inline
Comparison_result
compare_angle_with_x_axis(const Direction_2<K>& d1,
                          const Direction_2<K>& d2)
{
  return CGALi::compare_angle_with_x_axis(d1, d2, K());
}

template <class K >
inline
Comparison_result
compare_distance_to_point(const Point_2<K>& p,
                          const Point_2<K>& q,
                          const Point_2<K>& r)
{
  return CGALi::compare_distance_to_point(p, q, r, K());
}

template <class K>
inline
Comparison_result
compare_signed_distance_to_line(const Point_2<K>& p,
				const Point_2<K>& q,
				const Point_2<K>& r,
				const Point_2<K>& s)
{
  return CGALi::compare_signed_distance_to_line(p, q, r, s, K());
}

template <class K>
inline
Comparison_result
compare_signed_distance_to_line(const Line_2<K>& l,
				const Point_2<K>& p,
				const Point_2<K>& q)
{
  return CGALi::compare_signed_distance_to_line(l, p, q, K());
}

/* FIXME : Undocumented, obsolete...
template < class K >
inline
Comparison_result
compare_lexicographically_xy(const Point_2<K> &p,
                             const Point_2<K> &q)
{
  return K().compare_xy_2_object()(p, q);
}
*/

template < class K >
inline
Comparison_result
compare_slopes(const Line_2<K> &l1, const Line_2<K> &l2)
{
  return CGALi::compare_slopes(l1, l2, K());
}

template < class K >
inline
Comparison_result
compare_slopes(const Segment_2<K> &s1, const Segment_2<K> &s2)
{
  return CGALi::compare_slopes(s1, s2, K());
}

template < class K >
inline
Comparison_result
compare_x(const Point_2<K> &p, const Point_2<K> &q)
{
  return CGALi::compare_x(p, q, K());
}

template < class K >
inline
Comparison_result
compare_x(const Point_2<K>& p,
          const Line_2<K>& l1,
          const Line_2<K>& l2)
{
  return CGALi::compare_x(p, l1, l2, K());
}

template < class K >
inline
Comparison_result
compare_x(const Line_2<K> &l,
          const Line_2<K> &h1,
          const Line_2<K> &h2)
{
  return CGALi::compare_x(l, h1, h2, K());
}

template < class K >
inline
Comparison_result
compare_x(const Line_2<K> &l1,
          const Line_2<K> &h1,
          const Line_2<K> &l2,
          const Line_2<K> &h2)
{
  return CGALi::compare_x(l1, h1, l2, h2, K());
}

template < class K >
inline
Comparison_result
compare_x_at_y(const Point_2<K>& p, const Line_2<K>& h)
{
  return CGALi::compare_x_at_y(p, h, K());
}

/* Undocumented
template < class K >
inline
Comparison_result
compare_x_at_y(const Point_2<K>& p, const Segment_2<K>& s)
{
  return CGALi::compare_x_at_y(p, s, K());
}
*/

template < class K >
inline
Comparison_result
compare_x_at_y(const Point_2<K> &p,
               const Line_2<K> &h1,
               const Line_2<K> &h2)
{
  return CGALi::compare_x_at_y(p, h1, h2, K());
}

template < class K >
inline
Comparison_result
compare_x_at_y(const Line_2<K> &l1,
               const Line_2<K> &l2,
               const Line_2<K> &h)
{
  return CGALi::compare_x_at_y(l1, l2, h, K());
}

template < class K >
inline
Comparison_result
compare_x_at_y(const Line_2<K> &l1,
               const Line_2<K> &l2,
               const Line_2<K> &h1,
               const Line_2<K> &h2)
{
  return CGALi::compare_x_at_y(l1, l2, h1, h2, K());
}

template < class K >
inline
Comparison_result
compare_xy(const Point_2<K> &p, const Point_2<K> &q)
{
  return CGALi::compare_xy(p, q, K());
}

template < class K >
inline
Comparison_result
compare_y(const Point_2<K> &p, const Point_2<K> &q)
{
  return CGALi::compare_y(p, q, K());
}

template < class K >
inline
Comparison_result
compare_y(const Point_2<K> &p,
          const Line_2<K> &l1,
          const Line_2<K> &l2)
{
  return CGALi::compare_y(p, l1, l2, K());
}

template < class K >
inline
Comparison_result
compare_y(const Line_2<K> &l1,
          const Line_2<K> &l2,
          const Line_2<K> &h1,
          const Line_2<K> &h2)
{
  return CGALi::compare_y(l1, l2, h1, h2, K());
}

template < class K >
inline
Comparison_result
compare_y(const Line_2<K> &l,
          const Line_2<K> &h1,
          const Line_2<K> &h2)
{
  return CGALi::compare_y(l, h1, h2, K());
}

template < class K >
inline
Comparison_result
compare_y_at_x(const Point_2<K> &p, const Segment_2<K> &s)
{
  return CGALi::compare_y_at_x(p, s, K());
}

template < class K >
inline
Comparison_result
compare_y_at_x(const Point_2<K> &p,
               const Segment_2<K> &s1,
               const Segment_2<K> &s2)
{
  return CGALi::compare_y_at_x(p, s1, s2, K());
}

template < class K >
inline
Comparison_result
compare_y_at_x(const Point_2<K> &p, const Line_2<K> &h)
{
  return CGALi::compare_y_at_x(p, h, K());
}  

template < class K >
inline
Comparison_result
compare_y_at_x(const Point_2<K> &p,
               const Line_2<K> &h1,
               const Line_2<K> &h2)
{
  return CGALi::compare_y_at_x(p, h1, h2, K());
}

template < class K >
inline
Comparison_result
compare_y_at_x(const Line_2<K> &l1,
               const Line_2<K> &l2,
               const Line_2<K> &h)
{
  return CGALi::compare_y_at_x(l1, l2, h, K());
}

template < class K >
inline
Comparison_result
compare_y_at_x(const Line_2<K> &l1,
               const Line_2<K> &l2,
               const Line_2<K> &h1,
               const Line_2<K> &h2)
{
  return CGALi::compare_y_at_x(l1, l2, h1, h2, K());
}

template <class K>
inline
bool
has_larger_distance_to_point(const Point_2<K>& p,
			     const Point_2<K>& q,
			     const Point_2<K>& r)
{
  return CGALi::has_larger_distance_to_point(p, q, r, K());
}

template <class K>
inline
bool
has_smaller_distance_to_point(const Point_2<K>& p,
                              const Point_2<K>& q,
                              const Point_2<K>& r)
{
  return CGALi::has_smaller_distance_to_point(p, q, r, K());
}

template <class K>
inline
bool
has_smaller_signed_distance_to_line(const Line_2<K>& l,
                                    const Point_2<K>& p,
                                    const Point_2<K>& q)
{
  return CGALi::has_smaller_signed_distance_to_line(l, p, q, K());
}

template <class K>
inline
bool
has_larger_signed_distance_to_line(const Line_2<K>& l,
				   const Point_2<K>& p,
				   const Point_2<K>& q)
{
  return CGALi::has_larger_signed_distance_to_line(l, p, q, K());
}

template <class K>
inline
bool
has_larger_signed_distance_to_line(const Point_2<K>& p,
				   const Point_2<K>& q,
				   const Point_2<K>& r,
				   const Point_2<K>& s)
{
  return CGALi::has_larger_signed_distance_to_line(p, q, r, s, K());
}

template <class K>
inline
bool
has_smaller_signed_distance_to_line(const Point_2<K>& p,
                                    const Point_2<K>& q,
                                    const Point_2<K>& r,
                                    const Point_2<K>& s)
{
  return CGALi::has_smaller_signed_distance_to_line(p, q, r, s, K());
}

template < class K >
inline
bool
left_turn(const Point_2<K> &p, const Point_2<K> &q, const Point_2<K> &r)
{
  return CGALi::left_turn(p, q, r, K());
}

template < class K >
inline
bool
less_x(const Point_2<K> &p, const Point_2<K> &q)
{
  return CGALi::less_x(p, q, K());
}

template < class K >
inline
bool
less_y(const Point_2<K> &p, const Point_2<K> &q)
{
  return CGALi::less_y(p, q, K());
}

template < class K >
inline
bool
lexicographically_xy_larger(const Point_2<K> &p, const Point_2<K> &q)
{
  return CGALi::lexicographically_xy_larger(p, q, K());
}

template < class K >
inline
bool
lexicographically_xy_larger_or_equal(const Point_2<K> &p, const Point_2<K> &q)
{
  return CGALi::lexicographically_xy_larger_or_equal(p, q, K());
}

template < class K >
inline
bool
lexicographically_xy_smaller(const Point_2<K> &p, const Point_2<K> &q)
{
  return CGALi::lexicographically_xy_smaller(p, q, K());
}

template < class K >
inline
bool
lexicographically_xy_smaller_or_equal(const Point_2<K> &p,
                                      const Point_2<K> &q)
{
  return CGALi::lexicographically_xy_smaller_or_equal(p, q, K());
}

template < class K >
inline
bool
lexicographically_yx_smaller(const Point_2<K> &p, const Point_2<K> &q)
{
  return CGALi::lexicographically_yx_smaller(p, q, K());
}

template < class K >
inline
bool
lexicographically_yx_smaller_or_equal(const Point_2<K> &p, const Point_2<K> &q)
{
  return CGALi::lexicographically_yx_smaller_or_equal(p, q, K());
}

// FIXME : Undocumented
template < class K >
inline
bool
lexicographically_yx_larger(const Point_2<K> &p, const Point_2<K> &q)
{
  return CGALi::lexicographically_yx_larger(p, q, K());
}

// FIXME : Undocumented
template < class K >
inline
bool
lexicographically_yx_larger_or_equal(const Point_2<K> &p, const Point_2<K> &q)
{
  return CGALi::lexicographically_yx_larger_or_equal(p, q, K());
}

template < class K >
inline
typename K::Point_2
midpoint(const Point_2<K> &p, const Point_2<K> &q)
{
  return CGALi::midpoint(p, q, K());
}

template < class K >
inline
typename K::Point_2
max_vertex(const Iso_rectangle_2<K> &ir)
{
  return CGALi::max_vertex(ir, K());
}

template < class K >
inline
typename K::Point_2
min_vertex(const Iso_rectangle_2<K> &ir)
{
  return CGALi::min_vertex(ir, K());
}

// FIXME TODO : What do we do with the operators ?
// They have no counter part with the kernel argument...
template < class K >
inline
bool
operator<(const Direction_2<K>& d1, const Direction_2<K>& d2)
{ return compare_angle_with_x_axis(d1, d2) == SMALLER; }

template < class K >
inline
bool
operator>(const Direction_2<K>& d1, const Direction_2<K>& d2)
{ return compare_angle_with_x_axis(d1, d2) == LARGER; }

template < class K >
inline
bool
operator>=(const Direction_2<K>& d1, const Direction_2<K>& d2)
{ return compare_angle_with_x_axis(d1, d2) != SMALLER; }

template < class K >
inline
bool
operator<=(const Direction_2<K>& d1, const Direction_2<K>& d2)
{ return compare_angle_with_x_axis(d1, d2) != LARGER; }

template < class K >
inline
bool
operator<(const Point_2<K>& p, const Point_2<K>& q)
{ return K().less_xy_2_object()(p, q); }

template < class K >
inline
bool
operator>(const Point_2<K>& p, const Point_2<K>& q)
{ return K().less_xy_2_object()(q, p); }

template < class K >
inline
bool
operator<=(const Point_2<K>& p, const Point_2<K>& q)
{ return ! K().less_xy_2_object()(q, p); }

template < class K >
inline
bool
operator>=(const Point_2<K>& p, const Point_2<K>& q)
{ return ! K().less_xy_2_object()(p, q); }

template < class K >
inline
typename K::Vector_2
operator*(const typename CGAL_WRAP(K)::FT &c, const Vector_2<K> &w)
{
  return K().construct_scaled_vector_2_object()(w, c);
}

template < class K >
inline
typename K::Vector_2
operator*(const Vector_2<K> &w, const typename CGAL_WRAP(K)::FT &c)
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
Orientation
orientation(const Point_2<K> &p, const Point_2<K> &q, const Point_2<K> &r)
{
  return CGALi::orientation(p, q, r, K());
}

// parallel() functions are in global_functions.h

template <typename K>
inline
bool
right_turn(const Point_2<K> &p, const Point_2<K> &q, const Point_2<K> &r)
{
  return CGALi::right_turn(p, q, r, K());
}

template <class K>
inline
Bounded_side
side_of_bounded_circle(const Point_2<K> &p,
                       const Point_2<K> &q,
                       const Point_2<K> &r,
                       const Point_2<K> &t)
{
  return CGALi::side_of_bounded_circle(p, q, r, t, K());
}

template <class K>
inline
Bounded_side
side_of_bounded_circle(const Point_2<K> &p,
                       const Point_2<K> &q,
                       const Point_2<K> &r)
{
  return CGALi::side_of_bounded_circle(p, q, r, K());
}

template <class K>
inline
Oriented_side
side_of_oriented_circle(const Point_2<K> &p,
                        const Point_2<K> &q,
                        const Point_2<K> &r,
                        const Point_2<K> &t)
{
  return CGALi::side_of_oriented_circle(p, q, r, t, K());
}

template < class K >
inline
typename K::FT
squared_radius(const Point_2<K> &p, const Point_2<K> &q)
{
  return CGALi::squared_radius(p, q, K());
}

template < class K >
CGAL_KERNEL_INLINE
typename K::FT
squared_radius(const Point_2<K>& p, const Point_2<K>& q, const Point_2<K>& r)
{
  return CGALi::squared_radius(p, q, r, K());
}

template < class K >
inline
bool
x_equal(const Point_2<K> &p, const Point_2<K> &q)
{
  return CGALi::x_equal(p, q, K());
}

template < class K >
inline
bool
y_equal(const Point_2<K> &p, const Point_2<K> &q)
{
  return CGALi::y_equal(p, q, K());
}

CGAL_END_NAMESPACE

#endif  // CGAL_KERNEL_GLOBAL_FUNCTIONS_2_H
