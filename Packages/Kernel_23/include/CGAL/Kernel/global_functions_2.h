// Copyright (c) 2003  Utrecht University (The Netherlands),
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

// Generic functions calling the kernel functor.

#include <CGAL/user_classes.h>

CGAL_BEGIN_NAMESPACE

template < class K >
inline
Angle
angle(const Point_2<K> &p,
      const Point_2<K> &q,
      const Point_2<K> &r)
{
  return K().angle_2_object()(p, q, r);
}

template < class K >
inline
bool
are_ordered_along_line(const Point_2<K> &p,
                       const Point_2<K> &q,
                       const Point_2<K> &r)
{
  return K().are_ordered_along_line_2_object()(p, q, r);
}

template < class K >
inline
bool
are_strictly_ordered_along_line(const Point_2<K> &p,
                                const Point_2<K> &q,
                                const Point_2<K> &r)
{
  return K().are_strictly_ordered_along_line_2_object()(p, q, r);
}

template <typename K>
inline
typename K::Line_2
bisector(const typename K::Point_2 &p,
         const typename K::Point_2 &q, const K &k)
{
  return k.construct_bisector_2_object()(p, q);
}

template <typename K>
inline
typename K::Line_2
bisector(const typename K::Line_2 &l1,
         const typename K::Line_2 &l2, const K &k)
{
  return k.construct_bisector_2_object()(l1, l2);
}

template <typename K>
inline
typename K::Line_2
bisector(const Point_2<K> &p, const Point_2<K> &q)
{
  return K().construct_bisector_2_object()(p, q);
}

template <typename K>
inline
typename K::Line_2
bisector(const Line_2<K> &l1, const Line_2<K> &l2)
{
  return K().construct_bisector_2_object()(l1, l2);
}

template < class K >
inline
typename K::Point_2
centroid(const Point_2<K> &p,
         const Point_2<K> &q,
         const Point_2<K> &r)
{
  return K().construct_centroid_2_object()(p, q, r);
}

template < class K >
inline
typename K::Point_2
centroid(const Point_2<K> &p,
         const Point_2<K> &q,
         const Point_2<K> &r,
         const Point_2<K> &s)
{
  return K().construct_centroid_2_object()(p, q, r, s);
}

template < class K >
inline
typename K::Point_2
circumcenter(const Point_2<K> &p,
             const Point_2<K> &q,
             const Point_2<K> &r)
{
  return K().construct_circumcenter_2_object()(p, q, r);
}

template < class K >
inline
bool
collinear(const Point_2<K> &p, const Point_2<K> &q, const Point_2<K> &r)
{
  return K().collinear_2_object()(p, q, r);
}

template < class K >
inline
bool
collinear_are_ordered_along_line(const Point_2<K> &p,
                                 const Point_2<K> &q,
                                 const Point_2<K> &r)
{
  return K().collinear_are_ordered_along_line_2_object()(p, q, r);
}

template < class K >
inline
bool
collinear_are_strictly_ordered_along_line(const Point_2<K> &p,
                                          const Point_2<K> &q,
                                          const Point_2<K> &r)
{
  return K().collinear_are_strictly_ordered_along_line_2_object()(p, q, r);
}

template < typename K >
inline
Comparison_result
compare_angle_with_x_axis(const typename CGAL_WRAP(K)::Direction_2& d1,
                          const typename CGAL_WRAP(K)::Direction_2& d2,
                          const K& k)
{
  return k.compare_angle_with_x_axis_2_object()(d1, d2);
}

template < typename K >
inline
Comparison_result
compare_angle_with_x_axis(const Direction_2<K>& d1,
                          const Direction_2<K>& d2)
{
  return K().compare_angle_with_x_axis_2_object()(d1, d2);
}

template <class K >
inline
Comparison_result
compare_distance_to_point(const Point_2<K>& p,
                          const Point_2<K>& q,
                          const Point_2<K>& r)
{
  return K().compare_distance_2_object()(p, q, r);
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
compare_x(const Point_2<K> &p, const Point_2<K> &q)
{
  return K().compare_x_2_object()(p, q);
}

template < class K >
inline
Comparison_result
compare_x(const Point_2<K>& p,
          const Line_2<K>& l1,
          const Line_2<K>& l2)
{
  return K().compare_x_2_object()(p, l1, l2);
}

template < class K >
inline
Comparison_result
compare_x(const Line_2<K> &l,
          const Line_2<K> &h1,
          const Line_2<K> &h2)
{
  return K().compare_x_2_object()(l, h1, h2);
}

template < class K >
inline
Comparison_result
compare_x(const Line_2<K> &l1,
          const Line_2<K> &h1,
          const Line_2<K> &l2,
          const Line_2<K> &h2)
{
  return K().compare_x_2_object()(l1, h1, l2, h2);
}

/* Conflicts...
template < class K >
inline
Comparison_result
compare_x_at_y(const typename CGAL_WRAP(K)::Point_2& p,
               const typename CGAL_WRAP(K)::Line_2& h, const K& k)
{
  return K().compare_x_at_y_2_object()(p, h);
}
*/

template < class K >
inline
Comparison_result
compare_x_at_y(const Point_2<K>& p, const Line_2<K>& h)
{
  return K().compare_x_at_y_2_object()(p, h);
}

template < class K >
inline
Comparison_result
compare_x_at_y(const Point_2<K> &p,
               const Line_2<K> &h1,
               const Line_2<K> &h2)
{
  return K().compare_x_at_y_2_object()(p, h1, h2);
}

template < class K >
inline
Comparison_result
compare_x_at_y(const Line_2<K> &l1,
               const Line_2<K> &l2,
               const Line_2<K> &h)
{
  return K().compare_x_at_y_2_object()(l1, l2, h);
}

template < class K >
inline
Comparison_result
compare_x_at_y(const Line_2<K> &l1,
               const Line_2<K> &l2,
               const Line_2<K> &h1,
               const Line_2<K> &h2)
{
  return K().compare_x_at_y_2_object()(l1, l2, h1, h2);
}

template < class K >
inline
Comparison_result
compare_xy(const Point_2<K> &p, const Point_2<K> &q)
{
  return K().compare_xy_2_object()(p, q);
}

template < class K >
inline
Comparison_result
compare_y(const Point_2<K> &p, const Point_2<K> &q)
{
  return K().compare_y_2_object()(p, q);
}

template < class K >
inline
Comparison_result
compare_y(const Point_2<K> &p,
          const Line_2<K> &l1,
          const Line_2<K> &l2)
{
  return K().compare_y_2_object()(p, l1, l2);
}

template < class K >
inline
Comparison_result
compare_y(const Line_2<K> &l1,
          const Line_2<K> &l2,
          const Line_2<K> &h1,
          const Line_2<K> &h2)
{
  return K().compare_y_2_object()(l1, l2, h1, h2);
}

template < class K >
inline
Comparison_result
compare_y(const Line_2<K> &l,
          const Line_2<K> &h1,
          const Line_2<K> &h2)
{
  return K().compare_y_2_object()(l, h1, h2);
}

/*  This one is conflicting with (Point, Segment, Segment)...
    So disabled for now.
template < class K >
inline
Comparison_result
compare_y_at_x(const typename CGAL_WRAP(K)::Point_2 &p,
               const typename CGAL_WRAP(K)::Segment_2 &s, const K& k)
{
  return k.compare_y_at_x_2_object()(p, s);
}
*/

template < class K >
inline
Comparison_result
compare_y_at_x(const Point_2<K> &p, const Segment_2<K> &s)
{
  return K().compare_y_at_x_2_object()(p, s);
}

/* Conflicts.
template < class K >
inline
Comparison_result
compare_y_at_x(const typename CGAL_WRAP(K)::Point_2 &p,
               const typename CGAL_WRAP(K)::Segment_2 &s1,
               const typename CGAL_WRAP(K)::Segment_2 &s2, const K& k)
{
  return k.compare_y_at_x_2_object()(p, s1, s2);
}
*/

template < class K >
inline
Comparison_result
compare_y_at_x(const Point_2<K> &p,
               const Segment_2<K> &s1,
               const Segment_2<K> &s2)
{
  return K().compare_y_at_x_2_object()(p, s1, s2);
}

template < class K >
inline
Comparison_result
compare_y_at_x(const Point_2<K> &p, const Line_2<K> &h)
{
  return K().compare_y_at_x_2_object()(p, h);
}  

template < class K >
inline
Comparison_result
compare_y_at_x(const Point_2<K> &p,
               const Line_2<K> &h1,
               const Line_2<K> &h2)
{
  return K().compare_y_at_x_2_object()(p, h1, h2);
}

template < class K >
inline
Comparison_result
compare_y_at_x(const Line_2<K> &l1,
               const Line_2<K> &l2,
               const Line_2<K> &h)
{
  return K().compare_y_at_x_2_object()(l1, l2, h);
}

template < class K >
inline
Comparison_result
compare_y_at_x(const Line_2<K> &l1,
               const Line_2<K> &l2,
               const Line_2<K> &h1,
               const Line_2<K> &h2)
{
  return K().compare_y_at_x_2_object()(l1, l2, h1, h2);
}

template <class K>
inline
bool
has_smaller_distance_to_point(const Point_2<K>& p,
                              const Point_2<K>& q,
                              const Point_2<K>& r)
{
  return K().less_distance_to_point_2_object()(p, q, r);
}

template < class K >
inline
bool
left_turn(const Point_2<K> &p, const Point_2<K> &q, const Point_2<K> &r)
{
  return K().left_turn_2_object()(p, q, r);
}

template < class K >
inline
bool
less_x(const Point_2<K> &p, const Point_2<K> &q)
{
  return K().less_x_2_object()(p, q);
}

template < class K >
inline
bool
less_y(const Point_2<K> &p, const Point_2<K> &q)
{
  return K().less_y_2_object()(p, q);
}

template < class K >
inline
bool
lexicographically_xy_larger(const Point_2<K> &p, const Point_2<K> &q)
{
  return K().compare_xy_2_object()(p, q) == LARGER;
}

template < class K >
inline
bool
lexicographically_xy_larger_or_equal(const Point_2<K> &p, const Point_2<K> &q)
{
  return K().compare_xy_2_object()(p, q) != SMALLER;
}

template < class K >
inline
bool
lexicographically_xy_smaller(const Point_2<K> &p, const Point_2<K> &q)
{
  return K().less_xy_2_object()(p, q);
}

template < class K >
inline
bool
lexicographically_xy_smaller_or_equal(const Point_2<K> &p,
                                      const Point_2<K> &q)
{
  return K().compare_xy_2_object()(p, q) != LARGER;
}

template < class K >
inline
bool
lexicographically_yx_smaller(const Point_2<K> &p, const Point_2<K> &q)
{
  return K().less_yx_2_object()(p, q);
}

// FIXME : Undocumented
template < class K >
inline
bool
lexicographically_yx_larger(const Point_2<K> &p, const Point_2<K> &q)
{
  return lexicographically_yx_smaller(q, p);
}

// FIXME : Undocumented
template < class K >
inline
bool
lexicographically_yx_larger_or_equal(const Point_2<K> &p, const Point_2<K> &q)
{
  return !lexicographically_yx_smaller(p, q);
}

template < class K >
inline
typename K::Point_2
midpoint(const typename CGAL_WRAP(K)::Point_2 &p,
         const typename CGAL_WRAP(K)::Point_2 &q, const K &k)
{
  return k.construct_midpoint_2_object()(p, q);
}

template < class K >
inline
typename K::Point_2
midpoint(const Point_2<K> &p, const Point_2<K> &q)
{
  return K().construct_midpoint_2_object()(p, q);
}

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

template <typename K>
inline
Orientation
orientation(const typename CGAL_WRAP(K)::Point_2 &p,
            const typename CGAL_WRAP(K)::Point_2 &q,
            const typename CGAL_WRAP(K)::Point_2 &r, const K &k)
{
  return k.orientation_2_object()(p, q, r);
}

template <typename K>
inline
Orientation
orientation(const Point_2<K> &p,
            const Point_2<K> &q,
            const Point_2<K> &r)
{
  return K().orientation_2_object()(p, q, r);
}

template <typename K>
inline
bool
parallel(const typename CGAL_WRAP(K)::Line_2 &l1,
         const typename CGAL_WRAP(K)::Line_2 &l2, const K &k)
{
  return k.are_parallel_2_object()(l1, l2);
}

template <typename K>
inline
bool
parallel(const typename CGAL_WRAP(K)::Ray_2 &r1,
         const typename CGAL_WRAP(K)::Ray_2 &r2, const K &k)
{
  return k.are_parallel_2_object()(r1, r2);
}

template <typename K>
inline
bool
parallel(const typename CGAL_WRAP(K)::Segment_2 &s1,
         const typename CGAL_WRAP(K)::Segment_2 &s2, const K &k)
{
  return k.are_parallel_2_object()(s1, s2);
}

template <typename K>
inline
bool
right_turn(const typename CGAL_WRAP(K)::Point_2 &p,
           const typename CGAL_WRAP(K)::Point_2 &q,
           const typename CGAL_WRAP(K)::Point_2 &r, const K &k)
{
  return orientation(p, q, r, k) == RIGHT_TURN;
}

template <typename K>
inline
bool
right_turn(const Point_2<K> &p,
           const Point_2<K> &q,
           const Point_2<K> &r)
{
  return right_turn(p, q, r, K());
}

template <class K>
inline
Bounded_side
side_of_bounded_circle(const typename CGAL_WRAP(K)::Point_2 &p,
                       const typename CGAL_WRAP(K)::Point_2 &q,
                       const typename CGAL_WRAP(K)::Point_2 &r,
                       const typename CGAL_WRAP(K)::Point_2 &t, const K &k)
{
  return k.side_of_bounded_circle_2_object()(p, q, r, t);
}

template <class K>
inline
Bounded_side
side_of_bounded_circle(const Point_2<K> &p,
                       const Point_2<K> &q,
                       const Point_2<K> &r,
                       const Point_2<K> &t)
{
  return K().side_of_bounded_circle_2_object()(p, q, r, t);
}

template <class K>
inline
Bounded_side
side_of_bounded_circle(const typename CGAL_WRAP(K)::Point_2 &p,
                       const typename CGAL_WRAP(K)::Point_2 &q,
                       const typename CGAL_WRAP(K)::Point_2 &r, const K &k)
{
  return k.side_of_bounded_circle_2_object()(p, q, r);
}

template <class K>
inline
Bounded_side
side_of_bounded_circle(const Point_2<K> &p,
                       const Point_2<K> &q,
                       const Point_2<K> &r)
{
  return K().side_of_bounded_circle_2_object()(p, q, r);
}

template <class K>
inline
Oriented_side
side_of_oriented_circle(const typename CGAL_WRAP(K)::Point_2 &p,
                        const typename CGAL_WRAP(K)::Point_2 &q,
                        const typename CGAL_WRAP(K)::Point_2 &r,
                        const typename CGAL_WRAP(K)::Point_2 &t, const K &k)
{
  return k.side_of_oriented_circle_2_object()(p, q, r, t);
}

template <class K>
inline
Oriented_side
side_of_oriented_circle(const Point_2<K> &p,
                        const Point_2<K> &q,
                        const Point_2<K> &r,
                        const Point_2<K> &t)
{
  return K().side_of_oriented_circle_2_object()(p, q, r, t);
}

template < class K >
inline
typename K::FT
squared_radius(const Point_2<K> &p, const Point_2<K> &q)
{
  return K().compute_squared_radius_2_object()(p, q);
}

template < class K >
CGAL_KERNEL_INLINE
typename K::FT
squared_radius( const Point_2<K>& p,
                const Point_2<K>& q,
                const Point_2<K>& r )
{
  return K().compute_squared_radius_2_object()(p, q, r);
}

template < class K >
inline
bool
x_equal(const Point_2<K> &p, const Point_2<K> &q)
{
  return K().equal_x_2_object()(p, q);
}

template < class K >
inline
bool
y_equal(const Point_2<K> &p, const Point_2<K> &q)
{
  return K().equal_y_2_object()(p, q);
}

CGAL_END_NAMESPACE

#endif  // CGAL_KERNEL_GLOBAL_FUNCTIONS_2_H
