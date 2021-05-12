// Copyright (c) 2003-2010
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),,
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_KERNEL_GLOBAL_FUNCTIONS_INTERNAL_2_H
#define CGAL_KERNEL_GLOBAL_FUNCTIONS_INTERNAL_2_H

// Generic functions calling the kernel functor, taking the kernel as
// parameter.

// These functions are not documented for now, but could be as some point.

#include <CGAL/basic.h>
#include <CGAL/Dimension.h>
#include <boost/utility/enable_if.hpp>
#include "boost/mpl/equal_to.hpp"
#include <boost/mpl/int.hpp>
#include <boost/mpl/integral_c.hpp>

namespace CGAL {

namespace internal {

template < class K >
inline
typename K::Angle
angle(const typename K::Vector_2 &u,
      const typename K::Vector_2 &v, const K& k)
{
  return k.angle_2_object()(u, v);
}

template < class K >
inline
typename K::Angle
angle(const typename K::Point_2 &p,
      const typename K::Point_2 &q,
      const typename K::Point_2 &r, const K& k)
{
  return k.angle_2_object()(p, q, r);
}

template < class K >
inline
typename K::Angle
angle(const typename K::Point_2 &p,
      const typename K::Point_2 &q,
      const typename K::Point_2 &r,
      const typename K::Point_2 &s, const K& k)
{
  return k.angle_2_object()(p, q, r, s);
}

template < class K >
inline
typename K::Boolean
are_ordered_along_line(const typename K::Point_2 &p,
                       const typename K::Point_2 &q,
                       const typename K::Point_2 &r, const K& k)
{
  return k.are_ordered_along_line_2_object()(p, q, r);
}

template < class K >
inline
typename K::Boolean
are_strictly_ordered_along_line(const typename K::Point_2 &p,
                                const typename K::Point_2 &q,
                                const typename K::Point_2 &r,
                                const K& k)
{
  return k.are_strictly_ordered_along_line_2_object()(p, q, r);
}

template < class K >
inline
typename K::FT
area(const typename K::Point_2 &p,
     const typename K::Point_2 &q,
     const typename K::Point_2 &r,
     const K& k)
{
  return k.compute_area_2_object()(p, q, r);
}

template < class K >
inline
typename K::Point_2
barycenter(const typename K::Point_2 &p1, const typename K::FT& w1,
           const typename K::Point_2 &p2, const K& k)
{
  return k.construct_barycenter_2_object()(p1, w1, p2);
}

template < class K >
inline
typename K::Point_2
barycenter(const typename K::Point_2 &p1, const typename K::FT& w1,
           const typename K::Point_2 &p2, const typename K::FT& w2, const K& k)
{
  return k.construct_barycenter_2_object()(p1, w1, p2, w2);
}

template < class K >
inline
typename K::Point_2
barycenter(const typename K::Point_2 &p1, const typename K::FT& w1,
           const typename K::Point_2 &p2, const typename K::FT& w2,
           const typename K::Point_2 &p3, const K& k)
{
  return k.construct_barycenter_2_object()(p1, w1, p2, w2, p3);
}

template < class K >
inline
typename K::Point_2
barycenter(const typename K::Point_2 &p1, const typename K::FT& w1,
           const typename K::Point_2 &p2, const typename K::FT& w2,
           const typename K::Point_2 &p3, const typename K::FT& w3, const K& k)
{
  return k.construct_barycenter_2_object()(p1, w1, p2, w2, p3, w3);
}

template < class K >
inline
typename K::Point_2
barycenter(const typename K::Point_2 &p1, const typename K::FT& w1,
           const typename K::Point_2 &p2, const typename K::FT& w2,
           const typename K::Point_2 &p3, const typename K::FT& w3,
           const typename K::Point_2 &p4, const K& k)
{
  return k.construct_barycenter_2_object()(p1, w1, p2, w2, p3, w3, p4);
}

template < class K >
inline
typename K::Point_2
barycenter(const typename K::Point_2 &p1, const typename K::FT& w1,
           const typename K::Point_2 &p2, const typename K::FT& w2,
           const typename K::Point_2 &p3, const typename K::FT& w3,
           const typename K::Point_2 &p4, const typename K::FT& w4, const K& k)
{
  return k.construct_barycenter_2_object()(p1, w1, p2, w2, p3, w3, p4, w4);
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

template < class K >
inline
typename K::Point_2
centroid(const typename K::Point_2 &p,
         const typename K::Point_2 &q,
         const typename K::Point_2 &r, const K& k)
{
  return k.construct_centroid_2_object()(p, q, r);
}

template < class K >
inline
typename K::Point_2
centroid(const typename K::Point_2 &p,
         const typename K::Point_2 &q,
         const typename K::Point_2 &r,
         const typename K::Point_2 &s, const K& k)
{
  return k.construct_centroid_2_object()(p, q, r, s);
}

template < class K >
inline
typename K::Point_2
centroid(const typename K::Triangle_2 &t, const K& k)
{
  return k.construct_centroid_2_object()(t);
}

template < class K >
inline
typename K::Point_2
circumcenter(const typename K::Point_2 &p,
             const typename K::Point_2 &q, const K& k)
{
  return k.construct_circumcenter_2_object()(p, q);
}

template < class K >
inline
typename K::Point_2
circumcenter(const typename K::Point_2 &p,
             const typename K::Point_2 &q,
             const typename K::Point_2 &r, const K& k)
{
  return k.construct_circumcenter_2_object()(p, q, r);
}

template < class K >
inline
typename K::Point_2
circumcenter(const typename K::Triangle_2 &t, const K& k)
{
  return k.construct_circumcenter_2_object()(t);
}

template < class K >
inline
typename K::Boolean
collinear(const typename K::Point_2 &p,
          const typename K::Point_2 &q,
          const typename K::Point_2 &r, const K& k)
{
  return k.collinear_2_object()(p, q, r);
}

template < class K >
inline
typename K::Boolean
collinear_are_ordered_along_line(const typename K::Point_2 &p,
                                 const typename K::Point_2 &q,
                                 const typename K::Point_2 &r,
                                 const K& k)
{
  return k.collinear_are_ordered_along_line_2_object()(p, q, r);
}

template < class K >
inline
typename K::Boolean
collinear_are_strictly_ordered_along_line(
             const typename K::Point_2 &p,
             const typename K::Point_2 &q,
             const typename K::Point_2 &r, const K& k)
{
  return k.collinear_are_strictly_ordered_along_line_2_object()(p, q, r);
}

template < typename K >
inline
typename K::Comparison_result
compare_angle_with_x_axis(const typename K::Direction_2& d1,
                          const typename K::Direction_2& d2,
                          const K& k)
{
  return k.compare_angle_with_x_axis_2_object()(d1, d2);
}

template <class K, class T1, class T2, class T3>
inline
typename boost::enable_if<
  boost::mpl::equal_to<boost::mpl::integral_c<int,
                                              Ambient_dimension<T1>::type::value>,
                       boost::mpl::integral_c<int, 2> >,
  typename K::Comparison_result>
::type
compare_distance(const T1 &o1,
                 const T2 &o2,
                 const T3 &o3, const K& k)
{
  return k.compare_distance_2_object()(o1, o2, o3);
}

template <class K, class T1, class T2, class T3, class T4>
inline
typename boost::enable_if<
  boost::mpl::equal_to<boost::mpl::integral_c<int,
                                              Ambient_dimension<T1>::type::value>,
                       boost::mpl::integral_c<int, 2> >,
  typename K::Comparison_result>
::type
compare_distance(const T1 &o1,
                 const T2 &o2,
                 const T3 &o3,
                 const T4 &o4, const K& k)
{
  return k.compare_distance_2_object()(o1, o2, o3, o4);
}

template <class K >
inline
typename K::Comparison_result
compare_distance_to_point(const typename K::Point_2 &p,
                          const typename K::Point_2 &q,
                          const typename K::Point_2 &r, const K& k)
{
  return k.compare_distance_2_object()(p, q, r);
}

template <class K >
inline
typename K::Comparison_result
compare_power_distance(const typename K::Point_2 &r,
                       const typename K::Weighted_point_2 &p,
                       const typename K::Weighted_point_2 &q, const K& k)
{
  return k.compare_power_distance_2_object()(r, p, q);
}

template <class K >
inline
typename K::Comparison_result
compare_squared_distance(const typename K::Point_2 &p,
                         const typename K::Point_2 &q,
                         const typename K::FT &d2, const K& k)
{
  return k.compare_squared_distance_2_object()(p, q, d2);
}

template <class K>
inline
typename K::Comparison_result
compare_signed_distance_to_line(const typename K::Point_2& p,
                                const typename K::Point_2& q,
                                const typename K::Point_2& r,
                                const typename K::Point_2& s,
                                const K& k)
{
  return k.compare_signed_distance_to_line_2_object()(p, q, r, s);
}

template <class K>
inline
typename K::Comparison_result
compare_signed_distance_to_line(const typename K::Line_2& l,
                                const typename K::Point_2& p,
                                const typename K::Point_2& q,
                                const K& k)
{
  return k.compare_signed_distance_to_line_2_object()(l, p, q);
}

template < class K >
inline
typename K::Comparison_result
compare_slope(const typename K::Line_2 &l1,
              const typename K::Line_2 &l2, const K& k)
{
  return k.compare_slope_2_object()(l1, l2);
}

template < class K >
inline
typename K::Comparison_result
compare_slope(const typename K::Segment_2 &s1,
              const typename K::Segment_2 &s2, const K& k)
{
  return k.compare_slope_2_object()(s1, s2);
}

template < class K >
inline
typename K::Comparison_result
compare_x(const typename K::Point_2 &p,
          const typename K::Point_2 &q, const K& k)
{
  return k.compare_x_2_object()(p, q);
}

template < class K >
inline
typename K::Comparison_result
compare_x(const typename K::Point_2 &p,
          const typename K::Line_2 &l1,
          const typename K::Line_2 &l2, const K& k)
{
  return k.compare_x_2_object()(p, l1, l2);
}

template < class K >
inline
typename K::Comparison_result
compare_x(const typename K::Line_2 &l,
          const typename K::Line_2 &h1,
          const typename K::Line_2 &h2, const K& k)
{
  return k.compare_x_2_object()(l, h1, h2);
}

template < class K >
inline
typename K::Comparison_result
compare_x(const typename K::Line_2 &l1,
          const typename K::Line_2 &h1,
          const typename K::Line_2 &l2,
          const typename K::Line_2 &h2, const K& k)
{
  return k.compare_x_2_object()(l1, h1, l2, h2);
}

template < class K >
inline
typename K::Comparison_result
compare_x_at_y(const typename K::Point_2& p,
               const typename K::Line_2& h, const K& k)
{
  return k.compare_x_at_y_2_object()(p, h);
}

/* Undocumented
template < class K >
inline
typename K::Comparison_result
compare_y_at_x(const typename K::Point_2 &p,
               const typename K::Segment_2 &s, const K& k)
{
  return k.compare_y_at_x_2_object()(p, s);
}
*/

template < class K >
inline
typename K::Comparison_result
compare_x_at_y(const typename K::Point_2 &p,
               const typename K::Line_2 &h1,
               const typename K::Line_2 &h2, const K& k)
{
  return k.compare_x_at_y_2_object()(p, h1, h2);
}

template < class K >
inline
typename K::Comparison_result
compare_x_at_y(const typename K::Line_2 &l1,
               const typename K::Line_2 &l2,
               const typename K::Line_2 &h, const K& k)
{
  return k.compare_x_at_y_2_object()(l1, l2, h);
}

template < class K >
inline
typename K::Comparison_result
compare_x_at_y(const typename K::Line_2 &l1,
               const typename K::Line_2 &l2,
               const typename K::Line_2 &h1,
               const typename K::Line_2 &h2, const K& k)
{
  return k.compare_x_at_y_2_object()(l1, l2, h1, h2);
}

template < class K >
inline
typename K::Comparison_result
compare_xy(const typename K::Point_2 &p,
           const typename K::Point_2 &q, const K& k)
{
  return k.compare_xy_2_object()(p, q);
}

template < class K >
inline
typename K::Comparison_result
compare_yx(const typename K::Point_2 &p,
           const typename K::Point_2 &q, const K& k)
{
  return k.compare_yx_2_object()(p, q);
}

template < class K >
inline
typename K::Comparison_result
compare_y(const typename K::Point_2 &p,
          const typename K::Point_2 &q, const K& k)
{
  return k.compare_y_2_object()(p, q);
}

template < class K >
inline
typename K::Comparison_result
compare_y(const typename K::Point_2 &p,
          const typename K::Line_2 &l1,
          const typename K::Line_2 &l2, const K& k)
{
  return k.compare_y_2_object()(p, l1, l2);
}

template < class K >
inline
typename K::Comparison_result
compare_y(const typename K::Line_2 &l1,
          const typename K::Line_2 &l2,
          const typename K::Line_2 &h1,
          const typename K::Line_2 &h2, const K& k)
{
  return k.compare_y_2_object()(l1, l2, h1, h2);
}

template < class K >
inline
typename K::Comparison_result
compare_y(const typename K::Line_2 &l,
          const typename K::Line_2 &h1,
          const typename K::Line_2 &h2, const K& k)
{
  return k.compare_y_2_object()(l, h1, h2);
}

template < class K >
inline
typename K::Comparison_result
compare_y_at_x(const typename K::Point_2 &p,
               const typename K::Segment_2 &s, const K& k)
{
  return k.compare_y_at_x_2_object()(p, s);
}

template < class K >
inline
typename K::Comparison_result
compare_y_at_x(const typename K::Point_2 &p,
               const typename K::Segment_2 &s1,
               const typename K::Segment_2 &s2, const K& k)
{
  return k.compare_y_at_x_2_object()(p, s1, s2);
}

template < class K >
inline
typename K::Comparison_result
compare_y_at_x(const typename K::Point_2 &p,
               const typename K::Line_2 &l, const K& k)
{
  return k.compare_y_at_x_2_object()(p, l);
}

template < class K >
inline
typename K::Comparison_result
compare_y_at_x(const typename K::Point_2 &p,
               const typename K::Line_2 &h1,
               const typename K::Line_2 &h2, const K& k)
{
  return k.compare_y_at_x_2_object()(p, h1, h2);
}

template < class K >
inline
typename K::Comparison_result
compare_y_at_x(const typename K::Line_2 &l1,
               const typename K::Line_2 &l2,
               const typename K::Line_2 &h, const K& k)
{
  return k.compare_y_at_x_2_object()(l1, l2, h);
}

template < class K >
inline
typename K::Comparison_result
compare_y_at_x(const typename K::Line_2 &l1,
               const typename K::Line_2 &l2,
               const typename K::Line_2 &h1,
               const typename K::Line_2 &h2, const K& k)
{
  return k.compare_y_at_x_2_object()(l1, l2, h1, h2);
}

template < class K >
inline
typename K::FT
determinant(const typename K::Vector_2 &v0,
            const typename K::Vector_2 &v1, const K &k)
{
  return k.compute_determinant_2_object()(v0, v1);
}

template <class K>
inline
typename K::Boolean
has_larger_distance_to_point(const typename K::Point_2 &p,
                             const typename K::Point_2 &q,
                             const typename K::Point_2 &r,
                             const K& k)
{
  return k.less_distance_to_point_2_object()(p, r, q);
}

template <class K>
inline
typename K::Boolean
has_smaller_distance_to_point(const typename K::Point_2 &p,
                              const typename K::Point_2 &q,
                              const typename K::Point_2 &r,
                              const K& k)
{
  return k.less_distance_to_point_2_object()(p, q, r);
}

template <class K>
inline
typename K::Boolean
has_smaller_signed_distance_to_line(const typename K::Line_2& l,
                                    const typename K::Point_2& p,
                                    const typename K::Point_2& q,
                                    const K& k)
{
  return k.less_signed_distance_to_line_2_object()(l, p, q);
}

template <class K>
inline
typename K::Boolean
has_larger_signed_distance_to_line(const typename K::Line_2& l,
                                   const typename K::Point_2& p,
                                   const typename K::Point_2& q,
                                   const K& k)
{
  return k.less_signed_distance_to_line_2_object()(l, q, p);
}

template <class K>
inline
typename K::Boolean
has_larger_signed_distance_to_line(const typename K::Point_2& p,
                                   const typename K::Point_2& q,
                                   const typename K::Point_2& r,
                                   const typename K::Point_2& s,
                                   const K& k)
{
  return k.less_signed_distance_to_line_2_object()(p, q, s, r);
}

template <class K>
inline
typename K::Boolean
has_smaller_signed_distance_to_line(const typename K::Point_2& p,
                                    const typename K::Point_2& q,
                                    const typename K::Point_2& r,
                                    const typename K::Point_2& s,
                                    const K& k)
{
  return k.less_signed_distance_to_line_2_object()(p, q, r, s);
}

template < class K >
inline
typename K::Boolean
left_turn(const typename K::Point_2 &p,
          const typename K::Point_2 &q,
          const typename K::Point_2 &r, const K& k)
{
  return k.left_turn_2_object()(p, q, r);
}

template < class K >
inline
typename K::Boolean
less_x(const typename K::Point_2 &p,
       const typename K::Point_2 &q, const K& k)
{
  return k.less_x_2_object()(p, q);
}

template < class K >
inline
typename K::Boolean
less_y(const typename K::Point_2 &p,
       const typename K::Point_2 &q, const K& k)
{
  return k.less_y_2_object()(p, q);
}

template < class K >
inline
typename K::Boolean
lexicographically_xy_larger(const typename K::Point_2 &p,
                            const typename K::Point_2 &q,
                            const K& k)
{
  return k.compare_xy_2_object()(p, q) == LARGER;
}

template < class K >
inline
typename K::Boolean
lexicographically_xy_larger_or_equal(const typename K::Point_2 &p,
                                     const typename K::Point_2 &q,
                                     const K& k)
{
  return k.compare_xy_2_object()(p, q) != SMALLER;
}

template < class K >
inline
typename K::Boolean
lexicographically_xy_smaller(const typename K::Point_2 &p,
                             const typename K::Point_2 &q,
                             const K& k)
{
  return k.less_xy_2_object()(p, q);
}

template < class K >
inline
typename K::Boolean
lexicographically_xy_smaller_or_equal(const typename K::Point_2 &p,
                                      const typename K::Point_2 &q,
                                      const K& k)
{
  return k.compare_xy_2_object()(p, q) != LARGER;
}

template < class K >
inline
typename K::Boolean
lexicographically_yx_smaller(const typename K::Point_2 &p,
                             const typename K::Point_2 &q,
                             const K& k)
{
  return k.less_yx_2_object()(p, q);
}

template < class K >
inline
typename K::Boolean
lexicographically_yx_smaller_or_equal(const typename K::Point_2 &p,
                                      const typename K::Point_2 &q,
                                      const K& k)
{
  return !k.less_yx_2_object()(q, p);
}

// FIXME : Undocumented
template < class K >
inline
typename K::Boolean
lexicographically_yx_larger(const typename K::Point_2 &p,
                            const typename K::Point_2 &q,
                            const K& k)
{
  return k.less_yx_2_object()(q, p);
}

// FIXME : Undocumented
template < class K >
inline
typename K::Boolean
lexicographically_yx_larger_or_equal(const typename K::Point_2 &p,
                                     const typename K::Point_2 &q,
                                     const K& k)
{
  return !k.less_yx_2_object()(p, q);
}

template < class K >
inline
typename K::FT
l_infinity_distance(const typename K::Point_2 &p,
                    const typename K::Point_2 &q,
                    const K& k)
{
  return k.compute_L_infinity_distance_2_object()(p, q);
}


template < class K >
inline
typename K::Point_2
midpoint(const typename K::Point_2 &p,
         const typename K::Point_2 &q, const K &k)
{
  return k.construct_midpoint_2_object()(p, q);
}

template < class K >
inline
typename K::Point_2
max_vertex(const typename K::Iso_rectangle_2 &ir, const K &k)
{
  return k.construct_max_vertex_2_object()(ir);
}

template < class K >
inline
typename K::Point_2
min_vertex(const typename K::Iso_rectangle_2 &ir, const K &k)
{
  return k.construct_min_vertex_2_object()(ir);
}

template <typename K>
inline
typename K::Orientation
orientation(const typename K::Point_2 &p,
            const typename K::Point_2 &q,
            const typename K::Point_2 &r, const K &k)
{
  return k.orientation_2_object()(p, q, r);
}

template <typename K>
inline
typename K::Orientation
orientation(const typename K::Vector_2 &u,
            const typename K::Vector_2 &v, const K &k)
{
  return k.orientation_2_object()(u, v);
}

template <typename K>
inline
typename K::Boolean
parallel(const typename K::Line_2 &l1,
         const typename K::Line_2 &l2, const K &k)
{
  return k.are_parallel_2_object()(l1, l2);
}

template <typename K>
inline
typename K::Boolean
parallel(const typename K::Ray_2 &r1,
         const typename K::Ray_2 &r2, const K &k)
{
  return k.are_parallel_2_object()(r1, r2);
}

template <typename K>
inline
typename K::Boolean
parallel(const typename K::Segment_2 &s1,
         const typename K::Segment_2 &s2, const K &k)
{
  return k.are_parallel_2_object()(s1, s2);
}

template <class K >
inline
typename K::FT
power_product(const typename K::Weighted_point_2 &p,
              const typename K::Weighted_point_2 &q, const K &k)
{
  return k.compute_power_product_2_object()(p, q);
}

template <class K >
inline
typename K::Bounded_side
power_side_of_bounded_power_circle(const typename K::Weighted_point_2 &p,
                                   const typename K::Weighted_point_2 &q, const K &k)
{
  return k.power_side_of_bounded_power_circle_2_object()(p, q);
}

template <class K >
inline
typename K::Bounded_side
power_side_of_bounded_power_circle(const typename K::Weighted_point_2 &p,
                                   const typename K::Weighted_point_2 &q,
                                   const typename K::Weighted_point_2 &r, const K &k)
{
  return k.power_side_of_bounded_power_circle_2_object()(p, q, r);
}

template <class K >
inline
typename K::Bounded_side
power_side_of_bounded_power_circle(const typename K::Weighted_point_2 &p,
                                   const typename K::Weighted_point_2 &q,
                                   const typename K::Weighted_point_2 &r,
                                   const typename K::Weighted_point_2 &s, const K &k)
{
  return k.power_side_of_bounded_power_circle_2_object()(p, q, r, s);
}

template <class K>
inline
typename K::Oriented_side
power_side_of_oriented_power_circle(const typename K::Weighted_point_2 &p,
                                    const typename K::Weighted_point_2 &q, const K &k)
{
  return k.power_side_of_oriented_power_circle_2_object()(p, q);
}

template <class K>
inline
typename K::Oriented_side
power_side_of_oriented_power_circle(const typename K::Weighted_point_2 &p,
                                    const typename K::Weighted_point_2 &q,
                                    const typename K::Weighted_point_2 &r, const K &k)
{
  return k.power_side_of_oriented_power_circle_2_object()(p, q, r);
}

template <class K>
inline
typename K::Oriented_side
power_side_of_oriented_power_circle(const typename K::Weighted_point_2 &p,
                                    const typename K::Weighted_point_2 &q,
                                    const typename K::Weighted_point_2 &r,
                                    const typename K::Weighted_point_2 &t, const K &k)
{
  return k.power_side_of_oriented_power_circle_2_object()(p, q, r, t);
}

template <typename K>
inline
typename K::Line_2
radical_axis(const typename K::Weighted_point_2 &p,
             const typename K::Weighted_point_2 &q, const K &k)
{
  return k.construct_radical_axis_2_object()(p, q);
}

template <typename K>
inline
typename K::Boolean
right_turn(const typename K::Point_2 &p,
           const typename K::Point_2 &q,
           const typename K::Point_2 &r, const K &k)
{
  return internal::orientation(p, q, r, k) == RIGHT_TURN;
}

template <class K>
inline
typename K::Bounded_side
side_of_bounded_circle(const typename K::Point_2 &p,
                       const typename K::Point_2 &q,
                       const typename K::Point_2 &r,
                       const typename K::Point_2 &t, const K &k)
{
  return k.side_of_bounded_circle_2_object()(p, q, r, t);
}

template <class K>
inline
typename K::Bounded_side
side_of_bounded_circle(const typename K::Point_2 &p,
                       const typename K::Point_2 &q,
                       const typename K::Point_2 &r, const K &k)
{
  return k.side_of_bounded_circle_2_object()(p, q, r);
}

template <class K>
inline
typename K::Oriented_side
side_of_oriented_circle(const typename K::Point_2 &p,
                        const typename K::Point_2 &q,
                        const typename K::Point_2 &r,
                        const typename K::Point_2 &t, const K &k)
{
  return k.side_of_oriented_circle_2_object()(p, q, r, t);
}

template < class K >
inline
typename K::FT
squared_radius(const typename K::Point_2 &p, const K &k)
{
  return k.compute_squared_radius_2_object()(p);
}

template < class K >
inline
typename K::FT
squared_radius(const typename K::Point_2 &p,
               const typename K::Point_2 &q, const K &k)
{
  return k.compute_squared_radius_2_object()(p, q);
}

template < class K >
inline
typename K::FT
squared_radius(const typename K::Point_2 &p,
               const typename K::Point_2 &q,
               const typename K::Point_2 &r, const K &k)
{
  return k.compute_squared_radius_2_object()(p, q, r);
}

template < class K >
inline
typename K::FT
squared_radius_smallest_orthogonal_circle(const typename K::Weighted_point_2 &p,
                                          const K &k)
{
  return k.compute_squared_radius_smallest_orthogonal_circle_2_object()(p);
}

template < class K >
inline
typename K::FT
squared_radius_smallest_orthogonal_circle(const typename K::Weighted_point_2 &p,
                                          const typename K::Weighted_point_2 &q,
                                          const K &k)
{
  return k.compute_squared_radius_smallest_orthogonal_circle_2_object()(p, q);
}

template < class K >
inline
typename K::FT
squared_radius_smallest_orthogonal_circle(const typename K::Weighted_point_2 &p,
                                          const typename K::Weighted_point_2 &q,
                                          const typename K::Weighted_point_2 &r,
                                          const K &k)
{
  return k.compute_squared_radius_smallest_orthogonal_circle_2_object()(p, q, r);
}

template < class K >
inline
typename K::Point_2
weighted_circumcenter(const typename K::Weighted_point_2 &p,
                      const typename K::Weighted_point_2 &q,
                      const typename K::Weighted_point_2 &r, const K &k)
{
  return k.construct_weighted_circumcenter_2_object()(p, q, r);
}

template < class K >
inline
typename K::Boolean
x_equal(const typename K::Point_2 &p,
        const typename K::Point_2 &q, const K &k)
{
  return k.equal_x_2_object()(p, q);
}

template < class K >
inline
typename K::Boolean
y_equal(const typename K::Point_2 &p,
        const typename K::Point_2 &q, const K &k)
{
  return k.equal_y_2_object()(p, q);
}

} // namespace internal

} //namespace CGAL

#endif  // CGAL_KERNEL_GLOBAL_FUNCTIONS_INTERNAL_2_H
