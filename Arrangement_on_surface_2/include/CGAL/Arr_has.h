// Copyright (c) 2005,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Efi Fogel      <efif@post.tau.ac.il>
//            Eric Berberich <ericb@post.tau.ac.il>
//            Shepard Liu    <shepard0liu@gmail.com>

#ifndef CGAL_ARR_HAS_H
#define CGAL_ARR_HAS_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <type_traits>

#include <CGAL/Bbox_2.h>
#include <CGAL/Arr_enums.h>

namespace CGAL {

namespace detail {
  // A type implicitly convertible to anything
  struct any_type {
    template <typename T>
    operator T() const;
  };
}

  // `Compare_x_2`
// Fallback selected if `Compare_x_2` is not defined in the traits `T` below.
template <typename, typename = std::void_t<>>
struct has_compare_x_2 : std::false_type {};

// Partial specialization selected if `T::Compare_x_2` is defined.
template <typename T>
struct has_compare_x_2<T, std::void_t<typename T::Compare_x_2>> : std::true_type {};

// `Compare_xy_2`
// Fallback selected if `Compare_xy_2` is not defined in the traits `T` below.
template <typename, typename = std::void_t<>>
struct has_compare_xy_2 : std::false_type {};

// Partial specialization selected if `T::Compare_xy_2` is defined.
template <typename T>
struct has_compare_xy_2<T, std::void_t<typename T::Compare_xy_2>> : std::true_type {};

// `Construct_min_vertex_2`
// Fallback selected if `Construct_min_vertex_2` is not defined in the traits `T` below.
template <typename, typename = std::void_t<>>
struct has_construct_min_vertex_2 : std::false_type {};

// Partial specialization selected if `T::Construct_min_vertex_2` is defined.
template <typename T>
struct has_construct_min_vertex_2<T, std::void_t<typename T::Construct_min_vertex_2>> : std::true_type {};

// `Construct_max_vertex_2`
// Fallback selected if `Construct_max_vertex_2` is not defined in the traits `T` below.
template <typename, typename = std::void_t<>>
struct has_construct_max_vertex_2 : std::false_type {};

// Partial specialization selected if `T::Construct_max_vertex_2` is defined.
template <typename T>
struct has_construct_max_vertex_2<T, std::void_t<typename T::Construct_max_vertex_2>> : std::true_type {};

// `Is_vertical_2`
// Fallback selected if `Is_vertical_2` is not defined in the traits `T` below.
template <typename, typename = std::void_t<>>
struct has_is_vertical_2 : std::false_type {};

// Partial specialization selected if `T::Is_vertical_2` is defined.
template <typename T>
struct has_is_vertical_2<T, std::void_t<typename T::Is_vertical_2>> : std::true_type {};

// `Compare_y_at_x_2`
// Fallback selected if `Compare_y_at_x_2` is not defined in the traits `T` below.
template <typename, typename = std::void_t<>>
struct has_compare_y_at_x_2 : std::false_type {};

// Partial specialization selected if `T::Compare_y_at_x_2` is defined.
template <typename T>
struct has_compare_y_at_x_2<T, std::void_t<typename T::Compare_y_at_x_2>> : std::true_type {};

// `Equal_2`
// Fallback selected if `Equal_2` is not defined in the traits `T` below.
template <typename, typename = std::void_t<>>
struct has_equal_2 : std::false_type {};

// Partial specialization selected if `T::Equal_2` is defined.
template <typename T>
struct has_equal_2<T, std::void_t<typename T::Equal>> : std::true_type {};

// `Compare_y_at_x_left_2`
// Fallback selected if `Compare_y_at_x_left_2` is not defined in the traits `T` below.
template <typename, typename = std::void_t<>>
struct has_compare_y_at_x_left_2 : std::false_type {};

// Partial specialization selected if `T::Compare_y_at_x_left_2` is defined.
template <typename T>
struct has_compare_y_at_x_left_2<T, std::void_t<typename T::Compare_y_at_x_left_2>> : std::true_type {};

// `Compare_y_at_x_right_2`
// Fallback selected if `Compare_y_at_x_right_2` is not defined in the traits `T` below.
template <typename, typename = std::void_t<>>
struct has_compare_y_at_x_right_2 : std::false_type {};

// Partial specialization selected if `T::Compare_y_at_x_right_2` is defined.
template <typename T>
struct has_compare_y_at_x_right_2<T, std::void_t<typename T::Compare_y_at_x_right_2>> : std::true_type {};

// `Make_x_monotone_2`
// Fallback selected if `Make_x_monotone_2` is not defined in the traits `T` below.
template <typename, typename = std::void_t<>>
struct has_make_x_monotone_2 : std::false_type {};

// Partial specialization selected if `T::Make_x_monotone_2` is defined.
template <typename T>
struct has_make_x_monotone_2<T, std::void_t<typename T::Make_x_monotone_2>> : std::true_type {};

// `Split_2`
// Fallback selected if `Split_2` is not defined in the traits `T` below.
template <typename, typename = std::void_t<>>
struct has_split_2 : std::false_type {};

// Partial specialization selected if `T::Split_2` is defined.
template <typename T>
struct has_split_2<T, std::void_t<typename T::Split_2>> : std::true_type {};

// `Intersect_2`
// Fallback selected if `Intersect_2` is not defined in the traits `T` below.
template <typename, typename = std::void_t<>>
struct has_intersect_2 : std::false_type {};

// Partial specialization selected if `T::Intersect_2` is defined.
template <typename T>
struct has_intersect_2<T, std::void_t<typename T::Intersect_2>> : std::true_type {};

// `Do_intersect_2`
// Fallback selected if `Do_intersect_2` is not defined in the traits `T` below.
template <typename, typename = std::void_t<>>
struct has_do_intersect_2 : std::false_type {};

// Partial specialization selected if `T::Do_intersect_2` is defined.
template <typename T>
struct has_do_intersect_2<T, std::void_t<typename T::Do_intersect_2>> : std::true_type {};

// `Are_mergeable_2`
// Fallback selected if `Are_mergeable_2` is not defined in the traits `T` below.
template <typename, typename = std::void_t<>>
struct has_are_mergeable_2 : std::false_type {};

// Partial specialization selected if `T::Are_mergeable_2` is defined.
template <typename T>
struct has_are_mergeable_2<T, std::void_t<typename T::Are_mergeable_2>> : std::true_type {};

// `Merge_2`
// Fallback selected if `Merge_2` is not defined in the traits `T` below.
template <typename, typename = std::void_t<>>
struct has_merge_2 : std::false_type {};

// Partial specialization selected if `T::Merge_2` is defined.
template <typename T>
struct has_merge_2<T, std::void_t<typename T::Merge_2>> : std::true_type {};

// `Construct_opposite_2`
// Fallback selected if `Construct_opposite_2` is not defined in the traits `T` below.
template <typename, typename = std::void_t<>>
struct has_construct_opposite_2 : std::false_type {};

// Partial specialization selected if `T::Construct_opposite_2` is defined.
template <typename T>
struct has_construct_opposite_2<T, std::void_t<typename T::Construct_opposite_2>> : std::true_type {};

// `Construct_point_2`
// Fallback selected if `Construct_point_2` is not defined in the traits `T` below..
template <typename, typename = std::void_t<>>
struct has_construct_point_2 : std::false_type {};

// Partial specialization selected `T::Construct_point_2` is defined.
template <typename T>
struct has_construct_point_2<T, std::void_t<typename T::Construct_point_2>> : std::true_type {};

/* Fallback selected if `Construct_point_2`, nested in the traits `T` below,
 * does not define an operator that accepts two coordinates.
 */
template <typename, typename = std::void_t<>>
struct has_construct_point_2_xy : std::false_type {};

// Partial specialization selected if `T::Construct_point_2` defines an operator that accepts two coordinates.
template <typename T>
struct has_construct_point_2_xy<T,
                                std::void_t<decltype(std::declval<typename T::Construct_point_2&>()
                                                     (std::declval<detail::any_type>(),
                                                      std::declval<detail::any_type>()))>> : std::true_type {};

// `Construct_x_monotone_curve_2`
// Fallback selected if `Construct_x_monotone_curve_2` is not defined in the traits `T` below.
template <typename, typename = std::void_t<>>
struct has_construct_x_monotone_curve_2 : std::false_type {};

// Partial specialization selected if `T::Construct_x_monotone_curve_2` is defined.
template <typename T>
struct has_construct_x_monotone_curve_2<T, std::void_t<typename T::Construct_x_monotone_curve_2>> : std::true_type {};

// `Construct_curve_2`
// Fallback selected if `Construct_curve_2` is not defined in the traits `T` below.
template <typename, typename = std::void_t<>>
struct has_construct_curve_2 : std::false_type {};

// Partial specialization selected if `T::Construct_curve_2` is defined.
template <typename T>
struct has_construct_curve_2<T, std::void_t<typename T::Construct_curve_2>> : std::true_type {};

// `Compare_endpoints_xy_2`
// Fallback selected if `Compare_endpoints_xy_2` is not defined in the traits `T` below.
template <typename, typename = std::void_t<>>
struct has_compare_endpoints_xy_2 : std::false_type {};

// Partial specialization selected if `T::Compare_endpoints_xy_2` is defined.
template <typename T>
struct has_compare_endpoints_xy_2<T, std::void_t<typename T::Compare_endpoints_xy_2>> : std::true_type {};

// `Approximate_2`
// Fallback selected if `Approximate_2` is not defined in the traits `T` below.
template <typename, typename = std::void_t<>>
struct has_approximate_2 : std::false_type {};

// Partial specialization selected if `T::Approximate_2` is defined.
template <typename T>
struct has_approximate_2<T, std::void_t<typename T::Approximate_2>> : std::true_type {};

/* Fallback selected if `Approximate_2`, nested in the traits `T` below,
 * does not define an operator that approximates a point.
 */
template <typename, typename = std::void_t<>>
struct has_approximate_2_point : std::false_type {};

// Partial specialization selected if `T::Approximate_2` defines an operator that approximates a point.
template <typename T>
struct has_approximate_2_point<T,
                               std::void_t<decltype(std::declval<typename T::Approximate_2&>()
                                                    (std::declval<const typename T::Point_2&>()))>> : std::true_type {};

/* Fallback selected if `Approximate_2`, nested in the traits `T` below,
 * does not define an operator that approximates a bounded \f$x\f$-monotone curve.
 */
template <typename, typename = std::void_t<>>
struct has_approximate_2_xcv : std::false_type {};

/* Partial specialization selected if `T::Approximate_2` defines an operator
 * that approximates a bounded \f$x\f$-monotone curve.
 */
template <typename T>
struct has_approximate_2_xcv<T,
                             std::void_t<decltype(std::declval<typename T::Approximate_2&>()
                                                  (std::declval<const typename T::X_monotone_curve_2&>(),
                                                   std::declval<double>(),
                                                   std::declval<void*>(),
                                                   std::declval<bool>()))>> : std::true_type {};

/* Fallback selected if `Approximate_2`, nested in the traits `T` below,
 * does not define an operator that approximates an ubounded \f$x\f$-monotone curve.
 */
template <typename, typename = std::void_t<>>
struct has_approximate_2_xcv_bounds : std::false_type {};

/* Partial specialization selected if `T::Approximate_2` defines an operator
 * that approximates an unbounded \f$x\f$-monotone curve.
 */
template <typename T>
struct has_approximate_2_xcv_bounds<T,
                                    std::void_t<decltype(std::declval<typename T::Approximate_2&>()
                                                         (std::declval<const typename T::X_monotone_curve_2&>(),
                                                          std::declval<double>(),
                                                          std::declval<void*>(),
                                                          std::declval<Bbox_2>(),
                                                          std::declval<bool>()))>> : std::true_type {};

// `Parameter_space_in_x_2`
// Fallback selected if `Parameter_space_in_x_2` is not defined in the traits `T` below.
template <typename, typename = std::void_t<>>
struct has_parameter_space_in_x_2 : std::false_type {};

// Partial specialization selected if `T::Parameter_space_in_x_2` is defined.
template <typename T>
struct has_parameter_space_in_x_2<T, std::void_t<typename T::Parameter_space_in_x_2>> : std::true_type {};

// `Is_on_x_identification_2`
// Fallback selected if `Is_on_x_identification_2` is not defined in the traits `T` below.
template <typename, typename = std::void_t<>>
struct has_is_on_x_identification_2 : std::false_type {};

// Partial specialization selected if `T::Is_on_x_identification_2` is defined.
template <typename T>
struct has_is_on_x_identification_2<T, std::void_t<typename T::Is_on_x_identification_2>> : std::true_type {};

// `Compare_y_on_boundary_2`
// Fallback selected if `Compare_y_on_boundary_2` is not defined in the traits `T` below.
template <typename, typename = std::void_t<>>
struct has_compare_y_on_boundary_2 : std::false_type {};

// Partial specialization selected if `T::Compare_y_on_boundary_2` is defined.
template <typename T>
struct has_compare_y_on_boundary_2<T, std::void_t<typename T::Compare_y_on_boundary_2>> : std::true_type {};

// `Compare_y_near_boundary_2`
// Fallback selected if `Compare_y_near_boundary_2` is not defined in the traits `T` below.
template <typename, typename = std::void_t<>>
struct has_compare_y_near_boundary_2 : std::false_type {};

// Partial specialization selected if `T::Compare_y_near_boundary_2` is defined.
template <typename T>
struct has_compare_y_near_boundary_2<T, std::void_t<typename T::Compare_y_near_boundary_2>> : std::true_type {};

// `Parameter_space_in_y_2`
// Fallback selected if `Parameter_space_in_y_2` is not defined in the traits `T` below.
template <typename, typename = std::void_t<>>
struct has_parameter_space_in_y_2 : std::false_type {};

// Partial specialization selected if `T::Parameter_space_in_y_2` is defined.
template <typename T>
struct has_parameter_space_in_y_2<T, std::void_t<typename T::Parameter_space_in_y_2>> : std::true_type {};

// `Is_on_y_identification_2`
// Fallback selected if `Is_on_y_identification_2` is not defined in the traits `T` below.
template <typename, typename = std::void_t<>>
struct has_is_on_y_identification_2 : std::false_type {};

// Partial specialization selected if `T::Is_on_y_identification_2` is defined.
template <typename T>
struct has_is_on_y_identification_2<T, std::void_t<typename T::Is_on_y_identification_2>> : std::true_type {};

// `Compare_x_on_boundary_2`
// Fallback selected if `Compare_x_on_boundary_2` is not defined in the traits `T` below.
template <typename, typename = std::void_t<>>
struct has_compare_x_on_boundary_2 : std::false_type {};

// Partial specialization selected if `T::Compare_x_on_boundary_2` is defined.
template <typename T>
struct has_compare_x_on_boundary_2<T, std::void_t<typename T::Compare_x_on_boundary_2>> : std::true_type {};

/* Fallback selected if `Compare_x_on_boundary_2`, nested in the traits `T` below,
 * does not define an operator that accepts two points.
 */
template <typename, typename = std::void_t<>>
struct has_compare_x_on_boundary_2_points : std::false_type {};

// Partial specialization selected if `T::Compare_x_on_boundary_2` defines an operator that accepts two points.
template <typename T>
struct has_compare_x_on_boundary_2_points<T,
                                          std::void_t<decltype(std::declval<typename T::Compare_x_on_boundary_2>()
                                                               (std::declval<const typename T::Point_2&>(),
                                                                std::declval<const typename T::Point_2&>()))>> :
    std::true_type {};

/* Fallback selected if `Compare_x_on_boundary_2`, nested in the traits `T` below,
 * does not define an operator that accepts a point and a curve end.
 */
template <typename, typename = std::void_t<>>
struct has_compare_x_on_boundary_2_point_curve_end : std::false_type {};

/* Partial specialization selected if `T::Compare_x_on_boundary_2` defines an
 * operator that accepts a point and a curve end.
 */
template <typename T>
struct has_compare_x_on_boundary_2_point_curve_end
<T, std::void_t<decltype(std::declval<typename T::Compare_x_on_boundary_2>()(
     std::declval<const typename T::Point_2&>(),
     std::declval<const typename T::X_monotone_curve_2&>(),
     std::declval<Arr_curve_end>()))>> : std::true_type {};

/* Fallback selected if `Compare_x_on_boundary_2`, nested in the traits `T` below,
 * does not define an operator that accepts two curve ends.
 */
template <typename, typename = std::void_t<>>
struct has_compare_x_on_boundary_2_curve_ends : std::false_type {};

/* Partial specialization selected if `T::Compare_x_on_boundary_2` defines an
 * operator that accepts two curve ends.
 */
template <typename T>
struct has_compare_x_on_boundary_2_curve_ends
<T, std::void_t<decltype(std::declval<typename T::Compare_x_on_boundary_2>()(
     std::declval<const typename T::X_monotone_curve_2&>(),
     std::declval<Arr_curve_end>(),
     std::declval<const typename T::X_monotone_curve_2&>(),
     std::declval<Arr_curve_end>()))>> : std::true_type {};

// `Compare_x_near_boundary_2`
// Fallback selected if `Compare_x_near_boundary_2` is not defined in the traits `T` below.
template <typename, typename = std::void_t<>>
struct has_compare_x_near_boundary_2 : std::false_type {};

// Partial specialization selected if `T::Compare_x_near_boundary_2` is defined.
template <typename T>
struct has_compare_x_near_boundary_2<T, std::void_t<typename T::Compare_x_near_boundary_2>> : std::true_type {};

}

#endif
