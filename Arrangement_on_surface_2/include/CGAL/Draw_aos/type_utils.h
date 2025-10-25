// Copyright (c) 2025
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
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Shepard Liu       <shepard0liu@gmail.com>

#ifndef CGAL_DRAW_AOS_TYPE_UTILS_H
#define CGAL_DRAW_AOS_TYPE_UTILS_H

#include <limits>
#include <vector>
#include <type_traits>
#include <cmath>

#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>

namespace CGAL {
namespace draw_aos {

enum class Boundary_side {
  Top = 0,
  Left = 1,
  Bottom = 2,
  Right = 3,
  None = -1,
};

template <typename, typename = std::void_t<>>
struct has_approximate_2_object : std::false_type {};

template <typename Gt>
struct has_approximate_2_object<Gt, std::void_t<decltype(std::declval<Gt>().approximate_2_object())>> : std::true_type
{};

// Detect whether Gt has defined a member function approximate_2_object()
template <typename Gt>
inline constexpr bool has_approximate_2_object_v = has_approximate_2_object<Gt>::value;

template <typename, typename, typename = std::void_t<>>
struct has_approximate_point : std::false_type {};

template <typename Gt, typename A>
struct has_approximate_point<Gt,
                             A,
                             std::void_t<decltype(std::declval<A>()(std::declval<const typename Gt::Point_2&>()))>> :
    std::true_type
{};

// Detect whether A has operator()(const Gt::Point_2&)
template <typename Gt, typename A>
inline constexpr bool has_approximate_point_v = has_approximate_point<Gt, A>::value;

template <typename, typename, typename, typename = std::void_t<>>
struct has_approximate_xcv : std::false_type {};

template <typename Gt, typename A, typename O>
struct has_approximate_xcv
<Gt,
 A,
 O,
 std::void_t<decltype(std::declval<A&>()(std::declval<const typename Gt::X_monotone_curve_2&>(),
                                         std::declval<double>(),
                                         std::declval<O>(),
                                         std::declval<bool>()))>> : std::true_type
{};

// Detect whether A has operator()(const Gt::X_monotone_curve_2&, double, OutputIterator, bool)?
template <typename Gt, typename A>
constexpr bool has_approximate_xcv_v = has_approximate_xcv<Gt, A, void*>::value;

template <typename, typename, typename, typename = std::void_t<>>
struct has_approximate_xcv_with_bounds : std::false_type
{};

template <typename Gt, typename A, typename O>
struct has_approximate_xcv_with_bounds
<Gt,
 A,
 O,
 std::void_t<decltype(std::declval<A&>()(std::declval<const typename Gt::X_monotone_curve_2&>(),
                                         std::declval<double>(),
                                         std::declval<O>(),
                                         std::declval<Bbox_2>(),
                                         std::declval<bool>()))>> : std::true_type
{};

// Detect whether A has operator()(const X_monotone_curve&, double, OutputIterator, Bbox_2, bool)
template <typename Gt, typename A>
inline constexpr bool has_approximate_xcv_with_bounds_v = has_approximate_xcv_with_bounds<Gt, A, void*>::value;

// Detect whether a geometry traits has all the necessary types and functions for approximation
template <typename Gt>
constexpr bool has_approximate_traits_v =
  has_approximate_2_object_v<Gt> && has_approximate_point_v<Gt, typename Gt::Approximate_2> &&
  (has_approximate_xcv_v<Gt, typename Gt::Approximate_2> ||
   has_approximate_xcv_with_bounds_v<Gt, typename Gt::Approximate_2>);

template <typename Gt, typename = std::void_t<>>
struct has_is_in_x_range : std::false_type
{};

template <typename Gt>
struct has_is_in_x_range<Gt,
                         std::void_t<decltype(std::declval<typename Gt::X_monotone_curve_2>().is_in_x_range
                                              (std::declval<const typename Gt::Point_2&>()))>> : std::true_type
{};

// Detect whether Gt::X_monotone_curve_2 has a member function bool is_in_x_range(const Gt::Point_2&)
template <typename Gt>
inline constexpr bool has_is_in_x_range_v = has_is_in_x_range<Gt>::value;

// Detect whether Gt is or derives from Arr_geodesic_arc_on_sphere_traits_2<*, *, *>
template <typename Gt>
struct is_or_derived_from_agas {
private:
  template <typename Kernel_, int AtanX, int AtanY>
  static std::true_type test(const Arr_geodesic_arc_on_sphere_traits_2<Kernel_, AtanX, AtanY>*);

  static std::false_type test(...);

public:
  static constexpr bool value = decltype(test(static_cast<const Gt*>(nullptr)))::value;
};

template <typename Gt>
inline constexpr bool is_or_derived_from_agas_v = is_or_derived_from_agas<Gt>::value;

// Detect whether T is or derives from a geometry traits on curved surfaces
template <typename Gt>
inline constexpr bool is_or_derived_from_curved_surf_traits_v = is_or_derived_from_agas_v<Gt>;

// Static helpers to get template arguments from a geometry traits
template <typename Gt>
struct tmpl_args {};

template <typename Kernel_, int AtanX, int AtanY>
struct tmpl_args<Arr_geodesic_arc_on_sphere_traits_2<Kernel_, AtanX, AtanY>> {
  using Kernel = Kernel_;
  static constexpr int atan_x = AtanX;
  static constexpr int atan_y = AtanY;
};

/*! \brief Approximation data types
 *
 * \tparam Gt Geometry traits
 */
template <typename Gt>
class Arr_approximate_traits {
  using Geom_traits = Gt;

  template <typename P, typename I>
  struct Triangle_soup_ {
    using Index = I;
    using Triangle = std::array<Index, 3>;
    using Point = P;

    std::vector<Point> points;
    std::vector<Triangle> triangles;
  };

public:
  using Approx_point = typename Geom_traits::Approximate_point_2;
  using Approx_nt = typename Geom_traits::Approximate_number_type;
  using Approx_kernel = typename Geom_traits::Approximate_kernel;

  // 2D parameter space point
  using Point = typename Approx_kernel::Point_2;
  using Polyline = std::vector<Point>;

  using Triangle_soup = Triangle_soup_<Point, int>;

  // A null point with NaN coordinates. Use ::is_null(pt) to check if a point is null.
  inline static const Point Null_point =
      Point(std::numeric_limits<Approx_nt>::signaling_NaN(), std::numeric_limits<Approx_nt>::signaling_NaN());
  static bool is_null(Point pt) { return std::isnan(pt.x()) || std::isnan(pt.y()); }
};

/*!
 * \brief Functor to construct a Point_2 from an Approximate_point_2.
 *
 * \tparam Gt Geometry traits
 */
template <typename Gt>
class Construct_gt_point_2 {
  using Approx_traits = Arr_approximate_traits<Gt>;
  using Approx_point = typename Approx_traits::Approx_point;
  using Gt_point = typename Gt::Point_2;

public:
  Gt_point operator()(const Approx_point& pt) const { return Gt_point(pt.x(), pt.y()); }
};

// Specialization for Arr_geodesic_arc_on_sphere_traits_2
template <typename Kernel, int AtanX, int AtanY>
class Construct_gt_point_2<Arr_geodesic_arc_on_sphere_traits_2<Kernel, AtanX, AtanY>> {
  using Geom_traits = Arr_geodesic_arc_on_sphere_traits_2<Kernel, AtanX, AtanY>;
  using Approx_traits = Arr_approximate_traits<Arr_geodesic_arc_on_sphere_traits_2<Kernel, AtanX, AtanY>>;
  using Approx_point = typename Approx_traits::Approx_point;
  using Gt_point = typename Geom_traits::Point_2;

public:
  Gt_point operator()(const Approx_point& pt) const {
    using Direction_3 = typename Kernel::Direction_3;
    return Gt_point(Direction_3(pt.dx(), pt.dy(), pt.dz()),
                    static_cast<typename Gt_point::Location_type>(pt.location()));
  }
};

} // namespace draw_aos
} // namespace CGAL

#endif // CGAL_DRAW_AOS_TYPE_UTILS_H
