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

#include <CGAL/license/Arrangement_on_surface_2.h>

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

// Obtain the traits that a traits class is built upon, if any:
// - A traits decorator, e.g., `Arr_tracing_traits_2` or `Arr_counting_traits_2`,
//   exposes the traits it decorates as the nested type `Decorated_traits_2`.
// - A polycurve traits, e.g., `Arr_polycurve_traits_2` or
//   `Arr_polyline_traits_2`, exposes the traits of its subcurves as the nested
//   type `Subcurve_traits_2`.
// If Gt defines both types, `Decorated_traits_2` takes precedence. If Gt
// defines neither, it is not built upon another traits.
template <typename Gt, typename = void>
struct subcurve_traits_of {};

template <typename Gt>
struct subcurve_traits_of<Gt, std::void_t<typename Gt::Subcurve_traits_2>>
{ using type = typename Gt::Subcurve_traits_2; };

template <typename Gt, typename = void>
struct underlying_traits_of : subcurve_traits_of<Gt> {};

template <typename Gt>
struct underlying_traits_of<Gt, std::void_t<typename Gt::Decorated_traits_2>>
{ using type = typename Gt::Decorated_traits_2; };

// Obtain the last (innermost) traits in a chain of traits classes, each built
// upon the next one (see `underlying_traits_of`), e.g.,
// Arr_counting_traits_2<Arr_polyline_traits_2<Arr_polycurve_traits_2<...>>>.
// If Gt is not built upon another traits, the result is Gt itself.
template <typename Gt, typename = void>
struct innermost_traits { using type = Gt; };

template <typename Gt>
struct innermost_traits<Gt, std::void_t<typename underlying_traits_of<Gt>::type>> :
    innermost_traits<typename underlying_traits_of<Gt>::type>
{};

template <typename Gt>
using innermost_traits_t = typename innermost_traits<Gt>::type;

// Detect whether Gt is a traits decorator, i.e., it defines `Decorated_traits_2`.
template <typename Gt, typename = void>
struct is_traits_decorator : std::false_type {};

template <typename Gt>
struct is_traits_decorator<Gt, std::void_t<typename Gt::Decorated_traits_2>> : std::true_type {};

template <typename Gt>
inline constexpr bool is_traits_decorator_v = is_traits_decorator<Gt>::value;

// Detect whether Gt is built upon another traits (see `underlying_traits_of`).
template <typename Gt, typename = void>
struct has_underlying_traits : std::false_type {};

template <typename Gt>
struct has_underlying_traits<Gt, std::void_t<typename underlying_traits_of<Gt>::type>> : std::true_type {};

template <typename Gt>
inline constexpr bool has_underlying_traits_v = has_underlying_traits<Gt>::value;

// Obtain the traits object that the given traits object is built upon: the
// decorated traits of a decorator or the subcurve traits of a polycurve traits.
template <typename Gt, std::enable_if_t<has_underlying_traits_v<Gt>, int> = 0>
const typename underlying_traits_of<Gt>::type& underlying_traits(const Gt& traits) {
  if constexpr (is_traits_decorator_v<Gt>) return traits.traits();
  else return *traits.subcurve_traits_2();
}

// Detect whether the innermost traits of Gt (see `innermost_traits`) is or
// derives from
// Arr_geodesic_arc_on_sphere_traits_2<*, *, *>
template <typename Gt>
struct is_or_derived_from_agas {
private:
  template <typename Kernel_, int AtanX, int AtanY>
  static std::true_type test(const Arr_geodesic_arc_on_sphere_traits_2<Kernel_, AtanX, AtanY>*);

  static std::false_type test(...);

public:
  static constexpr bool value = decltype(test(static_cast<const innermost_traits_t<Gt>*>(nullptr)))::value;
};

template <typename Gt>
inline constexpr bool is_or_derived_from_agas_v = is_or_derived_from_agas<Gt>::value;

// Detect whether T is or derives from a geometry traits on curved surfaces
template <typename Gt>
inline constexpr bool is_or_derived_from_curved_surf_traits_v = is_or_derived_from_agas_v<Gt>;

// Static helpers to get template arguments from a geometry traits
template <typename Kernel_, int AtanX, int AtanY>
struct Agas_tmpl_args {
  using Kernel = Kernel_;
  static constexpr int atan_x = AtanX;
  static constexpr int atan_y = AtanY;
};

struct No_tmpl_args {};

template <typename Kernel_, int AtanX, int AtanY>
Agas_tmpl_args<Kernel_, AtanX, AtanY> agas_tmpl_args(const Arr_geodesic_arc_on_sphere_traits_2<Kernel_, AtanX, AtanY>*);

No_tmpl_args agas_tmpl_args(...);

// Obtain the template arguments of the innermost traits of Gt (see
// `innermost_traits`), provided that it is or derives from
// Arr_geodesic_arc_on_sphere_traits_2<*, *, *>.
template <typename Gt>
struct tmpl_args : decltype(agas_tmpl_args(static_cast<const innermost_traits_t<Gt>*>(nullptr))) {};

/*! constructs the identification curve of an arrangement on the sphere, that
 * is, the x-monotone curve of the given traits that lies on the identification
 * meridian from the south pole to the north pole.
 * The arc is constructed by the innermost traits, which is (or derives from)
 * `Arr_geodesic_arc_on_sphere_traits_2`, as the construction of a curve from
 * its endpoints and the normal of its plane is specific to this traits. Then,
 * it is lifted through the chain of traits: a decorator uses the curve type
 * of the traits it decorates, and a polycurve traits constructs an x-monotone
 * polycurve that consists of a single subcurve.
 */
template <typename Gt>
typename Gt::X_monotone_curve_2 construct_identification_curve(const Gt& traits) {
  static_assert(is_or_derived_from_agas_v<Gt>,
                "The innermost traits must be (or derive from) Arr_geodesic_arc_on_sphere_traits_2.");
  if constexpr (has_underlying_traits_v<Gt>) {
    auto xcv = construct_identification_curve(underlying_traits(traits));
    if constexpr (is_traits_decorator_v<Gt>) return xcv;
    else return traits.construct_x_monotone_curve_2_object()(xcv);
  }
  else {
    using Direction_3 = typename Gt::Direction_3;
    using Point_2 = typename Gt::Point_2;
    using Args = tmpl_args<Gt>;
    Point_2 src(Direction_3(0, 0, -1), Point_2::MIN_BOUNDARY_LOC);
    Point_2 trg(Direction_3(0, 0, 1), Point_2::MAX_BOUNDARY_LOC);
    return traits.construct_x_monotone_curve_2_object()(src, trg, Direction_3(Args::atan_y, -Args::atan_x, 0));
  }
}

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
  using Gt_approximate_2 = typename Geom_traits::Approximate_2;
  using Approx_point = typename Gt_approximate_2::Approximate_point_2;
  using Approx_nt = typename Gt_approximate_2::Approximate_number_type;
  using Approx_kernel = typename Gt_approximate_2::Approximate_kernel;

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
template <typename Gt, typename = void>
class Construct_gt_point_2 {
  using Approx_traits = Arr_approximate_traits<Gt>;
  using Approx_point = typename Approx_traits::Approx_point;
  using Gt_point = typename Gt::Point_2;

public:
  Gt_point operator()(const Approx_point& pt) const { return Gt_point(pt.x(), pt.y()); }
};

// Specialization for traits whose innermost traits (see `innermost_traits`) is
// or derives from Arr_geodesic_arc_on_sphere_traits_2
template <typename Gt>
class Construct_gt_point_2<Gt, std::enable_if_t<is_or_derived_from_agas_v<Gt>>> {
  using Geom_traits = Gt;
  using Approx_traits = Arr_approximate_traits<Geom_traits>;
  using Approx_point = typename Approx_traits::Approx_point;
  using Gt_point = typename Geom_traits::Point_2;

public:
  Gt_point operator()(const Approx_point& pt) const {
    using Direction_3 = typename tmpl_args<Geom_traits>::Kernel::Direction_3;
    return Gt_point(Direction_3(pt.dx(), pt.dy(), pt.dz()),
                    static_cast<typename Gt_point::Location_type>(pt.location()));
  }
};

} // namespace draw_aos
} // namespace CGAL

#endif
