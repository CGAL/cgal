// Copyright (c) 2020 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Dmitry Anisimov
//

#ifndef CGAL_WEIGHTS_INTERNAL_UTILS_H
#define CGAL_WEIGHTS_INTERNAL_UTILS_H

// STL includes.
#include <cmath>
#include <string>
#include <memory>
#include <vector>
#include <utility>
#include <iterator>
#include <iostream>
#include <cassert>

// Boost headers.
#include <boost/mpl/has_xxx.hpp>
#include <CGAL/boost/graph/helpers.h>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/property_map.h>
#include <CGAL/Point_2.h>
#include <CGAL/Point_3.h>

namespace CGAL {
namespace Weights {
namespace internal {

  // Sqrt helpers.
  template<typename GeomTraits>
  class Default_sqrt {

  private:
    using Traits = GeomTraits;
    using FT = typename Traits::FT;

  public:
    FT operator()(const FT value) const {
      return static_cast<FT>(
        CGAL::sqrt(CGAL::to_double(CGAL::abs(value))));
    }
  };

  BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_nested_type_Sqrt, Sqrt, false)

  // Case: do_not_use_default = false.
  template<typename GeomTraits,
  bool do_not_use_default = Has_nested_type_Sqrt<GeomTraits>::value>
  class Get_sqrt {

  public:
    using Traits = GeomTraits;
    using Sqrt = Default_sqrt<Traits>;

    static Sqrt sqrt_object(const Traits& ) {
      return Sqrt();
    }
  };

  // Case: do_not_use_default = true.
  template<typename GeomTraits>
  class Get_sqrt<GeomTraits, true> {

  public:
    using Traits = GeomTraits;
    using Sqrt = typename Traits::Sqrt;

    static Sqrt sqrt_object(const Traits& traits) {
      return traits.sqrt_object();
    }
  };

  // Normalize values.
  template<typename FT>
  void normalize(std::vector<FT>& values) {

    FT sum = FT(0);
    for (const FT& value : values) {
      sum += value;
    }

    CGAL_assertion(sum != FT(0));
    if (sum == FT(0)) {
      return;
    }

    const FT inv_sum = FT(1) / sum;
    for (FT& value : values) {
      value *= inv_sum;
    }
  }

  // Raises value to the power.
  template<typename GeomTraits>
  typename GeomTraits::FT power(
    const GeomTraits&,
    const typename GeomTraits::FT value,
    const typename GeomTraits::FT p) {

    using FT = typename GeomTraits::FT;
    const double base = CGAL::to_double(value);
    const double exp = CGAL::to_double(p);
    return static_cast<FT>(std::pow(base, exp));
  }

  // Computes distance between two 2D points.
  template<typename GeomTraits>
  typename GeomTraits::FT distance_2(
    const GeomTraits& traits,
    const typename GeomTraits::Point_2& p,
    const typename GeomTraits::Point_2& q) {

    using Get_sqrt = Get_sqrt<GeomTraits>;
    const auto sqrt = Get_sqrt::sqrt_object(traits);

    const auto squared_distance_2 =
      traits.compute_squared_distance_2_object();
    return sqrt(squared_distance_2(p, q));
  }

  // Computes length of a 2D vector.
  template<typename GeomTraits>
  typename GeomTraits::FT length_2(
    const GeomTraits& traits,
    const typename GeomTraits::Vector_2& v) {

    using Get_sqrt = Get_sqrt<GeomTraits>;
    const auto sqrt = Get_sqrt::sqrt_object(traits);

    const auto squared_length_2 =
      traits.compute_squared_length_2_object();
    return sqrt(squared_length_2(v));
  }

  // Normalizes a 2D vector.
  template<typename GeomTraits>
  void normalize_2(
    const GeomTraits& traits,
    typename GeomTraits::Vector_2& v) {

    using FT = typename GeomTraits::FT;
    const FT length = length_2(traits, v);
    CGAL_assertion(length != FT(0));
    if (length == FT(0)) {
      return;
    }
    v /= length;
  }

  // Computes cotanget between two 2D vectors.
  template<typename GeomTraits>
  typename GeomTraits::FT cotangent_2(
    const GeomTraits& traits,
    const typename GeomTraits::Point_2& p,
    const typename GeomTraits::Point_2& q,
    const typename GeomTraits::Point_2& r) {

    using FT = typename GeomTraits::FT;
    const auto dot_product_2 =
      traits.compute_scalar_product_2_object();
    const auto cross_product_2 =
      traits.compute_determinant_2_object();
    const auto construct_vector_2 =
      traits.construct_vector_2_object();

    const auto v1 = construct_vector_2(q, r);
    const auto v2 = construct_vector_2(q, p);

    const FT dot = dot_product_2(v1, v2);
    const FT cross = cross_product_2(v1, v2);

    const FT length = CGAL::abs(cross);
    // CGAL_assertion(length != FT(0)); not really necessary
    if (length != FT(0)) {
      return dot / length;
    } else {
      return FT(0); // undefined
    }
  }

  // Computes tanget between two 2D vectors.
  template<typename GeomTraits>
  typename GeomTraits::FT tangent_2(
    const GeomTraits& traits,
    const typename GeomTraits::Point_2& p,
    const typename GeomTraits::Point_2& q,
    const typename GeomTraits::Point_2& r) {

    using FT = typename GeomTraits::FT;
    const auto dot_product_2 =
      traits.compute_scalar_product_2_object();
    const auto cross_product_2 =
      traits.compute_determinant_2_object();
    const auto construct_vector_2 =
      traits.construct_vector_2_object();

    const auto v1 = construct_vector_2(q, r);
    const auto v2 = construct_vector_2(q, p);

    const FT dot = dot_product_2(v1, v2);
    const FT cross = cross_product_2(v1, v2);

    const FT length = CGAL::abs(cross);
    // CGAL_assertion(dot != FT(0)); not really necessary
    if (dot != FT(0)) {
      return length / dot;
    } else {
      return FT(0); // undefined
    }
  }

  // Computes distance between two 3D points.
  template<typename GeomTraits>
  typename GeomTraits::FT distance_3(
    const GeomTraits& traits,
    const typename GeomTraits::Point_3& p,
    const typename GeomTraits::Point_3& q) {

    using Get_sqrt = Get_sqrt<GeomTraits>;
    const auto sqrt = Get_sqrt::sqrt_object(traits);

    const auto squared_distance_3 =
      traits.compute_squared_distance_3_object();
    return sqrt(squared_distance_3(p, q));
  }

  // Computes length of a 3D vector.
  template<typename GeomTraits>
  typename GeomTraits::FT length_3(
    const GeomTraits& traits,
    const typename GeomTraits::Vector_3& v) {

    using Get_sqrt = Get_sqrt<GeomTraits>;
    const auto sqrt = Get_sqrt::sqrt_object(traits);

    const auto squared_length_3 =
      traits.compute_squared_length_3_object();
    return sqrt(squared_length_3(v));
  }

  // Normalizes a 3D vector.
  template<typename GeomTraits>
  void normalize_3(
    const GeomTraits& traits,
    typename GeomTraits::Vector_3& v) {

    using FT = typename GeomTraits::FT;
    const FT length = length_3(traits, v);
    CGAL_assertion(length != FT(0));
    if (length == FT(0)) {
      return;
    }
    v /= length;
  }

  // Computes cotanget between two 3D vectors.
  template<typename GeomTraits>
  typename GeomTraits::FT cotangent_3(
    const GeomTraits& traits,
    const typename GeomTraits::Point_3& p,
    const typename GeomTraits::Point_3& q,
    const typename GeomTraits::Point_3& r) {

    using FT = typename GeomTraits::FT;
    const auto dot_product_3 =
      traits.compute_scalar_product_3_object();
    const auto cross_product_3 =
      traits.construct_cross_product_vector_3_object();
    const auto construct_vector_3 =
      traits.construct_vector_3_object();

    const auto v1 = construct_vector_3(q, r);
    const auto v2 = construct_vector_3(q, p);

    const FT dot = dot_product_3(v1, v2);
    const auto cross = cross_product_3(v1, v2);

    const FT length = length_3(traits, cross);
    // TODO:
    // Not really necessary: since we handle case length = 0. Does this case happen?
    // Yes, e.g. in Surface Parameterization tests. Does it affect the results?
    // In current applications, not really.
    // CGAL_assertion(length != FT(0));
    if (length != FT(0)) {
      return dot / length;
    } else {
      return FT(0); // undefined
    }
  }

  // Computes tanget between two 3D vectors.
  template<typename GeomTraits>
  typename GeomTraits::FT tangent_3(
    const GeomTraits& traits,
    const typename GeomTraits::Point_3& p,
    const typename GeomTraits::Point_3& q,
    const typename GeomTraits::Point_3& r) {

    using FT = typename GeomTraits::FT;
    const auto dot_product_3 =
      traits.compute_scalar_product_3_object();
    const auto cross_product_3 =
      traits.construct_cross_product_vector_3_object();
    const auto construct_vector_3 =
      traits.construct_vector_3_object();

    const auto v1 = construct_vector_3(q, r);
    const auto v2 = construct_vector_3(q, p);

    const FT dot = dot_product_3(v1, v2);
    const auto cross = cross_product_3(v1, v2);

    const FT length = length_3(traits, cross);
    // CGAL_assertion(dot != FT(0)); not really necessary
    if (dot != FT(0)) {
      return length / dot;
    } else {
      return FT(0); // undefined
    }
  }

  // Computes 3D angle between two vectors.
  template<typename GeomTraits>
  double angle_3(
    const GeomTraits& traits,
    const typename GeomTraits::Vector_3& v1,
    const typename GeomTraits::Vector_3& v2) {

    const auto dot_product_3 =
      traits.compute_scalar_product_3_object();
    const double dot =
      CGAL::to_double(dot_product_3(v1, v2));

    double angle_rad = 0.0;
    if (dot < -1.0) {
      angle_rad = std::acos(-1.0);
    } else if (dot > 1.0) {
      angle_rad = std::acos(+1.0);
    } else {
      angle_rad = std::acos(dot);
    }
    return angle_rad;
  }

  // Rotates a 3D point around axis.
  template<typename GeomTraits>
  typename GeomTraits::Point_3 rotate_point_3(
    const GeomTraits&,
    const double angle_rad,
    const typename GeomTraits::Vector_3& axis,
    const typename GeomTraits::Point_3& query) {

    using FT = typename GeomTraits::FT;
    using Point_3 = typename GeomTraits::Point_3;

    const FT c = static_cast<FT>(std::cos(angle_rad));
    const FT s = static_cast<FT>(std::sin(angle_rad));
    const FT C = FT(1) - c;

    const auto x = axis.x();
    const auto y = axis.y();
    const auto z = axis.z();

    return Point_3(
      (x * x * C + c)     * query.x() +
      (x * y * C - z * s) * query.y() +
      (x * z * C + y * s) * query.z(),
      (y * x * C + z * s) * query.x() +
      (y * y * C + c)     * query.y() +
      (y * z * C - x * s) * query.z(),
      (z * x * C - y * s) * query.x() +
      (z * y * C + x * s) * query.y() +
      (z * z * C + c)     * query.z());
  }

  // Computes two 3D orthogonal base vectors wrt a given normal.
  template<typename GeomTraits>
  void orthogonal_bases_3(
    const GeomTraits& traits,
    const typename GeomTraits::Vector_3& normal,
    typename GeomTraits::Vector_3& b1,
    typename GeomTraits::Vector_3& b2) {

    using Vector_3 = typename GeomTraits::Vector_3;
    const auto cross_product_3 =
      traits.construct_cross_product_vector_3_object();

    const auto nx = normal.x();
    const auto ny = normal.y();
    const auto nz = normal.z();

    if (CGAL::abs(nz) >= CGAL::abs(ny)) {
      b1 = Vector_3(nz, 0, -nx);
    } else {
      b1 = Vector_3(ny, -nx, 0);
    }
    b2 = cross_product_3(normal, b1);

    normalize_3(traits, b1);
    normalize_3(traits, b2);
  }

  // Converts a 3D point into a 2D point wrt to a given plane.
  template<typename GeomTraits>
  typename GeomTraits::Point_2 to_2d(
    const GeomTraits& traits,
    const typename GeomTraits::Vector_3& b1,
    const typename GeomTraits::Vector_3& b2,
    const typename GeomTraits::Point_3& origin,
    const typename GeomTraits::Point_3& query) {

    using Point_2  = typename GeomTraits::Point_2;
    const auto dot_product_3 =
      traits.compute_scalar_product_3_object();
    const auto construct_vector_3 =
      traits.construct_vector_3_object();

    const auto v = construct_vector_3(origin, query);
    const auto x = dot_product_3(b1, v);
    const auto y = dot_product_3(b2, v);
    return Point_2(x, y);
  }

  // Flattening.

  // \cgalFigureBegin{flattening, flattening.svg}
  //   The non-planar configuration (top) is flattened to the planar configuration (bottom).
  // \cgalFigureEnd

  // When computing weights for a query point \f$q\f$ with respect to its neighbors
  // \f$p_0\f$, \f$p_1\f$, and \f$p_2\f$, the local configuration is a quadrilateral
  // [\f$p_0\f$, \f$p_1\f$, \f$p_2\f$, \f$q\f$] or two connected triangles [\f$q\f$, \f$p_0\f$, \f$p_1\f$]
  // and [\f$q\f$, \f$p_1\f$, \f$p_2\f$]. When working in 3D, these triangles are not
  // necessarily coplanar, in other words, they do not belong to the same common plane.
  // When they are not coplanar, they can be made coplanar through the process called *flattening* (see the Figure above),
  // however the latter introduces a distortion because the weights are computed with respect to the
  // flattened configuration rather than to the original non-flat configuration.

  // \subsection Weights_Examples_ProjectionTraits Computing 2D Weights in 3D

  // If you have a 2D polygon in 3D plane that is not an XY plane, you can still compute
  // the 2D weights, however you need to provide a special projection traits class.
  // The common plane that is used in this example is projectable to the XY plane. We first
  // compute `Mean_value_weights_2` for a 3D polygon in this plane. We then also show how to use
  // the projection traits to compute the \ref PkgWeightsRefWachspressWeights "2D Wachspress weight"
  // for 3D points which are not strictly coplanar.

  // \cgalExample{Weights/projection_traits.cpp}

  // Example of flattening:

  // 3D configuration.
  // const Point_3 p0(0, 1, 1);
  // const Point_3 p1(2, 0, 1);
  // const Point_3 p2(7, 1, 1);
  // const Point_3 q0(3, 1, 1);

  // Choose a type of the weight:
  // e.g. 0 - Wachspress (WP) weight.
  // const FT wp = FT(0);

  // Compute WP weights for q1 which is not on the plane [p0, p1, p2].

  // Point_3 q1(3, 1, 2);
  // std::cout << "3D wachspress (WP, q1): ";
  // std::cout << CGAL::Weights::three_point_family_weight(p0, p1, p2, q1, wp) << std::endl;

  // Converge q1 towards q0 that is we flatten the configuration.
  // We also compare the result with the authalic weight.

  // std::cout << "Converge q1 to q0: " << std::endl;
  // for (FT x = FT(0); x <= FT(1); x += step) {
  //   std::cout << "3D wachspress/authalic: ";
  //   q1 = Point_3(3, 1, FT(2) - x);
  //   std::cout << CGAL::Weights::three_point_family_weight(p0, p1, p2, q1, wp) << "/";
  //   std::cout << CGAL::Weights::authalic_weight(p0, p1, p2, q1) << std::endl;
  // }

  // Flattens an arbitrary quad into a planar quad.
  template<typename GeomTraits>
  void flatten(
    const GeomTraits& traits,
    const typename GeomTraits::Point_3& t, // prev neighbor/vertex/point
    const typename GeomTraits::Point_3& r, // curr neighbor/vertex/point
    const typename GeomTraits::Point_3& p, // next neighbor/vertex/point
    const typename GeomTraits::Point_3& q, // query point
    typename GeomTraits::Point_2& tf,
    typename GeomTraits::Point_2& rf,
    typename GeomTraits::Point_2& pf,
    typename GeomTraits::Point_2& qf) {

    // std::cout << std::endl;
    using Point_3 = typename GeomTraits::Point_3;
    using Vector_3 = typename GeomTraits::Vector_3;

    const auto cross_product_3 =
      traits.construct_cross_product_vector_3_object();
    const auto construct_vector_3 =
      traits.construct_vector_3_object();
    const auto centroid_3 =
      traits.construct_centroid_3_object();

    // Compute centroid.
    const auto center = centroid_3(t, r, p, q);
    // std::cout << "centroid: " << center << std::endl;

    // Translate.
    const Point_3 t1 = Point_3(
      t.x() - center.x(), t.y() - center.y(), t.z() - center.z());
    const Point_3 r1 = Point_3(
      r.x() - center.x(), r.y() - center.y(), r.z() - center.z());
    const Point_3 p1 = Point_3(
      p.x() - center.x(), p.y() - center.y(), p.z() - center.z());
    const Point_3 q1 = Point_3(
      q.x() - center.x(), q.y() - center.y(), q.z() - center.z());

    // std::cout << "translated t1: " << t1 << std::endl;
    // std::cout << "translated r1: " << r1 << std::endl;
    // std::cout << "translated p1: " << p1 << std::endl;
    // std::cout << "translated q1: " << q1 << std::endl;

    // Middle axis.
    auto ax = construct_vector_3(q1, r1);
    normalize_3(traits, ax);

    // Prev and next vectors.
    auto v1 = construct_vector_3(q1, t1);
    auto v2 = construct_vector_3(q1, p1);

    normalize_3(traits, v1);
    normalize_3(traits, v2);

    // Two triangle normals.
    auto n1 = cross_product_3(v1, ax);
    auto n2 = cross_product_3(ax, v2);

    normalize_3(traits, n1);
    normalize_3(traits, n2);

    // std::cout << "normal n1: " << n1 << std::endl;
    // std::cout << "normal n2: " << n2 << std::endl;

    // Angle between two normals.
    const double angle_rad = angle_3(traits, n1, n2);
    // std::cout << "angle deg n1 <-> n2: " << angle_rad * 180.0 / CGAL_PI << std::endl;

    // Rotate p1 around ax so that it lands onto the plane [q1, t1, r1].
    const auto& t2 = t1;
    const auto& r2 = r1;
    const auto  p2 = rotate_point_3(traits, angle_rad, ax, p1);
    const auto& q2 = q1;
    // std::cout << "rotated p2: " << p2 << std::endl;

    // Compute orthogonal base vectors.
    Vector_3 b1, b2;
    const auto& normal = n1;
    orthogonal_bases_3(traits, normal, b1, b2);

    // const auto angle12 = angle_3(traits, b1, b2);
    // std::cout << "angle deg b1 <-> b2: " << angle12 * 180.0 / CGAL_PI << std::endl;

    // Flatten a quad.
    const auto& origin = q2;
    tf = to_2d(traits, b1, b2, origin, t2);
    rf = to_2d(traits, b1, b2, origin, r2);
    pf = to_2d(traits, b1, b2, origin, p2);
    qf = to_2d(traits, b1, b2, origin, q2);

    // std::cout << "flattened qf: " << qf << std::endl;
    // std::cout << "flattened tf: " << tf << std::endl;
    // std::cout << "flattened rf: " << rf << std::endl;
    // std::cout << "flattened pf: " << pf << std::endl;

    // std::cout << "A1: " << area_2(traits, rf, qf, pf) << std::endl;
    // std::cout << "A2: " << area_2(traits, pf, qf, rf) << std::endl;
    // std::cout << "C: "  << area_2(traits, tf, rf, pf) << std::endl;
    // std::cout << "B: "  << area_2(traits, pf, qf, tf) << std::endl;
  }

  // Computes area of a 2D triangle.
  template<typename GeomTraits>
  typename GeomTraits::FT area_2(
    const GeomTraits& traits,
    const typename GeomTraits::Point_2& p,
    const typename GeomTraits::Point_2& q,
    const typename GeomTraits::Point_2& r) {

    const auto area_2 = traits.compute_area_2_object();
    return area_2(p, q, r);
  }

  // Computes positive area of a 2D triangle.
  template<typename GeomTraits>
  typename GeomTraits::FT positive_area_2(
    const GeomTraits& traits,
    const typename GeomTraits::Point_2& p,
    const typename GeomTraits::Point_2& q,
    const typename GeomTraits::Point_2& r) {

    return CGAL::abs(area_2(traits, p, q, r));
  }

  // Computes area of a 3D triangle.
  template<typename GeomTraits>
  typename GeomTraits::FT area_3(
    const GeomTraits& traits,
    const typename GeomTraits::Point_3& p,
    const typename GeomTraits::Point_3& q,
    const typename GeomTraits::Point_3& r) {

    using FT = typename GeomTraits::FT;
    using Point_3 = typename GeomTraits::Point_3;
    using Vector_3 = typename GeomTraits::Vector_3;

    const auto cross_product_3 =
      traits.construct_cross_product_vector_3_object();
    const auto construct_vector_3 =
      traits.construct_vector_3_object();
    const auto centroid_3 =
      traits.construct_centroid_3_object();

    // Compute centroid.
    const auto center = centroid_3(p, q, r);

    // Translate.
    const Point_3 a = Point_3(
      p.x() - center.x(), p.y() - center.y(), p.z() - center.z());
    const Point_3 b = Point_3(
      q.x() - center.x(), q.y() - center.y(), q.z() - center.z());
    const Point_3 c = Point_3(
      r.x() - center.x(), r.y() - center.y(), r.z() - center.z());

    // Prev and next vectors.
    auto v1 = construct_vector_3(b, a);
    auto v2 = construct_vector_3(b, c);
    normalize_3(traits, v1);
    normalize_3(traits, v2);

    // Compute normal.
    auto normal = cross_product_3(v1, v2);
    normalize_3(traits, normal);

    // Compute orthogonal base vectors.
    Vector_3 b1, b2;
    orthogonal_bases_3(traits, normal, b1, b2);

    // Compute area.
    const auto& origin = b;
    const auto pf = to_2d(traits, b1, b2, origin, a);
    const auto qf = to_2d(traits, b1, b2, origin, b);
    const auto rf = to_2d(traits, b1, b2, origin, c);

    const FT A = area_2(traits, pf, qf, rf);
    return A;
  }

  // Computes positive area of a 3D triangle.
  template<typename GeomTraits>
  typename GeomTraits::FT positive_area_3(
    const GeomTraits& traits,
    const typename GeomTraits::Point_3& p,
    const typename GeomTraits::Point_3& q,
    const typename GeomTraits::Point_3& r) {

    using FT = typename GeomTraits::FT;
    const auto construct_vector_3 =
      traits.construct_vector_3_object();

    const auto v1 = construct_vector_3(q, r);
    const auto v2 = construct_vector_3(q, p);

    const auto cross_product_3 =
      traits.construct_cross_product_vector_3_object();
    const auto cross = cross_product_3(v1, v2);
    const FT half = FT(1) / FT(2);
    const FT A = half * length_3(traits, cross);
    return A;
  }

  // Computes a clamped cotangent between two 3D vectors.
  // In the old version of weights in PMP, it has been called secure.
  // See Weights/internal/pmp_weights_deprecated.h for more information.
  template<typename GeomTraits>
  typename GeomTraits::FT cotangent_3_clamped(
    const GeomTraits& traits,
    const typename GeomTraits::Point_3& p,
    const typename GeomTraits::Point_3& q,
    const typename GeomTraits::Point_3& r) {

    using FT = typename GeomTraits::FT;
    using Get_sqrt = Get_sqrt<GeomTraits>;
    const auto sqrt = Get_sqrt::sqrt_object(traits);

    const auto dot_product_3 =
      traits.compute_scalar_product_3_object();
    const auto construct_vector_3 =
      traits.construct_vector_3_object();

    const auto v1 = construct_vector_3(q, r);
    const auto v2 = construct_vector_3(q, p);

    const FT dot = dot_product_3(v1, v2);
    const FT length_v1 = length_3(traits, v1);
    const FT length_v2 = length_3(traits, v2);

    const FT lb = -FT(999) / FT(1000), ub = FT(999) / FT(1000);
    FT cosine = dot / length_v1 / length_v2;
    cosine = (cosine < lb) ? lb : cosine;
    cosine = (cosine > ub) ? ub : cosine;
    const FT sine = sqrt(FT(1) - cosine * cosine);

    CGAL_assertion(sine != FT(0));
    if (sine != FT(0)) {
      return cosine / sine;
    }
    return FT(0); // undefined
  }

} // namespace internal
} // namespace Weights
} // namespace CGAL

#endif // CGAL_WEIGHTS_INTERNAL_UTILS_H
