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

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/number_utils.h>
#include <CGAL/property_map.h>
#include <CGAL/utils.h>

#include <boost/mpl/has_xxx.hpp>

#include <cmath>
#include <memory>
#include <iostream>
#include <iterator>
#include <vector>

namespace CGAL {
namespace Weights {
namespace internal {

// Sqrt helpers.
template<typename GeomTraits>
class Default_sqrt
{
private:
  using Traits = GeomTraits;
  using FT = typename Traits::FT;

public:
  FT operator()(const FT value) const
  {
    return static_cast<FT>(CGAL::sqrt(CGAL::to_double(CGAL::abs(value))));
  }
};

BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_nested_type_Sqrt, Sqrt, false)

// Case: do_not_use_default = false.
template<typename GeomTraits,
         bool do_not_use_default = Has_nested_type_Sqrt<GeomTraits>::value>
class Get_sqrt
{
public:
  using Traits = GeomTraits;
  using Sqrt = Default_sqrt<Traits>;

  static Sqrt sqrt_object(const Traits&)
  {
    return Sqrt();
  }
};

// Case: do_not_use_default = true.
template<typename GeomTraits>
class Get_sqrt<GeomTraits, true>
{
public:
  using Traits = GeomTraits;
  using Sqrt = typename Traits::Sqrt;

  static Sqrt sqrt_object(const Traits& traits)
  {
    return traits.sqrt_object();
  }
};

template<typename FT>
void normalize(std::vector<FT>& values)
{
  FT sum = FT(0);
  for (const FT& value : values)
    sum += value;

  CGAL_assertion(!is_zero(sum));
  if (is_zero(sum))
    return;

  const FT inv_sum = FT(1) / sum;
  for (FT& value : values)
    value *= inv_sum;
}

template<typename FT>
FT power(const FT value,
         const FT p)
{
  const double base = CGAL::to_double(value);
  const double exp = CGAL::to_double(p);

  return static_cast<FT>(std::pow(base, exp));
}

// 2D ==============================================================================================

template<typename GeomTraits>
typename GeomTraits::FT distance_2(const typename GeomTraits::Point_2& p,
                                   const typename GeomTraits::Point_2& q,
                                   const GeomTraits& traits)
{
  using Get_sqrt = Get_sqrt<GeomTraits>;
  auto sqrt = Get_sqrt::sqrt_object(traits);

  auto squared_distance_2 = traits.compute_squared_distance_2_object();
  return sqrt(squared_distance_2(p, q));
}

template<typename GeomTraits>
typename GeomTraits::FT length_2(const typename GeomTraits::Vector_2& v,
                                 const GeomTraits& traits)
{
  using Get_sqrt = Get_sqrt<GeomTraits>;
  auto sqrt = Get_sqrt::sqrt_object(traits);

  auto squared_length_2 = traits.compute_squared_length_2_object();
  return sqrt(squared_length_2(v));
}

template<typename GeomTraits>
void normalize_2(typename GeomTraits::Vector_2& v,
                 const GeomTraits& traits)
{
  using FT = typename GeomTraits::FT;

  const FT length = length_2(v, traits);
  CGAL_assertion(!is_zero(length));
  if (is_zero(length))
    return;

  v /= length;
}

// 3D ==============================================================================================

template<typename GeomTraits>
typename GeomTraits::FT distance_3(const typename GeomTraits::Point_3& p,
                                   const typename GeomTraits::Point_3& q,
                                   const GeomTraits& traits)
{
  auto squared_distance_3 = traits.compute_squared_distance_3_object();

  using Get_sqrt = Get_sqrt<GeomTraits>;
  auto sqrt = Get_sqrt::sqrt_object(traits);

  return sqrt(squared_distance_3(p, q));
}

template<typename GeomTraits>
typename GeomTraits::FT length_3(const typename GeomTraits::Vector_3& v,
                                 const GeomTraits& traits)
{
  using Get_sqrt = Get_sqrt<GeomTraits>;
  auto sqrt = Get_sqrt::sqrt_object(traits);

  auto squared_length_3 = traits.compute_squared_length_3_object();
  return sqrt(squared_length_3(v));
}

template<typename GeomTraits>
void normalize_3(typename GeomTraits::Vector_3& v,
                 const GeomTraits& traits)
{
  using FT = typename GeomTraits::FT;
  const FT length = length_3(v, traits);
  CGAL_assertion(!is_zero(length));
  if (is_zero(length))
    return;

  v /= length;
}

// the angle is in radians
template<typename GeomTraits>
double angle_3(const typename GeomTraits::Vector_3& v1,
               const typename GeomTraits::Vector_3& v2,
               const GeomTraits& traits)
{
  auto dot_product_3 = traits.compute_scalar_product_3_object();

  const double product = CGAL::sqrt(to_double(dot_product_3(v1,v1) * dot_product_3(v2,v2)));
  if(product == 0.)
    return 0.;

  const double dot = CGAL::to_double(dot_product_3(v1, v2));
  const double cosine = dot / product;

  if (cosine < -1.0)
    return std::acos(-1.0);
  else if (cosine > 1.0)
    return std::acos(+1.0);
  else
    return std::acos(cosine);
}

// Rotates a 3D point around axis.
template<typename GeomTraits>
typename GeomTraits::Point_3 rotate_point_3(const double angle_rad,
                                            const typename GeomTraits::Vector_3& axis,
                                            const typename GeomTraits::Point_3& query,
                                            const GeomTraits& traits)
{
  using FT = typename GeomTraits::FT;

  auto point_3 = traits.construct_point_3_object();

  const FT c = static_cast<FT>(std::cos(angle_rad));
  const FT s = static_cast<FT>(std::sin(angle_rad));
  const FT C = FT(1) - c;

  const FT& x = axis.x();
  const FT& y = axis.y();
  const FT& z = axis.z();

  return point_3(
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
void orthogonal_bases_3(const typename GeomTraits::Vector_3& normal,
                        typename GeomTraits::Vector_3& b1,
                        typename GeomTraits::Vector_3& b2,
                        const GeomTraits& traits)
{
  using FT = typename GeomTraits::FT;
  using Vector_3 = typename GeomTraits::Vector_3;

  auto cross_product_3 = traits.construct_cross_product_vector_3_object();

  const FT& nx = normal.x();
  const FT& ny = normal.y();
  const FT& nz = normal.z();

  if (CGAL::abs(nz) >= CGAL::abs(ny))
    b1 = Vector_3(nz, 0, -nx);
  else
    b1 = Vector_3(ny, -nx, 0);

  b2 = cross_product_3(normal, b1);

  normalize_3(b1, traits);
  normalize_3(b2, traits);
}

// Converts a 3D point into a 2D point wrt to a given plane.
template<typename GeomTraits>
typename GeomTraits::Point_2 to_2d(const typename GeomTraits::Vector_3& b1,
                                   const typename GeomTraits::Vector_3& b2,
                                   const typename GeomTraits::Point_3& origin,
                                   const typename GeomTraits::Point_3& query,
                                   const GeomTraits& traits)
{
  using FT = typename GeomTraits::FT;
  using Vector_3 = typename GeomTraits::Vector_3;

  auto dot_product_3 = traits.compute_scalar_product_3_object();
  auto vector_3 = traits.construct_vector_3_object();
  auto point_2 = traits.construct_point_2_object();

  const Vector_3 v = vector_3(origin, query);
  const FT x = dot_product_3(b1, v);
  const FT y = dot_product_3(b2, v);

  return point_2(x, y);
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
void flatten(const typename GeomTraits::Point_3& t, // prev neighbor/vertex/point
             const typename GeomTraits::Point_3& r, // curr neighbor/vertex/point
             const typename GeomTraits::Point_3& p, // next neighbor/vertex/point
             const typename GeomTraits::Point_3& q, // query point
             typename GeomTraits::Point_2& tf,
             typename GeomTraits::Point_2& rf,
             typename GeomTraits::Point_2& pf,
             typename GeomTraits::Point_2& qf,
             const GeomTraits& traits)
{
  // std::cout << std::endl;
  using Point_3 = typename GeomTraits::Point_3;
  using Vector_3 = typename GeomTraits::Vector_3;

  auto point_3 = traits.construct_point_3_object();
  auto cross_product_3 = traits.construct_cross_product_vector_3_object();
  auto vector_3 = traits.construct_vector_3_object();
  auto centroid_3 = traits.construct_centroid_3_object();

  // Compute centroid.
  const Point_3 center = centroid_3(t, r, p, q);
  // std::cout << "centroid: " << center << std::endl;

  // Translate.
  const Point_3 t1 = point_3(t.x() - center.x(), t.y() - center.y(), t.z() - center.z());
  const Point_3 r1 = point_3(r.x() - center.x(), r.y() - center.y(), r.z() - center.z());
  const Point_3 p1 = point_3(p.x() - center.x(), p.y() - center.y(), p.z() - center.z());
  const Point_3 q1 = point_3(q.x() - center.x(), q.y() - center.y(), q.z() - center.z());

  // std::cout << "translated t1: " << t1 << std::endl;
  // std::cout << "translated r1: " << r1 << std::endl;
  // std::cout << "translated p1: " << p1 << std::endl;
  // std::cout << "translated q1: " << q1 << std::endl;

  // Middle axis.
  Vector_3 ax = vector_3(q1, r1);
  normalize_3(ax, traits);

  // Prev and next vectors.
  Vector_3 v1 = vector_3(q1, t1);
  Vector_3 v2 = vector_3(q1, p1);

  // Two triangle normals.
  Vector_3 n1 = cross_product_3(v1, ax);
  Vector_3 n2 = cross_product_3(ax, v2);

  // std::cout << "normal n1: " << n1 << std::endl;
  // std::cout << "normal n2: " << n2 << std::endl;

  // Angle between two normals.
  const double angle_rad = angle_3(n1, n2, traits);
  // std::cout << "angle deg n1 <-> n2: " << angle_rad * 180.0 / CGAL_PI << std::endl;

  // Rotate p1 around ax so that it lands onto the plane [q1, t1, r1].
  const Point_3& t2 = t1;
  const Point_3& r2 = r1;
  const Point_3  p2 = rotate_point_3(angle_rad, ax, p1, traits);
  const Point_3& q2 = q1;
  // std::cout << "rotated p2: " << p2 << std::endl;

  // Compute orthogonal base vectors.
  Vector_3 b1, b2;
  const Vector_3& normal = n1;
  orthogonal_bases_3(normal, b1, b2, traits);

  // const Angle angle12 = angle_3(b1, b2, traits);
  // std::cout << "angle deg b1 <-> b2: " << angle12 * 180.0 / CGAL_PI << std::endl;

  // Flatten a quad.
  const Point_3& origin = q2;
  tf = to_2d(b1, b2, origin, t2, traits);
  rf = to_2d(b1, b2, origin, r2, traits);
  pf = to_2d(b1, b2, origin, p2, traits);
  qf = to_2d(b1, b2, origin, q2, traits);

  // std::cout << "flattened qf: " << qf << std::endl;
  // std::cout << "flattened tf: " << tf << std::endl;
  // std::cout << "flattened rf: " << rf << std::endl;
  // std::cout << "flattened pf: " << pf << std::endl;

  // std::cout << "A1: " << area_2(rf, qf, pf, traits) << std::endl;
  // std::cout << "A2: " << area_2(pf, qf, rf, traits) << std::endl;
  // std::cout << "C: "  << area_2(tf, rf, pf, traits) << std::endl;
  // std::cout << "B: "  << area_2(pf, qf, tf, traits) << std::endl;
}

template<typename GeomTraits>
typename GeomTraits::FT positive_area_2(const typename GeomTraits::Point_2& p,
                                        const typename GeomTraits::Point_2& q,
                                        const typename GeomTraits::Point_2& r,
                                        const GeomTraits& traits)
{
  auto area_2 = traits.compute_area_2_object();
  return CGAL::abs(area_2(p, q, r));
}

template<typename GeomTraits>
typename GeomTraits::FT area_3(const typename GeomTraits::Point_3& p,
                               const typename GeomTraits::Point_3& q,
                               const typename GeomTraits::Point_3& r,
                               const GeomTraits& traits)
{
  using Point_2 = typename GeomTraits::Point_2;
  using Point_3 = typename GeomTraits::Point_3;
  using Vector_3 = typename GeomTraits::Vector_3;

  auto area_2 = traits.compute_area_2_object();
  auto point_3 = traits.construct_point_3_object();
  auto cross_product_3 = traits.construct_cross_product_vector_3_object();
  auto vector_3 = traits.construct_vector_3_object();
  auto centroid_3 = traits.construct_centroid_3_object();

  // Compute centroid.
  const Point_3 center = centroid_3(p, q, r);

  // Translate.
  const Point_3 a = point_3(p.x() - center.x(), p.y() - center.y(), p.z() - center.z());
  const Point_3 b = point_3(q.x() - center.x(), q.y() - center.y(), q.z() - center.z());
  const Point_3 c = point_3(r.x() - center.x(), r.y() - center.y(), r.z() - center.z());

  // Prev and next vectors.
  Vector_3 v1 = vector_3(b, a);
  Vector_3 v2 = vector_3(b, c);

  // Compute normal.
  Vector_3 normal = cross_product_3(v1, v2);

  // Compute orthogonal base vectors.
  Vector_3 b1, b2;
  orthogonal_bases_3(normal, b1, b2, traits);

  // Compute area.
  const Point_3& origin = b;
  const Point_2 pf = to_2d(b1, b2, origin, a, traits);
  const Point_2 qf = to_2d(b1, b2, origin, b, traits);
  const Point_2 rf = to_2d(b1, b2, origin, c, traits);

  return area_2(pf, qf, rf);
}

// Computes positive area of a 3D triangle.
template<typename GeomTraits>
typename GeomTraits::FT positive_area_3(const typename GeomTraits::Point_3& p,
                                        const typename GeomTraits::Point_3& q,
                                        const typename GeomTraits::Point_3& r,
                                        const GeomTraits& traits)
{
  using Get_sqrt = Get_sqrt<GeomTraits>;
  auto sqrt = Get_sqrt::sqrt_object(traits);

  auto squared_area_3 = traits.compute_squared_area_3_object();

  return sqrt(squared_area_3(p, q, r));
}

} // namespace internal
} // namespace Weights
} // namespace CGAL

#endif // CGAL_WEIGHTS_INTERNAL_UTILS_H
