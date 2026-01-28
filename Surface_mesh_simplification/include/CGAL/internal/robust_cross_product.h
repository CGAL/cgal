// Copyright (c) 2025  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//  Author(s)     : Mael Rouxel-Labbé
//
#ifndef CGAL_SMS_INTERNAL_ROBUST_CROSS_PRODUCT_H
#define CGAL_SMS_INTERNAL_ROBUST_CROSS_PRODUCT_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/Surface_mesh_simplification/internal/Common.h>

#include <optional>

namespace CGAL {
namespace Surface_mesh_simplification {
namespace internal {

// a*b - c*d
// The next two functions are from https://stackoverflow.com/questions/63665010/accurate-floating-point-computation-of-the-sum-and-difference-of-two-products
inline double diff_of_products_kahan(const double a, const double b, const double c, const double d)
{
  double w = d * c;
  double e = std::fma(c, -d, w);
  double f = std::fma(a, b, -w);
  return f + e;
}

inline double diff_of_products_cht(const double a, const double b, const double c, const double d)
{
  double p1 = a * b;
  double p2 = c * d;
  double e1 = std::fma (a, b, -p1);
  double e2 = std::fma (c, -d, p2);
  double r = p1 - p2;
  double e = e1 + e2;
  return r + e;
}

inline double diff_of_products(const double a, const double b, const double c, const double d)
{
#if 0
  // this can create large errors with inexact constructions
  return a*b - c*d;
  // the next two are equivalent in results and speed
#elif 1
  return diff_of_products_kahan(a, b, c, d);
#elif 0
  return diff_of_products_cht(a, b, c, d);
#endif
}

template <typename OFT>
inline OFT diff_of_products(const OFT& a, const OFT& b, const OFT& c, const OFT& d)
{
  return a*b - c*d;
}

// balanced solution based on abusing the fact that here we expect u and v to have similar coordinates
template<class GeomTraits>
typename GeomTraits::Vector_3 similar_coordinates_cross_product(const typename GeomTraits::Vector_3& u,
                                                                const typename GeomTraits::Vector_3& v)
{
  using FT = typename GeomTraits::FT;
  using Vector = typename GeomTraits::Vector_3;

  const FT& ux = u.x();
  const FT& uy = u.y();
  const FT& uz = u.z();
  const FT& vx = v.x();
  const FT& vy = v.y();
  const FT& vz = v.z();

  auto minor = [](const FT& ui, const FT& vi, const FT& uj, const FT& vj)
  {
    // The main idea is that we expect ai and bi (and aj and bj) to have roughly the same magnitude
    // since this function is used to compute the cross product of two vectors that are defined
    // as (ORIGIN, pa) and (ORIGIN, pb) and pa and pb are part of the same triangle.
    //
    // We can abuse this fact to trade 2 extra subtractions to lower the error.
    return ui * (vj - uj) + uj * (ui - vi);
  };

  // ay*
  FT x = minor(uy, vy, uz, vz);
  FT y = minor(uz, vz, ux, vx);
  FT z = minor(ux, vx, uy, vy);

  return Vector(x, y, z);
}

#if 0
template<class GeomTraits>
typename GeomTraits::Vector_3 exact_cross_product(const typename GeomTraits::Vector_3& a,
                                                  const typename GeomTraits::Vector_3& b)
{
  CGAL::Cartesian_converter<GeomTraits, CGAL::Exact_predicates_exact_constructions_kernel> to_exact;
  CGAL::Cartesian_converter<CGAL::Exact_predicates_exact_constructions_kernel, GeomTraits> to_approx;
  auto exv = cross_product(to_exact(a), to_exact(b));
  exv.exact();
  return to_approx(exv);
}
#endif

template<class GeomTraits>
typename GeomTraits::Vector_3 robust_cross_product(const typename GeomTraits::Vector_3& u,
                                                   const typename GeomTraits::Vector_3& v)
{
#if 0
  // this can create large errors and spiky meshes for kernels with inexact constructions
  return CGAL::cross_product(u,v);
#elif 0
  // improves the problem mentioned above a bit, but not enough
  return { std::fma(u.y(), v.z(), -u.z()*v.y()),
           std::fma(u.z(), v.x(), -u.x()*v.z()),
           std::fma(u.x(), v.y(), -u.y()*v.x()) };
#elif 0
  // this is the best without resorting to exact, but it inflicts a 20% slowdown
  return { diff_of_products(u.y(), v.z(), u.z(), v.y()),
           diff_of_products(u.z(), v.x(), u.x(), v.z()),
           diff_of_products(u.x(), v.y(), u.y(), v.x()) };
#elif 0
  // obviously too slow
  return exact_cross_product(u, v);
#elif 1
  // balanced solution based on abusing the fact that in this package, we usually have that
  // u and v to have similar coordinates
  return similar_coordinates_cross_product<GeomTraits>(u, v);
#endif
}

} // namespace CGAL
} // namespace Surface_mesh_simplification
} // namespace internal

#endif // CGAL_SMS_INTERNAL_ROBUST_CROSS_PRODUCT_H


