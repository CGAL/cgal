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

namespace CGAL::Surface_mesh_simplification::internal{

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
  template<class Geom_traits>
  typename Geom_traits::Vector_3 robust_cross_product(const typename Geom_traits::Vector_3& a, const typename Geom_traits::Vector_3& b)
  {
    using FT = typename Geom_traits::FT;
    using Vector = typename Geom_traits::Vector_3;

    const FT& ax=a.x();
    const FT& ay=a.y();
    const FT& az=a.z();
    const FT& bx=b.x();
    const FT& by=b.y();
    const FT& bz=b.z();

    auto minor = [](const FT& ai, const FT& bi, const FT& aj, const FT& bj)
    {
      // The main idea is that we expect ai and bi (and aj and bj) to have roughly the same magnitude
      // since this function is used to compute the cross product of two vectors that are defined
      // as (ORIGIN, pa) and (ORIGIN, pb) and pa and pb are part of the same triangle.
      //
      // We can abuse this fact to trade 2 extra subtractions to lower the error.
      return ai * (bj - aj) + aj * (ai - bi);
    };

    // ay*
    FT x = minor(ay, by, az, bz);
    FT y = minor(az, bz, ax, bx);
    FT z = minor(ax, bx, ay, by);

    return Vector(x, y, z);
  }

#if 0
  template<class Geom_traits>
  typename Geom_traits::Vector_3 exact_cross_product(const typename Geom_traits::Vector_3& a, const typename Geom_traits::Vector_3& b)
  {
    CGAL::Cartesian_converter<Geom_traits, CGAL::Exact_predicates_exact_constructions_kernel> to_exact;
    CGAL::Cartesian_converter<CGAL::Exact_predicates_exact_constructions_kernel, Geom_traits> to_approx;
    auto exv = cross_product(to_exact(a), to_exact(b));
    exv.exact();
    return to_approx(exv);
  }



  template<class Geom_traits>
  typename Geom_traits::Vector_3 experimental_cross_product(const typename Geom_traits::Vector_3& a, const typename Geom_traits::Vector_3& b)
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
    // balanced solution based on abusing the fact that here we expect u and v to have similar coordinates
    return robust_cross_product(u, v);
#endif
  }

#endif


} // namespace CGAL

#endif // CGAL_SMS_INTERNAL_ROBUST_CROSS_PRODUCT_H


