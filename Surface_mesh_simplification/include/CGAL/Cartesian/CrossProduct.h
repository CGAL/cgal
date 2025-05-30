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
#ifndef CGAL_CARTESIAN_CROSSPRODUCT_H
#define CGAL_CARTESIAN_CROSSPRODUCT_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <optional>

namespace CGAL {

// a*b - c*d
  // The next two functions are from https://stackoverflow.com/questions/63665010/accurate-floating-point-computation-of-the-sum-and-difference-of-two-products
  static double diff_of_products_kahan(const double a, const double b, const double c, const double d)
  {
    double w = d * c;
    double e = std::fma(c, -d, w);
    double f = std::fma(a, b, -w);
    return f + e;
  }

  static double diff_of_products_cht(const double a, const double b, const double c, const double d)
  {
    double p1 = a * b;
    double p2 = c * d;
    double e1 = std::fma (a, b, -p1);
    double e2 = std::fma (c, -d, p2);
    double r = p1 - p2;
    double e = e1 + e2;
    return r + e;
  }

  static double diff_of_products(const double a, const double b, const double c, const double d)
  {
    // return a*b - c*d;
    // the next two are equivalent in results and speed
    return diff_of_products_kahan(a, b, c, d);
    // return diff_of_products_cht(a, b, c, d);
  }

  template <typename OFT>
  static OFT diff_of_products(const OFT& a, const OFT& b, const OFT& c, const OFT& d)
  {
    return a*b - c*d;
  }


} // namespace CGAL

#endif // CGAL_CARTESIAN_CROSSPRODUCT_H //
// EOF //


