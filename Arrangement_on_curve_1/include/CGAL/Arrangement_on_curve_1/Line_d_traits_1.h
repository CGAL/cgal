// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_LINE_D_TRAITS_1_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_LINE_D_TRAITS_1_H

#include <vector>

#include <CGAL/enum.h>

namespace CGAL {
namespace Arrangement_on_curve_1 {

template <typename Kernel_>
class Line_d_traits_1 {
public:
  using Kernel = Kernel_;
  using Point_1 = typename Kernel::Point_d;
  using Line_d = typename Kernel::Line_d;
  using Vector_d = typename Kernel::Vector_d;
  using FT = typename Kernel::FT;

  Line_d_traits_1(const Line_d& line) : m_line(line) {}

  class Compare_x_1 {
  private:
    Line_d m_line;

    // Helper function to calculate a dot product for arbitrary d-dimensional vectors
    FT dot_product(const Vector_d& v1, const Vector_d& v2) const {
      FT result = 0;
      // Accessing coordinate dimensions via Cartesian iterator representation
      auto it1 = v1.cartesian_begin();
      auto it2 = v2.cartesian_begin();
      while (it1 != v1.cartesian_end() && it2 != v2.cartesian_end()) {
        result += (*it1) * (*it2);
        ++it1;
        ++it2;
      }
      return result;
    }

  public:
    Compare_x_1(const Line_d& line) : m_line(line) {}

    Comparison_result operator()(const Point_1& p1, const Point_1& p2) const {
      if (p1 == p2) return EQUAL;

      // Parameterize points by projecting onto the d-dimensional directional vector
      Point_1 base = m_line.point(0);
      Vector_d dir = m_line.direction().vector();

      Vector_d v1 = p1 - base;
      Vector_d v2 = p2 - base;

      FT t1 = dot_product(v1, dir);
      FT t2 = dot_product(v2, dir);

      if (t1 < t2) return SMALLER;
      if (t1 > t2) return LARGER;
      return EQUAL;
    }
  };

  Compare_x_1 compare_x_1_object() const { return Compare_x_1(m_line); }

private:
  Line_d m_line;
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
