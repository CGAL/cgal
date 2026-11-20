// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_LINE_D_TRAITS_1_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_LINE_D_TRAITS_1_H

#include <vector>

#include <CGAL/enum.h>
#include <CGAL/Arrangement_on_curve_1/Side_tags.h>
#include <CGAL/Arrangement_on_curve_1/Parameter_space.h>

namespace CGAL {
namespace Arrangement_on_curve_1 {

// Geometry traits for a 1D arrangement of points on a line in
// \f$\mathbb{R}^d\f$. Points are ordered by their signed arc-length parameter
// along the master line.  The line extends to \f$\pm\infty\f$, so both sides
// are open.  Pair with `Open_topology_traits`.
template <typename Kernel_>
class Line_d_traits_1 {
public:
  using Kernel = Kernel_;
  using Point_1 = typename Kernel::Point_d;
  using Line_d = typename Kernel::Line_d;
  using Vector_d = typename Kernel::Vector_d;
  using FT = typename Kernel::FT;

  // Boundary-category tags (geometry traits owns these).  The line extends to
  // \f$\pm\infty\f$; no boundary vertex is stored at either end.
  using Left_side_category = Open_side_tag;
  using Right_side_category = Open_side_tag;

  Line_d_traits_1(const Line_d& line) : m_line(line) {}

  class Compare_x_1 {
  private:
    Line_d m_line;

    // Dot product for arbitrary d-dimensional vectors.
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

      // Parameterize points by projecting onto the d-dimensional direction.
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

  // Both sides are open; every finite point is in the interior.
  struct Parameter_space_in_x_1 {
    Parameter_space operator()(const Point_1&) const { return Parameter_space::INTERIOR; }
  };

  Parameter_space_in_x_1 parameter_space_in_x_1_object() const { return Parameter_space_in_x_1{}; }

private:
  Line_d m_line;
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
