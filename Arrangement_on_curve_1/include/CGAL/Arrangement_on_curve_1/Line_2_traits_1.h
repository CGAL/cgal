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

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_LINE_2_TRAITS_1_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_LINE_2_TRAITS_1_H

#include <CGAL/basic.h>
#include <CGAL/enum.h>

namespace CGAL {
namespace Arrangement_on_curve_1 {

template <typename Kernel_>
class Line_2_traits_1 {
public:
  using Kernel = Kernel_;
  using Point_1 = typename Kernel::Point_2;
  using Line_2 = typename Kernel::Line_2;
  using Vector_2 = typename Kernel::Vector_2;

  Line_2_traits_1(const Line_2& line) : m_line(line) {}

  class Compare_x_1 {
  private:
    Line_2 m_line;
  public:
    Compare_x_1(const Line_2& line) : m_line(line) {}

    Comparison_result operator()(const Point_1& p1, const Point_1& p2) const {
      if (p1 == p2) return EQUAL;

      Point_1 base = m_line.point();
      Vector_2 dir = m_line.to_vector();

      auto t1 = Vector_2(base, p1) * dir;
      auto t2 = Vector_2(base, p2) * dir;

      if (t1 < t2) return SMALLER;
      if (t1 > t2) return LARGER;
      return EQUAL;
    }
  };

  Compare_x_1 compare_x_1_object() const { return Compare_x_1(m_line); }

private:
  Line_2 m_line;
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
