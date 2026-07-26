// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_FT_TRAITS_1_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_FT_TRAITS_1_H

namespace CGAL {
namespace Arrangement_on_curve_1 {

// Minimal geometry traits model utilizing direct coordinate FT parameters
template <typename Kernel>
class Ft_traits_1 {
public:
  using Point_1 = typename Kernel::FT;

  class Compare_x_1 {
  public:
    CGAL::Comparison_result operator()(const Point_1& p1, const Point_1& p2) const {
      if (p1 < p2) return CGAL::SMALLER;
      if (p1 > p2) return CGAL::LARGER;
      return CGAL::EQUAL;
    }
  };

  Compare_x_1 compare_x_1_object() const { return Compare_x_1(); }
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
