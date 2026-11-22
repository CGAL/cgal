// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_FT_TRAITS_1_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_FT_TRAITS_1_H

#include <CGAL/enum.h>
#include <CGAL/Arrangement_on_curve_1/Side_tags.h>

namespace CGAL {
namespace Arrangement_on_curve_1 {

// Base class for FT-scalar geometry traits.
//
// Factors out the types and predicates that are common to all FT-based
// geometry traits: Point_1 = Kernel::FT, and Compare_x_1.  Derived classes
// add Left_side_category / Right_side_category and Parameter_space_in_x_.
//
// Not intended to be instantiated directly; use `Line_traits_1`, `Segment_traits_1`,
// `Left_ray_traits_1`, or `Right_ray_traits_1`.
template <typename Kernel_>
class Ft_traits_1 {
public:
  using Kernel = Kernel_;
  using Point_1 = typename Kernel::FT;

  // Compare_x_1---total order on FT scalars.
  class Compare_x_1 {
  public:
    Comparison_result operator()(const Point_1& p1, const Point_1& p2) const {
      if (p1 < p2) return SMALLER;
      if (p1 > p2) return LARGER;
      return EQUAL;
    }
  };

  Compare_x_1 compare_x_1_object() const { return Compare_x_1(); }
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
