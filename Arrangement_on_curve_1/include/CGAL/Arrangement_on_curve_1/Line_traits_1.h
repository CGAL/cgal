// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_LINE_TRAITS_1_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_LINE_TRAITS_1_H

#include <CGAL/Arrangement_on_curve_1/Ft_traits_1.h>
#include <CGAL/Arrangement_on_curve_1/Parameter_space.h>

namespace CGAL {
namespace Arrangement_on_curve_1 {

/*! Geometry traits for a 1D arrangement on the full real line
 * \f$(-\infty,+\infty)\f$.  The parameter space is open on both sides.  Points
 * are represented as `Kernel::FT`. Pair with `Open_topology_traits`.
 */
template <typename Kernel_>
class Line_traits_1 : public Ft_traits_1<Kernel_> {
  using Kernel = Kernel_;
  using Base = Ft_traits_1<Kernel>;

public:
  using Point_1 = typename Base::Point_1;
  using Left_side_category = Open_side_tag;
  using Right_side_category = Open_side_tag;

  // Parameter_space_in_x_1---every finite point is in the interior.
  class Parameter_space_in_x_1 {
    Parameter_space operator()(const Point_1&) const { return Parameter_space::INTERIOR; }
  };

  Parameter_space_in_x_1 parameter_space_in_x_1_object() const { return Parameter_space_in_x_1{}; }
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
