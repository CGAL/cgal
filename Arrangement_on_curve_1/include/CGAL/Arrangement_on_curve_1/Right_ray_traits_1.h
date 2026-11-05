// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_RIGHT_RAY_TRAITS_1_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_RIGHT_RAY_TRAITS_1_H

#include <CGAL/Arrangement_on_curve_1/Ft_traits_1.h>
#include <CGAL/Arrangement_on_curve_1/Parameter_space.h>

namespace CGAL {
namespace Arrangement_on_curve_1 {

// Geometry traits for a 1D arrangement on a riht-closed ray \f$(-\infty, b]\f$
// The right side is closed; the left side is open.  Points are represented as
// `Kernel::FT`.  Pair with `Right_closed_topology_traits`.
template <typename Kernel_>
class Right_ray_traits_1 : public Ft_traits_1<Kernel_> {
  using Kernel = Kernel_;
  using Base = Ft_traits_1<Kernel>;

public:
  using Point_1 = typename Base::Point_1;

  using Left_side_category = Open_side_tag;
  using Right_side_category = Closed_side_tag;

  explicit Right_ray_traits_1(const Point_1& right_bound) : m_right_bound(right_bound) {}

  const Point_1& right_bound() const { return m_right_bound; }

  // `Parameter_space_in_x_1` obtains the location of a point in the parameter space.
  class Parameter_space_in_x_1 {
    const Right_ray_traits_1& m_tr;

  public:
    explicit Parameter_space_in_x_1(const Right_ray_traits_1& tr) : m_tr(tr) {}

    Parameter_space operator()(const Point_1& p) const
    { return (p == m_tr.m_right_bound) ? Parameter_space::RIGHT_BOUNDARY : Parameter_space::INTERIOR; }
  };

  Parameter_space_in_x_1 parameter_space_in_x_1_object() const { return Parameter_space_in_x_1(*this); }

private:
  Point_1 m_right_bound;
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
