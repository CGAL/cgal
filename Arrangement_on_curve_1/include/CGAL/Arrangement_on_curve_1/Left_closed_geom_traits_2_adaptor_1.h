// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_LEFT_CLOSED_GEOM_TRAITS_2_ADAPTOR_1_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_LEFT_CLOSED_GEOM_TRAITS_2_ADAPTOR_1_H

#include <CGAL/Arrangement_on_curve_1/Side_tags.h>
#include <CGAL/Arrangement_on_curve_1/Parameter_space.h>
#include <CGAL/Arrangement_on_curve_1/Geom_traits_2_adaptor_1_base.h>

namespace CGAL {
namespace Arrangement_on_curve_1 {

// Adapts a 2D geometry traits class into a 1D geometry traits class whose
// parameter space is closed on the left side and open on the right side.  The
// left boundary point is the lexicographically smallest point of the master
// curve.  Pair with `Left_bounded_topology_traits`.
template <typename GeomTraits_2>
class Left_closed_geom_traits_2_adaptor_1 : public internal::Geom_traits_2_adaptor_1_base<GeomTraits_2> {
  using Base = internal::Geom_traits_2_adaptor_1_base<GeomTraits_2>;

public:
  using typename Base::Geom_traits_2;
  using typename Base::Point_1;
  using typename Base::Shared_geom_traits_2;

  using Left_side_category = Closed_side_tag;
  using Right_side_category = Open_side_tag;

  /*! constructs from a shared 2D traits object and the left boundary point.
   */
  Left_closed_geom_traits_2_adaptor_1(Shared_geom_traits_2 traits_2, const Point_1& left_bound) :
    Base(std::move(traits_2)),
    m_left_bound(left_bound)
  {}

  /*! constructs from a copy of a 2D traits object and the left boundary point.
   */
  Left_closed_geom_traits_2_adaptor_1(const Geom_traits_2& traits_2, const Point_1& left_bound) :
    Left_closed_geom_traits_2_adaptor_1(std::make_shared<const Geom_traits_2>(traits_2), left_bound)
  {}

  /*! obtains the left boundary point.
   */
  const Point_1& left_bound() const { return m_left_bound; }

  /*! determines the location of a point in the parameter space:
   * `LEFT_BOUNDARY` if it is the left boundary point and `INTERIOR` otherwise.
   */
  class Parameter_space_in_x_1 {
  public:
    explicit Parameter_space_in_x_1(const Left_closed_geom_traits_2_adaptor_1& tr) : m_tr(tr) {}

    Parameter_space operator()(const Point_1& p) const {
      auto cmp_xy = m_tr.geom_traits_2().compare_xy_2_object();
      return (cmp_xy(p, m_tr.m_left_bound) == EQUAL) ? Parameter_space::LEFT_BOUNDARY : Parameter_space::INTERIOR;
    }

  private:
    const Left_closed_geom_traits_2_adaptor_1& m_tr;
  };

  /*! obtains a `Parameter_space_in_x_1` functor.
   */
  Parameter_space_in_x_1 parameter_space_in_x_1_object() const { return Parameter_space_in_x_1(*this); }

private:
  Point_1 m_left_bound;
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
