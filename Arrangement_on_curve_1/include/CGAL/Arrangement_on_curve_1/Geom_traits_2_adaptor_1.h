// Copyright (c) 2026 Tel-Aviv University (Israel). All rights reserved.
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_GEOM_TRAITS_2_ADAPTOR_1_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_GEOM_TRAITS_2_ADAPTOR_1_H

#include <memory>

#include <CGAL/enum.h>

namespace CGAL {
namespace Arrangement_on_curve_1 {

// Adapts any 2D geometry traits class of the 2D Arrangements package
// (i.e., any model of the ArrangementTraits_2 concept, or a refinement
// thereof) into a 1D geometry traits class, suitable for use with
// Arrangement_on_curve_1. The points of the 1D arrangement are the
// Point_2 objects of the 2D traits, all assumed to lie on a single
// curve -- the "master arc" -- which is either x-monotone or a single
// vertical segment.
template <typename GeomTraits_2>
class Geom_traits_2_adaptor_1 {
public:
  using Geom_traits_2 = GeomTraits_2;
  using Point_1 = typename Geom_traits_2::Point_2;
  using Curve_2 = typename Geom_traits_2::Curve_2;
  using X_monotone_curve_2 = typename Geom_traits_2::X_monotone_curve_2;

  using Geom_traits_2_ptr = std::shared_ptr<const Geom_traits_2>;

  // Constructor accepts a shared pointer to the 2D geometry traits
  // and the specific x-monotone curve (or vertical segment) along which
  // all 1D points are assumed to lie.
  Geom_traits_2_adaptor_1(Geom_traits_2_ptr traits_2_ptr, const X_monotone_curve_2& master_arc)
    : m_traits_2(traits_2_ptr), m_master_arc(master_arc) {}

  class Compare_x_1 {
  private:
    Geom_traits_2_ptr m_traits_2;
    X_monotone_curve_2 m_arc;
  public:
    Compare_x_1(Geom_traits_2_ptr t2, const X_monotone_curve_2& arc) : m_traits_2(t2), m_arc(arc) {}

    Comparison_result operator()(const Point_1& p1, const Point_1& p2) const {
      // Lexicographic (x, then y) comparison.
      auto cmp_xy_obj = m_traits_2->compare_xy_2_object();
      auto cmp_xy = cmp_xy_obj(p1, p2);
      if (cmp_xy == EQUAL) return EQUAL;

      // The master arc is either x-monotone or a single vertical segment.
      // If it is x-monotone, the points are ordered by their x-coordinates.
      auto cmp_x_obj = m_traits_2->compare_x_2_object();
      auto cmp_x = cmp_x_obj(p1, p2);
      if (cmp_x != EQUAL) return cmp_x;

      // The two points share the same x-coordinate, i.e., the master arc
      // is a vertical segment. In this case, fall back to the (non-equal)
      // result of the lexicographic comparison, which orders the points
      // by their y-coordinates.
      return cmp_xy;
    }
  };

  Compare_x_1 compare_x_1_object() const { return Compare_x_1(m_traits_2, m_master_arc); }

private:
  Geom_traits_2_ptr m_traits_2;
  X_monotone_curve_2 m_master_arc;
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
