// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_GEOM_TRAITS_2_ADAPTOR_1_BASE_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_GEOM_TRAITS_2_ADAPTOR_1_BASE_H

#include <memory>

#include <CGAL/enum.h>
#include <CGAL/assertions.h>

namespace CGAL {
namespace Arrangement_on_curve_1 {
namespace internal {

// Base of all adaptors of a 2D geometry traits class to a 1D geometry traits
// class.  It holds the shared 2D traits and provides `Compare_x_1`.
//
// The points of the 1D arrangement are `Point_2` objects of the 2D traits that
// lie on a single master curve, which is either `x`-monotone or vertical.
// Along such a curve the 1D order is the lexicographic (`x`, then `y`) order,
// so `Compare_x_1` is `Compare_xy_2`.
template <typename GeomTraits_2>
class Geom_traits_2_adaptor_1_base {
public:
  using Geom_traits_2 = GeomTraits_2;
  using Point_1 = typename Geom_traits_2::Point_2;
  using Shared_geom_traits_2 = std::shared_ptr<const Geom_traits_2>;

  /*! obtains the underlying 2D traits object.
   */
  const Geom_traits_2& geom_traits_2() const { return *m_traits_2; }

  /*! obtains the shared pointer to the underlying 2D traits object.
   */
  const Shared_geom_traits_2& shared_geom_traits_2() const { return m_traits_2; }

  /*! compares two points along the master curve lexicographically.
   */
  class Compare_x_1 {
  public:
    explicit Compare_x_1(const Geom_traits_2& traits_2) : m_traits_2(traits_2) {}

    Comparison_result operator()(const Point_1& p1, const Point_1& p2) const
    { return m_traits_2.compare_xy_2_object()(p1, p2); }

  private:
    const Geom_traits_2& m_traits_2;
  };

  /*! obtains a `Compare_x_1` functor.
   */
  Compare_x_1 compare_x_1_object() const { return Compare_x_1(*m_traits_2); }

protected:
  explicit Geom_traits_2_adaptor_1_base(Shared_geom_traits_2 traits_2) : m_traits_2(std::move(traits_2))
  { CGAL_precondition(m_traits_2 != nullptr); }

  Shared_geom_traits_2 m_traits_2;
};

} // namespace internal
} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
