// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_DIRECTION_GEOM_TRAITS_2_ADAPTOR_1_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_DIRECTION_GEOM_TRAITS_2_ADAPTOR_1_H

#include <memory>
#include <utility>

#include <CGAL/assertions.h>
#include <CGAL/enum.h>
#include <CGAL/Arrangement_on_curve_1/Side_tags.h>

namespace CGAL {
namespace Arrangement_on_curve_1 {

// Adapts a 2D geometry traits class that provides directions in the plane and
// the predicates `Counterclockwise_in_between_2` and `Equal_2` (e.g., any
// kernel) into a 1D geometry traits class of the circle of directions, whose
// parameter space is identified.  The identification point is a given
// reference direction s, and the directions are ordered counterclockwise,
// starting right after s.
//
// The order is determined by `Counterclockwise_in_between_2`, which refers to
// no axis: a direction a precedes a direction b if and only if a is
// counterclockwise strictly in between s and b.  The equality of directions,
// and in particular the test whether a direction is the identification point,
// is determined by `Equal_2`.  Pair with `Identified_topology_traits`.
template <typename GeomTraits_2>
class Direction_geom_traits_2_adaptor_1 {
public:
  using Geom_traits_2 = GeomTraits_2;
  using Point_1 = typename Geom_traits_2::Direction_2;
  using Shared_geom_traits_2 = std::shared_ptr<const Geom_traits_2>;

  using Left_side_category = Identified_side_tag;
  using Right_side_category = Identified_side_tag;

  /*! constructs from a shared 2D traits object and a reference direction,
   * which is the identification point.
   */
  Direction_geom_traits_2_adaptor_1(Shared_geom_traits_2 traits_2, const Point_1& reference) :
    m_traits_2(std::move(traits_2)),
    m_reference(reference)
  { CGAL_precondition(m_traits_2 != nullptr); }

  /*! constructs from a copy of a 2D traits object and a reference direction,
   * which is the identification point.
   */
  Direction_geom_traits_2_adaptor_1(const Geom_traits_2& traits_2, const Point_1& reference) :
    m_traits_2(std::make_shared<const Geom_traits_2>(traits_2)),
    m_reference(reference)
  {}

  /*! obtains the underlying 2D traits object.
   */
  const Geom_traits_2& geom_traits_2() const { return *m_traits_2; }

  /*! obtains the shared pointer to the underlying 2D traits object.
   */
  const Shared_geom_traits_2& shared_geom_traits_2() const { return m_traits_2; }

  /*! obtains the reference direction, which is the identification point.
   */
  const Point_1& reference_direction() const { return m_reference; }

  /*! determines whether a direction is the identification point, that is,
   * whether it is equal to the reference direction.
   */
  class Is_on_x_identification_1 {
  public:
    explicit Is_on_x_identification_1(const Direction_geom_traits_2_adaptor_1& traits) : m_traits(traits) {}

    bool operator()(const Point_1& d) const
    { return m_traits.geom_traits_2().equal_2_object()(d, m_traits.reference_direction()); }

  private:
    const Direction_geom_traits_2_adaptor_1& m_traits;
  };

  /*! obtains an `Is_on_x_identification_1` functor.
   */
  Is_on_x_identification_1 is_on_x_identification_1_object() const { return Is_on_x_identification_1(*this); }

  /*! compares two directions in the counterclockwise order that starts right
   * after the reference direction.
   * \pre Neither direction is the identification point.
   */
  class Compare_x_1 {
  public:
    explicit Compare_x_1(const Direction_geom_traits_2_adaptor_1& traits) : m_traits(traits) {}

    Comparison_result operator()(const Point_1& d1, const Point_1& d2) const {
      CGAL_precondition(! m_traits.is_on_x_identification_1_object()(d1));
      CGAL_precondition(! m_traits.is_on_x_identification_1_object()(d2));
      const auto& ref = m_traits.reference_direction();
      auto ccw_in_between = m_traits.geom_traits_2().counterclockwise_in_between_2_object();
      if (ccw_in_between(d1, ref, d2)) return SMALLER;
      return ccw_in_between(d2, ref, d1) ? LARGER : EQUAL;
    }

  private:
    const Direction_geom_traits_2_adaptor_1& m_traits;
  };

  /*! obtains a `Compare_x_1` functor.
   */
  Compare_x_1 compare_x_1_object() const { return Compare_x_1(*this); }

  /*! determines whether two directions are equal.
   */
  class Equal_1 {
  public:
    explicit Equal_1(const Direction_geom_traits_2_adaptor_1& traits) : m_traits(traits) {}

    bool operator()(const Point_1& d1, const Point_1& d2) const
    { return m_traits.geom_traits_2().equal_2_object()(d1, d2); }

  private:
    const Direction_geom_traits_2_adaptor_1& m_traits;
  };

  /*! obtains an `Equal_1` functor.
   */
  Equal_1 equal_1_object() const { return Equal_1(*this); }

private:
  Shared_geom_traits_2 m_traits_2;
  Point_1 m_reference;          // the identification point
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
