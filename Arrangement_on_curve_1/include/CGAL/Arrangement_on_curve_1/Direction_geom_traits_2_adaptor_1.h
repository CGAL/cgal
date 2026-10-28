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
// the predicate `Counterclockwise_in_between_2` (e.g., any kernel) into a 1D
// geometry traits class of the circle of directions, whose parameter space is
// identified.  The identification point is a given reference direction s, and
// the directions are ordered counterclockwise, starting right after s.
//
// The order is determined by `Counterclockwise_in_between_2` alone, which
// refers to no axis: a direction a precedes a direction b if and only if a is
// counterclockwise strictly in between s and b; a and b are equal if and only
// if neither precedes the other, and a is equal to s if and only if a is not
// counterclockwise strictly in between s and s (the open range from s to s is
// the entire circle except for s).  Pair with `Identified_topology_traits`.
template <typename GeomTraits_2>
class Direction_geom_traits_2_adaptor_1 {
public:
  using Geom_traits_2 = GeomTraits_2;
  using Point_1 = typename Geom_traits_2::Direction_2;
  using Shared_geom_traits_2 = std::shared_ptr<const Geom_traits_2>;

  using Left_side_category = Identified_side_tag;
  using Right_side_category = Identified_side_tag;

private:
  using Ccw_in_between_2 = typename Geom_traits_2::Counterclockwise_in_between_2;

public:
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
    explicit Is_on_x_identification_1(const Direction_geom_traits_2_adaptor_1& traits) :
      m_ccw_in_between(traits.m_traits_2->counterclockwise_in_between_2_object()),
      m_reference(traits.m_reference)
    {}

    bool operator()(const Point_1& d) const { return ! m_ccw_in_between(d, m_reference, m_reference); }

  private:
    Ccw_in_between_2 m_ccw_in_between;
    const Point_1& m_reference;
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
    explicit Compare_x_1(const Direction_geom_traits_2_adaptor_1& traits) :
      m_ccw_in_between(traits.m_traits_2->counterclockwise_in_between_2_object()),
      m_reference(traits.m_reference)
    {}

    Comparison_result operator()(const Point_1& d1, const Point_1& d2) const {
      CGAL_precondition(m_ccw_in_between(d1, m_reference, m_reference));
      CGAL_precondition(m_ccw_in_between(d2, m_reference, m_reference));
      if (m_ccw_in_between(d1, m_reference, d2)) return SMALLER;
      return m_ccw_in_between(d2, m_reference, d1) ? LARGER : EQUAL;
    }

  private:
    Ccw_in_between_2 m_ccw_in_between;
    const Point_1& m_reference;
  };

  /*! obtains a `Compare_x_1` functor.
   */
  Compare_x_1 compare_x_1_object() const { return Compare_x_1(*this); }

  /*! determines whether two directions are equal.
   */
  class Equal_1 {
  public:
    explicit Equal_1(const Direction_geom_traits_2_adaptor_1& traits) :
      m_ccw_in_between(traits.m_traits_2->counterclockwise_in_between_2_object())
    {}

    bool operator()(const Point_1& d1, const Point_1& d2) const { return ! m_ccw_in_between(d1, d2, d2); }

  private:
    Ccw_in_between_2 m_ccw_in_between;
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
