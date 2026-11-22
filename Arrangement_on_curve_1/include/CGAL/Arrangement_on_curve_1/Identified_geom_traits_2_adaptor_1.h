// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_IDENTIFIED_GEOM_TRAITS_2_ADAPTOR_1_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_IDENTIFIED_GEOM_TRAITS_2_ADAPTOR_1_H

#include <type_traits>

#include <CGAL/Arr_tags.h>
#include <CGAL/Arrangement_on_curve_1/Side_tags.h>
#include <CGAL/Arrangement_on_curve_1/Geom_traits_2_adaptor_1_base.h>

namespace CGAL {
namespace Arrangement_on_curve_1 {

// Adapts a 2D geometry traits class whose parameter space has identified left
// and right sides (e.g., `Arr_geodesic_arc_on_sphere_traits_2`) into a 1D
// geometry traits class whose parameter space is identified.
//
// The points of the 1D arrangement lie on a single `x`-monotone curve whose
// two ends lie on the 2D identification curve: its minimal end on the left
// boundary and its maximal end on the right boundary of the 2D parameter
// space.  Both ends are the same point of the surface, which is the
// identification point of the 1D arrangement.  The other points are ordered
// lexicographically, which is the order along the curve.  Pair with
// `Identified_topology_traits`.
template <typename GeomTraits_2>
class Identified_geom_traits_2_adaptor_1 : public internal::Geom_traits_2_adaptor_1_base<GeomTraits_2> {
  using Base = internal::Geom_traits_2_adaptor_1_base<GeomTraits_2>;

  static_assert(std::is_same_v<typename GeomTraits_2::Left_side_category, Arr_identified_side_tag> &&
                std::is_same_v<typename GeomTraits_2::Right_side_category, Arr_identified_side_tag>,
                "The 2D traits must have identified left and right sides.");

public:
  using typename Base::Geom_traits_2;
  using typename Base::Point_1;
  using typename Base::Shared_geom_traits_2;

  using Left_side_category = Identified_side_tag;
  using Right_side_category = Identified_side_tag;

  /*! constructs from a default-constructed 2D traits object.
   */
  Identified_geom_traits_2_adaptor_1() : Base(std::make_shared<const Geom_traits_2>()) {}

  /*! constructs from a shared 2D traits object.
   */
  explicit Identified_geom_traits_2_adaptor_1(Shared_geom_traits_2 traits_2) : Base(std::move(traits_2)) {}

  /*! constructs from a copy of a 2D traits object.
   */
  explicit Identified_geom_traits_2_adaptor_1(const Geom_traits_2& traits_2) :
    Base(std::make_shared<const Geom_traits_2>(traits_2))
  {}

  /*! determines whether a point is the identification point, that is,
   * whether it lies on the 2D identification curve.
   */
  class Is_on_x_identification_1 {
  public:
    explicit Is_on_x_identification_1(const Geom_traits_2& traits_2) : m_traits_2(traits_2) {}

    bool operator()(const Point_1& p) const { return m_traits_2.is_on_y_identification_2_object()(p); }

  private:
    const Geom_traits_2& m_traits_2;
  };

  /*! obtains an `Is_on_x_identification_1` functor.
   */
  Is_on_x_identification_1 is_on_x_identification_1_object() const
  { return Is_on_x_identification_1(*this->m_traits_2); }
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
