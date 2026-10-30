// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_ARRANGEMENT_ON_CURVE_1_OPEN_GEOM_TRAITS_2_ADAPTOR_1_H
#define CGAL_ARRANGEMENT_ON_CURVE_1_OPEN_GEOM_TRAITS_2_ADAPTOR_1_H

#include <CGAL/Arrangement_on_curve_1/Side_tags.h>
#include <CGAL/Arrangement_on_curve_1/Geom_traits_2_adaptor_1_base.h>

namespace CGAL {
namespace Arrangement_on_curve_1 {

// Adapts a 2D geometry traits class into a 1D geometry traits class whose
// parameter space is open on both sides.  The master curve is an `x`-monotone
// curve or a vertical segment, and none of its endpoints is a point of the 1D
// arrangement.  Pair with `Unbounded_topology_traits`.
template <typename GeomTraits_2>
class Open_geom_traits_2_adaptor_1 : public internal::Geom_traits_2_adaptor_1_base<GeomTraits_2> {
  using Base = internal::Geom_traits_2_adaptor_1_base<GeomTraits_2>;

public:
  using typename Base::Geom_traits_2;
  using typename Base::Shared_geom_traits_2;

  using Left_side_category = Open_side_tag;
  using Right_side_category = Open_side_tag;

  /*! constructs from a default-constructed 2D traits object.
   */
  Open_geom_traits_2_adaptor_1() : Base(std::make_shared<const Geom_traits_2>()) {}

  /*! constructs from a shared 2D traits object.
   */
  explicit Open_geom_traits_2_adaptor_1(Shared_geom_traits_2 traits_2) : Base(std::move(traits_2)) {}

  /*! constructs from a copy of a 2D traits object.
   */
  explicit Open_geom_traits_2_adaptor_1(const Geom_traits_2& traits_2) :
    Base(std::make_shared<const Geom_traits_2>(traits_2))
  {}
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

#endif
