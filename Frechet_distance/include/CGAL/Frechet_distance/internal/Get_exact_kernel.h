// Copyright (c) 2024 Max-Planck-Institute Saarbruecken (Germany), GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : André Nusser <anusser@mpi-inf.mpg.de>
//                 Marvin Künnemann <marvin@mpi-inf.mpg.de>
//                 Karl Bringmann <kbringma@mpi-inf.mpg.de>
//                 Andreas Fabri
// =============================================================================

#pragma once

#include <CGAL/license/Frechet_distance.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_rational.h>

#include <CGAL/STL_Extension/internal/Has_nested_type_Has_filtered_predicates_tag.h>


namespace CGAL {
namespace Frechet_distance_ {
namespace internal {

    template <typename T>
struct Get_exact_kernel {

    using K = typename T::Kernel;

    static auto get_is_filtered()
    {
      if constexpr (::CGAL::internal::Has_nested_type_Has_filtered_predicates_tag<K>::value)
      {
        if constexpr (K::Has_filtered_predicates_tag::value) {
            return std::true_type();
        }
        else
          return std::false_type();
      }
      else
        return std::false_type();
    }

    static constexpr bool is_filtered = decltype(get_is_filtered())::value;
    static constexpr bool  is_floating_point = std::is_floating_point_v<typename K::FT>;

    static auto get_type()
    {
      if constexpr (is_filtered)
      {
        return typename K::Exact_kernel{};
      }
      else
      {
        if constexpr (is_floating_point)
            return CGAL::Simple_cartesian<CGAL::Exact_rational>{};
        else
            return K{};
      }
    }

    using type = decltype(get_type());

};

} } } // namespace CGAL::Frechet_distance_::internal