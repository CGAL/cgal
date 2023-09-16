// Copyright (c) 2023  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Jackson Campolattaro

#ifndef ORTHTREE_EXAMPLES_ORTHTREE_TRAITS_D_BASE_H
#define ORTHTREE_EXAMPLES_ORTHTREE_TRAITS_D_BASE_H

#include <CGAL/license/Orthtree.h>

#include <CGAL/Dimension.h>
#include <CGAL/Orthtree/Cartesian_ranges.h>
#include <CGAL/Kernel_d/Point_d.h>
#include <CGAL/Kernel_d/Sphere_d.h>

namespace CGAL {

/*!
  \ingroup PkgOrthtreeTraits

  The class `Orthtree_traits_d_base` can be subclassed for easier implementation of a dd OrthtreeTraits concept.

  \tparam GeomTraits model of `Kernel`.

  \cgalModels{OrthtreeTraits}
  \sa `CGAL::Quadtree`
  \sa `CGAL::Orthtree_traits_point_d`
*/
template <typename K, typename DimensionTag>
struct Orthtree_traits_d_base {
public:

  /// \name Types
  /// @{

  using Dimension = DimensionTag;
  using FT = typename K::FT;
  using Point_d = typename K::Point_d;
  using Bbox_d = typename K::Iso_box_d;
  using Sphere_d = typename K::Sphere_d;
  using Cartesian_const_iterator_d = typename K::Cartesian_const_iterator_d;

  /*!
    Adjacency type.

    \note This type is used to identify adjacency directions with
    easily understandable keywords (left, right, up, etc.) and is thus
    mainly useful for `Orthtree_traits_2` and `Orthtree_traits_3`. In
    higher dimensions, such keywords do not exist and this type is
    simply an integer. Conversions from this integer to bitsets still
    work but do not provide any easier API for adjacency selection.
  */
  using Adjacency = int;

  /// @}


  /// \name Operations
  /// @{

  auto construct_point_d_object() const {
    return [](auto... Args) -> Point_d {
      std::initializer_list<FT> args_list{Args...};
      return Point_d{args_list.size(), args_list.begin(), args_list.end()};
    };
  }

  /// @}

};

}

#endif //ORTHTREE_EXAMPLES_ORTHTREE_TRAITS_2_BASE_H
