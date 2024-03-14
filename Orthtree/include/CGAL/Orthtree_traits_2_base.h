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

#ifndef ORTHTREE_EXAMPLES_ORTHTREE_TRAITS_2_BASE_H
#define ORTHTREE_EXAMPLES_ORTHTREE_TRAITS_2_BASE_H

#include <CGAL/license/Orthtree.h>

#include <CGAL/Dimension.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Point_set_2.h>
#include <CGAL/Orthtree/Cartesian_ranges.h>

namespace CGAL {

/*!
  \ingroup PkgOrthtreeTraits

  The class `Orthtree_traits_2_base` can be subclassed for easier implementation of a 2d OrthtreeTraits concept.

  \tparam GeomTraits model of `Kernel`.

  \cgalModels{OrthtreeTraits}
  \sa `CGAL::Quadtree`
  \sa `CGAL::Orthtree_traits_point_2`
*/
template <typename K>
struct Orthtree_traits_2_base {
public:

  /// \name Types
  /// @{

  using Dimension = Dimension_tag<2>;
  using FT = typename K::FT;
  using Point_d = typename K::Point_2;
  using Bbox_d = typename K::Iso_rectangle_2;
  using Sphere_d = typename K::Circle_2;
  using Cartesian_const_iterator_d = typename K::Cartesian_const_iterator_2;

  /*!
   * \brief Two directions along each axis in Cartesian space, relative to a node.
   *
   * Directions are mapped to numbers as 2-bit integers.
   *
   * The first bit indicates the axis (0 = x, 1 = y),
   * the second bit indicates the direction along that axis (0 = -, 1 = +).
   *
   * The following diagram may be a useful reference:
   *
   *            3 *
   *              |
   *              |                    y+
   *              |                    *
   *     0 *------+------* 1           |
   *              |                    |
   *              |                    +-----* x+
   *              |
   *              * 2
   *
   * This lookup table may also be helpful:
   *
   * | Direction | bitset | number | Enum  |
   * | --------- | ------ | ------ | ----- |
   * | `-x`      |  00    | 0      | LEFT  |
   * | `+x`      |  01    | 1      | RIGHT |
   * | `-y`      |  10    | 2      | DOWN  |
   * | `+y`      |  11    | 3      | UP    |
   */
  enum Adjacency
  {
    LEFT,
    RIGHT,
    DOWN,
    UP
  };

  /// @}

  /// \name Operations
  /// @{

  auto construct_point_d_object() const {
    return [](const FT& x, const FT& y) -> Point_d {
      return {x, y};
    };
  }

  /// @}

};

}

#endif //ORTHTREE_EXAMPLES_ORTHTREE_TRAITS_2_BASE_H
