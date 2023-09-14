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

#ifndef ORTHTREE_EXAMPLES_ORTHTREE_TRAITS_3_BASE_H
#define ORTHTREE_EXAMPLES_ORTHTREE_TRAITS_3_BASE_H

#include <CGAL/license/Orthtree.h>

#include <CGAL/Dimension.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Orthtree/Cartesian_ranges.h>

namespace CGAL {

/*!
  \ingroup PkgOrthtreeTraits

  The class `Orthtree_traits_3_base` can be subclassed for easier implementation of a 3d OrthtreeTraits concept.

  \tparam GeomTraits model of `Kernel`.

  \cgalModels{OrthtreeTraits}
  \sa `CGAL::Quadtree`
  \sa `CGAL::Orthtree_traits_point_3`
*/
template <typename K>
struct Orthtree_traits_3_base {
public:

  /// \name Types
  /// @{

  using GeomTraits = K;
  using Dimension = Dimension_tag<3>;
  using FT = typename K::FT;
  using Point_d = typename K::Point_3;
  using Bbox_d = typename K::Iso_cuboid_3;
  using Sphere_d = typename K::Sphere_3;
  using Cartesian_const_iterator_d = typename K::Cartesian_const_iterator_3;

  /*!
   * \brief Two directions along each axis in Cartesian space, relative to a node.
   *
   * Directions are mapped to numbers as 3-bit integers,
   * though the numbers 6 and 7 are not used because there are only 6 different directions.
   *
   * The first two bits indicate the axis (00 = x, 01 = y, 10 = z),
   * the third bit indicates the direction along that axis (0 = -, 1 = +).
   *
   * The following diagram may be a useful reference:
   *
   *            3 *
   *              |  * 5
   *              | /                  y+
   *              |/                   *  z+
   *     0 *------+------* 1           | *
   *             /|                    |/
   *            / |                    +-----* x+
   *         4 *  |
   *              * 2
   *
   * This lookup table may also be helpful:
   *
   * | Direction | bitset | number | Enum  |
   * | --------- | ------ | ------ | ----- |
   * | `-x`      | 000    | 0      | LEFT  |
   * | `+x`      | 001    | 1      | RIGHT |
   * | `-y`      | 010    | 2      | DOWN  |
   * | `+y`      | 011    | 3      | UP    |
   * | `-z`      | 100    | 4      | BACK  |
   * | `+z`      | 101    | 5      | FRONT |
   */
  enum Adjacency {
    LEFT,
    RIGHT,
    DOWN,
    UP,
    BACK,
    FRONT
  };

  /// \cond SKIP_IN_MANUAL
  enum Child {
    LEFT_BOTTOM_BACK,
    RIGHT_BOTTOM_BACK,
    LEFT_TOP_BACK,
    RIGHT_TOP_BACK,
    LEFT_BOTTOM_FRONT,
    RIGHT_BOTTOM_FRONT,
    LEFT_TOP_FRONT,
    RIGHT_TOP_FRONT
  };
  /// \endcond

  /// @}


  /// \name Operations
  /// @{

  auto construct_point_d_object() const {
    return [](const FT& x, const FT& y, const FT& z) -> Point_d {
      return {x, y, z};
    };
  }

  /// @}

};

}

#endif //ORTHTREE_EXAMPLES_ORTHTREE_TRAITS_2_BASE_H
