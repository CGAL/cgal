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

#ifndef ORTHTREE_ORTHTREE_TRAITS_BASE_FOR_DIMENSION_H
#define ORTHTREE_ORTHTREE_TRAITS_BASE_FOR_DIMENSION_H

#include <CGAL/license/Orthtree.h>

#include <CGAL/Dimension.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Point_set_2.h>
#include <CGAL/Orthtree/Cartesian_ranges.h>
#include <CGAL/Orthtree.h>

namespace CGAL {

template <typename K, typename DimensionTag>
struct Orthtree_traits_base_for_dimension {
  /// \name Types
  /// @{
  using Kernel = K;
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
      return Point_d{static_cast<int>(args_list.size()), args_list.begin(), args_list.end()};
    };
  }

  auto locate_halfspace_object() const {
    return [](const FT& a, const FT& b) -> bool {
      return a < b;
      };
  }
  /// @}
};

template <typename K>
struct Orthtree_traits_base_for_dimension<K, Dimension_tag<2>> {
  /// \name Types
  /// @{
  using Kernel = K;
  using Dimension = Dimension_tag<2>;
  using FT = typename K::FT;
  using Point_d = typename K::Point_2;
  using Bbox_d = typename K::Iso_rectangle_2;
  using Sphere_d = typename K::Circle_2;
  using Cartesian_const_iterator_d = typename K::Cartesian_const_iterator_2;
  /*!
   * \brief Two directions along each axis in %Cartesian space, relative to a node.
   *
   * Directions are mapped to numbers as 2-bit integers.
   *
   * The first bit indicates the axis (0 = x, 1 = y),
   * the second bit indicates the direction along that axis (0 = -, 1 = +).
   *
   * The following diagram may be a useful reference:
   *
   *
   *                * 3
   *               /
   *              /                    y+
   *    0 *------+------* 1           *
   *            /                    /
   *           /                    +-----* x+
   *        2 *
   *
   *
   * This lookup table may also be helpful:
   *
   * | Direction | bitset | number | Enum  |
   * | --------- | ------ | ------ | ----- |
   * | `-x`      |  00    | 0      | LEFT  |
   * | `+x`      |  01    | 1      | RIGHT |
   * | `-y`      |  10    | 2      | FRONT |
   * | `+y`      |  11    | 3      | BACK  |
   */
  enum Adjacency {
    LEFT,
    RIGHT,
    FRONT,
    BACK
  };
  /// @}

  /// \name Operations
  /// @{
  auto construct_point_d_object() const {
    return [](const FT& x, const FT& y) -> Point_d {
      return {x, y};
    };
  }

  auto locate_halfspace_object() const {
    return [](const FT& a, const FT& b) -> bool {
      return a < b;
      };
  }
  /// @}
};

template <typename K>
struct Orthtree_traits_base_for_dimension<K, Dimension_tag<3>> {
  /// \name Types
  /// @{
  using Kernel = K;
  using Dimension = Dimension_tag<3>;
  using FT = typename K::FT;
  using Point_d = typename K::Point_3;
  using Bbox_d = typename K::Iso_cuboid_3;
  using Sphere_d = typename K::Sphere_3;
  using Cartesian_const_iterator_d = typename K::Cartesian_const_iterator_3;
  /*!
   * \brief Two directions along each axis in %Cartesian space, relative to a node.
   *
   * Directions are mapped to numbers as 3-bit integers,
   * though the numbers 6 and 7 are not used because there are only 6 different directions.
   *
   * The first two bits indicate the axis (00 = x, 01 = y, 10 = z),
   * the third bit indicates the direction along that axis (0 = -, 1 = +).
   *
   * The following diagram may be a useful reference:
   *
   *            5 *
   *              |  * 3
   *              | /                  z+
   *              |/                   *  y+
   *     0 *------+------* 1           | *
   *             /|                    |/
   *            / |                    +-----* x+
   *         2 *  |
   *              * 4
   *
   * This lookup table may also be helpful:
   *
   * | Direction | bitset | number | Enum  |
   * | --------- | ------ | ------ | ----- |
   * | `-x`      | 000    | 0      | LEFT  |
   * | `+x`      | 001    | 1      | RIGHT |
   * | `-y`      | 010    | 2      | FRONT |
   * | `+y`      | 011    | 3      | BACK  |
   * | `-z`      | 100    | 4      | DOWN  |
   * | `+z`      | 101    | 5      | UP    |
   */
  enum Adjacency {
    LEFT,
    RIGHT,
    FRONT,
    BACK,
    DOWN,
    UP
  };
  /// \cond SKIP_IN_MANUAL
  enum Child {
    LEFT_BOTTOM_FRONT,
    RIGHT_BOTTOM_FRONT,
    LEFT_BOTTOM_BACK,
    RIGHT_BOTTOM_BACK,
    LEFT_TOP_FRONT,
    RIGHT_TOP_FRONT,
    LEFT_TOP_BACK,
    RIGHT_TOP_BACK
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

  auto locate_halfspace_object() const {
    return [](const FT& a, const FT& b) -> bool {
      return a < b;
      };
  }
  /// @}
};

}

#endif //ORTHTREE_ORTHTREE_TRAITS_BASE_FOR_DIMENSION_H
