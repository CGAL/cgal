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

#ifndef ORTHTREE_ORTHTREE_TRAITS_BASE_H
#define ORTHTREE_ORTHTREE_TRAITS_BASE_H

#include <CGAL/license/Orthtree.h>

#include <CGAL/Bbox_2.h>
#include <CGAL/Point_set_2.h>
#include <CGAL/Orthtree/Cartesian_ranges.h>
#include <CGAL/Orthtree.h>

namespace CGAL {

/*!
  \ingroup PkgOrthtreeTraits

  The class `Orthtree_traits_base` is a base class providing common choices for types and functors.
  The base class is extended by `CGAL::Orthtree_traits_point<GeomTraits, PointRange, PointMap, dimension>` and by `CGAL::Orthtree_traits_face_graph<PolygonMesh, VertexPointMap>`.

  \tparam GeomTraits a model of `Kernel`.
  \tparam dim dimension of the ambient Euclidean space.

  \sa `CGAL::Orthtree_traits_point<GeomTraits, PointRange, PointMap, dim>`
  \sa `CGAL::Orthtree_traits_face_graph<PolygonMesh, VertexPointMap>`
*/

template <typename GeomTraits, int dim>
struct Orthtree_traits_base {
  /// \name Types
  /// @{
  using Node_index = std::size_t;
  using Kernel = GeomTraits;
  static constexpr int dimension = dim;
  using FT = typename GeomTraits::FT;
  using Point_d = typename GeomTraits::Point_d;
  using Bbox_d = typename GeomTraits::Iso_box_d;
  using Sphere_d = typename GeomTraits::Sphere_d;
  using Cartesian_const_iterator_d = typename GeomTraits::Cartesian_const_iterator_d;
  /*!
   * Adjacency type.
   *
   * \note This type is used to identify adjacency directions with
   * easily understandable keywords (left, right, up, down, ...) and is thus
   * mainly useful in 2D and 3D. In
   * higher dimensions, such keywords do not exist and this type is
   * simply an integer. Conversions from this integer to bitsets still
   * work but do not provide any user-friendly API for adjacency selection.
   *
   * Two directions along each axis in %Cartesian space, relative to a node.
   *
   * Directions are mapped to numbers as 3-bit integers in the 3D case or as 2-bit integers in the 2D case.
   * In the 3d case the numbers 6 and 7 are not used because there are only 6 different directions.
   *
   * The first two bits indicate the axis (00 = x, 01 = y, 10 = z),
   * the third bit indicates the direction along that axis (0 = -, 1 = +).
   *
   * The following diagram and table showing the 3D case may be a useful reference (2D case is identical with one dimension less):
   *
   *            3 *
   *              |  * 4
   *              | /                  y+
   *              |/                   *
   *     0 *------+------* 1           |
   *             /|                    |
   *            / |                    +-----* x+
   *         5 *  |                   /
   *              * 2                /
   *                                * z+
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
  using Adjacency = int;
  /// @}

  auto construct_point_d_object() const {
    return [](auto... Args) -> Point_d {
      std::initializer_list<FT> args_list{Args...};
      return Point_d{static_cast<int>(args_list.size()), args_list.begin(), args_list.end()};
    };
  }
};

template <typename GeomTraits>
struct Orthtree_traits_base<GeomTraits, 2> {
  using Node_index = std::size_t;
  using Kernel = GeomTraits;
  static constexpr int dimension = 2;
  using FT = typename GeomTraits::FT;
  using Point_d = typename GeomTraits::Point_2;
  using Bbox_d = typename GeomTraits::Iso_rectangle_2;
  using Sphere_d = typename GeomTraits::Circle_2;
  using Cartesian_const_iterator_d = typename GeomTraits::Cartesian_const_iterator_2;

  enum Adjacency {
    LEFT,
    RIGHT,
    DOWN,
    UP
  };

  auto construct_point_d_object() const {
    return [](const FT& x, const FT& y) -> Point_d {
      return {x, y};
    };
  }
};

template <typename GeomTraits>
struct Orthtree_traits_base<GeomTraits, 3> {
  using Node_index = std::size_t;
  using Kernel = GeomTraits;
  static constexpr int dimension = 3;
  using FT = typename GeomTraits::FT;
  using Point_d = typename GeomTraits::Point_3;
  using Bbox_d = typename GeomTraits::Iso_cuboid_3;
  using Sphere_d = typename GeomTraits::Sphere_3;
  using Cartesian_const_iterator_d = typename GeomTraits::Cartesian_const_iterator_3;

  enum Adjacency {
    LEFT,
    RIGHT,
    DOWN,
    UP,
    BACK,
    FRONT
  };

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

  auto construct_point_d_object() const {
    return [](const FT& x, const FT& y, const FT& z) -> Point_d {
      return {x, y, z};
    };
  }
};

}

#endif //ORTHTREE_ORTHTREE_TRAITS_BASE_H
