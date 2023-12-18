// Copyright (c) 2020  GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_ORTHTREE_TRAITS_3_H
#define CGAL_ORTHTREE_TRAITS_3_H

#include <CGAL/license/Orthtree.h>

#include <CGAL/Dimension.h>
#include <CGAL/Bbox_3.h>

namespace CGAL
{

/*!
  \ingroup PkgOrthtreeTraits

  The class `Orthtree_traits_3` can be used as a template parameter of
  the `Orthtree` class.

  \tparam GeomTraits model of `Kernel`.

  \cgalModels{OrthtreeTraits}
  \sa `CGAL::Octree`
  \sa `CGAL::Orthtree_traits_2`
  \sa `CGAL::Orthtree_traits_d`
*/
template <typename GeomTraits>
struct Orthtree_traits_3
{
public:

  /// \name Types
  /// @{

  typedef Dimension_tag<3> Dimension; ///< Dimension type.
  typedef Bbox_3 Bbox_d; ///< Bounding box type.
  typedef typename GeomTraits::FT FT; ///< Number type.
  typedef typename GeomTraits::Point_3 Point_d; ///< Point type.
  typedef typename GeomTraits::Sphere_3 Sphere_d; ///< Sphere type.
  typedef typename GeomTraits::Cartesian_const_iterator_3 Cartesian_const_iterator_d; ///< An iterator over the %Cartesian coordinates.
  typedef std::array<FT, Dimension::value> Array; ///< Array type.

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
  enum Adjacency
  {
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

#ifdef DOXYGEN_RUNNING
  /*!
    Functor with an operator to construct a `Point_d` from an `Array` object.
  */
  typedef unspecified_type Construct_point_d_from_array;
#else
  struct Construct_point_d_from_array
  {
    Point_d operator() (const Array& array) const
    {
      return Point_d (array[0], array[1], array[2]);
    }
  };
#endif

#ifdef DOXYGEN_RUNNING
  /*!
    Functor with an operator to construct a `Bbox_d` from two `Array` objects (coordinates of minimum and maximum points).
  */
  typedef unspecified_type Construct_bbox_d;
#else
  struct Construct_bbox_d
  {
    Bbox_d operator() (const Array& min,
                       const Array& max) const
    {
      return Bbox_d (min[0], min[1], min[2], max[0], max[1], max[2]);
    }
  };
#endif

  /// @}

  /// \name Operations
  /// @{

  /*!
    Function used to construct an object of type `Construct_point_d_from_array`.
  */
  Construct_point_d_from_array construct_point_d_from_array_object() const
  { return Construct_point_d_from_array(); }

  /*!
    Function used to construct an object of type `Construct_bbox_d`.
  */
  Construct_bbox_d construct_bbox_d_object() const
  { return Construct_bbox_d(); }

  /// @}
};

}

#endif // CGAL_ORTHTREE_TRAITS_3_H
