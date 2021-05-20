// Copyright (c) 2014 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Dmitry Anisimov, David Bommes, Kai Hormann, Pierre Alliez
//

#ifndef CGAL_BARYCENTRIC_TRIANGLE_COORDINATES_2_H
#define CGAL_BARYCENTRIC_TRIANGLE_COORDINATES_2_H

#include <CGAL/license/Barycentric_coordinates_2.h>

// Internal includes.
#include <CGAL/Barycentric_coordinates_2/internal/utils_2.h>

namespace CGAL {
namespace Barycentric_coordinates {

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes triangle coordinates.

    This function computes barycentric coordinates at a given `query` point
    with respect to the points `p0`, `p1`, and `p2`, which form a triangle, that is one
    coordinate per point. The coordinates are stored in a destination range
    beginning at `c_begin`.

    After the coordinates \f$b_0\f$, \f$b_1\f$, and \f$b_2\f$ are computed, the query
    point \f$q\f$ can be obtained as \f$q = b_0p_0 + b_1p_1 + b_2p_2\f$. See more details
    in the user manual \ref compute_tri_coord "here".

    \tparam OutIterator
    a model of `OutputIterator` that accepts values of type `GeomTraits::FT`

    \tparam GeomTraits
    a model of `BarycentricTraits_2`

    \param p0
    the first vertex of a triangle

    \param p1
    the second vertex of a triangle

    \param p2
    the third vertex of a triangle

    \param query
    a query point

    \param c_begin
    the beginning of the destination range with the computed coordinates

    \param traits
    a traits class with geometric objects, predicates, and constructions;
    this parameter can be omitted if the traits class can be deduced from the point type

    \return an output iterator to the element in the destination range,
    one past the last coordinate stored

    \pre area_2(p0, p1, p2) != 0
  */
  template<
  typename OutIterator,
  typename GeomTraits>
  OutIterator triangle_coordinates_2(
    const typename GeomTraits::Point_2& p0,
    const typename GeomTraits::Point_2& p1,
    const typename GeomTraits::Point_2& p2,
    const typename GeomTraits::Point_2& query,
    OutIterator c_begin,
    const GeomTraits& traits) {

    return internal::planar_coordinates_2(
      p0, p1, p2, query, c_begin, traits);
  }

  /// \cond SKIP_IN_MANUAL
  template<
  typename Point_2,
  typename OutIterator>
  OutIterator triangle_coordinates_2(
    const Point_2& p0,
    const Point_2& p1,
    const Point_2& p2,
    const Point_2& query,
    OutIterator c_begin) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    const GeomTraits traits;
    return triangle_coordinates_2(
      p0, p1, p2, query, c_begin, traits);
  }
  /// \endcond

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes triangle coordinates.

    This function computes barycentric coordinates at a given `query` point
    with respect to the points `p0`, `p1`, and `p2`, which form a triangle, that is one
    coordinate per point. The coordinates are returned in a tuple.

    After the coordinates \f$b_0\f$, \f$b_1\f$, and \f$b_2\f$ are computed, the query
    point \f$q\f$ can be obtained as \f$q = b_0p_0 + b_1p_1 + b_2p_2\f$. See more details
    in the user manual \ref compute_tri_coord "here".

    \tparam GeomTraits
    a model of `BarycentricTraits_2`

    \param p0
    the first vertex of a triangle

    \param p1
    the second vertex of a triangle

    \param p2
    the third vertex of a triangle

    \param query
    a query point

    \param traits
    a traits class with geometric objects, predicates, and constructions;
    this parameter can be omitted if the traits class can be deduced from the point type

    \return a tuple `std::tuple<GeomTraits::FT, GeomTraits::FT, GeomTraits::FT>`
    with the computed coordinates

    \pre area_2(p0, p1, p2) != 0
  */
  template<typename GeomTraits>
  std::tuple<
  typename GeomTraits::FT,
  typename GeomTraits::FT,
  typename GeomTraits::FT>
  triangle_coordinates_in_tuple_2(
    const typename GeomTraits::Point_2& p0,
    const typename GeomTraits::Point_2& p1,
    const typename GeomTraits::Point_2& p2,
    const typename GeomTraits::Point_2& query,
    const GeomTraits& traits) {

    using FT = typename GeomTraits::FT;
    std::vector<FT> coordinates;
    coordinates.reserve(3);
    internal::planar_coordinates_2(
      p0, p1, p2, query, std::back_inserter(coordinates), traits);
    CGAL_assertion(coordinates.size() == 3);
    return std::make_tuple(coordinates[0], coordinates[1], coordinates[2]);
  }

  /// \cond SKIP_IN_MANUAL
  template<typename Point_2>
  std::tuple<
  typename Kernel_traits<Point_2>::Kernel::FT,
  typename Kernel_traits<Point_2>::Kernel::FT,
  typename Kernel_traits<Point_2>::Kernel::FT>
  triangle_coordinates_in_tuple_2(
    const Point_2& p0,
    const Point_2& p1,
    const Point_2& p2,
    const Point_2& query) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    const GeomTraits traits;
    return triangle_coordinates_in_tuple_2(
      p0, p1, p2, query, traits);
  }
  /// \endcond

} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_TRIANGLE_COORDINATES_2_H
