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

#ifndef CGAL_BARYCENTRIC_SEGMENT_COORDINATES_2_H
#define CGAL_BARYCENTRIC_SEGMENT_COORDINATES_2_H

#include <CGAL/license/Barycentric_coordinates_2.h>

// Internal includes.
#include <CGAL/Barycentric_coordinates_2/internal/utils_2.h>

namespace CGAL {
namespace Barycentric_coordinates {

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes segment coordinates.

    This function computes barycentric coordinates at a given `query` point
    with respect to the end points `p0` and `p1` of a segment that is one
    coordinate per end point. The coordinates are stored in a destination range
    beginning at `c_begin`.

    After the coordinates \f$b_0\f$ and \f$b_1\f$ are computed, the query point \f$q\f$ can be
    obtained as \f$q = b_0p_0 + b_1p_1\f$. If \f$q\f$ does not belong to the line through \f$p_0\f$
    and \f$p_1\f$, it is projected onto this line, and only then the coordinates are
    computed. See more details in the user manual \ref compute_seg_coord "here".

    \tparam OutIterator
    a model of `OutputIterator` that accepts values of type `GeomTraits::FT`

    \tparam GeomTraits
    a model of `BarycentricTraits_2`

    \param p0
    the first vertex of a segment

    \param p1
    the second vertex of a segment

    \param query
    a query point

    \param c_begin
    the beginning of the destination range with the computed coordinates

    \param traits
    a traits class with geometric objects, predicates, and constructions;
    this parameter can be omitted if the traits class can be deduced from the point type

    \return an output iterator to the element in the destination range,
    one past the last coordinate stored

    \pre p0 != p1
  */
  template<
  typename OutIterator,
  typename GeomTraits>
  OutIterator segment_coordinates_2(
    const typename GeomTraits::Point_2& p0,
    const typename GeomTraits::Point_2& p1,
    const typename GeomTraits::Point_2& query,
    OutIterator c_begin,
    const GeomTraits& traits) {

    return internal::linear_coordinates_2(
      p0, p1, query, c_begin, traits);
  }

  /// \cond SKIP_IN_MANUAL
  template<
  typename Point_2,
  typename OutIterator>
  OutIterator segment_coordinates_2(
    const Point_2& p0,
    const Point_2& p1,
    const Point_2& query,
    OutIterator c_begin) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    const GeomTraits traits;
    return segment_coordinates_2(
      p0, p1, query, c_begin, traits);
  }
  /// \endcond

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes segment coordinates.

    This function computes barycentric coordinates at a given `query` point
    with respect to the end points `p0` and `p1` of a segment that is one
    coordinate per end point. The coordinates are returned in a pair.

    After the coordinates \f$b_0\f$ and \f$b_1\f$ are computed, the query point \f$q\f$ can be
    obtained as \f$q = b_0p_0 + b_1p_1\f$. If \f$q\f$ does not belong to the line through \f$p_0\f$
    and \f$p_1\f$, it is projected onto this line, and only then the coordinates are
    computed. See more details in the user manual \ref compute_seg_coord "here".

    \tparam GeomTraits
    a model of `BarycentricTraits_2`

    \param p0
    the first vertex of a segment

    \param p1
    the second vertex of a segment

    \param query
    a query point

    \param traits
    a traits class with geometric objects, predicates, and constructions;
    this parameter can be omitted if the traits class can be deduced from the point type

    \return a pair `std::pair<GeomTraits::FT, GeomTraits::FT>`
    with the computed coordinates

    \pre p0 != p1
  */
  template<typename GeomTraits>
  std::pair<
  typename GeomTraits::FT,
  typename GeomTraits::FT>
  segment_coordinates_in_pair_2(
    const typename GeomTraits::Point_2& p0,
    const typename GeomTraits::Point_2& p1,
    const typename GeomTraits::Point_2& query,
    const GeomTraits& traits) {

    using FT = typename GeomTraits::FT;
    std::vector<FT> coordinates;
    coordinates.reserve(2);
    internal::linear_coordinates_2(
      p0, p1, query, std::back_inserter(coordinates), traits);
    CGAL_assertion(coordinates.size() == 2);
    return std::make_pair(coordinates[0], coordinates[1]);
  }

  /// \cond SKIP_IN_MANUAL
  template<typename Point_2>
  std::pair<
  typename Kernel_traits<Point_2>::Kernel::FT,
  typename Kernel_traits<Point_2>::Kernel::FT>
  segment_coordinates_in_pair_2(
    const Point_2& p0,
    const Point_2& p1,
    const Point_2& query) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    const GeomTraits traits;
    return segment_coordinates_in_pair_2(
      p0, p1, query, traits);
  }
  /// \endcond

} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_SEGMENT_COORDINATES_2_H
