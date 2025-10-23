// Copyright (c) 2021 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Antonio Gomes, Dmitry Anisimov
//

#ifndef CGAL_BARYCENTRIC_TRIANGLE_COORDINATES_3_H
#define CGAL_BARYCENTRIC_TRIANGLE_COORDINATES_3_H

#include <CGAL/license/Barycentric_coordinates_3.h>

#include <CGAL/Barycentric_coordinates_3/internal/utils_3.h>

namespace CGAL{
namespace Barycentric_coordinates{

/*!
  \ingroup PkgBarycentricCoordinates3RefFunctions

  \brief computes tetrahedron coordinates.

  This function computes barycentric coordinates at a given `query` point
  with respect to the points `p0`, `p1`, `p2`, and `p3`, which form a tetrahedron, that is one
  coordinate per point. The coordinates are stored in a destination range
  beginning at `oi`.

  After the coordinates \f$b_0\f$, \f$b_1\f$, \f$b_2\f$, and \f$b_2\f$ are computed, the query
  point \f$q\f$ can be obtained as \f$q = b_0p_0 + b_1p_1 + b_2p_2 + b_3p_3\f$.

  \tparam OutputIterator
  a model of `OutputIterator` that accepts values of type `GeomTraits::FT`

  \tparam GeomTraits
  a model of `BarycentricTraits_3`

  \param p0
  the first vertex of a tetrahedron

  \param p1
  the second vertex of a tetrahedron

  \param p2
  the third vertex of a tetrahedron

  \param p3
  the fourth vertex of a tetrahedron

  \param query
  a query point

  \param oi
  the beginning of the destination range with the computed coordinates

  \param traits
  a traits class with geometric objects, predicates, and constructions;
  this parameter can be omitted if the traits class can be deduced from the point type

  \return an output iterator to the element in the destination range,
  one past the last coordinate stored

  \pre `traits.compute_volume_3_object()(p0, p1, p2, p3) != 0`
*/
template<
typename OutputIterator,
typename GeomTraits>
OutputIterator tetrahedron_coordinates(
  const typename GeomTraits::Point_3& p0,
  const typename GeomTraits::Point_3& p1,
  const typename GeomTraits::Point_3& p2,
  const typename GeomTraits::Point_3& p3,
  const typename GeomTraits::Point_3& query,
  OutputIterator oi,
  const GeomTraits& traits) {

  return internal::tetrahedron_coordinates_impl(
    p0, p1, p2, p3, query, oi, traits);
}

//return iterator(infer from point_3)
template<
typename Point_3,
typename OutputIterator>
OutputIterator tetrahedron_coordinates(
  const Point_3& p0,
  const Point_3& p1,
  const Point_3& p2,
  const Point_3& p3,
  const Point_3& query,
  OutputIterator oi) {

  using GeomTraits = typename Kernel_traits<Point_3>::Kernel;
  const GeomTraits traits;
  return tetrahedron_coordinates(
    p0, p1, p2, p3, query, oi, traits);
}

/*!
  \ingroup PkgBarycentricCoordinates3RefFunctions

  \brief computes tetrahedron coordinates.

  This function computes barycentric coordinates at a given `query` point
  with respect to the points `p0`, `p1`, `p2`, and `p3`, which form a tetrahedron, that is one
  coordinate per point. The coordinates are returned in an array.

  After the coordinates \f$b_0\f$, \f$b_1\f$, \f$b_2\f$, and \f$b_3\f$ are computed, the query
  point \f$q\f$ can be obtained as \f$q = b_0p_0 + b_1p_1 + b_2p_2 + b_3p_3\f$.

  \tparam GeomTraits
  a model of `BarycentricTraits_3`

  \param p0
  the first vertex of a tetrahedron

  \param p1
  the second vertex of a tetrahedron

  \param p2
  the third vertex of a tetrahedron

  \param p3
  the fourth vertex of a tetrahedron

  \param query
  a query point

  \param traits
  a traits class with geometric objects, predicates, and constructions;
  this parameter can be omitted if the traits class can be deduced from the point type

  \return an array `std::array<GeomTraits::FT, 4>`
  with the computed coordinates

  \pre `traits.compute_volume_3_object()(p0, p1, p2, p3) != 0`
*/
template<typename GeomTraits>
std::array<typename GeomTraits::FT, 4>
tetrahedron_coordinates_in_array(
  const typename GeomTraits::Point_3& p0,
  const typename GeomTraits::Point_3& p1,
  const typename GeomTraits::Point_3& p2,
  const typename GeomTraits::Point_3& p3,
  const typename GeomTraits::Point_3& query,
  const GeomTraits& traits) {

  using FT = typename GeomTraits::FT;
  std::vector<FT> coordinates;
  coordinates.reserve(4);
  internal::tetrahedron_coordinates_impl(
    p0, p1, p2, p3, query, std::back_inserter(coordinates), traits);
  CGAL_assertion(coordinates.size() == 4);
  return {coordinates[0], coordinates[1], coordinates[2], coordinates[3]};
}

//return array (infer from point_3)
template<typename Point_3>
std::array<typename Kernel_traits<Point_3>::Kernel::FT, 4>
tetrahedron_coordinates_in_array(
  const Point_3& p0,
  const Point_3& p1,
  const Point_3& p2,
  const Point_3& p3,
  const Point_3& query) {

  using GeomTraits = typename Kernel_traits<Point_3>::Kernel;
  const GeomTraits traits;
  return tetrahedron_coordinates_in_array(
    p0, p1, p2, p3, query, traits);
}

}
}

#endif
