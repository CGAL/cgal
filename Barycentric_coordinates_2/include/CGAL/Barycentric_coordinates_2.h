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

#ifndef CGAL_BARYCENTRIC_COORDINATES_2_H
#define CGAL_BARYCENTRIC_COORDINATES_2_H

#include <CGAL/license/Barycentric_coordinates_2.h>

// Deprecated headers.
// #include <CGAL/Barycentric_coordinates_2/Deprecated_headers_2.h>

// Internal includes.
#include <CGAL/Barycentric_coordinates_2/segment_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/triangle_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/boundary_coordinates_2.h>

#include <CGAL/Barycentric_coordinates_2/Delaunay_domain_2.h>
#include <CGAL/Barycentric_coordinates_2/Wachspress_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/Discrete_harmonic_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/Mean_value_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/Harmonic_coordinates_2.h>

namespace CGAL {
namespace Barycentric_coordinates {

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes 2D harmonic coordinates.

    This function first creates a triangulation of the polygon interior domain given
    the user specified `max_edge_length` parameter and several `seed` points. It then
    computes 2D harmonic coordinates at each vertex of this triangulation with respect
    to the `n` vertices of a simple `polygon`, that is one coordinate per polygon vertex.
    The coordinates are stored in a destination range beginning at `c_begin`, where each
    range element is a vector with `n` coordinates, and the size of range equals to the
    number of triangulation vertices.

    Internally, the classes `Delaunay_domain_2` and `Harmonic_coordinates_2` are used.
    If one wants to evaluate harmonic coordinates at multiple query points, which are
    not the vertices of the created triangulation, one needs to refer to those classes.

    \tparam PointRange
    a model of `ConstRange` whose iterator type is `RandomAccessIterator`
    and value type is `GeomTraits::Point_2`

    \tparam OutIterator
    a model of `OutputIterator` that accepts elements of type `std::vector<GeomTraits::FT>`

    \tparam GeomTraits
    a model of `BarycentricTraits_2`

    \param polygon
    an instance of `PointRange` with 2D points, which form a simple polygon

    \param seeds
    an instance of `PointRange`, which contains seed points indicating, which parts
    of the `polygon` should be partitioned and subdivided

    \param c_begin
    the beginning of the destination range with the computed coordinates

    \param traits
    a traits class with geometric objects, predicates, and constructions;
    this parameter can be omitted if the traits class can be deduced from the
    value type of the `PointRange`

    \param max_edge_length
    an upper bound on the length of the longest edge; the default is `0.01`

    \return an output iterator to the element in the destination range,
    one past the last vector of coordinates stored

    \pre polygon.size() >= 3
    \pre polygon is simple
  */
  template<
  typename PointRange,
  typename OutIterator,
  typename GeomTraits>
  OutIterator harmonic_coordinates_2(
    const PointRange& polygon,
    const PointRange& seeds,
    OutIterator c_begin,
    const GeomTraits& traits,
    const typename GeomTraits::FT max_edge_length =
    typename GeomTraits::FT(1) / typename GeomTraits::FT(100)) {

    using Domain =
      Delaunay_domain_2<PointRange, GeomTraits>;
    using Harmonic_coordinates_2 =
      Harmonic_coordinates_2<PointRange, Domain, GeomTraits>;
    using FT = typename GeomTraits::FT;

    Domain domain(polygon, traits);
    domain.create(max_edge_length, seeds);

    Harmonic_coordinates_2 harmonic_coordinates_2(
      polygon, domain, traits);
    harmonic_coordinates_2.compute();

    std::vector<FT> coordinates;
    coordinates.reserve(polygon.size());
    for (std::size_t k = 0; k < domain.number_of_vertices(); ++k) {
      coordinates.clear();
      harmonic_coordinates_2(k, std::back_inserter(coordinates));
      *(c_begin++) = coordinates;
    }
    return c_begin;
  }

  /// \cond SKIP_IN_MANUAL
  template<
  typename PointRange,
  typename OutIterator>
  OutIterator harmonic_coordinates_2(
    const PointRange& polygon,
    const PointRange& seeds,
    OutIterator c_begin,
    const typename Kernel_traits<typename PointRange::value_type>::Kernel::FT
    max_edge_length =
    typename Kernel_traits<typename PointRange::value_type>::Kernel::FT(1)  /
    typename Kernel_traits<typename PointRange::value_type>::Kernel::FT(100)) {

    using GeomTraits = typename Kernel_traits<
      typename PointRange::value_type>::Kernel;
    const GeomTraits traits;
    return harmonic_coordinates_2(
      polygon, seeds, c_begin, traits, max_edge_length);
  }
  /// \endcond

} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_COORDINATES_2_H
