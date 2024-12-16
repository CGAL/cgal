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

  /*
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

#include <CGAL/disable_warnings.h>

namespace CGAL {
namespace Barycentric_coordinates {

#if !defined(CGAL_NO_DEPRECATED_CODE) || defined(DOXYGEN_RUNNING)

  /*!
  * \ingroup PkgBarycentricCoordinates2RefDeprecated
  * The class `Triangle_coordinates_2` implements barycentric coordinates ( <a href="https://mathworld.wolfram.com/BarycentricCoordinates.html" target=blanc>[1]</a>,
  * <a href="https://en.wikipedia.org/wiki/Barycentric_coordinate_system" target=blanc>[2]</a> ) with respect to an arbitrary non-degenerate triangle in the plane.
  * This class is parameterized by a traits class `Traits`.

  * \deprecated This part of the package is deprecated since the version 5.4 of \cgal.

  \tparam Traits must be a model of the concept `BarycentricTraits_2`.

  */
  template<class Traits>
  class
  Triangle_coordinates_2
  {

  public:

    /// \name Types
    /// @{

    /// Number type.
    typedef typename Traits::FT      FT;

    /// Point type.
    typedef typename Traits::Point_2 Point_2;

    #ifdef DOXYGEN_RUNNING
      /// Range of vertices in a triangle.
      /// This type is a model of the concept `Range`. Its iterator type is `RandomAccessIterator`, and its value type is `Traits::Point_2`.
      typedef unspecified_type Vertex_range;
    #else
      typedef std::vector<Point_2> Vertex_range;
    #endif

    /// @}

    /// \name Creation
    /// @{

    /// Creates the class `Triangle_coordinates_2` that implements triangle coordinates with respect to an arbitrary non-degenerate triangle in the plane.
    /// The triangle is given by its three vertices.
    /// \pre Triangle is not degenerate.

    #ifndef DOXYGEN_RUNNING
      CGAL_DEPRECATED_MSG("This part of the package is deprecated since the version 5.4 of CGAL!")
    #endif
    Triangle_coordinates_2(
      const Point_2 &first_vertex,
      const Point_2 &second_vertex,
      const Point_2 &third_vertex,
      const Traits &b_traits = Traits()) :
    vertex(),
    barycentric_traits(b_traits),
    area_2(barycentric_traits.compute_area_2_object()),
    collinear_2(barycentric_traits.collinear_2_object())
    {
      CGAL_precondition( !collinear_2(first_vertex, second_vertex, third_vertex) );

      vertex.resize(3);
      vertex[0] = first_vertex;
      vertex[1] = second_vertex;
      vertex[2] = third_vertex;
    }

    /// @}

    /// \name Computation
    /// @{

    /// Computes triangle barycentric coordinates for a chosen query point with respect to all three vertices of the triangle.
    /// Computed coordinates are stored in the output iterator `output`.
    template<class OutputIterator>
    inline std::optional<OutputIterator> operator()(
      const Point_2 &query_point, OutputIterator output)
    {
      return triangle_coordinates_2(query_point, output);
    }

    /// @}

    /// \name Endpoint Accessors
    /// @{

    /// Returns all the vertices of the triangle.
    inline const Vertex_range& vertices() const
    {
      return vertex;
    }

    /// Returns the first vertex of the triangle.
    inline const Point_2& first_vertex() const
    {
      return vertex[0];
    }

    /// Returns the second vertex of the triangle.
    inline const Point_2& second_vertex() const
    {
      return vertex[1];
    }

    /// Returns the third vertex of the triangle.
    inline const Point_2& third_vertex() const
    {
      return vertex[2];
    }

    /// @}

    // Computes triangle barycentric coordinates for a chosen query point with respect to all three vertices of the triangle.
    // This function accepts a container of the type <a href="https://en.cppreference.com/w/cpp/container/vector">`std::vector`</a>
    // and returns an iterator of the type <a href="https://en.cppreference.com/w/cpp/iterator/back_insert_iterator">`std::back_insert_iterator`</a>
    // that is placed past-the-end of the resulting sequence of coordinate values.
    inline std::optional<std::back_insert_iterator<std::vector<FT> > > operator()(
      const Point_2 &query_point, std::vector<FT> &output_vector)
    {
      output_vector.reserve(output_vector.size() + 3);
      typedef typename std::back_insert_iterator<std::vector<FT> > OutputIterator;
      OutputIterator output = std::back_inserter(output_vector);
      return triangle_coordinates_2(query_point, output);
    }

    // Information Functions

    // This function prints some information about the used triangle and triangle coordinates.
    void print_information(std::ostream &output_stream = std::cout) const
    {
      output_stream << std::endl << "INFORMATION: " << std::endl;

      output_stream << std::endl << "DATA STRUCTURE: " << std::endl << std::endl;
      output_stream << "The internal data structure is triangle." << std::endl;

      output_stream << std::endl << "DEGENERACY: " << std::endl << std::endl;
      if(!collinear_2(vertex[0], vertex[1], vertex[2])) output_stream << "This triangle is not degenerate." << std::endl;
      else std::cout << "This triangle is degenerate. The correct computation is not expected!" << std::endl;

      output_stream << std::endl << "TYPE OF COORDINATES: " << std::endl << std::endl;
      output_stream << "The coordinate functions to be computed are triangle coordinates." << std::endl;

      output_stream << std::endl << "INFORMATION ABOUT COORDINATES: " << std::endl << std::endl;
      output_stream << "Triangle coordinates can be computed exactly for an arbitrary point in the plane." << std::endl;
    }

  private:

    // Internal global variables.
    Vertex_range vertex;

    const Traits &barycentric_traits;

    FT area_second;
    FT area_third;
    FT inverted_total_area;

    FT b_first;
    FT b_second;

    typename Traits::Compute_area_2 area_2;
    typename Traits::Collinear_2 collinear_2;

    // Compute triangle coordinates.
    template<class OutputIterator>
    std::optional<OutputIterator> triangle_coordinates_2(
      const Point_2 &query_point, OutputIterator &output) {

      // Compute some related sub-areas.
      area_second = area_2(vertex[1], vertex[2], query_point);
      area_third  = area_2(vertex[2], vertex[0], query_point);

      // Compute the total inverted area of the triangle.
      inverted_total_area = FT(1) / area_2(vertex[0], vertex[1], vertex[2]);

      // Compute the first and second coordinate functions.
      b_first  = area_second * inverted_total_area;
      b_second = area_third  * inverted_total_area;

      *output = b_first;
      ++output;

      *output = b_second;
      ++output;

      // Compute the last = third coordinate, using the partition of unity property.
      *output = FT(1) - b_first - b_second;
      ++output;

      // Output all coordinates.
      return std::optional<OutputIterator>(output);
    }
  };

  /*!
  * \relates Triangle_coordinates_2
  * This is a global function that takes three vertices of a triangle and computes triangle coordinates at a given query point with respect to these vertices.

  * \deprecated This part of the package is deprecated since the version 5.4 of \cgal.

  \tparam Traits must be a model of the concept `BarycentricTraits_2`.

  */
  template<class Traits>
  #ifndef DOXYGEN_RUNNING
  CGAL_DEPRECATED_MSG("This part of the package is deprecated since the version 5.4 of CGAL!")
  #endif
  inline std::array<typename Traits::FT,3> compute_triangle_coordinates_2(
    const typename Traits::Point_2 &first_vertex,
    const typename Traits::Point_2 &second_vertex,
    const typename Traits::Point_2 &third_vertex,
    const typename Traits::Point_2 &query_point,
    const Traits &barycentric_traits = Traits()) {

    // Some predefined functions.
    typename Traits::Compute_area_2 area_2 = barycentric_traits.compute_area_2_object();

    // Compute some related sub-areas.
    const typename Traits::FT area_second = area_2(second_vertex, third_vertex, query_point);
    const typename Traits::FT area_third  = area_2(third_vertex , first_vertex, query_point);

    // Compute the total inverted area of the triangle.
    const typename Traits::FT inverted_total_area = typename Traits::FT(1) / area_2(first_vertex, second_vertex, third_vertex);

    // Compute the first and second coordinate functions.
    const typename Traits::FT b_first  = area_second * inverted_total_area;
    const typename Traits::FT b_second = area_third  * inverted_total_area;

    // Return the std::array<FT,3> type of coordinates.
    return CGAL::make_array(b_first, b_second, typename Traits::FT(1) - b_first - b_second);
  }

#endif // CGAL_NO_DEPRECATED_CODE

} // namespace Barycentric_coordinates
} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_BARYCENTRIC_TRIANGLE_COORDINATES_2_H
