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

  /*
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

#include <CGAL/disable_warnings.h>

namespace CGAL {
namespace Barycentric_coordinates {

#if !defined(CGAL_NO_DEPRECATED_CODE) || defined(DOXYGEN_RUNNING)

  /*!
  \ingroup PkgBarycentricCoordinates2RefDeprecated
  * The class `Segment_coordinates_2` implements barycentric coordinates with respect to an arbitrary non-degenerate segment along an arbitrary line in the plane.
  * This class is parameterized by a traits class `Traits`.

  * \deprecated This part of the package is deprecated since the version 5.4 of \cgal.

  \tparam Traits must be a model of the concept `BarycentricTraits_2`.

  */
  template<class Traits>
  class
  Segment_coordinates_2
  {

  public:

    /// \name Types
    /// @{

    /// Number type.
    typedef typename Traits::FT      FT;

    /// Point type.
    typedef typename Traits::Point_2 Point_2;

    #ifdef DOXYGEN_RUNNING
      /// Range of vertices in a segment.
      /// This type is a model of the concept `Range`. Its iterator type is `RandomAccessIterator`, and its value type is `Traits::Point_2`.
      typedef unspecified_type Vertex_range;
    #else
      typedef std::vector<Point_2> Vertex_range;
    #endif

    /// @}

    /// \name Creation
    /// @{

    /// Creates the class `Segment_coordinates_2` that implements segment coordinates with respect to an arbitrary non-degenerate segment along an arbitrary line in the plane.
    /// The segment is given by its two vertices.
    /// \pre Segment is not degenerate.

    #ifndef DOXYGEN_RUNNING
      CGAL_DEPRECATED_MSG("This part of the package is deprecated since the version 5.4 of CGAL!")
    #endif
    Segment_coordinates_2(
      const Point_2 &first_vertex,
      const Point_2 &second_vertex,
      const Traits &b_traits = Traits()) :
    vertex(),
    barycentric_traits(b_traits),
    scalar_product_2(barycentric_traits.compute_scalar_product_2_object()),
    squared_distance_2(barycentric_traits.compute_squared_distance_2_object()),
    equal_2(barycentric_traits.equal_2_object())
    {
      CGAL_precondition( !equal_2(first_vertex, second_vertex) );

      vertex.resize(2);
      vertex[0] = first_vertex;
      vertex[1] = second_vertex;
    }

    /// @}

    /// \name Computation
    /// @{

    /// Computes segment barycentric coordinates for a chosen query point with respect to both vertices of the segment.
    /// Computed coordinates are stored in the output iterator `output`.
    template<class OutputIterator>
    inline std::optional<OutputIterator> operator()(
      const Point_2 &query_point, OutputIterator output)
    {
      return segment_coordinates_2(query_point, output);
    }

    /// @}

    /// \name Endpoint Accessors
    /// @{

    /// Returns both vertices of the segment.
    inline const Vertex_range& vertices() const
    {
      return vertex;
    }

    /// Returns the first vertex of the segment.
    inline const Point_2& first_vertex() const
    {
      return vertex[0];
    }

    /// Returns the second vertex of the segment.
    inline const Point_2& second_vertex() const
    {
      return vertex[1];
    }

    /// @}

    // Computes segment barycentric coordinates for a chosen query point with respect to both vertices of the segment.
    // This function accepts a container of the type <a href="https://en.cppreference.com/w/cpp/container/vector">`std::vector`</a>
    // and returns an iterator of the type <a href="https://en.cppreference.com/w/cpp/iterator/back_insert_iterator">`std::back_insert_iterator`</a>
    // that is placed past-the-end of the resulting sequence of coordinate values.
    inline std::optional<std::back_insert_iterator<std::vector<FT> > > operator()(
      const Point_2 &query_point, std::vector<FT> &output_vector)
    {
      output_vector.reserve(output_vector.size() + 2);
      typedef typename std::back_insert_iterator<std::vector<FT> > OutputIterator;
      OutputIterator output = std::back_inserter(output_vector);
      return segment_coordinates_2(query_point, output);
    }

    // Information Functions

    // This function prints some information about the used segment and segment coordinates.
    void print_information(std::ostream &output_stream = std::cout) const
    {
      output_stream << std::endl << "INFORMATION: " << std::endl;

      output_stream << std::endl << "DATA STRUCTURE: " << std::endl << std::endl;
      output_stream << "The internal data structure is segment." << std::endl;

      output_stream << std::endl << "DEGENERACY: " << std::endl << std::endl;
      if(!equal_2(vertex[0], vertex[1])) output_stream << "This segment is not degenerate." << std::endl;
      else output_stream << "This segment is degenerate. The correct computation is not expected!" << std::endl;

      output_stream << std::endl << "TYPE OF COORDINATES: " << std::endl << std::endl;
      output_stream << "The coordinate functions to be computed are segment coordinates." << std::endl;

      output_stream << std::endl << "INFORMATION ABOUT COORDINATES: " << std::endl << std::endl;
      output_stream << "Segment coordinates can be computed exactly for an arbitrary point along the line supporting the used segment." << std::endl;
      output_stream << "A slight offset from the line is allowed." << std::endl;
    }

  private:

    // Internal global variables.
    Vertex_range vertex;

    const Traits &barycentric_traits;

    FT b_first;
    FT opposite_scalar_product;

    typename Traits::Compute_scalar_product_2 scalar_product_2;
    typename Traits::Compute_squared_distance_2 squared_distance_2;
    typename Traits::Equal_2 equal_2;

    // Compute segment coordinates.
    template<class OutputIterator>
    std::optional<OutputIterator> segment_coordinates_2(
      const Point_2 &query_point, OutputIterator &output) {

      // Project point on the segment and compute the first coordinate.
      opposite_scalar_product = scalar_product_2(query_point - vertex[1], vertex[0] - vertex[1]);
      b_first = opposite_scalar_product / squared_distance_2(vertex[0], vertex[1]);

      // Compute the second coordinate, using the partition of unity property.
      *output = b_first;
      ++output;
      *output = FT(1) - b_first;
      ++output;

      // Output both coordinates.
      return std::optional<OutputIterator>(output);
    }
  };

  /*!
  * \relates Segment_coordinates_2
  * This is a global function that takes both vertices of a segment and computes segment coordinates at a given query point with respect to these vertices.

  * \deprecated This part of the package is deprecated since the version 5.4 of \cgal.

  \tparam Traits must be a model of the concept `BarycentricTraits_2`.

  */
  template<class Traits>
  #ifndef DOXYGEN_RUNNING
  CGAL_DEPRECATED_MSG("This part of the package is deprecated since the version 5.4 of CGAL!")
  #endif
  inline std::array<typename Traits::FT,2> compute_segment_coordinates_2(
    const typename Traits::Point_2 &first_vertex,
    const typename Traits::Point_2 &second_vertex,
    const typename Traits::Point_2 &query_point,
    const Traits &barycentric_traits = Traits()) {

    // Some predefined functions.
    typename Traits::Compute_scalar_product_2 scalar_product_2 = barycentric_traits.compute_scalar_product_2_object();
    typename Traits::Compute_squared_distance_2 squared_distance_2 = barycentric_traits.compute_squared_distance_2_object();

    // Project point on the segment and compute the first coordinate.
    const typename Traits::FT opposite_scalar_product = scalar_product_2(query_point - second_vertex, first_vertex - second_vertex);
    const typename Traits::FT b_first = opposite_scalar_product / squared_distance_2(first_vertex, second_vertex);

    // Return the std::array<FT,2> type of coordinates.
    return CGAL::make_array(b_first, typename Traits::FT(1) - b_first);
  }

#endif // CGAL_NO_DEPRECATED_CODE

} // namespace Barycentric_coordinates
} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_BARYCENTRIC_SEGMENT_COORDINATES_2_H
