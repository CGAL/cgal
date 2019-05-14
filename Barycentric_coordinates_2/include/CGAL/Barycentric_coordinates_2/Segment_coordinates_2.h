// Copyright (c) 2014 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is a part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY, AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s) : Dmitry Anisimov, David Bommes, Kai Hormann, and Pierre Alliez.

/*!
  \file Segment_coordinates_2.h
*/

#ifndef CGAL_SEGMENT_COORDINATES_2_H
#define CGAL_SEGMENT_COORDINATES_2_H

#include <CGAL/license/Barycentric_coordinates_2.h>

#include <CGAL/disable_warnings.h>

// STL headers.  
#include <vector>

// CGAL headers.
#include <CGAL/array.h>  
#include <CGAL/assertions.h>

// Boost headers.
#include <boost/optional.hpp>

// CGAL namespace.
namespace CGAL {

// Barycentric coordinates namespace.
namespace Barycentric_coordinates {

// Examples: see the User Manual - https://doc.cgal.org/latest/Manual/index.html.

/*!
 * \ingroup PkgBarycentric_coordinates_2
 * The class `Segment_coordinates_2` implements barycentric coordinates with respect to an arbitrary non-degenerate segment along an arbitrary line in the plane.
 * This class is parameterized by a traits class `Traits`.

\tparam Traits must be a model of the concept `BarycentricTraits_2`.

*/

template<class Traits> 
    class Segment_coordinates_2
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
    Segment_coordinates_2(const Point_2 &first_vertex, const Point_2 &second_vertex, const Traits &b_traits = Traits()) :
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
        inline boost::optional<OutputIterator> operator()(const Point_2 &query_point, OutputIterator output)
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
    // This function accepts a container of the type <a href="http://en.cppreference.com/w/cpp/container/vector">`std::vector`</a> 
    // and returns an iterator of the type <a href="http://en.cppreference.com/w/cpp/iterator/back_insert_iterator">`std::back_insert_iterator`</a>
    // that is placed past-the-end of the resulting sequence of coordinate values.
    inline boost::optional<std::back_insert_iterator<std::vector<FT> > > operator()(const Point_2 &query_point, std::vector<FT> &output_vector)
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
        boost::optional<OutputIterator> segment_coordinates_2(const Point_2 &query_point, OutputIterator &output)
    {   
        // Project point on the segment and compute the first coordinate.
        opposite_scalar_product = scalar_product_2(query_point - vertex[1], vertex[0] - vertex[1]);
        b_first = opposite_scalar_product / squared_distance_2(vertex[0], vertex[1]);

        // Compute the second coordinate, using the partition of unity property.
        *output = b_first;
        ++output;
        *output = FT(1) - b_first;

        // Output both coordinates.
        return boost::optional<OutputIterator>(output);
    }
};

// Global functions

/*!
   \anchor seg_coord_global
 * \relates Segment_coordinates_2
 * This is a global function that takes both vertices of a segment and computes segment coordinates at a given query point with respect to these vertices.
 
\tparam Traits must be a model of the concept `BarycentricTraits_2`.

*/

template<class Traits>
    inline CGAL::cpp11::array<typename Traits::FT,2> compute_segment_coordinates_2(const typename Traits::Point_2 &first_vertex, const typename Traits::Point_2 &second_vertex, const typename Traits::Point_2 &query_point, const Traits &barycentric_traits = Traits())
{
    // Some predefined functions.
    typename Traits::Compute_scalar_product_2 scalar_product_2 = barycentric_traits.compute_scalar_product_2_object();
    typename Traits::Compute_squared_distance_2 squared_distance_2 = barycentric_traits.compute_squared_distance_2_object();

    // Number type.
    typedef typename Traits::FT FT;

    // Project point on the segment and compute the first coordinate.
    const FT opposite_scalar_product = scalar_product_2(query_point - second_vertex, first_vertex - second_vertex);
    const FT b_first = opposite_scalar_product / squared_distance_2(first_vertex, second_vertex);

    // Return the CGAL::cpp11::array<FT,2> type of coordinates.
    return CGAL::make_array(b_first, FT(1) - b_first);
}

} // namespace Barycentric_coordinates

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_SEGMENT_COORDINATES_2_H
