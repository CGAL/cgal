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
//
// Author(s) : Dmitry Anisimov, David Bommes, Kai Hormann, and Pierre Alliez.

/*!
  \file Triangle_coordinates_2.h
*/

#ifndef CGAL_TRIANGLE_COORDINATES_2_H
#define CGAL_TRIANGLE_COORDINATES_2_H

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

// Examples: see the User Manual here - http://doc.cgal.org/latest/Manual/index.html.
// [1] Reference: Weisstein, Eric W. "Barycentric Coordinates." From MathWorld--A Wolfram Web Resource. http://mathworld.wolfram.com/BarycentricCoordinates.html.

/*!
 * \ingroup PkgBarycentric_coordinates_2
 * The class `Triangle_coordinates_2` implements barycentric coordinates ( <a href="http://mathworld.wolfram.com/BarycentricCoordinates.html" target=blanc>[1]</a>,
 * <a href="http://en.wikipedia.org/wiki/Barycentric_coordinate_system" target=blanc>[2]</a> ) with respect to an arbitrary non-degenerate triangle in the plane.
 * This class is parameterized by a traits class `Traits`.

\tparam Traits must be a model of the concept `BarycentricTraits_2`.

*/

template<class Traits> 
    class Triangle_coordinates_2
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
    Triangle_coordinates_2(const Point_2 &first_vertex, const Point_2 &second_vertex, const Point_2 &third_vertex, const Traits &b_traits = Traits()) :
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
        inline boost::optional<OutputIterator> operator()(const Point_2 &query_point, OutputIterator output)
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
    // This function accepts a container of the type <a href="http://en.cppreference.com/w/cpp/container/vector">`std::vector`</a> 
    // and returns an iterator of the type <a href="http://en.cppreference.com/w/cpp/iterator/back_insert_iterator">`std::back_insert_iterator`</a>
    // that is placed past-the-end of the resulting sequence of coordinate values.
    inline boost::optional<std::back_insert_iterator<std::vector<FT> > > operator()(const Point_2 &query_point, std::vector<FT> &output_vector)
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
        boost::optional<OutputIterator> triangle_coordinates_2(const Point_2 &query_point, OutputIterator &output)
    {
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

        // Output all coordinates.
        return boost::optional<OutputIterator>(output);
   }
};

// Global functions

/*!
 * \relates Triangle_coordinates_2
 * This is a global function that takes three vertices of a triangle and computes triangle coordinates at a given query point with respect to these vertices.
 
\tparam Traits must be a model of the concept `BarycentricTraits_2`.

*/

template<class Traits>
    inline CGAL::cpp11::array<typename Traits::FT,3> compute_triangle_coordinates_2(const typename Traits::Point_2 &first_vertex, const typename Traits::Point_2 &second_vertex, const typename Traits::Point_2 &third_vertex, const typename Traits::Point_2 &query_point, const Traits &barycentric_traits = Traits())
{
    // Some predefined functions.
    typename Traits::Compute_area_2 area_2 = barycentric_traits.compute_area_2_object();

    // Number type.
    typedef typename Traits::FT FT;

    // Compute some related sub-areas.
    const FT area_second = area_2(second_vertex, third_vertex, query_point);
    const FT area_third  = area_2(third_vertex , first_vertex, query_point);

    // Compute the total inverted area of the triangle.
    const FT inverted_total_area = FT(1) / area_2(first_vertex, second_vertex, third_vertex);

    // Compute the first and second coordinate functions.
    const FT b_first  = area_second * inverted_total_area;
    const FT b_second = area_third  * inverted_total_area;

    // Return the CGAL::cpp11::array<FT,3> type of coordinates.
    return CGAL::make_array(b_first, b_second, FT(1) - b_first - b_second);
}

} // namespace Barycentric_coordinates

} // namespace CGAL

#endif // CGAL_TRIANGLE_COORDINATES_2_H
