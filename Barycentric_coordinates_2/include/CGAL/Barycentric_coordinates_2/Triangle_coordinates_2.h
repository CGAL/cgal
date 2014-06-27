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

\cgalHeading{Template parameters}

\tparam Traits must be a model of the concept `BarycentricTraits_2`. In particular, it must provide the functions `Kernel::Compute_area_2` and `Kernel::Collinear_2`.

*/

template<class Traits> 
    class Triangle_coordinates_2
{

public:

    /// \name Types
    /// @{

    /// Number type.
    typedef typename Traits::FT      Scalar;

    /// Point type.
    typedef typename Traits::Point_2 Point;

    /// @}

    /// \name Creation
    /// @{

    /// Creates the class `Triangle_coordinates_2` that implements the behaviour of triangle coordinates with respect to an arbitrary non-degenerate triangle in the plane.
    /// The triangle is given by its three vertices.
    /// \pre Triangle is not degenerate.
    Triangle_coordinates_2(const Point &first_vertex, const Point &second_vertex, const Point &third_vertex, const Traits &b_traits = Traits()) :
        barycentric_traits(b_traits),
        vertex_0(first_vertex),
        vertex_1(second_vertex),
        vertex_2(third_vertex),
        area_2(barycentric_traits.compute_area_2_object()),
        collinear_2(barycentric_traits.collinear_2_object())
    {
        CGAL_precondition( !collinear_2(vertex_0, vertex_1, vertex_2) );
    }

    /// @}

    /// \name Computation of Basis Functions
    /// @{

    /// Computes triangle barycentric coordinates for a chosen query point with respect to all three vertices of the triangle.
    /// This function accepts any STL like iterator, which complies with the `Iterator` concept.
    template<class Iterator>
        inline std::pair<Iterator, bool> compute(const Point &query_point, Iterator output)
    {
        return triangle_coordinates_2(query_point, output);
    }

    /// Computes triangle barycentric coordinates for a chosen query point with respect to all three vertices of the triangle.
    /// This function accepts a container of the type <a href="http://en.cppreference.com/w/cpp/container/vector">`std::vector`</a> 
    /// and returns an iterator of the type <a href="http://en.cppreference.com/w/cpp/iterator/back_insert_iterator">`std::back_insert_iterator`</a>
    /// that is placed past-the-end of the resulting sequence of coordinate values.
    inline std::pair<std::back_insert_iterator<std::vector<Scalar> >, bool> compute(const Point &query_point, std::vector<Scalar> &output_vector)
    {
        output_vector.reserve(output_vector.size() + 3);
        typedef typename std::back_insert_iterator<std::vector<Scalar> > Iterator;
        Iterator output = std::back_inserter(output_vector);
        return triangle_coordinates_2(query_point, output);
    }

    /// This is a static function that takes three vertices of a triangle and computes triangle coordinates at a given query point with respect to these vertices.
    /// The coordinate values are returned as `CGAL::cpp11::array<FT,3>`.
    /// The function also requires a traits class of the concept `BarycentricTraits_2`.
    static inline CGAL::cpp11::array<typename Traits::FT,3> static_compute(const typename Traits::Point_2 &first_vertex, const typename Traits::Point_2 &second_vertex, const typename Traits::Point_2 &third_vertex, const typename Traits::Point_2 &query_point, const Traits &barycentric_traits = Traits())
    {
        return static_triangle_coordinates_2(first_vertex, second_vertex, third_vertex, query_point, barycentric_traits);
    }

    /// @}

    /// \name Information Functions
    /// @{

    /// This function prints some information about the used triangle and triangle coordinates.
    void print_information(std::ostream &output_stream = std::cout) const
    {
        output_stream << std::endl << "INFORMATION: " << std::endl;

        output_stream << std::endl << "DATA STRUCTURE: " << std::endl << std::endl;
        output_stream << "The internal data structure is triangle." << std::endl;

        output_stream << std::endl << "DEGENERACY: " << std::endl << std::endl;
        if(!collinear_2(vertex_0, vertex_1, vertex_2)) output_stream << "This triangle is not degenerate." << std::endl;
        else std::cout << "This triangle is degenerate. The correct computation is not expected!" << std::endl;

        output_stream << std::endl << "TYPE OF COORDINATES: " << std::endl << std::endl;
        output_stream << "The coordinate functions to be computed are triangle coordinates." << std::endl;

        output_stream << std::endl << "INFORMATION ABOUT COORDINATES: " << std::endl << std::endl;
        output_stream << "Triangle coordinates can be computed exactly for an arbitrary point in the plane." << std::endl;
    }

    /// @}

private:

    // Internal global variables.
    const Traits &barycentric_traits;

    const Point &vertex_0;
    const Point &vertex_1;
    const Point &vertex_2;

    Scalar area_second;
    Scalar area_third;
    Scalar inverted_total_area;

    Scalar b_first;
    Scalar b_second;

    typename Traits::Compute_area_2 area_2;
    typename Traits::Collinear_2 collinear_2;

    // Compute triangle coordinates.
    template<class Iterator>
        std::pair<Iterator, bool> triangle_coordinates_2(const Point &query_point, Iterator &output)
    {
        // Compute some related sub-areas.
        area_second = area_2(vertex_1, vertex_2, query_point);
        area_third  = area_2(vertex_2, vertex_0, query_point);

        // Compute the total inverted area of the triangle.
        inverted_total_area = Scalar(1) / area_2(vertex_0, vertex_1, vertex_2);

        // Compute the first and second coordinate functions.
        b_first  = area_second * inverted_total_area;
        b_second = area_third  * inverted_total_area;

        *output = b_first;
        ++output;

        *output = b_second;
        ++output;

        // Compute the last = third coordinate, using the partition of unity property.
        *output = Scalar(1) - b_first - b_second;

        // Output all coordinates.
        return std::make_pair(output, true);
   }

    // This is a static function that takes three vertices of a triangle and computes triangle coordinates at a given query point with respect to these vertices.
    static CGAL::cpp11::array<typename Traits::FT,3> static_triangle_coordinates_2(const typename Traits::Point_2 &vertex_0, const typename Traits::Point_2 &vertex_1, const typename Traits::Point_2 &vertex_2, const typename Traits::Point_2 &query_point, const Traits &barycentric_traits)
    {
        // Some predefined functions.
        typename Traits::Compute_area_2 area_2 = barycentric_traits.compute_area_2_object();

        // Number type.
        typedef typename Traits::FT Scalar;

        // Compute some related sub-areas.
        const Scalar area_second = area_2(vertex_1, vertex_2, query_point);
        const Scalar area_third  = area_2(vertex_2, vertex_0, query_point);

        // Compute the total inverted area of the triangle.
        const Scalar inverted_total_area = Scalar(1) / area_2(vertex_0, vertex_1, vertex_2);

        // Compute the first and second coordinate functions.
        const Scalar b_first  = area_second * inverted_total_area;
        const Scalar b_second = area_third  * inverted_total_area;

        // Return the CGAL::cpp11::array<FT,3> type of coordinates.
        return CGAL::make_array(b_first, b_second, Scalar(1) - b_first - b_second);
    }
};

} // namespace Barycentric_coordinates

} // namespace CGAL

#endif // CGAL_TRIANGLE_COORDINATES_2_H