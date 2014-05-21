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

// CGAL headers.
#include <CGAL/assertions.h>
#include <CGAL/Kernel/global_functions_2.h>

// CGAL namespace.
namespace CGAL {

// Barycentric coordinates namespace.
namespace Barycentric_coordinates {

// Examples: see the User Manual here - http://doc.cgal.org/latest/Manual/index.html.
// [1] Reference: Weisstein, Eric W. "Barycentric Coordinates." From MathWorld--A Wolfram Web Resource. http://mathworld.wolfram.com/BarycentricCoordinates.html.

/*!
 * \ingroup PkgBarycentric_coordinates_2
 * The class Triangle_coordinates_2 implements barycentric coordinates ( <a href="http://mathworld.wolfram.com/BarycentricCoordinates.html" target=blanc>[1]</a>,
 * <a href="http://en.wikipedia.org/wiki/Barycentric_coordinate_system" target=blanc>[2]</a> ) with respect to an arbitrary non-degenerate triangle in the plane.
 * This class is parameterized by a traits class `Traits`.
 *
 * \pre The triangle's vertices must be ordered.

 \cgalHeading{Requirements}

 <OL>
 <LI> `Traits` class must contain the following subset of types:
 <UL>
 <LI> `Traits::Point_2` - the type of a point used internally in the class, which is equivalent to the type `CGAL::Point_2`;
 <LI> `Traits::Point_3` - the type of a point used in the static function static_compute() to return coordinate values.
 <LI> `Traits::Point_d` - the user-defined type of a point;
 <LI> `Traits::FT` - the type of a coordinate value;
 </UL>
 </OL>

 */

template<class Traits> 
    class Triangle_coordinates_2
{

public:

    /// \name Types
    /// @{

    /// Number type.
    typedef typename Traits::FT      Scalar;

    /// Type of a general point.
    typedef typename Traits::Point_d Point_d;

    /// Type of 3D point.
    typedef typename Traits::Point_3 Point_3;

    /// Type of 2D point.
    typedef typename Traits::Point_2 Point_2;

    /// @}

    /// \name Creation
    /// @{

    /// Creates the class `Triangle_coordinates_2` that implements the behaviour of triangle coordinates with respect to an arbitrary non-degenerate triangle in the plane.
    /// The triangle is given by its three vertices.
    /// \pre Triangle is not degenerate.
    Triangle_coordinates_2(const Point_d &first_vertex, const Point_d &second_vertex, const Point_d &third_vertex) :
        barycentric_traits(Traits()),
        vertex_0(barycentric_traits.project(first_vertex)),
        vertex_1(barycentric_traits.project(second_vertex)),
        vertex_2(barycentric_traits.project(third_vertex)) 
    {
        CGAL_precondition( CGAL::area(vertex_0, vertex_1, vertex_2) != Scalar(0) );
    }

    /// @}

    /// \name Computation of basis functions
    /// @{

    /// Computes triangle barycentric coordinates for a chosen query point with respect to all three vertices of the triangle.
    /// This function accepts any STL like iterator, which complies with the `Iterator` concept.
    template<class Iterator>
        inline std::pair<Iterator, bool> compute(const Point_d &query_point, Iterator output)
    {
        return triangle_coordinates_2(barycentric_traits.project(query_point), output);
    }

    /// Computes triangle barycentric coordinates for a chosen query point with respect to all three vertices of the triangle.
    /// This function accepts a container of the type <a href="http://en.cppreference.com/w/cpp/container/vector">`std::vector`</a> 
    /// and returns an iterator of the type <a href="http://en.cppreference.com/w/cpp/iterator/back_insert_iterator">`std::back_insert_iterator`</a>
    /// that is placed past-the-end of the resulting sequence of coordinate values.
    inline std::pair<std::back_insert_iterator<std::vector<Scalar> >, bool> compute(const Point_d &query_point, std::vector<Scalar> &output_vector)
    {
        output_vector.reserve(output_vector.size() + 3);
        typedef typename std::back_insert_iterator<std::vector<Scalar> > Iterator;
        Iterator output = std::back_inserter(output_vector);
        return triangle_coordinates_2(barycentric_traits.project(query_point), output);
    }

    /// This is a static function which takes three vertices of a triangle and computes triangle coordinates at a given query point with respect to these vertices.
    /// These three coordinate values are returned as a point of the type `CGAL::Point_3`.
    /// The function also requires a traits class that converts a user-defined type `Traits::Point_d` of the triangle's vertices and query point to the type `CGAL::Point_2` used internally.
    static inline typename Traits::Point_3 static_compute(const typename Traits::Point_d &first_vertex, const typename Traits::Point_d &second_vertex, const typename Traits::Point_d &third_vertex, const typename Traits::Point_d &query_point, const Traits &traits_class)
    {
        return static_triangle_coordinates_2(traits_class.project(first_vertex), traits_class.project(second_vertex), traits_class.project(third_vertex), traits_class.project(query_point));
    }

    /// @}

    /// \name Information functions
    /// @{

    /// Print some information about the used triangle and triangle coordinates.
    void print_information(std::ostream &output_stream = std::cout) const
    {
        output_stream << std::endl << "INFORMATION: " << std::endl;

        output_stream << std::endl << "DATA STRUCTURE: " << std::endl << std::endl;
        output_stream << "The internal data structure is triangle." << std::endl;

        output_stream << std::endl << "DEGENERACY: " << std::endl << std::endl;
        if(CGAL::area(vertex_0, vertex_1, vertex_2) != Scalar(0)) output_stream << "This triangle is not degenerate." << std::endl;
        else std::cout << "This triangle is degenerate. The correct computation is not expected!" << std::endl;

        output_stream << std::endl << "TYPE OF COORDINATES: " << std::endl << std::endl;
        output_stream << "The coordinate functions to compute are triangle coordinates." << std::endl;

        output_stream << std::endl << "INFORMATION ABOUT COORDINATES: " << std::endl << std::endl;
        output_stream << "Triangle coordinates can be computed exactly for an arbitrary point in the plane." << std::endl;
    }

    /// @}

private:

    // Internal global variables.
    const Traits barycentric_traits;

    const Point_2 vertex_0;
    const Point_2 vertex_1;
    const Point_2 vertex_2;

    Scalar area_second;
    Scalar area_third;
    Scalar inverted_total_area;

    Scalar b_first;
    Scalar b_second;

    // Compute triangle coordinates.
    template<class Iterator>
        std::pair<Iterator, bool> triangle_coordinates_2(const Point_2 &query_point, Iterator &output)
    {
        // Compute some related sub-areas.
        area_second = CGAL::area(vertex_1, vertex_2, query_point);
        area_third  = CGAL::area(vertex_2, vertex_0, query_point);

        // Compute the total inverted area of the triangle.
        inverted_total_area = Scalar(1) / CGAL::area(vertex_0, vertex_1, vertex_2);

        // Compute the first and second coordinate functions.
        b_first  = area_second * inverted_total_area;
        b_second = area_third  * inverted_total_area;

        *output = b_first;
        ++output;

        *output = b_second;
        ++output;

        // Compute the last - third coordinate using the partition of unity property.
        *output = Scalar(1) - b_first - b_second;

        // Output all coordinates.
        return std::make_pair(output, true);
   }

    // ...detailed description...
    static typename Traits::Point_3 static_triangle_coordinates_2(const typename Traits::Point_2 &vertex_0, const typename Traits::Point_2 &vertex_1, const typename Traits::Point_2 &vertex_2, const typename Traits::Point_2 &query_point)
    {
        // Number type.
        typedef typename Traits::FT      Scalar;

        // Point type.
        typedef typename Traits::Point_3 Point_3;

        // Compute some related sub-areas.
        const Scalar area_second = CGAL::area(vertex_1, vertex_2, query_point);
        const Scalar area_third  = CGAL::area(vertex_2, vertex_0, query_point);

        // Compute the total inverted area of the triangle.
        const Scalar inverted_total_area = Scalar(1) / CGAL::area(vertex_0, vertex_1, vertex_2);

        // Compute the first and second coordinate functions.
        const Scalar b_first  = area_second * inverted_total_area;
        const Scalar b_second = area_third  * inverted_total_area;

        // Return the CGAL::Point_3 type of coordinates.
        return Point_3(b_first, b_second, Scalar(1) - b_first - b_second);
    }
};

} // namespace Barycentric_coordinates

} // namespace CGAL

#endif // CGAL_TRIANGLE_COORDINATES_2_H