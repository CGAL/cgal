// Copyright (c) 2013 INRIA Sophia-Antipolis (France).
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

// Barycentric coordinates headers.
#include <CGAL/barycentric_enum.h>

// CGAL namespace.
namespace CGAL {

// Barycentric coordinates namespace.
namespace Barycentric_coordinates {

// Examples: See User Manual here - http://doc.cgal.org/latest/Manual/index.html.
// [1] Reference: Weisstein, Eric W. "Barycentric Coordinates." From MathWorld--A Wolfram Web Resource. http://mathworld.wolfram.com/BarycentricCoordinates.html.

/*!
 * \ingroup PkgBarycentric_coordinates_2
 * The class Triangle_coordinates_2 implements barycentric coordinates ( <a href="http://mathworld.wolfram.com/BarycentricCoordinates.html" target=blanc>[1]</a>,
 * <a href="http://en.wikipedia.org/wiki/Barycentric_coordinate_system" target=blanc>[2]</a> ) with respect to an arbitrary non-degenerate triangle in the plane.
 * This class is parameterized by `CGAL::Triangle_2` class and `Iterator` class.
 * The latter can be any class that fulfills the requirements for an STL iterator.
 *
 * \sa `Iterator`
 *
 */

// This class does not allow the user to use different iterators to output coordinates after
// it has been created. In order to do so, we need to move template parameter Iterator in declaration of the function compute().
// Can we also use reference in front of iterator like OutputIterator &output when passing it to the function? 
template<typename Triangle_2, typename Iterator> 
    class Triangle_coordinates_2
{

public:

    // Creation.

    /// \name Types
    /// @{

    /// Type of 2D triangle.
	typedef Triangle_2               Triangle;

    /// Type of the used kernel.
	typedef typename Triangle::R     Kernel;

    /// Type of 2D point.
	typedef typename Kernel::Point_2 Point;

    /// Number type.
	typedef typename Kernel::FT      Scalar;

    /// Type of the used output iterator.
    typedef Iterator                 OutputIterator;

    /// @}

    /// \name Creation
    /// @{

    /// Creates an instance of Triangle_coordinates_2 class for a provided triangle passed as a reference.
    /// \pre `!_triangle.is_degenerate()`
    Triangle_coordinates_2(const Triangle &_triangle) : triangle(_triangle)
    {
        CGAL_precondition( !_triangle.is_degenerate() );
    }

    /// @}

    // Coordinates.

    /// \name Computation of the basis functions
    /// @{

    /// Computes Triangle barycentric coordinates for current query point with
    /// respect to all three vertices of the triangle.
    inline std::pair<OutputIterator, bool> compute(const Point &query_point, OutputIterator output)
    {
        return triangle_coordinates_2(query_point, output);
    }

    /// @}

    // Information about computed coordinates.

    /// \name Information functions
    /// @{

    /// Print some information about currently used triangle and Triangle coordinates.
    void print_info() const
    {
        std::cout << std::endl << "INFORMATION: " << std::endl;

        std::cout << std::endl << "DATA STRUCTURE: " << std::endl << std::endl;
        std::cout << "The used data structure is triangle." << std::endl;

        std::cout << std::endl << "DEGENERACY: " << std::endl << std::endl;
        if(!triangle.is_degenerate()) std::cout << "Current triangle is not degenerate." << std::endl;
        else std::cout << "Current triangle is degenerate. The correct computation is not expected!" << std::endl;

        std::cout << std::endl << "TYPE OF COORDINATES: " << std::endl << std::endl;
        std::cout << "Currently computed coordinate functions are Triangle coordinates." << std::endl;

        std::cout << std::endl << "INFORMATION ABOUT COORDINATES: " << std::endl << std::endl;
        std::cout << "Triangle coordinates can be computed exactly for an arbitrary point in the plane." << std::endl;
    }

    /// @}

private:

    // Internal global variables.
    const Triangle &triangle;

    Scalar area_second;
    Scalar area_third;
    Scalar inverted_total_area;

    Scalar b_first;
    Scalar b_second;

    // Compute Triangle coordinates.
    // Can we somehow use a referenced output: &output?
    std::pair<OutputIterator, bool> triangle_coordinates_2(const Point &query_point, OutputIterator output)
    {
        // Compute some related sub-areas.
        area_second = CGAL::area(triangle.vertex(1), triangle.vertex(2), query_point);
        area_third  = CGAL::area(triangle.vertex(2), triangle.vertex(0), query_point);

        // Compute the total inverted area of the triangle.
        inverted_total_area = Scalar(1) / CGAL::area(triangle.vertex(0), triangle.vertex(1), triangle.vertex(2));

        // Compute the first and second coordinate functions.
        b_first  = area_second * inverted_total_area;
        b_second = area_third  * inverted_total_area;

        *output = b_first;
        ++output;

        *output = b_second;
        ++output;

        // Compute the last - third coordinate using the partition of unity property.
        *output = Scalar(1) - b_first - b_second;

        // Output coordinates.
        return std::make_pair(output, true);
   }
};

// Class Triangle_coordinates_2 with particular std::back_insert_iterator instead of a general one.

/*!
 * \ingroup PkgBarycentric_coordinates_2
 * The class Tri_coordinates_2 implements barycentric coordinates ( <a href="http://mathworld.wolfram.com/BarycentricCoordinates.html" target=blanc>[1]</a>,
 * <a href="http://en.wikipedia.org/wiki/Barycentric_coordinate_system" target=blanc>[2]</a> ) with respect to an arbitrary non-degenerate triangle in the plane.
 * This class is parameterized by `CGAL::Triangle_2` class and a Container class.
 * The latter can be any class that fulfills the requirements for <a href="http://en.cppreference.com/w/cpp/iterator/back_insert_iterator">`std::back_insert_iterator`</a>. 
 * It defaults to <a href="http://en.cppreference.com/w/cpp/container/vector">`std::vector`</a> container.
 *
 * \sa <a href="http://en.cppreference.com/w/cpp/iterator/back_insert_iterator">`std::back_insert_iterator`</a>
 *
 */

template<typename Triangle_2, typename InputContainer = std::vector<typename Triangle_2::R::FT> > 
    class Tri_coordinates_2 : public Triangle_coordinates_2<Triangle_2, std::back_insert_iterator<InputContainer> >
{

public:

    // Creation.

    /// \name Types
    /// @{

    /// Type of the used container.
    typedef InputContainer Container;

    /// Type of the base class.
    typedef Triangle_coordinates_2<Triangle_2, std::back_insert_iterator<Container> > Base;

    /// @}

    /// \name Creation
    /// @{

    /// Creates an instance of Triangle_coordinates_2 class for a provided triangle passed as a reference.
    /// The used iterator is <a href="http://en.cppreference.com/w/cpp/iterator/back_insert_iterator">`std::back_insert_iterator`</a>.
    Tri_coordinates_2(const Triangle_2 &_triangle) : Base(_triangle) { }

    /// @}

    // Coordinates.

    /// \name Computation of the basis functions
    /// @{

    /// Computes Triangle barycentric coordinates for current query point with
    /// respect to all three vertices of the triangle.
    inline std::pair<std::back_insert_iterator<Container>, bool> compute(const typename Triangle_2::R::Point_2 &query_point, Container &container)
    {
        return Base::compute(query_point, std::back_inserter(container));
    }

    /// Static function, which gets a `CGAL::Triangle_2` and returns a `CGAL::Point_3` type of coordinates.
    static inline typename Triangle_2::R::Point_3 Compute(const Triangle_2 &_triangle, const typename Triangle_2::R::Point_2 &query_point)
    {
        return static_triangle_coordinates_2(_triangle, query_point);
    }

    /// @}

private:

    // Static function, which gets a CGAL::Triangle_2 and returns a CGAL::Point_3 with coordinates.
    static typename Triangle_2::R::Point_3 static_triangle_coordinates_2(const Triangle_2 &_triangle, const typename Triangle_2::R::Point_2 &query_point)
    {
        // Point type.
        typedef typename Triangle_2::R::Point_3 Point;

        // Triangle coordinates type.
        typedef typename CGAL::Barycentric_coordinates::Tri_coordinates_2<Triangle_2, Container> Triangle_coordinates;

        // Instantiate Triangle coordinates class for the triangle defined above.
        Triangle_coordinates triangle_coordinates(_triangle);

        // Create a container to store coordinates.
        Container coordinates;

        // Compute coordinates.
        triangle_coordinates.compute(query_point, coordinates);

        // Return CGAL::Point_3 type of coordinates.
        return Point(coordinates[0], coordinates[1], coordinates[2]);
    }
};

} // namespace Barycentric_coordinates

} // namespace CGAL

#endif // CGAL_TRIANGLE_COORDINATES_2_H