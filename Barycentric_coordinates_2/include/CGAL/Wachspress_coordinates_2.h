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
  \file Wachspress_coordinates_2.h
*/

#ifndef CGAL_WACHSPRESS_COORDINATES_2_H
#define CGAL_WACHSPRESS_COORDINATES_2_H

// CGAL headers.
#include <CGAL/assertions.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Kernel/global_functions_2.h> 

// Barycentric coordinates headers.
#include <CGAL/barycentric_enum.h>

// CGAL namespace.
namespace CGAL {

// Barycentric coordinates namespace.
namespace Barycentric_coordinates {

// Examples: see the User Manual here - http://doc.cgal.org/latest/Manual/index.html.
// [1] Reference: "M. S. Floater, K. Hormann, and G. Kos. A general construction of barycentric coordinates over convex polygons. Advances in Computational Mathematics, 24(1-4):311-331, 2006.".

/*!
 * \ingroup PkgBarycentric_coordinates_2
 * The class Wachspress_coordinates_2 implements 2D Wachspress coordinates ( \cite cgal:bc:fhk-gcbcocp-06, \cite cgal:bc:mlbd-gbcip-02, \cite cgal:bc:w-rfeb-75 ).
 * This class is parameterized by a traits class `Traits`, and it is used as a coordinate class to complete the base class `Barycentric_coordinates_2`.
 * For a polygon with three vertices it is better to use the class `CGAL::Barycentric_coordinates::Triangle_coordinates_2`.
 * Wachspress coordinates can be computed exactly, and they are always positive in the closure of a strictly convex polygon.
 *
 * \pre The polygon's vertices must be ordered.
 
 \cgalHeading{Requirements}

 <OL>
 <LI> `Traits` class must contain the following subset of types:
 <UL>
 <LI> `Traits::Point_2` - the type of a point used internally in the class, which is equivalent to the type `CGAL::Point_2`;
 <LI> `Traits::K` - the used kernel;
 <LI> `Traits::FT` - the type of a coordinate value;
 </UL>
 </OL>

 */
 
template<class Traits> 
    class Wachspress_coordinates_2
{

public:

    /// \name Types
    /// @{

    /// Number type.
    typedef typename Traits::FT      Scalar;

    /// Type of the used kernel.
    typedef typename Traits::K       Kernel;

    /// Type of 2D point.
    typedef typename Traits::Point_2 Point_2;

    /// @}

    /// \name Creation
    /// @{

    /// Creates the class `Wachspress_coordinates_2` that implements the behaviour of Wachspress coordinates for any query point that does not belong to the polygon's boundary.
    /// The polygon is given by a range of vertices of the type `CGAL::Point_2` stored in a container of the type <a href="http://en.cppreference.com/w/cpp/container/vector">`std::vector`</a>.
    Wachspress_coordinates_2(const std::vector<Point_2> &vertices) :
        vertex(vertices),
        number_of_vertices(vertex.size())
    {
        // Resize all the internal containers.
        A.resize(number_of_vertices);
        C.resize(number_of_vertices);

        weight.resize(number_of_vertices);
    }

    /// @}

    /// \name Computation of Wachspress weight functions
    /// @{

    /// This function is intended to compute Wachspress weights for a chosen query point.
    template<typename Iterator>
        inline std::pair<Iterator, bool> weights(const Point_2 &query_point, Iterator &output)
    {
        return weights_2(query_point, output);
    }

    /// @}

    /// \name Computation of Wachspress basis functions
    /// @{

    /// This function is intended to compute Wachspress barycentric coordinates for a chosen query point on the bounded side of a strictly convex polygon with the O(n^2) precise algorithm.
    /// \pre The provided polygon is strictly convex that is it complies with the constant `CGAL::Barycentric_coordinates::STRICTLY_CONVEX`. 
    template<typename Iterator>
        inline std::pair<Iterator, bool> coordinates_on_bounded_side_precise(const Point_2 &query_point, Iterator &output)
    {   
        return coordinates_on_bounded_side_precise_2(query_point, output);
    }

    /// This function is intended to compute Wachspress barycentric coordinates for a chosen query point on the bounded side of a strictly convex polygon with the O(n) fast algorithm.
    /// \pre The provided polygon is strictly convex that is it complies with the constant `CGAL::Barycentric_coordinates::STRICTLY_CONVEX`. 
    template<typename Iterator>
        inline std::pair<Iterator, bool> coordinates_on_bounded_side_fast(const Point_2 &query_point, Iterator &output)
    {   
        return coordinates_on_bounded_side_fast_2(query_point, output);
    }

    /// This function is intended to compute Wachspress barycentric coordinates for a chosen query point on the unbounded side of a strictly convex polygon with the O(n^2) precise algorithm.
    /// \pre The provided polygon is strictly convex that is it complies with the constant `CGAL::Barycentric_coordinates::STRICTLY_CONVEX`. 
    template<typename Iterator>
        inline std::pair<Iterator, bool> coordinates_on_unbounded_side_precise(const Point_2 &query_point, Iterator &output)
    {   
        return coordinates_on_unbounded_side_precise_2(query_point, output);
    }

    /// This function is intended to compute Wachspress barycentric coordinates for a chosen query point on the unbounded side of a strictly convex polygon with the O(n) fast algorithm.
    /// \pre The provided polygon is strictly convex that is it complies with the constant `CGAL::Barycentric_coordinates::STRICTLY_CONVEX`. 
    template<typename Iterator>
        inline std::pair<Iterator, bool> coordinates_on_unbounded_side_fast(const Point_2 &query_point, Iterator &output)
    {   
        return coordinates_on_unbounded_side_fast_2(query_point, output);
    }

    /// @}

    /// \name Information functions
    /// @{

    /// Print some information about 2D Wachspress coordinates.
    void print_coordinates_information(std::ostream &output_stream) const
    {
        return print_coordinates_information_2(output_stream);
    }

    /// @}

private:

    // Some convenient typedefs.
    typedef typename std::vector<Scalar> Scalar_vector;
    typedef typename std::vector<Point_2> Point_vector;

    // Internal global variables.
    const Point_vector &vertex;

    const size_t number_of_vertices;

    Scalar_vector A, C, weight;

    Scalar wp_denominator, inverted_wp_denominator;

    // WEIGHTS.

    // Compute 2D Wachspress weights without normalization.
    template<typename Iterator> 
        std::pair<Iterator, bool> weights_2(const Point_2 &query_point, Iterator &output)
    {
        // Get the number of vertices in the polygon.
        const int n = int(number_of_vertices);

        // Compute areas A and C following the area notation from [1]. Split the loop to make this computation faster.
        A[0] = CGAL::area(vertex[0]  , vertex[1], query_point);
        C[0] = CGAL::area(vertex[n-1], vertex[0], vertex[1]  );

        for(int i = 1; i < n-1; ++i) {
            A[i] = CGAL::area(vertex[i]  , vertex[i+1], query_point);
            C[i] = CGAL::area(vertex[i-1], vertex[i]  , vertex[i+1]);
        }

        A[n-1] = CGAL::area(vertex[n-1], vertex[0]  , query_point);
        C[n-1] = CGAL::area(vertex[n-2], vertex[n-1], vertex[0]  );

        // Compute unnormalized weights following the formula (28) from [1].
        CGAL_precondition( A[n-1] != Scalar(0) && A[0] != Scalar(0) );
        *output = C[0] / (A[n-1] * A[0]);
        ++output;

        for(int i = 1; i < n-1; ++i) {
            CGAL_precondition( A[i-1] != Scalar(0) && A[i] != Scalar(0) );
            *output = C[i] / (A[i-1] * A[i]);
            ++output;
        }

        CGAL_precondition( A[n-2] != Scalar(0) && A[n-1] != Scalar(0) );
        *output = C[n-1] / (A[n-2] * A[n-1]);

        // Return weights.
        return std::make_pair(output, true);
    }

    // COORDINATES ON BOUNDED SIDE.

    // Compute 2D Wachspress coordinates on the bounded side of the polygon with the slow O(n^2) but precise algorithm.
    // Here, n - is the number of the polygon's vertices.
    template<typename Iterator>
        std::pair<Iterator, bool> coordinates_on_bounded_side_precise_2(const Point_2 &query_point, Iterator &output)
    {
        CGAL_precondition( type_of_polygon() == STRICTLY_CONVEX );

        // Get the number of vertices in the polygon.
        const int n = int(number_of_vertices);

        // Compute areas A following the area notation from [1]. Split the loop to make this computation faster.
        A[0] = CGAL::area(vertex[0], vertex[1], query_point);
        for(int i = 1; i < n-1; ++i) A[i] = CGAL::area(vertex[i], vertex[i+1], query_point);
        A[n-1] = CGAL::area(vertex[n-1], vertex[0], query_point);

        // Initialize weights with areas C following the area notation from [1].
        // Then we multiply them by areas A as in the formula (5) from [1]. We also split the loop.
        weight[0] = CGAL::area(vertex[n-1], vertex[0], vertex[1]);
        for(int j = 1; j < n-1; ++j) weight[0] *= A[j];

        for(int i = 1; i < n-1; ++i) {
            weight[i] = CGAL::area(vertex[i-1], vertex[i], vertex[i+1]);
            for(int j = 0; j < i-1; ++j) weight[i] *= A[j];
            for(int j = i+1; j < n; ++j) weight[i] *= A[j];
        }

        weight[n-1] = CGAL::area(vertex[n-2], vertex[n-1], vertex[0]);
        for(int j = 0; j < n-2; ++j) weight[n-1] *= A[j];

        // Compute the sum of all weights - denominator of Wachspress coordinates.
        wp_denominator = weight[0];
        for(int i = 1; i < n; ++i) wp_denominator += weight[i];

        // Invert this denominator.
        CGAL_precondition( wp_denominator != Scalar(0) );
        inverted_wp_denominator = Scalar(1) / wp_denominator;

        // Normalize weights and save them as resulting Wachspress coordinates.
        for(int i = 0; i < n-1; ++i) {
            *output = weight[i] * inverted_wp_denominator;
            ++output;
        }
        *output = weight[n-1] * inverted_wp_denominator;

        // Return coordinates.
        return std::make_pair(output, true);
    }

    // Compute 2D Wachspress coordinates on the bounded side of the polygon with the fast O(n) but less precise algorithm.
    // Here, n - is the number of the polygon's vertices. Precision is lost near the boundary (~ 1.0e-10 and closer).
    template<typename Iterator>
        std::pair<Iterator, bool> coordinates_on_bounded_side_fast_2(const Point_2 &query_point, Iterator &output)
    {
        CGAL_precondition( type_of_polygon() == STRICTLY_CONVEX );

        // Get the number of vertices in the polygon.
        const int n = int(number_of_vertices);

        // Compute areas A and C following the area notation from [1]. Split the loop to make this computation faster.
        A[0] = CGAL::area(vertex[0]  , vertex[1], query_point);
        C[0] = CGAL::area(vertex[n-1], vertex[0], vertex[1]  );

        for(int i = 1; i < n-1; ++i) {
            A[i] = CGAL::area(vertex[i]  , vertex[i+1], query_point);
            C[i] = CGAL::area(vertex[i-1], vertex[i]  , vertex[i+1]);
        }

        A[n-1] = CGAL::area(vertex[n-1], vertex[0]  , query_point);
        C[n-1] = CGAL::area(vertex[n-2], vertex[n-1], vertex[0]  );

        // Compute the unnormalized weights following the formula (28) from [1].
        CGAL_precondition( A[n-1] != Scalar(0) && A[0] != Scalar(0) );
        weight[0] = C[0] / (A[n-1] * A[0]);

        for(int i = 1; i < n-1; ++i) {
            CGAL_precondition( A[i-1] != Scalar(0) && A[i] != Scalar(0) );
            weight[i] = C[i] / (A[i-1] * A[i]);
        }

        CGAL_precondition( A[n-2] != Scalar(0) && A[n-1] != Scalar(0) );
        weight[n-1] = C[n-1] / (A[n-2] * A[n-1]);

        // Compute the sum of all weights - denominator of Wachspress coordinates.
        wp_denominator = weight[0];
        for(int i = 1; i < n; ++i) wp_denominator += weight[i];

        // Invert this denominator.
        CGAL_precondition( wp_denominator != Scalar(0) );
        inverted_wp_denominator = Scalar(1) / wp_denominator;

        // Normalize weights and save them as resulting Wachspress coordinates.
        for(int i = 0; i < n-1; ++i) {
            *output = weight[i] * inverted_wp_denominator;
            ++output;
        }
        *output = weight[n-1] * inverted_wp_denominator;

        // Return coordinates.
        return std::make_pair(output, true);
    }

    // COORDINATES ON UNBOUNDED SIDE.

    // Compute 2D Wachspress coordinates on the unbounded side of the polygon with the slow O(n^2) but precise algorithm.
    // Here, n - is the number of the polygon's vertices.
    template<typename Iterator>
        std::pair<Iterator, bool> coordinates_on_unbounded_side_precise_2(const Point_2 &query_point, Iterator &output)
    {
        std::cout << std::endl << "WARNING: Wachspress coordinates might be not well-defined outside the polygon!" << std::endl;

        // Use the same formulas as for the bounded side since they are also valid on the unbounded side.
        return coordinates_on_bounded_side_precise_2(query_point, output);
    }

    // Compute 2D Wachspress coordinates on the unbounded side of the polygon with the fast O(n) but less precise algorithm.
    // Here, n - is the number of the polygon's vertices. Precision is lost near the boundary (~ 1.0e-10 and closer).
    template<typename Iterator>
        std::pair<Iterator, bool> coordinates_on_unbounded_side_fast_2(const Point_2 &query_point, Iterator &output)
    {
        std::cout << std::endl << "WARNING: Wachspress coordinates might be not well-defined outside the polygon!" << std::endl;

        // Use the same formulas as for the bounded side since they are also valid on the unbounded side.
        return coordinates_on_bounded_side_fast_2(query_point, output);
    }

    // OTHER FUNCTIONS.

    // Print some information about 2D Wachspress coordinates.
    void print_coordinates_information_2(std::ostream &output_stream) const
    {
        output_stream << std::endl << "CONVEXITY: " << std::endl << std::endl;

        if(type_of_polygon() == CONCAVE)         output_stream << "This polygon is not convex. The correct computation is not expected!" << std::endl;
        if(type_of_polygon() == WEAKLY_CONVEX)   output_stream << "This polygon is weakly convex. The correct computation is not expected!" << std::endl;
        if(type_of_polygon() == STRICTLY_CONVEX) output_stream << "This polygon is strictly convex." << std::endl;

        output_stream << std::endl << "TYPE OF COORDINATES: " << std::endl << std::endl;
        output_stream << "The coordinate functions to compute are Wachspress coordinates." << std::endl;

        output_stream << std::endl << "INFORMATION ABOUT COORDINATES: " << std::endl << std::endl;
        output_stream << "Wachspress coordinates are well-defined in the closure of an arbitrary strictly convex polygon and can be computed exactly." << std::endl;

        output_stream << std::endl;
        output_stream << "They satisfy the following properties: " << std::endl;
        output_stream << "1. Partition of unity or constant precision;" << std::endl;
        output_stream << "2. Homogeneity or linear precision;" << std::endl;
        output_stream << "3. Lagrange property;" << std::endl;
        output_stream << "4. Non-negativity;" << std::endl;
        output_stream << "5. Boundedness between 0 and 1;" << std::endl;
        output_stream << "6. Linearity along edges;" << std::endl;
        output_stream << "7. Smoothness;" << std::endl;
        output_stream << "8. Similarity invariance;" << std::endl;

        output_stream << std::endl;
        output_stream << "For polygons whose vertices lie on a common circle, they coincide with Discrete Harmonic coordinates." << std::endl;

        output_stream << std::endl << "REFERENCE: " << std::endl << std::endl;
        output_stream << "M. S. Floater, K. Hormann, and G. Kos. A general construction of barycentric coordinates over convex polygons. Advances in Computational Mathematics, 24(1-4):311-331, 2006." << std::endl;
    }

    // Check the type of the provided polygon - CONVEX, STRICTLY_CONVEX, or CONCAVE.
    Type_of_polygon type_of_polygon() const
    {
        // First, test the polygon on convexity.
        if(CGAL::is_convex_2(vertex.begin(), vertex.end(), Kernel())) {

            // Index of the last polygon's vertex.
            const int last = int(number_of_vertices) - 1;

            // Test all the consequent triplets of the polygon's vertices on collinearity. 
            // In case we find at least one, return WEAKLY_CONVEX polygon.
            if(CGAL::collinear(vertex[last], vertex[0], vertex[1]))
                return WEAKLY_CONVEX;

            for(int i = 1; i < last; ++i) {
                if(CGAL::collinear(vertex[i-1], vertex[i], vertex[i+1]))
                    return WEAKLY_CONVEX;
            }

            if(CGAL::collinear(vertex[last-1], vertex[last], vertex[0]))
                return WEAKLY_CONVEX;

            // Otherwise, return STRICTLY_CONVEX polygon.
            return STRICTLY_CONVEX;
        }

        // Otherwise, return CONCAVE polygon.
        return CONCAVE;
    }
};

} // namespace Barycentric_coordinates

} // namespace CGAL

#endif // CGAL_WACHSPRESS_COORDINATES_2_H