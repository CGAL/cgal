// Copyright (c) 2014 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is a part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Dmitry Anisimov, David Bommes, Kai Hormann, and Pierre Alliez.

/*!
  \file Discrete_harmonic_2.h
*/

#ifndef CGAL_DISCRETE_HARMONIC_2_H
#define CGAL_DISCRETE_HARMONIC_2_H

#include <CGAL/license/Barycentric_coordinates_2.h>

#include <CGAL/disable_warnings.h>

// CGAL headers.
#include <CGAL/assertions.h>
#include <CGAL/Polygon_2_algorithms.h>

// Barycentric coordinates headers.
#include <CGAL/Barycentric_coordinates_2/barycentric_enum_2.h>

// Boost headers.
#include <boost/optional/optional.hpp>

// CGAL namespace.
namespace CGAL {

// Barycentric coordinates namespace.
namespace Barycentric_coordinates {

// Examples: see the User Manual here - https://doc.cgal.org/latest/Manual/index.html.
// [1] Reference: "M. S. Floater, K. Hormann, and G. Kos. A general construction of barycentric coordinates over convex polygons. Advances in Computational Mathematics, 24(1-4):311-331, 2006.".

/*!
 * \ingroup PkgBarycentricCoordinates2Ref
 * The class `Discrete_harmonic_2` implements 2D discrete harmonic coordinates ( \cite cgal:bc:fhk-gcbcocp-06, \cite cgal:pp-cdmsc-93, \cite cgal:bc:eddhls-maam-95 ).
 * This class is parameterized by a traits class `Traits`, and it is used as a coordinate class to complete the class `Generalized_barycentric_coordinates_2`.
 * For a polygon with three vertices (triangle) it is better to use the class `Triangle_coordinates_2`.
 * Discrete harmonic coordinates can be computed exactly. By definition, they do not necesserily give positive values.

\tparam Traits must be a model of the concepts `BarycentricTraits_2` and `PolygonTraits_2`.

\cgalModels `BarycentricCoordinates_2`

\pre The provided polygon is strictly convex.

*/

template<class Traits>
    class Discrete_harmonic_2
{

public:

    /// \name Types
    /// @{

    /// Number type.
    typedef typename Traits::FT      FT;

    /// Point type.
    typedef typename Traits::Point_2 Point_2;

    /// @}

    // \name Creation

    // Creates the class `Discrete_harmonic_2` that implements the behaviour of discrete harmonic coordinates for any query point that does not belong to the polygon's boundary.
    // The polygon is given by a range of vertices of the type `Traits::Point_2` stored in a container of the type <a href="http://en.cppreference.com/w/cpp/container/vector">`std::vector`</a>.
    Discrete_harmonic_2(const std::vector<typename Traits::Point_2> &vertices, const Traits &b_traits) :
        vertex(vertices),
        barycentric_traits(b_traits),
        number_of_vertices(vertex.size()),
        area_2(barycentric_traits.compute_area_2_object()),
        squared_distance_2(barycentric_traits.compute_squared_distance_2_object()),
        collinear_2(barycentric_traits.collinear_2_object())
    {
        // Resize all the internal containers.
        r.resize(number_of_vertices);
        A.resize(number_of_vertices);
        B.resize(number_of_vertices);

        weight.resize(number_of_vertices);
    }

    // Computation of Discrete Harmonic Weight Functions

    // This function computes discrete harmonic weights (unnormalized coordinates) for a chosen query point.
    template<class OutputIterator>
        inline boost::optional<OutputIterator> weights(const Point_2 &query_point, OutputIterator &output)
    {
        return weights_2(query_point, output);
    }

    // Computation of Discrete Harmonic Basis Functions

    // This function computes discrete harmonic barycentric coordinates for a chosen query point on the bounded side of a strictly convex polygon.
    // \pre The provided polygon is strictly convex.
    template<class OutputIterator>
        inline boost::optional<OutputIterator> coordinates_on_bounded_side(const Point_2 &query_point, OutputIterator &output, const Type_of_algorithm type_of_algorithm)
    {
        switch(type_of_algorithm)
        {
            case PRECISE:
            return coordinates_on_bounded_side_precise_2(query_point, output);
            break;

            case FAST:
            return coordinates_on_bounded_side_fast_2(query_point, output);
            break;
        }

        // Pointer cannot be here. Something went wrong.
        const bool type_of_algorithm_failure = true;
        CGAL_postcondition( !type_of_algorithm_failure );
        if(!type_of_algorithm_failure) return boost::optional<OutputIterator>(output);
        else return boost::optional<OutputIterator>();
    }

    // This function computes discrete harmonic barycentric coordinates for a chosen query point on the unbounded side of a strictly convex polygon.
    // \pre The provided polygon is strictly convex.
    template<class OutputIterator>
        inline boost::optional<OutputIterator> coordinates_on_unbounded_side(const Point_2 &query_point, OutputIterator &output, const Type_of_algorithm type_of_algorithm, const bool warning_tag = true)
    {
        switch(type_of_algorithm)
        {
            case PRECISE:
            return coordinates_on_unbounded_side_precise_2(query_point, output, warning_tag);
            break;

            case FAST:
            return coordinates_on_unbounded_side_fast_2(query_point, output, warning_tag);
            break;
        }

        // Pointer cannot be here. Something went wrong.
        const bool type_of_algorithm_failure = true;
        CGAL_postcondition( !type_of_algorithm_failure );
        if(!type_of_algorithm_failure) return boost::optional<OutputIterator>(output);
        else return boost::optional<OutputIterator>();
    }

    // Information Functions

    // This function prints some information about discrete harmonic coordinates.
    void print_coordinates_information(std::ostream &output_stream) const
    {
        return print_coordinates_information_2(output_stream);
    }

private:

    // Some convenient typedefs.
    typedef typename std::vector<FT>      FT_vector;
    typedef typename std::vector<Point_2> Point_vector;

    // Internal global variables.
    const Point_vector &vertex;

    const Traits &barycentric_traits;

    const size_t number_of_vertices;

    FT_vector r, A, B, weight;

    FT dh_denominator, inverted_dh_denominator;

    typename Traits::Compute_area_2 area_2;
    typename Traits::Compute_squared_distance_2 squared_distance_2;
    typename Traits::Collinear_2 collinear_2;

    // WEIGHTS.

    // Compute discrete harmonic weights without normalization.
    template<class OutputIterator>
        boost::optional<OutputIterator> weights_2(const Point_2 &query_point, OutputIterator &output)
    {
        // Get the number of vertices in the polygon.
        const int n = int(number_of_vertices);

        // Compute areas A, B, and distances r following the notation from [1]. Split the loop to make this computation faster.
        r[0] = squared_distance_2(vertex[0] , query_point);
        A[0] = area_2(vertex[0]  , vertex[1], query_point);
        B[0] = area_2(vertex[n-1], vertex[1], query_point);

        for(int i = 1; i < n-1; ++i) {
            r[i] = squared_distance_2(vertex[i]   , query_point);
            A[i] = area_2(vertex[i]  , vertex[i+1], query_point);
            B[i] = area_2(vertex[i-1], vertex[i+1], query_point);
        }

        r[n-1] = squared_distance_2(vertex[n-1], query_point);
        A[n-1] = area_2(vertex[n-1], vertex[0] , query_point);
        B[n-1] = area_2(vertex[n-2], vertex[0] , query_point);

        // Compute unnormalized weights following the formula (25) with p = 2 from [1].
        CGAL_precondition( A[n-1] != FT(0) && A[0] != FT(0) );
        *output = (r[1]*A[n-1] - r[0]*B[0] + r[n-1]*A[0]) / (A[n-1] * A[0]);
        ++output;

        for(int i = 1; i < n-1; ++i) {
            CGAL_precondition( A[i-1] != FT(0) && A[i] != FT(0) );
            *output = (r[i+1]*A[i-1] - r[i]*B[i] + r[i-1]*A[i]) / (A[i-1] * A[i]);
            ++output;
        }

        CGAL_precondition( A[n-2] != FT(0) && A[n-1] != FT(0) );
        *output = (r[0]*A[n-2] - r[n-1]*B[n-1] + r[n-2]*A[n-1]) / (A[n-2] * A[n-1]);

        // Return weights.
        return boost::optional<OutputIterator>(output);
    }

    // COORDINATES ON BOUNDED SIDE.

    // Compute discrete harmonic coordinates on the bounded side of the polygon with the slow O(n^2) but precise algorithm.
    // Here, n - is the number of the polygon's vertices.
    template<class OutputIterator>
        boost::optional<OutputIterator> coordinates_on_bounded_side_precise_2(const Point_2 &query_point, OutputIterator &output)
    {
        CGAL_precondition( type_of_polygon() == STRICTLY_CONVEX );

        // Get the number of vertices in the polygon.
        const int n = int(number_of_vertices);

        // Compute areas A, B, and distances r following the notation from [1]. Split the loop to make this computation faster.
        r[0] = squared_distance_2(vertex[0] , query_point);
        A[0] = area_2(vertex[0]  , vertex[1], query_point);
        B[0] = area_2(vertex[n-1], vertex[1], query_point);

        for(int i = 1; i < n-1; ++i) {
            r[i] = squared_distance_2(vertex[i]   , query_point);
            A[i] = area_2(vertex[i]  , vertex[i+1], query_point);
            B[i] = area_2(vertex[i-1], vertex[i+1], query_point);
        }

        r[n-1] = squared_distance_2(vertex[n-1], query_point);
        A[n-1] = area_2(vertex[n-1], vertex[0] , query_point);
        B[n-1] = area_2(vertex[n-2], vertex[0] , query_point);

        // Initialize weights with the numerator of the formula (25) with p = 2 from [1].
        // Then we multiply them by areas A as in the formula (5) in [1]. We also split the loop.
        weight[0] = r[1]*A[n-1] - r[0]*B[0] + r[n-1]*A[0];
        for(int j = 1; j < n-1; ++j) weight[0] *= A[j];

        for(int i = 1; i < n-1; ++i) {
            weight[i] = r[i+1]*A[i-1] - r[i]*B[i] + r[i-1]*A[i];
            for(int j = 0; j < i-1; ++j) weight[i] *= A[j];
            for(int j = i+1; j < n; ++j) weight[i] *= A[j];
        }

        weight[n-1] = r[0]*A[n-2] - r[n-1]*B[n-1] + r[n-2]*A[n-1];
        for(int j = 0; j < n-2; ++j) weight[n-1] *= A[j];

        // Compute the sum of all weights - denominator of discrete harmonic coordinates.
        dh_denominator = weight[0];
        for(int i = 1; i < n; ++i) dh_denominator += weight[i];

        // Invert this denominator.
        CGAL_precondition( dh_denominator != FT(0) );
        inverted_dh_denominator = FT(1) / dh_denominator;

        // Normalize weights and save them as resulting discrete harmonic coordinates.
        for(int i = 0; i < n-1; ++i) {
            *output = weight[i] * inverted_dh_denominator;
            ++output;
        }
        *output = weight[n-1] * inverted_dh_denominator;

        // Return coordinates.
        return boost::optional<OutputIterator>(output);
    }

    // Compute discrete harmonic coordinates on the bounded side of the polygon with the fast O(n) but less precise algorithm.
    // Here, n - is the number of the polygon's vertices. Precision is lost near the boundary (~ 1.0e-10 and closer).
    template<class OutputIterator>
        boost::optional<OutputIterator> coordinates_on_bounded_side_fast_2(const Point_2 &query_point, OutputIterator &output)
    {
        CGAL_precondition( type_of_polygon() == STRICTLY_CONVEX );

        // Get the number of vertices in the polygon.
        const int n = int(number_of_vertices);

        // Compute areas A, B, and distances r following the notation from [1]. Split the loop to make this computation faster.
        r[0] = squared_distance_2(vertex[0] , query_point);
        A[0] = area_2(vertex[0]  , vertex[1], query_point);
        B[0] = area_2(vertex[n-1], vertex[1], query_point);

        for(int i = 1; i < n-1; ++i) {
            r[i] = squared_distance_2(vertex[i]   , query_point);
            A[i] = area_2(vertex[i]  , vertex[i+1], query_point);
            B[i] = area_2(vertex[i-1], vertex[i+1], query_point);
        }

        r[n-1] = squared_distance_2(vertex[n-1], query_point);
        A[n-1] = area_2(vertex[n-1], vertex[0] , query_point);
        B[n-1] = area_2(vertex[n-2], vertex[0] , query_point);

        // Compute unnormalized weights following the formula (25) with p = 2 from [1].
        CGAL_precondition( A[n-1] != FT(0) && A[0] != FT(0) );
        weight[0] = (r[1]*A[n-1] - r[0]*B[0] + r[n-1]*A[0]) / (A[n-1] * A[0]);

        for(int i = 1; i < n-1; ++i) {
            CGAL_precondition( A[i-1] != FT(0) && A[i] != FT(0) );
            weight[i] = (r[i+1]*A[i-1] - r[i]*B[i] + r[i-1]*A[i]) / (A[i-1] * A[i]);
        }

        CGAL_precondition( A[n-2] != FT(0) && A[n-1] != FT(0) );
        weight[n-1] = (r[0]*A[n-2] - r[n-1]*B[n-1] + r[n-2]*A[n-1]) / (A[n-2] * A[n-1]);

        // Compute the sum of all weights - denominator of discrete harmonic coordinates.
        dh_denominator = weight[0];
        for(int i = 1; i < n; ++i) dh_denominator += weight[i];

        // Invert this denominator.
        CGAL_precondition( dh_denominator != FT(0) );
        inverted_dh_denominator = FT(1) / dh_denominator;

        // Normalize weights and save them as resulting discrete harmonic coordinates.
        for(int i = 0; i < n-1; ++i) {
            *output = weight[i] * inverted_dh_denominator;
            ++output;
        }
        *output = weight[n-1] * inverted_dh_denominator;

        // Return coordinates.
        return boost::optional<OutputIterator>(output);
    }

    // COORDINATES ON UNBOUNDED SIDE.

    // Compute discrete harmonic coordinates on the unbounded side of the polygon with the slow O(n^2) but precise algorithm.
    // Here, n - is the number of the polygon's vertices.
    template<class OutputIterator>
        boost::optional<OutputIterator> coordinates_on_unbounded_side_precise_2(const Point_2 &query_point, OutputIterator &output, bool warning_tag)
    {
        if(warning_tag)
            std::cout << std::endl << "ATTENTION: Discrete harmonic coordinates might be not well-defined outside the polygon!" << std::endl;

        // Use the same formulas as for the bounded side since they are also valid on the unbounded side.
        return coordinates_on_bounded_side_precise_2(query_point, output);
    }

    // Compute discrete harmonic coordinates on the unbounded side of the polygon with the fast O(n) but less precise algorithm.
    // Here, n - is the number of the polygon's vertices. Precision is lost near the boundary (~ 1.0e-10 and closer).
    template<class OutputIterator>
        boost::optional<OutputIterator> coordinates_on_unbounded_side_fast_2(const Point_2 &query_point, OutputIterator &output, bool warning_tag)
    {
        if(warning_tag)
            std::cout << std::endl << "ATTENTION: Discrete harmonic coordinates might be not well-defined outside the polygon!" << std::endl;

        // Use the same formulas as for the bounded side since they are also valid on the unbounded side.
        return coordinates_on_bounded_side_fast_2(query_point, output);
    }

    // OTHER FUNCTIONS.

    // Print some information about discrete harmonic coordinates.
    void print_coordinates_information_2(std::ostream &output_stream) const
    {
        output_stream << std::endl << "CONVEXITY: " << std::endl << std::endl;

        if(type_of_polygon() == STRICTLY_CONVEX) {
            output_stream << "This polygon is strictly convex." << std::endl;
        } else if(type_of_polygon() == WEAKLY_CONVEX) {
            output_stream << "This polygon is weakly convex. The correct computation is not expected!" << std::endl;
        } else if(type_of_polygon() == CONCAVE) {
            output_stream << "This polygon polygon is not convex. The correct computation is not expected!" << std::endl;
        }

        output_stream << std::endl << "TYPE OF COORDINATES: " << std::endl << std::endl;
        output_stream << "The coordinate functions to be computed are discrete harmonic coordinates." << std::endl;

        output_stream << std::endl << "INFORMATION ABOUT COORDINATES: " << std::endl << std::endl;
        output_stream << "Discrete harmonic coordinates are well-defined in the closure of an arbitrary strictly convex polygon and can be computed exactly." << std::endl;

        output_stream << std::endl;
        output_stream << "They satisfy the following properties: " << std::endl;
        output_stream << "1. Partition of unity or constant precision;" << std::endl;
        output_stream << "2. Homogeneity or linear precision;" << std::endl;
        output_stream << "3. Lagrange property;" << std::endl;
        output_stream << "4. Linearity along edges;" << std::endl;
        output_stream << "5. Smoothness;" << std::endl;
        output_stream << "6. Similarity invariance;" << std::endl;

        output_stream << std::endl;
        output_stream << "For polygons whose vertices lie on a common circle, they coincide with Wachspress coordinates." << std::endl;

        output_stream << std::endl << "REFERENCE: " << std::endl << std::endl;
        output_stream << "M. S. Floater, K. Hormann, and G. Kos. A general construction of barycentric coordinates over convex polygons. Advances in Computational Mathematics, 24(1-4):311-331, 2006." << std::endl;
    }

    // Check the type of the provided polygon - CONVEX, STRICTLY_CONVEX, or CONCAVE.
    Type_of_polygon type_of_polygon() const
    {
        // First, test the polygon on convexity.
        if(CGAL::is_convex_2(vertex.begin(), vertex.end(), barycentric_traits)) {

            // Index of the last polygon's vertex.
            const int last = int(number_of_vertices) - 1;

            // Test all the consequent triplets of the polygon's vertices on collinearity.
            // In case we find at least one, return WEAKLY_CONVEX polygon.
            if(collinear_2(vertex[last], vertex[0], vertex[1]))
                return WEAKLY_CONVEX;

            for(int i = 1; i < last; ++i) {
                if(collinear_2(vertex[i-1], vertex[i], vertex[i+1]))
                    return WEAKLY_CONVEX;
            }

            if(collinear_2(vertex[last-1], vertex[last], vertex[0]))
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

#include <CGAL/enable_warnings.h>

#endif // CGAL_DISCRETE_HARMONIC_2_H
