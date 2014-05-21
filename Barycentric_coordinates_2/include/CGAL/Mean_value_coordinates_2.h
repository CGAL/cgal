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
  \file Mean_value_coordinates_2.h
*/

#ifndef CGAL_MEAN_VALUE_COORDINATES_2_H
#define CGAL_MEAN_VALUE_COORDINATES_2_H

// CGAL headers.
#include <CGAL/Vector_2.h>
#include <CGAL/assertions.h>
#include <CGAL/Kernel/global_functions_2.h>

// CGAL namespace.
namespace CGAL {

// Barycentric coordinates namespace.
namespace Barycentric_coordinates {

// Examples: see the User Manual here - http://doc.cgal.org/latest/Manual/index.html.
// [1] Reference: "K. Hormann and M. Floater. Mean value coordinates for arbitrary planar polygons. ACM Transactions on Graphics, 25(4):1424-1441, 2006.".
// [2] Reference: "M. S. Floater, Wachspress and mean value coordinates, to appear in the Proceedings of the 14th International Conference on Approximation Theory, G. Fasshauer and L. L. Schumaker (eds.)."

/*!
 * \ingroup PkgBarycentric_coordinates_2
 * The class Mean_value_coordinates_2 implements 2D mean value coordinates ( \cite cgal:bc:hf-mvcapp-06, \cite cgal:bc:fhk-gcbcocp-06, \cite cgal:f-mvc-03 ).
 * This class is parameterized by a traits class `Traits`, and it is used as a coordinate class to complete the base class `Barycentric_coordinates_2`.
 * For a polygon with three vertices (triangle) it is better to use the class `CGAL::Barycentric_coordinates::Triangle_coordinates_2`.
 * Mean value coordinates can be computed only approximately due to an inevitable square root operation, and they are necesserily positive only inside the kernel of a star-shaped polygon and inside any quadrilateral.
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
    class Mean_value_coordinates_2
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

    /// Creates the class `Mean_value_coordinates_2` that implements the behaviour of mean value coordinates for any query point that does not belong to the polygon's boundary.
    /// The polygon is given by a range of vertices of the type `CGAL::Point_2` stored in a container of the type <a href="http://en.cppreference.com/w/cpp/container/vector">`std::vector`</a>.
    Mean_value_coordinates_2(const std::vector<Point_2> &vertices) :
        vertex(vertices),
        number_of_vertices(vertex.size())
    {
        // Resize all the internal containers.
        s.resize(number_of_vertices);

        r.resize(number_of_vertices);
        A.resize(number_of_vertices);
        B.resize(number_of_vertices);
        D.resize(number_of_vertices);
        P.resize(number_of_vertices);
        t.resize(number_of_vertices);

        weight.resize(number_of_vertices);
    }

    /// @}

    /// \name Computation of mean value weight functions
    /// @{

    /// This function is intended to compute mean value weights for a chosen query point.
    template<class Iterator>
        inline std::pair<Iterator, bool> weights(const Point_2 &query_point, Iterator &output)
    {
        return weights_2(query_point, output);
    }

    /// @}

    /// \name Computation of mean value basis functions
    /// @{

    /// This function is intended to compute mean value barycentric coordinates for a chosen query point on the bounded side of a simple polygon with the O(n^2) precise algorithm.
    template<class Iterator>
        inline std::pair<Iterator, bool> coordinates_on_bounded_side_precise(const Point_2 &query_point, Iterator &output)
    {   
        return coordinates_on_bounded_side_precise_2(query_point, output);
    }

    /// This function is intended to compute mean value barycentric coordinates for a chosen query point on the bounded side of a simple polygon with the O(n) fast algorithm.
    template<class Iterator>
        inline std::pair<Iterator, bool> coordinates_on_bounded_side_fast(const Point_2 &query_point, Iterator &output)
    {   
        return coordinates_on_bounded_side_fast_2(query_point, output);
    }

    /// This function is intended to compute mean value barycentric coordinates for a chosen query point on the unbounded side of a simple polygon with the O(n^2) precise algorithm.
    template<class Iterator>
        inline std::pair<Iterator, bool> coordinates_on_unbounded_side_precise(const Point_2 &query_point, Iterator &output)
    {   
        return coordinates_on_unbounded_side_precise_2(query_point, output);
    }

    /// This function is intended to compute mean value barycentric coordinates for a chosen query point on the unbounded side of a simple polygon with the O(n) fast algorithm.
    template<class Iterator>
        inline std::pair<Iterator, bool> coordinates_on_unbounded_side_fast(const Point_2 &query_point, Iterator &output)
    {   
        return coordinates_on_unbounded_side_fast_2(query_point, output);
    }

    /// @}

    /// \name Information functions
    /// @{

    /// Print some information about 2D mean value coordinates.
    void print_coordinates_information(std::ostream &output_stream) const
    {
        return print_coordinates_information_2(output_stream);
    }

    /// @}

private:

    // Some convenient typedefs.
    typedef typename std::vector<Scalar> Scalar_vector;
    typedef typename std::vector<Point_2> Point_vector;
    typedef Vector_2<Kernel>                    Vector;
    typedef typename std::vector<Vector> Vector_vector;

    // Internal global variables.
    const Point_vector &vertex;

    const size_t number_of_vertices;

    Vector_vector s;

    Scalar_vector r, A, B, D, P, t, weight;

    Scalar mv_denominator, inverted_mv_denominator;

    // WEIGHTS.

    // Compute 2D mean value weights without normalization.
    template<class Iterator>
        std::pair<Iterator, bool> weights_2(const Point_2 &query_point, Iterator &output)
    {
        // Get the number of vertices in the polygon.
        const int n = int(number_of_vertices);

        // Compute vector s following the pseudo-code in the Figure 10 from [1].
        for(int i = 0; i < n; ++i) s[i] = vertex[i] - query_point;

        // Compute length r, area A, and dot product D following the pseudo-code in the Figure 10 from [1].
        // Split the loop to make this computation faster.
        r[0] = Scalar(CGAL::sqrt(CGAL::to_double(s[0].squared_length())));
        A[0] = CGAL::area(vertex[0], vertex[1], query_point);
        D[0] = s[0]*s[1];

        for(int i = 1; i < n-1; ++i) {
            r[i] = Scalar(CGAL::sqrt(CGAL::to_double(s[i].squared_length())));
            A[i] = CGAL::area(vertex[i], vertex[i+1], query_point);
            D[i] = s[i]*s[i+1];
        }

        r[n-1] = Scalar(CGAL::sqrt(CGAL::to_double(s[n-1].squared_length())));
        A[n-1] = CGAL::area(vertex[n-1], vertex[0], query_point);
        D[n-1] = s[n-1]*s[0];

        // Compute intermediate values t using the formulas from slide 19 here
        // - http://www.inf.usi.ch/hormann/nsfworkshop/presentations/Hormann.pdf 
        for(int i = 0; i < n-1; ++i) {
            CGAL_precondition( (r[i]*r[i+1] + D[i]) != Scalar(0) );
            t[i] = A[i] / (r[i]*r[i+1] + D[i]);
        }

        CGAL_precondition( (r[n-1]*r[0] + D[n-1]) != Scalar(0) );
        t[n-1] = A[n-1] / (r[n-1]*r[0] + D[n-1]);

        // Compute mean value weights using the same pseudo-code as before.
        CGAL_precondition( r[0] != Scalar(0) );
        *output = (t[n-1] + t[0]) / r[0];
        ++output;

        for(int i = 1; i < n-1; ++i) {
            CGAL_precondition( r[i] != Scalar(0) );
            *output = (t[i-1] + t[i]) / r[i];
            ++output;
        }

        CGAL_precondition( r[n-1] != Scalar(0) );
        *output = (t[n-2] + t[n-1]) / r[n-1];

        // Return weights.
        return std::make_pair(output, true);
    }

    // COORDINATES ON BOUNDED SIDE.

    // Compute 2D mean value coordinates on the bounded side of the polygon with the slow O(n^2) but precise algorithm.
    // Here, n - is the number of the polygon's vertices.
    template<class Iterator>
        std::pair<Iterator, bool> coordinates_on_bounded_side_precise_2(const Point_2 &query_point, Iterator &output)
    {
        // Get the number of vertices in the polygon.
        const int n = int(number_of_vertices);

        // Compute vector s and its length r following the pseudo-code in the Figure 10 from [1].
        s[0] = vertex[0] - query_point;
        r[0] = Scalar(CGAL::sqrt(CGAL::to_double(s[0].squared_length())));

        // Compute areas A and B following notations from [1] (see Figure 2). Split the loop to make this computation faster.
        A[0] = CGAL::area(vertex[0]  , vertex[1], query_point);
        B[0] = CGAL::area(vertex[n-1], vertex[1], query_point);

        for(int i = 1; i < n-1; ++i) {
            s[i] = vertex[i] - query_point;
            r[i] = Scalar(CGAL::sqrt(CGAL::to_double(s[i].squared_length())));

            A[i] = CGAL::area(vertex[i]  , vertex[i+1], query_point);
            B[i] = CGAL::area(vertex[i-1], vertex[i+1], query_point);
        }

        s[n-1] = vertex[n-1] - query_point;
        r[n-1] = Scalar(CGAL::sqrt(CGAL::to_double(s[n-1].squared_length())));

        A[n-1] = CGAL::area(vertex[n-1], vertex[0], query_point);
        B[n-1] = CGAL::area(vertex[n-2], vertex[0], query_point);

        // Following section 4.2 from [2] we denote P_j = r_j*r_{j+1} + dot_product(d_j, d_{j+1}).
        // Vector s_i from [1] corresponds that one with the name d_i in [2].
        for(int j = 0; j < n-1; ++j) P[j] = r[j]*r[j+1] + s[j]*s[j+1];
        P[n-1] = r[n-1]*r[0] + s[n-1]*s[0];

        // Compute mean value weights using formula (16) from [2].
        // Since formula (16) always gives positive values, we have to add a proper sign to all the weight functions.
        weight[0] = r[n-1]*r[1] - s[n-1]*s[1];
        for(int j = 1; j < n-1; ++j) weight[0] *= P[j];
        weight[0] = sign_of_weight(A[n-1], A[0], B[0])*Scalar(CGAL::sqrt(CGAL::to_double(weight[0])));

        for(int i = 1; i < n-1; ++i) {
            weight[i] = r[i-1]*r[i+1] - s[i-1]*s[i+1];
            for(int j = 0; j < i-1; ++j) weight[i] *= P[j];
            for(int j = i+1; j < n; ++j) weight[i] *= P[j];
            weight[i] = sign_of_weight(A[i-1], A[i], B[i])*Scalar(CGAL::sqrt(CGAL::to_double(weight[i])));
        }

        weight[n-1] = r[n-2]*r[0] - s[n-2]*s[0];
        for(int j = 0; j < n-2; ++j) weight[n-1] *= P[j];
        weight[n-1] = sign_of_weight(A[n-2], A[n-1], B[n-1])*Scalar(CGAL::sqrt(CGAL::to_double(weight[n-1])));

        // Compute the sum of all weights - denominator of mean value coordinates.
        mv_denominator = weight[0];
        for(int i = 1; i < n; ++i) mv_denominator += weight[i];

        // Invert this denominator.
        CGAL_precondition( mv_denominator != Scalar(0) );
        inverted_mv_denominator = Scalar(1) / mv_denominator;

        // Normalize weights and save them as resulting mean value coordinates.
        for(int i = 0; i < n-1; ++i) {
            *output = weight[i] * inverted_mv_denominator;
            ++output;
        }
        *output = weight[n-1] * inverted_mv_denominator;

        // Return coordinates.
        return std::make_pair(output, true);
    }

    // Compute 2D mean value coordinates on the bounded side of the polygon with the fast O(n) but less precise algorithm.
    // Here, n - is the number of the polygon's vertices. Precision is lost near the boundary (~ 1.0e-10 and closer).
    template<class Iterator>
        std::pair<Iterator, bool> coordinates_on_bounded_side_fast_2(const Point_2 &query_point, Iterator &output)
    {
        // Get the number of vertices in the polygon.
        const int n = int(number_of_vertices);

        // Compute vector s following the pseudo-code in the Figure 10 from [1].
        for(int i = 0; i < n; ++i) s[i] = vertex[i] - query_point;

        // Compute length r, area A, and dot product D following the pseudo-code in the Figure 10 from [1].
        // Split the loop to make this computation faster.
        r[0] = Scalar(CGAL::sqrt(CGAL::to_double(s[0].squared_length())));
        A[0] = CGAL::area(vertex[0], vertex[1], query_point);
        D[0] = s[0]*s[1];

        for(int i = 1; i < n-1; ++i) {
            r[i] = Scalar(CGAL::sqrt(CGAL::to_double(s[i].squared_length())));
            A[i] = CGAL::area(vertex[i], vertex[i+1], query_point);
            D[i] = s[i]*s[i+1];
        }

        r[n-1] = Scalar(CGAL::sqrt(CGAL::to_double(s[n-1].squared_length())));
        A[n-1] = CGAL::area(vertex[n-1], vertex[0], query_point);
        D[n-1] = s[n-1]*s[0];

        // Compute intermediate values t using the formulas from slide 19 here
        // - http://www.inf.usi.ch/hormann/nsfworkshop/presentations/Hormann.pdf 
        for(int i = 0; i < n-1; ++i) {
            CGAL_precondition( (r[i]*r[i+1] + D[i]) != Scalar(0) );
            t[i] = A[i] / (r[i]*r[i+1] + D[i]);
        }

        CGAL_precondition( (r[n-1]*r[0] + D[n-1]) != Scalar(0) );
        t[n-1] = A[n-1] / (r[n-1]*r[0] + D[n-1]);

        // Compute mean value weights using the same pseudo-code as before.
        CGAL_precondition( r[0] != Scalar(0) );
        weight[0] = (t[n-1] + t[0]) / r[0];

        for(int i = 1; i < n-1; ++i) {
            CGAL_precondition( r[i] != Scalar(0) );
            weight[i] = (t[i-1] + t[i]) / r[i];
        }

        CGAL_precondition( r[n-1] != Scalar(0) );
        weight[n-1] = (t[n-2] + t[n-1]) / r[n-1];

        // Compute the sum of all weights - denominator of mean value coordinates.
        mv_denominator = weight[0];
        for(int i = 1; i < n; ++i) mv_denominator += weight[i];

        // Invert this denominator.
        CGAL_precondition( mv_denominator != Scalar(0) );
        inverted_mv_denominator = Scalar(1) / mv_denominator;

        // Normalize weights and save them as resulting mean value coordinates.
        for(int i = 0; i < n-1; ++i) {
            *output = weight[i] * inverted_mv_denominator;
            ++output;
        }
        *output = weight[n-1] * inverted_mv_denominator;

        // Return coordinates.
        return std::make_pair(output, true);
    }

    // COORDINATES ON UNBOUNDED SIDE.

    // Compute 2D mean value coordinates on the unbounded side of the polygon with the slow O(n^2) but precise algorithm.
    // Here, n - is the number of the polygon's vertices.
    template<class Iterator>
        std::pair<Iterator, bool> coordinates_on_unbounded_side_precise_2(const Point_2 &query_point, Iterator &output)
    {
        // Use the same formulas as for the bounded side since they are also valid on the unbounded side.
        return coordinates_on_bounded_side_precise(query_point, output);
    }

    // Compute 2D mean value coordinates on the unbounded side of the polygon with the fast O(n) but less precise algorithm.
    // Here, n - is the number of the polygon's vertices. Precision is lost near the boundary (~ 1.0e-10 and closer).
    template<class Iterator>
        std::pair<Iterator, bool> coordinates_on_unbounded_side_fast_2(const Point_2 &query_point, Iterator &output)
    {
        // Use the same formulas as for the bounded side since they are also valid on the unbounded side.
        return coordinates_on_bounded_side_fast(query_point, output);
    }

    // OTHER FUNCTIONS.

    // Return sign of current mean value weight function.
    // We can have 3 different values: 0 if the weight = 0, -1 if the weight is negative, and +1 if the weight is positive.
    inline Scalar sign_of_weight(const Scalar &A_prev, const Scalar &A, const Scalar &B) const
    {
        if(A_prev > Scalar(0) && A > Scalar(0) && B <= Scalar(0)) return  Scalar(1);
        if(A_prev < Scalar(0) && A < Scalar(0) && B >= Scalar(0)) return Scalar(-1);
        if(B > Scalar(0)) return  Scalar(1);
        if(B < Scalar(0)) return Scalar(-1);

        return Scalar(0);
    }

    // Print some information about 2D mean value coordinates.
    void print_coordinates_information_2(std::ostream &output_stream) const
    {
        output_stream << std::endl << "TYPE OF COORDINATES: " << std::endl << std::endl;
        output_stream << "The coordinate functions to compute are mean value coordinates." << std::endl;

        output_stream << std::endl << "INFORMATION ABOUT COORDINATES: " << std::endl << std::endl;
        output_stream << "Mean value coordinates are well-defined for an arbitrary simple polygon and can be computed only approximately due to the involved computation of the square root." << std::endl;

        output_stream << std::endl;
        output_stream << "They satisfy the following properties: " << std::endl;
        output_stream << "1. Partition of unity or constant precision;" << std::endl;
        output_stream << "2. Homogeneity or linear precision;" << std::endl;
        output_stream << "3. Lagrange property;" << std::endl;
        output_stream << "4. Linearity along edges;" << std::endl;
        output_stream << "5. Smoothness;" << std::endl;
        output_stream << "6. Similarity invariance;" << std::endl;
        output_stream << "7. Linear independence;" << std::endl;
        output_stream << "8. Refinability;" << std::endl;

        output_stream << std::endl;
        output_stream << "Mean value coordinates satisfy the non-negativity and boundedness between 0 and 1 properties inside the kernel of an arbitrary star-shaped polygon, too." << std::endl;

        output_stream << std::endl << "REFERENCES: " << std::endl << std::endl;
        output_stream << "K. Hormann and M. Floater. Mean value coordinates for arbitrary planar polygons. ACM Transactions on Graphics, 25(4):1424-1441, 2006." << std::endl;
        output_stream << "M. S. Floater, Wachspress and mean value coordinates, to appear in the Proceedings of the 14th International Conference on Approximation Theory, G. Fasshauer and L. L. Schumaker (eds.)." << std::endl;
    }
};

} // namespace Barycentric_coordinates

} // namespace CGAL

#endif // CGAL_MEAN_VALUE_COORDINATES_2_H