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

// Barycentric coordinates headers.
#include <CGAL/Barycentric_coordinates_base_2.h>

// CGAL namespace.
namespace CGAL {

// Barycentric coordinates namespace.
namespace Barycentric_coordinates {

// Examples: See User Manual here - http://doc.cgal.org/latest/Manual/index.html.
// [1] Reference: "K. Hormann and M. Floater. Mean value coordinates for arbitrary planar polygons. ACM Transactions on Graphics, 25(4):1424-1441, 2006.".
// [2] Reference: "M. S. Floater, Wachspress and mean value coordinates, to appear in the Proceedings of the 14th International Conference on Approximation Theory, G. Fasshauer and L. L. Schumaker (eds.)."

/*!
 * \ingroup PkgBarycentric_coordinates_2
 * The class Mean_value_coordinates_2 implements 2D Mean Value coordinates ( \cite cgal:bc:hf-mvcapp-06, \cite cgal:bc:fhk-gcbcocp-06, \cite cgal:f-mvc-03 )
 * with respect to an arbitrary simple polygon. This class is parameterized by the `CGAL::Polygon_2` and the `Iterator` class.
 * The latter can be any class that fulfills the requirements for an STL iterator. This class is derived
 * from the class `CGAL::Barycentric_coordinates::Barycentric_coordinates_base_2`.
 * For a polygon with three vertices (triangle) it is better to use the class `CGAL::Barycentric_coordinates::Triangle_coordinates_2`.
 * Mean Value coordinates can be computed only approximately due to an inevitable square root operation,
 * and they are necesserily positive only inside the kernel of a star-shaped polygon and inside any quadrilateral.
 *
 * \sa `Iterator`
 *
 */

// This class does not allow the user to use different iterators to output coordinates after
// it has been created. In order to do so, we need to move template parameter Iterator in declaration of the class's functions.
// Can we also use reference in front of iterator like OutputIterator &output when passing it to the function?  
template<typename Polygon_2, typename Iterator> 
    class Mean_value_coordinates_2 : public Barycentric_coordinates_base_2<Polygon_2, Iterator>
{

public:

    // Creation.

    /// \name Types
    /// @{

    /// Type of 2D polygon.
    typedef Polygon_2                   Polygon;

    /// Type of the used kernel.
    typedef typename Polygon::Traits    Kernel;

    /// Type of 2D segment.
    typedef typename Polygon::Segment_2 Segment;

    /// Type of 2D point.
    typedef typename Polygon::Point_2   Point;

    /// Number type.
    typedef typename Polygon::FT        Scalar;

    /// Type of the used output iterator.
    typedef Iterator                    OutputIterator;

    /// @}

    /// \name Creation
    /// @{

    /// Creates an instance of the class Mean_value_coordinates_2 for a provided polygon passed as a reference.
    /// For preconditions and functions to compute weights or coordinates
    /// see the class `CGAL::Barycentric_coordinates::Barycentric_coordinates_base_2`.
    Mean_value_coordinates_2(const Polygon &_polygon) :
        Barycentric_coordinates_base_2<Polygon, OutputIterator>(_polygon)
    {
        const int number_of_polygon_vertices = _polygon.size();

        // Resize all internal containers.
        s.resize(number_of_polygon_vertices);

        r.resize(number_of_polygon_vertices);
        A.resize(number_of_polygon_vertices);
        B.resize(number_of_polygon_vertices);
        D.resize(number_of_polygon_vertices);
        P.resize(number_of_polygon_vertices);
        t.resize(number_of_polygon_vertices);

        weight.resize(number_of_polygon_vertices);
    }

    /// @}

private:

    // Some convenient typedefs.
    typedef Barycentric_coordinates_base_2<Polygon, OutputIterator> Base;
    typedef typename std::vector<Scalar>                   Scalar_vector;
    typedef Vector_2<Kernel>                                      Vector;
    typedef typename std::vector<Vector>                   Vector_vector;

    // Internal global variables.
    Vector_vector s;

    Scalar_vector r, A, B, D, P, t, weight;

    Scalar mv_denominator, inverted_mv_denominator;

    // COMPUTATION.

    // Can we somehow use a referenced output: &output?

    // WEIGHTS.

    // Compute Mean Value weights without normalization. 
    std::pair<OutputIterator, bool> weights(const Point &query_point, OutputIterator output)
    {
        // Get number of vertices in the polygon.
        const int n = Base::number_of_polygon_vertices;

        // Compute vectors s following the pseudo-code in the Figure 10 from [1].
        for(int i = 0; i < n; ++i) s[i] = Base::polygon.vertex(i) - query_point;

        // Compute length r, area A, and dot product D following the pseudo-code in the Figure 10 from [1].
        // Split the loop to make computation faster.
        r[0] = Scalar(CGAL::sqrt(CGAL::to_double(s[0].squared_length())));
        A[0] = CGAL::area(Base::polygon.vertex(0), Base::polygon.vertex(1), query_point);
        D[0] = s[0]*s[1];

        for(int i = 1; i < n-1; ++i) {
            r[i] = Scalar(CGAL::sqrt(CGAL::to_double(s[i].squared_length())));
            A[i] = CGAL::area(Base::polygon.vertex(i), Base::polygon.vertex(i+1), query_point);
            D[i] = s[i]*s[i+1];
        }

        r[n-1] = Scalar(CGAL::sqrt(CGAL::to_double(s[n-1].squared_length())));
        A[n-1] = CGAL::area(Base::polygon.vertex(n-1), Base::polygon.vertex(0), query_point);
        D[n-1] = s[n-1]*s[0];

        // Compute intermediate values t using formulas from slide 19 here
        // - http://www.inf.usi.ch/hormann/nsfworkshop/presentations/Hormann.pdf 
        for(int i = 0; i < n-1; ++i) {
            CGAL_precondition( (r[i]*r[i+1] + D[i]) != Scalar(0) );
            t[i] = A[i] / (r[i]*r[i+1] + D[i]);
        }

        CGAL_precondition( (r[n-1]*r[0] + D[n-1]) != Scalar(0) );
        t[n-1] = A[n-1] / (r[n-1]*r[0] + D[n-1]);

        // Compute Mean Value weights using the same pseudo-code as before.
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

    // Compute Mean Value coordinates on the bounded side of the polygon with slow O(n^2) but precise algorithm.
    // Here, n - is a number of polygon's vertices.
    std::pair<OutputIterator, bool> coordinates_on_bounded_side_precise(const Point &query_point, OutputIterator output)
    {
        // Get number of vertices in the polygon.
        const int n = Base::number_of_polygon_vertices;

        // Compute vector s and its length r following the pseudo-code in the Figure 10 from [1].
        s[0] = Base::polygon.vertex(0) - query_point;
        r[0] = Scalar(CGAL::sqrt(CGAL::to_double(s[0].squared_length())));

        // Compute areas A and B following notations from [1] (see Figure 2). Split the loop to make computation faster.
        A[0] = CGAL::area(Base::polygon.vertex(0)  , Base::polygon.vertex(1), query_point);
        B[0] = CGAL::area(Base::polygon.vertex(n-1), Base::polygon.vertex(1), query_point);

        for(int i = 1; i < n-1; ++i) {
            s[i] = Base::polygon.vertex(i) - query_point;
            r[i] = Scalar(CGAL::sqrt(CGAL::to_double(s[i].squared_length())));

            A[i] = CGAL::area(Base::polygon.vertex(i)  , Base::polygon.vertex(i+1), query_point);
            B[i] = CGAL::area(Base::polygon.vertex(i-1), Base::polygon.vertex(i+1), query_point);
        }

        s[n-1] = Base::polygon.vertex(n-1) - query_point;
        r[n-1] = Scalar(CGAL::sqrt(CGAL::to_double(s[n-1].squared_length())));

        A[n-1] = CGAL::area(Base::polygon.vertex(n-1), Base::polygon.vertex(0), query_point);
        B[n-1] = CGAL::area(Base::polygon.vertex(n-2), Base::polygon.vertex(0), query_point);

        // Following section 4.2 from [2] we denote P_j = r_j*r_{j+1} + dot_product(d_j, d_{j+1}).
        // Vector s_i from [1] corresponds that one with name d_i in [2].
        for(int j = 0; j < n-1; ++j) P[j] = r[j]*r[j+1] + s[j]*s[j+1];
        P[n-1] = r[n-1]*r[0] + s[n-1]*s[0];

        // Compute Mean Value weights using formula (16) from [2].
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

        // Compute the sum of all weights - denominator of Mean Value coordinates.
        mv_denominator = weight[0];
        for(int i = 1; i < n; ++i) mv_denominator += weight[i];

        // Invert this denominator.
        CGAL_precondition( mv_denominator != Scalar(0) );
        inverted_mv_denominator = Scalar(1) / mv_denominator;

        // Normalize weights and save them as resulting Mean Value coordinates.
        for(int i = 0; i < n-1; ++i) {
            *output = weight[i] * inverted_mv_denominator;
            ++output;
        }
        *output = weight[n-1] * inverted_mv_denominator;

        // Return coordinates.
        return std::make_pair(output, true);
    }

    // Compute Mean Value coordinates on the bounded side of the polygon with fast O(n) but less precise algorithm.
    // Here, n - is number of polygon's vertices. Precision is lost next to the boundary.
    // Can we somehow reuse here the computations from the function weights()?
    std::pair<OutputIterator, bool> coordinates_on_bounded_side_fast(const Point &query_point, OutputIterator output)
    {
        // Get number of vertices in the polygon.
        const int n = Base::number_of_polygon_vertices;

        // Compute vectors s following the pseudo-code in the Figure 10 from [1].
        for(int i = 0; i < n; ++i) s[i] = Base::polygon.vertex(i) - query_point;

        // Compute length r, area A, and dot product D following the pseudo-code in the Figure 10 from [1].
        // Split the loop to make computation faster.
        r[0] = Scalar(CGAL::sqrt(CGAL::to_double(s[0].squared_length())));
        A[0] = CGAL::area(Base::polygon.vertex(0), Base::polygon.vertex(1), query_point);
        D[0] = s[0]*s[1];

        for(int i = 1; i < n-1; ++i) {
            r[i] = Scalar(CGAL::sqrt(CGAL::to_double(s[i].squared_length())));
            A[i] = CGAL::area(Base::polygon.vertex(i), Base::polygon.vertex(i+1), query_point);
            D[i] = s[i]*s[i+1];
        }

        r[n-1] = Scalar(CGAL::sqrt(CGAL::to_double(s[n-1].squared_length())));
        A[n-1] = CGAL::area(Base::polygon.vertex(n-1), Base::polygon.vertex(0), query_point);
        D[n-1] = s[n-1]*s[0];

        // Compute intermediate values t using formulas from slide 19 here
        // - http://www.inf.usi.ch/hormann/nsfworkshop/presentations/Hormann.pdf 
        for(int i = 0; i < n-1; ++i) {
            CGAL_precondition( (r[i]*r[i+1] + D[i]) != Scalar(0) );
            t[i] = A[i] / (r[i]*r[i+1] + D[i]);
        }

        CGAL_precondition( (r[n-1]*r[0] + D[n-1]) != Scalar(0) );
        t[n-1] = A[n-1] / (r[n-1]*r[0] + D[n-1]);

        // Compute Mean Value weights using the same pseudo-code as before.
        CGAL_precondition( r[0] != Scalar(0) );
        weight[0] = (t[n-1] + t[0]) / r[0];

        for(int i = 1; i < n-1; ++i) {
            CGAL_precondition( r[i] != Scalar(0) );
            weight[i] = (t[i-1] + t[i]) / r[i];
        }

        CGAL_precondition( r[n-1] != Scalar(0) );
        weight[n-1] = (t[n-2] + t[n-1]) / r[n-1];

        // Compute the sum of all weights - denominator of Mean Value coordinates.
        mv_denominator = weight[0];
        for(int i = 1; i < n; ++i) mv_denominator += weight[i];

        // Invert this denominator.
        CGAL_precondition( mv_denominator != Scalar(0) );
        inverted_mv_denominator = Scalar(1) / mv_denominator;

        // Normalize weights and save them as resulting Mean Value coordinates.
        for(int i = 0; i < n-1; ++i) {
            *output = weight[i] * inverted_mv_denominator;
            ++output;
        }
        *output = weight[n-1] * inverted_mv_denominator;

        // Return coordinates.
        return std::make_pair(output, true);
    }

    // COORDINATES ON UNBOUNDED SIDE.

    // Compute Mean Value coordinates on the unbounded side of the polygon with slow O(n^2) but precise algorithm.
    // Here, n - is a number of polygon's vertices.
    std::pair<OutputIterator, bool> coordinates_on_unbounded_side_precise(const Point &query_point, OutputIterator output)
    {
        // Use the same formulas as for bounded side since they are also valid on unbounded side.
        return coordinates_on_bounded_side_precise(query_point, output);
    }

    // Compute Mean Value coordinates on the unbounded side of the polygon with fast O(n) but less precise algorithm.
    // Here, n - is number of polygon's vertices. Precision is lost next to the boundary.
    std::pair<OutputIterator, bool> coordinates_on_unbounded_side_fast(const Point &query_point, OutputIterator output)
    {
        // Use the same formulas as for bounded side since they are also valid on unbounded side.
        return coordinates_on_bounded_side_fast(query_point, output);
    }

    // OTHER FUNCTIONS.

    // Return sign of current Mean Value weight function.
    // We can have 3 different values: 0 if weight = 0, -1 if weight is negative, and +1 if weight is positive.
    inline Scalar sign_of_weight(const Scalar &A_prev, const Scalar &A, const Scalar &B) const
    {
        if(A_prev > Scalar(0) && A > Scalar(0) && B <= Scalar(0)) return  Scalar(1);
        if(A_prev < Scalar(0) && A < Scalar(0) && B >= Scalar(0)) return Scalar(-1);
        if(B > Scalar(0)) return  Scalar(1);
        if(B < Scalar(0)) return Scalar(-1);

        return Scalar(0);
    }

    // Print some info about Mean Value coordinates.
    void print_coordinates_info() const
    {
        std::cout << std::endl << "TYPE OF COORDINATES: " << std::endl << std::endl;
        std::cout << "Currently computed coordinate functions are Mean Value coordinates." << std::endl;

        std::cout << std::endl << "INFORMATION ABOUT COORDINATES: " << std::endl << std::endl;
        std::cout << "Mean Value coordinates are well-defined for an arbitrary simple polygon, and they can be computed only approximately due to the involved computation of the square root." << std::endl;

        std::cout << std::endl;
        std::cout << "For these polygons, they satisfy the following properties: " << std::endl;
        std::cout << "1. Partition of unity or constant precision;" << std::endl;
        std::cout << "2. Homogeneity or linear precision;" << std::endl;
        std::cout << "3. Lagrange property;" << std::endl;
        std::cout << "4. Linearity along edges;" << std::endl;
        std::cout << "5. Smoothness;" << std::endl;
        std::cout << "6. Similarity invariance;" << std::endl;
        std::cout << "7. Linear independence;" << std::endl;
        std::cout << "8. Refinability;" << std::endl;

        std::cout << std::endl;
        std::cout << "Mean Value coordinates satisfy the Non-negativity and Boundedness between 0 and 1 properties inside the kernel of an arbitrary star-shaped polygon, too." << std::endl;

        std::cout << std::endl << "REFERENCES: " << std::endl << std::endl;
        std::cout << "K. Hormann and M. Floater. Mean value coordinates for arbitrary planar polygons. ACM Transactions on Graphics, 25(4):1424-1441, 2006." << std::endl;
        std::cout << "M. S. Floater, Wachspress and mean value coordinates, to appear in the Proceedings of the 14th International Conference on Approximation Theory, G. Fasshauer and L. L. Schumaker (eds.)." << std::endl;
    }
};

// Class Mean_value_coordinates_2 with particular std::back_insert_iterator instead of a general one.

/*!
 * \ingroup PkgBarycentric_coordinates_2
 * The class MV_coordinates_2 implements 2D Mean Value coordinates ( \cite cgal:bc:hf-mvcapp-06, \cite cgal:bc:fhk-gcbcocp-06, \cite cgal:f-mvc-03 )
 * with respect to an arbitrary simple polygon. This class is parameterized by the `CGAL::Polygon_2` and a Container class.
 * The latter can be any class that fulfills the requirements for the <a href="http://en.cppreference.com/w/cpp/iterator/back_insert_iterator">`std::back_insert_iterator`</a>. 
 * It defaults to <a href="http://en.cppreference.com/w/cpp/container/vector">`std::vector`</a> container. 
 * This class is derived from the class `CGAL::Barycentric_coordinates::Barycentric_coordinates_base_2`.
 * For a polygon with three vertices (triangle) it is better to use the class `CGAL::Barycentric_coordinates::Tri_coordinates_2`.
 * Mean Value coordinates can be computed only approximately due to an inevitable square root operation,
 * and they are necesserily positive only inside the kernel of a star-shaped polygon and inside any quadrilateral.
 *
 * \sa <a href="http://en.cppreference.com/w/cpp/iterator/back_insert_iterator">`std::back_insert_iterator`</a>
 *
 */

template<typename Polygon_2, typename InputContainer = std::vector<typename Polygon_2::FT> > 
    class MV_coordinates_2 : public Mean_value_coordinates_2<Polygon_2, std::back_insert_iterator<InputContainer> >
{

public:

    // Creation.

    /// \name Types
    /// @{

    /// Type of the used container.
    typedef InputContainer Container;

    /// Type of the base class.
    typedef Mean_value_coordinates_2<Polygon_2, std::back_insert_iterator<Container> > Base;

    /// @}

    /// \name Creation
    /// @{

    /// Creates an instance of the class Mean_value_coordinates_2 for a provided polygon passed as a reference.
    /// The used iterator is <a href="http://en.cppreference.com/w/cpp/iterator/back_insert_iterator">`std::back_insert_iterator`</a>.
    MV_coordinates_2(const Polygon_2 &_polygon) : Base(_polygon) { }

    /// @}

    // Weights.

    /// \name Computation of weight functions
    /// @{

    /// Computes Mean Value weights for any strictly interior query point with respect to all the vertices of the polygon.
    /// \pre `_polygon.bounded_side(query_point) == CGAL::ON_BOUNDED_SIDE`
    inline std::pair<std::back_insert_iterator<Container>, bool> compute_weights(const typename Polygon_2::Point_2 &query_point, Container &container)
    {
        return Base::compute_weights(query_point, std::back_inserter(container));
    }

    /// @}

    // Coordinates.

    /// \name Computation of basis functions at the vertices (with index)
    /// @{

    /// Computes Mean Value coordinates for a query point, which coincides with one of the polygon's vertices, with known index.
    /// \pre `(0 <= index) && (index < number_of_polygon_vertices)`
    inline std::pair<std::back_insert_iterator<Container>, bool> compute_at_vertex(const int index, Container &container) const
    {
        return Base::compute_at_vertex(index, std::back_inserter(container));
    }

    /// @}

    /// \name Computation of basis functions along edges (with index)
    /// @{

    /// Computes Mean Value coordinates for a query point on the polygon's boundary with known index of the edge to which this point belongs.
    /// \pre `_polygon.bounded_side(query_point) == CGAL::ON_BOUNDARY`
    /// \pre `(0 <= index) && (index < number_of_polygon_vertices)`
    inline std::pair<std::back_insert_iterator<Container>, bool> compute_on_edge(const typename Polygon_2::Point_2 &query_point, const int index, Container &container) const
    {    
        return Base::compute_on_edge(query_point, index, std::back_inserter(container));
    }

    /// @}

    /// \name Computation of basis functions at an arbitrary point
    /// @{

    /// Computes Mean Value coordinates for any query point in the plane with respect to all the vertices of the polygon.
    inline std::pair<std::back_insert_iterator<Container>, bool> compute(const typename Polygon_2::Point_2 &query_point, Container &container, Query_point_location query_point_location = UNSPECIFIED_LOCATION, Type_of_algorithm type_of_algorithm = PRECISE)
    {   
        return Base::compute(query_point, std::back_inserter(container), query_point_location, type_of_algorithm);
    }

    /// @}
};

} // namespace Barycentric_coordinates

} // namespace CGAL

#endif // CGAL_MEAN_VALUE_COORDINATES_2_H