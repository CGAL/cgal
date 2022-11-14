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
// Author(s) : Dmitry Anisimov, David Bommes, Kai Hormann, and Pierre Alliez.

#ifndef CGAL_BARYCENTRIC_MEAN_VALUE_2_H
#define CGAL_BARYCENTRIC_MEAN_VALUE_2_H

#include <CGAL/license/Barycentric_coordinates_2.h>
#include <CGAL/disable_warnings.h>

// STL headers.
#include <vector>

// CGAL headers.
#include <CGAL/utils.h>
#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>

// Boost headers.
#include <boost/mpl/has_xxx.hpp>
#include <boost/optional/optional.hpp>

// Barycentric coordinates headers.
#include <CGAL/Barycentric_coordinates_2/barycentric_enum_2.h>

// CGAL namespace.
namespace CGAL {

// Barycentric coordinates namespace.
namespace Barycentric_coordinates {

#if !defined(CGAL_NO_DEPRECATED_CODE) || defined(DOXYGEN_RUNNING)

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Try to find a square root object in the provided `Traits` class. If not, then use the default square root from CGAL.

// Finds a square root of the provided value of the type `Kernel::FT` by first converting it to the double type and then taking the square root using the `CGAL::sqrt()` function.
template<class Traits>
    class Default_sqrt
{
    typedef typename Traits::FT FT;

public:
    FT operator()(const FT &value) const
    {
        return FT(CGAL::sqrt(CGAL::to_double(value)));
    }
};

BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_nested_type_Sqrt, Sqrt, false)

// Case: do_not_use_default = false.
template<class Traits, bool do_not_use_default = Has_nested_type_Sqrt<Traits>::value>
    class Get_sqrt
{
public:
    typedef Default_sqrt<Traits> Sqrt;

    static Sqrt sqrt_object(const Traits&)
    {
        return Sqrt();
    }
};

// Case: do_not_use_default = true.
template<class Traits>
    class Get_sqrt<Traits, true>
{
public:
    typedef typename Traits::Sqrt Sqrt;

    static Sqrt sqrt_object(const Traits &traits)
    {
        return traits.sqrt_object();
    }
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Examples: see the User Manual here - https://doc.cgal.org/latest/Manual/index.html.
// [1] Reference: "K. Hormann and M. Floater. Mean value coordinates for arbitrary planar polygons. ACM Transactions on Graphics, 25(4):1424-1441, 2006.".
// [2] Reference: "M. S. Floater, Wachspress and mean value coordinates, to appear in the Proceedings of the 14th International Conference on Approximation Theory, G. Fasshauer and L. L. Schumaker (eds.)."

/*!
 * \ingroup PkgBarycentricCoordinates2RefDeprecated
 * The class `Mean_value_2` implements 2D mean value coordinates ( \cite cgal:bc:hf-mvcapp-06, \cite cgal:bc:fhk-gcbcocp-06, \cite cgal:f-mvc-03 ).
 * This class is parameterized by a traits class `Traits`, and it is used as a coordinate class to complete the class `Generalized_barycentric_coordinates_2`.
 * For a polygon with three vertices (triangle) it is better to use the class `Triangle_coordinates_2`.
 * Mean value coordinates can be computed only approximately due to an inevitable square root operation, and they are necessarily positive only inside the kernel of a star-shaped polygon and inside any quadrilateral.

 * \deprecated This part of the package is deprecated since the version 5.4 of \cgal.

\tparam Traits must be a model of the concept `BarycentricTraits_2`.

\cgalModels `BarycentricCoordinates_2`

*/
template<class Traits>
class
#ifndef DOXYGEN_RUNNING
CGAL_DEPRECATED_MSG("This part of the package is deprecated since the version 5.4 of CGAL!")
#endif
Mean_value_2
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

    // Creates the class `Mean_value_2` that implements the behavior of mean value coordinates for any query point that does not belong to the polygon's boundary.
    // The polygon is given by a range of vertices of the type `Traits::Point_2` stored in a container of the type <a href="https://en.cppreference.com/w/cpp/container/vector">`std::vector`</a>.
    Mean_value_2(const std::vector<typename Traits::Point_2> &vertices, const Traits &b_traits) :
        vertex(vertices),
        barycentric_traits(b_traits),
        number_of_vertices(vertex.size()),
        area_2(barycentric_traits.compute_area_2_object()),
        squared_length_2(barycentric_traits.compute_squared_length_2_object()),
        scalar_product_2(barycentric_traits.compute_scalar_product_2_object()),
        sqrt(Get_sqrt<Traits>::sqrt_object(barycentric_traits))
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

    // Computation of Mean Value Weight Functions

    // This function computes mean value weights (unnormalized coordinates) for a chosen query point.
    template<class OutputIterator>
        inline boost::optional<OutputIterator> weights(const Point_2 &query_point, OutputIterator &output)
    {
        return weights_2(query_point, output);
    }

    // Computation of Mean Value Basis Functions

    // This function computes mean value barycentric coordinates for a chosen query point on the bounded side of a simple polygon.
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

    // This function computes mean value barycentric coordinates for a chosen query point on the unbounded side of a simple polygon.
    template<class OutputIterator>
        inline boost::optional<OutputIterator> coordinates_on_unbounded_side(const Point_2 &query_point, OutputIterator &output, const Type_of_algorithm type_of_algorithm)
    {
        switch(type_of_algorithm)
        {
            case PRECISE:
            return coordinates_on_unbounded_side_precise_2(query_point, output);
            break;

            case FAST:
            return coordinates_on_unbounded_side_fast_2(query_point, output);
            break;
        }

        // Pointer cannot be here. Something went wrong.
        const bool type_of_algorithm_failure = true;
        CGAL_postcondition( !type_of_algorithm_failure );
        if(!type_of_algorithm_failure) return boost::optional<OutputIterator>(output);
        else return boost::optional<OutputIterator>();
    }

    // Information Functions

    // This function prints some information about mean value coordinates.
    void print_coordinates_information(std::ostream &output_stream) const
    {
        return print_coordinates_information_2(output_stream);
    }

private:

    // Some convenient typedefs.
    typedef typename std::vector<FT>      FT_vector;
    typedef typename std::vector<Point_2> Point_vector;
    typedef typename Traits::Vector_2     Vector;
    typedef typename std::vector<Vector>  Vector_vector;

    // Internal global variables.
    const Point_vector &vertex;

    const Traits &barycentric_traits;

    const size_t number_of_vertices;

    Vector_vector s;

    FT_vector r, A, B, D, P, t, weight;

    FT mv_denominator, inverted_mv_denominator;

    typename Traits::Compute_area_2 area_2;
    typename Traits::Compute_squared_length_2 squared_length_2;
    typename Traits::Compute_scalar_product_2 scalar_product_2;
    typename Get_sqrt<Traits>::Sqrt sqrt;

    // WEIGHTS.

    // Compute mean value weights without normalization.
    template<class OutputIterator>
        boost::optional<OutputIterator> weights_2(const Point_2 &query_point, OutputIterator &output)
    {
        // Get the number of vertices in the polygon.
        const int n = int(number_of_vertices);

        // Compute vectors s following the pseudo-code in the Figure 10 from [1].
        for(int i = 0; i < n; ++i) s[i] = vertex[i] - query_point;

        // Compute lengths r, areas A, and dot products D following the pseudo-code in the Figure 10 from [1].
        // Split the loop to make this computation faster.
        r[0] = sqrt(squared_length_2(s[0]));
        A[0] = area_2(vertex[0], vertex[1], query_point);
        D[0] = scalar_product_2(s[0], s[1]);

        for(int i = 1; i < n-1; ++i) {
            r[i] = sqrt(squared_length_2(s[i]));
            A[i] = area_2(vertex[i], vertex[i+1], query_point);
            D[i] = scalar_product_2(s[i], s[i+1]);
        }

        r[n-1] = sqrt(squared_length_2(s[n-1]));
        A[n-1] = area_2(vertex[n-1], vertex[0], query_point);
        D[n-1] = scalar_product_2(s[n-1], s[0]);

        // Compute intermediate values t using the formulas from slide 19 here
        // - https://www.inf.usi.ch/hormann/nsfworkshop/presentations/Hormann.pdf
        for(int i = 0; i < n-1; ++i) {
            CGAL_precondition( (r[i]*r[i+1] + D[i]) != FT(0) );
            t[i] = A[i] / (r[i]*r[i+1] + D[i]);
        }

        CGAL_precondition( (r[n-1]*r[0] + D[n-1]) != FT(0) );
        t[n-1] = A[n-1] / (r[n-1]*r[0] + D[n-1]);

        // Compute mean value weights using the same pseudo-code as before.
        CGAL_precondition( r[0] != FT(0) );
        *output = (t[n-1] + t[0]) / r[0];
        ++output;

        for(int i = 1; i < n-1; ++i) {
            CGAL_precondition( r[i] != FT(0) );
            *output = (t[i-1] + t[i]) / r[i];
            ++output;
        }

        CGAL_precondition( r[n-1] != FT(0) );
        *output = (t[n-2] + t[n-1]) / r[n-1];
        ++output;

        // Return weights.
        return boost::optional<OutputIterator>(output);
    }

    // COORDINATES ON BOUNDED SIDE.

    // Compute mean value coordinates on the bounded side of the polygon with the slow O(n^2) but precise algorithm.
    // Here, n - is the number of the polygon's vertices.
    template<class OutputIterator>
        boost::optional<OutputIterator> coordinates_on_bounded_side_precise_2(const Point_2 &query_point, OutputIterator &output)
    {
        // Get the number of vertices in the polygon.
        const int n = int(number_of_vertices);

        // Compute vectors s and its lengths r following the pseudo-code in the Figure 10 from [1].
        s[0] = vertex[0] - query_point;
        r[0] = sqrt(squared_length_2(s[0]));

        // Compute areas A and B following the notation from [1] (see Figure 2). Split the loop to make this computation faster.
        A[0] = area_2(vertex[0]  , vertex[1], query_point);
        B[0] = area_2(vertex[n-1], vertex[1], query_point);

        for(int i = 1; i < n-1; ++i) {
            s[i] = vertex[i] - query_point;
            r[i] = sqrt(squared_length_2(s[i]));

            A[i] = area_2(vertex[i]  , vertex[i+1], query_point);
            B[i] = area_2(vertex[i-1], vertex[i+1], query_point);
        }

        s[n-1] = vertex[n-1] - query_point;
        r[n-1] = sqrt(squared_length_2(s[n-1]));

        A[n-1] = area_2(vertex[n-1], vertex[0], query_point);
        B[n-1] = area_2(vertex[n-2], vertex[0], query_point);

        // Following section 4.2 from [2] we denote P_j = r_j*r_{j+1} + dot_product(d_j, d_{j+1}).
        // Vector s_i from [1] corresponds to that one with the name d_i in [2].
        for(int j = 0; j < n-1; ++j)
            P[j] = CGAL::max<FT>(r[j]*r[j+1] + scalar_product_2(s[j], s[j+1]), 0);
        P[n-1] = CGAL::max<FT>(r[n-1]*r[0] + scalar_product_2(s[n-1], s[0]), 0);

        // Compute mean value weights using the formula (16) from [2].
        // Since the formula (16) always gives positive values, we have to add a proper sign to all the weight functions.
        weight[0] = r[n-1]*r[1] - scalar_product_2(s[n-1], s[1]);
        for(int j = 1; j < n-1; ++j) weight[0] *= P[j];
        weight[0] = sign_of_weight(A[n-1], A[0], B[0])*sqrt(weight[0]);

        for(int i = 1; i < n-1; ++i) {
            weight[i] = r[i-1]*r[i+1] - scalar_product_2(s[i-1], s[i+1]);
            for(int j = 0; j < i-1; ++j) weight[i] *= P[j];
            for(int j = i+1; j < n; ++j) weight[i] *= P[j];
            weight[i] = sign_of_weight(A[i-1], A[i], B[i])*sqrt(weight[i]);
        }

        weight[n-1] = r[n-2]*r[0] - scalar_product_2(s[n-2], s[0]);
        for(int j = 0; j < n-2; ++j) weight[n-1] *= P[j];
        weight[n-1] = sign_of_weight(A[n-2], A[n-1], B[n-1])*sqrt(weight[n-1]);

        // Compute the sum of all weights - denominator of mean value coordinates.
        mv_denominator = weight[0];
        for(int i = 1; i < n; ++i) mv_denominator += weight[i];

        // Invert this denominator.
        CGAL_precondition( mv_denominator != FT(0) );
        inverted_mv_denominator = FT(1) / mv_denominator;

        // Normalize weights and save them as resulting mean value coordinates.
        for(int i = 0; i < n-1; ++i) {
            *output = weight[i] * inverted_mv_denominator;
            ++output;
        }
        *output = weight[n-1] * inverted_mv_denominator;
        ++output;

        // Return coordinates.
        return boost::optional<OutputIterator>(output);
    }

    // Compute mean value coordinates on the bounded side of the polygon with the fast O(n) but less precise algorithm.
    // Here, n - is the number of the polygon's vertices. Precision is lost near the boundary (~ 1.0e-10 and closer).
    template<class OutputIterator>
        boost::optional<OutputIterator> coordinates_on_bounded_side_fast_2(const Point_2 &query_point, OutputIterator &output)
    {
        // Get the number of vertices in the polygon.
        const int n = int(number_of_vertices);

        // Compute vectors s following the pseudo-code in the Figure 10 from [1].
        for(int i = 0; i < n; ++i) s[i] = vertex[i] - query_point;

        // Compute lengths r, areas A, and dot products D following the pseudo-code in the Figure 10 from [1].
        // Split the loop to make this computation faster.
        r[0] = sqrt(squared_length_2(s[0]));
        A[0] = area_2(vertex[0], vertex[1], query_point);
        D[0] = scalar_product_2(s[0], s[1]);

        for(int i = 1; i < n-1; ++i) {
            r[i] = sqrt(squared_length_2(s[i]));
            A[i] = area_2(vertex[i], vertex[i+1], query_point);
            D[i] = scalar_product_2(s[i], s[i+1]);
        }

        r[n-1] = sqrt(squared_length_2(s[n-1]));
        A[n-1] = area_2(vertex[n-1], vertex[0], query_point);
        D[n-1] = scalar_product_2(s[n-1], s[0]);

        // Compute intermediate values t using the formulas from slide 19 here
        // - https://www.inf.usi.ch/hormann/nsfworkshop/presentations/Hormann.pdf
        for(int i = 0; i < n-1; ++i) {
            CGAL_precondition( (r[i]*r[i+1] + D[i]) != FT(0) );
            t[i] = A[i] / (r[i]*r[i+1] + D[i]);
        }

        CGAL_precondition( (r[n-1]*r[0] + D[n-1]) != FT(0) );
        t[n-1] = A[n-1] / (r[n-1]*r[0] + D[n-1]);

        // Compute mean value weights using the same pseudo-code as before.
        CGAL_precondition( r[0] != FT(0) );
        weight[0] = (t[n-1] + t[0]) / r[0];

        for(int i = 1; i < n-1; ++i) {
            CGAL_precondition( r[i] != FT(0) );
            weight[i] = (t[i-1] + t[i]) / r[i];
        }

        CGAL_precondition( r[n-1] != FT(0) );
        weight[n-1] = (t[n-2] + t[n-1]) / r[n-1];

        // Compute the sum of all weights - denominator of mean value coordinates.
        mv_denominator = weight[0];
        for(int i = 1; i < n; ++i) mv_denominator += weight[i];

        // Invert this denominator.
        CGAL_precondition( mv_denominator != FT(0) );
        inverted_mv_denominator = FT(1) / mv_denominator;

        // Normalize weights and save them as resulting mean value coordinates.
        for(int i = 0; i < n-1; ++i) {
            *output = weight[i] * inverted_mv_denominator;
            ++output;
        }
        *output = weight[n-1] * inverted_mv_denominator;
        ++output;

        // Return coordinates.
        return boost::optional<OutputIterator>(output);
    }

    // COORDINATES ON UNBOUNDED SIDE.

    // Compute mean value coordinates on the unbounded side of the polygon with the slow O(n^2) but precise algorithm.
    // Here, n - is the number of the polygon's vertices.
    template<class OutputIterator>
        boost::optional<OutputIterator> coordinates_on_unbounded_side_precise_2(const Point_2 &query_point, OutputIterator &output)
    {
        // Use the same formulas as for the bounded side since they are also valid on the unbounded side.
        return coordinates_on_bounded_side_precise_2(query_point, output);
    }

    // Compute mean value coordinates on the unbounded side of the polygon with the fast O(n) but less precise algorithm.
    // Here, n - is the number of the polygon's vertices. Precision is lost near the boundary (~ 1.0e-10 and closer).
    template<class OutputIterator>
        boost::optional<OutputIterator> coordinates_on_unbounded_side_fast_2(const Point_2 &query_point, OutputIterator &output)
    {
        // Use the same formulas as for the bounded side since they are also valid on the unbounded side.
        return coordinates_on_bounded_side_fast_2(query_point, output);
    }

    // OTHER FUNCTIONS.

    // Return the sign of a mean value weight function.
    // We can have 3 different values: 0 if the weight = 0, -1 if the weight is negative, and +1 if the weight is positive.
    inline FT sign_of_weight(const FT &A_prev, const FT &A, const FT &B) const
    {
        if(A_prev > FT(0) && A > FT(0) && B <= FT(0)) return FT(1);
        if(A_prev < FT(0) && A < FT(0) && B >= FT(0)) return FT(-1);
        if(B > FT(0)) return FT(1);
        if(B < FT(0)) return FT(-1);

        return FT(0);
    }

    // Print some information about mean value coordinates.
    void print_coordinates_information_2(std::ostream &output_stream) const
    {
        output_stream << std::endl << "TYPE OF COORDINATES: " << std::endl << std::endl;
        output_stream << "The coordinate functions to be computed are mean value coordinates." << std::endl;

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
        output_stream << "8. Refinability." << std::endl;

        output_stream << std::endl;
        output_stream << "Mean value coordinates satisfy the non-negativity and boundedness between 0 and 1 properties inside the kernel of an arbitrary star-shaped polygon, too." << std::endl;

        output_stream << std::endl << "REFERENCES: " << std::endl << std::endl;
        output_stream << "K. Hormann and M. Floater. Mean value coordinates for arbitrary planar polygons. ACM Transactions on Graphics, 25(4):1424-1441, 2006." << std::endl;
        output_stream << "M. S. Floater, Wachspress and mean value coordinates, to appear in the Proceedings of the 14th International Conference on Approximation Theory, G. Fasshauer and L. L. Schumaker (eds.)." << std::endl;
    }
};

#endif // CGAL_NO_DEPRECATED_CODE

} // namespace Barycentric_coordinates

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_BARYCENTRIC_MEAN_VALUE_2_H
