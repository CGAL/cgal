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
  \file Segment_coordinates_2.h
*/

#ifndef CGAL_SEGMENT_COORDINATES_2_H
#define CGAL_SEGMENT_COORDINATES_2_H

// STL headers.  
#include <vector>

// CGAL headers.
#include <CGAL/array.h>  
#include <CGAL/assertions.h>

// CGAL namespace.
namespace CGAL {

// Barycentric coordinates namespace.
namespace Barycentric_coordinates {

// Examples: see the User Manual - http://doc.cgal.org/latest/Manual/index.html.

/*!
 * \ingroup PkgBarycentric_coordinates_2
 * The class `Segment_coordinates_2` implements barycentric coordinates with respect to an arbitrary non-degenerate segment along an arbitrary line in the plane.
 * This class is parameterized by a traits class `Traits`.

\cgalHeading{Template parameters}

\tparam Traits must be a model of the concept `BarycentricTraits_2`. In particular, it must provide the functions `Kernel::Compute_scalar_product_2`, `Kernel::Compute_squared_distance_2`, and `Kernel::Equal_2`. 

*/

template<class Traits> 
    class Segment_coordinates_2
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

    /// Creates the class `Segment_coordinates_2` that implements the behaviour of segment coordinates with respect to an arbitrary non-degenerate segment along an arbitrary line in the plane.
    /// The segment is given by its two vertices.
    /// \pre Segment is not degenerate.
    Segment_coordinates_2(const Point &first_vertex, const Point &second_vertex, const Traits &b_traits = Traits()) :
        barycentric_traits(b_traits),
        vertex_0(first_vertex),
        vertex_1(second_vertex),
        scalar_product_2(barycentric_traits.compute_scalar_product_2_object()),
        squared_distance_2(barycentric_traits.compute_squared_distance_2_object()),
        equal_2(barycentric_traits.equal_2_object()) 
    {
        CGAL_precondition( !equal_2(vertex_0, vertex_1) );
    }

    /// @}

    /// \name Computation of Basis Functions
    /// @{

    /// Computes segment barycentric coordinates for a chosen query point with respect to both vertices of the segment.
    /// This function accepts any STL like iterator, which complies with the `Iterator` concept.
    template<class Iterator> 
        inline std::pair<Iterator, bool> compute(const Point &query_point, Iterator output)
    {
        return segment_coordinates_2(query_point, output);
    }

    /// Computes segment barycentric coordinates for a chosen query point with respect to both vertices of the segment.
    /// This function accepts a container of the type <a href="http://en.cppreference.com/w/cpp/container/vector">`std::vector`</a> 
    /// and returns an iterator of the type <a href="http://en.cppreference.com/w/cpp/iterator/back_insert_iterator">`std::back_insert_iterator`</a>
    /// that is placed past-the-end of the resulting sequence of coordinate values.
    inline std::pair<std::back_insert_iterator<std::vector<Scalar> >, bool> compute(const Point &query_point, std::vector<Scalar> &output_vector)
    {
        output_vector.reserve(output_vector.size() + 2); 
        typedef typename std::back_insert_iterator<std::vector<Scalar> > Iterator;
        Iterator output = std::back_inserter(output_vector);
        return segment_coordinates_2(query_point, output);
    }

    /// This is a static function that takes both vertices of a segment and computes segment coordinates at a given query point with respect to these vertices.
    /// The coordinate values are returned as `CGAL::cpp11::array<FT,2>`.
    /// The function also requires a traits class of the concept `BarycentricTraits_2`.
    static inline CGAL::cpp11::array<typename Traits::FT,2> static_compute(const typename Traits::Point_2 &first_vertex, const typename Traits::Point_2 &second_vertex, const typename Traits::Point_2 &query_point, const Traits &barycentric_traits = Traits())
    {
        return static_segment_coordinates_2(first_vertex, second_vertex, query_point, barycentric_traits);
    }

    /// @}

    /// \name Information Functions
    /// @{

    /// This function prints some information about the used segment and segment coordinates.
    void print_information(std::ostream &output_stream = std::cout) const
    {
        output_stream << std::endl << "INFORMATION: " << std::endl;

        output_stream << std::endl << "DATA STRUCTURE: " << std::endl << std::endl;
        output_stream << "The internal data structure is segment." << std::endl;

        output_stream << std::endl << "DEGENERACY: " << std::endl << std::endl;
        if(!equal_2(vertex_0, vertex_1)) output_stream << "This segment is not degenerate." << std::endl;
        else output_stream << "This segment is degenerate. The correct computation is not expected!" << std::endl;

        output_stream << std::endl << "TYPE OF COORDINATES: " << std::endl << std::endl;
        output_stream << "The coordinate functions to be computed are segment coordinates." << std::endl;

        output_stream << std::endl << "INFORMATION ABOUT COORDINATES: " << std::endl << std::endl;
        output_stream << "Segment coordinates can be computed exactly for an arbitrary point along the line supporting the used segment." << std::endl;
        output_stream << "A slight offset from the line is allowed." << std::endl;
    }

    /// @}

private:

    // Internal global variables.
    const Traits &barycentric_traits;

    const Point &vertex_0;
    const Point &vertex_1;

    Scalar b_first;
    Scalar opposite_scalar_product;

    typename Traits::Compute_scalar_product_2 scalar_product_2;
    typename Traits::Compute_squared_distance_2 squared_distance_2;
    typename Traits::Equal_2 equal_2;

    // Compute segment coordinates.
    template<class Iterator> 
        std::pair<Iterator, bool> segment_coordinates_2(const Point &query_point, Iterator &output)
    {   
        // Project point on the segment and compute the first coordinate.
        opposite_scalar_product = scalar_product_2(query_point - vertex_1, vertex_0 - vertex_1);
        b_first = opposite_scalar_product / squared_distance_2(vertex_0, vertex_1);

        // Compute the second coordinate, using the partition of unity property.
        *output = b_first;
        ++output;
        *output = Scalar(1) - b_first;

        // Output both coordinates.
        return std::make_pair(output, true);
    }

    // This is a static function that takes both vertices of a segment and computes segment coordinates at a given query point with respect to these vertices.
    static typename CGAL::cpp11::array<typename Traits::FT,2> static_segment_coordinates_2(const typename Traits::Point_2 &vertex_0, const typename Traits::Point_2 &vertex_1, const typename Traits::Point_2 &query_point, const Traits &barycentric_traits)
    {
        // Some predefined functions.
        typename Traits::Compute_scalar_product_2 scalar_product_2 = barycentric_traits.compute_scalar_product_2_object();
        typename Traits::Compute_squared_distance_2 squared_distance_2 = barycentric_traits.compute_squared_distance_2_object();

        // Number type.
        typedef typename Traits::FT Scalar;

        // Project point on the segment and compute the first coordinate.
        const Scalar opposite_scalar_product = scalar_product_2(query_point - vertex_1, vertex_0 - vertex_1);
        const Scalar b_first = opposite_scalar_product / squared_distance_2(vertex_0, vertex_1);

        // Return the CGAL::cpp11::array<FT,2> type of coordinates.
        return CGAL::make_array(b_first, Scalar(1) - b_first);
    }
};

} // namespace Barycentric_coordinates

} // namespace CGAL

#endif // CGAL_SEGMENT_COORDINATES_2_H