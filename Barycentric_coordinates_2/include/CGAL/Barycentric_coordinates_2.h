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
  \file Barycentric_coordinates_2.h
*/

#ifndef CGAL_BARYCENTRIC_COORDINATES_2_H
#define CGAL_BARYCENTRIC_COORDINATES_2_H

// CGAL headers.
#include <CGAL/Segment_2.h>
#include <CGAL/Polygon_2_algorithms.h>

// Barycentric coordinates headers.
#include <CGAL/Segment_coordinates_2.h> 
#include <CGAL/barycentric_enum.h>

// CGAL namespace.
namespace CGAL {

// Barycentric coordinates namespace.
namespace Barycentric_coordinates {

// Examples: see the User Manual here - http://doc.cgal.org/latest/Manual/index.html.

/*!
 * \ingroup PkgBarycentric_coordinates_2
 * The class Barycentric_coordinates_2 is the base class for all two-dimensional generalized barycentric coordinates. 
 * This class implements the behaviour of coordinate functions along the polygon's boundary and provides a common interface for all coordinate classes.
 * This class is parameterized by the the class `InputIterator`, a coordinate class `Coordinate_2`, and a traits class `Traits`.
 *
 * \pre The vertices of a polygon must be ordered.

 \cgalHeading{Requirements}

 <OL>
 <LI> `Traits` class must contain the following subset of types:
 <UL>
 <LI> `Traits::Point_2` - the type of a point used internally in the class `Barycentric_coordinates_2`;
 <LI> `Traits::Point_d` - the user-defined type of a point;
 <LI> `Traits::K` - the used kernel;
 <LI> `Traits::FT` - the type of a coordinate value;
 </UL>
 <LI> The value type of `InputIterator` is equivalent to the type `Traits::Point_d` that can be any user-defined type of a polygon's vertex and query point.
 <LI> The coordinate function of the type `Coordinate_2` must be one of the following classes: `Wachspress_coordinates_2`, `Mean_value_coordinates_2`, or `Discrete_harmonic_coordinates_2`.
 <LI> The traits class must provide a function `project()` that projects a point and sequence of points from the type `Traits::Point_d` to the type `Traits::Point_2`, which is equivalent to the type `CGAL::Point_2`, used internally in all coordinate classes.
 </OL>

 */

template<class InputIterator, class Coordinate_2, class Traits> 
    class Barycentric_coordinates_2
{

public:

    /// \name Types
    /// @{

    /// Number type.
    typedef typename Traits::FT      Scalar;

    /// Type of the used kernel.
    typedef typename Traits::K       Kernel;

    /// Type of a general point.
    typedef typename Traits::Point_d Point_d;

    /// Type of 2D point.
    typedef typename Traits::Point_2 Point_2;

    /// @}

    /// \name Creation
    /// @{

    /// Creates the base class `Barycentric_coordinates_2` that implements the behaviour along the polygon's boundary for a chosen coordinate function with respect to the polygon given by a range of vertices `[first_vertex, last_vertex]`.
    /// \pre Number of the polygon's vertices > 2.
    /// \pre The provided polygon is simple.
    Barycentric_coordinates_2(const InputIterator &first_vertex, const InputIterator &last_vertex) :
        barycentric_traits(Traits()),
        vertex(barycentric_traits.project(first_vertex, last_vertex)),
        number_of_vertices(vertex.size()),
        coordinate(Coordinate_2(vertex))
    {
        CGAL_precondition( int(number_of_vertices) > 2 );

        CGAL_precondition( CGAL::is_simple_2(vertex.begin(), vertex.end(), Kernel()) );
    }

    /// @}

    /// \name Computation of weight functions
    /// @{

    /// Computes generalized barycentric weights for any strictly interior query point with respect to all the vertices of the polygon.
    /// This function accepts any STL like iterator, which complies with the `Iterator` concept.
    ///
    /// \pre The provided query point belongs to the polygon's interior excluding the boundary.
    template<typename Iterator>
        inline std::pair<Iterator, bool> compute_weights(const Point_d &query_point, Iterator output)
    {
        return weights_2(barycentric_traits.project(query_point), output);
    }

    /// Computes generalized barycentric weights for any strictly interior query point with respect to all the vertices of the polygon.
    /// This function accepts a container of the type <a href="http://en.cppreference.com/w/cpp/container/vector">`std::vector`</a> 
    /// and returns an iterator of the type <a href="http://en.cppreference.com/w/cpp/iterator/back_insert_iterator">`std::back_insert_iterator`</a>
    /// that is placed past-the-end of the resulting sequence of weight values.
    ///
    /// \pre The provided query point belongs to the polygon's interior excluding the boundary.
    inline std::pair<std::back_insert_iterator<std::vector<Scalar> >, bool> compute_weights(const Point_d &query_point, std::vector<Scalar> &output_vector)
    {
        output_vector.reserve(output_vector.size() + number_of_vertices); 
        typedef typename std::back_insert_iterator<std::vector<Scalar> > Iterator;
        Iterator output = std::back_inserter(output_vector);
        return weights_2(barycentric_traits.project(query_point), output);
    }

    /// @}

    /// \name Computation of basis functions at the vertices (with index)
    /// @{

    /// Computes generalized barycentric coordinates for a query point, which coincides with one of the polygon's vertices, with known index.
    /// This function accepts any STL like iterator, which complies with the `Iterator` concept.
    ///
    /// \pre (0 <= index) && (index < number of the polygon's vertices).
    template<typename Iterator>
        inline std::pair<Iterator, bool> compute_on_vertex(const int index, Iterator output) const
    {
        return coordinates_on_vertex_2(index, output);
    }

    /// Computes generalized barycentric coordinates for a query point, which coincides with one of the polygon's vertices, with known index.
    /// This function accepts a container of the type <a href="http://en.cppreference.com/w/cpp/container/vector">`std::vector`</a> 
    /// and returns an iterator of the type <a href="http://en.cppreference.com/w/cpp/iterator/back_insert_iterator">`std::back_insert_iterator`</a>
    /// that is placed past-the-end of the resulting sequence of coordinate values.
    ///
    /// \pre (0 <= index) && (index < number of the polygon's vertices).
    inline std::pair<std::back_insert_iterator<std::vector<Scalar> >, bool> compute_on_vertex(const int index, std::vector<Scalar> &output_vector) const
    {
        output_vector.reserve(output_vector.size() + number_of_vertices); 
        typedef typename std::back_insert_iterator<std::vector<Scalar> > Iterator;
        Iterator output = std::back_inserter(output_vector);
        return coordinates_on_vertex_2(index, output);
    }

    /// @}

    /// \name Computation of basis functions along edges (with index)
    /// @{

    /// Computes generalized barycentric coordinates for a query point on the polygon's boundary with known index of the edge to which this point belongs.
    /// This function accepts any STL like iterator, which complies with the `Iterator` concept.
    ///
    /// \pre The provided query point belongs to the polygon's boundary.
    /// \pre (0 <= index) && (index < number of the polygon's vertices).
    template<typename Iterator>
        inline std::pair<Iterator, bool> compute_on_edge(const Point_d &query_point, const int index, Iterator output) const
    {    
        return coordinates_on_boundary_2(barycentric_traits.project(query_point), index, output);
    }

    /// Computes generalized barycentric coordinates for a query point on the polygon's boundary with known index of the edge to which this point belongs.
    /// This function accepts a container of the type <a href="http://en.cppreference.com/w/cpp/container/vector">`std::vector`</a> 
    /// and returns an iterator of the type <a href="http://en.cppreference.com/w/cpp/iterator/back_insert_iterator">`std::back_insert_iterator`</a>
    /// that is placed past-the-end of the resulting sequence of coordinate values.
        ///
    /// \pre The provided query point belongs to the polygon's boundary.
    /// \pre (0 <= index) && (index < number of the polygon's vertices).
    inline std::pair<std::back_insert_iterator<std::vector<Scalar> >, bool> compute_on_edge(const Point_d &query_point, const int index, std::vector<Scalar> &output_vector) const
    {
        output_vector.reserve(output_vector.size() + number_of_vertices); 
        typedef typename std::back_insert_iterator<std::vector<Scalar> > Iterator;
        Iterator output = std::back_inserter(output_vector);
        return coordinates_on_boundary_2(barycentric_traits.project(query_point), index, output);
    }

    /// @}

    /// \name Computation of basis functions at an arbitrary point
    /// @{

    /// Computes generalized barycentric coordinates for any query point in the plane with respect to all the vertices of the polygon.
    /// This function accepts any STL like iterator, which complies with the `Iterator` concept.
    ///
    /// Different choices of the parameter `query_point_location` are possible: 
    /// `CGAL::Barycentric_coordinates::UNSPECIFIED_LOCATION` - default constant with automatic check for a location, 
    /// `CGAL::Barycentric_coordinates::ON_BOUNDED_SIDE` - for a strictly interior query point, 
    /// `CGAL::Barycentric_coordinates::ON_BOUNDARY` - for a query point on the boundary of the polygon,
    /// `CGAL::Barycentric_coordinates::ON_VERTEX` - for a query point at one of the polygon's vertices, and
    /// `CGAL::Barycentric_coordinates::ON_UNBOUNDED_SIDE` - for a strictly exterior query point.
    ///
    /// Another parameter is `type_of_algorithm` with the following possible constants: 
    /// `CGAL::Barycentric_coordinates::PRECISE` - default slow algorithm, which is as precise as possible and
    /// `CGAL::Barycentric_coordinates::FAST` - fast algorithm, which is less precise but much faster.
    template<typename Iterator>
        inline std::pair<Iterator, bool> compute(const Point_d &query_point, Iterator output, Query_point_location query_point_location = UNSPECIFIED_LOCATION, Type_of_algorithm type_of_algorithm = PRECISE)
    {   
        return coordinates_2(barycentric_traits.project(query_point), output, query_point_location, type_of_algorithm);
    }

    /// Computes generalized barycentric coordinates for any query point in the plane with respect to all the vertices of the polygon.
    /// This function accepts a container of the type <a href="http://en.cppreference.com/w/cpp/container/vector">`std::vector`</a> 
    /// and returns an iterator of the type <a href="http://en.cppreference.com/w/cpp/iterator/back_insert_iterator">`std::back_insert_iterator`</a>
    /// that is placed past-the-end of the resulting sequence of coordinate values.
    inline std::pair<std::back_insert_iterator<std::vector<Scalar> >, bool> compute(const Point_d &query_point, std::vector<Scalar> &output_vector, Query_point_location query_point_location = UNSPECIFIED_LOCATION, Type_of_algorithm type_of_algorithm = PRECISE)
    {
        output_vector.reserve(output_vector.size() + number_of_vertices); 
        typedef typename std::back_insert_iterator<std::vector<Scalar> > Iterator;
        Iterator output = std::back_inserter(output_vector);
        return coordinates_2(barycentric_traits.project(query_point), output, query_point_location, type_of_algorithm);
    }

    /// @}

    /// \name Information functions
    /// @{

    /// Print some information about the used polygon and coordinate functions.
    void print_information(std::ostream &output_stream = std::cout) const
    {
        output_stream << std::endl << "INFORMATION: " << std::endl;

        output_stream << std::endl << "DATA STRUCTURE: " << std::endl << std::endl;
        output_stream << "The internal data structure is polygon." << std::endl;

        output_stream << std::endl << "NUMBER OF VERTICES: " << std::endl << std::endl;
        if(int(number_of_vertices) < 3) {
            output_stream << "This polygon has " << number_of_vertices << " vertices." << std::endl;
            output_stream << "Since number of vertices is less than 3, generalized barycentric coordinates cannot be computed!" << std::endl;
            output_stream << "Please use the class CGAL::Barycentric_coordinates::Segment_coordinates_2!" << std::endl;
        } else {
            if(int(number_of_vertices) == 3) {
                 output_stream << "This polygon has " << number_of_vertices << " vertices." << std::endl;
                 output_stream << "For triangles it is better to use the class CGAL::Barycentric_coordinates::Triangle_coordinates_2!" << std::endl;
            }
            else output_stream << "This polygon has " << number_of_vertices << " vertices." << std::endl;
        }

        output_stream << std::endl << "SIMPLICITY: " << std::endl << std::endl;
        if(CGAL::is_simple_2(vertex.begin(), vertex.end(), Kernel())) output_stream << "This polygon is simple." << std::endl;
        else output_stream << "This polygon is not simple. The correct computation is not expected!" << std::endl;

        coordinate.print_coordinates_information(output_stream);
    }

    /// @}

private:

    // Some convenient typedefs.
    typedef typename std::vector<Point_2> Point_vector;

    // Internal global variables.
    const Traits barycentric_traits;

    const Point_vector vertex;

    const size_t number_of_vertices;

    Coordinate_2 coordinate; 

    // WEIGHTS.

    // Compute weights on the bounded side of the polygon - see precondition.
    template<typename Iterator>
        inline std::pair<Iterator, bool> weights_2(const Point_2 &query_point, Iterator &output)
    {
        // This is the only global precondition on the computation of weights.
        CGAL_precondition( CGAL::bounded_side_2(vertex.begin(), vertex.end(), query_point, Kernel()) == CGAL::ON_BOUNDED_SIDE );

        return coordinate.weights(query_point, output);
    }

    // COORDINATES EVERYWHERE.

    // Compute coordinates at any point in the plane.
    template<typename Iterator>
        std::pair<Iterator, bool> coordinates_2(const Point_2 &query_point, Iterator &output, const Query_point_location query_point_location, const Type_of_algorithm type_of_algorithm)
    {
        // Determine a location of the current query point provided by the user.
        switch(query_point_location)
        {
        case UNSPECIFIED_LOCATION:
            return coordinates_unspecified_2(query_point, output, type_of_algorithm);
            break;

        case ON_BOUNDED_SIDE:
            return coordinates_on_bounded_side_2(query_point, output, type_of_algorithm);
            break;

        case ON_BOUNDARY:
            return coordinates_on_boundary_2(query_point, output);
            break;

        case ON_VERTEX:
            return coordinates_on_vertex_2(query_point, output);
            break;

        case ON_UNBOUNDED_SIDE:
            return coordinates_on_unbounded_side_2(query_point, output, type_of_algorithm);
            break;
        }

        // Pointer cannot be here. Something went wrong.
        const bool query_point_location_failure = true;
        CGAL_postcondition( !query_point_location_failure );
        return std::make_pair(output, !query_point_location_failure);
    }

    // Compute coordinates at any point in the plane with unspecified location.
    template<typename Iterator>
        std::pair<Iterator, bool> coordinates_unspecified_2(const Point_2 &query_point, Iterator &output, const Type_of_algorithm type_of_algorithm)
    {
        // Determine a global location of the current query point.
        switch(CGAL::bounded_side_2(vertex.begin(), vertex.end(), query_point, Kernel()))
        {
        case CGAL::ON_BOUNDED_SIDE:
            return coordinates_on_bounded_side_2(query_point, output, type_of_algorithm);
            break;

        case CGAL::ON_BOUNDARY:
        {
            int index = -1;
            if(is_query_point_at_vertex(query_point, index)) return coordinates_on_vertex_2(index, output);
            else return coordinates_on_boundary_2(query_point, output);
        }
            break;

        case CGAL::ON_UNBOUNDED_SIDE:
            return coordinates_on_unbounded_side_2(query_point, output, type_of_algorithm);
            break;
        }

        // Pointer cannot be here. Something went wrong.
        const bool query_point_location_failure = true;
        CGAL_postcondition( !query_point_location_failure );
        return std::make_pair(output, !query_point_location_failure);
    }

    // COORDINATES ON BOUNDED SIDE.

    // Compute coordinates on the bounded side of the polygon - precise or fast.
    template<typename Iterator>
        inline std::pair<Iterator, bool> coordinates_on_bounded_side_2(const Point_2 &query_point, Iterator &output, const Type_of_algorithm type_of_algorithm)
    {
        CGAL_precondition( CGAL::bounded_side_2(vertex.begin(), vertex.end(), query_point, Kernel()) == CGAL::ON_BOUNDED_SIDE );

        // Choose an algorithm to compute coordinates on the bounded side of the polygon.
        switch(type_of_algorithm)
        {
        case PRECISE:
            return coordinate.coordinates_on_bounded_side_precise(query_point, output);
            break;

        case FAST:
            return coordinate.coordinates_on_bounded_side_fast(query_point, output);
            break;
        }

        // Pointer cannot be here. Something went wrong.
        const bool type_of_algorithm_failure = true;
        CGAL_postcondition( !type_of_algorithm_failure );
        return std::make_pair(output, !type_of_algorithm_failure);
    }

    // COORDINATES ON BOUNDARY.

    // Compute coordinates along the boundary of the polygon with beforehand known index of the edge to which the query point belongs.
    template<typename Iterator>
        std::pair<Iterator, bool> coordinates_on_boundary_2(const Point_2 &query_point, const int index, Iterator &output) const
    {
        CGAL_precondition( CGAL::bounded_side_2(vertex.begin(), vertex.end(), query_point, Kernel()) == CGAL::ON_BOUNDARY ); 
        CGAL_precondition( (0 <= index) && (index < int(number_of_vertices)) );

        // Index of the last polygon's vertex.
        const int last = int(number_of_vertices) - 1;

        // If the query point is on the last edge call the function coordinates_on_last_edge_2().
        if(index == last) return coordinates_on_last_edge_2(query_point, last, output);
        else {
            // Otherwise, all the coordinates are zero apart those for the chosen edge with index = `index`.
            for(int i = 0; i < index; ++i) {
                *output = Scalar(0);
                ++output;
            }

            // Compute segment coordinates along the chosen edge with index = `index`.
            typedef CGAL::Barycentric_coordinates::Barycentric_traits_2<Kernel> B_traits;
            Segment_coordinates_2<B_traits> segment_coordinates(vertex[index], vertex[index+1]);
            std::pair<Iterator, bool> result = segment_coordinates.compute(query_point, output);
            ++output;

            for(int i = index + 1; i < last; ++i) {
                *output = Scalar(0);
                ++output;
            }

            // Return coordinates.
            return std::make_pair(output, result.second);
        }

        // Pointer cannot be here. Something went wrong.
        const bool coordinates_on_boundary_failure = true;
        CGAL_postcondition( !coordinates_on_boundary_failure );
        return std::make_pair(output, !coordinates_on_boundary_failure);
    }

    // Compute coordinates along the boundary of the polygon without beforehand known index of the edge to which current query point belongs.
    template<typename Iterator>
        std::pair<Iterator, bool> coordinates_on_boundary_2(const Point_2 &query_point, Iterator &output) const
    {
        CGAL_precondition( CGAL::bounded_side_2(vertex.begin(), vertex.end(), query_point, Kernel()) == CGAL::ON_BOUNDARY );

        // Some convenient typedefs.
        typedef CGAL::Segment_2<Kernel> Segment;

        // Index of the last polygon's vertex.
        const int last = int(number_of_vertices) - 1;

        // If the query point is on the last edge call the function coordinates_on_last_edge_2().
        if(Segment(vertex[last],vertex[0]).has_on(query_point)) return coordinates_on_last_edge_2(query_point, last, output);
        else {
            // Otherwise, all the coordinates are zeros apart from those for the edge with the query point.
            int index;
            bool status = false;
            for(index = 0; index < last; ++index) {
                if(Segment(vertex[index],vertex[index+1]).has_on(query_point)) {

                    // Compute segment coordinates along the edge with the query point.
                    typedef CGAL::Barycentric_coordinates::Barycentric_traits_2<Kernel> B_traits;
                    Segment_coordinates_2<B_traits> segment_coordinates(vertex[index], vertex[index+1]);
                    std::pair<Iterator, bool> result = segment_coordinates.compute(query_point, output);
                    status = result.second;
                    ++output;
                    break;

                } else {
                    *output = Scalar(0);
                    ++output;
                }
            }

            for(int i = index + 1; i < last; ++i) {
                *output = Scalar(0);
                ++output;
            }

            // Return coordinates.
            return std::make_pair(output, status);
        }

        // Pointer cannot be here. Something went wrong.
        const bool coordinates_on_boundary_failure = true;
        CGAL_postcondition( !coordinates_on_boundary_failure );
        return std::make_pair(output, !coordinates_on_boundary_failure);
    }

    // Compute coordinates for a query point lying on the last edge of the polygon.
    template<typename Iterator>
        std::pair<Iterator, bool> coordinates_on_last_edge_2(const Point_2 &query_point, const int last, Iterator &output) const
    {
        // Some convenient typedefs.
        typedef CGAL::Barycentric_coordinates::Barycentric_traits_2<Kernel> B_traits;

        typedef std::vector<Scalar> Coordinate_vector;

        typedef std::back_insert_iterator<Coordinate_vector> Vector_insert_iterator;

        // Create a coordinate vector.
        Coordinate_vector coordinate;
        coordinate.reserve(2);

        // Compute segment coordinates along the last edge of the polygon.
        Segment_coordinates_2<B_traits> segment_coordinates(vertex[last], vertex[0]);
        std::pair<Vector_insert_iterator, bool> result = segment_coordinates.compute(query_point, std::back_inserter(coordinate));

        // Store all the coordinate values. 
        // All the values are zeros apart from those corresponding to the first and last vertices of the polygon.
        *output = coordinate[1];
        ++output;
        for(int i = 1; i < last; ++i) {
            *output = Scalar(0);
            ++output;
        }
        *output = coordinate[0];

        // Return computed coordinates.
        return std::make_pair(output, result.second);
    }

    // COORDINATES AT VERTEX.

    // Compute coordinates for a query point lying at one of the polygon's vertices with beforehand known vertex's index.
    template<typename Iterator>
        std::pair<Iterator, bool> coordinates_on_vertex_2(const int index, Iterator &output) const
    {
        CGAL_precondition( (0 <= index) && (index < int(number_of_vertices)) );

        // All the coordinate values are zeros apart from that one with the index = `index`, which is one.
        for(int i = 0; i < index; ++i) {
            *output = Scalar(0);
            ++output;
        }

        *output = Scalar(1);
        ++output;

        for(int i = index + 1; i < int(number_of_vertices); ++i) {
            *output = Scalar(0);
            ++output;
        }

        // Return coordinates.
        return std::make_pair(output, true);
    }

    // Compute coordinates for a query point lying at one of the polygon's vertices without beforehand known vertex's index.
    template<typename Iterator>
        std::pair<Iterator, bool> coordinates_on_vertex_2(const Point_2 &query_point, Iterator &output) const
    {
        int index = -1;
        CGAL_precondition( is_query_point_at_vertex(query_point, index) == true );

        // All the coordinate values are zeros apart from that one corresponding to the vertex coinciding with the query point, which is one.
        bool coordinates_on_vertex_failure = true;
        for(index = 0; index < int(number_of_vertices); ++index) {
            if(is_query_point_equal_to_vertex(query_point, index)) {

                // We have just found vertex coinciding with the query point.
                *output = Scalar(1);
                ++output;
                coordinates_on_vertex_failure = false;
                break;

            }
            else {
                *output = Scalar(0);
                ++output;
            }
        }

        for(int i = index + 1; i < int(number_of_vertices); ++i) {
            *output = Scalar(0);
            ++output;
        }

        // Return coordinates.
        CGAL_postcondition( !coordinates_on_vertex_failure );
        return std::make_pair(output, !coordinates_on_vertex_failure);
    }

    // COORDINATES ON UNBOUNDED SIDE.

    // Compute coordinates on the unbounded side of the polygon - precise or fast.
    template<typename Iterator>
        inline std::pair<Iterator, bool> coordinates_on_unbounded_side_2(const Point_2 &query_point, Iterator &output, const Type_of_algorithm type_of_algorithm)
    {
        CGAL_precondition( CGAL::bounded_side_2(vertex.begin(), vertex.end(), query_point, Kernel()) == CGAL::ON_UNBOUNDED_SIDE );

        // Choose an algorithm to compute coordinates on the unbounded side of the polygon.
        switch(type_of_algorithm)
        {
        case PRECISE:
            return coordinate.coordinates_on_unbounded_side_precise(query_point, output);
            break;

        case FAST:
            return coordinate.coordinates_on_unbounded_side_fast(query_point, output);
            break;
        }

        // Pointer cannot be here. Something went wrong.
        const bool type_of_algorithm_failure = true;
        CGAL_postcondition( !type_of_algorithm_failure );
        return std::make_pair(output, !type_of_algorithm_failure);
    }

    // HELP FUNCTIONS.

    // Test if the query point is at vertex or not.
    inline bool is_query_point_at_vertex(const Point_2 &query_point, int &index) const
    {
        for(index = 0; index < int(number_of_vertices); ++index)
            if(is_query_point_equal_to_vertex(query_point, index)) return true;

        return false;
    }

    // Test if the query point is equal to the vertex with the index = `index` or not.
    inline bool is_query_point_equal_to_vertex(const Point_2 &query_point, const int index) const
    {
        if(query_point == vertex[index]) return true;
        return false;
    }
};

} // namespace Barycentric_coordinates

} // namespace CGAL

#endif // CGAL_BARYCENTRIC_COORDINATES_BASE_2_H