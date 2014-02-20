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
  \file Barycentric_coordinates_base_2.h
*/

#ifndef CGAL_BARYCENTRIC_COORDINATES_BASE_2_H
#define CGAL_BARYCENTRIC_COORDINATES_BASE_2_H

// CGAL headers.
#include <CGAL/Kernel/global_functions_2.h>

// Barycentric coordinates headers.
#include <CGAL/Segment_coordinates_2.h>

// CGAL namespace.
namespace CGAL {

// Barycentric coordinates namespace.
namespace Barycentric_coordinates {

// Examples: See User Manual here - http://doc.cgal.org/latest/Manual/index.html.

/*!
 * \ingroup PkgBarycentric_coordinates_2
 * The Barycentric_coordinates_base_2 is an abstract base class for
 * all the types of two-dimensional generalized barycentric coordinates.
 * This class is parameterized by `CGAL::Polygon_2` class and `Iterator` class.
 * The latter can be any class that fulfills the requirements for an STL iterator.
 *
 * \sa `Iterator`
 *
 */

// This class does not allow the user to use different iterators to output coordinates after
// it has been created. In order to do so, we need to move template parameter Iterator
// in declaration of the class's functions but it is impossible because we use virtual functions.
// Can we also use reference in front of iterator like Iterator &output when passing it to the function?
template<typename Polygon_2, typename Iterator> 
    class Barycentric_coordinates_base_2
{

public:

    /// \name Creation
    /// @{

    /// Creates an instance of Barycentric_coordinates_base_2 class for a provided polygon passed as a reference.
    /// \pre `number_of_polygon_vertices > 2`
    /// \pre `_polygon.is_simple()`
    Barycentric_coordinates_base_2(const Polygon_2 &_polygon) : 
        number_of_polygon_vertices(_polygon.size()),
        polygon(_polygon)
    {
        CGAL_precondition( number_of_polygon_vertices > 2 );

        CGAL_precondition( _polygon.is_simple() );
    }

    /// Distructor of the class.
    virtual ~Barycentric_coordinates_base_2() { }

    /// @}

    // Compute weights.

    /// \name Computation of the weight functions
    /// @{

    /// Computes generalized barycentric weights for any strictly interior
    /// query point with respect to all the vertices of the polygon.
    /// This function is overloaded within each of the derived
    /// classes for a particular type of the coordinate function.
    /// \pre `_polygon.bounded_side(query_point) == CGAL::ON_BOUNDED_SIDE`
    inline std::pair<Iterator, bool> compute_weights(const typename Polygon_2::Point_2 &query_point, Iterator output)
    {
        return weights_2(query_point, output);
    }

    /// @}

    // Compute coordinates only at vertex.

    /// \name Computation of the basis functions at the vertices (with index)
    /// @{

    /// Computes generalized barycentric coordinates for a query point, which coincides with one
    /// of the polygon's vertices, with beforehand known index.
    /// \pre `(0 <= index) && (index < number_of_polygon_vertices)`
    inline std::pair<Iterator, bool> compute_at_vertex(const int index, Iterator output) const
    {
        return coordinates_at_vertex_2(index, output);
    }

    /// @}

    // Compute coordinates only along edge.

    /// \name Computation of the basis functions along edges (with index)
    /// @{

    /// Computes generalized barycentric coordinates for a query point
    /// on the polygon's boundary with beforehand known index of the edge to which this point belongs.
    /// \pre `_polygon.bounded_side(query_point) == CGAL::ON_BOUNDARY`
    /// \pre `(0 <= index) && (index < number_of_polygon_vertices)`
    inline std::pair<Iterator, bool> compute_on_edge(const typename Polygon_2::Point_2 &query_point, const int index, Iterator output) const
    {    
        return coordinates_on_boundary_2(query_point, index, output);
    }

    /// @}

    // Compute coordinates.

    /// \name Computation of the basis functions at an arbitrary point
    /// @{

    /// Computes generalized barycentric coordinates for any query point with respect to all the vertices of the polygon.
    /// This function is overloaded within each of the derived classes for a particular type of the coordinate function.
    /// Different choices of `query_point_location` parameter are possible: `CGAL::BC::UNSPECIFIED_LOCATION` -
    /// default constant with automatic check for a location, `CGAL::BC::ON_BOUNDED_SIDE` - for a strictly
    /// interior query point, `CGAL::BC::ON_BOUNDARY` - for a query point on the boundary of the polygon, and
    /// `CGAL::ON_UNBOUNDED_SIDE` - for a strictly exterior query point. Another parameter is `type_of_algorithm`
    /// with the following possible constants: `CGAL::BC::PRECISE` - default slow algorithm, which is precise as
    /// much as possible and `CGAL::BC::FAST` - fast algorithm, which is less precise but much faster. For some
    /// coordinate functions both algorithms can be the same.
    inline std::pair<Iterator, bool> compute(const typename Polygon_2::Point_2 &query_point, Iterator output, Query_point_location query_point_location = UNSPECIFIED_LOCATION, Type_of_algorithm type_of_algorithm = PRECISE)
    {   
        return coordinates_2(query_point, output, query_point_location, type_of_algorithm);
    }

    /// @}

    // Information about computed coordinates.

    /// \name Information functions
    /// @{

    /// Print some information about currently used polygon and coordinate functions.
    /// This function is partly overloaded within each of the
    /// derived classes for a particular type of the coordinate function.
    void print_info() const
    {
        std::cout << std::endl << "INFORMATION: " << std::endl;

        std::cout << std::endl << "DATA STRUCTURE: " << std::endl << std::endl;
        std::cout << "The used data structure is polygon." << std::endl;

        std::cout << std::endl << "NUMBER OF VERTICES: " << std::endl << std::endl;
        if(number_of_polygon_vertices < 3) {
            std::cout << "Provided polygon has " << number_of_polygon_vertices << " vertices." << std::endl;
            std::cout << "Since number of vertices is less than 3, generalized barycentric coordinates cannot be computed!" << std::endl;
            std::cout << "Please use the CGAL::Barycentric_coordinates::Segment_coordinates_2 class!" << std::endl;
        } else {
            if(number_of_polygon_vertices == 3) {
                 std::cout << "Provided polygon has " << number_of_polygon_vertices << " vertices." << std::endl;
                 std::cout << "For triangles it is better to use CGAL::Barycentric_coordinates::Triangle_coordinates_2 class!" << std::endl;
            }
            else std::cout << "Provided polygon has " << number_of_polygon_vertices << " vertices." << std::endl;
        }

        std::cout << std::endl << "SIMPLICITY: " << std::endl << std::endl;
        if(polygon.is_simple()) std::cout << "Current polygon is simple." << std::endl;
        else std::cout << "Current polygon is not simple. The correct computation is not expected!" << std::endl;

        print_coordinates_info();
    }

    /// @}

protected:

    // Internal global variables.
    const int number_of_polygon_vertices;
    const Polygon_2 &polygon;

    // Check the type of the provided polygon - CONVEX, STRICTLY_CONVEX, or CONCAVE.
    Type_of_polygon type_of_polygon() const
    {
        // First, test polygon on convexity.
        if(polygon.is_convex()) {
            // Index of the last polygon's vertex.
            const int last = number_of_polygon_vertices - 1;

            // Test all the consequent triplets of the polygon's vertices on collinearity. 
            // In case we find at least one, return WEAKLY_CONVEX polygon.
            if(CGAL::collinear(polygon.vertex(last), polygon.vertex(0), polygon.vertex(1)))
                return WEAKLY_CONVEX;

            for(int i = 1; i < last; ++i) {
                if(CGAL::collinear(polygon.vertex(i - 1), polygon.vertex(i), polygon.vertex(i + 1)))
                    return WEAKLY_CONVEX;
            }

            if(CGAL::collinear(polygon.vertex(last - 1), polygon.vertex(last), polygon.vertex(0)))
                return WEAKLY_CONVEX;

            // Otherwise, return STRICTLY_CONVEX polygon.
            return STRICTLY_CONVEX;
        }

        // Otherwise, return CONCAVE polygon.
        return CONCAVE;
    }

private:

    // Compute weights.

    // WEIGHTS.

    // Compute weights on bounded side of the polygon - see precondition.
    inline std::pair<Iterator, bool> weights_2(const typename Polygon_2::Point_2 &query_point, Iterator output)
    {
        // This is the only global precondition on the computation of weights.
        CGAL_precondition( polygon.bounded_side(query_point) == CGAL::ON_BOUNDED_SIDE );

        return weights(query_point, output);
    }

    // This function for weights must be overloaded within each of the derived classes.
    virtual std::pair<Iterator, bool> weights(const typename Polygon_2::Point_2 &query_point, Iterator output) = 0;

    // Compute coordinates.

    // COORDINATES EVERYWHERE.

    // Compute coordinates at any point in the plane with beforehand specified location.
    std::pair<Iterator, bool> coordinates_2(const typename Polygon_2::Point_2 &query_point, Iterator output, const Query_point_location query_point_location, const Type_of_algorithm type_of_algorithm)
    {
        // Determine location of the current query point provided by the user.
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

        case AT_VERTEX:
            return coordinates_at_vertex_2(query_point, output);
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
    std::pair<Iterator, bool> coordinates_unspecified_2(const typename Polygon_2::Point_2 &query_point, Iterator output, const Type_of_algorithm type_of_algorithm)
    {
        // Determine global location of the current query point.
        switch(polygon.bounded_side(query_point))
        {
        case CGAL::ON_BOUNDED_SIDE:
            return coordinates_on_bounded_side_2(query_point, output, type_of_algorithm);
            break;

        case CGAL::ON_BOUNDARY:
        {
            int index = -1;
            if(is_query_point_at_vertex(query_point, index)) return coordinates_at_vertex_2(index, output);
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

    // Compute coordinates on bounded side of the polygon - precise or fast.
    inline std::pair<Iterator, bool> coordinates_on_bounded_side_2(const typename Polygon_2::Point_2 &query_point, Iterator output, const Type_of_algorithm type_of_algorithm)
    {
        CGAL_precondition( polygon.bounded_side(query_point) == CGAL::ON_BOUNDED_SIDE );

        // Choose algorithm to compute coordinates on bounded side of the polygon.
        switch(type_of_algorithm)
        {
        case PRECISE:
            return coordinates_on_bounded_side_precise(query_point, output);
            break;

        case FAST:
            return coordinates_on_bounded_side_fast(query_point, output);
            break;
        }

        // Pointer cannot be here. Something went wrong.
        const bool type_of_algorithm_failure = true;
        CGAL_postcondition( !type_of_algorithm_failure );
        return std::make_pair(output, !type_of_algorithm_failure);
    }

    // This is a function to compute coordinates on bounded side of the polygon with maximum possible precision. 
    // It must be overloaded within each of the derived classes.
    virtual std::pair<Iterator, bool> coordinates_on_bounded_side_precise(const typename Polygon_2::Point_2 &query_point, Iterator output) = 0;

    // This is a function to compute coordinates on bounded side of the polygon with maximum possible speed. 
    // It must be overloaded within each of the derived classes.
    virtual std::pair<Iterator, bool> coordinates_on_bounded_side_fast(const typename Polygon_2::Point_2 &query_point, Iterator output) = 0;

    // COORDINATES ON BOUNDARY.

    // Compute coordinates along boundary of the polygon with beforehand known index of the edge to which current query point belongs.
    std::pair<Iterator, bool> coordinates_on_boundary_2(const typename Polygon_2::Point_2 &query_point, const int index, Iterator output) const
    {
        CGAL_precondition( polygon.bounded_side(query_point) == CGAL::ON_BOUNDARY ); 
        CGAL_precondition(  (0 <= index) && (index < number_of_polygon_vertices)  );

        // Some convenient typedefs.
        typedef typename Polygon_2::Segment_2 Segment;
        typedef typename Polygon_2::FT        Scalar;

        // Index of the last polygon's vertex.
        const int last = number_of_polygon_vertices - 1;

        // If the query point is on the last edge call coordinates_on_last_edge_2() function.
        if(index == last) return coordinates_on_last_edge_2(query_point, last, output);
        else {
            // Otherwise, all the coordinates are zero apart those for the chosen edge with index = `index`.
            int i;
            for(i = 0; i < index; ++i) {
                *output = Scalar(0);
                ++output;
            }

            // Compute Segment coordinates along the chosen edge with index = `index`.
            Segment segment = polygon.edge(index);
            Segment_coordinates_2<Segment, Iterator> segment_coordinates(segment);
            std::pair<Iterator, bool> result = segment_coordinates.compute(query_point, output);
            ++output;

            for(i = index + 1; i < last; ++i) {
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

    // Compute coordinates along boundary of the polygon without beforehand known index of the edge to which current query point belongs.
    std::pair<Iterator, bool> coordinates_on_boundary_2(const typename Polygon_2::Point_2 &query_point, Iterator output) const
    {
        CGAL_precondition( polygon.bounded_side(query_point) == CGAL::ON_BOUNDARY );

        // Some convenient typedefs.
        typedef typename Polygon_2::Segment_2 Segment;
        typedef typename Polygon_2::FT        Scalar;

        // Index of the last polygon's vertex.
        const int last = number_of_polygon_vertices - 1;

        // If the query point is on the last edge call coordinates_on_last_edge_2() function.
        if(polygon.edge(last).has_on(query_point)) return coordinates_on_last_edge_2(query_point, last, output);
        else {
            // Otherwise, all the coordinates are zero apart those for the edge with query point.
            int index;
            bool status = false;
            for(index = 0; index < last; ++index) {
                if(polygon.edge(index).has_on(query_point)) {

                    // Compute Segment coordinates along the edge with query point.
                    Segment segment = polygon.edge(index);
                    Segment_coordinates_2<Segment, Iterator> segment_coordinates(segment);
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

    // Compute coordinates for query point lying on the last edge of the polygon.
    std::pair<Iterator, bool> coordinates_on_last_edge_2(const typename Polygon_2::Point_2 &query_point, const int last, Iterator output) const
    {
        // Some convenient typedefs.
        typedef typename Polygon_2::Segment_2 Segment;
        typedef typename Polygon_2::FT        Scalar;

        typedef std::vector<Scalar> Coordinate_vector;

        typedef std::back_insert_iterator<Coordinate_vector> Vector_insert_iterator;

        // Create a coordinate vector.
        Coordinate_vector coordinates;

        // Compute Segment coordinates along the last edge of the polygon.
        Segment segment = polygon.edge(last);
        Segment_coordinates_2<Segment, Vector_insert_iterator> segment_coordinates(segment);
        std::pair<Vector_insert_iterator, bool> result = segment_coordinates.compute(query_point, std::back_inserter(coordinates));

        // Store all the coordinate values. 
        // All values are zeros apart those corresponding to the first and last vertices of the polygon.
        *output = coordinates[1];
        ++output;
        for(int i = 1; i < last; ++i) {
            *output = Scalar(0);
            ++output;
        }
        *output = coordinates[0];

        // Return computed coordinates.
        return std::make_pair(output, result.second);
    }

    // COORDINATES AT VERTEX.

    // Compute coordinates for a query point lying at on of the polygon's vertices with beforehand known vertex's index.
    std::pair<Iterator, bool> coordinates_at_vertex_2(const int index, Iterator output) const
    {
        CGAL_precondition( (0 <= index) && (index < number_of_polygon_vertices) );

        // Some convenient typedefs.
        typedef typename Polygon_2::FT Scalar;

        // All coordinate values are zeros apart that one with index = `index`, which is one.
        int i;
        for(i = 0; i < index; ++i) {
            *output = Scalar(0);
            ++output;
        }

        *output = Scalar(1);
        ++output;

        for(i = index + 1; i < number_of_polygon_vertices; ++i) {
            *output = Scalar(0);
            ++output;
        }

        // Return coordinates.
        return std::make_pair(output, true);
    }

    // Compute coordinates for a query point lying at on of the polygon's vertices without beforehand known vertex's index.
    std::pair<Iterator, bool> coordinates_at_vertex_2(const typename Polygon_2::Point_2 &query_point, Iterator output) const
    {
        int index = -1;
        CGAL_precondition( is_query_point_at_vertex(query_point, index) == true );

        // Some convenient typedefs.
        typedef typename Polygon_2::FT Scalar;

        // All coordinate values are zeros apart that one corresponding to the vertex coinciding with current query point, which is one.
        bool coordinates_at_vertex_failure = true;
        for(index = 0; index < number_of_polygon_vertices; ++index) {
            if(is_query_point_equal_to_vertex(query_point, index)) {

                // We have just found vertex coinciding with current query point.
                *output = Scalar(1);
                ++output;
                coordinates_at_vertex_failure = false;
                break;

            }
            else {
                *output = Scalar(0);
                ++output;
            }
        }

        for(int i = index + 1; i < number_of_polygon_vertices; ++i) {
            *output = Scalar(0);
            ++output;
        }

        // Return coordinates.
        CGAL_postcondition( !coordinates_at_vertex_failure );
        return std::make_pair(output, !coordinates_at_vertex_failure);
    }

    // COORDINATES ON UNBOUNDED SIDE.

    // Compute coordinates on unbounded side of the polygon - precise or fast.
    inline std::pair<Iterator, bool> coordinates_on_unbounded_side_2(const typename Polygon_2::Point_2 &query_point, Iterator output, const Type_of_algorithm type_of_algorithm)
    {
        CGAL_precondition( polygon.bounded_side(query_point) == CGAL::ON_UNBOUNDED_SIDE );

        // Choose algorithm to compute coordinates on unbounded side of the polygon.
        switch(type_of_algorithm)
        {
        case PRECISE:
            return coordinates_on_unbounded_side_precise(query_point, output);
            break;

        case FAST:
            return coordinates_on_unbounded_side_fast(query_point, output);
            break;
        }

        // Pointer cannot be here. Something went wrong.
        const bool type_of_algorithm_failure = true;
        CGAL_postcondition( !type_of_algorithm_failure );
        return std::make_pair(output, !type_of_algorithm_failure);
    }

    // This is a function to compute coordinates on unbounded side of the polygon with maximum possible precision. 
    // It must be overloaded within each of the derived classes.
    virtual std::pair<Iterator, bool> coordinates_on_unbounded_side_precise(const typename Polygon_2::Point_2 &query_point, Iterator output) = 0;

    // This is a function to compute coordinates on unbounded side of the polygon with maximum possible speed. 
    // It must be overloaded within each of the derived classes.
    virtual std::pair<Iterator, bool> coordinates_on_unbounded_side_fast(const typename Polygon_2::Point_2 &query_point, Iterator output) = 0;

    // HELP FUNCTIONS.

    // Test if the current query point is at vertex or not.
    inline bool is_query_point_at_vertex(const typename Polygon_2::Point_2 &query_point, int &index) const
    {
        for(index = 0; index < number_of_polygon_vertices; ++index)
            if(is_query_point_equal_to_vertex(query_point, index)) return true;

        return false;
    }

    // Test if the current query point is equal to the vertex with index = `index` or not.
    inline bool is_query_point_equal_to_vertex(const typename Polygon_2::Point_2 &query_point, const int index) const
    {
        if(query_point == polygon.vertex(index)) return true;
        return false;
    }

    // This is a function to print some info about currently computed coordinate functions.
    // It must be overloaded within each of the derived classes.
    virtual void print_coordinates_info() const = 0;
};

} // namespace Barycentric_coordinates

} // namespace CGAL

#endif // CGAL_BARYCENTRIC_COORDINATES_BASE_2_H