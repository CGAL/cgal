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
  \file Segment_coordinates_2.h
*/

#ifndef CGAL_SEGMENT_COORDINATES_2_H
#define CGAL_SEGMENT_COORDINATES_2_H

// CGAL headers.
#include <CGAL/assertions.h>
#include <CGAL/squared_distance_2.h>

// Barycentric coordinates headers.
#include <CGAL/barycentric_enum.h>

// CGAL namespace.
namespace CGAL {

// Barycentric coordinates namespace.
namespace Barycentric_coordinates {

// Examples: See User Manual here - http://doc.cgal.org/latest/Manual/index.html.

/*!
 * \ingroup PkgBarycentric_coordinates_2
 * The class Segment_coordinates_2 implements barycentric coordinates with respect to an arbitrary
 * non-degenerate segment along an arbitrary line in the plane. This class is parameterized by
 * `CGAL::Segment_2` class and `Iterator` class. The latter can be any class that fulfills
 * the requirements for an STL iterator.
 *
 * \sa `Iterator`
 *
 */

// This class does not allow the user to use different iterators to output coordinates after
// it has been created. In order to do so, we need to move template parameter Iterator in declaration of the function compute().
// Can we also use reference in front of iterator like OutputIterator &output when passing it to the function? 
template<typename Segment_2, typename Iterator> 
    class Segment_coordinates_2
{

public:

    // Creation.

    /// \name Types
    /// @{

    /// Type of 2D segment.
    typedef Segment_2                Segment;

    /// Type of the used kernel.
    typedef typename Segment::R      Kernel;

    /// Type of 2D point.
    typedef typename Kernel::Point_2 Point;

    /// Number type.
    typedef typename Kernel::FT      Scalar;

    /// Type of the used output iterator.
    typedef Iterator                 OutputIterator;

    /// @}

    /// \name Creation
    /// @{

    /// Creates an instance of Segment_coordinates_2 class for a provided segment passed as a reference.
    /// \pre `!_segment.is_degenerate()`
    Segment_coordinates_2(const Segment &_segment) : segment(_segment)
    {
        CGAL_precondition( !_segment.is_degenerate() );
    }

    /// @}

    // Coordinates.

    /// \name Computation of the basis functions
    /// @{

    /// Computes Segment barycentric coordinates for current query point with
    /// respect to both vertices of the segment.
    inline std::pair<OutputIterator, bool> compute(const Point &query_point, OutputIterator output)
    {
        return segment_coordinates_2(query_point, output);
    }

    /// @}

    // Information about computed coordinates.

    /// \name Information functions
    /// @{

    /// Print some information about currently used segment and Segment coordinates.
    void print_info() const
    {
        std::cout << std::endl << "INFORMATION: " << std::endl;

        std::cout << std::endl << "DATA STRUCTURE: " << std::endl << std::endl;
        std::cout << "The used data structure is segment." << std::endl;

        std::cout << std::endl << "DEGENERACY: " << std::endl << std::endl;
        if(!segment.is_degenerate()) std::cout << "Current segment is not degenerate." << std::endl;
        else std::cout << "Current segment is degenerate. The correct computation is not expected!" << std::endl;

        std::cout << std::endl << "TYPE OF COORDINATES: " << std::endl << std::endl;
        std::cout << "Currently computed coordinate functions are Segment coordinates." << std::endl;

        std::cout << std::endl << "INFORMATION ABOUT COORDINATES: " << std::endl << std::endl;
        std::cout << "Segment coordinates can be computed exactly for an arbitrary point along the line supporting the used segment." << std::endl;
        std::cout << "The slight offset from the line is allowed." << std::endl;
    }

    /// @}

private:

    // Internal global variables.
    const Segment &segment;

    Scalar b_first;
    Scalar opposite_scalar_product;

    // Compute Segment coordinates.
    // Can we somehow use a referenced output: &output?
    std::pair<OutputIterator, bool> segment_coordinates_2(const Point &query_point, OutputIterator output)
    {
        // Project point on the segment and compute the first coordinate.
        opposite_scalar_product = (query_point - segment.vertex(1))*(segment.vertex(0) - segment.vertex(1));
        b_first = opposite_scalar_product / CGAL::squared_distance(segment.vertex(0), segment.vertex(1));

        // Compute second coordinate using the partition of unity property.
        *output = b_first;
        ++output;
        *output = Scalar(1) - b_first;

        // Output coordinates.
        return std::make_pair(output, true);
    }
};

// Class Segment_coordinates_2 with particular std::back_insert_iterator instead of a general one.

/*!
 * \ingroup PkgBarycentric_coordinates_2
 * The class Seg_coordinates_2 implements barycentric coordinates with respect to an arbitrary
 * non-degenerate segment along an arbitrary line in the plane. This class is parameterized by
 * `CGAL::Segment_2` class and a Container class. The latter can be any class that fulfills
 * the requirements for <a href="http://en.cppreference.com/w/cpp/iterator/back_insert_iterator">`std::back_insert_iterator`</a>. 
 * It defaults to <a href="http://en.cppreference.com/w/cpp/container/vector">`std::vector`</a> container.
 *
 * \sa <a href="http://en.cppreference.com/w/cpp/iterator/back_insert_iterator">`std::back_insert_iterator`</a>
 *
 */

template<typename Segment_2, typename InputContainer = std::vector<typename Segment_2::R::FT> > 
    class Seg_coordinates_2 : public Segment_coordinates_2<Segment_2, std::back_insert_iterator<InputContainer> >
{

public:

    // Creation.

    /// \name Types
    /// @{

    /// Type of the used container.
    typedef InputContainer Container;

    /// Type of the base class.
    typedef Segment_coordinates_2<Segment_2, std::back_insert_iterator<Container> > Base;

    /// @}

    /// \name Creation
    /// @{

    /// Creates an instance of Segment_coordinates_2 class for a provided segment passed as a reference.
    /// The used iterator is <a href="http://en.cppreference.com/w/cpp/iterator/back_insert_iterator">`std::back_insert_iterator`</a>.
    Seg_coordinates_2(const Segment_2 &_segment) : Base(_segment) { }

    /// @}

    // Coordinates.

    /// \name Computation of the basis functions
    /// @{

    /// Computes Segment barycentric coordinates for current query point with
    /// respect to both vertices of the segment.
    inline std::pair<std::back_insert_iterator<Container>, bool> compute(const typename Segment_2::R::Point_2 &query_point, Container &container)
    {
        return Base::compute(query_point, std::back_inserter(container));
    }

    /// Static function, which gets a `CGAL::Segment_2` and returns a `CGAL::Point_2` type of coordinates.
    static inline typename Segment_2::R::Point_2 Compute(const Segment_2 &_segment, const typename Segment_2::R::Point_2 &query_point)
    {
        return static_segment_coordinates_2(_segment, query_point);
    }

    /// @}

private:

    // Static function, which gets a CGAL::Segment_2 and returns a CGAL::Point_2 with coordinates.
    static typename Segment_2::R::Point_2 static_segment_coordinates_2(const Segment_2 &_segment, const typename Segment_2::R::Point_2 &query_point)
    {
        // Point type.
        typedef typename Segment_2::R::Point_2 Point;

        // Segment coordinates type.
        typedef typename CGAL::Barycentric_coordinates::Seg_coordinates_2<Segment_2, Container> Segment_coordinates;

        // Instantiate Segment coordinates class for the segment defined above.
        Segment_coordinates segment_coordinates(_segment);

        // Create a container to store coordinates.
        Container coordinates;

        // Compute coordinates.
        segment_coordinates.compute(query_point, coordinates);

        // Return CGAL::Point_2 type of coordinates.
        return Point(coordinates[0], coordinates[1]);
    }
};

} // namespace Barycentric_coordinates

} // namespace CGAL

#endif // CGAL_SEGMENT_COORDINATES_2_H