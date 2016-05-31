// Copyright (c) 2013-2015  The University of Western Sydney, Australia.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Authors: Weisheng Si, Quincy Tse

#ifndef CGAL_LESS_BY_DIRECTION_2_H
#define CGAL_LESS_BY_DIRECTION_2_H

#include <iostream>
#include <cstdlib>
#include <utility>
#include <CGAL/Polynomial.h>
#include <CGAL/number_utils.h>
#include <CGAL/enum.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Aff_transformation_2.h>

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>

namespace CGAL {

/*  Function object that orders 2D graph vertex_descriptors based on the order
 *  induced by the direction D described in the book:
 *
 *  Giri Narasimhan and Michiel Smid, Chapter 4: Spanners based on the Theta graph, Geometric Spanner Networks,
 *  Cambridge University Press, 2007.
 *
 *  The ties are broken according to the direction of ccw90(D).
 *  The way of breaking ties in this functor is intended to prevent the overlapping of cone boundaries.
 *  As a result, a vertex on the cw boundary will be considered to be inside this cone,
 *  while a vertex on the ccw boundary will not.
 *
 *  This function object is implemented using the function `CGAL::compare_signed_distance_to_line_2()`,
 *  which orders two points according to their signed distance to a base line.
 *
 */
template <typename Kernel_, typename Graph_>
class  Less_by_direction_2 : public std::binary_function <typename Graph_::vertex_descriptor,
        typename Graph_::vertex_descriptor, bool> {

public:
    // typedef for C++11 - doesn't hurt to also have for C++98
    typedef typename Graph_::vertex_descriptor first_argument_type;
    typedef typename Graph_::vertex_descriptor second_argument_type;
    typedef bool     result_type;

    // typedef for Direction_2 and Line_2
    typedef typename Kernel_::Direction_2 Direction_2;
    typedef typename Kernel_::Line_2 Line_2;
    typedef typename Kernel_::Point_2 Point_2;
    typedef typename Kernel_::Aff_transformation_2 Transformation;

    // constructor
    Less_by_direction_2(const Graph_& g, const Direction_2& d)
        : graph(g), base_line(Point_2(0,0), d) {};

    bool operator() (const typename Graph_::vertex_descriptor& p,
                     const typename Graph_::vertex_descriptor& q) const {
        Comparison_result outcome;
        outcome = compare_signed_distance_to_line(base_line, graph[p], graph[q]);
        if (outcome == SMALLER)
            return true;
        else {
            if (outcome == LARGER)
                return false;
        }

        /* otherwise, outcome == CGAL::EQUAL, ties will be broken by a second order
         * according to the ccw90(base_line) direction.
         */
        // define a rotation of counter clockwise 90
        Transformation ccw90(0, -1, 1,  0);
        // rotate
        Line_2 ccw90_line = ccw90(base_line);
        outcome = compare_signed_distance_to_line(ccw90_line, graph[p], graph[q]);
        if (outcome == SMALLER)
            return true;
        else
            return false;
    }

private:
    const Graph_& graph;
    const Line_2 base_line;

};      // class Less_by_direction_2

}  // namespace CGAL


#endif
