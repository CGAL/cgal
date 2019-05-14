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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Authors: Weisheng Si, Quincy Tse, Frédérik Paradis

/*! \file Construct_theta_graph_2.h
 *
 * This header implements the functor for constructing Theta graphs.
 */

#ifndef CGAL_CONSTRUCT_THETA_GRAPH_2_H
#define CGAL_CONSTRUCT_THETA_GRAPH_2_H

#include <CGAL/license/Cone_spanners_2.h>


#include <iostream>
#include <cstdlib>
#include <utility>
#include <CGAL/Aff_transformation_2.h>
#include <CGAL/Compute_cone_boundaries_2.h>
#include <CGAL/Cone_spanners_2/Less_by_direction_2.h>
#include <CGAL/Cone_spanners_2/Plane_scan_tree.h>
#include <CGAL/Cone_spanners_enum_2.h>
#include <CGAL/tss.h>

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>

namespace CGAL {

/*! \ingroup PkgConeBasedSpanners

 \brief A template functor for constructing Theta graphs with a given set of 2D points and
         a given initial direction for the cone boundaries.

 \tparam Traits_  Must be either `CGAL::Exact_predicates_exact_constructions_kernel_with_root_of`
                  or `CGAL::Exact_predicates_inexact_constructions_kernel`.

 \tparam Graph_   The graph type to store the constructed cone based spanner.
                  It must be <A HREF="http://www.boost.org/libs/graph/doc/adjacency_list.html">`boost::adjacency_list`</A>
                  with `Traits_::Point_2` as `VertexProperties`.
 */
template <typename Traits_, typename Graph_>
class Construct_theta_graph_2 {

public:

    /*! the geometric traits class.  */
    typedef Traits_ Traits;

    /*! the specific type of `boost::adjacency_list`. */
    typedef Graph_                           Graph;

    /*! the point type */
    typedef typename Traits::Point_2                 Point_2;

    /*! the direction type */
    typedef typename Traits::Direction_2             Direction_2;

private:


    typedef typename Traits::Line_2                  Line_2;
    typedef typename Traits::Aff_transformation_2    Transformation;
    typedef Less_by_direction_2<Traits, Graph_>      Less_by_direction;

    /* Store the number of cones.  */
    unsigned int  cone_number;

    /* Store whether even, odd or all cones are selected to construct graph. */
    Cones_selected cones_choice;

    /* Store the directions of the rays dividing the plane. The initial direction will be
     * stored in rays[0].
     */
    std::vector<Direction_2>   rays;

public:
    /*! \brief Constructor.

     \param k     Number of cones to divide space into
     \param initial_direction  A direction denoting one of the rays dividing the
                   cones. This allows arbitary rotations of the rays that divide
                   the plane.  (default: positive x-axis)
     \param cones_selected  Indicates whether even, odd or all cones are
                   selected to construct graph.
     */
    Construct_theta_graph_2 (unsigned int k,
                             Direction_2 initial_direction = Direction_2(1,0),
                             Cones_selected cones_selected = ALL_CONES
                            ): cone_number(k), cones_choice(cones_selected), rays(std::vector<Direction_2>(k))

    {
        if (k<2) {
            std::cout << "The number of cones must be larger than 1!" << std::endl;
            CGAL_assertion(false);
        }

        /* Initialize a functor, specialization will happen here depending on the kernel type to
         compute the cone boundaries either exactly or inexactly */
        Compute_cone_boundaries_2<Traits> compute_cones;
        // compute the rays using the functor
        compute_cones(k, initial_direction, rays.begin());
    }

    /*!
     \brief Function operator to construct a Theta graph.

     \details For the details of this algorithm, please refer to the User Manual.

     \tparam  PointInputIterator an `InputIterator` with value type `Point_2`.
     \param[in] start An iterator pointing to the first vertex of the input.
     \param[in] end   An iterator pointing to the past-the-end location of the input.
     \param[out] g    The constructed graph object.
     */
    template <typename PointInputIterator>
    Graph_& operator()(const PointInputIterator& start,
                       const PointInputIterator& end,
                       Graph_& g) {
        // add vertices into the graph
        for (PointInputIterator curr = start; curr != end; ++curr) {
            g[boost::add_vertex(g)] = *curr;
        }

        unsigned int i;   // ray index of the cw ray
        unsigned int j;   // index of the ccw ray

        // add edges into the graph for every cone
        int new_start = cones_choice != ALL_CONES ? cones_choice : 0;
        int increment = cones_choice != ALL_CONES ? 2 : 1;
        for (i = new_start; i < cone_number; i += increment) {
            j = (i+1) % cone_number;
            add_edges_in_cone(rays[i], rays[j], g);
        }

        return g;
    }

    /*! \brief returns the number of cones.
     */
    unsigned int number_of_cones() const {
        return cone_number;
    }

    /*! \brief outputs the set of directions to the iterator `result`.

      \tparam DirectionOutputIterator  an `OutputIterator` with value type `Direction_2`.
       \return `result`
     */
    template<class DirectionOutputIterator>
    DirectionOutputIterator directions(DirectionOutputIterator result) {
        typename std::vector<Direction_2>::iterator it;
        for (it=rays.begin(); it!=rays.end(); it++) {
            *result++ = *it;
        }
        return result;
    }

protected:

    /* Construct edges in one cone bounded by two directions.

     \param cwBound      The direction of the clockwise boundary of the cone.
     \param ccwBound     The direction of the counter-clockwise boundary.
     \param g            The Theta graph to be built.
    */
    void add_edges_in_cone(const Direction_2& cwBound, const Direction_2& ccwBound, Graph_& g) {
        if (ccwBound == cwBound) {
            // Degenerate case,  not allowed.
            throw std::out_of_range("The cw boundary and the ccw boundary shouldn't be same!");
        }

        // Find angle bisector (requiring sqrt(), not exact)
        Line_2 cwLine(ORIGIN, cwBound);
        Line_2 ccwLine(ORIGIN, ccwBound);
        Direction_2 bisector_direction = bisector(cwLine, ccwLine).direction();

        // Rotational transformation of cw 90 degree
        CGAL_STATIC_THREAD_LOCAL_VARIABLE_4(Transformation, cw90, 0, 1, -1,  0);
        
        // Ordering
        // here D1 is the reverse of D1 in the book, we find this is easier to implement
        const Less_by_direction  orderD1 (g, ccwBound);
        const Less_by_direction  orderD2 (g, cwBound);
        const Less_by_direction  orderMid(g, cw90(bisector_direction));

        typename Graph_::vertex_iterator vit, ve;
        boost::tie(vit, ve) = boost::vertices(g);

        // Step 1: Sort S according to order induced by D1
        std::vector<typename Graph_::vertex_descriptor> S(vit, ve);
        std::sort(S.begin (), S.end (), orderD1);

        // Step 2: Initialise an empty set to store vertices sorted by orderD2
        typedef CGAL::ThetaDetail::Plane_scan_tree<typename Graph_::vertex_descriptor,
                typename Graph_::vertex_descriptor,
                Less_by_direction,
                Less_by_direction > PSTree;
        PSTree pst(orderD2, orderMid);

        // Step 3: visit S in orderD1
        //         insert '*it' into T
        //         find ri = T.minAbove(*it)
        //         add an edge
        for (typename std::vector<typename Graph_::vertex_descriptor>::const_iterator
                it = S.begin(); it != S.end(); ++it) {
            pst.add(*it, *it);
            const typename Graph_::vertex_descriptor *const ri = pst.minAbove(*it);
            if ( ri != NULL ) {
                typename Graph_::edge_descriptor existing_e;
                bool                    existing;
                // check whether the edge already exists
                boost::tie(existing_e, existing)=boost::edge(*it, *ri, g);
                if (!existing)
                    boost::add_edge(*it, *ri, g);
            }

        }  // end of for
    };     // end of add edges in cone

};      // class Construct_theta_graph_2


}  // namespace CGAL


#endif
