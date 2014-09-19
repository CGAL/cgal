// Copyright (c) 2013-2014  The University of Western Sydney, Australia.
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
// Authors: Quincy Tse, Weisheng Si

/** @file Theta_graph_2.h
 *
 * This header implements the class Theta_graph_2, the constructor of which 
 * builds a Theta graph on a set of given vertices.
 */

#ifndef CGAL_THETA_GRAPH_2_H
#define CGAL_THETA_GRAPH_2_H

#include <CGAL/Cone_spanners_2/_cxx0x_hack.h>

#include <iostream>
#include <cstdlib>
#include <set>
#ifndef GXX11
#include <functional>
#endif

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <CGAL/Cone_spanners_2.h>
#include <CGAL/Cone_spanners_2/Plane_Scan_Tree.h>

namespace CGAL {

	/**
        * \ingroup PkgConeBasedSpanners
        * @brief A derived class for constructing Theta graphs with a given set of 2D points.
	*
	*  Its base class is Cone_spanners_2. 
	*  Directed,undirected and bidirectional graphs are supported. For differences among these
	*  three types of graphs, please see the documentation of BGL.
	*/
	template <typename Kernel, 
		      typename Directedness=boost::undirectedS, 
			  typename EdgeProperty=boost::no_property>
	class Theta_graph_2 : public Cone_spanners_2<Kernel, Directedness, EdgeProperty>
	{
	public:
        typedef typename Kernel::Direction_2             Direction_2;
	    typedef typename Kernel::Point_2                 Point_2;
		typedef typename Kernel::Line_2                  Line_2;
		typedef typename Kernel::Aff_transformation_2    Transformation;

		typedef Cone_spanners_2<Kernel, Directedness, EdgeProperty>  base_class;
		typedef typename base_class::Graph  Graph;
		typedef typename base_class::vertex_smaller_2  vertex_smaller_2;

		/** @brief constructs a Theta graph object.
		*
		* @param k     Number of cones to divide space into
		* @param start An iterator pointing to the first point (vertex) in the graph.
		*              (default: nullptr)
		* @param end   An iterator pointing to the place that passes the last point.  (default: nullptr)
		* @param ray0  A direction denoting one of the rays deviding the
		*              cones. This allows arbitary rotations of the way space is
		*              divided.  (default: positive x-axis) 
		*/
#ifdef GXX11
	    template <typename PointInputIterator=Point_2*>
#else
	    template <typename PointInputIterator>  
#endif
		Theta_graph_2(const unsigned int k,
			        const PointInputIterator& start=nullptr, 
					const PointInputIterator& end=nullptr,
			        const Direction_2& ray0 = Direction_2(1,0)
					)
					: Cone_spanners_2<Kernel, Directedness, EdgeProperty>(k, start, end, ray0)
		{
			build_edges();
		}

		/** @brief copy constructor.
		*  @param x  another Theta_graph_2 object to copy from.
		*/
		Theta_graph_2 (const Theta_graph_2& x) 
			: Cone_spanners_2<Kernel, Directedness, EdgeProperty>(x) {}

		/** @brief This function implements the algorithm for adding edges to build the Theta graph.
		 * The algorithm implemented is described in 
		 * Giri Narasimhan and Michiel Smid, Chapter 4: Spanners based on the Theta graph, Geometric Spanner Networks,
		 * Cambridge University Press, 2007.
		 * This algorithm has the complexity of O(n*log(n)), which is optimal.
		 *
		 *  @return   the constructed graph object.
		 */
		virtual Graph& build_edges() {
			unsigned int i;   // ray index of the cw ray
		    unsigned int j;   // index of the ccw ray

			for (i = 0; i < this->num_cones; i++) {
			    j = (i+1) % this->num_cones;   
				add_edges_in_cone(this->rays[i], this->rays[j]);
			}
			return this->g;
		}

		/** @brief Construct edges bounded by two directions. 
		*
		* @param cwBound      The direction that bounds the cone on the clockwise
		*                      direction.
		* @param ccwBound     The direction that bounds the cone on the counter-clockwise
		*                      direction.
		*/
		void add_edges_in_cone(const Direction_2& cwBound, const Direction_2& ccwBound) {
			if (ccwBound == cwBound) {
				// Degenerate case - k = 1
				// not allowed.
				throw std::out_of_range("k should be >= 2");
			}

    		// Find angle bisector (requiring sqrt(), not exact) 
		    Line_2 cwLine(ORIGIN, cwBound);
		    Line_2 ccwLine(ORIGIN, ccwBound);
		    Direction_2 bisectorDir = bisector(cwLine, ccwLine).direction();
		   
			// Rotational transformation of cw 90 degree
			static const Transformation cw90( 0, 1, -1,  0);

			// Ordering
			Graph& g = this->g;
			// here D1 is the reverse of D1 in the book, we find this is easier to implement
			const vertex_smaller_2 orderD1 (g, ccwBound);
			const vertex_smaller_2 orderD2 (g, cwBound);
            const vertex_smaller_2 orderMid(g, cw90(bisectorDir));

			typename Graph::vertex_iterator vit, ve;
			boost::tie(vit, ve) = boost::vertices(g);

			// Step 1: Sort S according to order induced by D1
			std::vector<typename Graph::vertex_descriptor> S(vit, ve);
			std::sort(S.begin (), S.end (), orderD1);

			// Step 2: Initialise an empty set to store vertices sorted by orderD2 
			typedef CGAL::ThetaDetail::Plane_Scan_Tree<typename Graph::vertex_descriptor,
												 typename Graph::vertex_descriptor,
												 vertex_smaller_2, vertex_smaller_2> PSTree;
			PSTree pst(orderD2, orderMid);
#ifndef NDEBUG
#ifdef REALLY_VERBOSE_TREE_STATE_AFTER_EVERY_TREE_UPDATE__SAFE_TO_REMOVE_FOR_PRODUCTION
    		int i = 0;
#endif
#endif
			// Step 3: visit S in orderD1
            //     * insert pi into T
            //     * ri = T.minAbove(pi)
			for (typename std::vector<typename Graph::vertex_descriptor>::const_iterator
				it = S.begin(); it != S.end(); ++it) {
				pst.add(*it, *it);
				const typename Graph::vertex_descriptor *const ri = pst.minAbove(*it);
				if (nullptr != ri)
				  boost::add_edge(*it, *ri, g);

#ifndef NDEBUG
#ifdef REALLY_VERBOSE_TREE_STATE_AFTER_EVERY_TREE_UPDATE__SAFE_TO_REMOVE_FOR_PRODUCTION
// Prints the current tree
// To see the tree, pipe output to dot. eg
//    ./a.out <whatever arguments...> | dot -Tpng -O
// You'll get a sequence of png files:
// noname.dot.png
// noname.dot.2.png
// noname.dot.3.png
// ...etc...
//
// The tree output shades the new value added, and states what action was taken.
std::cout << "graph Plane_Scan_Tree {" << std::endl <<
  pst << std::endl << std::endl;
int j = 1;
for (auto rit = S.rbegin(); rit <= it; ++rit) {
  auto p = g[*rit];
  std::cout << "\t\"" << *rit << "\"[label=\"" << j++ << "\"";
  if (rit == it)
    std::cout << ",style=filled";
  std::cout << "];" << std::endl;
}

if (pst.size() > 1) {
  std::cout << "\t{rank=same;" << std::endl;
  std::cout << "\"" << pst.begin()->first << "\"";
  for (auto pit = ++(pst.begin()); pit != pst.end(); ++pit) {
    std::cout << "--\"" << pit->first << "\"";
  }
  std::cout << "[color=white];" << std::endl;
  std::cout << "rankdir=LR;" << std::endl;
  std::cout << "}" << std::endl;
}

std::cout << "\tlabel=\"" << ++i << ": Added (" << g[*it].x().to_double() << "," << g[*it].y().to_double() << ").";
if (nullptr != ri)
  std::cout << " -- (" << g[*ri].x().to_double() << "," << g[*ri].y().to_double() << ").";
std::cout << "\";" << std::endl;
std::cout << "\ttableloc=\"b\";" << std:: endl;
std::cout << "}" << std::endl << std::endl;
#endif
#endif
           }  // end of for
		};     // end of buildcone
		
};  // class theta_graph


/* serialization, to implement in future

	template < typename Kernel, typename Directedness, typename EdgeProperty >
	std::istream& operator>> (std::istream& is, Theta_graph_2<Kernel, Directedness, EdgeProperty>& theta_graph);

	template < typename Kernel, typename Directedness, typename EdgeProperty > 
	std::ostream& operator<< (std::ostream& os, const Theta_graph_2<Kernel, Directedness, EdgeProperty>& theta_graph);
*/

}  // namespace CGAL

#ifdef GXX11
#undef GXX11
#endif

#endif
