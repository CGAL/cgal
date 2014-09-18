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

/** @file Yao_graph_2.h
 * 
 *  This header implements the class Yao_graph_2, the constructor of which 
 *  builds a Yao graph on a set of given vertices.
 * 
 */

#ifndef CGAL_YAO_GRAPH_2_H
#define CGAL_YAO_GRAPH_2_H

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

namespace CGAL {

	/** @brief A derived class for constructing Yao graphs with a given set of 2D points.
	*  
	*  Its base class is Cone_spanners_2.
	*  Directed,undirected and bidirectional graphs are supported. For differences among these
	*  three types of graphs, please see the documentation of BGL.
	*
	*/
	template <typename Kernel, 
		      typename Directedness=boost::undirectedS, 
			  typename EdgeProperty=boost::no_property>
	class Yao_graph_2 : public Cone_spanners_2<Kernel, Directedness, EdgeProperty>
	{
	public:
        typedef typename Kernel::Direction_2             Direction_2;
	    typedef typename Kernel::Point_2                 Point_2;
		typedef typename Kernel::Aff_transformation_2    Transformation;

		typedef Cone_spanners_2<Kernel, Directedness, EdgeProperty>  base_class;
		typedef typename base_class::Graph  Graph;
		typedef typename base_class::vertex_smaller_2  vertex_smaller_2;

		// a type for the set to store vertices sorted by a direction 
    	typedef std::set<typename Graph::vertex_descriptor, vertex_smaller_2> pointSet;

		/** @brief Constructor
		*  Constructs a Yao_graph_2 Graph object.
		*
		* @param k     Number of cones to divide space into
		* @param start An iterator pointing to the first point (vertex) in the graph.
		*              (default: nullptr)
		* @param end   An iterator pointing to the place that passes the last point. (default: nullptr)
		* @param ray0  A direction denoting one of the rays deviding the
		*              cones. This allows arbitary rotations of the way space is
		*              divided.  (default: positive x-axis) 
		*/
#ifdef GXX11
	    template <typename PointInputIterator=Point_2*>
#else
	    template <typename PointInputIterator>  
#endif
		Yao_graph_2(const unsigned int k,
			        const PointInputIterator& start=nullptr, 
					const PointInputIterator& end=nullptr,
			        const Direction_2& ray0 = Direction_2(1,0)
					)
					: Cone_spanners_2<Kernel, Directedness, EdgeProperty>(k, start, end, ray0)
		{
			build_edges();
		}

		/** @brief Copy Constructor
		*  @param x  another Yao_graph_2 object to copy from.
		*/
		Yao_graph_2 (const Yao_graph_2& x) 
			: Cone_spanners_2<Kernel, Directedness, EdgeProperty>(x) {}

		/** @brief This function implements the algorithm for adding edges to build the Yao graph.
		  * The algorithm implemented is a slight adaptation to the algorithm for Theta graph described in 
		  * Giri Narasimhan and Michiel Smid, Chapter 4: Spanners based on the Theta graph, Geometric Spanner Networks,
		  * Cambridge University Press, 2007.
		  * The adaptation lies in the way how the 'closest' node is searched. 
		  * A binary tree search is not possible here, so the search here has complexity O(n), 
		  * giving rise to the complexity of O(n^2) of the entire algorithm.
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

		/**  @brief Construct edges bounded by two directions.
		*
		* @param cwBound      The direction that bounds the cone on the clockwise
		*                      direction.
		* @param ccwBound     The direction that bounds the cone on the counter-clockwise
		*                      direction.
		* @return The updated underlying graph.
        *
		*  @see  G. Narasimhan and M. Smid, Geometric Spanner Networks: Cambridge
        *        University Press, 2007,
		*/
		void add_edges_in_cone(const Direction_2& cwBound, const Direction_2& ccwBound) {
			if (ccwBound == cwBound) {
				// Degenerate case - k = 1
				// not allowed.
				throw std::out_of_range("k should be >= 2");
			}

			// Ordering
			Graph& g = this->g;
			// here D1 is the reverse of D1 in the book, we find this is easier to implement
			const vertex_smaller_2 orderD1 (g, ccwBound);
			const vertex_smaller_2 orderD2 (g, cwBound);

			typename Graph::vertex_iterator vit, ve;
			boost::tie(vit, ve) = boost::vertices(g);

			// Step 1: Sort S according to order induced by D1
			std::vector<typename Graph::vertex_descriptor> S(vit, ve);
			std::sort(S.begin (), S.end (), orderD1);

			// Step 2: Initialise an empty set to store vertices sorted by orderD2 
			pointSet pst(orderD2);

			// Step 3: visit S in orderD1
			//         * insert 'it' into T
			//         search the min in T
            for (typename std::vector<typename Graph::vertex_descriptor>::const_iterator
                it = S.begin(); it != S.end(); ++it) 
			{
					Compare_Euclidean_Distance comp(g[*it], g);

					pst.insert(*it);
					// Find the last added node - O(logn)
					typename pointSet::iterator it2 = pst.find(*it);
					// Find minimum in tree from last ended node - O(n)
					typename pointSet::iterator min = std::min_element(++it2, pst.end(), comp);
					if (min != pst.end())
						boost::add_edge(*it, *min, g);
			}

		}

		/** functor for comparing Euclidean distances of two vertices to a given vertex */
		struct Compare_Euclidean_Distance {
			const Point_2& p;
			const Graph& g;

			Compare_Euclidean_Distance(const Point_2&p, const Graph& g) : p(p), g(g) {}
			bool operator() (const typename pointSet::iterator::value_type& x, const typename pointSet::iterator::value_type& y) {
				const Point_2& px = g[x];
				const Point_2& py = g[y];
				return (px.x()-p.x())*(px.x()-p.x()) + (px.y()-p.y())*(px.y()-p.y()) < (py.x()-p.x())*(py.x()-p.x()) + (py.y()-p.y())*(py.y()-p.y());
			}
		};

};  // class yao_graph

/*  serialization, to implement in future

	template < typename Kernel, typename Directedness, typename EdgeProperty >
	std::istream& operator>> (std::istream& is, Yao_graph_2<Kernel, Directedness, EdgeProperty>& yao_graph);

	template < typename Kernel, typename Directedness, typename EdgeProperty > 
	std::ostream& operator<< (std::ostream& os, const Yao_graph_2<Kernel, Directedness, EdgeProperty>& yao_graph);

*/

}  // namespace CGAL

#ifdef GXX11
#undef GXX11
#endif

#endif
