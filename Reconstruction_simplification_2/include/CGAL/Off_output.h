#ifndef OFF_OUTPUT_H_
#define OFF_OUTPUT_H_


// Copyright (c) 2014  INRIA Sophia-Antipolis (France), INRIA Lorraine LORIA.
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
// Author(s)     : Ivo Vigan

//CGAL
#include <CGAL/basic.h>

//local
#include <CGAL/Reconstruction_triangulation_2.h>
#include <CGAL/List_output.h>

#include <ostream>

namespace CGAL {

/*!
\ingroup PkgReconstructionSimplification2Models


\brief The class `Off_output` is a model for the `OutputModule` concept.

\details It allows accessing the isolated vertices and the edges
of the reconstructed shape in the OFF format via an std::ostream object.


\tparam Kernel is the geometric kernel, used for the reconstruction and
					simplification task.
 */
template<class Kernel>
class Off_output {


	/// \cond SKIP_IN_MANUAL

public:
	typedef Reconstruction_triangulation_2<Kernel> Rt_2;
	typedef typename Rt_2::Triangulation_data_structure Tds_2;
	typedef typename Kernel::Point_2 Point;
	typedef typename Kernel::Segment_2 Segment;
	typedef typename Rt_2::Face_handle Face_handle;

private:
	typedef std::back_insert_iterator<std::vector<Point> >   Point_it;
	typedef std::back_insert_iterator<std::vector<Segment> > Edge_it;

	std::vector<Point> isolated_points;
	std::vector<Segment> edges;

	CGAL::List_output<Kernel, Point_it, Edge_it> list_output;


	void save_one_edge(std::ostream& os, const Segment& edge, std::set<Point>& edge_vertices) {
	    Point a = edge.source();
	    Point b = edge.target();

	    typename std::set<Point>::iterator it_a = edge_vertices.find(a);
	    typename std::set<Point>::iterator it_b = edge_vertices.find(b);

	    int pos_a = std::distance(edge_vertices.begin(), it_a);
	    int pos_b = std::distance(edge_vertices.begin(), it_b);

	    os << "2 "  << pos_a + isolated_points.size() << " "
	    		<< pos_b + isolated_points.size() << std::endl;
	}

	void vertices_of_edges(std::set<Point>& edge_vertices) {

	  	for (typename std::vector<Segment>::iterator it = edges.begin();
				it != edges.end(); it++) {

			Point a = (*it).source();
		    Point b = (*it).target();

			edge_vertices.insert(a);
			edge_vertices.insert(b);
	    }
	}

	/// \endcond

public:

	Off_output() : list_output(Point_it(isolated_points), Edge_it(edges)) { }

	/*!
	Extracts the solid edges and vertices from the `Reconstruction_simplification_2` module.

	\param rt2 The `Reconstruction_triangulation_2` from which the solid edges and vertices are extracted.

	*/
	void store_marked_elements(Rt_2& rt2) {
		list_output.store_marked_elements(rt2);
	}

	/*!
	Writes the edges and vertices of the output simplex into an `std::ostream`
	in the OFF format.

	\param os The `std::ostream` where the OFF data will be written to.
	*/
	void get_os_output(std::ostream& os) {
		std::set<Point> edge_vertices;
		vertices_of_edges(edge_vertices);

		os << "OFF " << isolated_points.size() + edge_vertices.size() <<
				" 0 " << edges.size()  << std::endl;

		for (typename std::vector<Point>::iterator it = isolated_points.begin();
					it != isolated_points.end(); it++) {
			os << *it << std::endl;
	  	}

	  	for (typename std::set<Point>::iterator it = edge_vertices.begin();
				it != edge_vertices.end(); it++) {

	  		os << *it << std::endl;
	  	}

		for (int i = 0; i < isolated_points.size(); i++) {
			os << "1 " <<  i << std::endl;
		}

	  	for (typename std::vector<Segment>::iterator it = edges.begin();
				it != edges.end(); it++) {

	  		save_one_edge(os, *it,edge_vertices);
	    }
	}
};
}


#endif /* OFF_OUTPUT_H_ */
