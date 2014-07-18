/*
 * Off_output.h
 *
 *  Created on: Jul 10, 2014
 *      Author: ivovigan
 */

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
of the reconstructed shape via an std::ostream object.


\tparam Kernel is the geometric kernel, used for the reconstruction and
					simplification task.
 */
template<class Kernel>
class Off_output {
public:
	typedef Reconstruction_triangulation_2<Kernel> Rt_2;
	typedef typename Rt_2::Triangulation_data_structure Tds_2;
	typedef typename Rt_2::Edge      Edge;
	typedef typename Rt_2::Vertex    Vertex;
	typedef typename Kernel::Point_2 Point;

	typedef typename Rt_2::Face_handle Face_handle;

	typedef typename CGAL::List_output<Kernel>::Output_Vertex_Iterator
				Output_Vertex_Iterator;

	typedef typename CGAL::List_output<Kernel>::Output_Edge_Iterator
				Output_Edge_Iterator;
private:
	List_output<Kernel> list_output;


	void save_one_vertex(std::ostream& os, const Vertex& v) {
	    os << v << std::endl;
	}

	void save_one_edge(std::ostream& os, const Edge& edge, std::set<Point>& edge_vertices) {
	    int i = edge.second;
	    Face_handle face = edge.first;
	    Point a = face->vertex((i+1)%3)->point();
	    Point b = face->vertex((i+2)%3)->point();

	    typename std::set<Point>::iterator it_a = edge_vertices.find(a);
	    typename std::set<Point>::iterator it_b = edge_vertices.find(b);

	    int pos_a = std::distance(edge_vertices.begin(), it_a);
	    int pos_b = std::distance(edge_vertices.begin(), it_b);

	    os << "2 "  << pos_a + list_output.vertex_count() << " "
	    		<< pos_b + list_output.vertex_count() << std::endl;
	}

	void vertices_of_edges(std::set<Point>& edge_vertices) {

		for (Output_Edge_Iterator it = list_output.edges_start();
				it != list_output.edges_beyond(); it++) {

			int i = (*it).second;
			Face_handle face = (*it).first;
			Point a = face->vertex((i+1)%3)->point();
			Point b = face->vertex((i+2)%3)->point();

			edge_vertices.insert(a);
			edge_vertices.insert(b);
	    }
	}


public:
	void store_marked_elements(Rt_2& rt2, int nb_ignore) {
		list_output.store_marked_elements(rt2, nb_ignore);
	}

	/*!
	Writes the edges and vertices of the output simplex into an `std::ostream`
	in the OFF format.

	\param os The `std::ostream` where the OFF data will be written to.
	*/
	void get_os_output(std::ostream& os) {
		std::set<Point> edge_vertices;
		vertices_of_edges(edge_vertices);

		os << "OFF " << list_output.vertex_count() + edge_vertices.size() <<
				" 0 " << list_output.edge_count()  << std::endl;

	  	for (Output_Vertex_Iterator it = list_output.vertices_start();
				it != list_output.vertices_beyond(); it++) {
	  		save_one_vertex(os, *it);
	  	}


	  	for (typename std::set<Point>::iterator it = edge_vertices.begin();
				it != edge_vertices.end(); it++) {

	  		os << *it << std::endl;
	  	}

		for (int i = 0; i < list_output.vertex_count(); i++) {
			os << "1 " <<  i << std::endl;
		}

		for (Output_Edge_Iterator it = list_output.edges_start();
				it != list_output.edges_beyond(); it++) {

			save_one_edge(os, *it,edge_vertices);
	    }
	}
};
}


#endif /* OFF_OUTPUT_H_ */
