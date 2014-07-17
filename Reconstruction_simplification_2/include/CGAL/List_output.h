#ifndef LIST_OUTPUT_H_
#define LIST_OUTPUT_H_

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
#include <CGAL/Cost.h>


// boost
#include <boost/multi_index/mem_fun.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/member.hpp>

namespace CGAL {

/*!
\ingroup PkgReconstructionSimplification2Models


\brief The class `List_output` is a model for the `OutputModule` concept.

\details It returns Output-iterators which allow iterating over both the
isolated vertices and the edges of the reconstructed shape


\tparam Kernel is the geometric kernel, used for the reconstruction and
					simplification task.
 */
template<class Kernel>
class List_output {
public:
	typedef typename Kernel::FT                 				    FT;
	typedef Cost<FT>												Cost;
	typedef Reconstruction_triangulation_2<Kernel> Rt_2;
	typedef typename Kernel::Segment_2 Segment;
	typedef typename Kernel::Point_2 Point;

	typedef typename Rt_2::Vertex_handle Vertex_handle;

	typedef typename Rt_2::Edge_iterator   Edge_iterator;
	typedef typename Rt_2::Vertex_iterator Vertex_iterator;

	typedef typename Rt_2::Vertex	 Vertex;
	typedef typename Rt_2::Edge      Edge;

	typedef std::list<Vertex> Vertices;
	typedef std::list<Edge> Edges;
	typedef typename Vertices::const_iterator 		Output_Vertex_Iterator;
	typedef typename Edges::const_iterator 			Output_Edge_Iterator;

	typedef typename Rt_2::Reconstruction_edge_2 Reconstruction_edge_2;
	typedef typename Rt_2::MultiIndex MultiIndex;

	typedef typename Rt_2::Finite_edges_iterator Finite_edges_iterator;



private:

	Vertices vertices;
	Edges    edges;

	FT get_mass(const Edge& edge) const {
		return edge.first->mass(edge.second);
	}

	const Cost& get_cost(const Edge& edge) const {
		return edge.first->cost(edge.second);
	}

	Vertex_handle source_vertex(const Edge& edge) const {
			return edge.first->vertex(Rt_2::ccw(edge.second));
	}

	Vertex_handle target_vertex(const Edge& edge) const {
		return edge.first->vertex(Rt_2::cw(edge.second));
	}

	Segment get_segment(const Edge& edge) const {
		const Point& ps = source_vertex(edge)->point();
		const Point& pt = target_vertex(edge)->point();
		return Segment(ps, pt);
	}

	FT get_length(const Edge& edge) const {
		Segment segment = get_segment(edge);
		return std::sqrt(segment.squared_length());
	}

	FT get_edge_relevance(const Edge& edge) const {
		FT M = get_mass(edge);
		if (M == 0.0)
			return 0.0;

		FT L = get_length(edge);
		FT cost = get_cost(edge).finalize();
		return M * L * L / cost;
	}

	bool is_ghost(const Edge& edge) const {
		return edge.first->ghost(edge.second);
	}

public:
	  inline Output_Vertex_Iterator vertices_start() const {
	    return vertices.begin();
	  }

	  inline Output_Vertex_Iterator vertices_beyond() const {
	    return vertices.end();
	  }

	  inline Output_Edge_Iterator edges_start() const {
		 return edges.begin();
	  }

	  inline Output_Edge_Iterator edges_beyond() const {
	    return edges.end();
	  }


	  void clear() {
		  vertices.clear();
		  edges.clear();
	  }

	  int vertex_count() {
		  return vertices.size();
	  }

	  int edge_count() {
		  return edges.size();
	  }

	  void store_marked_vertices(Rt_2& rt2) {

		  for (Vertex_iterator vi = rt2.vertices_begin();
				  vi != rt2.vertices_end(); ++vi)
		  {
			  bool incident_edges_have_sample = false;
			  typename Rt_2::Edge_circulator start = rt2.incident_edges(vi);

			  typename Rt_2::Edge_circulator cur = start;

			  do {
				  if (!is_ghost(*cur)) {

					  incident_edges_have_sample = true;
					  break;
				  }
				  ++cur;
			  } while (cur != start);

			  if (!incident_edges_have_sample) {
				  if ((*vi).has_sample_assigned())
					  vertices.push_back(*vi);
			  }
		  }
	  }

	  void store_marked_edges(Rt_2& rt2, int nb_ignore) {
		 MultiIndex mindex;
		 for (Finite_edges_iterator ei = rt2.finite_edges_begin(); ei != rt2.finite_edges_end(); ++ei)
		 {
			Edge edge = *ei;
			if (is_ghost(edge)) continue;
			FT value = get_edge_relevance(edge); // >= 0
			mindex.insert(Reconstruction_edge_2(edge, value));
		 }


		int nb_remove = (std::min)(nb_ignore, int(mindex.size()));

		for (int i = 0; i < nb_remove; ++i)
		{
			Reconstruction_edge_2 pedge = *(mindex.template get<1>()).begin();
			(mindex.template get<0>()).erase(pedge);
		}

		while (!mindex.empty())
		{
			Reconstruction_edge_2 pedge = *(mindex.template get<1>()).begin();
			(mindex.template get<0>()).erase(pedge);
			edges.push_back(pedge.edge());

		}
	  }

	  void store_marked_elements(Rt_2& rt2, int nb_ignore) {
		  store_marked_vertices(rt2);
		  store_marked_edges(rt2, nb_ignore);
	}
};


} //namespace CGAL

#endif /* LIST_OUTPUT_H_ */
