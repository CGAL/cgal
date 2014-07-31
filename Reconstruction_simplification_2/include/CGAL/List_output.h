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

\details It takes two `Output-Iterators`, one for storing the
isolated points and one for storing the edges of the reconstructed shape.


\tparam Kernel is the geometric kernel, used for the reconstruction and
					simplification task.

\tparam Output_Vertex_Iterator The `Output-Iterator` type for storing the points

\tparam Output_Edge_Iterator The `Output-Iterator` type for storing the
										edges (as Segments).
 */
template<class Kernel, class Output_Vertex_Iterator, class Output_Edge_Iterator>
class List_output {
public:

	/// \cond SKIP_IN_MANUAL

	typedef typename Kernel::FT                 				    FT;
	typedef typename Kernel::Point_2                 				Point;
	typedef typename Kernel::Segment_2                 				Segment;

	typedef Reconstruction_triangulation_2<Kernel> Rt_2;

	typedef typename Rt_2::Vertex	 Vertex;
	typedef typename Rt_2::Edge      Edge;
	typedef typename Rt_2::Vertex_iterator Vertex_iterator;

	typedef std::list<Point> Vertices;
	typedef std::list<Segment> Edges;

	typedef typename Rt_2::Reconstruction_edge_2 Reconstruction_edge_2;
	typedef typename Rt_2::MultiIndex MultiIndex;

	typedef typename Rt_2::Finite_edges_iterator Finite_edges_iterator;

private:
	Output_Vertex_Iterator m_v_it;
	Output_Edge_Iterator m_e_it;

	void store_marked_vertices(Rt_2& rt2) {

		for (Vertex_iterator vi = rt2.vertices_begin();
				vi != rt2.vertices_end(); ++vi)
		{
			bool incident_edges_have_sample = false;
			typename Rt_2::Edge_circulator start = rt2.incident_edges(vi);

			typename Rt_2::Edge_circulator cur = start;

			do {
				if (!rt2.is_ghost(*cur)) {
					incident_edges_have_sample = true;
					break;
				}
				++cur;
			} while (cur != start);

			if (!incident_edges_have_sample) {
				if ((*vi).has_sample_assigned()) {
					Point p = (*vi).point();
					*m_v_it = p;
					m_v_it++;
				}
			}
		}
	}

	void store_marked_edges(Rt_2& rt2, int nb_ignore) {
		MultiIndex mindex;
		for (Finite_edges_iterator ei = rt2.finite_edges_begin(); ei != rt2.finite_edges_end(); ++ei)
		{
			Edge edge = *ei;
			if (rt2.is_ghost(edge)) continue;
			FT value = rt2.get_edge_relevance(edge); // >= 0
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
			Segment s(pedge.source()->point(), pedge.target()->point());
			//edges.push_back(s);
			*m_e_it = s;
			m_e_it++;

		}
	}

	/// \endcond
public:

	/// \name Creation
	  /// @{

	/*!
		 Instantiates a new List_output object
				  for two given Output_Iterators.

		 \param v_it An Output_Vertex_Iterator for storing the points.

		 \param e_it An Output_Edge_Iterator for storing the edges (as Segments).
	 */
	List_output(Output_Vertex_Iterator v_it, Output_Edge_Iterator e_it)  :
		m_v_it(v_it), m_e_it(e_it) { }
	  /// @}



	/*!
	Extracts the solid edges and vertices from the `Reconstruction_simplification_2` module.

	\param rt2 The `Reconstruction_triangulation_2` from which the solid edges and vertices are extracted.
	*/
	void store_marked_elements(Rt_2& rt2) {
		store_marked_vertices(rt2);
		store_marked_edges(rt2, 0); //TODO: IV do we want the nb_ignore parameter
	}
};


} //namespace CGAL

#endif /* LIST_OUTPUT_H_ */
