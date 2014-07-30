#ifndef TDS_OUTPUT_H_
#define TDS_OUTPUT_H_

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


namespace CGAL {


/*!
\ingroup PkgReconstructionSimplification2Models


\brief The class `Tds_output` is a model for the `OutputModule` concept.

\details It provides access to the `Tds_2` of the reconstruction simplex.


\tparam Kernel is the geometric kernel, used for the reconstruction and
					simplification task.
 */
template<class Kernel>
class Tds_output {

	/// \cond SKIP_IN_MANUAL


public:
	typedef Reconstruction_triangulation_2<Kernel> Rt_2;
	typedef typename Kernel::FT                 				    FT;

	typedef typename Rt_2::Vertex	 Vertex;
	typedef typename Rt_2::Edge      Edge;
	typedef typename Rt_2::Vertex_iterator Vertex_iterator;
	typedef typename Rt_2::Finite_edges_iterator Finite_edges_iterator;


private:
	Rt_2 m_rt2;

	void mark_vertices() {
		for (Vertex_iterator vi = m_rt2.vertices_begin();
			  vi != m_rt2.vertices_end(); ++vi)
		{

		  bool incident_edges_have_sample = false;
		  typename Rt_2::Edge_circulator start = m_rt2.incident_edges(vi);
		  typename Rt_2::Edge_circulator cur = start;

		  do {
			  if (!m_rt2.is_ghost(*cur)) {
				  incident_edges_have_sample = true;
				  break;
			  }
			  ++cur;
		  } while (cur != start);

		  if (!incident_edges_have_sample) {
			  if ((*vi).has_sample_assigned())
				  (*vi).set_relevance(1);
		  }
		}
	}

	void mark_edges() {

		for (Finite_edges_iterator ei = m_rt2.finite_edges_begin(); ei != m_rt2.finite_edges_end(); ++ei)
		{
			Edge edge = *ei;
			FT relevance = 0;
			if (!m_rt2.is_ghost(edge)) {
				relevance = m_rt2.get_edge_relevance(edge); // >= 0
			}
			edge.first->relevance(edge.second) = relevance;
		}
	}

	/// \endcond
public:

	/*!
	Extracts the solid edges and vertices from the `Reconstruction_simplification_2` module.

	\param rt2 The `Reconstruction_triangulation_2` from which the solid edges and vertices are extracted.
	\param nb_ignore The number of verticess to be ignored in the output.

	*/
	void store_marked_elements(Rt_2& rt2, int nb_ignore) {
		  m_rt2 = rt2;
		  mark_vertices();
		  mark_edges();
	}


	/*!
	Allows accessing the `Reconstruction_triangulation_2` of the `Reconstruction_simplification_2` module.

	\param rt2 The `Reconstruction_triangulation_2`.
	*/
	void extract_reconstruction_tds(Rt_2& rt2) {
		rt2 = m_rt2;
	}
};
}


#endif /* TDS_OUTPUT_H_ */
