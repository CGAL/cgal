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

//std
#include <utility>


//local
#include <CGAL/Reconstruction_triangulation_2.h>

namespace CGAL {

/*!
\ingroup PkgReconstructionSimplification2Models


\brief The class `List_output` is a model for the `ReconstructionSimplificationOutput_2` concept.

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

	typedef Reconstruction_vertex_base_2<Kernel> Reconstruction_vertex_base_2;

	typedef typename Rt_2::Vertex	 Vertex;
	typedef typename Rt_2::Edge      Edge;
	typedef typename Rt_2::Vertex_iterator Vertex_iterator;
	typedef typename Rt_2::Vertex_handle Vertex_handle;
	typedef typename Rt_2::Face_handle Face_handle;

	typedef std::list<Point> Vertices;
	typedef std::list<Segment> Edges;

	typedef typename Rt_2::Finite_edges_iterator Finite_edges_iterator;

private:
	Output_Vertex_Iterator m_v_it;
	Output_Edge_Iterator m_e_it;

	void store_solid_vertices(Rt_2& rt2) {

		for (Vertex_iterator vi = rt2.vertices_begin();
				vi != rt2.vertices_end(); ++vi)
		{
			bool incident_edges_have_sample = false;
			typename Rt_2::Edge_circulator start = rt2.incident_edges(vi);
			typename Rt_2::Edge_circulator cur   = start;

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

	void store_solid_edges(Rt_2& rt2) {
		for (Finite_edges_iterator ei = rt2.finite_edges_begin(); ei != rt2.finite_edges_end(); ++ei)
		{
			Edge edge = *ei;
			if (rt2.is_ghost(edge)) continue;

	        int index = edge.second;
	        Vertex_handle source = edge.first->vertex( (index+1)%3 );
	        Vertex_handle target = edge.first->vertex( (index+2)%3 );

			Segment s(source->point(), target->point());
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
	void store_solid_elements(Rt_2& rt2) {
		store_solid_vertices(rt2);
		store_solid_edges(rt2);
	}
};

} //namespace CGAL

#endif /* LIST_OUTPUT_H_ */
