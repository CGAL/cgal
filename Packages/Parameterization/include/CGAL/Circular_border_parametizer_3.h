// Copyright (c) 2005  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$
// $Name$
//
// Author(s)     : Laurent Saboret, Bruno Levy, Pierre Alliez


#ifndef CIRCULARBORDERPARAMETIZER_3_H
#define CIRCULARBORDERPARAMETIZER_3_H

#include "CGAL/parameterization_assertions.h"

CGAL_BEGIN_NAMESPACE


//
// Declaration 
//

// Class Circular_border_parametizer_3
// Model of BorderParametizer_3
// This class parameterizes the border of a 3D surface onto a circle.
// 
// Implementation note: to simplify the implementation, BorderParametizer_3 knows only the MeshAdaptor_3 class. It doesn't 
// know the parameterization algorithm requirements nor the kind of sparse linear system used.
template <class MeshAdaptor_3>			// 3D surface
class Circular_border_parametizer_3 
{
// Public types
public:
				typedef MeshAdaptor_3									Mesh_adaptor_3;

// Public operations
public:
				// Default constructor, copy constructor and operator =() are fine

				// Assign to mesh's border vertices a 2D position (ie a (u,v) pair) on border's shape
				// Mark them as "parameterized"
				// Return false on error
				bool  parameterize_border (MeshAdaptor_3* mesh);

				// Indicate if border's shape is convex
				bool  is_border_convex () { return true; }

// Private types
private:
				typedef typename MeshAdaptor_3::Border_vertex_iterator	Border_vertex_iterator;
				typedef typename MeshAdaptor_3::Point_2					Point_2;
				typedef typename MeshAdaptor_3::Vector_3				Vector_3;

// Private operations
private:
				// compute  total length of boundary
				double compute_boundary_length(MeshAdaptor_3* mesh);
};


//
// Implementation 
//

// compute  total length of boundary
template <class MeshAdaptor_3>
inline 
double Circular_border_parametizer_3<MeshAdaptor_3>::compute_boundary_length(MeshAdaptor_3* mesh)
{
//std::cerr << "  compute boundary length: boundary vertices are...";
	double len = 0.0;
	for(Border_vertex_iterator it = mesh->mesh_border_vertices_begin(); it != mesh->mesh_border_vertices_end(); it++)
	{
		CGAL_parameterization_assertion(mesh->is_vertex_on_border(*it));
//std::cerr << " (" << mesh->get_vertex_position(*it).x() << "," << mesh->get_vertex_position(*it).y() << "," << mesh->get_vertex_position(*it).z() << ")";

		// Get next iterator (looping)
		Border_vertex_iterator next = it; 
		next++;
		if(next == mesh->mesh_border_vertices_end())
			next = mesh->mesh_border_vertices_begin();

		// Add length of it -> next vector to 'len'
		Vector_3 v = mesh->get_vertex_position(*next) - mesh->get_vertex_position(*it);
		len += std::sqrt(v*v);
	}
//std::cerr << " done" << std::endl;
	return len;
}

// Assign to mesh's border vertices a 2D position (ie a (u,v) pair) on border's shape
// Mark them as "parameterized"
// Return false on error
template <class MeshAdaptor_3>
inline 
bool Circular_border_parametizer_3<MeshAdaptor_3>::parameterize_border (MeshAdaptor_3* mesh)
{
	CGAL_parameterization_assertion(mesh != NULL);

	// Nothing to do if no boundary
	if (mesh->mesh_border_vertices_begin() == mesh->mesh_border_vertices_end())
		return false;

	// compute the total boundary length	
	double total_len = compute_boundary_length(mesh);
	std::cerr << "  total boundary len: " << total_len << std::endl;
	CGAL_parameterization_assertion(total_len != 0);

	std::cerr << "  map on a circle...";
	const double PI = 3.14159265359;
	const double tmp = 2*PI/total_len;
	double len = 0.0;						// current position on the circle in [0, total_len]
//	Border_vertex_iterator first = mesh->mesh_border_vertices_begin();
//	bool border = (*first)->is_border();
	for(Border_vertex_iterator it = mesh->mesh_border_vertices_begin(); it != mesh->mesh_border_vertices_end(); it++)
	{
		CGAL_parameterization_assertion(mesh->is_vertex_on_border(*it));

		double angle = len*tmp;				// current position on the circle in radians

		// map vertex on unit circle
		Point_2 uv;
		//if(border)
			uv = Point_2(0.5+0.5*cos(-angle),0.5+0.5*sin(-angle));
		//else
		//	uv = Point_2(0.5+0.5*cos(angle),0.5+0.5*sin(angle));
//		std::cerr << "(" << uv.x() << "," << uv.y() << ") ";
		mesh->set_vertex_uv(&*it, uv);

		// Mark vertex as "parameterized"
		mesh->set_vertex_parameterized(&*it, true);

		// Get next iterator (looping)
		Border_vertex_iterator next = it; 
		next++;
		if(next == mesh->mesh_border_vertices_end())
			next = mesh->mesh_border_vertices_begin();

		// Add length of it -> next vector to 'len'
		Vector_3 v = mesh->get_vertex_position(*next) - mesh->get_vertex_position(*it);
		len += std::sqrt(v*v);
	}

  std::cerr << "done" << std::endl;

  return true;
}


CGAL_END_NAMESPACE

#endif //CIRCULARBORDERPARAMETIZER_3_H

