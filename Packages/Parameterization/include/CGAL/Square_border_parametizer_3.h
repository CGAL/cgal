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


#ifndef CGAL_SQUAREBORDERPARAMETIZER_3_H
#define CGAL_SQUAREBORDERPARAMETIZER_3_H

#include <CGAL/parameterization_assertions.h>

#include <cfloat>
#include <climits>
#include <vector>

CGAL_BEGIN_NAMESPACE


//
// Declaration 
//

// Class Square_border_parametizer_3
// Model of BorderParametizer_3
// This class parameterizes the border of a 3D surface onto a square.
//
// Design pattern: 
// Square_border_parametizer_3 is an Strategy (see [GOF95]): it implements 
// a strategy of boundary parameterization for models of Parametizer_3
//
// Implementation note: 
// To simplify the implementation, BorderParametizer_3 models know only the 
// MeshAdaptor_3 class. They don't know the parameterization algorithm 
// requirements nor the kind of sparse linear system used.

template <class MeshAdaptor_3>			// 3D surface
class Square_border_parametizer_3 
{
// Public types
public:
	// Export Mesh_Adaptor_3 type and subtypes
	typedef MeshAdaptor_3					Adaptor;
	typedef typename Parametizer_3<Adaptor>::ErrorCode	
											ErrorCode;
	typedef typename Adaptor::NT			NT;
	typedef typename Adaptor::Face_handle	Face_handle;
	typedef typename Adaptor::Face_const_handle	
											Face_const_handle;
	typedef typename Adaptor::Vertex_handle	Vertex_handle;
	typedef typename Adaptor::Vertex_const_handle		
											Vertex_const_handle;
	typedef typename Adaptor::Point_3		Point_3;
	typedef typename Adaptor::Point_2		Point_2;
	typedef typename Adaptor::Vector_3		Vector_3;
	typedef typename Adaptor::Vector_2		Vector_2;
	typedef typename Adaptor::Face_iterator	Face_iterator;
	typedef typename Adaptor::Face_const_iterator		
											Face_const_iterator;
	typedef typename Adaptor::Vertex_iterator Vertex_iterator;
	typedef typename Adaptor::Vertex_const_iterator		
											Vertex_const_iterator;
	typedef typename Adaptor::Border_vertex_iterator	
											Border_vertex_iterator;
	typedef typename Adaptor::Border_vertex_const_iterator		
											Border_vertex_const_iterator;
	typedef typename Adaptor::Vertex_around_face_circulator		
											Vertex_around_face_circulator;
	typedef typename Adaptor::Vertex_around_face_const_circulator	
											Vertex_around_face_const_circulator;
	typedef typename Adaptor::Vertex_around_vertex_circulator		
											Vertex_around_vertex_circulator;
	typedef typename Adaptor::Vertex_around_vertex_const_circulator	
											Vertex_around_vertex_const_circulator;

// Public operations
public:
	// Default constructor, copy constructor and operator =() are fine

	// Assign to mesh's border vertices a 2D position (ie a (u,v) pair) 
	// on border's shape. Mark them as "parameterized".
	// Return false on error
	bool  parameterize_border (Adaptor* mesh);

	// Indicate if border's shape is convex
	bool  is_border_convex () { return true; }

// Private types
private:
	typedef typename std::vector<double>	Offset_map; 

// Private operations
private:
	// compute  total length of boundary
	double compute_boundary_length(const Adaptor& mesh);

	// Compute mesh iterator whose offset is closest to 'value'
	Border_vertex_iterator closest_iterator(Adaptor* mesh, 
											const Offset_map& offsets, 
											double value);
};


//
// Implementation 
//

// compute  total length of boundary
template <class MeshAdaptor_3>
inline 
double Square_border_parametizer_3<MeshAdaptor_3>::compute_boundary_length(
														const MeshAdaptor_3& mesh)
{
	double len = 0.0;
	for(Border_vertex_const_iterator it = mesh.mesh_border_vertices_begin(); 
		it != mesh.mesh_border_vertices_end(); 
		it++)
	{
		CGAL_parameterization_assertion(mesh.is_vertex_on_border(it));

		// Get next iterator (looping)
		Border_vertex_const_iterator next = it; 
		next++;
		if(next == mesh.mesh_border_vertices_end())
			next = mesh.mesh_border_vertices_begin();

		// Add length of it -> next vector to 'len'
		Vector_3 v = mesh.get_vertex_position(next) 
			       - mesh.get_vertex_position(it);
		len += std::sqrt(v*v);
	}
	return len;
}

// Assign to mesh's border vertices a 2D position (ie a (u,v) pair) 
// on border's shape. Mark them as "parameterized".
// Return false on error
template <class MeshAdaptor_3>
inline 
bool Square_border_parametizer_3<MeshAdaptor_3>::parameterize_border (
														MeshAdaptor_3* mesh)
{
	CGAL_parameterization_assertion(mesh != NULL);

	// Nothing to do if no boundary
	if (mesh->mesh_border_vertices_begin() == mesh->mesh_border_vertices_end())
		return false;

	// compute the total boundary length	
	double total_len = compute_boundary_length(*mesh);
	std::cerr << "  total boundary len: " << total_len << std::endl;
	CGAL_parameterization_assertion(total_len != 0);

	// map to [0,4[
	std::cerr << "  map on a square...";
 	double len = 0.0;			// current position on square in [0, total_len[
	Offset_map offsets;			// vertex index -> offset map
	offsets.reserve(mesh->count_mesh_vertices());
	Border_vertex_iterator it;
	for(it = mesh->mesh_border_vertices_begin(); 
		it != mesh->mesh_border_vertices_end(); 
		it++)
	{
		CGAL_parameterization_assertion(mesh->is_vertex_on_border(it));

		offsets[mesh->get_vertex_index(it)] = 4.0f*len/total_len;
								// current position on square in [0,4[ 

		// Get next iterator (looping)
		Border_vertex_iterator next = it; 
		next++;
		if(next == mesh->mesh_border_vertices_end())
			next = mesh->mesh_border_vertices_begin();

		// Add length of it -> next vector to 'len'
		Vector_3 v = mesh->get_vertex_position(next) 
				   - mesh->get_vertex_position(it);
		len += std::sqrt(v*v);
	}

	// First square corner is mapped to first vertex. 
	// Find closest points for three other corners.
	Border_vertex_iterator it0 = mesh->mesh_border_vertices_begin();
 	Border_vertex_iterator it1 = closest_iterator(mesh, offsets, 1.0);
 	Border_vertex_iterator it2 = closest_iterator(mesh, offsets, 2.0);
 	Border_vertex_iterator it3 = closest_iterator(mesh, offsets, 3.0);
 	assert(it1 != it0);
 	assert(it1 != it2);
 	assert(it2 != it3);
 	assert(it1 != it3);
	//
	// Snap these vertices to corners
	offsets[mesh->get_vertex_index(it0)] = 0.0;
	offsets[mesh->get_vertex_index(it1)] = 1.0;
	offsets[mesh->get_vertex_index(it2)] = 2.0;
	offsets[mesh->get_vertex_index(it3)] = 3.0;

 	// Set vertices along square's sides and mark them as "parameterized"
	for(it = it0; it != it1; it++)	// 1st side
	{
		Point_2 uv(offsets[mesh->get_vertex_index(it)], 0.0);
		mesh->set_vertex_uv(it, uv); 
//		std::cerr << "(" << uv.x() << "," << uv.y() << ") ";
		mesh->set_vertex_parameterized(it, true);
	}
	for(it = it1; it != it2; it++)								// 2nd side
	{
		Point_2 uv(1.0, offsets[mesh->get_vertex_index(it)]-1);
		mesh->set_vertex_uv(it, uv); 
//		std::cerr << "(" << uv.x() << "," << uv.y() << ") ";
		mesh->set_vertex_parameterized(it, true);
	}
	for(it = it2; it != it3; it++)								// 3rd side
	{
		Point_2 uv(3-offsets[mesh->get_vertex_index(it)], 1.0);
		mesh->set_vertex_uv(it, uv); 
//		std::cerr << "(" << uv.x() << "," << uv.y() << ") ";
		mesh->set_vertex_parameterized(it, true);
	}
	for(it = it3; it != mesh->mesh_border_vertices_end(); it++)	// 4th side
	{
		Point_2 uv(0.0, 4-offsets[mesh->get_vertex_index(it)]);
		mesh->set_vertex_uv(it, uv); 
//		std::cerr << "(" << uv.x() << "," << uv.y() << ") ";
		mesh->set_vertex_parameterized(it, true);
	}

	std::cerr << "done" << std::endl;

	return true;
}

// Utility method for parameterize_border()
// Compute mesh iterator whose offset is closest to 'value'
template <class MeshAdaptor_3>
inline 
typename Square_border_parametizer_3<MeshAdaptor_3>::Border_vertex_iterator 
Square_border_parametizer_3<MeshAdaptor_3>::closest_iterator(MeshAdaptor_3* mesh, 
															 const Offset_map& offsets, 
															 double value)
{
	Border_vertex_iterator best;
	double min = DBL_MAX;			// distance for 'best'

	for (Border_vertex_iterator it = mesh->mesh_border_vertices_begin(); 
		 it != mesh->mesh_border_vertices_end(); 
		 it++)
	{
		double d = fabs(offsets[mesh->get_vertex_index(it)] - value);
		if (d < min)
		{
			best = it;
			min = d;
		}
	}

	return best;
}


CGAL_END_NAMESPACE

#endif //CGAL_SQUAREBORDERPARAMETIZER_3_H

