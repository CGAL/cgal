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


#ifndef CGAL_TWO_VERTICES_PARAMETIZER_3_H_INCLUDED
#define CGAL_TWO_VERTICES_PARAMETIZER_3_H_INCLUDED

#include <CGAL/parameterization_assertions.h>

#include <cfloat>
#include <climits>

CGAL_BEGIN_NAMESPACE


//
// Declaration 
//

// Class Two_vertices_parametizer_3
// Model of BorderParametizer_3
// This class parameterizes 2 extreme vertices of a 3D surface. 
// This kind of border parameterization is used by free border parameterizations.
//
// Design pattern: 
// Two_vertices_parametizer_3 is an Strategy (see [GOF95]): it implements 
// a strategy of boundary parameterization for models of Parametizer_3
//
// Implementation note: 
// To simplify the implementation, BorderParametizer_3 models know only the 
// MeshAdaptor_3 class. They don't know the parameterization algorithm 
// requirements nor the kind of sparse linear system used.

template <class MeshAdaptor_3>			// 3D surface
class Two_vertices_parametizer_3 
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

	// Map 2 extreme vertices of the 3D mesh and mark them as "parameterized"
	// Return false on error
	bool  parameterize_border (Adaptor* mesh);

	// Indicate if border's shape is convex. 
	// Meaningless for free border parameterization algorithms.
	bool  is_border_convex () { return false; }
};


//
// Implementation 
//

// Map 2 extreme vertices of the 3D mesh and mark them as "parameterized"
// Return false on error
//
// Implementation note: 
// This method is a copy of Bruno Levy's LSCM::project() method in lscm_with_generic_api.cpp
template <class Adaptor>
inline bool
Two_vertices_parametizer_3<Adaptor>::parameterize_border(Adaptor* mesh)
{
	Vertex_iterator it;

	CGAL_parameterization_assertion(mesh != NULL);

	// Nothing to do if no boundary
	if (mesh->mesh_border_vertices_begin() == mesh->mesh_border_vertices_end())
		return false;

	std::cerr << "  map 2 vertices...";

	// Get mesh's bounding box
    double xmin =  1e30 ;
    double ymin =  1e30 ;
    double zmin =  1e30 ;
    double xmax = -1e30 ;
    double ymax = -1e30 ;
    double zmax = -1e30 ;
	for (it = mesh->mesh_vertices_begin(); it != mesh->mesh_vertices_end(); it++)
	{
		Point_3 position = mesh->get_vertex_position(it);

		xmin = std::min(position.x(), xmin) ;
		ymin = std::min(position.y(), ymin) ;		// LS 04/2005: was ", xmin)"
		zmin = std::min(position.z(), zmin) ;		// LS 04/2005: was ", xmin)"

		xmax = std::max(position.x(), xmax) ;		// LS 04/2005: was ", xmin)"
		ymax = std::max(position.y(), ymax) ;		// LS 04/2005: was ", xmin)"
		zmax = std::max(position.z(), zmax) ;		// LS 04/2005: was ", xmin)"
	}

	// Find longest bounding box axes
	double dx = xmax - xmin ;
	double dy = ymax - ymin ;
	double dz = zmax - zmin ;
	enum { X_AXIS, Y_AXIS, Z_AXIS } longest_axis, second_longest_axis;
	if(dx < dy && dx < dz) {
		if(dy > dz) {
			longest_axis        = Y_AXIS; 
			second_longest_axis = Z_AXIS;
		} else {
			longest_axis        = Z_AXIS; 
			second_longest_axis = Y_AXIS;
		}
	} else if(dy < dx && dy < dz) {
		if(dx > dz) {
			longest_axis        = X_AXIS; 
			second_longest_axis = Z_AXIS;
		} else {
			longest_axis        = Z_AXIS; 
			second_longest_axis = X_AXIS;
		}
		} else { // (dz < dx && dz < dy)
		if(dx > dy) {
			longest_axis        = X_AXIS; 
			second_longest_axis = Y_AXIS;
		} else {
			longest_axis        = Y_AXIS; 
			second_longest_axis = X_AXIS;
		}
	}
	Vector_3 V1,				// bounding box' longest axis
		     V2 ;				// bounding box' 2nd longest axis
	double V1_min, V1_max;		// bounding box' dimensions along V1
	double V2_min, V2_max;		// bounding box' dimensions along V2
	switch (longest_axis)
	{
	case X_AXIS:
		V1 = Vector_3(1,0,0) ;
		V1_min = xmin; 
		V1_max = xmax; 
		break;
	case Y_AXIS:
		V1 = Vector_3(0,1,0) ;
		V1_min = ymin; 
		V1_max = ymax; 
		break;
	case Z_AXIS:
		V1 = Vector_3(0,0,1) ;
		V1_min = zmin; 
		V1_max = zmax; 
		break;
	default:
		assert(false);
	}
	switch (second_longest_axis)
	{
	case X_AXIS:
		V2 = Vector_3(1,0,0) ;
		V2_min = xmin; 
		V2_max = xmax; 
		break;
	case Y_AXIS:
		V2 = Vector_3(0,1,0) ;
		V2_min = ymin; 
		V2_max = ymax; 
		break;
	case Z_AXIS:
		V2 = Vector_3(0,0,1) ;
		V2_min = zmin; 
		V2_max = zmax; 
		break;
	default:
		assert(false);
	}

	// Project onto longest bounding box axes,
	// Set extrema vertices' (u,v) in unit square and mark them as "parameterized"
	Vertex_handle vxmin = NULL ;
	double  umin  = DBL_MAX ;
	Vertex_handle vxmax = NULL ;
	double  umax  = DBL_MIN ;
	for (it = mesh->mesh_vertices_begin(); it != mesh->mesh_vertices_end(); it++)
	{
		Point_3  position = mesh->get_vertex_position(it);
		Vector_3 position_as_vector = position - Point_3(0,0,0);

		// coordinate along the bounding box' main axes
		double u = position_as_vector * V1 ;	
		double v = position_as_vector * V2 ;

		// LS 04/2005: convert to unit square coordinates
		assert(V1_max > V1_min);
		assert(V2_max > V2_min);
		u = (u - V1_min) / (V1_max - V1_min);		
		v = (v - V2_min) / (V2_max - V2_min);	

		mesh->set_vertex_uv(it, Point_2(u,v)) ;	// useful only for vxmin and vxmax

		if(u < umin) {
			vxmin = it ;
			umin = u ;
		} 
		if(u > umax) {
			vxmax = it ;
			umax = u ;
		} 
	}
	mesh->set_vertex_parameterized(vxmin, true) ;
	mesh->set_vertex_parameterized(vxmax, true) ;

	std::cerr << "done" << std::endl;

	return true;
}


CGAL_END_NAMESPACE

#endif //CGAL_TWO_VERTICES_PARAMETIZER_3_H_INCLUDED

