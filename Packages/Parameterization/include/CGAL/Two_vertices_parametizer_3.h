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


#ifndef Two_vertices_parametizer_3_H_INCLUDED
#define Two_vertices_parametizer_3_H_INCLUDED

#include "CGAL/parameterization_assertions.h"

#include <cfloat>
#include <climits>

CGAL_BEGIN_NAMESPACE


//
// Declaration 
//

// Class Two_vertices_parametizer_3
// Model of BorderParametizer_3
// This class parameterizes 2 extreme vertices of a 3D surface. This kind of border parameterization is used by free border parameterization algorithms.
//
// Design pattern: Two_vertices_parametizer_3 is an Strategy (see [GOF95]): it implements a strategy of boundary parameterization for models of Parametizer_3
//
// Implementation note: to simplify the implementation, BorderParametizer_3 models know only the MeshAdaptor_3 class. They don't 
//                      know the parameterization algorithm requirements nor the kind of sparse linear system used.
template <class MeshAdaptor_3>			// 3D surface
class Two_vertices_parametizer_3 
{
// Public types
public:
				// Export Mesh_Adaptor_3 type and subtypes
				typedef MeshAdaptor_3													Mesh_adaptor_3;
				typedef typename Parametizer_3<MeshAdaptor_3>::ErrorCode				ErrorCode;
				typedef typename MeshAdaptor_3::NT										NT;
				typedef typename MeshAdaptor_3::Face_handle								Face_handle;
				typedef typename MeshAdaptor_3::Face_const_handle						Face_const_handle;
				typedef typename MeshAdaptor_3::Vertex_handle							Vertex_handle;
				typedef typename MeshAdaptor_3::Vertex_const_handle						Vertex_const_handle;
				typedef typename MeshAdaptor_3::Point_3									Point_3;
				typedef typename MeshAdaptor_3::Point_2									Point_2;
				typedef typename MeshAdaptor_3::Vector_3								Vector_3;
				typedef typename MeshAdaptor_3::Face_iterator							Face_iterator;
				typedef typename MeshAdaptor_3::Face_const_iterator						Face_const_iterator;
				typedef typename MeshAdaptor_3::Vertex_iterator							Vertex_iterator;
				typedef typename MeshAdaptor_3::Vertex_const_iterator					Vertex_const_iterator;
				typedef typename MeshAdaptor_3::Border_vertex_iterator					Border_vertex_iterator;
				typedef typename MeshAdaptor_3::Border_vertex_const_iterator			Border_vertex_const_iterator;
				typedef typename MeshAdaptor_3::Vertex_around_face_circulator			Vertex_around_face_circulator;
				typedef typename MeshAdaptor_3::Vertex_around_face_const_circulator		Vertex_around_face_const_circulator;
				typedef typename MeshAdaptor_3::Vertex_around_vertex_circulator			Vertex_around_vertex_circulator;
				typedef typename MeshAdaptor_3::Vertex_around_vertex_const_circulator	Vertex_around_vertex_const_circulator;

// Public operations
public:
				// Default constructor, copy constructor and operator =() are fine

				// Map 2 extreme vertices of the 3D mesh and mark them as "parameterized"
				// Return false on error
				bool  parameterize_border (MeshAdaptor_3* mesh);

				// Indicate if border's shape is convex. Meaningless for free border parameterization algorithms.
				bool  is_border_convex () { return false; }
};


//
// Implementation 
//

// Map 2 extreme vertices of the 3D mesh and mark them as "parameterized"
// Return false on error
//
// Note: this method is a copy of Bruno Levy's LSCM::project() method in lscm_with_generic_api.cpp
template <class MeshAdaptor_3>
inline 
bool Two_vertices_parametizer_3<MeshAdaptor_3>::parameterize_border (MeshAdaptor_3* mesh)
{
	Vertex_iterator it;

	CGAL_parameterization_assertion(mesh != NULL);

	// Nothing to do if no boundary
	if (mesh->mesh_border_vertices_begin() == mesh->mesh_border_vertices_end())
		return false;

	std::cerr << "  map 2 vertices...";

	// Get surface's bounding box
	double xmin = DBL_MAX ;
	double ymin = DBL_MAX ;
	double zmin = DBL_MAX ;
	double xmax = DBL_MIN ;
	double ymax = DBL_MIN ;
	double zmax = DBL_MIN ;
	for (it = mesh->mesh_vertices_begin(); it != mesh->mesh_vertices_end(); it++)
	{
		Point_3 position = mesh->get_vertex_position(it);

		xmin = std::min(position.x(), xmin) ;
		ymin = std::min(position.y(), xmin) ;
		zmin = std::min(position.z(), xmin) ;

		xmax = std::max(position.x(), xmin) ;
		ymax = std::max(position.y(), xmin) ;
		zmax = std::max(position.z(), xmin) ;
	}

	// Find shortest bounding box axis
	Vector_3 V1,V2 ;
	double dx = xmax - xmin ;
	double dy = ymax - ymin ;
	double dz = zmax - zmin ;
	if(dx < dy && dx < dz) {
		if(dy > dz) {
			V1 = Vector_3(0,1,0) ;
			V2 = Vector_3(0,0,1) ;
		} else {
			V2 = Vector_3(0,1,0) ;
			V1 = Vector_3(0,0,1) ;
		}
	} else if(dy < dx && dy < dz) {
		if(dx > dz) {
			V1 = Vector_3(1,0,0) ;
			V2 = Vector_3(0,0,1) ;
		} else {
			V2 = Vector_3(1,0,0) ;
			V1 = Vector_3(0,0,1) ;
		}
	} else if(dz < dx && dz < dy) {
		if(dx > dy) {
			V1 = Vector_3(1,0,0) ;
			V2 = Vector_3(0,1,0) ;
		} else {
			V2 = Vector_3(1,0,0) ;
			V1 = Vector_3(0,1,0) ;
		}
	}

	// Project onto shortest bounding box axis,
	// and mark extrema vertices as "parameterized"
	Vertex_handle vxmin = NULL ;
	double  umin  = DBL_MAX ;
	Vertex_handle vxmax = NULL ;
	double  umax  = DBL_MIN ;
	for (it = mesh->mesh_vertices_begin(); it != mesh->mesh_vertices_end(); it++)
	{
		Point_3  position = mesh->get_vertex_position(it);
		Vector_3 position_as_vector = position - Point_3(0,0,0);

		double u = position_as_vector * V1 ;	/* dot product */
		double v = position_as_vector * V2 ;	/* dot product */
		mesh->set_vertex_uv(it, Point_2(u,v)) ;	// LS 02/05: I guess that this is useful only for *vxmin and *vxmax

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

#endif //Two_vertices_parametizer_3_H_INCLUDED

