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


#ifndef MESHADAPTORPOLYHEDRONEX_H
#define MESHADAPTORPOLYHEDRONEX_H

#include <CGAL/iterator.h>
#include <CGAL/circulator.h>
#include <CGAL/Iterator_project.h>
#include <CGAL/Circulator_project.h>
#include <CGAL/function_objects.h>

#include "CGAL/parameterization_assertions.h"

#include "cgal_types.h"		
#include "Mesh_feature_extractor.h"

#include <cassert>


namespace CGAL 
{
	// Utility class for Mesh_adaptor_polyhedron_ex
	// This class is used to generate the Border_vertex_iterator type
	template <>	struct Project_vertex<Polyhedron_ex::Halfedge_handle> {
		typedef Polyhedron_ex::Halfedge_handle	argument_type;
		typedef Polyhedron_ex::Vertex			Vertex;
		typedef Vertex							result_type;
		typedef Arity_tag<1>					Arity;
		// Get the vertex inside a Halfedge handle
		Vertex&       operator()( Polyhedron_ex::Halfedge_handle& x)       const { return *(x->vertex()); }
		const Vertex& operator()( const Polyhedron_ex::Halfedge_handle& x) const { return *(x->vertex()); }
	};

	// Utility class for Mesh_adaptor_polyhedron_ex
	// This class is used to generate the Vertex_around_face_circulator type
	template <>	struct Project_vertex<Polyhedron_ex::Halfedge> {
		typedef Polyhedron_ex::Halfedge			argument_type;
		typedef Polyhedron_ex::Vertex			Vertex;
		typedef Vertex							result_type;
		typedef Arity_tag<1>					Arity;
		// Get the vertex of a Halfedge
		Vertex&       operator()( Polyhedron_ex::Halfedge& x)       const { return *(x.vertex()); }
		const Vertex& operator()( const Polyhedron_ex::Halfedge& x) const { return *(x.vertex()); }
	};

	// Utility class for Mesh_adaptor_polyhedron_ex
	// This class is used to generate the Vertex_around_vertex_circulator type
	template <class Node>	
	struct Project_opposite_vertex {
		typedef Node							argument_type;
		typedef typename Node::Vertex			Vertex;
		typedef Vertex							result_type;
		typedef Arity_tag<1>					Arity;
		// Get the opposite vertex 
		Vertex&       operator()( Node& x)       const { return *(x.opposite()->vertex()); }
		const Vertex& operator()( const Node& x) const { return *(x.opposite()->vertex()); }
	};
}; // mamespace CGAL


// Utility class for Mesh_adaptor_polyhedron_ex
// This class adds a conversion to handle/const handle to an iterator class
template <class Iter,
          class ConstHandle,
          class Handle = void*>
class Convertible_iterator : public Iter {
	typedef Iter											Base; 
	typedef Convertible_iterator<Iter, ConstHandle, Handle>	Self;
public:
	Convertible_iterator() {}
	Convertible_iterator(Base b) : Base(b) {}
	Self & operator++()		{ Base::operator++(); return *this; }
	Self & operator--()		{ Base::operator--(); return *this; }
	Self operator++(int)	{ Self tmp(*this); ++(*this); return tmp; }
	Self operator--(int)	{ Self tmp(*this); --(*this); return tmp; }
	operator Handle()				{ return Base::ptr(); }
	operator ConstHandle() const	{ return Base::ptr(); }
};


// Class Mesh_adaptor_polyhedron_ex 
// Model of MeshAdaptor_3 concept
//
// Adaptor class to access to a Polyhedron_ex 3D mesh on an uniform manner. 
// The input mesh can be of any genus, but it has to come with a description of the boundary of a topological disc. If no boundary is given, 
// we assume that it coincides with the unique boundary already in the input mesh. 
class Mesh_adaptor_polyhedron_ex
{
// Private types
private:
				typedef Polyhedron_ex::Facet								Face;
				typedef Polyhedron_ex::Vertex								Vertex;
				typedef Polyhedron_ex::Halfedge								Halfedge;
				typedef Polyhedron_ex::Halfedge_handle						Halfedge_handle;
				typedef Feature_backbone<Polyhedron_ex::Vertex_handle, Halfedge_handle>	
																			Backbone;
				// Iterator over mesh boundary vertices
				typedef CGAL::Iterator_project<std::list<Halfedge_handle>::iterator, 
											   CGAL::Project_vertex<Halfedge_handle> >	
																			Border_vertex_iterator_base;
				typedef CGAL::Iterator_project<std::list<Halfedge_handle>::const_iterator, 
											   CGAL::Project_vertex<Halfedge_handle> >	
																			Border_vertex_const_iterator_base;
				// Circulator over a face's vertices
				typedef CGAL::Circulator_project<Polyhedron_ex::Halfedge_around_facet_circulator, 
												 CGAL::Project_vertex<Halfedge>, Vertex&, Vertex*>					
																			Vertex_around_face_circulator_base;
				typedef CGAL::Circulator_project<Polyhedron_ex::Halfedge_around_facet_const_circulator, 
											     CGAL::Project_vertex<Halfedge>, const Vertex&, const Vertex*>					
																			Vertex_around_face_const_circulator_base;
				// Circulator over the vertices incident to a vertex
				typedef CGAL::Circulator_project<Polyhedron_ex::Halfedge_around_vertex_circulator, 
												 CGAL::Project_opposite_vertex<Halfedge>, Vertex&, Vertex*>					
																			Vertex_around_vertex_circulator_base;
				typedef CGAL::Circulator_project<Polyhedron_ex::Halfedge_around_vertex_const_circulator, 
												 CGAL::Project_opposite_vertex<Halfedge>, const Vertex&, const Vertex*>					
																			Vertex_around_vertex_const_circulator_base;
	
// Public types
public:
				typedef Polyhedron_ex::Traits::FT							NT;
				typedef Polyhedron_ex::Facet_handle							Face_handle;
				typedef Polyhedron_ex::Facet_const_handle					Face_const_handle;
				typedef Polyhedron_ex::Vertex_handle						Vertex_handle;
				typedef Polyhedron_ex::Vertex_const_handle					Vertex_const_handle;
				typedef Polyhedron_ex::Traits::Point_2						Point_2;
				typedef Polyhedron_ex::Traits::Point_3						Point_3;
				typedef Polyhedron_ex::Traits::Vector_2						Vector_2;
				typedef Polyhedron_ex::Traits::Vector_3						Vector_3;
				// Iterator over all mesh faces
				typedef Polyhedron_ex::Facet_iterator						Face_iterator;
				typedef Polyhedron_ex::Facet_const_iterator					Face_const_iterator;
				// Iterator over all mesh vertices
				typedef Polyhedron_ex::Vertex_iterator						Vertex_iterator;
				typedef Polyhedron_ex::Vertex_const_iterator				Vertex_const_iterator;
				// Iterator over mesh boundary vertices
				typedef Convertible_iterator<Border_vertex_iterator_base, Vertex_const_handle, Vertex_handle>
																			Border_vertex_iterator;		
				typedef Convertible_iterator<Border_vertex_const_iterator_base, Vertex_const_handle>
																			Border_vertex_const_iterator;	
				// Circulator over a face's vertices
				typedef Convertible_iterator<Vertex_around_face_circulator_base, Vertex_const_handle, Vertex_handle>
																			Vertex_around_face_circulator;	
				typedef Convertible_iterator<Vertex_around_face_const_circulator_base, Vertex_const_handle>
																			Vertex_around_face_const_circulator; 
				// Circulator over the vertices incident to a vertex
				typedef Convertible_iterator<Vertex_around_vertex_circulator_base, Vertex_const_handle, Vertex_handle>
																			Vertex_around_vertex_circulator;	
				typedef Convertible_iterator<Vertex_around_vertex_const_circulator_base, Vertex_const_handle>
																			Vertex_around_vertex_const_circulator; 

// Public operations
public:
				//
				// LIFE CYCLE
				//

				// Create an adaptator for an existing Polyhedron_ex mesh
				// The input mesh must be a topological disc (thus have a unique boundary)
				Mesh_adaptor_polyhedron_ex(Polyhedron_ex* mesh) 
				{
					assert(mesh != NULL);
					m_mesh = mesh;

					// Extract mesh boundary
					m_boundary = extract_unique_boundary();
				}
				// Create an adaptator for an existing Polyhedron_ex mesh
				// The input mesh can be of any genus, but it has to come with a description of the boundary of a topological disc 
				template <class InputIterator>
				Mesh_adaptor_polyhedron_ex(Polyhedron_ex* mesh,
										   InputIterator first_boundary_halfedge, InputIterator last_boundary_halfedge) 
				{
					assert(mesh != NULL);
					m_mesh = mesh;

					// Copy input boundary
					std::copy(first_boundary_halfedge, last_boundary_halfedge, std::back_inserter(m_boundary));
				}

				// Default copy constructor and operator =() are fine

				// Update the hidden Polyhedron_ex mesh when parameterization is over
				~Mesh_adaptor_polyhedron_ex() 
				{
					// Update the halfedges' (u,v) from the vertices
					fprintf(stderr,"  update texture coordinates per corner...");
					Polyhedron_ex::Halfedge_iterator pHalfedge;
					for(pHalfedge = m_mesh->halfedges_begin(); pHalfedge != m_mesh->halfedges_end(); pHalfedge++)
					{
  						Polyhedron_ex::Vertex_handle pVertex =  pHalfedge->vertex();
						if (is_vertex_parameterized(pVertex))
							pHalfedge->uv(pVertex->u(),pVertex->v());	
					}
					fprintf(stderr,"ok\n");
				}

				//
				// MESH INTERFACE
				//

				// Get iterator over first vertex of mesh
				Vertex_iterator  mesh_vertices_begin () {
					return m_mesh->vertices_begin();
				}
				Vertex_const_iterator  mesh_vertices_begin () const {
					return m_mesh->vertices_begin();
				}

				// Get iterator over past-the-end vertex of mesh
				Vertex_iterator  mesh_vertices_end () {
					return m_mesh->vertices_end();
				}
				Vertex_const_iterator  mesh_vertices_end () const {
					return m_mesh->vertices_end();
				}

				// Count the number of vertices of the mesh
				int  count_mesh_vertices () const {
					int index = 0;
					for (Vertex_const_iterator it = mesh_vertices_begin(); it != mesh_vertices_end(); it++) 
						index++;
					return index;
				}

				// Index vertices of the mesh for 0 to count_mesh_vertices()-1
				void  index_mesh_vertices () {
					int index = 0;
					for (Vertex_iterator it = mesh_vertices_begin(); it != mesh_vertices_end(); it++)
						set_vertex_index(it, index++);
				}

				// Return true of all mesh's faces are triangles
				bool  is_mesh_triangular () const {
					for (Face_const_iterator it = mesh_faces_begin(); it != mesh_faces_end(); it++) 
						if (count_face_vertices(it) != 3)
							return false;
					return true;						// mesh is triangular if we reach this point
				}

				// Compute the genus of the mesh
				int  get_mesh_genus () const {
					return m_mesh->genus();
				}

				// Count the number of boundaries of the mesh
				int  count_mesh_boundaries () const {
					return m_mesh->nb_boundaries();
				}

				// Get iterator over first vertex of mesh's border.
				Border_vertex_iterator  mesh_border_vertices_begin () {
					return m_boundary.begin();
				}
				Border_vertex_const_iterator  mesh_border_vertices_begin () const {
					return m_boundary.begin();
				}

				// Get iterator over past-the-end vertex of mesh's border.
				Border_vertex_iterator  mesh_border_vertices_end () {
					return m_boundary.end();
				}
				Border_vertex_const_iterator  mesh_border_vertices_end () const {
					return m_boundary.end();
				}

				// Get iterator over first face of mesh
				Face_iterator  mesh_faces_begin () {
					return m_mesh->facets_begin();
				}
				Face_const_iterator  mesh_faces_begin () const {
					return m_mesh->facets_begin();
				}

				// Get iterator over past-the-end face of mesh
				Face_iterator  mesh_faces_end () {
					return m_mesh->facets_end();
				}
				Face_const_iterator  mesh_faces_end () const {
					return m_mesh->facets_end();
				}

				// Count the number of faces of the mesh
				int  count_mesh_faces () const {
					int index = 0;
					for (Face_const_iterator it = mesh_faces_begin(); it != mesh_faces_end(); it++) 
						index++;
					return index;
				}

				//
				// FACE INTERFACE
				//

				// Get circulator over face's vertices
				Vertex_around_face_circulator  face_vertices_begin (Face_handle face) {
					return face->facet_begin();
				}
				Vertex_around_face_const_circulator  face_vertices_begin (Face_const_handle face) const {
					return face->facet_begin();
				}

				// Count the number of vertices of a face
				int  count_face_vertices (Face_const_handle face) const {
					int index = 0;
					Vertex_around_face_const_circulator cir_begin = face_vertices_begin(face), 
														cir_end   = cir_begin;
					CGAL_For_all(cir_begin, cir_end) { 
						index++;
					}
					return index;
				}

				//
				// VERTEX INTERFACE
				//

				// Get the 3D position of a vertex
				Point_3  get_vertex_position (Vertex_const_handle vertex) const {
					return vertex->point();	
				}

				// Get the 2D position of a vertex
				Point_2  get_vertex_uv (Vertex_const_handle vertex) const {
					return Point_2(vertex->u(), vertex->v());		
				}
				// Set the 2D position of a vertex
				void  set_vertex_uv (Vertex_handle vertex, const Point_2& uv) {
					vertex->u(uv.x());
					vertex->v(uv.y());
				}

				// Set "is parameterized" field of vertex
				void  set_vertex_parameterized (Vertex_handle vertex, bool parameterized) {
					vertex->halfedge()->is_parameterized(parameterized);
				}
				// Indicate if a vertex is already parameterized
				bool  is_vertex_parameterized (Vertex_const_handle vertex) const {
					return vertex->halfedge()->is_parameterized();
				}

				// Set vertex index
				void  set_vertex_index (Vertex_handle vertex, int new_index) {
					vertex->index(new_index);
				}
				// Get the index of a vertex
				int  get_vertex_index (Vertex_const_handle vertex) const {
					return vertex->index();		
				}

				// Return true if a vertex belongs to the mesh's boundary
				bool  is_vertex_on_border (Vertex_const_handle vertex) const {
					return m_mesh->is_border(vertex);
				}

				// Get circulator over the vertices incident to 'vertex'
				Vertex_around_vertex_circulator  vertices_around_vertex_begin (Vertex_handle vertex) {
					return vertex->vertex_begin();
				}
				Vertex_around_vertex_const_circulator  vertices_around_vertex_begin (Vertex_const_handle vertex) const {
					return vertex->vertex_begin();
				}

private:
				// Extract mesh UNIQUE boundary
				std::list<Halfedge_handle> extract_unique_boundary()
				{
					std::list<Halfedge_handle> boundary;
					m_mesh->free_skeleton();
					Mesh_feature_extractor feature_extractor(m_mesh);
					int nb_boundaries = feature_extractor.extract_boundaries(false,true);
					assert(nb_boundaries == 1);
					Backbone* pBackbone = (*m_mesh->get_skeleton()->backbones())[0];
					boundary = *(pBackbone->halfedges());
					m_mesh->free_skeleton();
					return boundary;
				}

// Fields
private:
				// The actual Polyhedron_ex mesh hidden inside *this
				Polyhedron_ex*	m_mesh;

				// Boundary of a topological disc inside the actual mesh
				std::list<Halfedge_handle> m_boundary;
};


#endif //MESHADAPTORPOLYHEDRONEX_H

