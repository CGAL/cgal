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
#include "Convertible_iterator.h"

#include <cassert>


// Class Mesh_adaptor_polyhedron_ex 
// Model of MeshAdaptor_3 concept, whose purpose is to allow the parameterization package to access meshes on an uniform manner
//
// Mesh_adaptor_polyhedron_ex is an adaptor class to access to a Polyhedron_ex 3D mesh using MeshAdaptor_3 interface. 
// The input mesh can be of any genus, but it has to come with a description of the boundary of a topological disc. If no boundary is given, 
// we assume that it coincides with the unique boundary already in the input mesh. 
// Mesh_adaptor_polyhedron_ex exports points and facets of a Polyhedron_ex mesh. Points and facets outside of the input boundary
// are not exported, i.e. are hidden to the parameterization package.
// Mesh_adaptor_polyhedron_ex stores 1 u/v pair per inner vertex (all halfedges of a given vertex share the same u/v) + 1 u/v pair per border halfedge
//
// Design pattern: 
// Mesh_adaptor_polyhedron_ex is an Adaptor (see [GOF95]): it changes the Polyhedron_ex interface to match the MeshAdaptor_3 concept
// 
// Implementation notes:
// * Mesh_adaptor_polyhedron_ex::Face	 = Polyhedron_ex::Facet (as you would expect)
// * Mesh_adaptor_polyhedron_ex::Vertex = Polyhedron_ex::Halfedge. Mesh_adaptor_polyhedron_ex represents inner vertices by their first incident halfedge and
//   border vertices by their halfedge on the inner side of the border. The fields expected by the MeshAdaptor_3 concept in vertices (u/v pair, 
//   index and "is_parameterized" boolean) are stored in halfedges.

class Mesh_adaptor_polyhedron_ex
{
// Private types
private:
				// Forward references
				typedef Polyhedron_ex::Halfedge								Vertex;	
				struct														Project_halfedge_vertex;
				struct														Project_halfedge_handle_vertex;
				struct														Project_opposite_halfedge_vertex;
				struct														Ignored;
				// Halfedge
				typedef Polyhedron_ex::Halfedge								Halfedge;
				typedef Polyhedron_ex::Halfedge_handle						Halfedge_handle;
				typedef Polyhedron_ex::Halfedge_const_handle				Halfedge_const_handle;
				typedef Polyhedron_ex::Halfedge_iterator					Halfedge_iterator;
				typedef Polyhedron_ex::Halfedge_const_iterator				Halfedge_const_iterator;
				typedef Polyhedron_ex::Halfedge_around_vertex_circulator	Halfedge_around_vertex_circulator;
				// Iterator over all mesh vertices
				typedef CGAL::Filter_iterator<Halfedge_iterator, Ignored>	Vertex_iterator_base;
				typedef CGAL::Filter_iterator<Halfedge_const_iterator, Ignored> Vertex_const_iterator_base;
				// Iterator over mesh boundary vertices
				typedef CGAL::Iterator_project<std::list<Halfedge_handle>::iterator, 
											   Project_halfedge_handle_vertex>	
																			Border_vertex_iterator_base;
				typedef CGAL::Iterator_project<std::list<Halfedge_handle>::const_iterator, 
											   Project_halfedge_handle_vertex>	
																			Border_vertex_const_iterator_base;
				// Circulator over a face's vertices
				typedef CGAL::Circulator_project<Polyhedron_ex::Halfedge_around_facet_circulator, 
												 Project_halfedge_vertex, Vertex&, Vertex*>					
																			Vertex_around_face_circulator_base;
				typedef CGAL::Circulator_project<Polyhedron_ex::Halfedge_around_facet_const_circulator, 
											     Project_halfedge_vertex, const Vertex&, const Vertex*>					
																			Vertex_around_face_const_circulator_base;
				// Circulator over the vertices incident to a vertex
				// @@@ INONDATION
				typedef CGAL::Circulator_project<Polyhedron_ex::Halfedge_around_vertex_circulator, 
												 Project_opposite_halfedge_vertex, Vertex&, Vertex*>					
																			Vertex_around_vertex_circulator_base;
				typedef CGAL::Circulator_project<Polyhedron_ex::Halfedge_around_vertex_const_circulator, 
												 Project_opposite_halfedge_vertex, const Vertex&, const Vertex*>					
																			Vertex_around_vertex_const_circulator_base;
	
// Public types
public:
				typedef Polyhedron_ex::Traits::FT							NT;
				// Points and vectors
				typedef Polyhedron_ex::Traits::Point_2						Point_2;
				typedef Polyhedron_ex::Traits::Point_3						Point_3;
				typedef Polyhedron_ex::Traits::Vector_2						Vector_2;
				typedef Polyhedron_ex::Traits::Vector_3						Vector_3;
				// Face
				typedef Polyhedron_ex::Facet								Face;
				typedef Polyhedron_ex::Facet_handle							Face_handle;
				typedef Polyhedron_ex::Facet_const_handle					Face_const_handle;
				// Iterator over all mesh faces
				typedef Polyhedron_ex::Facet_iterator						Face_iterator;
				typedef Polyhedron_ex::Facet_const_iterator					Face_const_iterator;
				// Mesh_adaptor_polyhedron_ex::Vertex = Polyhedron_ex::Halfedge
				typedef Polyhedron_ex::Halfedge								Vertex;	
				typedef Polyhedron_ex::Halfedge_handle						Vertex_handle;
				typedef Polyhedron_ex::Halfedge_const_handle				Vertex_const_handle;
				// Iterator over all mesh vertices
				typedef Convertible_iterator<Vertex_iterator_base, Vertex_const_handle, Vertex_handle>
																			Vertex_iterator;		
				typedef Convertible_iterator<Vertex_const_iterator_base, Vertex_const_handle>
																			Vertex_const_iterator;	
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
				// @@@ INONDATION
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
					// Extract mesh boundary
					std::list<Halfedge_handle> boundary = extract_unique_boundary(mesh);

					// Call common part of the constructors
					init(mesh, boundary.begin(), boundary.end());
				}
				// Create an adaptator for an existing Polyhedron_ex mesh
				// The input mesh can be of any genus, but it has to come with a description of the boundary of a topological disc 
				template <class InputIterator>
				Mesh_adaptor_polyhedron_ex(Polyhedron_ex* mesh,
										   InputIterator first_boundary_halfedge, InputIterator last_boundary_halfedge) 
				{
					// Call common part of the constructors
					init(mesh, first_boundary_halfedge, last_boundary_halfedge);
				}

				// Default copy constructor and operator =() are fine

				//
				// MESH INTERFACE
				//

				// Get iterator over first vertex of mesh
				Vertex_iterator  mesh_vertices_begin () {
					return Vertex_iterator_base(m_mesh->halfedges_end(), Ignored(), m_mesh->halfedges_begin());
				}
				Vertex_const_iterator  mesh_vertices_begin () const {
					return Vertex_const_iterator_base(m_mesh->halfedges_end(), Ignored(), m_mesh->halfedges_begin());
				}

				// Get iterator over past-the-end vertex of mesh
				Vertex_iterator  mesh_vertices_end () {
					return Vertex_iterator_base(m_mesh->halfedges_end(), Ignored());
				}
				Vertex_const_iterator  mesh_vertices_end () const {
					return Vertex_const_iterator_base(m_mesh->halfedges_end(), Ignored());
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
					return (Border_vertex_iterator_base) m_boundary.begin();
				}
				Border_vertex_const_iterator  mesh_border_vertices_begin () const {
					return (Border_vertex_const_iterator_base) m_boundary.begin();
				}

				// Get iterator over past-the-end vertex of mesh's border.
				Border_vertex_iterator  mesh_border_vertices_end () {
					return (Border_vertex_iterator_base) m_boundary.end();
				}
				Border_vertex_const_iterator  mesh_border_vertices_end () const {
					return (Border_vertex_const_iterator_base) m_boundary.end();
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
					assert(is_valid(face));
					return (Vertex_around_face_circulator_base) face->facet_begin();
				}
				Vertex_around_face_const_circulator  face_vertices_begin (Face_const_handle face) const {
					assert(is_valid(face));
					return (Vertex_around_face_const_circulator_base) face->facet_begin();
				}

				// Count the number of vertices of a face
				int  count_face_vertices (Face_const_handle face) const {
					assert(is_valid(face));
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

				// Get the 3D position of a vertex (stored in Polyhedron vertex)
				Point_3  get_vertex_position (Vertex_const_handle adaptor_vertex) const {
					assert(is_valid(adaptor_vertex));
					return adaptor_vertex->vertex()->point();	
				}

				// Get/set the 2D position (u/v pair) of a vertex (stored in halfedge)
				Point_2  get_vertex_uv (Vertex_const_handle adaptor_vertex) const {
					assert(is_valid(adaptor_vertex));
					return Point_2(adaptor_vertex->u(), adaptor_vertex->v());		
				}
				void  set_vertex_uv (Vertex_handle adaptor_vertex, const Point_2& uv) 
				{
					assert(is_valid(adaptor_vertex));
					// Update all halfedges sharing the same vertex (except outer halfedges)
					Polyhedron_ex::Vertex_handle polyhedron_vertex = adaptor_vertex->vertex();
					Halfedge_around_vertex_circulator cir     = polyhedron_vertex->vertex_begin(), 
													  cir_end = cir;
					CGAL_For_all(cir, cir_end) { 
						if ( ! cir->is_border() ) {	// skip outer halfedges @@@ SEAM
							cir->u(uv.x());	
							cir->v(uv.y());
						}
					}
				}

				// Get/set "is parameterized" field of vertex (stored in halfedge)
				bool  is_vertex_parameterized (Vertex_const_handle adaptor_vertex) const {
					assert(is_valid(adaptor_vertex));
					return adaptor_vertex->is_parameterized();
				}
				void  set_vertex_parameterized (Vertex_handle adaptor_vertex, bool parameterized) 
				{
					assert(is_valid(adaptor_vertex));
					// Update all halfedges sharing the same vertex (except outer halfedges)
					Polyhedron_ex::Vertex_handle polyhedron_vertex = adaptor_vertex->vertex();
					Halfedge_around_vertex_circulator cir     = polyhedron_vertex->vertex_begin(), 
													  cir_end = cir;
					CGAL_For_all(cir, cir_end) { 
						if ( ! cir->is_border() )	// skip outer halfedges @@@ SEAM
							cir->is_parameterized(parameterized);
					}
				}

				// Get/set vertex index (stored in halfedge)
				int  get_vertex_index (Vertex_const_handle adaptor_vertex) const {
					assert(is_valid(adaptor_vertex));
					return adaptor_vertex->index();		
				}
				void  set_vertex_index (Vertex_handle adaptor_vertex, int new_index) 
				{
					assert(is_valid(adaptor_vertex));
					// Update all halfedges sharing the same vertex (except outer halfedges)
					Polyhedron_ex::Vertex_handle polyhedron_vertex = adaptor_vertex->vertex();
					Halfedge_around_vertex_circulator cir     = polyhedron_vertex->vertex_begin(), 
													  cir_end = cir;
					CGAL_For_all(cir, cir_end) { 
						if ( ! cir->is_border() )	// skip outer halfedges @@@ SEAM
							cir->index(new_index);
					}
				}

				// Return true if a vertex belongs to the mesh's boundary
				bool  is_vertex_on_border (Vertex_const_handle adaptor_vertex) const {
					assert(is_valid(adaptor_vertex));
					return adaptor_vertex->opposite()->is_border();		// @@@ SEAM
				}

				// Get circulator over the vertices incident to 'adaptor_vertex'
				// @@@ INONDATION
				Vertex_around_vertex_circulator  vertices_around_vertex_begin (Vertex_handle adaptor_vertex) {
					assert(is_valid(adaptor_vertex));
					return (Vertex_around_vertex_circulator) adaptor_vertex->vertex_begin();
				}
				Vertex_around_vertex_const_circulator  vertices_around_vertex_begin (Vertex_const_handle adaptor_vertex) const {
					assert(is_valid(adaptor_vertex));
					return (Vertex_around_vertex_const_circulator) adaptor_vertex->vertex_begin();
				}

// Private types
private:
				// Utility class to generate the Vertex_around_face_circulator type
				struct Project_halfedge_vertex {
					typedef Halfedge							argument_type;
					typedef Mesh_adaptor_polyhedron_ex::Vertex	Vertex;
					typedef Vertex								result_type;
					typedef CGAL::Arity_tag<1>					Arity;
					// Get the adaptor vertex of the halfedge 'h'
					Vertex&       operator()(Halfedge& h)       const { return *(get_adaptor_vertex(&h)); }
					const Vertex& operator()(const Halfedge& h) const { return *(get_adaptor_vertex(&h)); }
				};

				// Utility class to generate the Border_vertex_iterator type
				struct Project_halfedge_handle_vertex {
					typedef Halfedge_handle						argument_type;
					typedef Mesh_adaptor_polyhedron_ex::Vertex	Vertex;
					typedef Vertex								result_type;
					typedef CGAL::Arity_tag<1>					Arity;
					// Get the adaptor vertex of the halfedge handle 'h'
					Vertex&       operator()(Halfedge_handle& h)       const { return *(get_adaptor_vertex(h)); }
					const Vertex& operator()(const Halfedge_handle& h) const { return *(get_adaptor_vertex(h)); }
				};

				// Utility class to generate the Vertex_around_vertex_circulator type
				struct Project_opposite_halfedge_vertex {
					typedef Halfedge							argument_type;
					typedef Mesh_adaptor_polyhedron_ex::Vertex	Vertex;
					typedef Vertex								result_type;
					typedef CGAL::Arity_tag<1>					Arity;
					// Get the adaptor vertex of the halfedge opposite to h
					Vertex&       operator()(Halfedge& h)       const { return *(get_adaptor_vertex(h.opposite())); }
					const Vertex& operator()(const Halfedge& h) const { return *(get_adaptor_vertex(h.opposite())); }
				};

				// Utility class to generate the Vertex_iterator type
				struct Ignored {
					// Return true <=> the object is not exported by Mesh_adaptor_polyhedron_ex
					bool operator()(const Halfedge_iterator& h) const		{ return get_adaptor_vertex(h) != h; }	// @@@ INONDATION
					bool operator()(const Halfedge_const_iterator& h) const { return get_adaptor_vertex(h) != h; }	// @@@ INONDATION
					bool operator()(const Face_iterator& f) const			{ return false; }						// @@@ INONDATION
					bool operator()(const Face_const_iterator& f) const		{ return false; }						// @@@ INONDATION
				};

 // Private operations
private:
				// Extract mesh UNIQUE boundary
				static std::list<Halfedge_handle> extract_unique_boundary(Polyhedron_ex* mesh)
				{
					typedef Feature_backbone<Polyhedron_ex::Vertex_handle, Halfedge_handle>	Backbone;
					assert(mesh != NULL);
					std::list<Halfedge_handle> boundary;	// returned list
					mesh->free_skeleton();
					Mesh_feature_extractor feature_extractor(mesh);
					int nb_boundaries = feature_extractor.extract_boundaries(true);
					assert(nb_boundaries == 1);
					Backbone *pBackbone = (*mesh->get_skeleton()->backbones())[0];
					boundary = *(pBackbone->halfedges());
					mesh->free_skeleton();
					return boundary;
				}

				// Utility method that contains the common part of the constructors:
				// Initialize an adaptator for an existing Polyhedron_ex mesh
				// The input mesh can be of any genus, but it has to come with a description of the boundary of a topological disc 
				template <class InputIterator>
				void init(Polyhedron_ex* mesh, InputIterator first_boundary_halfedge, InputIterator last_boundary_halfedge) 
				{
					assert(mesh != NULL);
					m_mesh = mesh;

#ifndef NDEBUG
					// Index Polyhedron_ex vertices to ease debugging
					int index = 0;
					for (Polyhedron_ex::Vertex_iterator it2 = mesh->vertices_begin(); it2 != mesh->vertices_end(); it2++)
						it2->index(index++);

					// Index all halfedges with negative numbers to ease debugging
					/*int*/ index = -1;
					for (Halfedge_iterator it3 = mesh->halfedges_begin(); it3 != mesh->halfedges_end(); it3++)
						it3->index(index--);
#endif

					// Copy input boundary
					assert(m_boundary.empty());
					for (InputIterator it = first_boundary_halfedge; it != last_boundary_halfedge; it++)
						m_boundary.push_back(*it);

					// @@@ INONDER LA SURFACE ICI
				}

				// Convert any polyhedron halfedge to its adaptor vertex 
				// Note: there is no polyhedron vertex to adaptor vertex conversion because a vertex may belong several times to the seam
				//
				// Implementation note: Mesh_adaptor_polyhedron_ex::Vertex = Polyhedron_ex::Halfedge. Mesh_adaptor_polyhedron_ex represents 
				//                      inner vertices by their first incident halfedge and border vertices by their halfedge on the inner side of the border.
				static Vertex_const_handle get_adaptor_vertex(Halfedge_const_handle halfedge)
				{
					assert(halfedge != NULL);

					Vertex_const_handle adaptor_vertex;					// returned value

					// If the polyhedron vertex is on the boundary
					Polyhedron_ex::Vertex_const_handle polyhedron_vertex = halfedge->vertex();
					if (Polyhedron_ex::is_border(polyhedron_vertex)) 	// @@@ SEAM
					{
						// do a counter clockwise search around the vertex until we reach the inner boundary halfedge
						adaptor_vertex = halfedge;
						while ( ! adaptor_vertex->opposite()->is_border() )
							adaptor_vertex = adaptor_vertex->opposite()->prev();
					} 
					else // if inner vertex, get its first incident halfedge
					{
						adaptor_vertex = polyhedron_vertex->vertex_begin();
					}

					assert(adaptor_vertex != NULL);
					assert(adaptor_vertex->vertex() == halfedge->vertex());
					return adaptor_vertex;
				}
				inline static Vertex_handle get_adaptor_vertex(Halfedge_handle halfedge) {
					Vertex_const_handle adaptor_vertex = get_adaptor_vertex( (Halfedge_const_handle)halfedge );
					return (Vertex*) (&*adaptor_vertex);
				}
				inline static Vertex_handle get_adaptor_vertex(Halfedge* halfedge) {
					return get_adaptor_vertex( (Halfedge_handle)halfedge );
				}

				// Check if variables are valid (for debug purpose)
				inline bool is_valid (Face_const_handle face) const {
					return (face != NULL);	
				}
				inline bool is_valid (Vertex_const_handle adaptor_vertex) const {
					if (adaptor_vertex == NULL)
						return false;
					if ( adaptor_vertex->is_border() )					// outer halfedges are not exported @@@ SEAM
						return false;
					return (get_adaptor_vertex(adaptor_vertex) == adaptor_vertex);	
				}

// Fields
private:
				Polyhedron_ex*	m_mesh;							// The adapted mesh 
				std::list<Halfedge_handle> m_boundary;			// Inner halfedges of the boundary of a topological disc inside m_mesh
};


#endif //MESHADAPTORPOLYHEDRONEX_H

