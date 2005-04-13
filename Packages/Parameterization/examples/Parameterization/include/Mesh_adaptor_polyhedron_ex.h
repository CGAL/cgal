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
#include <CGAL/Filter_circulator.h>
#include <CGAL/function_objects.h>
#include <CGAL/HalfedgeDS_iterator.h>

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
// * Mesh_adaptor_polyhedron_ex::Face	= Polyhedron_ex::Facet (as you would expect)
// * Mesh_adaptor_polyhedron_ex::Vertex = Polyhedron_ex::Halfedge. Mesh_adaptor_polyhedron_ex represents inner vertices by their first incident halfedge and
//   border vertices by their halfedge on the inner side of the border. The fields expected by the MeshAdaptor_3 concept in vertices (u/v pair, 
//   index and "is_parameterized" boolean) are stored in halfedges.

class Mesh_adaptor_polyhedron_ex
{
// Private types
private:
	// Forward references
	struct															Project_halfedge_vertex;
	struct															Project_halfedge_handle_vertex;
	struct															All_faces_filter;
	struct															All_vertices_filter;
	template <class It, class Ctg> class							Adaptor_vertex_around_vertex_circulator;
	// Halfedge
	typedef Polyhedron_ex::Halfedge									Halfedge;
	typedef Polyhedron_ex::Halfedge_handle							Halfedge_handle;
	typedef Polyhedron_ex::Halfedge_const_handle					Halfedge_const_handle;
	typedef Polyhedron_ex::Halfedge_iterator						Halfedge_iterator;
	typedef Polyhedron_ex::Halfedge_const_iterator					Halfedge_const_iterator;
	typedef Polyhedron_ex::Halfedge_around_vertex_circulator		Halfedge_around_vertex_circulator;
	typedef Polyhedron_ex::Halfedge_around_vertex_const_circulator	Halfedge_around_vertex_const_circulator;
	typedef Polyhedron_ex::Halfedge_around_facet_circulator			Halfedge_around_facet_circulator;
	typedef Polyhedron_ex::Halfedge_around_facet_const_circulator	Halfedge_around_facet_const_circulator;
    typedef Halfedge::Supports_halfedge_prev						Supports_prev;
	// Seaming flag for faces, vertices and halfedges
	enum Seaming_status	 { OUTER, INNER, BORDER /* = outer border for HE */};

// Public types
public:
	typedef Polyhedron_ex::Traits::FT								NT;
	// Points and vectors
	typedef Polyhedron_ex::Traits::Point_2							Point_2;
	typedef Polyhedron_ex::Traits::Point_3							Point_3;
	typedef Polyhedron_ex::Traits::Vector_2							Vector_2;
	typedef Polyhedron_ex::Traits::Vector_3							Vector_3;
	// Face
	typedef Polyhedron_ex::Facet									Face;
	typedef Polyhedron_ex::Facet_handle								Face_handle;
	typedef Polyhedron_ex::Facet_const_handle						Face_const_handle;
    // Circulator category
	typedef CGAL::HalfedgeDS_circulator_traits<Supports_prev>		Ctr;
    typedef Ctr::iterator_category									circulator_category;
	// Iterator over all mesh faces
	typedef CGAL::Filter_iterator<Polyhedron_ex::Facet_iterator, All_faces_filter>	
																	Face_iterator_base;
	typedef Convertible_iterator<Face_iterator_base, Face_const_handle, Face_handle>
																	Face_iterator;		
	typedef CGAL::Filter_iterator<Polyhedron_ex::Facet_const_iterator, All_faces_filter> 
																	Face_const_iterator_base;
	typedef Convertible_iterator<Face_const_iterator_base, Face_const_handle>
																	Face_const_iterator;	
	// Mesh_adaptor_polyhedron_ex::Vertex = Polyhedron_ex::Halfedge
	typedef Polyhedron_ex::Halfedge									Vertex;	
	typedef Polyhedron_ex::Halfedge_handle							Vertex_handle;
	typedef Polyhedron_ex::Halfedge_const_handle					Vertex_const_handle;
	// Iterator over all mesh vertices
	typedef CGAL::Filter_iterator<Halfedge_iterator, All_vertices_filter> 
																	Vertex_iterator_base;
	typedef Convertible_iterator<Vertex_iterator_base, Vertex_const_handle, Vertex_handle>
																	Vertex_iterator;		
	typedef CGAL::Filter_iterator<Halfedge_const_iterator, All_vertices_filter> 
																	Vertex_const_iterator_base;
	typedef Convertible_iterator<Vertex_const_iterator_base, Vertex_const_handle>
																	Vertex_const_iterator;	
	// Iterator over mesh boundary vertices
	typedef CGAL::Iterator_project<std::list<Halfedge_handle>::iterator, 
								   Project_halfedge_handle_vertex>	Border_vertex_iterator_base;
	typedef Convertible_iterator<Border_vertex_iterator_base, Vertex_const_handle, Vertex_handle>
																	Border_vertex_iterator;		
	typedef CGAL::Iterator_project<std::list<Halfedge_handle>::const_iterator, 
								   Project_halfedge_handle_vertex>	Border_vertex_const_iterator_base;
	typedef Convertible_iterator<Border_vertex_const_iterator_base, Vertex_const_handle>
																	Border_vertex_const_iterator;	
	// Circulator over a face's vertices
	typedef CGAL::Circulator_project<Halfedge_around_facet_circulator, Project_halfedge_vertex,  
									 Vertex&, Vertex*>				Vertex_around_face_circulator_base;
	typedef Convertible_iterator<Vertex_around_face_circulator_base, Vertex_const_handle, Vertex_handle>
																	Vertex_around_face_circulator;	
	typedef CGAL::Circulator_project<Halfedge_around_facet_const_circulator, Project_halfedge_vertex, 
									 const Vertex&, const Vertex*>	Vertex_around_face_const_circulator_base;
	typedef Convertible_iterator<Vertex_around_face_const_circulator_base, Vertex_const_handle>
																	Vertex_around_face_const_circulator; 
	// Circulator over the vertices incident to a vertex
	typedef Adaptor_vertex_around_vertex_circulator<Vertex_handle, circulator_category>
																	Vertex_around_vertex_circulator;	
	typedef Adaptor_vertex_around_vertex_circulator<Vertex_const_handle, circulator_category>
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
		// Extract mesh's longest boundary
		std::list<Halfedge_handle> boundary = extract_longest_boundary(mesh);

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
		return Vertex_iterator_base(m_mesh->halfedges_end(), All_vertices_filter(), m_mesh->halfedges_begin());
	}
	Vertex_const_iterator  mesh_vertices_begin () const {
		return Vertex_const_iterator_base(m_mesh->halfedges_end(), All_vertices_filter(), m_mesh->halfedges_begin());
	}

	// Get iterator over past-the-end vertex of mesh
	Vertex_iterator  mesh_vertices_end () {
		return Vertex_iterator_base(m_mesh->halfedges_end(), All_vertices_filter());
	}
	Vertex_const_iterator  mesh_vertices_end () const {
		return Vertex_const_iterator_base(m_mesh->halfedges_end(), All_vertices_filter());
	}

	// Get the number of vertices of the mesh
	int  count_mesh_vertices () const {
		return m_mesh_vertices_count;
	}

	// Index vertices of the mesh for 0 to count_mesh_vertices()-1
	void  index_mesh_vertices () 
	{
		fprintf(stderr,"  index Mesh_adaptor_polyhedron_ex vertices... ");
		int index = 0;
		for (Vertex_iterator it = mesh_vertices_begin(); it != mesh_vertices_end(); it++)
			set_vertex_index(it, index++);
		fprintf(stderr,"ok\n");

//#ifndef NDEBUG
//		// Print halfedges 
//		fprintf(stderr,"  print all halfedges... ");
//		for (Halfedge_iterator he_it = m_mesh->halfedges_begin(); he_it != m_mesh->halfedges_end(); he_it++)
//			fprintf(stderr, "HE %d=%d->%d ", (int)he_it->index(), (int)he_it->opposite()->vertex()->index(), (int)he_it->vertex()->index());
//		fprintf(stderr,"ok\n");
//#endif
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
		// the exported mesh is a topological disc:
		return 0; 
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
		return Face_iterator_base(m_mesh->facets_end(), All_faces_filter(), m_mesh->facets_begin());
	}
	Face_const_iterator  mesh_faces_begin () const {
		return Face_const_iterator_base(m_mesh->facets_end(), All_faces_filter(), m_mesh->facets_begin());
	}

	// Get iterator over past-the-end face of mesh
	Face_iterator  mesh_faces_end () {
		return Face_iterator_base(m_mesh->facets_end(), All_faces_filter());
	}
	Face_const_iterator  mesh_faces_end () const {
		return Face_const_iterator_base(m_mesh->facets_end(), All_faces_filter());
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

	// Get/set the 2D position (u/v pair) of a vertex (stored in halfedges sharing the same vertex)
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
			if (get_adaptor_vertex(cir) == adaptor_vertex) {	// if on the same side of a seam
				if (is_inner_or_inner_border(cir)) {			// skip outer halfedges
					assert(Halfedge_handle(cir) == adaptor_vertex || cir->index() < 0);
					cir->u(uv.x());	
					cir->v(uv.y());
				}
			}
		}
	}

	// Get/set "is parameterized" field of vertex (stored in halfedges sharing the same vertex)
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
			if (get_adaptor_vertex(cir) == adaptor_vertex) {	// if on the same side of a seam
				if (is_inner_or_inner_border(cir)) {			// skip outer halfedges
					assert(Halfedge_handle(cir) == adaptor_vertex || cir->index() < 0);
					cir->is_parameterized(parameterized);
				}
			}
		}
	}

	// Get/set vertex index (stored in "main" halfedge)
	int  get_vertex_index (Vertex_const_handle adaptor_vertex) const {
		assert(is_valid(adaptor_vertex));
		int index = adaptor_vertex->index();		
		assert(index >= 0 && index < m_mesh_vertices_count);
		return index;
	}
	void  set_vertex_index (Vertex_handle adaptor_vertex, int index)  {
		assert(is_valid(adaptor_vertex));
		assert(index >= 0 && index < m_mesh_vertices_count);
		adaptor_vertex->index(index);	
	}

	// Return true if a vertex belongs to the mesh's boundary
	bool  is_vertex_on_border (Vertex_const_handle adaptor_vertex) const {
		assert(is_valid(adaptor_vertex));
		// Border vertices are exported as inner border halfedges. Halfedges marked BORDER are outer border ones.
		return (adaptor_vertex->opposite()->seaming() == BORDER);
	}

	// Get circulator over the vertices incident to 'adaptor_vertex'
	Vertex_around_vertex_circulator  vertices_around_vertex_begin (Vertex_handle adaptor_vertex) {
		assert(is_valid(adaptor_vertex));
		return Vertex_around_vertex_circulator(adaptor_vertex->next(), adaptor_vertex);
	}
	Vertex_around_vertex_const_circulator  vertices_around_vertex_begin (Vertex_const_handle adaptor_vertex) const {
		assert(is_valid(adaptor_vertex));
		return Vertex_around_vertex_const_circulator(adaptor_vertex->next(), adaptor_vertex);
	}

// Private operations
private:
	// Extract mesh's longest boundary
	static std::list<Halfedge_handle> extract_longest_boundary(Polyhedron_ex* mesh)
	{
		typedef Feature_backbone<Polyhedron_ex::Vertex_handle, Halfedge_handle>	Backbone;
		assert(mesh != NULL);
		std::list<Halfedge_handle> boundary;	// returned list

		// Extract mesh boundaries
		mesh->free_skeleton();
		Mesh_feature_extractor feature_extractor(mesh);
		int nb_boundaries = feature_extractor.extract_boundaries(true);

		// Get longest one
		assert(nb_boundaries > 0);
		Backbone *pBackbone = (*mesh->get_skeleton()->backbones())[0];
		boundary = *(pBackbone->halfedges());
		mesh->free_skeleton();
		return boundary;
	}

	// Utility method that contains the common part of the constructors:
	// Initialize an adaptator for an existing Polyhedron_ex mesh
	// The input mesh can be of any genus, but it has to come with a description of the boundary of a topological disc 
	//
	// Implementation note: at this stage, we cannot assume that the mesh is a triangular surface
	template <class InputIterator>
	void init(Polyhedron_ex* mesh, InputIterator first_boundary_halfedge, InputIterator last_boundary_halfedge) 
	{
		assert(mesh != NULL);
		m_mesh = mesh;
		m_mesh_vertices_count = 0;	// not yet initialized

#ifndef NDEBUG
		// Index Polyhedron_ex vertices to ease debugging
		fprintf(stderr,"  index Polyhedron vertices: ");
		int vtx_index = 0;
		for (Polyhedron_ex::Vertex_iterator vtx_it = m_mesh->vertices_begin(); vtx_it != m_mesh->vertices_end(); vtx_it++)
		{
			//fprintf(stderr, "#%d=(%f,%f,%f) ", (int)vtx_index, (float)vtx_it->point().x(), (float)vtx_it->point().y(), (float)vtx_it->point().z());
			vtx_it->index(vtx_index++);
		}
		fprintf(stderr,"ok\n");

		// Index all halfedges with negative numbers to ease debugging
		fprintf(stderr,"  index all halfedges with negative numbers: ");
		int he_index = -1;
		for (Halfedge_iterator he_it = m_mesh->halfedges_begin(); he_it != m_mesh->halfedges_end(); he_it++)
		{
			//fprintf(stderr, "HE %d=%d->%d ", (int)he_index, (int)he_it->opposite()->vertex()->index(), (int)he_it->vertex()->index());
			he_it->index(he_index--);
		}
		fprintf(stderr,"ok\n");

		// Dump input outer boundary (for debug purpose)
		//fprintf(stderr,"  input boundary is: ");
		//for (InputIterator border_it = first_boundary_halfedge; border_it != last_boundary_halfedge; border_it++)
		//	fprintf(stderr, "%d->%d ", (int)(*border_it)->opposite()->vertex()->index(), (int)(*border_it)->vertex()->index());
		//fprintf(stderr,"ok\n");
#endif

		// Check that the input boundary is a loop => it "cuts" a topological disc inside m_mesh
		bool input_boundary_is_valid = (first_boundary_halfedge != last_boundary_halfedge);
		assert(input_boundary_is_valid);
		for (InputIterator border_it = first_boundary_halfedge; border_it != last_boundary_halfedge; border_it++)
		{
			// Get next halfedge iterator 
			InputIterator next_boundary_it = border_it;
			next_boundary_it++;
			if (next_boundary_it == last_boundary_halfedge)
				next_boundary_it = first_boundary_halfedge;
			// Check that end of current HE == start of next one
			input_boundary_is_valid = ((*border_it)->vertex() == (*next_boundary_it)->opposite()->vertex());
			assert(input_boundary_is_valid);
		}
		// TO DO: check that the input boundary is not self-intersecting

		// Copy input outer boundary
		assert(m_boundary.empty());
		for (InputIterator border_it = first_boundary_halfedge; border_it != last_boundary_halfedge; border_it++)
			m_boundary.push_back(*border_it);

		// Set seaming flag of all halfedges to INNER, BORDER and OUTER wrt the boundary m_boundary
		flag_halfedges_seaming();

		// Count the number of vertices of the mesh (must be done after flagging the halfedges)
		m_mesh_vertices_count = 0;
		for (Vertex_iterator it = mesh_vertices_begin(); it != mesh_vertices_end(); it++) 
			m_mesh_vertices_count++;
	}

	// Set seaming flag of all halfedges to INNER, BORDER and OUTER wrt the boundary m_boundary
	//
	// Preconditions: 
	// * m_boundary must contain the outer boundary
	// Posconditions: 
	// * outer border halfedges are marked BORDER
	// * inner border halfedges are marked INNER, except if they are outer border for the other side of a seam
	// * outer halfedges are marked OUTER
	void flag_halfedges_seaming()
	{
		fprintf(stderr, "  tag topological disc...");

		// Initialize the seaming flag of all halfedges to OUTER
		m_mesh->flag_halfedges_seaming(OUTER);

		// Set seaming flag of outer boundary halfedges to BORDER
		assert( ! m_boundary.empty() );
		for (std::list<Halfedge_handle>::iterator it = m_boundary.begin(); it != m_boundary.end(); it++)
			(*it)->seaming(BORDER);

		// Set the seaming flag of the halfedges of the topological disc defined surrounded by m_boundary as INNER
		flag_surface_halfedges_seaming((*m_boundary.begin())->opposite()->facet());

		fprintf(stderr,"ok\n");
	}

	// Set the seaming flag of the topological disc defined by the boundary m_boundary as INNER
	//
	// Preconditions: 
	// * inner halfedges are marked as OUTER, outer border halfedges as BORDER
	// * pSeedFacet != NULL
	// Posconditions: 
	// * outer border halfedges are marked BORDER
	// * inner border halfedges are marked INNER, except if they are outer border for the other side of a seam
	// * outer halfedges are marked OUTER
	//
	// Implementation note: at this stage, we cannot assume that the mesh triangular
	void flag_surface_halfedges_seaming(Face_handle pSeedFacet)
	{
		if (pSeedFacet == NULL)
			return;							// Gloups... topological disc is empty!

		// List of triangles to flag = pSeedFacet initially
		std::list<Face_handle> facets;
		facets.push_front(pSeedFacet);

		// For each face in the list: pop it out, flag it as INNER and 
		// add its surrounding triangles to the list
		while(!facets.empty())
		{
			Face_handle pFacet = facets.front();
			facets.pop_front();
			assert(pFacet != NULL);

			// Flag this triangle's halfedges as INNER (except outer border ones)
			bool already_done = true;
			Halfedge_around_facet_circulator pHalfedge = pFacet->facet_begin(),
											 end       = pHalfedge;
			CGAL_For_all(pHalfedge, end)
			{
				assert(pHalfedge != NULL);
				if (pHalfedge->seaming() == OUTER)
				{
					pHalfedge->seaming(INNER);
					already_done = false;
				}
			}

			// Skip this triangle if it is done
			if (already_done)
				continue;

			// Add surrounding triangles to list without crossing the border
			pHalfedge = pFacet->facet_begin();
			end       = pHalfedge;
			CGAL_For_all(pHalfedge, end)
			{
				Halfedge_handle opposite_he = pHalfedge->opposite();
				assert(opposite_he != NULL);
				if (opposite_he->seaming() != BORDER)
				{
					// If the neighbor face is a hole, flag the neighbor outer vertex as INNER
					// Else, add the face to the list to flag it later
					Face_handle neighbor_face = opposite_he->facet();
					if (neighbor_face == NULL)
						opposite_he->seaming(INNER);
					else
						facets.push_front(neighbor_face);
				}
			}
		}
	}

	// Get face seaming status (INNER or OUTER)
	static Seaming_status seaming(Polyhedron_ex::Facet_const_handle face)
	{
		assert(face != NULL);
		Halfedge_const_handle halfedge = face->halfedge();
		assert(halfedge != NULL);
		return is_inner_or_inner_border(halfedge) ? INNER : OUTER;
	}

	// Get Polyhedron vertex seaming status
	// Implementation note: we could save time by adding a seaming flag to Polyhedron vertex class
	//                      instead of recomputing it over and over
	static Seaming_status seaming(Polyhedron_ex::Vertex_const_handle polyhedron_vertex)
	{
		Halfedge_around_vertex_const_circulator pHalfedge = polyhedron_vertex->vertex_begin(),
					                            end       = pHalfedge;
		
		// if isolated vertex
		if(pHalfedge == NULL) 
			return BORDER;
		
		// Parse the seaming status of the halfedges incident to 'polyhedron_vertex'. 
		// Priority is BORDER > INNER > OUTER.
		Seaming_status seaming = OUTER;
		CGAL_For_all(pHalfedge,end) 
		{
			Seaming_status current_he_seaming  = (Seaming_status) pHalfedge->seaming(),
				           opposite_he_seaming = (Seaming_status) pHalfedge->opposite()->seaming();

			// Valid tags are:
			assert( (current_he_seaming == INNER  && opposite_he_seaming == BORDER)
				 || (current_he_seaming == BORDER && opposite_he_seaming == INNER)
				 || (current_he_seaming == INNER  && opposite_he_seaming == INNER)
				 || (current_he_seaming == BORDER && opposite_he_seaming == BORDER)
				 || (current_he_seaming == OUTER  && opposite_he_seaming == OUTER) );

			if (current_he_seaming == INNER)
				seaming = INNER;
			if (opposite_he_seaming == BORDER)
				return BORDER;
		}
		return seaming;
	}

	// Indicate if an halfedge is in the topological disk (including its inner boundary)
	inline static bool is_inner_or_inner_border(Halfedge_const_handle halfedge) 
	{
		assert(halfedge != NULL);

		// Inner border halfedges are marked INNER (general case) or BORDER 
		// (if it is on a seam that virtually "cuts" the shape)
		if (halfedge->opposite()->seaming() == BORDER) {
			assert(halfedge->seaming() == INNER
			    || halfedge->seaming() == BORDER);
		}

		return (halfedge->seaming() == INNER || halfedge->opposite()->seaming() == BORDER);	
	}

	// Convert any polyhedron halfedge to its adaptor vertex 
	// Note: there is no polyhedron vertex to adaptor vertex conversion because a vertex may belong several times to the seam
	//
	// Implementation note: Mesh_adaptor_polyhedron_ex::Vertex = Polyhedron_ex::Halfedge. Mesh_adaptor_polyhedron_ex represents 
	//                      inner vertices by their first incident halfedge and border vertices by their halfedge on the inner side of the border.
	static Vertex_const_handle get_adaptor_vertex(Halfedge_const_handle halfedge)
	{
		assert(halfedge != NULL);

		Vertex_const_handle adaptor_vertex = NULL;				// returned value

		// If the polyhedron vertex is on the boundary
		Polyhedron_ex::Vertex_const_handle polyhedron_vertex = halfedge->vertex();
		if (seaming(polyhedron_vertex) == BORDER) 
		{
			// do a counterclockwise search around the vertex until we reach the inner boundary halfedge
			adaptor_vertex = halfedge;
			while (adaptor_vertex->opposite()->seaming() != BORDER)
				adaptor_vertex = adaptor_vertex->opposite()->prev();
		} 
		// else if inner vertex, get its first incident halfedge
		else if (halfedge->seaming() == INNER)
		{
			adaptor_vertex = polyhedron_vertex->vertex_begin();
		}
		else 
		{
			assert(halfedge->seaming() == OUTER);
		}

		assert( (adaptor_vertex == NULL) || (adaptor_vertex->vertex() == halfedge->vertex()) );
		return adaptor_vertex;
	}
	inline static Vertex_handle get_adaptor_vertex(Halfedge_handle halfedge) {
		Vertex_const_handle adaptor_vertex = get_adaptor_vertex( (Halfedge_const_handle)halfedge );
		return (Vertex*) (&*adaptor_vertex);
	}

	// Same as get_adaptor_vertex() + checks that result is a valid 
	// Mesh_adaptor_polyhedron_ex vertex (for debug purpose)
	inline static Vertex_const_handle get_safe_adaptor_vertex(Halfedge_const_handle halfedge) {
		Vertex_const_handle adaptor_vertex = get_adaptor_vertex(halfedge);
		assert(is_valid(adaptor_vertex));
		return adaptor_vertex;
	}
	inline static Vertex_handle get_safe_adaptor_vertex(Halfedge_handle halfedge) {
		Vertex_handle adaptor_vertex = get_adaptor_vertex(halfedge);
		assert(is_valid(adaptor_vertex));
		return adaptor_vertex;
	}

	// Check if a Mesh_adaptor_polyhedron_ex face is valid (for debug purpose)
	inline static bool is_valid (Face_const_handle face) {
		if (face == NULL)
			return false;
		// outer triangles are not exported
		if (seaming(face) != INNER)					
			return false;
		// eventually: ok
		return true;
	}

	// Check if a Mesh_adaptor_polyhedron_ex vertex is valid (for debug purpose)
	inline static bool is_valid (Vertex_const_handle adaptor_vertex) {
		if (adaptor_vertex == NULL)
			return false;
		// outer halfedges are not exported
		if ( ! is_inner_or_inner_border(adaptor_vertex)	)				
			return false;
		// only "main" halfedges are exported
		if (get_adaptor_vertex(adaptor_vertex) != adaptor_vertex)					
			return false;
		// eventually: ok
		return true;	
	}

// Fields
private:
	Polyhedron_ex*	m_mesh;							// The adapted mesh 
	std::list<Halfedge_handle> m_boundary;			// Inner halfedges of the boundary of a topological disc inside m_mesh
	int m_mesh_vertices_count;						// Number of adaptor vertices of the mesh


// Private types
private:
	// Utility class to generate the Face_iterator type
	struct All_faces_filter {
		// Return TRUE <=> the face IS NOT EXPORTED by Mesh_adaptor_polyhedron_ex, ie is out of the topological disc
		bool operator()(const Polyhedron_ex::Facet_iterator& f) const		{ return seaming(f) == OUTER; }
		bool operator()(const Polyhedron_ex::Facet_const_iterator& f) const	{ return seaming(f) == OUTER; }
	};

	// Utility class to generate the Vertex_iterator type
	struct All_vertices_filter {
		// Return TRUE <=> the halfedge IS NOT EXPORTED by Mesh_adaptor_polyhedron_ex
		//             <=> the halfedge is outer or is inner but not its vertex's "main" halfedge
		// Implementation note: we could save time by tagging halfedges to skip as INNER_IGNORED
		bool operator()(const Halfedge_iterator& h) const		{ return get_adaptor_vertex(h) != h; }
		bool operator()(const Halfedge_const_iterator& h) const { return get_adaptor_vertex(h) != h; }
	};

	// Utility class to generate the Vertex_around_face_circulator type
	struct Project_halfedge_vertex {
		typedef Halfedge							argument_type;
		typedef Mesh_adaptor_polyhedron_ex::Vertex	Vertex;
		typedef Vertex								result_type;
		typedef CGAL::Arity_tag<1>					Arity;

		// Get the adaptor vertex of the halfedge 'h'
		Vertex&       operator()(Halfedge& h)       const { return *(get_safe_adaptor_vertex((Halfedge_handle)&h)); }
		const Vertex& operator()(const Halfedge& h) const { return *(get_safe_adaptor_vertex(&h)); }
	};

	// Utility class to generate the Border_vertex_iterator type
	struct Project_halfedge_handle_vertex {
		typedef Halfedge_handle						argument_type;
		typedef Mesh_adaptor_polyhedron_ex::Vertex	Vertex;
		typedef Vertex								result_type;
		typedef CGAL::Arity_tag<1>					Arity;

		// Get the adaptor vertex of the halfedge handle 'h'
		Vertex&       operator()(Halfedge_handle& h)       const { return *(get_safe_adaptor_vertex(h)); }
		const Vertex& operator()(const Halfedge_handle& h) const { return *(get_safe_adaptor_vertex(h)); }
	};

	// Circulator over the vertices incident to a vertex (clockwise)
	// Utility class to generate the Vertex_around_vertex_circulator type
	template <class It,								// Mesh_adaptor_polyhedron_ex::Vertex_handle or Vertex_const_handle
			class Ctg>								// Iterator category
	class Adaptor_vertex_around_vertex_circulator : public It 
	{
	public:
		typedef  It													Iterator;
		typedef  Ctg												iterator_category;
		typedef  Adaptor_vertex_around_vertex_circulator<It,Ctg>	Self;
		typedef  std::iterator_traits<It>							Traits;
		typedef  typename Traits::value_type						value_type;
		typedef  typename Traits::difference_type					difference_type;
		typedef  std::size_t										size_type;
		typedef  typename Traits::reference							reference;
		typedef  typename Traits::pointer							pointer;

	// CREATION
	// --------

		// Circulator pointing to NULL
		Adaptor_vertex_around_vertex_circulator() {}
		// Constructor from an out vertex
		explicit Adaptor_vertex_around_vertex_circulator(It out_halfedge, It center_vertex) 
			: m_current_halfedge(out_halfedge), 
			m_center_vertex(center_vertex) 
		{
			assert(is_valid(m_center_vertex));

			// If m_center_vertex is a border vertex, then [ccw_inner_border, m_cw_inner_border]
			// will be the range of the inner halfedges in or out of m_center_vertex
			if (seaming(m_center_vertex->vertex()) == BORDER)
			{
				// do a counterclockwise search around the vertex until we reach the inner boundary 
				// halfedge incident to m_center_vertex
				m_ccw_inner_border = center_vertex;
				while (m_ccw_inner_border->opposite()->seaming() != BORDER)
					m_ccw_inner_border = m_ccw_inner_border->opposite()->prev();

				// do a clockwise search around the vertex until we reach the inner boundary 
				// halfedge out of m_center_vertex
				m_cw_inner_border = center_vertex->next();
				while (m_cw_inner_border->opposite()->seaming() != BORDER)
					m_cw_inner_border = m_cw_inner_border->opposite()->next();
			}

			// m_ccw_inner_border is incident to m_center_vertex. The other values of m_current_halfedge are out halfedges.
			assert(m_current_halfedge == m_ccw_inner_border || m_current_halfedge->opposite()->vertex() == center_vertex->vertex());

			// Update the inherited neighbor adaptor vertex pointer 
			update_adaptor_vertex();
		}
		// Copy constructor
		template <class It2, class Ctg2>
		Adaptor_vertex_around_vertex_circulator(const Adaptor_vertex_around_vertex_circulator<It2,Ctg2> &c)
			: It((const It2&)(c)), 
			m_current_halfedge(c.m_current_halfedge), 
			m_center_vertex(c.m_center_vertex)
		{}

		// Default operator=() is fine

	// OPERATIONS Forward Category
	// ---------------------------

		bool operator==(CGAL_NULL_TYPE p) const {
			CGAL_assertion(p == 0);
			return It::operator==(It());
		}
		bool operator!=(CGAL_NULL_TYPE p) const { return !(*this == p); }
		bool operator==(const Self& i)    const { return  It::operator==(i); }
		bool operator!=(const Self& i)    const { return !(*this == i); }

		//	operator*() and operator->() are inherited

		// Clockwise rotation
		Self& operator++() 
		{
			// Increment m_current_halfedge
			//
			// If non border vertex, m_current_halfedge simply circulates over out halfedges
			if (m_cw_inner_border == NULL)
			{
				m_current_halfedge = m_current_halfedge->opposite()->next();
			}
			else // if border vertex, circulates only from m_ccw_inner_border to m_cw_inner_border (included)
			{	 // CAUTION: m_ccw_inner_border is incident to m_center_vertex. The others are out halfedges.
				if (m_current_halfedge == m_cw_inner_border)
					m_current_halfedge = m_ccw_inner_border;
				else if (m_current_halfedge == m_ccw_inner_border)
					m_current_halfedge = m_current_halfedge->next();
				else
					m_current_halfedge = m_current_halfedge->opposite()->next();
			}

			// Update the inherited neighbor adaptor vertex pointer
			update_adaptor_vertex();

			return *this;
		}
		Self  operator++(int) {
			Self tmp = *this;
			++*this;
			return tmp;
		}

	// OPERATIONS Bidirectional Category
	// ---------------------------------

		// Counterclockwise rotation
		Self& operator--() 
		{
			// Decrement m_current_halfedge
			//
			// If non border vertex, m_current_halfedge simply circulates over out halfedges
			if (m_cw_inner_border == NULL)
			{
				m_current_halfedge = m_current_halfedge->prev()->opposite();
			}
			else // if border vertex, circulates only from m_cw_inner_border to m_ccw_inner_border (included)
			{	 // CAUTION: m_ccw_inner_border is incident to m_center_vertex. The others are out halfedges.
				if (m_current_halfedge == m_ccw_inner_border) {
					m_current_halfedge = m_cw_inner_border;
				} else {
					m_current_halfedge = m_current_halfedge->prev()->opposite();
					if (m_current_halfedge == m_ccw_inner_border->opposite())
						m_current_halfedge = m_ccw_inner_border;
				}
			}

			// Update the inherited neighbor adaptor vertex pointer
			update_adaptor_vertex();

			return *this;
		}
		Self  operator--(int) {
			Self tmp = *this;
			--*this;
			return tmp;
		}

	private:
		// Update the inherited neighbor adaptor vertex pointer
		// Precondition: m_current_halfedge and m_center_vertex are valid
		void update_adaptor_vertex()
		{
			// Get halfedge incident to current neighbor vertex
			It	neighbor_halfedge = (m_current_halfedge->vertex() != m_center_vertex->vertex()) ? m_current_halfedge
																								: m_current_halfedge->prev();
			// The neighbor adaptor vertex is:
			It::operator=(get_safe_adaptor_vertex(neighbor_halfedge));
		}

	private:
		It	m_center_vertex;		// Polyhedron vertex center of the circulation
		It	m_current_halfedge;		// Current position of the circulator among m_center_vertex in and out halfedges
		It	m_ccw_inner_border,		// For border vertices: - border halfedge incident to m_center_vertex when searching counterclockwise
			m_cw_inner_border;		//                      - border halfedge out of m_center_vertex when searching counterclockwise

	}; // Adaptor_vertex_around_vertex_circulator

}; // Mesh_adaptor_polyhedron_ex


#endif //MESHADAPTORPOLYHEDRONEX_H

