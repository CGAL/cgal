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
// Author(s)     : Laurent Saboret, Pierre Alliez


#ifndef CGAL_MESHADAPTORPOLYHEDRONEX_H
#define CGAL_MESHADAPTORPOLYHEDRONEX_H

#include <CGAL/iterator.h>
#include <CGAL/circulator.h>
#include <CGAL/Convertible_iterator_project.h>
#include <CGAL/Convertible_circulator_project.h>

#include "cgal_types.h"
#include "Mesh_feature_extractor.h"

#include <cassert>


// Class Mesh_adaptor_polyhedron_ex
// Model of MeshAdaptor_3 concept, whose purpose is to allow
// the parameterization package to access meshes on an uniform manner.
//
// Mesh_adaptor_polyhedron_ex is an adaptor class to access to a Polyhedron_ex
// 3D mesh using MeshAdaptor_3 interface.
//
// The input mesh can be of any genus.
// It can have have any number of boundaries. Its "main border"
// will be the mesh's longest boundary (if there is at least one boundary). 
//
// Design pattern:
// Mesh_adaptor_polyhedron_ex is an Adaptor (see [GOF95]): it changes the
// Polyhedron_ex interface to match the MeshAdaptor_3 concept

class Mesh_adaptor_polyhedron_ex
{
// Private types
private:

    // Forward references
    struct                                  Project_halfedge_vertex;
    struct                                  Project_vertex_handle_vertex;
    struct                                  Project_opposite_halfedge_vertex;
    // Halfedge
    typedef Polyhedron_ex::Halfedge         Halfedge;
    typedef Polyhedron_ex::Halfedge_handle  Halfedge_handle;
    typedef Polyhedron_ex::Halfedge_const_handle
                                            Halfedge_const_handle;
    typedef Polyhedron_ex::Halfedge_iterator Halfedge_iterator;
    typedef Polyhedron_ex::Halfedge_const_iterator
                                            Halfedge_const_iterator;
    typedef Polyhedron_ex::Halfedge_around_vertex_circulator
                                            Halfedge_around_vertex_circulator;
    typedef Polyhedron_ex::Halfedge_around_vertex_const_circulator
                                            Halfedge_around_vertex_const_circulator;
    typedef Polyhedron_ex::Halfedge_around_facet_circulator
                                            Halfedge_around_facet_circulator;
    typedef Polyhedron_ex::Halfedge_around_facet_const_circulator
                                            Halfedge_around_facet_const_circulator;

public:

    //******************************************************************
    // Public types
    //******************************************************************

    // Number type
    typedef Polyhedron_ex::Traits::FT       NT;
    // Points and vectors
    typedef Polyhedron_ex::Traits::Point_2  Point_2;
    typedef Polyhedron_ex::Traits::Point_3  Point_3;
    typedef Polyhedron_ex::Traits::Vector_2 Vector_2;
    typedef Polyhedron_ex::Traits::Vector_3 Vector_3;
    // Facet
    typedef Polyhedron_ex::Facet            Facet;
    typedef Polyhedron_ex::Facet_handle     Facet_handle;
    typedef Polyhedron_ex::Facet_const_handle Facet_const_handle;
    // Iterator over all mesh facets
    typedef Polyhedron_ex::Facet_iterator   Facet_iterator;
    typedef Polyhedron_ex::Facet_const_iterator 
                                            Facet_const_iterator;
    // Vertex
    typedef Polyhedron_ex::Vertex           Vertex;
    typedef Polyhedron_ex::Vertex_handle    Vertex_handle;
    typedef Polyhedron_ex::Vertex_const_handle Vertex_const_handle;
    // Iterator over all mesh vertices
    typedef Polyhedron_ex::Vertex_iterator  Vertex_iterator;
    typedef Polyhedron_ex::Vertex_const_iterator 
                                            Vertex_const_iterator;
    // Iterator over mesh boundary vertices
    typedef CGAL::Convertible_iterator_project<std::list<Vertex_handle>::iterator,
                                               Project_vertex_handle_vertex,
                                               Vertex_const_handle, 
                                               Vertex_handle>
                                            Border_vertex_iterator;
    typedef CGAL::Convertible_iterator_project<std::list<Vertex_handle>::const_iterator,
                                               Project_vertex_handle_vertex,
                                               Vertex_const_handle>
                                            Border_vertex_const_iterator;
    // Counter-clockwise circulator over a facet's vertices
    typedef CGAL::Convertible_circulator_project<Halfedge_around_facet_circulator,
                                                 Project_halfedge_vertex,
                                                 Vertex&, 
                                                 Vertex*,
                                                 Vertex_const_handle, 
                                                 Vertex_handle>
                                            Vertex_around_facet_circulator;
    typedef CGAL::Convertible_circulator_project<Halfedge_around_facet_const_circulator,
                                                 Project_halfedge_vertex,
                                                 const Vertex&, 
                                                 const Vertex*,
                                                 Vertex_const_handle>
                                            Vertex_around_facet_const_circulator;
    // Clockwise circulator over the vertices incident to a vertex
    typedef CGAL::Convertible_circulator_project<Halfedge_around_vertex_circulator,
                                                 Project_opposite_halfedge_vertex,
                                                 Vertex&, 
                                                 Vertex*,
                                                 Vertex_const_handle, 
                                                 Vertex_handle>
                                            Vertex_around_vertex_circulator;
    typedef CGAL::Convertible_circulator_project<Halfedge_around_vertex_const_circulator,
                                                 Project_opposite_halfedge_vertex,
                                                 const Vertex&, 
                                                 const Vertex*,
                                                 Vertex_const_handle>
                                            Vertex_around_vertex_const_circulator;
    // Seaming flag
    enum Seaming_status  { OUTER, INNER, BORDER };

// Public operations
public:

    //******************************************************************
    // LIFE CYCLE
    //******************************************************************

    // Create an adaptator for an existing Polyhedron_ex mesh.
    // The input mesh can be of any genus.
    // It can have have any number of boundaries. Its "main border"
    // will be the mesh's longest boundary (if there is at least one boundary). 
    Mesh_adaptor_polyhedron_ex(Polyhedron_ex* mesh)
    {
        assert(mesh != NULL);

        // Extract mesh's longest boundary
        std::list<Halfedge_handle> boundary = extract_longest_boundary(mesh);

        // Store adapted mesh pointer
        m_polyhedron = mesh;

        // Copy main boundary
        for (std::list<Halfedge_handle>::iterator border_it = boundary.begin();
             border_it != boundary.end();
             border_it++)
        {
            m_main_border.push_back((*border_it)->vertex());
        }

        // Initialize the seaming flag of all vertices and halfedges to INNER or BORDER
        //
        // 1) Initialize the seaming flag of all vertices to INNER
        set_vertices_seaming(INNER);
        //
        // 2) Initialize the seaming flag of all halfedges to INNER
        set_edges_seaming(INNER);
        //
        // 3) Set seaming flag of main border vertices and halfedges to BORDER
        for (std::list<Vertex_handle>::iterator border_it = m_main_border.begin();
             border_it != m_main_border.end();
             border_it++)
        {
            // Set vertex seaming flag
            set_vertex_seaming(*border_it, BORDER);

            // Get next iterator (looping)
            std::list<Vertex_handle>::iterator next_border_it = border_it; 
            next_border_it++;
            if (next_border_it == m_main_border.end())
                next_border_it = m_main_border.begin();

            // Set halfedge seaming flag
            set_edge_seaming(*border_it, *next_border_it, BORDER);
        }

#ifndef NDEBUG
        // Index vertices right away to ease debugging
        index_mesh_vertices();
#endif
    }

    // Default copy constructor and operator =() are fine

    //******************************************************************
    // LEVEL 1 INTERFACE:
    // for classes attempting to parameterize complete topological disks
    // and compute 1 (u,v) pair per vertex
    // Example: all parameterization methods
    //******************************************************************

    // MESH INTERFACE

    // Get iterator over first vertex of mesh
    Vertex_iterator  mesh_vertices_begin() {
        return m_polyhedron->vertices_begin();
    }
    Vertex_const_iterator  mesh_vertices_begin() const {
        return m_polyhedron->vertices_begin();
    }

    // Get iterator over past-the-end vertex of mesh
    Vertex_iterator  mesh_vertices_end() {
        return m_polyhedron->vertices_end();
    }
    Vertex_const_iterator  mesh_vertices_end() const {
        return m_polyhedron->vertices_end();
    }

    // Count the number of vertices of the mesh
    int  count_mesh_vertices() const {
        int index = 0;
        for (Vertex_const_iterator it=mesh_vertices_begin(); it!=mesh_vertices_end(); it++)
            index++;
        return index;
    }

    // Index vertices of the mesh from 0 to count_mesh_vertices()-1
    void  index_mesh_vertices ()
    {
        fprintf(stderr,"  index Mesh_adaptor_polyhedron_ex vertices:\n");
        int index = 0;
        for (Vertex_iterator it=mesh_vertices_begin(); it!=mesh_vertices_end(); it++)
        {
            Point_3 position = get_vertex_position(it);
#ifdef DEBUG_TRACE
            fprintf(stderr, "    %d=(%f,%f,%f)\n",
                            index,
                            (float)position.x(), 
                            (float)position.y(), 
                            (float)position.z());
#endif
            set_vertex_index(it, index++);
        }
        fprintf(stderr,"    ok\n");
    }

    // Return true of all mesh's facets are triangles
    bool  is_mesh_triangular() const {
        for (Facet_const_iterator it = mesh_facets_begin(); it != mesh_facets_end(); it++)
            if (count_facet_vertices(it) != 3)
                return false;
        return true;            // mesh is triangular if we reach this point
    }

    // Compute the genus of the mesh
    int  get_mesh_genus() const {
        return m_polyhedron->genus();
    }

    // Count the number of boundaries of the mesh
    int  count_mesh_boundaries() const {
        return m_polyhedron->nb_boundaries();
    }

    // Get iterator over first vertex of mesh's main border
    Border_vertex_iterator  mesh_main_border_vertices_begin() {
        return m_main_border.begin();
    }
    Border_vertex_const_iterator  mesh_main_border_vertices_begin() const {
        return m_main_border.begin();
    }

    // Get iterator over past-the-end vertex of mesh's main border
    Border_vertex_iterator  mesh_main_border_vertices_end() {
        return m_main_border.end();
    }
    Border_vertex_const_iterator  mesh_main_border_vertices_end() const {
        return m_main_border.end();
    }

    // Get iterator over first facet of mesh
    Facet_iterator  mesh_facets_begin() {
        return m_polyhedron->facets_begin();
    }
    Facet_const_iterator  mesh_facets_begin() const {
        return m_polyhedron->facets_begin();
    }

    // Get iterator over past-the-end facet of mesh
    Facet_iterator  mesh_facets_end() {
        return m_polyhedron->facets_end();
    }
    Facet_const_iterator  mesh_facets_end() const {
        return m_polyhedron->facets_end();
    }

    // Count the number of facets of the mesh
    int  count_mesh_facets() const {
        int index = 0;
        for (Facet_const_iterator it=mesh_facets_begin(); it!=mesh_facets_end(); it++)
            index++;
        return index;
    }

    // FACET INTERFACE

    // Get circulator over facet's vertices
    Vertex_around_facet_circulator  facet_vertices_begin(Facet_handle facet) {
        return facet->facet_begin();
    }
    Vertex_around_facet_const_circulator  facet_vertices_begin(Facet_const_handle facet) const {
        return facet->facet_begin();
    }

    // Count the number of vertices of a facet
    int  count_facet_vertices(Facet_const_handle facet) const {
        int index = 0;
        Vertex_around_facet_const_circulator cir     = facet_vertices_begin(facet),
                                             cir_end = cir;
        CGAL_For_all(cir, cir_end)
            index++;
        return index;
    }

    // VERTEX INTERFACE

    // Get the 3D position of a vertex
    Point_3 get_vertex_position(Vertex_const_handle vertex) const {
        return vertex->point();    
    }

    // Get/set the 2D position (u/v pair) of a vertex
    // (stored in halfedges sharing the same vertex)
    Point_2  get_vertex_uv(Vertex_const_handle vertex) const {
        return get_corners_uv(vertex, NULL, NULL);
    }
    void  set_vertex_uv(Vertex_handle vertex, const Point_2& uv) {
        set_corners_uv(vertex, NULL, NULL, uv);
    }

    // Get/set "is parameterized" field of vertex
    // (stored in halfedges sharing the same vertex)
    bool  is_vertex_parameterized(Vertex_const_handle vertex) const {
        return are_corners_parameterized(vertex, NULL, NULL);
    }
    void  set_vertex_parameterized(Vertex_handle vertex, bool parameterized) {
        set_corners_parameterized(vertex, NULL, NULL, parameterized);
    }

    // Get/set vertex index
    // (stored in Polyhedron_ex vertex for debugging purpose)
    int  get_vertex_index(Vertex_const_handle vertex) const {
        //return get_corners_index(vertex, NULL, NULL);
        return vertex->index();
    }
    void  set_vertex_index(Vertex_handle vertex, int index)  
    {
        //set_corners_index(vertex, NULL, NULL, index);
        vertex->index(index);
    }

    //// Return true if a vertex belongs to ANY mesh's boundary
    //bool  is_vertex_on_border(Vertex_const_handle vertex) const {
    //    return m_polyhedron->is_border(vertex);
    //}

    // Return true if a vertex belongs to the UNIQUE mesh's main boundary
    bool  is_vertex_on_main_border(Vertex_const_handle vertex) const {
        return get_vertex_seaming(vertex) == BORDER;
    }

    // Get circulator over the vertices incident to 'vertex'
    // 'start_position' defines the optional initial position of the circulator
    Vertex_around_vertex_circulator vertices_around_vertex_begin(
                            Vertex_handle vertex, 
                            Vertex_handle start_position = NULL) 
    {
        if (start_position == NULL) 
            return vertex->vertex_begin();
        else
            return Halfedge_around_vertex_circulator(get_halfedge(start_position, 
                                                                  vertex));
    }
    Vertex_around_vertex_const_circulator vertices_around_vertex_begin(
                            Vertex_const_handle vertex, 
                            Vertex_const_handle start_position = NULL) const
    {
        if (start_position == NULL) 
            return vertex->vertex_begin();
        else
            return Halfedge_around_vertex_const_circulator(get_halfedge(start_position, 
                                                                        vertex));
    }

    //******************************************************************
    // LEVEL 2 INTERFACE:
    // for classes attempting to parameterize (part of) 3D surfaces
    // of any genus and with any number of connected components.
    // They compute 1 (u,v) pair per corner.
    // Example: Mesh_adaptor_patch_3
    //******************************************************************

    // VERTEX INTERFACE

    // Get/set vertex seaming flag, ie position of the vertex 
    // wrt to the UNIQUE main boundary
    Seaming_status  get_vertex_seaming(Vertex_const_handle vertex) const {
        return (Seaming_status) vertex->seaming();
    }
    void set_vertex_seaming(Vertex_handle vertex, Seaming_status seaming) {
        return vertex->seaming(seaming);
    }

    // Set all vertices seaming flag, ie position of the vertices
    // wrt to the UNIQUE main boundary
    void set_vertices_seaming(Seaming_status seaming)
    {
        for (Polyhedron_ex::Vertex_iterator v = m_polyhedron->vertices_begin();
             v != m_polyhedron->vertices_end();
             v++)
        {
             v->seaming(seaming);
        }
    }

    // EDGE INTERFACE

    // Get/set oriented edge's seaming flag, ie position of the oriented edge 
    // wrt to the UNIQUE main boundary
    Seaming_status  get_edge_seaming(Vertex_const_handle source, Vertex_const_handle target) const {
        return (Seaming_status) get_halfedge(source, target)->seaming();
    }
    void set_edge_seaming(Vertex_handle source, Vertex_handle target, Seaming_status seaming) {
        return get_halfedge(source, target)->seaming(seaming);
    }

    // Set all edges seaming flag, ie position of the edges
    // wrt to the UNIQUE main boundary
    void set_edges_seaming(Seaming_status seaming)
    {
        for (Polyhedron_ex::Halfedge_iterator h = m_polyhedron->halfedges_begin();
             h != m_polyhedron->halfedges_end();
             h++)
        {
             h->seaming(seaming);
        }
    }

    // CORNER INTERFACE

    // Get/set the 2D position (= (u,v) pair) of corners at the "right" 
    // of the prev_vertex -> vertex -> next_vertex line
    // (stored in incident halfedges)
    Point_2 get_corners_uv(Vertex_const_handle vertex, 
                           Vertex_const_handle prev_vertex,
                           Vertex_const_handle next_vertex) const 
    {
        // if inner vertex
        if (prev_vertex == NULL && next_vertex == NULL)
        {
            // get (u,v) pair from any incident halfedge
            return Point_2(vertex->halfedge()->u(), vertex->halfedge()->v());
        }
        else // if seam vertex
        {
            assert(prev_vertex != NULL);
            assert(next_vertex != NULL);

            // get (u,v) pair from first inner halfedge (clockwise)
            Halfedge_around_vertex_const_circulator cir( 
                                get_halfedge(next_vertex, vertex) );
            return Point_2(cir->u(), cir->v());
        }
    }
    void set_corners_uv(Vertex_handle vertex, 
                        Vertex_const_handle prev_vertex,
                        Vertex_const_handle next_vertex,
                        const Point_2& uv)
    {
        // if inner vertex
        if (prev_vertex == NULL && next_vertex == NULL)
        {
            // Loop over all incident halfedges
            Halfedge_around_vertex_circulator cir     = vertex->vertex_begin(),
                                              cir_end = cir;
            CGAL_For_all(cir, cir_end)
                cir->uv(uv.x(), uv.y());
        }
        else // if seam vertex
        {
            assert(prev_vertex != NULL);
            assert(next_vertex != NULL);

            // first inner halfedge (for a clockwise rotation)
            Halfedge_around_vertex_circulator cir( 
                                get_halfedge((Vertex*)&*next_vertex, vertex) );

            // past-the-end inner halfedge (for a clockwise rotation)
            Halfedge_around_vertex_circulator cir_end( 
                                get_halfedge((Vertex*)&*prev_vertex, vertex) );

            // Loop over incident halfedges at the "right" 
            // of the prev_vertex -> vertex -> next_vertex line
            CGAL_For_all(cir, cir_end)
                cir->uv(uv.x(), uv.y());
        }
    }

    // Get/set "is parameterized" field of corners at the "right" 
    // of the prev_vertex -> vertex -> next_vertex line
    // (stored in incident halfedges)
    bool are_corners_parameterized(Vertex_const_handle vertex, 
                                   Vertex_const_handle prev_vertex,
                                   Vertex_const_handle next_vertex) const 
    {
        // if inner vertex
        if (prev_vertex == NULL && next_vertex == NULL)
        {
            // get "is parameterized" field from any incident halfedge
            return vertex->halfedge()->is_parameterized();
        }
        else // if seam vertex
        {
            assert(prev_vertex != NULL);
            assert(next_vertex != NULL);

            // get "is parameterized" field from first inner halfedge (clockwise)
            Halfedge_around_vertex_const_circulator cir( 
                                get_halfedge(next_vertex, vertex) );
            return cir->is_parameterized();
        }
    }
    void set_corners_parameterized(Vertex_handle vertex, 
                                   Vertex_const_handle prev_vertex,
                                   Vertex_const_handle next_vertex,
                                   bool parameterized)
    {
        // if inner vertex
        if (prev_vertex == NULL && next_vertex == NULL)
        {
            // Loop over all incident halfedges
            Halfedge_around_vertex_circulator cir     = vertex->vertex_begin(),
                                              cir_end = cir;
            CGAL_For_all(cir, cir_end)
                cir->is_parameterized(parameterized);
        }
        else // if seam vertex
        {
            assert(prev_vertex != NULL);
            assert(next_vertex != NULL);

            // first inner halfedge (for a clockwise rotation)
            Halfedge_around_vertex_circulator cir( 
                                get_halfedge((Vertex*)&*next_vertex, vertex) );

            // past-the-end inner halfedge (for a clockwise rotation)
            Halfedge_around_vertex_circulator cir_end( 
                                get_halfedge((Vertex*)&*prev_vertex, vertex) );

            // Loop over incident halfedges at the "right" 
            // of the prev_vertex -> vertex -> next_vertex line
            CGAL_For_all(cir, cir_end)
                cir->is_parameterized(parameterized);
        }
    }

    // Get/set index of corners at the "right" 
    // of the prev_vertex -> vertex -> next_vertex line
    // (stored in incident halfedges)
    int get_corners_index(Vertex_const_handle vertex, 
                          Vertex_const_handle prev_vertex,
                          Vertex_const_handle next_vertex) const 
    {
        // if inner vertex
        if (prev_vertex == NULL && next_vertex == NULL)
        {
            // get index from any incident halfedge
            return vertex->halfedge()->index();
        }
        else // if seam vertex
        {
            assert(prev_vertex != NULL);
            assert(next_vertex != NULL);

            // get index from first inner halfedge (clockwise)
            Halfedge_around_vertex_const_circulator cir( 
                                get_halfedge(next_vertex, vertex) );
            return cir->index();
        }
    }
    void set_corners_index(Vertex_handle vertex, 
                           Vertex_const_handle prev_vertex,
                           Vertex_const_handle next_vertex,
                           int index)  
    {
        // if inner vertex
        if (prev_vertex == NULL && next_vertex == NULL)
        {
            // Loop over all incident halfedges
            Halfedge_around_vertex_circulator cir     = vertex->vertex_begin(),
                                              cir_end = cir;
            CGAL_For_all(cir, cir_end)
                cir->index(index);
        }
        else // if seam vertex
        {
            assert(prev_vertex != NULL);
            assert(next_vertex != NULL);

            // first inner halfedge (for a clockwise rotation)
            Halfedge_around_vertex_circulator cir( 
                                get_halfedge((Vertex*)&*next_vertex, vertex) );

            // past-the-end inner halfedge (for a clockwise rotation)
            Halfedge_around_vertex_circulator cir_end( 
                                get_halfedge((Vertex*)&*prev_vertex, vertex) );

            // Loop over incident halfedges at the "right" 
            // of the prev_vertex -> vertex -> next_vertex line
            CGAL_For_all(cir, cir_end)
                cir->index(index);
        }
    }

// Private operations
private:

    // Extract mesh's longest boundary
    static std::list<Halfedge_handle> extract_longest_boundary(Polyhedron_ex* mesh)
    {
        typedef Feature_backbone<Vertex_handle, Halfedge_handle> Backbone;

        std::list<Halfedge_handle> boundary;    // returned list

        assert(mesh != NULL);

        // Extract mesh boundaries
        mesh->free_skeleton();
        Mesh_feature_extractor feature_extractor(mesh);
        int nb_boundaries = feature_extractor.extract_boundaries(true);

        // Get longest one
        if (nb_boundaries > 0) {
            Backbone *pBackbone = (*mesh->get_skeleton()->backbones())[0];
            boundary = *(pBackbone->halfedges());
        }

        // Cleanup
        mesh->free_skeleton();

        return boundary;
    }

    // Get halfedge from source and target vertices
    // Will assert if such an halfedge doesn't exist
    Halfedge_const_handle get_halfedge(Vertex_const_handle source, 
                                       Vertex_const_handle target) const
    {
        assert(source != NULL);
        assert(target != NULL);

        Halfedge_around_vertex_const_circulator cir     = target->vertex_begin(),
                                                cir_end = cir;
        CGAL_For_all(cir, cir_end)
            if (cir->opposite()->vertex() == source)
                return cir;

        assert(false);              // error if we reach this point
        return NULL;
    }
    inline Halfedge_handle get_halfedge(Vertex_handle source, 
                                        Vertex_handle target)
    {
        Halfedge_const_handle halfedge =
            get_halfedge((Vertex_const_handle)source, (Vertex_const_handle)target);
        return (Halfedge*) (&*halfedge);
    }

// Fields
private:

    // The adapted mesh (cannot be NULL)
    Polyhedron_ex*              m_polyhedron;

    // Main boundary of a topological disc inside m_polyhedron (may be empty)
    std::list<Vertex_handle>    m_main_border;

// Private types
private:

    // Utility class to generate the Vertex_around_facet_circulator type
    struct Project_halfedge_vertex {
        typedef Halfedge                            argument_type;
        typedef Mesh_adaptor_polyhedron_ex::Vertex  Vertex;
        typedef Vertex                              result_type;
        typedef CGAL::Arity_tag<1>                  Arity;

        // Get the target vertex of a halfedge
        Vertex&       operator()(Halfedge& h)       const {
            return *(h.vertex());
        }
        const Vertex& operator()(const Halfedge& h) const {
            return *(h.vertex());
        }
    };

    // Utility class to generate the Border_vertex_iterator type
    struct Project_vertex_handle_vertex {
        typedef Vertex_handle                       argument_type;
        typedef Mesh_adaptor_polyhedron_ex::Vertex  Vertex;
        typedef Vertex                              result_type;
        typedef CGAL::Arity_tag<1>                  Arity;

        // Convert Vertex_handle to Vertex
        Vertex&       operator()(Vertex_handle& vh)       const { return *vh; }
        const Vertex& operator()(const Vertex_handle& vh) const { return *vh; }
    };

    // This class is used to generate the Vertex_around_vertex_circulator type
    struct Project_opposite_halfedge_vertex {
        typedef Halfedge                            argument_type;
        typedef Mesh_adaptor_polyhedron_ex::Vertex  Vertex;
        typedef Vertex                              result_type;
        typedef CGAL::Arity_tag<1>                  Arity;

        // Get the source vertex of a halfedge
        Vertex&       operator()(Halfedge& h)       const {
            return *(h.opposite()->vertex()); 
        }
        const Vertex& operator()(const Halfedge& h) const {
            return *(h.opposite()->vertex()); 
        }
    };

}; // Mesh_adaptor_polyhedron_ex


#endif //CGAL_MESHADAPTORPOLYHEDRONEX_H

