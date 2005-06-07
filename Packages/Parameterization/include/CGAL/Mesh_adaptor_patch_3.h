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


#ifndef CGAL_MESH_ADAPTOR_PATCH_3_H
#define CGAL_MESH_ADAPTOR_PATCH_3_H

#include <CGAL/iterator.h>
#include <CGAL/circulator.h>
#include <CGAL/Convertible_iterator_project.h>
#include <CGAL/Convertible_circulator_project.h>
#include <CGAL/Convertible_filter_iterator.h>
#include <CGAL/Mesh_adaptor_patch_vertex.h>
#include <CGAL/Mesh_adaptor_patch_iterators.h>
#include <CGAL/Mesh_adaptor_patch_circulators.h>

#include <CGAL/parameterization_assertions.h>

CGAL_BEGIN_NAMESPACE


// Class Mesh_adaptor_patch_3
// Model of MeshAdaptor_3 concept, whose purpose is to allow
// the parameterization package to access meshes on an uniform manner.
//
// Mesh_adaptor_patch_3 is an decorator class to virtually "cut" a patch
// in a MeshAdaptor_3 3D surface. Only the patch is exported,
// making the 3D surface look like a topological disk.
//
// The input mesh can be of any genus, but it has to come with a "seam" that
// describes the boundary of a topological disc. This boundary may be an actual
// border of the mesh or a virtual border.
//
// Design pattern:
// Mesh_adaptor_patch_3 is an Decorator (see [GOF95]): it changes the
// behavior of a MeshAdaptor_3 3D surface without changing its interface.

template<class MeshAdaptor_3>
class Mesh_adaptor_patch_3
{
// Private types
private:

    // Forward references
    struct                                  Inner_facets_filter;

public:

    //******************************************************************
    // Public types
    //******************************************************************

    // Export template parameter
    typedef MeshAdaptor_3                   Adaptor;
    // Number type
    typedef typename Adaptor::NT            NT;
    // Points and vectors
    typedef typename Adaptor::Point_2       Point_2;
    typedef typename Adaptor::Point_3       Point_3;
    typedef typename Adaptor::Vector_2      Vector_2;
    typedef typename Adaptor::Vector_3      Vector_3;
    // Facet
    typedef typename Adaptor::Facet         Facet;
    typedef typename Adaptor::Facet_handle  Facet_handle;
    typedef typename Adaptor::Facet_const_handle
                                            Facet_const_handle;
    // Iterator over all mesh facets
    typedef Convertible_filter_iterator<typename Adaptor::Facet_iterator,
                                        Inner_facets_filter,
                                        Facet_const_handle,
                                        Facet_handle>
                                            Facet_iterator;
    typedef Convertible_filter_iterator<typename Adaptor::Facet_const_iterator,
                                        Inner_facets_filter,
                                        Facet_const_handle>
                                            Facet_const_iterator;
    // Vertex
    typedef Mesh_adaptor_patch_vertex<Adaptor> Vertex;
    typedef Mesh_adaptor_patch_vertex_handle<Adaptor>
                                            Vertex_handle;
    typedef Mesh_adaptor_patch_vertex_const_handle<Adaptor>
                                            Vertex_const_handle;
    // Iterator over all mesh vertices
    typedef Mesh_adaptor_patch_vertex_list_iterator<Adaptor>
                                            Vertex_iterator;
    typedef Mesh_adaptor_patch_vertex_list_const_iterator<Adaptor>
                                            Vertex_const_iterator;
    // Iterator over mesh boundary vertices
    typedef Vertex_iterator                 Border_vertex_iterator;
    typedef Vertex_const_iterator           Border_vertex_const_iterator;
    // Counter-clockwise circulator over a facet's vertices
    typedef Mesh_patch_vertex_around_facet_cir<Mesh_adaptor_patch_3*,
                                               Vertex_handle,
                                               typename Adaptor::Vertex_around_facet_circulator>
                                            Vertex_around_facet_circulator;
    typedef Mesh_patch_vertex_around_facet_cir<const Mesh_adaptor_patch_3*,
                                               Vertex_const_handle,
                                               typename Adaptor::Vertex_around_facet_const_circulator>
                                            Vertex_around_facet_const_circulator;
    // Clockwise circulator over the vertices incident to a vertex
    typedef Mesh_patch_vertex_around_vertex_cir<Mesh_adaptor_patch_3*,
                                                Vertex_handle,
                                                typename Adaptor::Vertex_around_vertex_circulator>
                                            Vertex_around_vertex_circulator;
    typedef Mesh_patch_vertex_around_vertex_cir<const Mesh_adaptor_patch_3*,
                                                Vertex_const_handle,
                                                typename Adaptor::Vertex_around_vertex_const_circulator>
                                            Vertex_around_vertex_const_circulator;
    // Seaming flag
    enum Seaming_status  { OUTER = Adaptor::OUTER, INNER = Adaptor::INNER, BORDER = Adaptor::BORDER };

    friend class Vertex;
    friend class Vertex_handle;
    friend class Vertex_const_handle;
    friend class Vertex_iterator;
    friend class Vertex_const_iterator;
    friend class Vertex_around_facet_circulator;
    friend class Vertex_around_facet_const_circulator;
    friend class Vertex_around_vertex_circulator;
    friend class Vertex_around_vertex_const_circulator;

// Public operations
public:

    //******************************************************************
    // LIFE CYCLE
    //******************************************************************

    // Create an decorator for an existing MeshAdaptor_3 mesh
    // The input mesh can be of any genus, but it has to come with a "seam" that
    // describes the boundary of a topological disc. This boundary may be an actual
    // border of the mesh or a virtual border.
    //
    // Preconditions:
    // * first_seam_vertex -> end_seam_vertex defines the outer seam,
    //   ie Mesh_adaptor_patch_3 will export the "right" of the seam
    // * the "seam" is given as a container of Adaptor::Vertex_handle elements.
    template<class InputIterator>
    Mesh_adaptor_patch_3(Adaptor* mesh,
                         InputIterator first_seam_vertex,
                         InputIterator end_seam_vertex)
    {
        CGAL_parameterization_assertion(mesh != NULL);
        m_mesh_adaptor = mesh;

#ifdef DEBUG_TRACE
        // Dump input boundary (for debug purpose)
        fprintf(stderr,"  input boundary is: ");
        for (InputIterator it = first_seam_vertex; it != end_seam_vertex; it++)
            fprintf(stderr, "%s ", get_vertex_index_as_string(*it).c_str());
        fprintf(stderr,"ok\n");
#endif

        // Set seaming flag of all vertices and edges to INNER, BORDER or OUTER
        // wrt the first_seam_vertex -> end_seam_vertex boundary
        set_mesh_seaming(first_seam_vertex, end_seam_vertex);

        // Construct the list of all exported vertices, ie INNER and BORDER vertices
        //
        // 1) add inner vertices
        for (typename Adaptor::Vertex_iterator it = mesh->mesh_vertices_begin();
             it != mesh->mesh_vertices_end();
             it++)
        {
            if (m_mesh_adaptor->get_vertex_seaming(it) == Adaptor::INNER)
                m_inner_and_border_vertices.push_back( Vertex(it) );
        }
        // 2) add seam vertices, wrt outer seam/boundary order
        InputIterator border_it      = first_seam_vertex;
        InputIterator prev_border_it = end_seam_vertex; prev_border_it--;;
        InputIterator next_border_it = first_seam_vertex; next_border_it++;
        while (border_it != end_seam_vertex)
        {
            // Get outer border vertex
            Vertex v;
            // if border vertex
            if (m_mesh_adaptor->get_edge_seaming(*border_it, *prev_border_it) != Adaptor::BORDER)
                v = Vertex(*border_it, *prev_border_it, *next_border_it);
            else // if seam vertex
                v = Vertex(*border_it, *next_border_it, *prev_border_it);   // order inverted!

            // Add vertex
            m_inner_and_border_vertices.push_back(v);

            // Increment iterators
            border_it++;
            //
            prev_border_it++;
            if (prev_border_it == end_seam_vertex)
                prev_border_it = first_seam_vertex;
            //
            next_border_it++;
            if (next_border_it == end_seam_vertex)
                next_border_it = first_seam_vertex;
        }

        // Initialize m_seam_begin = iterator to beginning of seam/main boundary
        // inside m_inner_and_border_vertices
        m_seam_begin = mesh_vertices_end();
        for (Vertex_iterator it=mesh_vertices_begin(); it!=mesh_vertices_end(); it++) {
            if (get_vertex_seaming(it) == BORDER) {
                m_seam_begin = it;
                break;
            }
        }

#ifndef NDEBUG
        // Index vertices right away to ease debugging
        index_mesh_vertices();

    #ifdef DEBUG_TRACE
        // Dump seam (for debug purpose)
        fprintf(stderr,"  seam is: ");
        for (Border_vertex_iterator border_it = mesh_main_border_vertices_begin();
             border_it != mesh_main_border_vertices_end();
             border_it++)
        {
            fprintf(stderr, "#%d ", get_vertex_index(border_it));
        }
        fprintf(stderr,"ok\n");
    #endif
#endif
    }

    //******************************************************************
    // LEVEL 1 INTERFACE:
    // for classes attempting to parameterize complete topological disks
    // and compute 1 (u,v) pair per vertex
    // Example: all parameterization methods
    //******************************************************************

    // MESH INTERFACE

    // Get iterator over first vertex of mesh
    Vertex_iterator  mesh_vertices_begin() {
        return m_inner_and_border_vertices.begin();
    }
    Vertex_const_iterator  mesh_vertices_begin() const {
        return m_inner_and_border_vertices.begin();
    }

    // Get iterator over past-the-end vertex of mesh
    Vertex_iterator  mesh_vertices_end() {
        return m_inner_and_border_vertices.end();
    }
    Vertex_const_iterator  mesh_vertices_end() const {
        return m_inner_and_border_vertices.end();
    }

    // Count the number of vertices of the mesh
    int  count_mesh_vertices() const {
        int index = 0;
        for (Vertex_const_iterator it=mesh_vertices_begin(); it!=mesh_vertices_end(); it++)
            index++;
        return index;
    }

    // Index vertices of the mesh from 0 to count_mesh_vertices()-1
    void  index_mesh_vertices()
    {
        fprintf(stderr,"  index Mesh_adaptor_patch_3 vertices:\n");
        int index = 0;
        for (Vertex_iterator it=mesh_vertices_begin(); it!=mesh_vertices_end(); it++)
        {
#ifdef DEBUG_TRACE
            fprintf(stderr, "    #%d = {%s,%s,%s}\n",
                            index,
                            get_vertex_index_as_string(it->get_adaptor_vertex()).c_str(),
                            get_vertex_index_as_string(it->get_prev_seam_vertex()).c_str(),
                            get_vertex_index_as_string(it->get_next_seam_vertex()).c_str());
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
    int  get_mesh_genus() const
    {
        // We assume that cutting created a topological disk
        // TODO: compute genus
        return 0;
    }

    // Count the number of boundaries of the mesh
    int  count_mesh_boundaries() const
    {
        // We have at least 1 boundary (the "main one")
        // TODO: compute the number of boundaries
        return 1;
    }

    // Get iterator over first vertex of mesh's main border (aka "seam")
    Border_vertex_iterator  mesh_main_border_vertices_begin() {
        return m_seam_begin;
    }
    Border_vertex_const_iterator  mesh_main_border_vertices_begin() const {
        return (Border_vertex_const_iterator) m_seam_begin;
    }

    // Get iterator over past-the-end vertex of mesh's main border (aka "seam")
    Border_vertex_iterator  mesh_main_border_vertices_end() {
        return mesh_vertices_end();
    }
    Border_vertex_const_iterator  mesh_main_border_vertices_end() const {
        return mesh_vertices_end();
    }

    // Get iterator over first facet of mesh
    Facet_iterator  mesh_facets_begin() {
        return Facet_iterator(m_mesh_adaptor->mesh_facets_end(),
                              Inner_facets_filter(this),
                              m_mesh_adaptor->mesh_facets_begin());
    }
    Facet_const_iterator  mesh_facets_begin() const {
        return Facet_const_iterator(m_mesh_adaptor->mesh_facets_end(),
                                    Inner_facets_filter(this),
                                    m_mesh_adaptor->mesh_facets_begin());
    }

    // Get iterator over past-the-end facet of mesh
    Facet_iterator  mesh_facets_end() {
        return Facet_iterator(m_mesh_adaptor->mesh_facets_end(),
                              Inner_facets_filter(this));
    }
    Facet_const_iterator  mesh_facets_end() const {
        return Facet_const_iterator(m_mesh_adaptor->mesh_facets_end(),
                                    Inner_facets_filter(this));
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
        CGAL_parameterization_assertion(is_valid(facet));
        return Vertex_around_facet_circulator(this, m_mesh_adaptor->facet_vertices_begin(facet));
    }
    Vertex_around_facet_const_circulator  facet_vertices_begin(Facet_const_handle facet) const {
        CGAL_parameterization_assertion(is_valid(facet));
        return Vertex_around_facet_const_circulator(this, m_mesh_adaptor->facet_vertices_begin(facet));
    }

    // Count the number of vertices of a facet
    int  count_facet_vertices(Facet_const_handle facet) const {
        CGAL_parameterization_assertion(is_valid(facet));
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
        CGAL_parameterization_assertion(is_valid(vertex));
        return m_mesh_adaptor->get_vertex_position(vertex->get_adaptor_vertex());
    }

    // Get/set the 2D position (u/v pair) of a vertex
    Point_2  get_vertex_uv(Vertex_const_handle vertex) const {
        CGAL_parameterization_assertion(is_valid(vertex));
        return m_mesh_adaptor->get_corners_uv(vertex->get_adaptor_vertex(),
                                              vertex->get_prev_seam_vertex(),
                                              vertex->get_next_seam_vertex());
    }
    void  set_vertex_uv(Vertex_handle vertex, const Point_2& uv)
    {
#ifdef DEBUG_TRACE
        std::cerr << "    #" << get_vertex_index(vertex) << "(" << vertex->get_adaptor_vertex()->index() << ") <- (" << uv.x() << "," << uv.y() << ")\n";
#endif
        CGAL_parameterization_assertion(is_valid(vertex));
        return m_mesh_adaptor->set_corners_uv(vertex->get_adaptor_vertex(),
                                              vertex->get_prev_seam_vertex(),
                                              vertex->get_next_seam_vertex(),
                                              uv);
    }

    // Get/set "is parameterized" field of vertex
    bool  is_vertex_parameterized(Vertex_const_handle vertex) const {
        CGAL_parameterization_assertion(is_valid(vertex));
        return m_mesh_adaptor->are_corners_parameterized(vertex->get_adaptor_vertex(),
                                                         vertex->get_prev_seam_vertex(),
                                                         vertex->get_next_seam_vertex());
    }
    void  set_vertex_parameterized(Vertex_handle vertex, bool parameterized)
    {
        CGAL_parameterization_assertion(is_valid(vertex));
        return m_mesh_adaptor->set_corners_parameterized(vertex->get_adaptor_vertex(),
                                                         vertex->get_prev_seam_vertex(),
                                                         vertex->get_next_seam_vertex(),
                                                         parameterized);
    }

    // Get/set vertex index
    int  get_vertex_index(Vertex_const_handle vertex) const {
        CGAL_parameterization_assertion(is_valid(vertex));
        return m_mesh_adaptor->get_corners_index(vertex->get_adaptor_vertex(),
                                                 vertex->get_prev_seam_vertex(),
                                                 vertex->get_next_seam_vertex());
    }
    void  set_vertex_index(Vertex_handle vertex, int index)  {
        CGAL_parameterization_assertion(is_valid(vertex));
        return m_mesh_adaptor->set_corners_index(vertex->get_adaptor_vertex(),
                                                 vertex->get_prev_seam_vertex(),
                                                 vertex->get_next_seam_vertex(),
                                                 index);
    }

    //// Return true if a vertex belongs to ANY mesh's boundary
    //bool  is_vertex_on_border(Vertex_const_handle vertex) const {
    //    CGAL_parameterization_assertion(is_valid(vertex));
    //    return m_mesh_adaptor->is_vertex_on_border(vertex->get_adaptor_vertex());
    //}

    // Return true if a vertex belongs to the UNIQUE mesh's main border
    bool  is_vertex_on_main_border(Vertex_const_handle vertex) const {
        CGAL_parameterization_assertion(is_valid(vertex));
        return get_vertex_seaming(vertex) == BORDER;
    }

    // Get circulator over the vertices incident to 'vertex'
    // 'start_position' defines the optional initial position of the circulator
    Vertex_around_vertex_circulator vertices_around_vertex_begin(
                            Vertex_handle vertex,
                            Vertex_handle start_position = NULL)
    {
        CGAL_parameterization_assertion(is_valid(vertex));

        // Construct an adaptor circulator over the vertices
        // incident to vertex->get_adaptor_vertex()
        typename Adaptor::Vertex_around_vertex_circulator adaptor_circulator;
        if (start_position != NULL)
        {
            CGAL_parameterization_assertion(is_valid(start_position));
            adaptor_circulator = m_mesh_adaptor->vertices_around_vertex_begin(
                                            vertex->get_adaptor_vertex(),
                                            start_position->get_adaptor_vertex());
        }
        else
        {
            // If 'vertex' is a seam vertex, then vertex->get_prev_seam_vertex()
            // is a valid neighbor; else it is NULL, which is
            // a valid parameter for Adaptor::vertices_around_vertex_begin()
            adaptor_circulator = m_mesh_adaptor->vertices_around_vertex_begin(
                                            vertex->get_adaptor_vertex(),
                                            vertex->get_prev_seam_vertex());
        }

        return Vertex_around_vertex_circulator(this, vertex, adaptor_circulator);
    }
    Vertex_around_vertex_const_circulator vertices_around_vertex_begin(
                                    Vertex_const_handle vertex,
                                    Vertex_const_handle start_position = NULL) const
    {
        CGAL_parameterization_assertion(is_valid(vertex));

        // Construct an adaptor circulator over the vertices
        // incident to vertex->get_adaptor_vertex()
        typename Adaptor::Vertex_around_vertex_const_circulator adaptor_circulator;
        if (start_position != NULL)
        {
            CGAL_parameterization_assertion(is_valid(start_position));
            adaptor_circulator = m_mesh_adaptor->vertices_around_vertex_begin(
                                            vertex->get_adaptor_vertex(),
                                            start_position->get_adaptor_vertex());
        }
        else
        {
            // If 'vertex' is a seam vertex, then vertex->get_prev_seam_vertex()
            // is a valid neighbor; else it is NULL, which is
            // a valid parameter for Adaptor::vertices_around_vertex_begin()
            adaptor_circulator = m_mesh_adaptor->vertices_around_vertex_begin(
                                            vertex->get_adaptor_vertex(),
                                            vertex->get_prev_seam_vertex());
        }

        return Vertex_around_vertex_const_circulator(this, vertex, adaptor_circulator);
    }

    //******************************************************************
    // LEVEL 2 INTERFACE:
    // for classes attempting to parameterize (part of) 3D surfaces
    // of any genus and with any number of connected components.
    // They compute 1 (u,v) pair per corner.
    // Example: Mesh_adaptor_patch_3
    //******************************************************************

    // Mesh_adaptor_patch_3 exports no level 2 interface

// Private operations
private:

    // Default copy constructor and operator =() are not (yet) implemented
    Mesh_adaptor_patch_3(const Mesh_adaptor_patch_3& toCopy);
    Mesh_adaptor_patch_3& operator =(const Mesh_adaptor_patch_3& toCopy);

    // Set seaming flag of all vertices and edges to INNER, BORDER or OUTER
    // wrt the first_seam_vertex -> end_seam_vertex boundary
    // (outer seam edges are marked BORDER)
    //
    // Preconditions:
    // * first_seam_vertex -> end_seam_vertex defines the outer seam,
    //   ie Mesh_adaptor_patch_3 will export the "right" of the seam
    // * the "seam" is given as a container of Adaptor::Vertex_handle elements.
    template<class InputIterator>
    void set_mesh_seaming(InputIterator first_seam_vertex,
                          InputIterator end_seam_vertex)
    {
        fprintf(stderr, "  tag topological disc...");

        // Initialize the seaming flag of all vertices to OUTER
        m_mesh_adaptor->set_vertices_seaming(Adaptor::OUTER);

        // Initialize the seaming flag of all edges to OUTER
        m_mesh_adaptor->set_edges_seaming(Adaptor::OUTER);

        // Set seaming flag of seam vertices to BORDER.
        // Set seaming flag of outer seam edges to BORDER
        // and inner seam vertices to INNER.
        for (InputIterator border_it = first_seam_vertex;
             border_it != end_seam_vertex;
             border_it++)
        {
            // Set vertex seaming flag
            m_mesh_adaptor->set_vertex_seaming(*border_it,
                                               Adaptor::BORDER);

            // Get next iterator (looping)
            InputIterator next_border_it = border_it;
            next_border_it++;
            if (next_border_it == end_seam_vertex)
                next_border_it = first_seam_vertex;

            // Set outer seam edge to BORDER
            m_mesh_adaptor->set_edge_seaming(*border_it, *next_border_it,
                                             Adaptor::BORDER);

            // Set inner seam edge to INNER (except if also BORDER)
            if (m_mesh_adaptor->get_edge_seaming(*next_border_it,
                                                 *border_it) != Adaptor::BORDER) {
                m_mesh_adaptor->set_edge_seaming(*next_border_it, *border_it,
                                                Adaptor::INNER);
            }
        }

        // Set the seaming flag of inner vertices and edges to INNER
        for (InputIterator border_it = first_seam_vertex;
             border_it != end_seam_vertex;
             border_it++)
        {
            // Get next iterator (looping)
            InputIterator next_border_it = border_it;
            next_border_it++;
            if (next_border_it == end_seam_vertex)
                next_border_it = first_seam_vertex;

            // Get inner point at the "right" of *border_it
            // by a counter-clockwise rotation around the next seam vertex
            typename Adaptor::Vertex_around_vertex_circulator cir =
                m_mesh_adaptor->vertices_around_vertex_begin(*next_border_it,
                                                             *border_it);
            cir--;

            // Fill topological disk
            if (m_mesh_adaptor->get_vertex_seaming(cir) != Adaptor::BORDER)
                set_inner_region_seaming(cir);
        }

        fprintf(stderr,"ok\n");
    }

    // Set the seaming flag of inner vertices and edges to INNER
    // by filling the topological disk
    //
    // Preconditions:
    // * Inner vertices are marked as OUTER, seam vertices as BORDER
    // * Inner edges are marked as OUTER,
    //   outer seam edges as BORDER, inner seam edges as INNER
    // * pSeedVertex is in the inner region
    // * pSeedVertex != NULL
    //
    // Implementation note:
    // The seaming status of inner edges is unused, thus this part is not tested
    //
    void set_inner_region_seaming(typename Adaptor::Vertex_handle pSeedVertex)
    {
        if (pSeedVertex == NULL)
            return;                 // Gloups... topological disc is empty!

        // List of vertices to flag = pSeedVertex initially
        std::list<typename Adaptor::Vertex_handle> vertices;
        vertices.push_front(pSeedVertex);

        // For each vertex in the list: pop it out, flag it as INNER and
        // add its surrounding vertices to the list
        while (!vertices.empty())
        {
            typename Adaptor::Vertex_handle pVertex = vertices.front();
            vertices.pop_front();
            CGAL_parameterization_assertion(pVertex != NULL);

            // Flag this vertex as INNER
            if (m_mesh_adaptor->get_vertex_seaming(pVertex) == Adaptor::OUTER)
                m_mesh_adaptor->set_vertex_seaming(pVertex, Adaptor::INNER);
            else
                continue;           // Skip this vertex if it is already done

            // For each neighbor vertex
            typename Adaptor::Vertex_around_vertex_circulator cir, cir_end;
            cir     = m_mesh_adaptor->vertices_around_vertex_begin(pVertex);
            cir_end = cir;
            CGAL_For_all(cir, cir_end)
            {
                // Flag both oriented edges pVertex <-> cir
                m_mesh_adaptor->set_edge_seaming(pVertex, cir, Adaptor::INNER);
                m_mesh_adaptor->set_edge_seaming(cir, pVertex, Adaptor::INNER);

                // Add surrounding vertices to list without crossing the border
                if (m_mesh_adaptor->get_vertex_seaming(cir) == Adaptor::OUTER)
                    vertices.push_front(cir);
            }
        }
    }

    // Get facet' seaming status (INNER or OUTER)
    Seaming_status get_facet_seaming(typename Adaptor::Facet_const_handle facet) const
    {
        // don't call is_valid() to avoid an infinite loop
        CGAL_parameterization_assertion(facet != NULL);

        typename Adaptor::Vertex_around_facet_const_circulator
                            cir = m_mesh_adaptor->facet_vertices_begin(facet);
        CGAL_parameterization_assertion(cir != NULL);
        return (m_mesh_adaptor->get_vertex_seaming(cir) == Adaptor::OUTER) ?
               OUTER :
               INNER;
    }

    // Get/set vertex seaming flag,
    // ie position of the vertex wrt to the UNIQUE main boundary
    Seaming_status  get_vertex_seaming(Vertex_const_handle vertex) const {
        // don't call is_valid() to avoid an infinite loop
        CGAL_parameterization_assertion(vertex != NULL);

        return (Seaming_status) m_mesh_adaptor->get_vertex_seaming(
                                                    vertex->get_adaptor_vertex());
    }
    void set_vertex_seaming(Vertex_handle vertex, Seaming_status seaming) {
        m_mesh_adaptor->set_vertex_seaming(vertex->get_adaptor_vertex(),
                                           (typename Adaptor::Seaming_status)seaming);
    }

    // Create a patch vertex from an adaptor vertex + one of its neighbors
    //
    // Preconditions:
    // * adaptor_neighbor is a neighbor of adaptor_vertex
    // * (adaptor_vertex, adaptor_neighbor) must NOT be a seam (non-oriented) edge
    Vertex_const_handle get_decorated_vertex(
                    typename Adaptor::Vertex_const_handle adaptor_vertex,
                    typename Adaptor::Vertex_const_handle adaptor_neighbor) const
    {
        Vertex vertex;                      // returned variable

        // We need at least an inner neighbor as input
        assert(m_mesh_adaptor->get_edge_seaming(adaptor_vertex,
                                                adaptor_neighbor) != Adaptor::BORDER
            || m_mesh_adaptor->get_edge_seaming(adaptor_neighbor,
                                                adaptor_vertex) != Adaptor::BORDER);

        // if inner vertex
        if ( ! m_mesh_adaptor->is_vertex_on_main_border(adaptor_vertex) )
        {
            // No extra information needed if inner vertex
            vertex = Vertex((typename Adaptor::Vertex*)&*adaptor_vertex);
        }
        else // if seam vertex
        {
            // find previous vertex on seam by a clockwise rotation
            typename Adaptor::Vertex_around_vertex_const_circulator prev_seam_vertex_cir
                = m_mesh_adaptor->vertices_around_vertex_begin(adaptor_vertex,
                                                               adaptor_neighbor);
            while (m_mesh_adaptor->get_edge_seaming(prev_seam_vertex_cir,
                                                    adaptor_vertex) != Adaptor::BORDER)
            {
                prev_seam_vertex_cir++;
            }

            // find next vertex on seam by a counter-clockwise rotation
            typename Adaptor::Vertex_around_vertex_const_circulator next_seam_vertex_cir
                = m_mesh_adaptor->vertices_around_vertex_begin(adaptor_vertex,
                                                               adaptor_neighbor);
            while (m_mesh_adaptor->get_edge_seaming(adaptor_vertex,
                                                    next_seam_vertex_cir) != Adaptor::BORDER)
            {
                next_seam_vertex_cir--;
            }

            // The decorated vertex is then:
            vertex = Vertex((typename Adaptor::Vertex*)&*adaptor_vertex,
                            (typename Adaptor::Vertex*)&*prev_seam_vertex_cir,
                            (typename Adaptor::Vertex*)&*next_seam_vertex_cir);
        }

        // Implementation note:
        // The next line seems to return a reference to a local Vertex variable.
        // In fact, Vertex_[const_]handle constructor copies the Vertex object.
        // The purpose is to save the time of searching the Vertex in
        // m_inner_and_border_vertices list.
        return &vertex;
    }
    Vertex_handle get_decorated_vertex(
                        typename Adaptor::Vertex_handle adaptor_vertex,
                        typename Adaptor::Vertex_const_handle adaptor_neighbor)
    {
        Vertex_const_handle vertex_hdl = get_decorated_vertex(
                        (typename Adaptor::Vertex_const_handle) &*adaptor_vertex,
                        adaptor_neighbor);
        return (Vertex*) &*vertex_ref;
    }

    // Debug utility: Check if a Mesh_adaptor_patch_3 facet is valid
    bool is_valid(Facet_const_handle facet) const
    {
        if (facet == NULL)
            return false;
        // outer facets are not exported
        if (get_facet_seaming(facet) != INNER)
            return false;
        // eventually: ok
        return true;
    }

    // Debug utility: Check if a Mesh_adaptor_patch_3 vertex is valid
    bool is_valid(Vertex_const_handle vertex) const
    {
        if (vertex == NULL)
            return false;
        // outer halfedges are not exported
        if (get_vertex_seaming(vertex) == OUTER)
            return false;
        // eventually: ok
        return true;
    }

    // Debug utility: get vertex index as string ("-" if NULL vertex)
    std::string get_vertex_index_as_string(typename Adaptor::Vertex_const_handle vertex) const
    {
        if (vertex == NULL) {
            return std::string("-");
        } else {
            char index_as_string[64];
            sprintf(index_as_string, "%d", (int)m_mesh_adaptor->get_vertex_index(vertex));
            return std::string(index_as_string);
        }
    }

// Fields
private:

    // The decorated mesh
    Adaptor* m_mesh_adaptor;

    // List of all exported vertices
    // Order is: inner vertices, then seam/main boundary ones
    Mesh_adaptor_patch_vertex_list<Adaptor> m_inner_and_border_vertices;

    // Iterator to first seam vertex inside m_inner_and_border_vertices
    Border_vertex_iterator m_seam_begin;

// Private types
private:

    // Utility class to generate the Facet_iterator type
    struct Inner_facets_filter
    {
        Inner_facets_filter(const Mesh_adaptor_patch_3* mesh) : m_mesh_patch(mesh) {}

        // Return TRUE <=> the facet IS NOT EXPORTED by Mesh_adaptor_patch_3,
        // ie is out of the topological disc
        bool operator()(const typename Adaptor::Facet_iterator& f) const       {
            return m_mesh_patch->get_facet_seaming(f) == OUTER;
        }
        bool operator()(const typename Adaptor::Facet_const_iterator& f) const {
            return m_mesh_patch->get_facet_seaming(f) == OUTER;
        }

    private:
        const Mesh_adaptor_patch_3* m_mesh_patch;
    };

}; // Mesh_adaptor_patch_3


CGAL_END_NAMESPACE

#endif //CGAL_MESH_ADAPTOR_PATCH_3_H

