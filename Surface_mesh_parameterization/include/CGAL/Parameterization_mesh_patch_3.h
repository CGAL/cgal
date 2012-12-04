// Copyright (c) 2005  INRIA (France).
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
//
// Author(s)     : Laurent Saboret, Pierre Alliez, Bruno Levy


#ifndef CGAL_PARAMETERIZATION_MESH_PATCH_3_H
#define CGAL_PARAMETERIZATION_MESH_PATCH_3_H

#include <cstdio>
#include <CGAL/iterator.h>
#include <CGAL/circulator.h>
#include <CGAL/Convertible_iterator_project.h>
#include <CGAL/Convertible_circulator_project.h>
#include <CGAL/Convertible_filter_iterator.h>
#include <CGAL/Param_mesh_patch_vertex.h>
#include <CGAL/Param_mesh_patch_iterators.h>
#include <CGAL/Param_mesh_patch_circulators.h>

#include <CGAL/surface_mesh_parameterization_assertions.h>

/// \file Parameterization_mesh_patch_3.h

namespace CGAL {


/// \ingroup  PkgSurfaceParameterizationMesh
///
/// Parameterization_mesh_patch_3 is a Decorator class to <i>virtually</i> cut a patch
/// in a ParameterizationPatchableMesh_3 3D surface. Only the patch is exported,
/// making the 3D surface look like a topological disk.
///
/// The input mesh can be of any genus, but it has to come with a <i>seam</i> that
/// describes the border of a topological disc. This border may be an actual
/// border of the mesh or a virtual border.
///
/// \cgalModels `ParameterizationMesh_3`
///
///
/// \tparam ParameterizationPatchableMesh_3 3D surface mesh.

template<class ParameterizationPatchableMesh_3>
class Parameterization_mesh_patch_3
{
// Private types
private:

    // Forward references
    struct                                  Inner_facets_filter;

    /// Seaming flag.
    enum Seaming_status  { OUTER, INNER, BORDER };

public:

    /// Export template parameter.
    typedef ParameterizationPatchableMesh_3 Adaptor;

    #ifndef DOXYGEN_RUNNING
    /// \name Types implementing the ParameterizationMesh_3 interface
    /// @{

    /// Number type to represent coordinates.
    typedef typename Adaptor::NT            NT;

    /// 2D point that represents `(u,v)` coordinates computed
    /// by parameterization methods. Must provide `x()` and `y()` methods.
    typedef typename Adaptor::Point_2       Point_2;
    /// 3D point that represents vertices coordinates. Must provide `x()` and `y()` methods.
    typedef typename Adaptor::Point_3       Point_3;
    /// 2D vector. Must provide `x()` and `y()` methods.
    typedef typename Adaptor::Vector_2      Vector_2;
    /// 3D vector. Must provide `x()` and `y()` methods.
    typedef typename Adaptor::Vector_3      Vector_3;

    /// Opaque type representing a facet of the 3D mesh. No methods are expected.
    typedef typename Adaptor::Facet         Facet;
    /// Handle to a facet. Model of the Handle concept.
    typedef typename Adaptor::Facet_handle  Facet_handle;
    typedef typename Adaptor::Facet_const_handle
                                            Facet_const_handle;
    /// Iterator over all mesh facets. Model of the `ForwardIterator` concept.
    typedef Convertible_filter_iterator<typename Adaptor::Facet_iterator,
                                        Inner_facets_filter,
                                        Facet_const_handle,
                                        Facet_handle>
                                            Facet_iterator;
    typedef Convertible_filter_iterator<typename Adaptor::Facet_const_iterator,
                                        Inner_facets_filter,
                                        Facet_const_handle>
                                            Facet_const_iterator;

    /// Opaque type representing a vertex of the 3D mesh. No methods are expected.
    typedef Param_mesh_patch_vertex<Adaptor> Vertex;
    /// Handle to a vertex. Model of the `Handle` concept.
    typedef Param_mesh_patch_vertex_handle<Adaptor>
                                            Vertex_handle;
    typedef Param_mesh_patch_vertex_const_handle<Adaptor>
                                            Vertex_const_handle;
    /// Iterator over all vertices of a mesh. Model of the `ForwardIterator` concept.
    typedef Param_mesh_patch_vertex_list_iterator<Adaptor>
                                            Vertex_iterator;
    typedef Param_mesh_patch_vertex_list_const_iterator<Adaptor>
                                            Vertex_const_iterator;
    /// Iterator over vertices of the mesh "main border".
    /// Model of the `ForwardIterator` concept.
    typedef Vertex_iterator                 Border_vertex_iterator;
    typedef Vertex_const_iterator           Border_vertex_const_iterator;
    /// Counter-clockwise circulator over a facet's vertices.
    /// Model of the `BidirectionalCirculator` concept.
    typedef Mesh_patch_vertex_around_facet_cir<Parameterization_mesh_patch_3,
                                               Vertex_handle,
                                               typename Adaptor::Vertex_around_facet_circulator>
                                            Vertex_around_facet_circulator;
    typedef Mesh_patch_vertex_around_facet_cir<const Parameterization_mesh_patch_3,
                                               Vertex_const_handle,
                                               typename Adaptor::Vertex_around_facet_const_circulator>
                                            Vertex_around_facet_const_circulator;
    /// Clockwise circulator over the vertices incident to a vertex.
    /// Model of the `BidirectionalCirculator` concept.
    typedef Mesh_patch_vertex_around_vertex_cir<Parameterization_mesh_patch_3,
                                                Vertex_handle,
                                                typename Adaptor::Vertex_around_vertex_circulator,
                                                typename Adaptor::Vertex_handle>
                                            Vertex_around_vertex_circulator;
    typedef Mesh_patch_vertex_around_vertex_cir<const Parameterization_mesh_patch_3,
                                                Vertex_const_handle,
                                                typename Adaptor::Vertex_around_vertex_const_circulator,
                                                typename Adaptor::Vertex_const_handle>
                                            Vertex_around_vertex_const_circulator;

    /// @} // end of Types implementing the ParameterizationMesh_3 interface
    #endif //DOXYGEN_RUNNING
// Public operations
public:

    /// Create a Decorator for an existing ParameterizationPatchableMesh_3 mesh.
    /// The input mesh can be of any genus, but it has to come with a <i>seam</i> that
    /// describes the border of a topological disc. This border may be an actual
    /// border of the mesh or a virtual border.
    ///
    /// \pre `first_seam_vertex -> end_seam_vertex` defines the outer seam,
    ///   i.e. Parameterization_mesh_patch_3 will export the <i>right</i> of the seam.
    /// \pre The <i>seam</i> is given as a container of `Adaptor::Vertex_handle` elements.
    /// \pre The <i>seam</i> is implicitely a loop. The first vertex must *not* be
    ///   duplicated at the end.
    template<class InputIterator>
    Parameterization_mesh_patch_3(Adaptor& mesh,
                                  InputIterator first_seam_vertex,
                                  InputIterator end_seam_vertex)
        // Store reference to adapted mesh
      : m_mesh_adaptor(mesh)
    {
// #ifdef DEBUG_TRACE
//         // Dump input border (for debug purpose)
//         fprintf(stderr,"  input border is: ");
//         for (InputIterator it = first_seam_vertex; it != end_seam_vertex; it++)
//             fprintf(stderr, "%s ", get_vertex_index_as_string(*it).c_str());
//         fprintf(stderr,"ok\n");
// #endif

        // Set seaming flag of all vertices and edges to INNER, BORDER or OUTER
        // w.r.t. the first_seam_vertex -> end_seam_vertex border
        set_mesh_seaming(first_seam_vertex, end_seam_vertex);

        // Check that the cut mesh is 2-manifold
        m_is_valid = mesh.is_valid() && check_seam(first_seam_vertex, end_seam_vertex);

        // Construct the list of all exported vertices, i.e. INNER and BORDER vertices
        //
        // 1) add inner vertices
        for (typename Adaptor::Vertex_iterator it = mesh.mesh_vertices_begin();
             it != mesh.mesh_vertices_end();
             it++)
        {
            if (m_mesh_adaptor.get_vertex_seaming(it) == INNER)
                m_inner_and_border_vertices.push_back( Vertex(it) );
        }
        // 2) add seam vertices, w.r.t. outer seam/border order
        InputIterator border_it      = first_seam_vertex;
        InputIterator prev_border_it = end_seam_vertex; prev_border_it--;
        InputIterator next_border_it = first_seam_vertex; next_border_it++;
        while (border_it != end_seam_vertex)
        {
            // Get outer border vertex
            Vertex v;
            // if non main border vertex
            if (m_mesh_adaptor.get_halfedge_seaming(*border_it, *prev_border_it) != BORDER)
                v = Vertex(*border_it, *prev_border_it, *next_border_it);
            else // if main border (aka seam) vertex
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

        // Initialize m_seam_begin = iterator to beginning of seam/main border
        // inside m_inner_and_border_vertices
        m_seam_begin = mesh_vertices_end();
        for (Vertex_iterator it=mesh_vertices_begin(); it!=mesh_vertices_end(); it++) {
            if (get_vertex_seaming(it) == BORDER) {
                m_seam_begin = it;
                break;
            }
        }

#ifndef CGAL_NDEBUG
        // Index vertices right away to ease debugging
        index_mesh_vertices();

/*    #ifdef DEBUG_TRACE
        // Dump seam (for debug purpose)
        fprintf(stderr,"  seam is: ");
        for (Border_vertex_iterator border_it = mesh_main_border_vertices_begin();
             border_it != mesh_main_border_vertices_end();
             border_it++)
        {
            fprintf(stderr, "#%d ", get_vertex_index(border_it));
        }
        fprintf(stderr,"ok\n");
    #endif*/
#endif
    }

    /// @return the decorated mesh.
    Adaptor&       get_decorated_mesh()       { return m_mesh_adaptor; }
    const Adaptor& get_decorated_mesh() const { return m_mesh_adaptor; }

    /// \name Methods implementing the ParameterizationMesh_3 interface
    /// @{

    // MESH INTERFACE

    /// Indicate if the mesh matches the ParameterizationMesh_3 concept.
    bool is_valid() const {
        return m_is_valid;
    }

    /// Get iterator over first vertex of mesh.
    Vertex_iterator  mesh_vertices_begin() {
        return m_inner_and_border_vertices.begin();
    }
    Vertex_const_iterator  mesh_vertices_begin() const {
        return m_inner_and_border_vertices.begin();
    }

    /// Get iterator over past-the-end vertex of mesh.
    Vertex_iterator  mesh_vertices_end() {
        return m_inner_and_border_vertices.end();
    }
    Vertex_const_iterator  mesh_vertices_end() const {
        return m_inner_and_border_vertices.end();
    }

    /// Count the number of vertices of the mesh.
    int  count_mesh_vertices() const {
        int index = 0;
        for (Vertex_const_iterator it=mesh_vertices_begin(); it!=mesh_vertices_end(); it++)
            index++;
        return index;
    }

    /// Index vertices of the mesh from 0 to count_mesh_vertices()-1.
    void  index_mesh_vertices ()
    {
#ifdef DEBUG_TRACE
        fprintf(stderr,"  index Parameterization_mesh_patch_3 vertices:\n");
#endif
        int index = 0;
        for (Vertex_iterator it=mesh_vertices_begin(); it!=mesh_vertices_end(); it++)
        {
/*#ifdef DEBUG_TRACE
            fprintf(stderr, "    #%d = {%s,%s,%s}\n",
                            index,
                            get_vertex_index_as_string(it->vertex()).c_str(),
                            get_vertex_index_as_string(it->last_cw_neighbor()).c_str(),
                            get_vertex_index_as_string(it->first_cw_neighbor()).c_str());
#endif*/
            set_vertex_index(it, index++);
        }
#ifdef DEBUG_TRACE
        fprintf(stderr,"    ok\n");
#endif
    }

    /// Get iterator over first vertex of mesh's main border (aka <i>seam</i>).
    Border_vertex_iterator  mesh_main_border_vertices_begin() {
        return m_seam_begin;
    }
    Border_vertex_const_iterator  mesh_main_border_vertices_begin() const {
        return (Border_vertex_const_iterator) m_seam_begin;
    }

    /// Get iterator over past-the-end vertex of mesh's main border (aka <i>seam</i>).
    Border_vertex_iterator  mesh_main_border_vertices_end() {
        return mesh_vertices_end();
    }
    Border_vertex_const_iterator  mesh_main_border_vertices_end() const {
        return mesh_vertices_end();
    }

    /// @return the border containing seed_vertex (or an empty list if not found).
    /// @param seed_vertex a border vertex.
    std::list<Vertex_handle> get_border(Vertex_handle seed_vertex)
    {
        std::list<Vertex_handle> border;    // returned list

        // If seam vertex, return the seam
        if (is_vertex_on_main_border(seed_vertex))
        {
            for (Border_vertex_iterator it = mesh_main_border_vertices_begin();
                 it != mesh_main_border_vertices_end();
                 it++)
            {
                CGAL_surface_mesh_parameterization_assertion(is_vertex_on_main_border(it));
                border.push_back(it);
            }
        }
        else // if vertex on the border of a hole
        {
            // Get list of vertices on this border
            std::list<typename Adaptor::Vertex_handle> adaptor_border =
                m_mesh_adaptor.get_border(seed_vertex->vertex());

            // Copy them into 'border'
            for (typename std::list<typename Adaptor::Vertex_handle>::iterator it = adaptor_border.begin();
                 it != adaptor_border.end();
                 it++)
            {
                CGAL_surface_mesh_parameterization_assertion(is_vertex_on_border( Vertex_handle(*it) ));
                border.push_back( Vertex_handle(*it) );
            }
        }

        return border;
    }

    /// Get iterator over first facet of mesh.
    Facet_iterator  mesh_facets_begin() {
        return Facet_iterator(m_mesh_adaptor.mesh_facets_end(),
                              Inner_facets_filter(*this),
                              m_mesh_adaptor.mesh_facets_begin());
    }
    Facet_const_iterator  mesh_facets_begin() const {
        return Facet_const_iterator(m_mesh_adaptor.mesh_facets_end(),
                                    Inner_facets_filter(*this),
                                    m_mesh_adaptor.mesh_facets_begin());
    }

    /// Get iterator over past-the-end facet of mesh.
    Facet_iterator  mesh_facets_end() {
        return Facet_iterator(m_mesh_adaptor.mesh_facets_end(),
                              Inner_facets_filter(*this));
    }
    Facet_const_iterator  mesh_facets_end() const {
        return Facet_const_iterator(m_mesh_adaptor.mesh_facets_end(),
                                    Inner_facets_filter(*this));
    }

    /// Count the number of facets of the mesh.
    int  count_mesh_facets() const {
        int index = 0;
        for (Facet_const_iterator it=mesh_facets_begin(); it!=mesh_facets_end(); it++)
            index++;
        return index;
    }

    /// Return true of all mesh's facets are triangles.
    bool  is_mesh_triangular() const {
        for (Facet_const_iterator it = mesh_facets_begin(); it != mesh_facets_end(); it++)
            if (count_facet_vertices(it) != 3)
                return false;
        return true;            // mesh is triangular if we reach this point
    }

    /// Count the number of halfedges of the mesh.
    int  count_mesh_halfedges() const {
        int index = 0;
        for (Vertex_const_iterator it=mesh_vertices_begin(); it!=mesh_vertices_end(); it++)
        {
            // Count each neighbor vertex
            Vertex_around_vertex_const_circulator cir = vertices_around_vertex_begin(it);
            Vertex_around_vertex_const_circulator cir_end = cir;
            CGAL_For_all(cir, cir_end)
                index++;
        }
        return index;
    }

    // FACET INTERFACE

    /// Get circulator over facet's vertices.
    Vertex_around_facet_circulator  facet_vertices_begin(Facet_handle facet) {
        CGAL_surface_mesh_parameterization_assertion(is_valid(facet));
        return Vertex_around_facet_circulator(*this, m_mesh_adaptor.facet_vertices_begin(facet));
    }
    Vertex_around_facet_const_circulator  facet_vertices_begin(Facet_const_handle facet) const {
        CGAL_surface_mesh_parameterization_assertion(is_valid(facet));
        return Vertex_around_facet_const_circulator(*this, m_mesh_adaptor.facet_vertices_begin(facet));
    }

    /// Count the number of vertices of a facet.
    int  count_facet_vertices(Facet_const_handle facet) const {
        CGAL_surface_mesh_parameterization_assertion(is_valid(facet));
        int index = 0;
        Vertex_around_facet_const_circulator cir     = facet_vertices_begin(facet),
                                             cir_end = cir;
        CGAL_For_all(cir, cir_end)
            index++;
        return index;
    }

    // VERTEX INTERFACE

    /// Get the 3D position of a vertex.
    Point_3 get_vertex_position(Vertex_const_handle vertex) const {
        CGAL_surface_mesh_parameterization_assertion(is_valid(vertex));
        return m_mesh_adaptor.get_vertex_position(vertex->vertex());
    }

    /// Get/set the 2D position (u/v pair) of a vertex. Default value is undefined.
    Point_2  get_vertex_uv(Vertex_const_handle vertex) const {
        CGAL_surface_mesh_parameterization_assertion(is_valid(vertex));
        return m_mesh_adaptor.get_corners_uv(vertex->vertex(),
                                             vertex->last_cw_neighbor(),
                                             vertex->first_cw_neighbor());
    }
    void  set_vertex_uv(Vertex_handle vertex, const Point_2& uv)
    {
        CGAL_surface_mesh_parameterization_assertion(is_valid(vertex));
        return m_mesh_adaptor.set_corners_uv(vertex->vertex(),
                                             vertex->last_cw_neighbor(),
                                             vertex->first_cw_neighbor(),
                                             uv);
    }

    /// Get/set "is parameterized" field of vertex. Default value is undefined.
    bool  is_vertex_parameterized(Vertex_const_handle vertex) const {
        CGAL_surface_mesh_parameterization_assertion(is_valid(vertex));
        return m_mesh_adaptor.are_corners_parameterized(vertex->vertex(),
                                                        vertex->last_cw_neighbor(),
                                                        vertex->first_cw_neighbor());
    }
    void  set_vertex_parameterized(Vertex_handle vertex, bool parameterized)
    {
        CGAL_surface_mesh_parameterization_assertion(is_valid(vertex));
        return m_mesh_adaptor.set_corners_parameterized(vertex->vertex(),
                                                        vertex->last_cw_neighbor(),
                                                        vertex->first_cw_neighbor(),
                                                        parameterized);
    }

    /// Get/set vertex index. Default value is undefined.
    int  get_vertex_index(Vertex_const_handle vertex) const {
        CGAL_surface_mesh_parameterization_assertion(is_valid(vertex));
        return m_mesh_adaptor.get_corners_index(vertex->vertex(),
                                                vertex->last_cw_neighbor(),
                                                vertex->first_cw_neighbor());
    }
    void  set_vertex_index(Vertex_handle vertex, int index)  {
        CGAL_surface_mesh_parameterization_assertion(is_valid(vertex));
        return m_mesh_adaptor.set_corners_index(vertex->vertex(),
                                                vertex->last_cw_neighbor(),
                                                vertex->first_cw_neighbor(),
                                                index);
    }

    /// Get/set vertex' all purpose tag. Default value is undefined.
    int  get_vertex_tag(Vertex_const_handle vertex) const {
        CGAL_surface_mesh_parameterization_assertion(is_valid(vertex));
        return m_mesh_adaptor.get_corners_tag(vertex->vertex(),
                                              vertex->last_cw_neighbor(),
                                              vertex->first_cw_neighbor());
    }
    void set_vertex_tag(Vertex_handle vertex, int tag) {
        return m_mesh_adaptor.set_corners_tag(vertex->vertex(),
                                              vertex->last_cw_neighbor(),
                                              vertex->first_cw_neighbor(),
                                              tag);
    }

    /// Return `true` if `vertex` belongs to the border of any mesh.
    bool  is_vertex_on_border(Vertex_const_handle vertex) const {
        CGAL_surface_mesh_parameterization_assertion(is_valid(vertex));
        return is_vertex_on_main_border(vertex) ||
               m_mesh_adaptor.is_vertex_on_border(vertex->vertex());
    }

    /// Return `true` if `vertex` belongs to the UNIQUE mesh's main border
    /// set by the constructor.
    bool  is_vertex_on_main_border(Vertex_const_handle vertex) const {
        CGAL_surface_mesh_parameterization_assertion(is_valid(vertex));
        return get_vertex_seaming(vertex) == BORDER;
    }

    /// Get circulator over the vertices incident to `vertex`.
    /// `start_position` defines the optional initial position of the circulator.
    Vertex_around_vertex_circulator vertices_around_vertex_begin(
                            Vertex_handle vertex,
                            Vertex_handle start_position = Vertex_handle())
    {
        CGAL_surface_mesh_parameterization_assertion(is_valid(vertex));
        CGAL_surface_mesh_parameterization_assertion(start_position == NULL ||
                                        is_valid(start_position));

        // If no start position provided, pick one
        if (start_position == NULL)
        {
            // If 'vertex' is an inner vertex, pick any neighbor
            if (vertex->last_cw_neighbor() == NULL)
            {
                typename Adaptor::Vertex_around_vertex_circulator adaptor_circulator
                    = m_mesh_adaptor.vertices_around_vertex_begin(vertex->vertex());
                start_position = get_decorated_vertex_from_inner_edge(adaptor_circulator,
                                                                      vertex->vertex());
            }
            else // If 'vertex' is a seam vertex, pick its last clockwise neighbor
            {
                start_position = get_decorated_vertex_from_border_edge(vertex->last_cw_neighbor(),
                                                                       NULL,
                                                                       vertex->vertex());
            }
        }

        return Vertex_around_vertex_circulator(*this, vertex, start_position);
    }
    Vertex_around_vertex_const_circulator vertices_around_vertex_begin(
                                Vertex_const_handle vertex,
                                Vertex_const_handle start_position = Vertex_const_handle()) const
    {
        CGAL_surface_mesh_parameterization_assertion(is_valid(vertex));
        CGAL_surface_mesh_parameterization_assertion(start_position == NULL ||
                                        is_valid(start_position));

        // If no start position provided, pick one
        if (start_position == NULL)
        {
            // If 'vertex' is an inner vertex, pick any neighbor
            if (vertex->last_cw_neighbor() == NULL)
            {
                typename Adaptor::Vertex_around_vertex_const_circulator adaptor_circulator
                    = m_mesh_adaptor.vertices_around_vertex_begin(vertex->vertex());
                start_position = get_decorated_vertex_from_inner_edge(adaptor_circulator,
                                                                      vertex->vertex());
            }
            else // If 'vertex' is a seam vertex, pick its last clockwise neighbor
            {
                start_position = get_decorated_vertex_from_border_edge(vertex->last_cw_neighbor(),
                                                                       NULL,
                                                                       vertex->vertex());
            }
        }

        return Vertex_around_vertex_const_circulator(*this, vertex, start_position);
    }

    /// @} // end of Methods implementing the ParameterizationMesh_3 interface


// Private operations
private:

    /// Copy constructor and operator =() are not implemented.
    Parameterization_mesh_patch_3(const Parameterization_mesh_patch_3& toCopy);
    Parameterization_mesh_patch_3& operator =(const Parameterization_mesh_patch_3& toCopy);

    /// Set seaming flag of all vertices and edges to INNER, BORDER or OUTER
    /// w.r.t. the first_seam_vertex -> end_seam_vertex border
    /// (outer seam edges are marked BORDER).
    ///
    /// \pre first_seam_vertex -> end_seam_vertex defines the outer seam,
    ///   i.e. Parameterization_mesh_patch_3 will export the <i>right</i> of the seam.
    /// \pre The <i>seam</i> is given as a container of Adaptor::Vertex_handle elements.
    /// \pre The <i>seam</i> is implicitely a loop. The first vertex should *not* be
    ///   duplicated at the end.
    template<class InputIterator>
    void set_mesh_seaming(InputIterator first_seam_vertex,
                          InputIterator end_seam_vertex)
    {
#ifdef DEBUG_TRACE
        fprintf(stderr, "  tag topological disc...");
#endif

        // Initialize the seaming flag of all vertices to OUTER
        for (typename Adaptor::Vertex_iterator it = m_mesh_adaptor.mesh_vertices_begin();
             it != m_mesh_adaptor.mesh_vertices_end();
             it++)
        {
             m_mesh_adaptor.set_vertex_seaming(it, OUTER);
        }

        // Initialize the seaming flag of all halfedges to OUTER
        for (typename Adaptor::Vertex_iterator it = m_mesh_adaptor.mesh_vertices_begin();
             it != m_mesh_adaptor.mesh_vertices_end();
             it++)
        {
            // For each neighbor vertex
            typename Adaptor::Vertex_around_vertex_circulator cir, cir_end;
            cir     = m_mesh_adaptor.vertices_around_vertex_begin(it);
            cir_end = cir;
            CGAL_For_all(cir, cir_end)
                m_mesh_adaptor.set_halfedge_seaming(it, cir, OUTER);
        }

        // Set seaming flag of seam vertices to BORDER.
        // Set seaming flag of outer seam edges to BORDER
        // and inner seam edges to INNER.
        for (InputIterator border_it = first_seam_vertex;
             border_it != end_seam_vertex;
             border_it++)
        {
            // Set vertex seaming flag
            m_mesh_adaptor.set_vertex_seaming(*border_it, BORDER);

            // Get next iterator (looping)
            InputIterator next_border_it = border_it;
            next_border_it++;
            if (next_border_it == end_seam_vertex)
                next_border_it = first_seam_vertex;

            // Set outer seam edge to BORDER
            m_mesh_adaptor.set_halfedge_seaming(*border_it, *next_border_it,
                                                BORDER);

            // Set inner seam edge to INNER (except if also BORDER)
            if (m_mesh_adaptor.get_halfedge_seaming(*next_border_it,
                                                    *border_it) != BORDER) {
                m_mesh_adaptor.set_halfedge_seaming(*next_border_it, *border_it,
                                                    INNER);
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
                m_mesh_adaptor.vertices_around_vertex_begin(*next_border_it,
                                                            *border_it);
            cir--;

            // Fill topological disk
            if (m_mesh_adaptor.get_vertex_seaming(cir) != BORDER)
                set_inner_region_seaming(cir);
        }

#ifdef DEBUG_TRACE
        fprintf(stderr,"ok\n");
#endif
    }

    /// Set the seaming flag of inner vertices and edges to INNER
    /// by filling the topological disk.
    ///
    /// \pre Inner vertices are marked as OUTER, seam vertices as BORDER.
    /// \pre Inner edges are marked as OUTER,
    ///   outer seam edges as BORDER, inner seam edges as INNER.
    /// \pre pSeedVertex is in the inner region.
    /// \pre pSeedVertex != NULL.
    ///
    /// Implementation note:
    /// The seaming status of inner edges is unused, thus this part is not tested.
    ///
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
            CGAL_surface_mesh_parameterization_assertion(pVertex != NULL);

            // Flag this vertex as INNER
            if (m_mesh_adaptor.get_vertex_seaming(pVertex) == OUTER)
                m_mesh_adaptor.set_vertex_seaming(pVertex, INNER);
            else
                continue;           // Skip this vertex if it is already done

            // For each neighbor vertex
            typename Adaptor::Vertex_around_vertex_circulator cir, cir_end;
            cir     = m_mesh_adaptor.vertices_around_vertex_begin(pVertex);
            cir_end = cir;
            CGAL_For_all(cir, cir_end)
            {
                // Flag both oriented edges pVertex <-> cir
                m_mesh_adaptor.set_halfedge_seaming(pVertex, cir, INNER);
                m_mesh_adaptor.set_halfedge_seaming(cir, pVertex, INNER);

                // Add surrounding vertices to list without crossing the border
                if (m_mesh_adaptor.get_vertex_seaming(cir) == OUTER)
                    vertices.push_front(cir);
            }
        }
    }

    // Check that the seam is valid, i.e. that the cut mesh is 2-manifold.
    ///
    /// \pre first_seam_vertex -> end_seam_vertex defines the outer seam,
    ///   i.e. Parameterization_mesh_patch_3 will export the <i>right</i> of the seam.
    /// \pre The <i>seam</i> is given as a container of Adaptor::Vertex_handle elements.
    /// \pre The <i>seam</i> is implicitely a loop. The first vertex should *not* be
    ///   duplicated at the end.
    /// \pre The seaming flag of all vertices and edges to INNER, BORDER or OUTER
    ///   w.r.t. the first_seam_vertex -> end_seam_vertex border is set
    ///   (outer seam edges are marked BORDER).
    template<class InputIterator>
    bool check_seam(InputIterator first_seam_vertex,
                    InputIterator end_seam_vertex) const
    {
        // The input vertices list can be either a "seam along a line"
        // (that virtually cut the mesh along a line) or a "cut-out seam"
        // (loop that cuts out a part of the mesh).
        // A "seam along a line" is given as a 2-ways list of vertices.
        InputIterator second_seam_vertex = first_seam_vertex; second_seam_vertex++;
        bool is_seam_along_a_line
            = (m_mesh_adaptor.get_halfedge_seaming(*second_seam_vertex,
                                                   *first_seam_vertex) == BORDER);

        // One cannot mix "seam along a line" and "cut-out seam"
        int seam_length = 0;
        for (InputIterator border_it = first_seam_vertex;
             border_it != end_seam_vertex;
             border_it++)
        {
            seam_length++;  // compute seam length

            // Get next iterator (looping)
            InputIterator next_border_it = border_it;
            next_border_it++;
            if (next_border_it == end_seam_vertex)
                next_border_it = first_seam_vertex;

            // Opposite halfedges are on seam iff this is a "seam along a line"
            if ( is_seam_along_a_line !=
                 (m_mesh_adaptor.get_halfedge_seaming(*next_border_it,
                                                      *border_it) == BORDER) )
            {
                return false;
            }

            // In a "cut-out seam", a vertex cannot belong twice to the seam
            // (see e.g. "8" shape seam)
            if (!is_seam_along_a_line)
            {
                int nb_occurences = count(first_seam_vertex, end_seam_vertex, *border_it);
                if (nb_occurences > 1)
                    return false;
            }
        }

        // A "seam along a line" must be at least 2 edges long (i.e. 4 vertices in list)
        if (is_seam_along_a_line) {
            if (seam_length < 4)
                return false;
        // A "cut-out seam" must contain at least a triangle, thus 3 vertices
        } else {
            if (seam_length < 3)
                return false;
        }

        // else: ok
        return true;
    }

    /// Get facet' seaming status (INNER or OUTER).
    Seaming_status get_facet_seaming(typename Adaptor::Facet_const_handle facet) const
    {
        // do not call is_valid() to avoid an infinite loop
        CGAL_surface_mesh_parameterization_assertion(facet != NULL);

        typename Adaptor::Vertex_around_facet_const_circulator
                            cir = m_mesh_adaptor.facet_vertices_begin(facet);
        CGAL_surface_mesh_parameterization_assertion(cir != NULL);

        bool is_inner_facet = true;
        typename Adaptor::Vertex_around_facet_const_circulator cir_end = cir;
        do
        {
            is_inner_facet &=
                (m_mesh_adaptor.get_vertex_seaming(cir) != OUTER);
            ++cir;
        } while ( cir != cir_end );

        return (is_inner_facet)?
               INNER :
               OUTER;
    }

    /// Get/set vertex seaming flag,
    /// i.e. position of the vertex w.r.t. to the UNIQUE main border.
    Seaming_status  get_vertex_seaming(Vertex_const_handle vertex) const {
        CGAL_surface_mesh_parameterization_assertion(is_valid(vertex));
        return (Seaming_status) m_mesh_adaptor.get_vertex_seaming(
                                                    vertex->vertex());
    }
    void set_vertex_seaming(Vertex_handle vertex, Seaming_status seaming) {
        m_mesh_adaptor.set_vertex_seaming(vertex->vertex(), seaming);
    }

    /// Create a patch vertex from an adaptor vertex + one of its neighbors.
    ///
    /// \pre adaptor_neighbor is a neighbor of adaptor_vertex.
    /// \pre (adaptor_vertex, adaptor_neighbor) must *not* be a seam (non-oriented) edge.
    Vertex_const_handle get_decorated_vertex_from_inner_edge(
                typename Adaptor::Vertex_const_handle adaptor_vertex,
                typename Adaptor::Vertex_const_handle adaptor_neighbor) const
    {
        // We need at least an inner neighbor as input
        CGAL_surface_mesh_parameterization_assertion(
               m_mesh_adaptor.get_halfedge_seaming(adaptor_vertex,
                                                   adaptor_neighbor) != BORDER
            || m_mesh_adaptor.get_halfedge_seaming(adaptor_neighbor,
                                                   adaptor_vertex) != BORDER);

        // if inner vertex
        if (m_mesh_adaptor.get_vertex_seaming(adaptor_vertex) != BORDER)
        {
            // No extra information needed if inner vertex
            return Vertex_const_handle(adaptor_vertex);
        }
        else // if seam vertex
        {
            // find last neighbor (on seam) for a clockwise rotation
            typename Adaptor::Vertex_around_vertex_const_circulator last_cw_neighbor_cir
                = m_mesh_adaptor.vertices_around_vertex_begin(adaptor_vertex,
                                                              adaptor_neighbor);
            while (m_mesh_adaptor.get_halfedge_seaming(last_cw_neighbor_cir,
                                                       adaptor_vertex) != BORDER)
            {
                last_cw_neighbor_cir++;
            }

            // find first clockwise neighbor (on seam) by a counter-clockwise rotation
            typename Adaptor::Vertex_around_vertex_const_circulator first_cw_neighbor_cir
                = m_mesh_adaptor.vertices_around_vertex_begin(adaptor_vertex,
                                                              adaptor_neighbor);
            while (m_mesh_adaptor.get_halfedge_seaming(adaptor_vertex,
                                                       first_cw_neighbor_cir) != BORDER)
            {
                first_cw_neighbor_cir--;
            }

            // The decorated vertex is then:
            return Vertex_const_handle(adaptor_vertex,
                                       last_cw_neighbor_cir,
                                       first_cw_neighbor_cir);
        }
    }
    Vertex_handle get_decorated_vertex_from_inner_edge(
                typename Adaptor::Vertex_handle adaptor_vertex,
                typename Adaptor::Vertex_handle adaptor_neighbor)
    {
        // Call the const version of get_decorated_vertex_from_inner_edge()
        Vertex_const_handle vertex_hdl = get_decorated_vertex_from_inner_edge(
                        (typename Adaptor::Vertex_const_handle)adaptor_vertex,
                        (typename Adaptor::Vertex_const_handle)adaptor_neighbor);
        return const_cast<Vertex*>(&*vertex_hdl);
    }

    /// Create a patch vertex from a border/seam adaptor vertex
    /// + one of its neighbors on the seam.
    ///
    /// \pre adaptor_vertex is a border/seam vertex.
    /// \pre [first_cw_neighbor, last_cw_neighbor] defines the range
    ///   of the valid neighbors of adaptor_vertex (included) or are NULL.
    /// \pre Either first_cw_neighbor or last_cw_neighbor are not NULL.
    Vertex_const_handle get_decorated_vertex_from_border_edge(
                typename Adaptor::Vertex_const_handle adaptor_vertex,
                typename Adaptor::Vertex_const_handle last_cw_neighbor,
                typename Adaptor::Vertex_const_handle first_cw_neighbor) const
    {
        CGAL_surface_mesh_parameterization_assertion(adaptor_vertex != NULL);
        CGAL_surface_mesh_parameterization_assertion(
                      last_cw_neighbor != NULL || first_cw_neighbor != NULL);

        CGAL_surface_mesh_parameterization_assertion(
               last_cw_neighbor == NULL
            || m_mesh_adaptor.get_halfedge_seaming(adaptor_vertex,
                                                   last_cw_neighbor) == BORDER
            || m_mesh_adaptor.get_halfedge_seaming(last_cw_neighbor,
                                                   adaptor_vertex) == BORDER);
        CGAL_surface_mesh_parameterization_assertion(
               first_cw_neighbor == NULL
            || m_mesh_adaptor.get_halfedge_seaming(adaptor_vertex,
                                                   first_cw_neighbor) == BORDER
            || m_mesh_adaptor.get_halfedge_seaming(first_cw_neighbor,
                                                   adaptor_vertex) == BORDER);

        // If both first_cw_neighbor and last_cw_neighbor are provided
        if (last_cw_neighbor != NULL && first_cw_neighbor != NULL)
        {
            // we're done (quick)
            return Vertex_const_handle(adaptor_vertex,
                                       last_cw_neighbor,
                                       first_cw_neighbor);
        }
        else // if either first_cw_neighbor or last_cw_neighbor is missing
        {
            // search in border vertices list (slow)
            for (Border_vertex_const_iterator it = mesh_main_border_vertices_begin();
                 it != mesh_main_border_vertices_end();
                 it++)
            {
                if (it->vertex() == adaptor_vertex
                 && (last_cw_neighbor == NULL  || it->last_cw_neighbor() == last_cw_neighbor)
                 && (first_cw_neighbor == NULL || it->first_cw_neighbor() == first_cw_neighbor))
                {
                    return it;
                }
            }

            // we should not get here
            CGAL_error();
            return NULL;
        }
    }
    Vertex_handle get_decorated_vertex_from_border_edge(
                        typename Adaptor::Vertex_handle adaptor_vertex,
                        typename Adaptor::Vertex_handle last_cw_neighbor,
                        typename Adaptor::Vertex_handle first_cw_neighbor)
    {
        // Call the const version of get_decorated_vertex_from_border_edge()
        Vertex_const_handle vertex_hdl = get_decorated_vertex_from_border_edge(
                        (typename Adaptor::Vertex_const_handle)adaptor_vertex,
                        (typename Adaptor::Vertex_const_handle)last_cw_neighbor,
                        (typename Adaptor::Vertex_const_handle)first_cw_neighbor);
        return const_cast<Vertex*>(&*vertex_hdl);
    }

    /// Debug utility: Check if a Parameterization_mesh_patch_3 facet is valid.
    bool is_valid(Facet_const_handle facet) const
    {
        if (facet == NULL)
            return false;
        // outer facets are not exported
        if (get_facet_seaming(facet) != INNER)
            return false;
        // else: ok
        return true;
    }

    /// Debug utility: Check if a Parameterization_mesh_patch_3 vertex is valid.
    bool is_valid(Vertex_const_handle vertex) const
    {
        if (vertex == NULL)
            return false;
        // outer vertices are not exported
        if (m_mesh_adaptor.get_vertex_seaming(vertex->vertex()) == OUTER)
            return false;
        // prev/next vertices must be on the main border
        if (vertex->last_cw_neighbor() != NULL &&
            m_mesh_adaptor.get_vertex_seaming(vertex->last_cw_neighbor()) != BORDER)
            return false;
        if (vertex->first_cw_neighbor() != NULL &&
            m_mesh_adaptor.get_vertex_seaming(vertex->first_cw_neighbor()) != BORDER)
            return false;
        // else: ok
        return true;
    }

    /// Debug utility: get vertex index as string ("-" if NULL vertex).
    std::string get_vertex_index_as_string(typename Adaptor::Vertex_const_handle vertex) const
    {
        if (vertex == NULL) {
            return std::string("-");
        } else {
            char index_as_string[64];
            std::sprintf(index_as_string, "%d", (int)m_mesh_adaptor.get_vertex_index(vertex));
            return std::string(index_as_string);
        }
    }

// Public fields
public:

    /// The decorated mesh.
    Adaptor& m_mesh_adaptor;

// Private fields
private:

    /// List of all exported vertices.
    /// Order is: inner vertices, then seam/main border ones.
    Param_mesh_patch_vertex_list<Adaptor> m_inner_and_border_vertices;

    /// Iterator to first seam vertex inside m_inner_and_border_vertices.
    Border_vertex_iterator m_seam_begin;

    /// Indicate if the mesh matches the ParameterizationMesh_3 concept.
    bool m_is_valid;

// Private types
private:

    /// Utility class to generate the Facet_iterator type.
    struct Inner_facets_filter
    {
        Inner_facets_filter(const Parameterization_mesh_patch_3& mesh) : m_mesh_patch(mesh) {}

        /// Return true <=> the facet is *not* exported by Parameterization_mesh_patch_3,
        /// i.e. is out of the topological disc.
        bool operator()(const typename Adaptor::Facet_iterator& f) const       {
            return m_mesh_patch.get_facet_seaming(f) == OUTER;
        }
        bool operator()(const typename Adaptor::Facet_const_iterator& f) const {
            return m_mesh_patch.get_facet_seaming(f) == OUTER;
        }

    private:
        const Parameterization_mesh_patch_3& m_mesh_patch;
    };

// Friends
    
    #ifndef DOXYGEN_RUNNING
    friend class Param_mesh_patch_vertex<Adaptor>;
    friend class Param_mesh_patch_vertex_handle<Adaptor>;
    friend class Param_mesh_patch_vertex_const_handle<Adaptor>;
    friend class Param_mesh_patch_vertex_list_iterator<Adaptor>;
    friend class Param_mesh_patch_vertex_list_const_iterator<Adaptor>;
    friend class Mesh_patch_vertex_around_facet_cir<Parameterization_mesh_patch_3,
                                                    Vertex_handle,
                                                    typename Adaptor::Vertex_around_facet_circulator>;
    friend class Mesh_patch_vertex_around_facet_cir<const Parameterization_mesh_patch_3,
                                                    Vertex_const_handle,
                                                    typename Adaptor::Vertex_around_facet_const_circulator>;
    friend class Mesh_patch_vertex_around_vertex_cir<Parameterization_mesh_patch_3,
                                                     Vertex_handle,
                                                     typename Adaptor::Vertex_around_vertex_circulator,
                                                     typename Adaptor::Vertex_handle>;
    friend class Mesh_patch_vertex_around_vertex_cir<const Parameterization_mesh_patch_3,
                                                     Vertex_const_handle,
                                                     typename Adaptor::Vertex_around_vertex_const_circulator,
                                                     typename Adaptor::Vertex_const_handle>;
    #endif
}; // Parameterization_mesh_patch_3


} //namespace CGAL

#endif //CGAL_SURFACE_MESH_PARAMETERIZATION_MESH_PATCH_3_H
