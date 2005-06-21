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

#ifndef MESH_ADAPTOR_FEATURE_EXTRACTOR_H
#define MESH_ADAPTOR_FEATURE_EXTRACTOR_H

#include <CGAL/basic.h>

#include <list>
#include <vector>

CGAL_BEGIN_NAMESPACE


// Class Mesh_adaptor_feature_extractor
// 
// This class computes features (genus, boundaries, ...)
// of a 3D surface model of the MeshAdaptor_3 concept.

template<class MeshAdaptor_3>           // 3D surface
class Mesh_adaptor_feature_extractor
{
// Public types
public:

    // Export Mesh_Adaptor_3 type and subtypes
    typedef MeshAdaptor_3                   Adaptor;
    typedef typename Adaptor::NT            NT;
    typedef typename Adaptor::Facet_handle  Facet_handle;
    typedef typename Adaptor::Facet_const_handle
                                            Facet_const_handle;
    typedef typename Adaptor::Vertex_handle Vertex_handle;
    typedef typename Adaptor::Vertex_const_handle
                                            Vertex_const_handle;
    typedef typename Adaptor::Point_3       Point_3;
    typedef typename Adaptor::Point_2       Point_2;
    typedef typename Adaptor::Vector_3      Vector_3;
    typedef typename Adaptor::Vector_2      Vector_2;
    typedef typename Adaptor::Facet_iterator Facet_iterator;
    typedef typename Adaptor::Facet_const_iterator
                                            Facet_const_iterator;
    typedef typename Adaptor::Vertex_iterator Vertex_iterator;
    typedef typename Adaptor::Vertex_const_iterator
                                            Vertex_const_iterator;
    typedef typename Adaptor::Border_vertex_iterator
                                            Border_vertex_iterator;
    typedef typename Adaptor::Border_vertex_const_iterator
                                            Border_vertex_const_iterator;
    typedef typename Adaptor::Vertex_around_facet_circulator
                                            Vertex_around_facet_circulator;
    typedef typename Adaptor::Vertex_around_facet_const_circulator
                                            Vertex_around_facet_const_circulator;
    typedef typename Adaptor::Vertex_around_vertex_circulator
                                            Vertex_around_vertex_circulator;
    typedef typename Adaptor::Vertex_around_vertex_const_circulator
                                            Vertex_around_vertex_const_circulator;

    // Mesh boundary
    typedef std::list<Vertex_handle>        Boundary;
    // List of all boundaries of a mesh
    typedef std::vector<Boundary*>          Skeleton;

// Public operations
public:

    // Constructor
    // 
    // CAUTION: Caller must NOT modify 'mesh' during the
    //          Mesh_adaptor_feature_extractor life cycle.
    Mesh_adaptor_feature_extractor(Adaptor *mesh)
    {
        m_mesh_adaptor = mesh;
        CGAL_parameterization_assertion(m_mesh_adaptor != NULL);

        // m_mesh_adaptor features are not yet computed
        m_nb_connex_components = -1;
        m_nb_boundaries = -1;
        m_genus = -1;
    }
    virtual ~Mesh_adaptor_feature_extractor() {}

    // Get number of boundaries 
    int get_nb_boundaries() 
    { 
        // At first call, extract boundaries and put longest one first
        if (m_nb_boundaries == -1)
            extract_boundaries();

        return m_nb_boundaries; 
    }
    // Get extracted boundaries
    // The longest boundary is the first one
    const Skeleton& get_boundaries() 
    { 
        // At first call, extract boundaries and put longest one first
        if (m_nb_boundaries == -1)
            extract_boundaries();

        return m_skeleton; 
    }
    // Get longest boundary 
    const Boundary* get_longest_boundary() 
    { 
        // At first call, extract boundaries and put longest one first
        if (m_nb_boundaries == -1)
            extract_boundaries();

        return m_skeleton[0]; 
    }

    // Get # of connected components
    int get_nb_connex_components() 
    { 
        // At first call, count the number of connex components
        if (m_nb_connex_components == -1)
            count_connex_components();

        return m_nb_connex_components;
    }

    // Get the genus
    int get_genus() 
    { 
        // At first call, compute the genus
        if (m_genus == -1)
            compute_genus();

        return m_genus;
    }

// Private operations
private:

    // Extract boundaries and put longest one first
    // Result is in m_nb_boundaries and m_skeleton
    void extract_boundaries()
    {
        assert(m_skeleton.size() == 0);

        // Tag all vertices as unprocessed
        const int tag_free = 0;
        const int tag_done = 1;
        for (Vertex_iterator it = m_mesh_adaptor->mesh_vertices_begin(); 
             it != m_mesh_adaptor->mesh_vertices_end(); 
             it++)
        {
             m_mesh_adaptor->set_vertex_tag(it, tag_free);
        }

        // find all closed boundaries
        while (add_boundary(tag_free,tag_done)) {}

        // #boundaries
        m_nb_boundaries = m_skeleton.size();

        // put longest boundary first if required
        if (m_nb_boundaries>1)
        {
            int index = get_index_longest_boundary();
            Boundary *tmp = m_skeleton[0];
            m_skeleton[0] = m_skeleton[index];
            m_skeleton[index] = tmp;
        }

        std::cerr << "  " << m_nb_boundaries << " boundary(ies) found" << std::endl;
    }

    // add closed boundary Boundary
    bool add_boundary(int tag_free, int tag_done)
    {
        // Find a boundary tagged as "free" and tag it as "processed"
        // Return an empty list if not found
        std::list<Vertex_handle> boundary = find_free_boundary(tag_free, tag_done);
        if(boundary.empty())
            return false;

        // add one boundary to list
        Boundary *pNewBoundary = new Boundary;
        *pNewBoundary = boundary;
        m_skeleton.push_back(pNewBoundary);

        return true;
    }

    // Find a boundary tagged as "free" and tag it as "processed"
    // Return an empty list if not found
    std::list<Vertex_handle> find_free_boundary(int tag_free, int tag_done)
    {
        std::list<Vertex_handle> boundary;  // returned list

        // get any border vertex with "free" tag
        Vertex_handle seed_vertex = NULL;
        for (Vertex_iterator pVertex = m_mesh_adaptor->mesh_vertices_begin();
             pVertex != m_mesh_adaptor->mesh_vertices_end();
             pVertex++)
        {
            if (m_mesh_adaptor->is_vertex_on_border(pVertex) && 
                m_mesh_adaptor->get_vertex_tag(pVertex) == tag_free) 
            {
                seed_vertex = pVertex;
                break;
            }
        }
        if (seed_vertex == NULL)
            return boundary;                // return empty list

        // Get the boundary containing seed_vertex
        boundary = m_mesh_adaptor->get_boundary(seed_vertex);

        // Tag boundary vertices as "processed"
        std::list<Vertex_handle>::iterator it;
        for(it = boundary.begin(); it != boundary.end(); it++)
            m_mesh_adaptor->set_vertex_tag(*it, tag_done);

        return boundary;
    }

    // get index of the longest boundary
    int get_index_longest_boundary() const
    {
        int index = 0;
        double max = 0.0;

        // #boundaries
        int nb = m_skeleton.size();

        for(int i=0;i<nb;i++)
        {
            const Boundary *pBoundary = m_skeleton[i];
            double length = len(pBoundary);
            if (length > max)
            {
                index = i;
                max = length;
            }
        }

        return index;
    }

    // compute  total len of a boundary
    double len(const Boundary* pBoundary) const
    {
        double len = 0.0;
        std::list<Adaptor::Vertex_handle>::const_iterator it;
        for(it = pBoundary->begin(); it != pBoundary->end(); it++)
        {
            // Get next iterator (looping)
            std::list<Adaptor::Vertex_handle>::const_iterator next = it;
            next++;
            if (next == pBoundary->end())
                next = pBoundary->begin();

            Vector_3 v = m_mesh_adaptor->get_vertex_position(*next)
                    - m_mesh_adaptor->get_vertex_position(*it);
            len += std::sqrt(v*v);
        }
        return len;
    }

    // Count # of connected components
    // Result is in m_nb_connex_components
    void count_connex_components()
    {
        m_nb_connex_components = 0;

        const int tag_free = 0;
        const int tag_done = 1;
        for (Vertex_iterator it = m_mesh_adaptor->mesh_vertices_begin(); 
             it != m_mesh_adaptor->mesh_vertices_end(); 
             it++)
        {
             m_mesh_adaptor->set_vertex_tag(it, tag_free);
        }

        Vertex_handle seed_vertex = NULL;
        while((seed_vertex = get_any_vertex_tag(tag_free)) != NULL)
        {
            m_nb_connex_components++;
            tag_component(seed_vertex, tag_free, tag_done);
        }
    }

    // get any vertex with tag
    Vertex_handle get_any_vertex_tag(int tag)
    {
        for (Vertex_iterator it = m_mesh_adaptor->mesh_vertices_begin(); 
             it != m_mesh_adaptor->mesh_vertices_end(); 
             it++)
        {
            if (m_mesh_adaptor->get_vertex_tag(it) == tag)
                return it;
        }

        return NULL;
    }

    // tag component
    void tag_component(Vertex_handle pSeedVertex,
                       const int tag_free,
                       const int tag_done)
    {
        assert(m_mesh_adaptor->get_vertex_tag(pSeedVertex) == tag_free);

        std::list<Vertex_handle> vertices;
        vertices.push_front(pSeedVertex);

        while (!vertices.empty())
        {
            Vertex_handle pVertex = vertices.front();
            vertices.pop_front();

            // Stop if already done
            if (m_mesh_adaptor->get_vertex_tag(pVertex) == tag_done)
                continue;

            m_mesh_adaptor->set_vertex_tag(pVertex, tag_done);

            Vertex_around_vertex_circulator cir, cir_end;
            cir     = m_mesh_adaptor->vertices_around_vertex_begin(pVertex);
            cir_end = cir;
            CGAL_For_all(cir,cir_end)
                if (m_mesh_adaptor->get_vertex_tag(cir) == tag_free)
                    vertices.push_front(cir);
        }
    }

    // Compute the genus
    // Result is in m_genus
    //
    // Implementation note:
    //  G = (2*C + E - B - F - V)/2 with
    //  G : genus
    //  C : # of connected components
    //  E : # of edges
    //  B : # of boundaries
    //  F : # of facets
    //  V : # of vertices
    void compute_genus()
    {
        int c = get_nb_connex_components();
        int b = get_nb_boundaries();
        int v = m_mesh_adaptor->count_mesh_vertices();
        int e = m_mesh_adaptor->count_mesh_halfedges()/2;
        int f = m_mesh_adaptor->count_mesh_facets();
        
        m_genus = (2*c+e-b-f-v)/2;
        std::cerr << "  " << v << " vertices, " << f << " facets, ";
        std::cerr << e << " edges, " << b << " boundary(ies), genus " << m_genus << std::endl;
    }

// Fields
private:

    // Pointer to mesh to parse
    Adaptor*    m_mesh_adaptor;

    // m_mesh_adaptor features:
    int         m_nb_boundaries;
    Skeleton    m_skeleton;             // List of boundaries of m_mesh_adaptor
    int         m_nb_connex_components;
    int         m_genus;

}; // Mesh_adaptor_feature_extractor


CGAL_END_NAMESPACE

#endif // MESH_ADAPTOR_FEATURE_EXTRACTOR_H
