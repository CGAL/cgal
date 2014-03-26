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

#ifndef CGAL_PARAMETERIZATION_MESH_FEATURE_EXTRACTOR_H
#define CGAL_PARAMETERIZATION_MESH_FEATURE_EXTRACTOR_H

#include <list>
#include <vector>
#include <cmath>
#include <CGAL/circulator.h>

/// \file Parameterization_mesh_feature_extractor.h

namespace CGAL {


/// \ingroup  PkgSurfaceParameterizationHelper
///
/// The class Parameterization_mesh_feature_extractor
/// computes features (genus, borders, ...) of a 3D surface,
/// model of the ParameterizationMesh_3 concept.
template<class ParameterizationMesh_3>      //< 3D surface
class Parameterization_mesh_feature_extractor
{
// Public types
public:
    /// Export ParameterizationMesh_3 template parameter.
    typedef ParameterizationMesh_3          Adaptor;

    /// Type representing a border = STL container of vertex handles.
    typedef std::list<typename Adaptor::Vertex_handle>
                                            Border;
    /// Type representing the list of all borders of the mesh
    /// = STL container of Border elements.
    typedef std::vector<Border*>            Skeleton;

// Private types
private:
    // Mesh_Adaptor_3 subtypes:
    typedef typename Adaptor::NT            NT;
    typedef typename Adaptor::Point_2       Point_2;
    typedef typename Adaptor::Point_3       Point_3;
    typedef typename Adaptor::Vector_2      Vector_2;
    typedef typename Adaptor::Vector_3      Vector_3;
    typedef typename Adaptor::Facet         Facet;
    typedef typename Adaptor::Facet_handle  Facet_handle;
    typedef typename Adaptor::Facet_const_handle
                                            Facet_const_handle;
    typedef typename Adaptor::Facet_iterator Facet_iterator;
    typedef typename Adaptor::Facet_const_iterator
                                            Facet_const_iterator;
    typedef typename Adaptor::Vertex        Vertex;
    typedef typename Adaptor::Vertex_handle Vertex_handle;
    typedef typename Adaptor::Vertex_const_handle
                                            Vertex_const_handle;
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

// Public operations
public:

    /// Constructor.
    ///
    /// \attention This class caches the result of feature extractions
    /// => The caller must *not* modify `mesh` during the
    /// Parameterization_mesh_feature_extractor life cycle.
    Parameterization_mesh_feature_extractor(Adaptor& mesh)
        // Store reference to adapted mesh
      : m_mesh_adaptor(mesh)
    {
        // m_mesh_adaptor features are not yet computed
        m_nb_connex_components = -1;
        m_nb_borders = -1;
        m_genus = -1;
    }

    virtual ~Parameterization_mesh_feature_extractor() {
      for (typename Skeleton::iterator iter = m_skeleton.begin();
           iter != m_skeleton.end(); ++iter)
          delete *iter;
    }

    /// Get number of borders.
    int get_nb_borders()
    {
        // At first call, extract borders and put longest one first
        if (m_nb_borders == -1)
            extract_borders();

        return m_nb_borders;
    }
    /// Get extracted borders.
    /// The longest border is the first one.
    const Skeleton& get_borders()
    {
        // At first call, extract borders and put longest one first
        if (m_nb_borders == -1)
            extract_borders();

        return m_skeleton;
    }
    /// Get longest border.
    const Border& get_longest_border()
    {
        // At first call, extract borders and put longest one first
        if (m_nb_borders == -1)
            extract_borders();

        return *(m_skeleton[0]);
    }

    /// Get # of connected components.
    int get_nb_connex_components()
    {
        // At first call, count the number of connex components
        if (m_nb_connex_components == -1)
            count_connex_components();

        return m_nb_connex_components;
    }

    /// Get the genus.
    int get_genus()
    {
        // At first call, compute the genus
        if (m_genus == -1)
            compute_genus();

        return m_genus;
    }

// Private operations
private:

    /// Extract borders and put longest one first.
    /// Result is in m_nb_borders and m_skeleton.
    void extract_borders()
    {
        CGAL_surface_mesh_parameterization_precondition(m_skeleton.size() == 0);

        // Tag all vertices as unprocessed
        const int tag_free = 0;
        const int tag_done = 1;
        for (Vertex_iterator it = m_mesh_adaptor.mesh_vertices_begin();
             it != m_mesh_adaptor.mesh_vertices_end();
             it++)
        {
             m_mesh_adaptor.set_vertex_tag(it, tag_free);
        }

        // find all closed borders
        while (add_border(tag_free,tag_done)) {}

        // #borders
        m_nb_borders = static_cast<int>(m_skeleton.size());

        // put longest border first
        if (m_nb_borders>1)
        {
            int index = get_index_longest_border();
            Border *tmp = m_skeleton[0];
            m_skeleton[0] = m_skeleton[index];
            m_skeleton[index] = tmp;
        }
    }

    /// Add closed border.
    bool add_border(int tag_free, int tag_done)
    {
        // Find a border tagged as "free" and tag it as "processed"
        // Return an empty list if not found
        std::list<Vertex_handle> border = find_free_border(tag_free, tag_done);
        if(border.empty())
            return false;

        // add one border to list
        Border *pNewBorder = new Border;
        *pNewBorder = border;
        m_skeleton.push_back(pNewBorder);

        return true;
    }

    /// Find a border tagged as <i>free</i> and tag it as <i>processed</i>.
    /// Return an empty list if not found.
    std::list<Vertex_handle> find_free_border(int tag_free, int tag_done)
    {
        std::list<Vertex_handle> border;    // returned list

        // get any border vertex with "free" tag
        Vertex_handle seed_vertex = NULL;
        for (Vertex_iterator pVertex = m_mesh_adaptor.mesh_vertices_begin();
             pVertex != m_mesh_adaptor.mesh_vertices_end();
             pVertex++)
        {
            if (m_mesh_adaptor.is_vertex_on_border(pVertex) &&
                m_mesh_adaptor.get_vertex_tag(pVertex) == tag_free)
            {
                seed_vertex = pVertex;
                break;
            }
        }
        if (seed_vertex == NULL)
            return border;                  // return empty list

        // Get the border containing seed_vertex
        border = m_mesh_adaptor.get_border(seed_vertex);

        // Tag border vertices as "processed"
        typename std::list<Vertex_handle>::iterator it;
        for(it = border.begin(); it != border.end(); it++)
            m_mesh_adaptor.set_vertex_tag(*it, tag_done);

        return border;
    }

    /// Get index of the longest border.
    int get_index_longest_border() const
    {
        int index = 0;
        double max = 0.0;

        // #borders
        int nb = static_cast<int>(m_skeleton.size());

        for(int i=0;i<nb;i++)
        {
            const Border *pBorder = m_skeleton[i];
            double length = len(*pBorder);
            if (length > max)
            {
                index = i;
                max = length;
            }
        }

        return index;
    }

    /// Compute  total len of a border.
    double len(const Border& border) const
    {
        double len = 0.0;
        typename std::list<typename Adaptor::Vertex_handle>::const_iterator it;
        for(it = border.begin(); it != border.end(); it++)
        {
            // Get next iterator (looping)
            typename std::list<typename Adaptor::Vertex_handle>::const_iterator next = it;
            next++;
            if (next == border.end())
                next = border.begin();

            Vector_3 v = m_mesh_adaptor.get_vertex_position(*next)
                       - m_mesh_adaptor.get_vertex_position(*it);
            len += std::sqrt(v*v);
        }
        return len;
    }

    /// Count # of connected components.
    /// Result is in m_nb_connex_components.
    void count_connex_components()
    {
        m_nb_connex_components = 0;

        const int tag_free = 0;
        const int tag_done = 1;
        for (Vertex_iterator it = m_mesh_adaptor.mesh_vertices_begin();
             it != m_mesh_adaptor.mesh_vertices_end();
             it++)
        {
             m_mesh_adaptor.set_vertex_tag(it, tag_free);
        }

        Vertex_handle seed_vertex = NULL;
        while((seed_vertex = get_any_vertex_tag(tag_free)) != NULL)
        {
            m_nb_connex_components++;
            tag_component(seed_vertex, tag_free, tag_done);
        }
    }

    /// Get any vertex with tag.
    Vertex_handle get_any_vertex_tag(int tag)
    {
        for (Vertex_iterator it = m_mesh_adaptor.mesh_vertices_begin();
             it != m_mesh_adaptor.mesh_vertices_end();
             it++)
        {
            if (m_mesh_adaptor.get_vertex_tag(it) == tag)
                return it;
        }

        return NULL;
    }

    /// Tag component.
    void tag_component(Vertex_handle pSeedVertex,
                       const int tag_free,
                       const int tag_done)
    {
        CGAL_surface_mesh_parameterization_precondition(
             m_mesh_adaptor.get_vertex_tag(pSeedVertex) == tag_free);

        std::list<Vertex_handle> vertices;
        vertices.push_front(pSeedVertex);

        while (!vertices.empty())
        {
            Vertex_handle pVertex = vertices.front();
            vertices.pop_front();

            // Stop if already done
            if (m_mesh_adaptor.get_vertex_tag(pVertex) == tag_done)
                continue;

            m_mesh_adaptor.set_vertex_tag(pVertex, tag_done);

            Vertex_around_vertex_circulator cir, cir_end;
            cir     = m_mesh_adaptor.vertices_around_vertex_begin(pVertex);
            cir_end = cir;
            CGAL_For_all(cir,cir_end)
                if (m_mesh_adaptor.get_vertex_tag(cir) == tag_free)
                    vertices.push_front(cir);
        }
    }

    /// Compute the genus. Result is in m_genus.
    ///
    /// Implementation note:
    ///  G = (2*C + E - B - F - V)/2 with
    ///  G : genus
    ///  C : # of connected components
    ///  E : # of edges
    ///  B : # of borders
    ///  F : # of facets
    ///  V : # of vertices
    void compute_genus()
    {
        int c = get_nb_connex_components();
        int b = get_nb_borders();
        int v = m_mesh_adaptor.count_mesh_vertices();
        int e = m_mesh_adaptor.count_mesh_halfedges()/2;
        int f = m_mesh_adaptor.count_mesh_facets();

        m_genus = (2*c+e-b-f-v)/2;
    }

// Fields
private:

    /// Pointer to mesh to parse
    Adaptor&    m_mesh_adaptor;

    /// m_mesh_adaptor features:
    int         m_nb_borders;
    Skeleton    m_skeleton;             ///< List of borders of m_mesh_adaptor
    int         m_nb_connex_components;
    int         m_genus;

}; // Parameterization_mesh_feature_extractor


} //namespace CGAL

#endif // CGAL_SURFACE_MESH_PARAMETERIZATION_MESH_FEATURE_EXTRACTOR_H
