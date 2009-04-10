// Copyright (c) 2005-2009  INRIA (France).
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
// $URL$
// $Id$
//
//
// Author(s)     : Laurent Saboret, Pierre Alliez

#ifndef CGAL_POLYHEDRON_CONNECTED_COMPONENTS_H
#define CGAL_POLYHEDRON_CONNECTED_COMPONENTS_H

#include <CGAL/Polyhedron_3.h>

#include <map>
#include <list>

CGAL_BEGIN_NAMESPACE


// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
namespace CGALi {


/// Get any vertex with tag == tag_value.
///
/// @commentheading Template Parameters:
/// @param Polyhedron an instance of CGAL::Polyhedron_3<>.
///
/// @return a list of pairs (component size, a vertex of the component), 
/// ordered by size.

template<class Polyhedron> 
typename Polyhedron::Vertex_handle 
get_any_vertex_tag(Polyhedron& polyhedron, 
                   std::map<typename Polyhedron::Vertex*, int>& tags, 
                   const int tag_value)
{
    typedef typename Polyhedron::Vertex_handle Vertex_handle;
    typedef typename Polyhedron::Vertex_iterator Vertex_iterator;

    for (Vertex_iterator it = polyhedron.vertices_begin();
         it != polyhedron.vertices_end();
         it++)
    {
        if (tags[&*it] == tag_value)
            return it;
    }

    return NULL;
}

/// Tag a "free" connected component as "done". 
///
/// @commentheading Template Parameters:
/// @param Polyhedron an instance of CGAL::Polyhedron_3<>.
///
/// @return the size (number of vertices) of the component.
template<class Polyhedron> 
unsigned int tag_component(Polyhedron& polyhedron, 
                           typename Polyhedron::Vertex_handle pSeedVertex,
                           std::map<typename Polyhedron::Vertex*, int>& tags, 
                           const int tag_free,
                           const int tag_done)
{
    typedef typename Polyhedron::Vertex_handle Vertex_handle;
    typedef typename Polyhedron::Vertex_iterator Vertex_iterator;
    typedef typename Polyhedron::Halfedge_around_vertex_circulator 
      Halfedge_around_vertex_circulator;

    unsigned int number_of_vertices = 0; // size (number of vertices) of the component
    
    std::list<Vertex_handle> vertices;
    vertices.push_front(pSeedVertex);
    while (!vertices.empty())
    {
        Vertex_handle pVertex = vertices.front();
        vertices.pop_front();

        // Skip vertex if already done
        if (tags[&*pVertex] == tag_done)
            continue;

        // Mark vertex done
        tags[&*pVertex] = tag_done;
        number_of_vertices++;

        // Add free neighbors to the list
        Halfedge_around_vertex_circulator cir, cir_end;
        cir     = pVertex->vertex_begin();
        cir_end = cir;
        CGAL_For_all(cir,cir_end)
        {
            Vertex_handle neighbor = cir->opposite()->vertex();
            if (tags[&*neighbor] == tag_free)
                vertices.push_front(neighbor);
        }
    }
    
    return number_of_vertices;
}


} /* namespace CGALi */


// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------


/// Compute the list of all connected components of a polyhedron.
///
/// @commentheading Template Parameters:
/// @param Polyhedron an instance of CGAL::Polyhedron_3<> that supports positions.
///
/// @return a list of pairs (component size, a halfedge of the component), 
/// ordered by size.

template<class Polyhedron> 
std::multimap<unsigned int, typename Polyhedron::Halfedge_handle>
get_polyhedron_connected_components(Polyhedron& polyhedron)
{
    // Implementation note: 
    // We tag vertices instead of halfedges to save a factor 6.
    // The drawback is that we require the Polyhedron_3<> to support positions.
    // TODO: replace std::map to tag vertices by a property map.
    Assert_compile_time_tag(Polyhedron::Supports_halfedge_vertex(), Tag_true());
    std::map<typename Polyhedron::Vertex*, int> tags;

    typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
    typedef typename Polyhedron::Vertex_handle Vertex_handle;
    typedef typename Polyhedron::Vertex_iterator Vertex_iterator;

    /// list of all connected components of a polyhedron, ordered by size.
    std::multimap<unsigned int, Halfedge_handle> components;

    // Tag all mesh vertices as "free".
    const int tag_free = 0;
    const int tag_done = 1;
    for (Vertex_iterator it = polyhedron.vertices_begin();
         it != polyhedron.vertices_end();
         it++)
    {
         tags[&*it] = tag_free;
    }

    // For each component
    Vertex_handle seed_vertex = NULL;
    while((seed_vertex = CGALi::get_any_vertex_tag(polyhedron, tags, tag_free)) != NULL)
    {
        // Tag it as "done" and compute its size
        unsigned int component_size = CGALi::tag_component(polyhedron, seed_vertex, tags, tag_free, tag_done);

        // Add component to list, ordered by size
        components.insert(std::make_pair(component_size, seed_vertex->halfedge()));
    }
    
    return components;
}


CGAL_END_NAMESPACE

#endif // CGAL_POLYHEDRON_CONNECTED_COMPONENTS_H
