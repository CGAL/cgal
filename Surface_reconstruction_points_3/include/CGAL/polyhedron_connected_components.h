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
#include <CGAL/point_set_processing_assertions.h>

#include <map>
#include <list>

CGAL_BEGIN_NAMESPACE


// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
namespace CGALi {


/// Possible values of a vertex tag.
enum { tag_free, tag_done };


/// Get any vertex with tag == tag_free.
///
/// @commentheading Template Parameters:
/// @param Polyhedron an instance of CGAL::Polyhedron_3<Traits>.
///
/// @return a list of pairs (component's size (number of vertices), a vertex of the component),
/// ordered by size.

template<class Polyhedron>
typename Polyhedron::Vertex_handle
get_any_free_vertex(Polyhedron& polyhedron,
                    std::map<typename Polyhedron::Vertex*, int>& tags)
{
    typedef typename Polyhedron::Vertex_handle Vertex_handle;
    typedef typename Polyhedron::Vertex_iterator Vertex_iterator;

    for (Vertex_iterator it = polyhedron.vertices_begin();
         it != polyhedron.vertices_end();
         it++)
    {
        if (tags[&*it] == tag_free)
            return it;
    }

    return NULL;
}

/// Tag a "free" connected component as "done".
///
/// @commentheading Template Parameters:
/// @param Polyhedron an instance of CGAL::Polyhedron_3<Traits>.
///
/// @return the size (number of vertices) of the component.
template<class Polyhedron>
unsigned int tag_component(Polyhedron& polyhedron,
                           typename Polyhedron::Vertex_handle pSeedVertex,
                           std::map<typename Polyhedron::Vertex*, int>& tags)
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

        // Add vertex's "free" neighbors to the list
        Halfedge_around_vertex_circulator neighbor_cir, neighbor_end;
        neighbor_cir = pVertex->vertex_begin();
        neighbor_end = neighbor_cir;
        CGAL_For_all(neighbor_cir,neighbor_end)
        {
            Vertex_handle neighbor = neighbor_cir->opposite()->vertex();
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
/// @param Polyhedron an instance of CGAL::Polyhedron_3<Traits> that supports vertices.
///
/// @return a list of components expressed as pairs (number of vertices, vertex),
/// ordered by size.
template<class Polyhedron>
std::multimap<unsigned int, typename Polyhedron::Vertex_handle>
get_polyhedron_connected_components(Polyhedron& polyhedron)
{
    // Implementation note:
    // We tag vertices instead of halfedges to save a factor 6.
    // The drawback is that we require the Polyhedron_3<Traits> to support vertices.
    // TODO: replace std::map to tag vertices by a property map.
    Assert_compile_time_tag(typename Polyhedron::Supports_halfedge_vertex(), Tag_true());
    std::map<typename Polyhedron::Vertex*, int> tags;

    typedef typename Polyhedron::Vertex_handle Vertex_handle;
    typedef typename Polyhedron::Vertex_iterator Vertex_iterator;

    // list of all connected components of a polyhedron, ordered by size.
    std::multimap<unsigned int, Vertex_handle> components;

    // Tag all mesh vertices as "free".
    for (Vertex_iterator it = polyhedron.vertices_begin();
         it != polyhedron.vertices_end();
         it++)
    {
         tags[&*it] = CGALi::tag_free;
    }

    // For each component
    Vertex_handle seed_vertex = NULL;
    while((seed_vertex = CGALi::get_any_free_vertex(polyhedron, tags)) != NULL)
    {
        // Tag it as "done" and compute its size (number of vertices)
        unsigned int number_of_vertices = CGALi::tag_component(polyhedron, seed_vertex, tags);

        // Add component to (ordered) list
          components.insert(std::make_pair(number_of_vertices, seed_vertex));
    }

    return components;
}


/// Erase small connected components of a polyhedron:
/// erase all connected components but the largest.
///
/// @commentheading Template Parameters:
/// @param Polyhedron an instance of CGAL::Polyhedron_3<Traits> that supports
/// vertices and removal operation.
///
/// @return the number of connected components erased.
template<class Polyhedron>
unsigned int
erase_small_polyhedron_connected_components(Polyhedron& polyhedron)
{
    Assert_compile_time_tag(typename Polyhedron::Supports_removal(), Tag_true());

    typedef typename Polyhedron::Vertex_handle Vertex_handle;

    unsigned int nb_erased_components = 0,
                 nb_isolated_vertices = 0;

    // Get list of connected components, ordered by size (number of vertices)
    std::multimap<unsigned int, Vertex_handle>
      components = CGAL::get_polyhedron_connected_components(polyhedron);

    // Erase all connected components but the largest
    while (components.size() > 1)
    {
      unsigned int number_of_vertices = components.begin()->first;
      Vertex_handle vertex = components.begin()->second;

      // Remove component from list
      components.erase(components.begin());

      if (vertex->halfedge() != NULL) // if not isolated vertex
      {
        CGAL_TRACE_STREAM << "  Erase connected component (" << number_of_vertices << " vertices)\n";
        polyhedron.erase_connected_component(vertex->halfedge());
        nb_erased_components++;
      }
      else // if isolated vertex
      {
        // TODO: erase isolated vertices?
        // Note: Polyhedron_3 does not export HalfedgeDS::vertices_erase(Vertex_handle v)

        nb_isolated_vertices++;
      }
    }

    if (nb_isolated_vertices > 0)
      CGAL_TRACE_STREAM << "  Skipped " << nb_isolated_vertices << " isolated vertices\n";

    return nb_erased_components;
}


CGAL_END_NAMESPACE

#endif // CGAL_POLYHEDRON_CONNECTED_COMPONENTS_H
