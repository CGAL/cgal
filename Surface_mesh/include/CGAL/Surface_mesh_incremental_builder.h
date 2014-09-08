// Copyright (c) 1997  ETH Zurich (Switzerland).
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
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>)

#ifndef CGAL_SURFACE_MESH_INCREMENTAL_BUILDER_H
#define CGAL_SURFACE_MESH_INCREMENTAL_BUILDER_H 1

#include <CGAL/basic.h>
#include <CGAL/Random_access_adaptor.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/IO/Verbose_ostream.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

#include <vector>
#include <cstddef>

namespace CGAL {

template < class HalfedgeDS_>
class Surface_mesh_incremental_builder {
public:
    typedef HalfedgeDS_                     HDS; // internal
    typedef HalfedgeDS_                     HalfedgeDS;
    typedef typename HDS::Vertex_descriptor     Vertex_descriptor;
    typedef typename HDS::Halfedge_descriptor   Halfedge_descriptor;
    typedef typename HDS::Face_descriptor       Face_descriptor;

    typedef typename HDS::Point          Point_3;
    typedef typename HDS::size_type         size_type;

protected:
    typedef typename HDS::Vertex_iterator           Vertex_iterator;
    typedef typename HDS::Halfedge_iterator         Halfedge_iterator;
  typedef Random_access_value_adaptor<Vertex_iterator,Vertex_descriptor>  Random_access_index;

    bool                      m_error;
    bool                      m_verbose;
    HDS&                      hds;
    size_type                 rollback_v;
    size_type                 rollback_f;
    size_type                 rollback_h;
    size_type                 new_vertices;
    size_type                 new_faces;
    size_type                 new_halfedges;
    Face_descriptor               current_face;
    Random_access_index       index_to_vertex_map;
    std::vector< Halfedge_descriptor>  vertex_to_edge_map;

    Halfedge_descriptor           g1;      // first halfedge, 0 denotes none.
    Halfedge_descriptor           gprime;
    Halfedge_descriptor           h1;      // current halfedge
    size_type                 w1;      // first vertex.
    size_type                 w2;      // second vertex.
    size_type                 v1;      // current vertex
    bool                      first_vertex;
    bool                      last_vertex;

    CGAL_assertion_code( int check_protocol;)  // use to check protocol.
    // states for checking: 0 = created, 1 = constructing, 2 = make face.

    // Implement the vertex_to_edge_map either with 
    // the halfedge pointer in the vertices
    // ----------------------------------------------------
    void initialize_vertex_to_edge_map( size_type  n, bool mode) {
        vertex_to_edge_map.clear();
        vertex_to_edge_map.resize(n);
        if ( mode) {
            // go through all halfedges and keep a halfedge for each
            // vertex found in a hashmap.
            size_type i = 0;
            for ( Vertex_iterator vi = hds.vertices_begin();
                  vi != hds.vertices_end();
                  ++vi) {
              set_vertex_to_edge_map( i, halfedge(*vi,hds));
                ++i;
            }
        }
    }
  
    void push_back_vertex_to_edge_map( Halfedge_descriptor h) {
        vertex_to_edge_map.push_back(h);
    }
    Halfedge_descriptor get_vertex_to_edge_map( size_type i) {
      // Use the halfedge pointer within the vertex.
      return halfedge(index_to_vertex_map[i], hds);
    }

    void set_vertex_to_edge_map( size_type i, Halfedge_descriptor h) {
      // Use the self-managed array vertex_to_edge_map.
        CGAL_assertion(i < vertex_to_edge_map.size());
        vertex_to_edge_map[i] = h;
        // Use the halfedge pointer within the vertex.
        set_halfedge(index_to_vertex_map[i], h, hds);
    }

// An Incremental Builder for Polyhedral Surfaces
// ----------------------------------------------
// DEFINITION
//
// Surface_mesh_incremental_builder<HDS> is an auxiliary class that
// supports the incremental construction of polyhedral surfaces. This is
// for example convinient when constructing polyhedral surfaces from
// files. The incremental construction starts with a list of all point
// coordinates and concludes with a list of all facet polygons. Edges are
// not explicitly specified. They are derived from the incidence
// information provided from the facet polygons. These are given as a
// sequence of vertex indices. The correct protocol of method calls to
// build a polyhedral surface can be stated as regular expression:
//
// `begin_surface (add_vertex | (begin_facet add_vertex_to_facet*
//  end_facet))* end_surface '
//
// PARAMETERS
//
// `HDS' is the halfedge data structure used to represent the
// polyhedral surface that is to be constructed.
//
// CREATION
public:
    bool error() const { return m_error; }

    Surface_mesh_incremental_builder( HDS& h, bool verbose = false)
        // stores a reference to the halfedge data structure `h' in the
        // internal state. The previous polyhedral surface in `h'
        // remains unchanged. The incremental builder adds the new
        // polyhedral surface to the old one.
      : m_error( false), m_verbose( verbose), hds(h) {
        CGAL_assertion_code(check_protocol = 0;)
    }

    ~Surface_mesh_incremental_builder() {
        CGAL_assertion( check_protocol == 0);
    }

// OPERATIONS
    enum { RELATIVE_INDEXING = 0, ABSOLUTE_INDEXING = 1};


    void begin_surface( std::size_t v, std::size_t f, std::size_t h = 0,
                        int mode = RELATIVE_INDEXING);
        // starts the construction. v is the number of new
        // vertices to expect, f the number of new facets, and h the number of
        // new halfedges. If h is unspecified (`== 0') it is estimated using
        // Euler equations (plus 5% for the so far unkown holes and genus
        // of the object). These values are used to reserve space in the
        // surface_mesh representation `HDS'. If the representation
        // supports insertion these values do not restrict the class of
        // readable surface_meshs. If the representation does not support
        // insertion the object must fit in the reserved sizes.
        //    If `mode' is set to ABSOLUTE_INDEXING the incremental builder
        // uses absolute indexing and the vertices of the old polyhedral 
        // surface can be used in new facets. Otherwise relative indexing is 
        // used starting with new indices for the new construction.


    Vertex_descriptor add_vertex( const Point_3& p) {
        // adds p to the vertex list.
        CGAL_assertion( check_protocol == 1);
        Vertex_descriptor v = CGAL::add_vertex(hds);
        set_halfedge( v, Halfedge_descriptor(),hds);
        put(vertex_point, hds, v, p);
        index_to_vertex_map.push_back(Vertex_iterator(v, &hds));
        push_back_vertex_to_edge_map( Halfedge_descriptor());
        ++new_vertices;
        return v;
    }

  
    // returns handle for the vertex of index i
    Vertex_descriptor vertex( std::size_t i) {
        if ( i < new_vertices)
            return index_to_vertex_map[i];
        return Vertex_descriptor();
    }
  
    Face_descriptor begin_facet() {
        // starts a facet.
        if ( m_error)
            return Face_descriptor();
        CGAL_assertion( check_protocol == 1);
        CGAL_assertion_code( check_protocol = 2;)
       
        // initialize all status variables.
        first_vertex = true;  // denotes 'no vertex yet'
        g1 =  Halfedge_descriptor();  // denotes 'no halfedge yet'
        last_vertex = false;

        return add_face(hds);
    }

    void add_vertex_to_facet( std::size_t i);
        // adds a vertex with index i to the current facet. The first
        // point added with `add_vertex()' has the index 0.

    Halfedge_descriptor end_facet() {
        // ends a facet.
        if ( m_error)
            return Halfedge_descriptor();
        CGAL_assertion( check_protocol == 2);
        CGAL_assertion( ! first_vertex);
        // cleanup all static status variables
        add_vertex_to_facet( w1);
        if ( m_error)
            return Halfedge_descriptor();
        last_vertex = true;
        add_vertex_to_facet( w2);
        if ( m_error)
            return Halfedge_descriptor();
        CGAL_assertion( check_protocol == 2);
        CGAL_assertion_code( check_protocol = 1;)
        Halfedge_descriptor h = get_vertex_to_edge_map(w1);
        set_halfedge(current_face, h, hds);
        ++new_faces;
        return h;
    }

    template <class InputIterator>
    Halfedge_descriptor add_facet( InputIterator first, InputIterator beyond) {
        // synonym for begin_facet(), a call to add_vertex_to_facet() for each iterator
        // value type, and end_facet().
        begin_facet();
        for ( ; ! m_error && first != beyond; ++first)
            add_vertex_to_facet( *first);
        if ( m_error)
            return Halfedge_descriptor();
        return end_facet();
    }

    template <class InputIterator>
    bool test_facet( InputIterator first, InputIterator beyond) {
        // tests if the facet described by the vertex indices in the 
        // range [first,beyond) can be inserted without creating a 
        // a non-manifold (and therefore invalid) situation.
        // First, create a copy of the indices and close it cyclically
        std::vector< std::size_t> indices( first, beyond);
        if ( indices.size() < 3)
            return false;
        indices.push_back( indices[0]);
        return test_facet_indices( indices);
    }

    bool test_facet_indices( std::vector< std::size_t> indices);

    void end_surface();
        // ends the construction.

    bool check_unconnected_vertices();
        // returns `true' if unconnected vertices are detected. If `verb'
        // is set to `true' debug information about the unconnected
        // vertices is printed.

    bool remove_unconnected_vertices();
  

    void rollback();

protected:
    Halfedge_descriptor lookup_hole( std::size_t w) {
        CGAL_assertion( w < new_vertices);
        return lookup_hole( get_vertex_to_edge_map( w));
    }

    size_type find_vertex( Vertex_descriptor v) {
        // Returns 0 if v == NULL.
        if ( v == Vertex_descriptor() )
            return 0;
        size_type n = 0;
        typename HDS::Vertex_iterator it, end;
        boost::tie(it,end) = hds.vertices();
        while ( *it != v) {
            CGAL_assertion( it != end);
            ++n;
            ++it;
        }
        n = n - rollback_v;
        return n;
    }

    size_type find_facet( Face_descriptor f) {
        // Returns 0 if f == NULL.
        if ( f == Face_descriptor())
            return 0;
        size_type n = 0;
        typename HDS::Face_iterator it,end;
        boost::tie(it,end) = hds.faces();
        while ( *it != f) {
            CGAL_assertion( it != end);
            ++n;
            ++it;
        }
        n = n - rollback_f;
        return n;
    }

    Halfedge_descriptor lookup_halfedge( size_type w, size_type v) {
        // Pre: 0 <= w,v < new_vertices
        // Case a: It exists an halfedge g from w to v:
        //     g must be a border halfedge and the facet of g->opposite()
        //     must be set and different from the current facet.
        //     Set the facet of g to the current facet. Return the
        //     halfedge pointing to g.
        // Case b: It exists no halfedge from w to v:
        //     Create a new pair of halfedges g and g->opposite().
        //     Set the facet of g to the current facet and g->opposite()
        //     to a border halfedge. Assign the vertex references.
        //     Set g->opposite()->next() to g. Return g->opposite().
       
        CGAL_assertion( w < new_vertices);
        CGAL_assertion( v < new_vertices);
        CGAL_assertion( ! last_vertex);

        Halfedge_descriptor e = get_vertex_to_edge_map( w);
        if ( e != Halfedge_descriptor()) {
          CGAL_assertion( target(e,hds) == index_to_vertex_map[w]);
            // check that the facet has no self intersections
            if ( current_face != Face_descriptor()
                 && current_face == face(e,hds)) {
                Verbose_ostream verr( m_verbose);
                verr << " " << std::endl;
                verr << "CGAL::Surface_mesh_incremental_builder<HDS>::"
                     << std::endl;
                verr << "lookup_halfedge(): input error: facet "
                     << new_faces << " has a self intersection at vertex "
                     << w << "." << std::endl;
                m_error = true;
                return Halfedge_descriptor();
            }
            Halfedge_descriptor start_edge(e);
            do {
              if ( target(next(e,hds),hds) == index_to_vertex_map[v]) {
                if ( ! is_border(next(e,hds),hds)) {
                        Verbose_ostream verr( m_verbose);
                        verr << " " << std::endl;
                        verr << "CGAL::Surface_mesh_incremental_builder"
                                "<HDS>::" << std::endl;
                        verr << "lookup_halfedge(): input error: facet "
                             << new_faces << " shares a halfedge from "
                                "vertex " <<  w << " to vertex " << v
                             << " with";
                        if (  m_verbose && current_face != Face_descriptor())
                            verr << " facet "
                                 << find_facet( face(next(e,hds), hds))
                                 << '.' << std::endl;
                        else
                            verr << " another facet." << std::endl;
                        m_error = true;
                        return Halfedge_descriptor();
                    }
                CGAL_assertion( ! is_border(opposite(next(e,hds),hds), hds) );
                    if ( current_face != Face_descriptor() && current_face ==
                         face(opposite(next(e,hds),hds),hds)) {
                        Verbose_ostream verr( m_verbose);
                        verr << " " << std::endl;
                        verr << "CGAL::Surface_mesh_incremental_builder"
                                "<HDS>::" << std::endl;
                        verr << "lookup_halfedge(): input error: facet "
                             << new_faces << " has a self intersection "
                                "at the halfedge from vertex " << w
                             << " to vertex " << v << "." << std::endl;
                        m_error = true;
                        return Halfedge_descriptor();
                    }
                    set_face( next(e,hds), current_face, hds);
                    // The following line prevents e->next() to be picked
                    // by get_vertex_to_edge_map(v) in an upcoming call
                    // of lookup_halfedge(v, *)
                    set_vertex_to_edge_map( v, opposite(next(next(e,hds),hds),hds));
                    return e;
                }
              e = opposite(next(e,hds),hds);
            } while ( e != start_edge);
        }
        // create a new halfedge
        
        e = halfedge(add_edge(hds),hds);
        new_halfedges++;
        new_halfedges++;
        set_face(e, current_face,hds);
        set_target(e, index_to_vertex_map[v], hds);
        set_next(e, Halfedge_descriptor(),hds);
        //??? decorator.set_prev( e, e->opposite());
        hds.set_prev_only(e, opposite(e,hds));
        e = opposite(e,hds);
        set_target(e, index_to_vertex_map[w], hds);
        set_next(e, opposite(e,hds),hds);
        return e;
    }

    Halfedge_descriptor lookup_hole( Halfedge_descriptor e) {
        // Halfedge e points to a vertex w. Walk around w to find a hole
        // in the facet structure. Report an error if none exist. Return
        // the halfedge at this hole that points to the vertex w.
        CGAL_assertion( e != Halfedge_descriptor());

        Halfedge_descriptor start_edge( e);
        do {
          if ( is_border(next(e,hds),hds)) {
                return e;
            }
          e = opposite(next(e,hds),hds);
        } while ( e != start_edge);

        Verbose_ostream verr( m_verbose);
        verr << " " << std::endl;
        verr << "CGAL::Surface_mesh_incremental_builder<HDS>::" << std::endl;
        verr << "lookup_hole(): input error: at vertex "
             << find_vertex( target(e,hds))
             << " a closed surface already exists and facet "
             << new_faces << " is nonetheless adjacent." << std::endl;
        if (  m_verbose && current_face != Face_descriptor()) {
            verr << "             The closed cycle of facets is:";
            do {
              if ( ! is_border(e,hds))
                verr << " " << find_facet( face(e,hds));
              e = opposite(next(e,hds),hds);
            } while ( e != start_edge);
            verr << '.' << std::endl;
        }
        m_error = true;
        return Halfedge_descriptor();
    }
};

template < class HDS>
void
Surface_mesh_incremental_builder<HDS>::
rollback() {
    CGAL_assertion( rollback_v <= num_vertices(hds));
    CGAL_assertion( rollback_h <= num_halfedges(hds));
    CGAL_assertion( rollback_f <= num_faces(hds));
    if ( rollback_v == 0 && rollback_h == 0 && rollback_f == 0) {
        hds.clear();
    } else {
      while ( rollback_v != (hds.num_vertices() - hds.num_removed_vertices()))
            hds.vertices_pop_back();

      while ( rollback_h != (num_halfedges() - hds.num_removed_halfedges()))
            hds.edges_pop_back();
      while ( rollback_f != (hds.num_faces() -  - hds.num_removed_faces()))
            hds.faces_pop_back();
    }
    m_error = false;
    CGAL_assertion_code( check_protocol = 0;)
      }

template < class HDS>  CGAL_MEDIUM_INLINE
void
Surface_mesh_incremental_builder<HDS>::
begin_surface( std::size_t v, std::size_t f, std::size_t h, int mode) {
    CGAL_assertion( check_protocol == 0);
    CGAL_assertion_code( check_protocol = 1;)
    CGAL_assertion( ! m_error);
    if ( mode == RELATIVE_INDEXING) {
        new_vertices  = 0;
        new_faces     = 0;
        new_halfedges = 0;
        rollback_v    = hds.num_vertices();
        rollback_f    = hds.num_faces();
        rollback_h    = hds.num_halfedges();
    } else {
        new_vertices  = hds.num_vertices();
        new_faces     = hds.num_faces();
        new_halfedges = hds.num_halfedges();
        rollback_v    = 0;
        rollback_f    = 0;
        rollback_h    = 0;
    }
    if ( h == 0) {
        // Use the Eulerian equation for connected planar graphs. We do
        // not know the number of facets that are holes and we do not
        // know the genus of the surface. So we add 12 and a factor of
        // 5 percent.
      h = (std::size_t)((double)(v + f - 2 + 12) * 2.1);
    }
    hds.reserve( hds.num_vertices()  + v,
                 hds.num_halfedges() + h,
                 hds.num_faces()     + f);
    if ( mode == RELATIVE_INDEXING) {
        index_to_vertex_map = Random_access_index( hds.vertices().second);
        index_to_vertex_map.reserve(v);
        initialize_vertex_to_edge_map( v, false);
    } else {
      Vertex_iterator b,e;
      boost::tie(b,e) = vertices(hds);
      index_to_vertex_map = Random_access_index(b,e);
        index_to_vertex_map.reserve( hds.num_vertices() + v);
        initialize_vertex_to_edge_map( hds.num_vertices() + v, true);
    }
}

template < class HDS>
void
Surface_mesh_incremental_builder<HDS>::
add_vertex_to_facet( std::size_t v2) {
    if ( m_error)
        return;
    CGAL_assertion( check_protocol == 2);
    if ( v2 >= new_vertices) {
        Verbose_ostream verr( m_verbose);
        verr << " " << std::endl;
        verr << "CGAL::Surface_mesh_incremental_builder<HDS>::"
             << std::endl;
        verr << "add_vertex_to_facet(): vertex index " << v2
             << " is out-of-range [0," << new_vertices-1 << "]."
             << std::endl;
        m_error = true;
        return;
    }

    if ( first_vertex) {
        CGAL_assertion( ! last_vertex);
        w1 = v2;
        first_vertex = false;
        return;
    }
    if ( g1 == Halfedge_descriptor()) {
        CGAL_assertion( ! last_vertex);
        gprime  = lookup_halfedge( w1, v2);
        if ( m_error)
            return;
        h1 = g1 = next(gprime,hds);
        v1 = w2 = v2;
        return;
    }
    // g1, h1, v1, w1, w2 are set. Insert halfedge.
    // Lookup v1-->v2
    Halfedge_descriptor hprime;
    if ( last_vertex)
        hprime = gprime;
    else {
        hprime = lookup_halfedge( v1, v2);
        if ( m_error)
            return;
    }
    Halfedge_descriptor h2 = next(hprime,hds);
    CGAL_assertion( ! last_vertex || h2 == g1);
    Halfedge_descriptor prev = next(h1,hds);
    set_next(h1, h2, hds);

    if ( get_vertex_to_edge_map( v1) == Halfedge_descriptor()) {  // case 1:
      set_next(opposite(h2,hds), opposite(h1,hds),hds);

    } else {                                                  // case 2:
      bool b1 = is_border(opposite(h1,hds), hds);
      bool b2 = is_border(opposite(h2,hds), hds);
        if ( b1 && b2) {
            Halfedge_descriptor hole = lookup_hole( v1);
            if ( m_error)
                return;
            CGAL_assertion( hole != Halfedge_descriptor());
            set_next(opposite(h2,hds),next(hole,hds), hds);
            set_next(hole, opposite(h1,hds), hds);
        } else if ( b2) {                                     // case 2.b:
          CGAL_assertion( is_border(prev,hds));
          set_next(opposite(h2,hds), prev, hds);
        } else if ( b1) {                                     // case 2.c:
          CGAL_assertion( is_border(hprime,hds));
          set_next(hprime, opposite(h1,hds),hds);

        } else if ( next(opposite(h2,hds),hds) == opposite(h1,hds)) {// case 2.d:
            // f1 == f2
          CGAL_assertion( face(opposite(h1,hds),hds) ==
                          face(opposite(h2,hds),hds) );
        } else {                                              // case 2.e:
            if ( prev == h2) {                                // case _i:
                // nothing to be done, hole is closed.
            } else {                                          // case _ii:
              CGAL_assertion( is_border(prev,hds));
              CGAL_assertion( is_border(hprime,hds));
              set_next(hprime, prev, hds);
                // Check whether the halfedges around v1 are connected.
                // It is sufficient to check it for h1 to prev.
                // Assert loop termination:
                CGAL_assertion_code( std::size_t k = 0;)
                // Look for a hole in the facet complex starting at h1.
                Halfedge_descriptor hole;
                Halfedge_descriptor e = h1;
                do {
                  if ( is_border(e,hds))
                        hole = e;
                  e = opposite(next(e,hds),hds);
                    CGAL_assertion( k++ < hds.num_halfedges());
                } while ( next(e,hds) != prev && e != h1);
                if ( e == h1) {
                    // disconnected facet complexes
                    if ( hole != Halfedge_descriptor()) {
                        // The complex can be connected with
                        // the hole at hprime.
                      set_next(hprime, next(hole,hds),hds);
                      set_next(hole, prev,hds);
                    } else {
                        Verbose_ostream verr( m_verbose);
                        verr << " " << std::endl;
                        verr << "CGAL::Surface_mesh_incremental_builder<"
                                "HDS>::" << std::endl;
                        verr << "add_vertex_to_facet(): input error: "
                                "disconnected facet complexes at vertex "
                             << v1 << ":" << std::endl;

                        if ( m_verbose && current_face != Face_descriptor()) {
                            verr << "           involved facets are:";
                            do {
                              if ( ! is_border(e,hds))
                                verr << " " << find_facet(face(e,hds));
                              e = opposite(next(e,hds),hds);
                            } while ( e != h1);
                            verr << " (closed cycle) and";
                            e = hprime;
                            do {
                              if ( ! is_border(e,hds))
                                verr << " " << find_facet(face(e,hds));
                            } while ( e != hprime);
                            verr << "." << std::endl;
                        }
                        m_error = true;
                        return;
                    }
                }
            }
        }
    }
    if ( target(h1,hds) == index_to_vertex_map[v1])
        set_vertex_to_edge_map( v1, h1);
    CGAL_assertion( target(h1,hds) == index_to_vertex_map[v1]);
    h1 = h2;
    v1 = v2;
}

template < class HDS>
bool
Surface_mesh_incremental_builder<HDS>::
test_facet_indices( std::vector< std::size_t> indices) {

    // tests if the facet described by the vertex indices can be inserted 
    // without creating a non-manifold (and therefore invalid) situation.
    // indices are cyclically closed once.
    std::size_t n = indices.size() - 1;
    // Test if a vertex is not twice in the indices
    for ( std::size_t i = 0; i < n; ++i) {
        CGAL_precondition( indices[i] < new_vertices);
        // check if vertex indices[i] is already in the sequence [0..i-1]
        for ( std::size_t k = 0; k < i; ++k) {
            if ( indices[k] == indices[i])
                return false;
        }
    }
    // Test non-manifold halfedges
    for ( std::size_t i = 0; i < n; ++i) {
        // halfedge goes from vertex indices[i] to indices[i+1]
        // we know already that the halfedge is only once in the sequence
        // (otherwise the end-vertices would be twice in the sequence too)
        // check if halfedge is already in the HDS and is not border halfedge
        Halfedge_descriptor v = get_vertex_to_edge_map(indices[i]);
        Vertex_descriptor   w = index_to_vertex_map[indices[i+1]];
        if ( v != Halfedge_descriptor()
             && get_vertex_to_edge_map(indices[i+1]) != Halfedge_descriptor()) {
            // cycle through halfedge-loop and find edge to indices[i+1]
            Halfedge_descriptor vstart = v;
            do {
              v = opposite(next(,hds),hds);
            } while ( target(next(v,hds),hds) != w && v != vstart);
            if ( target(next(v,hds),hds) == w && ! is_border(next(v,hds),hds))
                return false;
        }
    }
    // test non-manifold vertices
    for ( std::size_t i = 0; i < n; ++i) {
        // since we don't allow duplicates in indices[..] and we 
        // tested for non-manifold halfedges already, we just need to check
        // if the vertex indices[i] is not a closed manifold yet.
        Halfedge_descriptor v = get_vertex_to_edge_map(indices[i]);
        if ( v != Halfedge_descriptor()) {
            Halfedge_descriptor vstart = v;
            do {
                v = opposite(next(v,hds),hds);
            } while ( ! is_border(v,hds) && v != vstart);
            if ( ! is_border(v,hds))
                return false;
        }
    }
    
    //Test if all halfedges of the new face 
    //are possibly consecutive border halfedges in the HDS.
    //Possibly because it may be not directly encoded in the HDS
    //(using next() function ). This situation can occur when one or
    //more facets share only a vertex: For example, the new facet we try to add
    //would make the vertex indices[i] a manifold but this should be forbidden
    //if a facet only incident to that vertex has already been inserted.
    //We check this for each vertex of the sequence.
    for ( std::size_t i = 0; i < n; ++i) {
      std::size_t prev_index=indices[ (i-1+n)%n];
      std::size_t next_index=indices[ (i+1)%n];
      Vertex_descriptor   previous_vertex = index_to_vertex_map[ prev_index ];
      Vertex_descriptor   next_vertex     = index_to_vertex_map[ next_index ];
      
      Halfedge_descriptor v = get_vertex_to_edge_map(indices[i]);
      
      if ( v == Halfedge_descriptor() || 
           get_vertex_to_edge_map(prev_index) == Halfedge_descriptor() ||
           get_vertex_to_edge_map(next_index) == Halfedge_descriptor()
         ) continue;
      
      Halfedge_descriptor start=v;
      //halfedges pointing to/running out from vertex indices[i]
      //and that need to be possibly consecutive
      Halfedge_descriptor previous=Halfedge_descriptor(),next=Halfedge_descriptor();
      
      //look for a halfedge incident to vertex indices[i]
      //and which opposite is incident to previous_vertex
      do{
        if (target(opposite(v,hds),hds)==previous_vertex){
          previous=v;
          CGAL_precondition(is_border(previous,hds));
          break;
        }
        v = opposite(next(v,hds),hds);
      }
      while (v!=start);
      
      if (previous!=Halfedge_descriptor()){
        v = opposite(next(v,hds),hds);
        //previous and next are already consecutive in the HDS
        if (target(opposite(v,hds),hds)==next_vertex) continue;
        
        //look for a border halfedge which opposite is
        //incident to next_vertex: set next halfedge
        do
        {
          if (target(opposite(v,hds),hds)==next_vertex){
            next = opposite(v,hd);
            break;
          }
          v= opposite(next(v,hds));
        }
        while(v!=previous);
        if (next==Halfedge_descriptor()) continue;
        
        //check if no constraint prevents
        //previous and next to be adjacent: 
        do{
          v = opposite(next(v,hds),hds);
          if ( is_border(opposite(v,hds),hds) ) break;
        }
        while (v!=previous);
        if (v==previous) return false;
        start=v;
      }
    }
    
    return true;
}


template < class HDS>  CGAL_MEDIUM_INLINE
void
Surface_mesh_incremental_builder<HDS>::
end_surface() {
    if ( m_error)
        return;
    CGAL_assertion( check_protocol == 1);
    CGAL_assertion_code( check_protocol = 0;)
}

template < class HDS>
bool
Surface_mesh_incremental_builder<HDS>::
check_unconnected_vertices() {
    if ( m_error)
        return false;
    bool unconnected = false;
    Verbose_ostream verr( m_verbose);
    for ( std::size_t i = 0; i < new_vertices; i++) {
        if ( get_vertex_to_edge_map( i) == Halfedge_descriptor()) {
            verr << "CGAL::Surface_mesh_incremental_builder<HDS>::\n"
                 << "check_unconnected_vertices( verb = true): "
                 << "vertex " << i << " is unconnected." << std::endl;
            unconnected = true;
        }
    }
    return unconnected;
}

template < class HDS>
bool
Surface_mesh_incremental_builder<HDS>::
remove_unconnected_vertices() {
    if ( m_error)
        return true;
    for( std::size_t i = 0; i < new_vertices; i++) {
        if( get_vertex_to_edge_map( i) == Halfedge_descriptor()) {
          remove_vertex( index_to_vertex_map[i], hds);
        }
    }
    return true;
}

} //namespace CGAL

#endif // CGAL_SURFACE_MESH_INCREMENTAL_BUILDER_H //
// EOF //
