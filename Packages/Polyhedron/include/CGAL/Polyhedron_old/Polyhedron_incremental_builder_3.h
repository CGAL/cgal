// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : Polyhedron_incremental_builder_3.h
// chapter       : $CGAL_Chapter: 3D-Polyhedral Surfaces $
// package       : $CGAL_Package: Polyhedron 2.9 (13 Sep 2000) $
// source        : polyhedron_builder.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : MPI Saarbruecken (Stefan Schirra <stschirr@mpi-sb.mpg.de>)
//
// Incremental Construction of Polyhedral Surfaces.
// ============================================================================

#ifndef CGAL_POLYHEDRON_OLD_POLYHEDRON_INCREMENTAL_BUILDER_3_H
#define CGAL_POLYHEDRON_OLD_POLYHEDRON_INCREMENTAL_BUILDER_3_H 1

#ifndef CGAL_RANDOM_ACCESS_VALUE_ADAPTOR_H
#include <CGAL/Random_access_value_adaptor.h>
#endif // CGAL_RANDOM_ACCESS_VALUE_ADAPTOR_H
#ifndef CGAL_PROTECT_VECTOR
#include <vector>
#define CGAL_PROTECT_VECTOR
#endif // CGAL_PROTECT_VECTOR
#ifndef CGAL_HALFEDGE_DATA_STRUCTURE_DECORATOR_H
#include <CGAL/Halfedge_data_structure_decorator.h>
#endif // CGAL_HALFEDGE_DATA_STRUCTURE_DECORATOR_H
#ifndef CGAL_IO_VERBOSE_OSTREAM_H
#include <CGAL/IO/Verbose_ostream.h>
#endif // CGAL_IO_VERBOSE_OSTREAM_H

CGAL_BEGIN_NAMESPACE

template < class HDS>
class Polyhedron_incremental_builder_3 {
public:
    typedef HDS                      Halfedge_data_structure;
    typedef HDS                      HalfedgeDS;
    typedef typename HDS::Vertex     Vertex;
    typedef typename HDS::Halfedge   Halfedge;
    typedef typename HDS::Facet      Facet;
    typedef typename HDS::Point      Point;
    typedef typename HDS::Point      Point_3;
    typedef typename HDS::Size       Size;

protected:
    typedef typename HDS::Supports_vertex_halfedge
        Supports_vertex_halfedge;
    typedef typename HDS::Supports_removal          Supports_removal;
    typedef typename HDS::Vertex_iterator           Vertex_iterator;
    typedef Random_access_value_adaptor<Vertex_iterator,Vertex>
                                                    Random_access_index;

    bool                    m_error;
    bool                    m_verbose;
    HDS&                    hds;
    Size                    rollback_v;
    Size                    rollback_f;
    Size                    rollback_h;
    Size                    new_vertices;
    Size                    new_facets;
    Size                    new_halfedges;
    Facet*                  current_facet;
    Random_access_index     index_to_vertex_map;
    std::vector<Halfedge*>  vertex_to_edge_map;

    Halfedge*               g1;      // first halfedge, 0 denotes none.
    Halfedge*               gprime;
    Halfedge*               h1;      // current halfedge
    Size                    w1;      // first vertex.
    Size                    w2;      // second vertex.
    Size                    v1;      // current vertex
    bool                    first_vertex;
    bool                    last_vertex;

    CGAL_assertion_code( int check_protocoll;)  // use to check protocoll.
    // states for checking: 0 = created, 1 = constructing, 2 = make facet.

    // Implement the vertex_to_edge_map either with an array or
    // the halfedge pointer in the vertices (if supported).
    // ----------------------------------------------------
    void initialize_vertex_to_edge_map( Size  , Tag_true) {}
    void initialize_vertex_to_edge_map( Size n, Tag_false) {
        vertex_to_edge_map = std::vector<Halfedge*>();
        vertex_to_edge_map.reserve(n);
    }
    void initialize_vertex_to_edge_map( Size n) {
        initialize_vertex_to_edge_map( n, Supports_vertex_halfedge());
    }
    void push_back_vertex_to_edge_map( Halfedge*  , Tag_true) {}
    void push_back_vertex_to_edge_map( Halfedge* h, Tag_false) {
        vertex_to_edge_map.push_back(h);
    }
    void push_back_vertex_to_edge_map( Halfedge* h) {
        push_back_vertex_to_edge_map( h, Supports_vertex_halfedge());
    }
    Halfedge* get_vertex_to_edge_map( int i, Tag_true) {
        // Use the halfedge pointer within the vertex.
        return index_to_vertex_map[i].halfedge();
    }
    Halfedge* get_vertex_to_edge_map( int i, Tag_false) {
        // Use the self-managed array vertex_to_edge_map.
        return vertex_to_edge_map[i];
    }
    Halfedge* get_vertex_to_edge_map( int i) {
        return get_vertex_to_edge_map( i, Supports_vertex_halfedge());
    }
    void set_vertex_to_edge_map( int i, Halfedge* h, Tag_true) {
        // Use the halfedge pointer within the vertex.
        index_to_vertex_map[i].set_halfedge(h);
    }
    void set_vertex_to_edge_map( int i, Halfedge* h, Tag_false) {
        // Use the self-managed array vertex_to_edge_map.
        vertex_to_edge_map[i] = h;
    }
    void set_vertex_to_edge_map( int i, Halfedge* h) {
        set_vertex_to_edge_map( i, h, Supports_vertex_halfedge());
    }

// An Incremental Builder for Polyhedral Surfaces
// ----------------------------------------------
// DEFINITION
//
// Polyhedron_incremental_builder_3<HDS> is an auxiliary class that
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

    Polyhedron_incremental_builder_3(HDS& h, bool verbose = false)
        // stores a reference to the halfedge data structure `h' in the
        // internal state. The previous polyhedral surface in `h'
        // remains unchanged. The incremental builder adds the new
        // polyhedral surface to the old one.
      : m_error( false), m_verbose( verbose), hds(h) {
        CGAL_assertion_code(check_protocoll = 0;)
    }

    ~Polyhedron_incremental_builder_3() {
        CGAL_assertion( check_protocoll == 0);
    }

// OPERATIONS

    #ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
    void begin_surface( Size v, Size f, Size h = 0);
    #else
    void begin_surface( std::size_t v, std::size_t f, std::size_t h = 0);
    #endif
        // starts the construction. v is the number of
        // vertices to expect, f the number of facets, and h the number of
        // halfedges. If h is unspecified (`== 0') it is estimated using
        // Euler equations (plus 5% for the so far unkown holes and genus
        // of the object). These values are used to reserve space in the
        // polyhedron representation `HDS'. If the representation
        // supports insertion these
        // values do not restrict the class of readable polyhedrons.
        // If the representation does not support insertion the object
        // must fit in the reserved sizes.

#ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
    void add_vertex( const Point& p);
        // adds p to the vertex list.
#endif

    void begin_facet();
        // starts a facet.

#ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
    void add_vertex_to_facet( Size i);
#else
    void add_vertex_to_facet( std::size_t i);
#endif
        // adds a vertex with index i to the current facet. The first
        // point added with `add_vertex()' has the index 0.

    void end_facet();
        // ends a facet.

    void end_surface();
        // ends the construction.

    bool check_unconnected_vertices();
        // returns `true' if unconnected vertices are detected. If `verb'
        // is set to `true' debug information about the unconnected
        // vertices is printed.

    bool remove_unconnected_vertices( Tag_true);
    bool remove_unconnected_vertices( Tag_false) {
        return ! check_unconnected_vertices();
    }
    bool remove_unconnected_vertices() {
        // returns `true' if all unconnected vertices could be removed
        // succesfully.
        return remove_unconnected_vertices( Supports_removal());
    }

    void rollback();

protected:
    Halfedge* lookup_hole( Size w) {
        CGAL_assertion( w < new_vertices);
        return lookup_hole( get_vertex_to_edge_map( w));
    }

#ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
    Size find_vertex( Vertex* v);
        // Returns 0 if v == NULL.

    Size find_facet( Facet* f);
        // Returns 0 if f == NULL.

    Halfedge* lookup_halfedge( Size w, Size v);
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

    Halfedge* lookup_hole( Halfedge* e);
        // Halfedge e points to a vertex w. Walk around w to find a hole
        // in the facet structure. Report an error if none exist. Return
        // the halfedge at this hole that points to the vertex w.

#else // CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS //
    #ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
    template < class HDS>  CGAL_LARGE_INLINE
    typename HDS::Size
    Polyhedron_incremental_builder_3<HDS>::
    #else
    Size
    #endif
    find_vertex( Vertex* v) {
        if ( ! v)
            return 0;
        Size n = 0;
        typename HDS::Vertex_iterator it = hds.vertices_begin();
        while ( &(*it) != v) {
            CGAL_assertion( it != hds.vertices_end());
            ++n;
            ++it;
        }
        n = n - ( hds.size_of_vertices() - new_vertices);
        return n;
    }
    
    #ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
    template < class HDS>  CGAL_LARGE_INLINE
    typename HDS::Size
    Polyhedron_incremental_builder_3<HDS>::
    #else
    Size
    #endif
    find_facet( Facet* f) {
        if ( ! f)
            return 0;
        Size n = 0;
        typename HDS::Facet_iterator it = hds.facets_begin();
        while ( &(*it) != f) {
            CGAL_assertion( it != hds.facets_end());
            ++n;
            ++it;
        }
        n = n - ( hds.size_of_facets() - new_facets);
        return n;
    }
    
    #ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
    template < class HDS>  CGAL_LARGE_INLINE
    typename HDS::Halfedge*
    Polyhedron_incremental_builder_3<HDS>::
    #else
    Halfedge*
    #endif
    lookup_halfedge( Size w, Size v) {
        typedef typename HDS::Supports_halfedge_vertex
            Supports_halfedge_vertex;
        Assert_compile_time_tag( Supports_halfedge_vertex(), Tag_true());
        CGAL_assertion( w < new_vertices);
        CGAL_assertion( v < new_vertices);
        CGAL_assertion( ! last_vertex);
        Verbose_ostream verr( m_verbose);
        Halfedge_data_structure_decorator<HDS> decorator;
        Halfedge* e = get_vertex_to_edge_map( w);
        if ( e) {
            CGAL_assertion( e->vertex() == &(index_to_vertex_map[w]));
            // check that the facet has no self intersections
            if ( current_facet && current_facet == decorator.get_facet(e)) {
                verr << " " << std::endl;
                verr << "Polyhedron_incremental_builder_3<HDS>::"
                     << std::endl;
                verr << "lookup_halfedge(): input error: facet "
                     << new_facets << " has a self intersection at vertex "
                     << w << "." << std::endl;
                m_error = true;
                return 0;
            }
            Halfedge* start_edge = e;
            do {
                if ( e->next()->vertex() == &(index_to_vertex_map[v]) ) {
                    if ( ! e->next()->is_border()) {
                        verr << " " << std::endl;
                        verr << "Polyhedron_incremental_builder_3"
                                "<HDS>::" << std::endl;
                        verr << "lookup_halfedge(): input error: facet "
                             << new_facets << " shares a halfedge from "
                                "vertex " <<  w << " to vertex " << v
                             << " with";
                        if ( current_facet)
                            verr << " facet "
                                 << find_facet( decorator.get_facet(e->next()))
                                 << '.' << std::endl;
                        else
                            verr << " another facet." << std::endl;
                        m_error = true;
                        return 0;
                    }
                    CGAL_assertion( ! e->next()->opposite()->is_border());
                    if ( current_facet &&
                         current_facet ==
                         decorator.get_facet( e->next()->opposite())) {
                        verr << " " << std::endl;
                        verr << "Polyhedron_incremental_builder_3"
                                "<HDS>::" << std::endl;
                        verr << "lookup_halfedge(): input error: facet "
                             << new_facets << " has a self intersection "
                                "at the halfedge from vertex " << w
                             << " to vertex " << v << "." << std::endl;
                        m_error = true;
                        return 0;
                    }
                    decorator.set_facet( e->next(), current_facet);
                    return e;
                }
                e = e->next()->opposite();
            } while ( e != start_edge);
        }
        // create a new halfedge
        if ( hds.size_of_halfedges() >= hds.capacity_of_halfedges()) {
            verr << " " << std::endl;
            verr << "Polyhedron_incremental_builder_3<HDS>::" << std::endl;
            verr << "lookup_halfedge(): capacity error: more than "
                 << new_halfedges << " halfedges added while creating facet"
                 << new_facets << '.' << std::endl;
            m_error = true;
            return 0;
        }
        e = hds.new_edge();
        new_halfedges++;
        new_halfedges++;
        decorator.set_facet( e, current_facet);
        e->set_vertex( &(index_to_vertex_map[v]));
        e->set_next( 0);
        decorator.set_prev( e, e->opposite());
        e = e->opposite();
        e->set_vertex( &(index_to_vertex_map[w]));
        e->set_next( e->opposite());
        return e;
    }
    
    #ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
    template < class HDS>  CGAL_LARGE_INLINE
    typename HDS::Halfedge*
    Polyhedron_incremental_builder_3<HDS>::
    #else
    Halfedge*
    #endif
    lookup_hole( Halfedge* e) {
        CGAL_assertion( e != NULL);
        Halfedge_data_structure_decorator<HDS> decorator;
        Halfedge* start_edge = e;
        do {
            if ( e->next()->is_border()) {
                return e;
            }
            e = e->next()->opposite();
        } while ( e != start_edge);
    
        Verbose_ostream verr( m_verbose);
        verr << " " << std::endl;
        verr << "Polyhedron_incremental_builder_3<HDS>::" << std::endl;
        verr << "lookup_hole(): input error: at vertex "
            << find_vertex( e->vertex())
            << " a closed surface already exists and facet "
            << new_facets << " is nonetheless adjacent." << std::endl;
        if ( current_facet) {
            verr << "             The closed cycle of facets is:";
            do {
                if ( ! e->is_border())
                    verr << " " << find_facet( decorator.get_facet(e));
                e = e->next()->opposite();
            } while ( e != start_edge);
            verr << '.' << std::endl;
        }
        m_error = true;
        return 0;
    }
    
    #ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
    template < class HDS>  CGAL_MEDIUM_INLINE
    void
    Polyhedron_incremental_builder_3<HDS>::
    #else
    public:
    void
    #endif
    add_vertex( const Point& p) {
        CGAL_assertion( check_protocoll == 1);
        if ( hds.size_of_vertices() >= hds.capacity_of_vertices()) {
            Verbose_ostream verr( m_verbose);
            verr << " " << std::endl;
            verr << "Polyhedron_incremental_builder_3<HDS>::" << std::endl;
            verr << "add_vertex(): capacity error: more than " << new_vertices
                 << " vertices added." << std::endl;
            m_error = true;
            return;
        }
        Halfedge_data_structure_decorator<HDS> decorator;
        Vertex* v = decorator.new_vertex( hds, p);
        index_to_vertex_map.push_back( Vertex_iterator(v));
        decorator.set_vertex_halfedge( v, 0);
        push_back_vertex_to_edge_map( 0);
        ++new_vertices;
    }    
#endif // CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS //
};

#ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
#ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
template < class HDS>  CGAL_LARGE_INLINE
typename HDS::Size
Polyhedron_incremental_builder_3<HDS>::
#else
Size
#endif
find_vertex( Vertex* v) {
    if ( ! v)
        return 0;
    Size n = 0;
    typename HDS::Vertex_iterator it = hds.vertices_begin();
    while ( &(*it) != v) {
        CGAL_assertion( it != hds.vertices_end());
        ++n;
        ++it;
    }
    n = n - ( hds.size_of_vertices() - new_vertices);
    return n;
}

#ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
template < class HDS>  CGAL_LARGE_INLINE
typename HDS::Size
Polyhedron_incremental_builder_3<HDS>::
#else
Size
#endif
find_facet( Facet* f) {
    if ( ! f)
        return 0;
    Size n = 0;
    typename HDS::Facet_iterator it = hds.facets_begin();
    while ( &(*it) != f) {
        CGAL_assertion( it != hds.facets_end());
        ++n;
        ++it;
    }
    n = n - ( hds.size_of_facets() - new_facets);
    return n;
}

#ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
template < class HDS>  CGAL_LARGE_INLINE
typename HDS::Halfedge*
Polyhedron_incremental_builder_3<HDS>::
#else
Halfedge*
#endif
lookup_halfedge( Size w, Size v) {
    typedef typename HDS::Supports_halfedge_vertex
        Supports_halfedge_vertex;
    Assert_compile_time_tag( Supports_halfedge_vertex(), Tag_true());
    CGAL_assertion( w < new_vertices);
    CGAL_assertion( v < new_vertices);
    CGAL_assertion( ! last_vertex);
    Verbose_ostream verr( m_verbose);
    Halfedge_data_structure_decorator<HDS> decorator;
    Halfedge* e = get_vertex_to_edge_map( w);
    if ( e) {
        CGAL_assertion( e->vertex() == &(index_to_vertex_map[w]));
        // check that the facet has no self intersections
        if ( current_facet && current_facet == decorator.get_facet(e)) {
            verr << " " << std::endl;
            verr << "Polyhedron_incremental_builder_3<HDS>::"
                 << std::endl;
            verr << "lookup_halfedge(): input error: facet "
                 << new_facets << " has a self intersection at vertex "
                 << w << "." << std::endl;
            m_error = true;
            return 0;
        }
        Halfedge* start_edge = e;
        do {
            if ( e->next()->vertex() == &(index_to_vertex_map[v]) ) {
                if ( ! e->next()->is_border()) {
                    verr << " " << std::endl;
                    verr << "Polyhedron_incremental_builder_3"
                            "<HDS>::" << std::endl;
                    verr << "lookup_halfedge(): input error: facet "
                         << new_facets << " shares a halfedge from "
                            "vertex " <<  w << " to vertex " << v
                         << " with";
                    if ( current_facet)
                        verr << " facet "
                             << find_facet( decorator.get_facet(e->next()))
                             << '.' << std::endl;
                    else
                        verr << " another facet." << std::endl;
                    m_error = true;
                    return 0;
                }
                CGAL_assertion( ! e->next()->opposite()->is_border());
                if ( current_facet &&
                     current_facet ==
                     decorator.get_facet( e->next()->opposite())) {
                    verr << " " << std::endl;
                    verr << "Polyhedron_incremental_builder_3"
                            "<HDS>::" << std::endl;
                    verr << "lookup_halfedge(): input error: facet "
                         << new_facets << " has a self intersection "
                            "at the halfedge from vertex " << w
                         << " to vertex " << v << "." << std::endl;
                    m_error = true;
                    return 0;
                }
                decorator.set_facet( e->next(), current_facet);
                return e;
            }
            e = e->next()->opposite();
        } while ( e != start_edge);
    }
    // create a new halfedge
    if ( hds.size_of_halfedges() >= hds.capacity_of_halfedges()) {
        verr << " " << std::endl;
        verr << "Polyhedron_incremental_builder_3<HDS>::" << std::endl;
        verr << "lookup_halfedge(): capacity error: more than "
             << new_halfedges << " halfedges added while creating facet"
             << new_facets << '.' << std::endl;
        m_error = true;
        return 0;
    }
    e = hds.new_edge();
    new_halfedges++;
    new_halfedges++;
    decorator.set_facet( e, current_facet);
    e->set_vertex( &(index_to_vertex_map[v]));
    e->set_next( 0);
    decorator.set_prev( e, e->opposite());
    e = e->opposite();
    e->set_vertex( &(index_to_vertex_map[w]));
    e->set_next( e->opposite());
    return e;
}

#ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
template < class HDS>  CGAL_LARGE_INLINE
typename HDS::Halfedge*
Polyhedron_incremental_builder_3<HDS>::
#else
Halfedge*
#endif
lookup_hole( Halfedge* e) {
    CGAL_assertion( e != NULL);
    Halfedge_data_structure_decorator<HDS> decorator;
    Halfedge* start_edge = e;
    do {
        if ( e->next()->is_border()) {
            return e;
        }
        e = e->next()->opposite();
    } while ( e != start_edge);

    Verbose_ostream verr( m_verbose);
    verr << " " << std::endl;
    verr << "Polyhedron_incremental_builder_3<HDS>::" << std::endl;
    verr << "lookup_hole(): input error: at vertex "
        << find_vertex( e->vertex())
        << " a closed surface already exists and facet "
        << new_facets << " is nonetheless adjacent." << std::endl;
    if ( current_facet) {
        verr << "             The closed cycle of facets is:";
        do {
            if ( ! e->is_border())
                verr << " " << find_facet( decorator.get_facet(e));
            e = e->next()->opposite();
        } while ( e != start_edge);
        verr << '.' << std::endl;
    }
    m_error = true;
    return 0;
}

#ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
template < class HDS>  CGAL_MEDIUM_INLINE
void
Polyhedron_incremental_builder_3<HDS>::
#else
public:
void
#endif
add_vertex( const Point& p) {
    CGAL_assertion( check_protocoll == 1);
    if ( hds.size_of_vertices() >= hds.capacity_of_vertices()) {
        Verbose_ostream verr( m_verbose);
        verr << " " << std::endl;
        verr << "Polyhedron_incremental_builder_3<HDS>::" << std::endl;
        verr << "add_vertex(): capacity error: more than " << new_vertices
             << " vertices added." << std::endl;
        m_error = true;
        return;
    }
    Halfedge_data_structure_decorator<HDS> decorator;
    Vertex* v = decorator.new_vertex( hds, p);
    index_to_vertex_map.push_back( Vertex_iterator(v));
    decorator.set_vertex_halfedge( v, 0);
    push_back_vertex_to_edge_map( 0);
    ++new_vertices;
}    
#endif // CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS //
template < class HDS>
void
Polyhedron_incremental_builder_3<HDS>::
rollback() {
    CGAL_assertion( rollback_v <= hds.size_of_vertices());
    CGAL_assertion( rollback_h <= hds.size_of_halfedges());
    CGAL_assertion( rollback_f <= hds.size_of_facets());
    while ( rollback_v != hds.size_of_vertices())
        hds.vertex_pop_back();
    CGAL_assertion((( hds.size_of_halfedges() - rollback_h) & 1) == 0);
    while ( rollback_h != hds.size_of_halfedges())
        hds.edge_pop_back();
    while ( rollback_f != hds.size_of_facets())
        hds.facet_pop_back();
    m_error = false;
    CGAL_assertion_code( check_protocoll = 0;)
}

template < class HDS>  CGAL_MEDIUM_INLINE
void
Polyhedron_incremental_builder_3<HDS>::
#ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
begin_surface( Size v, Size f, Size h) {
#else
begin_surface( std::size_t v, std::size_t f, std::size_t h) {
#endif
    CGAL_assertion( check_protocoll == 0);
    CGAL_assertion_code( check_protocoll = 1;)
    CGAL_assertion( ! m_error);
    new_vertices  = 0;
    new_facets    = 0;
    new_halfedges = 0;
    rollback_v = hds.size_of_vertices();
    rollback_f = hds.size_of_facets();
    rollback_h = hds.size_of_halfedges();
    if ( h == 0) {
        // Use the Eulerian equation for connected planar graphs. We do
        // not know the number of facets that are holes and we do not
        // know the genus of the surface. So we add 12 and a factor of
        // 5 percent.
        h = int((v + f - 2 + 12) * 2.1);
    }
    hds.reserve( hds.size_of_vertices()  + v,
                 hds.size_of_halfedges() + h,
                 hds.size_of_facets()    + f);
    index_to_vertex_map = Random_access_index( hds.vertices_end());
    index_to_vertex_map.reserve(v);
    initialize_vertex_to_edge_map( v);
}

template < class HDS>  CGAL_MEDIUM_INLINE
void
Polyhedron_incremental_builder_3<HDS>::
begin_facet() {
    if ( m_error)
        return;
    CGAL_assertion( check_protocoll == 1);
    CGAL_assertion_code( check_protocoll = 2;)
    if ( hds.size_of_facets() >= hds.capacity_of_facets()) {
        Verbose_ostream verr( m_verbose);
        verr << " " << std::endl;
        verr << "Polyhedron_incremental_builder_3<HDS>::" << std::endl;
        verr << "begin_facet(): capacity error: more than " << new_vertices
             << " facets added." << std::endl;
        m_error = true;
        return;
    }
    // initialize all status variables.
    first_vertex = true;  // denotes 'no vertex yet'
    g1 =  0;  // denotes 'no halfedge yet'
    last_vertex = false;

    Halfedge_data_structure_decorator<HDS> decorator;
    current_facet = decorator.new_facet( hds);
}

template < class HDS>
void
Polyhedron_incremental_builder_3<HDS>::
#ifndef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
add_vertex_to_facet( Size v2) {
#else
add_vertex_to_facet( std::size_t v2) {
#endif
    if ( m_error)
        return;
    Verbose_ostream verr( m_verbose);
    CGAL_assertion( check_protocoll == 2);
    if ( v2 >= new_vertices) {
        verr << " " << std::endl;
        verr << "Polyhedron_incremental_builder_3<HDS>::" << std::endl;
        verr << "add_vertex_to_facet(): vertex index " << v2
             << " is out-of-range [0," << new_vertices-1 << "]."
             << std::endl;
        m_error = true;
        return;
    }
    Halfedge_data_structure_decorator<HDS> decorator;

    if ( first_vertex) {
        CGAL_assertion( ! last_vertex);
        w1 = v2;
        first_vertex = false;
        return;
    }
    if ( g1 == 0) {
        CGAL_assertion( ! last_vertex);
        gprime  = lookup_halfedge( w1, v2);
        if ( m_error)
            return;
        h1 = g1 = gprime->next();
        v1 = w2 = v2;
        return;
    }
    // g1, h1, v1, w1, w2 are set. Insert halfedge.
    // Lookup v1-->v2
    Halfedge* hprime;
    if ( last_vertex)
        hprime = gprime;
    else {
        hprime = lookup_halfedge( v1, v2);
        if ( m_error)
            return;
    }
    Halfedge* h2 = hprime->next();
    CGAL_assertion( ! last_vertex || h2 == g1);
    Halfedge* prev = h1->next();
    h1->set_next( h2);
    decorator.set_prev( h2, h1);

    if ( get_vertex_to_edge_map( v1) == 0) {                  // case 1:
        h2->opposite()->set_next( h1->opposite());
        decorator.set_prev( h1->opposite(), h2->opposite());
    } else {                                                  // case 2:
        bool b1 = h1->opposite()->is_border();
        bool b2 = h2->opposite()->is_border();
        if ( b1 && b2) {
            Halfedge* hole = lookup_hole( v1);
            if ( m_error)
                return;
            CGAL_assertion( hole != NULL);
            h2->opposite()->set_next( hole->next());
            decorator.set_prev( hole->next(), h2->opposite());
            hole->set_next( h1->opposite());
            decorator.set_prev( h1->opposite(), hole);
        } else if ( b2) {                                     // case 2.b:
            CGAL_assertion( prev->is_border());
            h2->opposite()->set_next( prev);
            decorator.set_prev( prev, h2->opposite());
        } else if ( b1) {                                     // case 2.c:
            CGAL_assertion( hprime->is_border());
            hprime->set_next( h1->opposite());
            decorator.set_prev( h1->opposite(), hprime);
        } else if ( h2->opposite()->next() == h1->opposite()) {// case 2.d:
            // f1 == f2
            CGAL_assertion( decorator.get_facet( h1->opposite()) ==
                       decorator.get_facet( h2->opposite()));
        } else {                                              // case 2.e:
            if ( prev == h2) {                                // case _i:
                // nothing to be done, hole is closed.
            } else {                                          // case _ii:
                CGAL_assertion( prev->is_border());
                CGAL_assertion( hprime->is_border());
                hprime->set_next( prev);
                decorator.set_prev( prev, hprime);
                // Check whether the halfedges around v1 are connected.
                // It is sufficient to check it for h1 to prev.
                // Assert loop termination.
                CGAL_assertion_code( std::size_t k = 0;)
                // Look for a hole in the facet complex starting at h1.
                Halfedge* hole = 0;
                Halfedge* e    = h1;
                do {
                    if ( e->is_border())
                        hole = e;
                    e = e->next()->opposite();
                    CGAL_assertion( k++ < hds.size_of_halfedges());
                } while ( e->next() != prev && e != h1);
                if ( e == h1) {
                    // disconnected facet complexes
                    if ( hole) {
                        // The complex can be connected with
                        // the hole at hprime.
                        hprime->set_next( hole->next());
                        decorator.set_prev( hole->next(), hprime);
                        hole->set_next( prev);
                        decorator.set_prev( prev, hole);
                    } else {
                        verr << " " << std::endl;
                        verr << "Polyhedron_incremental_builder_3<"
                                "HDS>::" << std::endl;
                        verr << "add_vertex_to_facet(): input error: "
                                "disconnected facet complexes at vertex "
                             << v1 << ":" << std::endl;

                        if ( current_facet && m_verbose) {
                            verr << "           involved facets are:";
                            do {
                                if ( ! e->is_border())
                                    verr << " " << find_facet(
                                                decorator.get_facet(e));
                                e = e->next()->opposite();
                            } while ( e != h1);
                            verr << " (closed cycle) and";
                            e = hprime;
                            do {
                                if ( ! e->is_border())
                                    verr << " " << find_facet(
                                                decorator.get_facet(e));
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
    if( h1->vertex() == &(index_to_vertex_map[v1]))
        set_vertex_to_edge_map( v1, h1);
    CGAL_assertion( h1->vertex() == &(index_to_vertex_map[v1]));
    h1 = h2;
    v1 = v2;
}

template < class HDS>
void
Polyhedron_incremental_builder_3<HDS>::
end_facet() {
    if ( m_error)
        return;
    CGAL_assertion( check_protocoll == 2);
    CGAL_assertion( ! first_vertex);
    // cleanup all static status variables
    add_vertex_to_facet( w1);
    last_vertex = true;
    add_vertex_to_facet( w2);
    CGAL_assertion( check_protocoll == 2);
    CGAL_assertion_code( check_protocoll = 1;)
    Halfedge_data_structure_decorator<HDS> decorator;
    decorator.set_facet_halfedge( current_facet,
                                  get_vertex_to_edge_map(w1));
    ++new_facets;
}

template < class HDS>  CGAL_MEDIUM_INLINE
void
Polyhedron_incremental_builder_3<HDS>::
end_surface() {
    if ( m_error)
        return;
    CGAL_assertion( check_protocoll == 1);
    CGAL_assertion_code( check_protocoll = 0;)
}

template < class HDS>
bool
Polyhedron_incremental_builder_3<HDS>::
check_unconnected_vertices() {
    if ( m_error)
        return false;
    bool unconnected = false;
    Verbose_ostream verr( m_verbose);
    for( std::size_t i = 0; i < new_vertices; i++) {
        if( get_vertex_to_edge_map( i) == NULL) {
            verr << "Polyhedron_incremental_builder_3<HDS>::\n"
                 << "check_unconnected_vertices( verb = true): "
                 << "vertex " << i << " is unconnected." << std::endl;
            unconnected = true;
        }
    }
    return unconnected;
}

template < class HDS>
bool
Polyhedron_incremental_builder_3<HDS>::
remove_unconnected_vertices( Tag_true) {
    if ( m_error)
        return true;
    for( std::size_t i = 0; i < new_vertices; i++) {
        if( get_vertex_to_edge_map( i) == NULL) {
            hds.delete_vertex( &index_to_vertex_map[i]);
        }
    }
    return true;
}

CGAL_END_NAMESPACE
#endif // CGAL_POLYHEDRON_OLD_POLYHEDRON_INCREMENTAL_BUILDER_3_H //
// EOF //
