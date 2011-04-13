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
// file          : Halfedge_data_structure_decorator.h
// chapter       : $CGAL_Chapter: Halfedge Data Structures $
// package       : $CGAL_Package: Halfedge_DS 2.8 (13 Sep 2000) $
// source        : hds_decorator.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : MPI Saarbruecken (Stefan Schirra <stschirr@mpi-sb.mpg.de>)
//
// Halfedge Data Structure Decorator.
// ============================================================================

#ifndef CGAL_HALFEDGE_DATA_STRUCTURE_DECORATOR_H
#define CGAL_HALFEDGE_DATA_STRUCTURE_DECORATOR_H 1

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif
#ifndef CGAL_PROTECT_VECTOR
#include <vector>
#define CGAL_PROTECT_VECTOR
#endif
#ifndef CGAL_IO_VERBOSE_OSTREAM_H
#include <CGAL/IO/Verbose_ostream.h>
#endif // CGAL_IO_VERBOSE_OSTREAM_H

CGAL_BEGIN_NAMESPACE

template < class _HDS >
class Halfedge_data_structure_decorator {
public:

// The class Halfedge_data_structure_decorator<HDS> provides
// auxiliary functions to examine and modify a halfedge data structure
// with respect to the different capabilities supported by the different
// representations. The functions evaluate the support type tags of
// the halfedge data structure to decide on the actions. If the
// feature is not supported nothing is done. Note that for example
// the creation of new halfedges is mandatory for all halfedge data
// structure and will not be repeated here.

// TYPES
// ----------------------------------
    typedef _HDS                     HDS;
    typedef _HDS                     Halfedge_data_structure;
    typedef typename _HDS::Vertex    Vertex;
    typedef typename _HDS::Halfedge  Halfedge;
    typedef typename _HDS::Facet     Facet;

    // Point needed for Vertex constructor for efficiency reasons.
    typedef typename _HDS::Point     Point;

// The following types are equal to either `Tag_true' or `Tag_false',
// dependant whether the named feature is supported or not.

    typedef typename _HDS::Supports_vertex_point
                                                Supports_vertex_point;
    typedef typename _HDS::Supports_vertex_halfedge
                                                Supports_vertex_halfedge;
    typedef typename _HDS::Supports_halfedge_prev
                                                Supports_halfedge_prev;
    typedef typename _HDS::Supports_halfedge_vertex
                                                Supports_halfedge_vertex;
    typedef typename _HDS::Supports_halfedge_facet
                                                Supports_halfedge_facet;
    typedef typename _HDS::Supports_facet_halfedge
                                                Supports_facet_halfedge;

    typedef typename _HDS::Supports_removal     Supports_removal;


// CREATION
// ----------------------------------

    // Halfedge_data_structure_decorator() {}

// Access Functions
// ----------------------------------

    Halfedge* get_vertex_halfedge( Vertex* v) const {
        // returns the incident halfedge of v if supported, `NULL'
        // otherwise.
        return get_vertex_halfedge( v, Supports_vertex_halfedge());
    }

    Vertex* get_vertex( Halfedge* h) const {
        // returns the incident vertex of h if supported, `NULL'
        // otherwise.
        return get_vertex( h, Supports_halfedge_vertex());
    }

    Halfedge* get_prev( Halfedge* h) const {
        // returns the previous halfedge of h if supported, `NULL'
        // otherwise.
        return get_prev( h, Supports_halfedge_prev());
    }

    Halfedge* find_prev( Halfedge* h) const {
        // returns the previous halfedge of h. Uses the `prev()' method if
        // supported or performs a search around the facet using `next()'.
        return find_prev( h, Supports_halfedge_prev());
    }

    Halfedge* find_prev_around_vertex( Halfedge* h) const {
        // returns the previous halfedge of h. Uses the `prev()' method if
        // supported or performs a search around the vertex using `next()'.
        return find_prev_around_vertex( h, Supports_halfedge_prev());
    }

    Facet* get_facet( Halfedge* h) const {
        // returns the incident facet of h if supported, `NULL' otherwise.
        return get_facet( h, Supports_halfedge_facet());
    }

    Halfedge* get_facet_halfedge( Facet* f) const {
        // returns the incident halfedge of f if supported, `NULL'
        // otherwise.
        return get_facet_halfedge( f, Supports_facet_halfedge());
    }

// Const Access Functions
// ----------------------------------

    const Halfedge* get_vertex_halfedge( const Vertex* v) const {
        // returns the incident halfedge of v if supported, `NULL'
        // otherwise.
        return get_vertex_halfedge( v, Supports_vertex_halfedge());
    }

    const Vertex* get_vertex( const Halfedge* h) const {
        // returns the incident vertex of h if supported, `NULL'
        // otherwise.
        return get_vertex( h, Supports_halfedge_vertex());
    }

    const Halfedge* get_prev( const Halfedge* h) const {
        // returns the previous halfedge of h if supported, `NULL'
        // otherwise.
        return get_prev( h, Supports_halfedge_prev());
    }

    const Halfedge* find_prev( const Halfedge* h) const {
        // returns the previous halfedge of h. Uses the `prev()' method if
        // supported or performs a search around the facet using `next()'.
        return find_prev( h, Supports_halfedge_prev());
    }

    const Halfedge* find_prev_around_vertex( const Halfedge* h) const {
        // returns the previous halfedge of h. Uses the `prev()' method if
        // supported or performs a search around the vertex using `next()'.
        return find_prev_around_vertex( h, Supports_halfedge_prev());
    }

    const Facet* get_facet( const Halfedge* h) const {
        // returns the incident facet of h if supported, `NULL' otherwise.
        return get_facet( h, Supports_halfedge_facet());
    }

    const Halfedge* get_facet_halfedge( const Facet* f) const {
        // returns the incident halfedge of f if supported, `NULL'
        // otherwise.
        return get_facet_halfedge( f, Supports_facet_halfedge());
    }

// Creation of New Elements
// ----------------------------------

    Vertex* new_vertex( HDS& hds) const {
        // returns a new vertex from `hds' if vertices are supported,
        // `NULL' otherwise.
        return new_vertex( hds, Supports_halfedge_vertex());
    }

    Vertex* new_vertex( HDS& hds, const Point& p) const {
        // returns a new vertex from `hds' initialized to p if vertices
        // are supported, `NULL' otherwise.
        return new_vertex( hds, p, Supports_halfedge_vertex());
    }

    Vertex* new_vertex( HDS& hds, const Vertex* v) const {
        // returns a copy of v from `hds' if vertices
        // are supported, `NULL' otherwise.
        return new_vertex( hds, v, Supports_halfedge_vertex());
    }

    Facet* new_facet( HDS& hds) const {
        // returns a new facet from `hds' if facets are supported, `NULL'
        // otherwise.
        return new_facet( hds, Supports_halfedge_facet());
    }

    Facet* new_facet( HDS& hds, const Facet* f) const {
        // returns a copy of f from `hds' if facets
        // are supported, `NULL' otherwise.
        return new_facet( hds, f, Supports_halfedge_facet());
    }

// Creation of New Composed Items

    Halfedge* create_loop(HDS& hds) {
        // returns a halfedge from a newly created loop in `hds'
        // consisting of a single closed edge, one vertex and two facets
        // (if supported respectively).
        Halfedge* h = hds.new_edge();
        h->set_next( h);
        h->opposite()->set_next( h->opposite());
        set_prev( h, h);
        set_prev( h->opposite(), h->opposite());
        set_vertex( h, new_vertex(hds));
        set_vertex( h->opposite(), get_vertex(h));
        set_facet( h, new_facet(hds));
        set_facet( h->opposite(), new_facet(hds));
        set_facet_halfedge( h);
        set_facet_halfedge( h->opposite());
        set_vertex_halfedge( h);
        return h;
    }

    Halfedge* create_segment(HDS& hds) {
        // returns a halfedge from a newly created segment in `hds'
        // consisting of a single open edge, two vertices and one facet
        // (if supported respectively).
        Halfedge* h = hds.new_edge();
        h->set_next( h->opposite());
        h->opposite()->set_next( h);
        set_prev( h, h->opposite());
        set_prev( h->opposite(), h);
        set_vertex( h, new_vertex(hds));
        set_vertex( h->opposite(), new_vertex(hds));
        set_facet( h, new_facet(hds));
        set_facet( h->opposite(), get_facet(h));
        set_facet_halfedge( h);
        set_vertex_halfedge( h);
        set_vertex_halfedge( h->opposite());
        return h;
    }

// Removal of Elements
// ----------------------------------

    void delete_vertex(HDS& hds, Vertex* v) {
        // removes the vertex v if vertices are supported.
        delete_vertex( hds, v, Supports_halfedge_vertex());
    }

    void delete_facet(HDS& hds, Facet* f) {
        // removes the facet f if facets are supported.
        delete_facet( hds, f, Supports_halfedge_facet());
    }

// Modifying Functions (Composed)
// ----------------------------------

    void close_tip( Halfedge* h) const {
        // makes `h->opposite()' the successor of h.
        h->set_next( h->opposite());
        set_prev( h->opposite(), h);
    }

    void close_tip( Halfedge* h, Vertex* v) const {
        // makes `h->opposite()' the successor of h and sets the vertex to
        // v.
        h->set_next( h->opposite());
        set_prev( h->opposite(), h);
        set_vertex( h, v);
    }

    void insert_tip( Halfedge* h, Halfedge* v) const {
        // inserts the tip of the edge h into the halfedges around the
        // vertex pointed to by v. `h->opposite()' is the new successor of
        // v and `h->next()' will be set to `v->next()'. The vertex of h
        // will be the one v points to if vertices are supported.
        h->set_next( v->next());
        v->set_next( h->opposite());
        set_prev( h->next(), h);
        set_prev( h->opposite(), v);
        set_vertex( h, get_vertex( v));
    }

    void remove_tip( Halfedge* h) const {
        // removes `h->next()->opposite()' from the halfedge circle around
        // the vertex pointed to by h. The new successor of h will be
        // `h->next()->opposite()->next()'.
        h->set_next( h->next()->opposite()->next());
        set_prev( h->next(), h);
    }

    void insert_halfedge( Halfedge* h, Halfedge* f) const {
        // inserts the halfedge h between f and `f->next()'. The facet of
        // h will be the one f points to if facets are supported.
        h->set_next( f->next());
        f->set_next( h);
        set_prev( h->next(), h);
        set_prev( h, f);
        set_facet( h, get_facet( f));
    }

    void remove_halfedge( Halfedge* h) const {
        // removes `h->next()' from the halfedge circle around the facet
        // pointed to by h. The new successor of h will be `h->next()
        // ->next()'.
        h->set_next( h->next()->next());
        set_prev( h->next(), h);
    }

    void set_vertex_in_vertex_loop( Halfedge*  , Vertex* ,
                                    Tag_false) {}
    void set_vertex_in_vertex_loop( Halfedge* h, Vertex* v,
                                    Tag_true) {
        CGAL_assertion_code( std::size_t termination_count = 0;)
        Halfedge* end = h;
        do {
            CGAL_assertion( ++termination_count != 0);
            h->set_vertex( v);
            h = h->next()->opposite();
        } while ( h != end);
    }

    void set_vertex_in_vertex_loop( Halfedge* h, Vertex* v) {
        // Loops around the vertex incident to h and sets all vertex
        // pointers to v. Precondition: `h != NULL'.
        CGAL_precondition( h != 0);
        set_vertex_in_vertex_loop( h, v, Supports_halfedge_vertex());
    }

    void set_facet_in_facet_loop( Halfedge*  , Facet*, Tag_false) {}
    void set_facet_in_facet_loop( Halfedge* h, Facet* f, Tag_true) {
        CGAL_assertion_code( std::size_t termination_count = 0;)
        Halfedge* end = h;
        do {
            CGAL_assertion( ++termination_count != 0);
            h->set_facet( f);
            h = h->next();
        } while ( h != end);
    }

    void set_facet_in_facet_loop( Halfedge* h, Facet* f) {
        // Loops around the facet incident to h and sets all facet
        // pointers to f. Precondition: `h != NULL'.
        CGAL_precondition( h != 0);
        set_facet_in_facet_loop( h, f, Supports_halfedge_facet());
    }

// Modifying Functions (Euler Operators)
// ----------------------------------


    Halfedge* split_facet( HDS& hds, Halfedge* h, Halfedge* g) {
        // split the facet incident to `h' and `g' into two facets with a
        // new diagonal between the two vertices denoted by `h' and `g'
        // respectively. The second (new) facet is a copy of the first
        // facet. It returns the new diagonal. The time is proportional to
        // the distance from `h' to `g' around the facet.
        Halfedge* hnew = hds.new_edge();
        Facet*    fnew = new_facet( hds, get_facet(h));
        insert_tip( hnew, g);
        insert_tip( hnew->opposite(), h);
        set_facet( hnew, get_facet(h));
        set_facet_in_facet_loop( hnew->opposite(), fnew);
        set_facet_halfedge( hnew);
        set_facet_halfedge( hnew->opposite());
        return hnew;
    }

    Halfedge* join_facet( HDS& hds, Halfedge* h) {
        // join the two facets incident to h. The facet incident to
        // `h->opposite()' gets removed. Both facets might be holes.
        // Returns the predecessor of h. The invariant `join_facet(
        // split_facet( h, g))' returns h and keeps `hds'
        // unchanged. The time is proportional to the size of the facet
        // removed and the time to compute `h.prev()'. Precondition:
        // `HDS' supports removal of facets.
        Assert_compile_time_tag( Tag_true(), Supports_removal());
        Halfedge* hprev = find_prev( h);
        Halfedge* gprev = find_prev( h->opposite());
        remove_tip( hprev);
        remove_tip( gprev);
        hds.delete_edge( h);
        if ( get_facet( gprev) != 0)
            delete_facet( hds, get_facet( gprev));
        h = hprev;
        // 'half' of the halfedges have their correct facets.
        // Here we do the remaining halfedges.
        CGAL_assertion_code( std::size_t termination_count = 0;)
        while ( h != gprev) {
            CGAL_assertion( ++termination_count != 0);
            h = h->next();
            set_facet( h, get_facet( hprev));
        }
        if ( get_facet( hprev) != 0)
            set_facet_halfedge( hprev);
        set_vertex_halfedge( hprev);
        set_vertex_halfedge( gprev);
        return hprev;
    }

    Halfedge* split_vertex( HDS& hds, Halfedge* h, Halfedge* g) {
        // split the vertex incident to `h' and `g' into two vertices and
        // connects them with a new edge. The second (new) vertex is a
        // copy of the first vertex. It returns the new edge. The time is
        // proportional to the distance from `h' to `g' around the vertex.
        Halfedge* hnew = hds.new_edge();
        Vertex*   vnew = new_vertex( hds, get_vertex(h));
        insert_halfedge( hnew, g);
        insert_halfedge( hnew->opposite(), h);
        set_vertex( hnew, get_vertex(h));
        set_vertex_in_vertex_loop( hnew->opposite(), vnew);
        set_vertex_halfedge( hnew);
        set_vertex_halfedge( hnew->opposite());
        return hnew;
    }

    Halfedge* join_vertex( HDS& hds, Halfedge* h) {
        // join the two vertices incident to h. The vertex denoted by
        // `h->opposite()' gets removed. Returns the predecessor of h. The
        // invariant `join_vertex( split_vertex( h, g))' returns h and
        // keeps `hds' unchanged. The time is proportional to the
        // degree of the vertex removed and the time to compute `h.prev()'
        // . Precondition: `HDS' supports removal of vertices.
        Assert_compile_time_tag( Tag_true(), Supports_removal());
        Halfedge* hprev = find_prev( h->opposite());
        Halfedge* gprev = find_prev( h);
        remove_halfedge( hprev);
        remove_halfedge( gprev);
        hds.delete_edge( h);
        delete_vertex( hds, get_vertex( gprev));
        // 'half' of the halfedges have their correct vertex.
        // Here we do the remaining halfedges.
        h = hprev;
        CGAL_assertion_code( std::size_t termination_count = 0;)
        while ( h != gprev) {
            CGAL_assertion( ++termination_count != 0);
            h = h->next()->opposite();
            set_vertex( h, get_vertex( hprev));
        }
        set_vertex_halfedge( hprev);
        if ( ! hprev->is_border())
            set_facet_halfedge(  hprev);
        if ( ! gprev->is_border())
            set_facet_halfedge(  gprev);
        return hprev;
    }

    Halfedge* split_loop( HDS& hds, Halfedge* h, Halfedge* i, Halfedge* j){
        // cuts `hds' into two parts along the cycle (h,i,j).
        // Three copies of the vertices, three copies of the halfedges,
        // and two new triangles will be created. h,i,j will be incident
        // to the first new triangle. The returnvalue will be an halfedge
        // iterator denoting the new halfedge of the second new triangle
        // which was h beforehand. Precondition: h,i,j are distinct,
        // consecutive vertices of `hds' and form a cycle: i.e.
        // `h->vertex() == i->opposite()->vertex()', ..., `j->vertex() ==
        // h->opposite()->vertex()'.
        CGAL_precondition( h != i);
        CGAL_precondition( h != j);
        CGAL_precondition( i != j);
        CGAL_precondition( get_vertex(h) == get_vertex(i->opposite()));
        CGAL_precondition( get_vertex(i) == get_vertex(j->opposite()));
        CGAL_precondition( get_vertex(j) == get_vertex(h->opposite()));
        // Create a copy of the triangle.
        Halfedge* hnew = hds.new_edge( h);
        Halfedge* inew = hds.new_edge( i);
        Halfedge* jnew = hds.new_edge( j);
        close_tip( hnew, new_vertex( hds, get_vertex( h)));
        close_tip( inew, new_vertex( hds, get_vertex( i)));
        close_tip( jnew, new_vertex( hds, get_vertex( j)));
        insert_tip( inew->opposite(), hnew);
        insert_tip( jnew->opposite(), inew);
        insert_tip( hnew->opposite(), jnew);
        // Make the new incidences with the old stucture.
        CGAL_assertion_code( std::size_t termination_count = 0;)
        if ( h->next() != i) {
            Halfedge* g = h->next();
            h->set_next( i);
            set_prev( i, h);
            hnew->set_next( g);
            set_prev( g, hnew);
            g = g->opposite();
            while ( g->next() != i) {
                CGAL_assertion( ++termination_count != 0);
                set_vertex( g, get_vertex( hnew));
                g = g->next()->opposite();
            }
            set_vertex( g, get_vertex( hnew));
            g->set_next( inew);
            set_prev( inew, g);
        }
        if ( i->next() != j) {
            Halfedge* g = i->next();
            i->set_next( j);
            set_prev( j, i);
            inew->set_next( g);
            set_prev( g, inew);
            g = g->opposite();
            while ( g->next() != j) {
                CGAL_assertion( ++termination_count != 0);
                set_vertex( g, get_vertex( inew));
                g = g->next()->opposite();
            }
            set_vertex( g, get_vertex( inew));
            g->set_next( jnew);
            set_prev( jnew, g);
        }
        if ( j->next() != h) {
            Halfedge* g = j->next();
            j->set_next( h);
            set_prev( h, j);
            jnew->set_next( g);
            set_prev( g, jnew);
            g = g->opposite();
            while ( g->next() != h) {
                CGAL_assertion( ++termination_count != 0);
                set_vertex( g, get_vertex( jnew));
                g = g->next()->opposite();
            }
            set_vertex( g, get_vertex( jnew));
            g->set_next( hnew);
            set_prev( hnew, g);
        }
        // Fill the holes with two new facets.
        Facet* f = new_facet(hds);
        set_facet( h, f);
        set_facet( i, f);
        set_facet( j, f);
        set_facet_halfedge( h);
        f = new_facet(hds);
        set_facet( hnew->opposite(), f);
        set_facet( inew->opposite(), f);
        set_facet( jnew->opposite(), f);
        set_facet_halfedge( hnew->opposite());
        // Take care of maybe changed halfedge pointers.
        set_facet_halfedge( hnew);
        set_facet_halfedge( inew);
        set_facet_halfedge( jnew);
        set_vertex_halfedge( hnew);
        set_vertex_halfedge( inew);
        set_vertex_halfedge( jnew);
        return hnew->opposite();
    }

    Halfedge* join_loop( HDS& hds, Halfedge* h, Halfedge* g) {
        // glues the boundary of two facets together. Both facets and the
        // vertices of the facet loop g gets removed. Returns h. The
        // invariant `join_loop( h, split_loop( h, i, j))' returns h and
        // keeps `hds' unchanged. Precondition: `HDS' supports
        // removal of vertices and facets. The facets denoted by h and g
        // are different and have have equal size.
        Assert_compile_time_tag( Tag_true(), Supports_removal());
        CGAL_precondition( get_facet(h) == 0
                           || get_facet(h) != get_facet(g));
        if ( get_facet(h) != 0)
            delete_facet( hds, get_facet(h));
        if ( get_facet(g) != 0)
            delete_facet( hds, get_facet(g));
        Halfedge* hi = h;
        Halfedge* gi = g;
        CGAL_assertion_code( std::size_t termination_count = 0;)
        do {
            CGAL_assertion( ++termination_count != 0);
            Halfedge* hii = hi;
            Halfedge* gii = gi;
            hi = hi->next();
            // gi = find_prev(gi); // Replaced by search around vertex.
            set_facet( hii, get_facet( gii->opposite()));
            set_facet_halfedge( hii);
            delete_vertex( hds, get_vertex( gii->opposite()));
            if ( gii->opposite()->next()->opposite()->next() == gii) {
                gi = gii->opposite()->next()->opposite();
            } else {
                hii->set_next( gii->opposite()->next());
                set_prev( hii->next(), hii);
                gii = gii->opposite()->next()->opposite();
                set_vertex( gii, get_vertex(hii));
                while ( gii->next()->opposite()->next() != gi) {
                    CGAL_assertion( ++termination_count != 0);
                    gii = gii->next()->opposite();
                    set_vertex( gii, get_vertex(hii));
                }
                gi = gii->next()->opposite();
                gii->set_next( hi);
                set_prev( gii->next(), gii);
            }
        } while ( hi != h);
        CGAL_assertion( gi == g);
        do {
            Halfedge* gii = gi;
            gi = gi->next();
            hds.delete_edge( gii);
        } while ( gi != g);
        return h;
    }

  protected:
    // supports facet or not.
    void make_hole( HDS&, Halfedge*, Tag_false) {}
    void fill_hole( HDS&, Halfedge*, Tag_false) {}

    void make_hole( HDS& hds, Halfedge* h, Tag_true) {
        Assert_compile_time_tag( Tag_true(), Supports_removal());
        CGAL_precondition( h != 0);
        CGAL_precondition( ! h->is_border());
        hds.delete_facet( h->facet());
        set_facet_in_facet_loop( h, 0);
    }

    void fill_hole( HDS& hds, Halfedge* h, Tag_true) {
        CGAL_precondition( h != 0);
        CGAL_precondition( h->is_border());
        Facet* f = new_facet(hds);
        set_facet_in_facet_loop( h, f);
        set_facet_halfedge( h);
    }

  public:

    Halfedge* make_hole( HDS& hds, Halfedge* h) {
        // removes incident facet and makes all halfedges incident to the
        // facet to border edges. Returns h. Precondition: `HDS'
        // supports removal of facets. `! h.is_border()'.
        make_hole( hds, h, Supports_halfedge_facet());
        return h;
    }

    Halfedge* fill_hole( HDS& hds, Halfedge* h) {
        // fill a hole with a new created facet. Makes all border
        // halfedges of the hole denoted by h incident to the new facet.
        // Returns h. Precondition: `h.is_border()'.
        fill_hole( hds, h, Supports_halfedge_facet());
        return h;
    }

    Halfedge* add_facet_to_border( HDS& hds, Halfedge* h, Halfedge* g) {
        // creates a new facet from `hds' within the hole incident to h
        // and g by connecting the tip of g with the tip of h with a new
        // halfedge from `hds' and filling this separated part of the hole
        // with a new facet from `hds'. Returns the new halfedge incident
        // to the new facet. Precondition: `h != NULL', `g != NULL',
        // `h->is_border()', `g->is_border()' and g can be reached along
        // the same hole starting with h.
        CGAL_precondition( h != 0);
        CGAL_precondition( g != 0);
        CGAL_precondition( h->is_border());
        CGAL_precondition( g->is_border());
        Halfedge* hh = hds.new_edge();
        insert_tip( hh, h);
        insert_tip( hh->opposite(), g);
        fill_hole( hds, g);
        return hh;
    }

    Halfedge* flip_edge( Halfedge* h) {
        // Left rotate of edge.
        // Precond: both incident facets are triangles.
        CGAL_precondition( h == h->next()->next()->next());
        CGAL_precondition( h->opposite() ==
                    h->opposite()->next()->next()->next());
        Halfedge* hprev = h->next()->next();
        Halfedge* gprev = h->opposite()->next()->next();
        remove_tip( hprev);
        remove_tip( gprev);
        set_facet_halfedge(  hprev);
        set_facet_halfedge(  gprev);
        set_vertex_halfedge( hprev);
        set_vertex_halfedge( gprev);
        set_facet( hprev->next(), hprev->facet());
        set_facet( gprev->next(), gprev->facet());
        hprev = hprev->next();
        gprev = gprev->next();
        insert_tip( h, gprev);
        insert_tip( h->opposite(), hprev);
        CGAL_postcondition( h == h->next()->next()->next());
        CGAL_postcondition( h->opposite() ==
                     h->opposite()->next()->next()->next());
        return h;
    }

// Erasing
// ----------------------------------
  protected:
    // supports facet or not.
    void erase_facet( HDS&,     Halfedge*,   Tag_false) {}
    void erase_facet( HDS& hds, Halfedge* h, Tag_true)  {
        Assert_compile_time_tag( Tag_true(), Supports_removal());
        CGAL_precondition( h != 0);
        CGAL_precondition( ! h->is_border());
        hds.delete_facet( h->facet());
        CGAL_assertion_code( std::size_t termination_count = 0;)
        Halfedge* end = h;
        do {
            CGAL_assertion( ++termination_count != 0);
            set_facet( h, 0);
            Halfedge* g = h->next();
            bool h_tag = h->opposite()->is_border();
            bool g_tag = g->opposite()->is_border();
            if ( h_tag && g_tag && g->opposite()->next() == h->opposite()){
                delete_vertex( hds, get_vertex(h));
                if ( h != end)
                    hds.delete_edge(h);
            } else {
                if ( g_tag) {
                    set_vertex_halfedge(g->opposite()->next()->opposite());
                    remove_tip(h);
                }
                if ( h_tag) {
                    set_vertex_halfedge(h->next()->opposite());
                    remove_tip( find_prev_around_vertex( h->opposite()));
                    if ( h != end)
                        hds.delete_edge(h);
                }
            }
            h = g;
        } while ( h != end);
        if ( h->opposite()->is_border())
            hds.delete_edge(h);
    }

  public:
    void erase_facet( HDS& hds, Halfedge* h) {
        // removes the facet incident to `h' from `hds' and changes all
        // halfedges incident to the facet into border edges or removes
        // them from the polyhedral surface if they were already border
        // edges. See `make_hole(h)' for a more specialized variant.
        // Precondition: `h != NULL'. If facets are supported,
        // `Supports_removal' must be equivalent to `Tag_true'.
        erase_facet( hds, h, Supports_halfedge_facet());
    }

  protected:                               // Supports_halfedge_vertices
      void erase_connected_component_vertex( HDS&, Halfedge*, Tag_false) {}
      void erase_connected_component_vertex( HDS& hds, Halfedge* h,
                                             Tag_true) {
          // Erases the the vertex incident to h and sets all references
          // from halfedges around this vertex to NULL,
          // if the incident vertex handle is not already equal to
          // NULL. It is used to erase vertices as soon
          // as an vertex is encountered in the graph traversal. At this
          // point of the graph traversal the halfedge cycle around the
          // vertex is still closed. Lateron it will be broken.
          if ( h->vertex() != NULL) {
              hds.delete_vertex( h->vertex());
              set_vertex_in_vertex_loop( h, NULL);
          }
      }
      void erase_connected_component_vertex( HDS& hds, Halfedge* h) {
          erase_connected_component_vertex( hds, h,
                                            Supports_halfedge_vertex());
      }

      typedef std::vector<Halfedge*> HVector;
      void erase_connected_component_face_cycle( HDS& hds, Halfedge* h,
                                                 HVector& stack) {
          // Delete incident facet and set all incidences to 0.
          if ( get_facet(h) != 0) {
              delete_facet( hds, get_facet(h));
              set_facet_in_facet_loop( h, 0);
          }
          // Cycle around facet, delete incident vertices, push new
          // edges on the stack and mark edges as visited.
          erase_connected_component_vertex( hds, h);
          Halfedge* g = h->next();
          h->set_next(0);
          while (g != h) {
              erase_connected_component_vertex( hds, g);
              if ( g->opposite()->next() != 0)
                  stack.push_back( g->opposite());
              Halfedge* gg = g->next();
              g->set_next(0);
              g = gg;
          }
      }

  public:
    void erase_connected_component( HDS& hds, Halfedge* h) {
        Assert_compile_time_tag( Tag_true(), Supports_removal());
        HVector stack;
        // Algorithm: The next() pointer is used as visited tag
        //     for a graph search. If the next pointer of an halfedge
        //     or its opposite halfedge is set to 0, this edge has already
        //     been visited and must not be put on the stack again.
        // Initializing: Cycle through the face-cycle of h and put
        //     all opposite halfedges on the stack. Put h->opposite()
        //     on the stack. Note that even if the face cycle of h looks
        //     ugly ( e.g. h->opposite() is also in the cycle), neither
        //     h nor h->opposite() will be put on the stack. If
        //     h->opposite() is in the cycle, when h will be popped from
        //     the stack it will be immediately deleted.
        // Loop invariant: For each edge h on the stack h->opposite()->
        //     next() == NULL.
        // Looping: For each edge h on the stack, if h->next() is
        //     not already equal to NULL, cycle through the face-cycle
        //     of h and put all opposite halfedges on the stack.
        //     Delete h.
        // Where: Cycle through a face means: If h->facet() != NULL
        //     delete h->facet() and set all facet handles to NULL.
        //     Loop through the halfedges g around the facet, call
        //     erase_connected_component_vertex for each g, push
        //     g->opposite() on the stack if g->opposite()->next()
        //     is not already NULL. This implies that h->opposite()
        //     is not put on the stack again.
        erase_connected_component_face_cycle( hds, h, stack);
        stack.push_back( h->opposite());
        while ( ! stack.empty()) {
            h = stack.back();
            stack.pop_back();
            CGAL_assertion( h->opposite()->next() == 0);
            if ( h->next() != 0)
                erase_connected_component_face_cycle( hds, h, stack);
            hds.delete_edge( h);
        }
    }

// Modifying Functions (Primitives)
// ----------------------------------

    void set_point( Vertex* v, const Point& p) const {
        set_point( v, p, Supports_vertex_point());
        // sets the point of v to p.
    }

    void set_point( Halfedge* h, const Point& p) const {
        // sets the point of the vertex incident to h to p.
        set_point( h, p, Supports_halfedge_vertex());
    }

    void set_vertex_halfedge( Vertex* v, Halfedge* g) const {
        // sets the incident halfedge of v to g.
        set_vertex_halfedge( v, g, Supports_vertex_halfedge());
    }

    void set_vertex_halfedge( Halfedge* h) const {
        // sets the incident halfedge of the vertex incident to h to h.
        set_vertex_halfedge( h, h, Supports_halfedge_vertex());
    }

    void set_vertex( Halfedge* h, Vertex* v) const {
        // sets the incident vertex of h to v.
        set_vertex(h, v, Supports_halfedge_vertex());
    }

    void set_prev( Halfedge* h, Halfedge* g) const {
        // sets the previous link of h to g.
        set_prev( h, g, Supports_halfedge_prev());
    }

    void set_facet( Halfedge* h, Facet* f) const {
        // sets the incident facet of h to f.
        set_facet(h, f, Supports_halfedge_facet());
    }

    void set_facet_halfedge( Facet* f, Halfedge* g) const {
        // sets the incident halfedge of f to g.
        set_facet_halfedge( f, g, Supports_facet_halfedge());
    }

    void set_facet_halfedge( Halfedge* h) const {
        // sets the incident halfedge of the facet incident to h to h.
        set_facet_halfedge( h, h, Supports_halfedge_facet());
    }

// Implementing These Functions.
// ====================================================
// Access Functions
// ----------------------------------

    Halfedge* get_vertex_halfedge( Vertex*  ,Tag_false) const {
        return NULL;
    }
    Halfedge* get_vertex_halfedge( Vertex* v,Tag_true) const {
        return v->halfedge();
    }

    Vertex* get_vertex( Halfedge*  , Tag_false) const { return NULL;}
    Vertex* get_vertex( Halfedge* h, Tag_true)  const {
        return h->vertex();
    }

    Halfedge* get_prev(Halfedge*  , Tag_false) const {return NULL;}
    Halfedge* get_prev(Halfedge* h, Tag_true)  const {return h->prev();}

    Halfedge* find_prev(Halfedge* h, Tag_true) const { return h->prev(); }
    Halfedge* find_prev(Halfedge* h, Tag_false) const {
        Halfedge* g = h;
        while ( g->next() != h)
            g = g->next();
        return g;
    }

    Halfedge* find_prev_around_vertex(Halfedge* h, Tag_true) const {
        return h->prev();
    }
    Halfedge* find_prev_around_vertex(Halfedge* h, Tag_false) const {
        Halfedge* g = h->opposite();
        while ( g->next() != h)
            g = g->next()->opposite();
        return g;
    }

    Facet* get_facet(Halfedge*  , Tag_false) const { return NULL;}
    Facet* get_facet(Halfedge* h, Tag_true)  const { return h->facet();}

    Halfedge* get_facet_halfedge( Facet*  , Tag_false) const {
        return NULL;
    }
    Halfedge* get_facet_halfedge( Facet* f, Tag_true) const {
        return f->halfedge();
    }

// Const Access Functions
// ----------------------------------

    const Halfedge*
    get_vertex_halfedge( const Vertex*  ,Tag_false) const { return NULL; }
    const Halfedge*
    get_vertex_halfedge( const Vertex* v,Tag_true) const {
        return v->halfedge();
    }

    const Vertex*
    get_vertex( const Halfedge*  , Tag_false) const { return NULL; }
    const Vertex*
    get_vertex( const Halfedge* h, Tag_true)  const { return h->vertex(); }

    const Halfedge*
    get_prev( const Halfedge*  , Tag_false) const { return NULL; }
    const Halfedge*
    get_prev( const Halfedge* h, Tag_true)  const { return h->prev(); }

    const Halfedge*
    find_prev( const Halfedge* h, Tag_true) const { return h->prev(); }
    const Halfedge*
    find_prev( const Halfedge* h, Tag_false) const {
        const Halfedge* g = h;
        while ( g->next() != h)
            g = g->next();
        return g;
    }

    const Halfedge* find_prev_around_vertex( const Halfedge* h,
                                             Tag_true) const {
        return h->prev();
    }
    const Halfedge* find_prev_around_vertex( const Halfedge* h,
                                             Tag_false) const {
        const Halfedge* g = h->opposite();
        while ( g->next() != h)
            g = g->next()->opposite();
        return g;
    }

    const Facet*
    get_facet( const Halfedge*  , Tag_false) const { return NULL;}
    const Facet*
    get_facet( const Halfedge* h, Tag_true)  const { return h->facet();}

    const Halfedge*
    get_facet_halfedge( const Facet*  , Tag_false) const {
        return NULL;
    }
    const Halfedge*
    get_facet_halfedge( const Facet* f, Tag_true) const {
        return f->halfedge();
    }

// Creation of New Elements
// ----------------------------------

    Vertex* new_vertex( HDS&    , Tag_false) const { return NULL;}
    Vertex* new_vertex( HDS& hds, Tag_true) const  {
        return hds.new_vertex();
    }

    Vertex* new_vertex( HDS&    , const Point&  , Tag_false) const {
        return NULL;
    }
    Vertex* new_vertex( HDS& hds, const Point& p, Tag_true)  const {
        return hds.new_vertex(p);
    }

    Vertex* new_vertex( HDS&    , const Vertex*  ,Tag_false) const {
        return NULL;
    }
    Vertex* new_vertex( HDS& hds, const Vertex* v,Tag_true)  const {
        return hds.new_vertex(v);
    }

    Facet* new_facet( HDS&    , Tag_false) const { return NULL;}
    Facet* new_facet( HDS& hds, Tag_true)  const {
        return hds.new_facet();
    }

    Facet* new_facet( HDS&    , const Facet*  , Tag_false) const {
        return NULL;
    }
    Facet* new_facet( HDS& hds, const Facet* f, Tag_true)  const {
        return hds.new_facet(f);
    }


// Removal of Elements
// ----------------------------------

    void delete_vertex(HDS&    , Vertex*  , Tag_false) {}
    void delete_vertex(HDS& hds, Vertex* v, Tag_true) {
        hds.delete_vertex( v);
    }

    void delete_facet(HDS&    , Facet*  , Tag_false) {}
    void delete_facet(HDS& hds, Facet* f, Tag_true) {
        hds.delete_facet( f);
    }
// Modifying Function Primitives
// ----------------------------------

    void set_point( Vertex*  , const Point&  , Tag_false) const {}
    void set_point( Vertex* v, const Point& p, Tag_true)  const {
        v->point() = p;
    }

    void set_point( Halfedge*  , const Point&  , Tag_false) const {}
    void set_point( Halfedge* h, const Point& p, Tag_true)  const {
        set_point( h->vertex(), p);
    }

    void set_vertex_halfedge( Vertex*  , Halfedge* , Tag_false) const {}
    void set_vertex_halfedge( Vertex* v, Halfedge* g, Tag_true)  const {
        v->set_halfedge(g);
    }

    void set_vertex_halfedge( Halfedge*, Halfedge*, Tag_false) const {}
    void set_vertex_halfedge( Halfedge* h, Halfedge* g,Tag_true) const {
        set_vertex_halfedge( h->vertex(), g);
    }

    void set_vertex( Halfedge*,   Vertex*  , Tag_false) const {}
    void set_vertex( Halfedge* h, Vertex* v, Tag_true)  const {
        h->set_vertex(v);
    }

    void set_prev( Halfedge*  , Halfedge*  , Tag_false) const {}
    void set_prev( Halfedge* h, Halfedge* g, Tag_true)  const {
        h->set_prev( g);
    }

    void set_facet( Halfedge*,   Facet*  , Tag_false) const {}
    void set_facet( Halfedge* h, Facet* f, Tag_true)  const {
        h->set_facet(f);
    }

    void set_facet_halfedge( Facet*  , Halfedge*  , Tag_false) const {}
    void set_facet_halfedge( Facet* f, Halfedge* g, Tag_true)  const {
        f->set_halfedge(g);
    }

    void set_facet_halfedge( Halfedge*  , Halfedge*, Tag_false) const {}
    void set_facet_halfedge( Halfedge* h, Halfedge* g, Tag_true) const {
        set_facet_halfedge( h->facet(), g);
    }

    bool normalized_border_is_valid( const _HDS& hds, bool verb = false)
        const;
        // checks whether all non-border edges precedes the border edges.

    bool is_valid( const _HDS& hds, bool verb = false, int level = 0)
        const;
        // checks the combinatorial consistency.

    void reorient_facet( Halfedge* first);
                  // Supports: halfedge pointer in facets.
    void inside_out( _HDS& hds, Tag_false);
    void inside_out( _HDS& hds, Tag_true);

    void inside_out(_HDS& hds ) {
        inside_out( hds, Supports_facet_halfedge());
    }
};

template < class _HDS >
bool
Halfedge_data_structure_decorator<_HDS>::
normalized_border_is_valid( const _HDS& hds, bool verb) const {
    // checks whether all non-border edges precedes the border edges.
    typedef typename _HDS::Halfedge_const_iterator Halfedge_const_iterator;
    typedef typename _HDS::Size                    Size;
    bool valid = true;
    Verbose_ostream verr(verb);
    verr << "begin Halfedge_data_structure_decorator<HDS>::"
            "normalized_border_is_valid( verb=true):" << std::endl;

    Halfedge_const_iterator e = hds.halfedges_begin();
    Size count = 0;
    while( e != hds.halfedges_end() && ! e->is_border() && !
           e->opposite()->is_border()) {
        ++e;
        ++e;
        ++count;
    }
    verr << "    non-border edges: " << count << std::endl;
    if ( e != hds.border_halfedges_begin()) {
        verr << "    first border edge does not start at "
                "border_halfedges_begin()" << std::endl;
        valid = false;
    } else {
        count = 0;
        while( valid && e != hds.halfedges_end() &&
               e->opposite()->is_border()) {
            ++e;
            ++e;
            ++count;
        }
        verr << "    border     edges: " << count << std::endl;
        verr << "    total      edges: " << hds.size_of_halfedges()/2
             << std::endl;
        if ( e != hds.halfedges_end()) {
            if ( e->is_border()) {
                verr << "    border edge " << count
                     << ": wrong orientation." << std::endl;
            }
            verr << "    the sum of full + border equals not total edges."
                 << std::endl;
            valid = false;
        }
    }
    verr << "end of Halfedge_data_structure_decorator<HDS>::"
            "normalized_border_is_valid(): structure is "
         << ( valid ? "valid." : "NOT VALID.") << std::endl;
    return valid;
}

template < class _HDS >
bool
Halfedge_data_structure_decorator<_HDS>::
is_valid( const _HDS& hds, bool verb, int level) const {
    // checks the combinatorial consistency.
    typedef typename _HDS::Halfedge_const_iterator Halfedge_const_iterator;
    typedef typename _HDS::Vertex_const_iterator   Vertex_const_iterator;
    typedef typename _HDS::Facet_const_iterator    Facet_const_iterator;
    typedef typename _HDS::Size                    Size;
    typedef typename _HDS::Supports_halfedge_prev  Supports_halfedge_prev;
    typedef typename _HDS::Supports_halfedge_vertex
        Supports_halfedge_vertex;
    typedef typename _HDS::Supports_vertex_halfedge
        Supports_vertex_halfedge;
    typedef typename _HDS::Supports_halfedge_facet
        Supports_halfedge_facet;
    typedef typename _HDS::Supports_facet_halfedge
        Supports_facet_halfedge;

    Verbose_ostream verr(verb);
    verr << "begin Halfedge_data_structure_decorator<HDS>::is_valid("
            " verb=true, level = " << level << "):" << std::endl;

    bool valid = ( 1 != (hds.size_of_halfedges() & 1));
    if ( ! valid)
        verr << "number of halfedges is odd." << std::endl;

    // All halfedges.
    Halfedge_const_iterator begin = hds.halfedges_begin();
    Halfedge_const_iterator end   = hds.halfedges_end();
    Size  n = 0;
    Size nb = 0;
    for( ; valid && (begin != end); begin++) {
        verr << "halfedge " << n << std::endl;
        if ( begin->is_border())
            verr << "    is border halfedge" << std::endl;
        // Pointer integrity.
        valid = valid && ( begin->next() != NULL);
        valid = valid && ( begin->opposite() != NULL);
        if ( ! valid) {
            verr << "    pointer integrity corrupted (ptr==NULL)."
                 << std::endl;
            break;
        }
        // opposite integrity.
        valid = valid && ( begin->opposite() != &*begin);
        valid = valid && ( begin->opposite()->opposite() == &*begin);
        if ( ! valid) {
            verr << "    opposite pointer integrity corrupted."
                 << std::endl;
            break;
        }
        // previous integrity.
        valid = valid && ( ! check_tag( Supports_halfedge_prev()) ||
                           get_prev(begin->next()) == &*begin);
        if ( ! valid) {
            verr << "    previous pointer integrity corrupted."
                 << std::endl;
            break;
        }
        if ( level > 0) {
            // vertex integrity.
            valid = valid && ( ! check_tag( Supports_halfedge_vertex())
                               || get_vertex( &*begin) != NULL);
            valid = valid && ( get_vertex( &*begin) ==
                               get_vertex( begin->next()->opposite()));
            if ( ! valid) {
                verr << "    vertex pointer integrity corrupted."
                     << std::endl;
                break;
            }
            // facet integrity.
            valid = valid && ( ! check_tag( Supports_halfedge_facet())
                     || begin->is_border() || get_facet( &*begin) != NULL);
            valid = valid && ( get_facet( &*begin) ==
                        get_facet( begin->next()));
            if ( ! valid) {
                verr << "    facet pointer integrity corrupted."
                     << std::endl;
                break;
            }
        }
        ++n;
        if ( begin->is_border())
            ++nb;
    }
    verr << "summe border halfedges (2*nb) = " << 2 * nb << std::endl;
    if ( valid && n != hds.size_of_halfedges())
        verr << "counting halfedges failed." << std::endl;
    if ( valid && level >= 4 && (nb != hds.size_of_border_halfedges()))
        verr << "counting border halfedges failed." << std::endl;
    valid = valid && ( n  == hds.size_of_halfedges());
    valid = valid && ( level < 4 ||
                       (nb == hds.size_of_border_halfedges()));

    // All vertices.
    Vertex_const_iterator vbegin = hds.vertices_begin();
    Vertex_const_iterator vend   = hds.vertices_end();
    Size v = 0;
    n = 0;
    for( ; valid && (vbegin != vend); ++vbegin) {
        verr << "vertex " << v << std::endl;
        // Pointer integrity.
        if ( get_vertex_halfedge( &*vbegin) != NULL)
            valid = valid && ( ! check_tag(
                Supports_halfedge_vertex()) ||
                get_vertex( get_vertex_halfedge( &*vbegin)) == &*vbegin);
        else
            valid = valid && (! check_tag(
                Supports_vertex_halfedge()));
        if ( ! valid) {
            verr << "    halfedge pointer in vertex corrupted."
                 << std::endl;
            break;
        }
        // cycle-around-vertex test.
        const Halfedge* h = get_vertex_halfedge( &*vbegin);
        if ( h) {
            const Halfedge* g = h;
            do {
                verr << "    halfedge " << n << std::endl;
                ++n;
                h = h->next()->opposite();
                valid = valid && ( n <= hds.size_of_halfedges() && n != 0);
            } while ( valid && (h != g));
        }
        ++v;
    }
    if ( valid && v != hds.size_of_vertices())
        verr << "counting vertices failed." << std::endl;
    if ( valid && level >= 2 && ( check_tag( Supports_vertex_halfedge())
         && n  != hds.size_of_halfedges()))
        verr << "counting halfedges via vertices failed." << std::endl;
    valid = valid && ( v == hds.size_of_vertices());
    valid = valid && ( level < 2 ||
                       ! check_tag( Supports_vertex_halfedge()) ||
                       n  == hds.size_of_halfedges());

    // All facets.
    Facet_const_iterator fbegin = hds.facets_begin();
    Facet_const_iterator fend   = hds.facets_end();
    Size f = 0;
    n = 0;
    for( ; valid && (fbegin != fend); ++fbegin) {
        verr << "facet " << f << std::endl;
        // Pointer integrity.
        if ( get_facet_halfedge( &*fbegin) != NULL)
            valid = valid && ( ! check_tag(
                Supports_halfedge_facet()) ||
                get_facet( get_facet_halfedge( &*fbegin)) == &*fbegin);
        else
            valid = valid && (! check_tag(
                Supports_facet_halfedge()) || begin->is_border());
        if ( ! valid) {
            verr << "    halfedge pointer in facet corrupted."
                 << std::endl;
            break;
        }
        // cycle-around-facet test.
        const Halfedge* h = get_facet_halfedge( &*fbegin);
        if ( h) {
            const Halfedge* g = h;
            do {
                verr << "    halfedge " << n << std::endl;
                ++n;
                h = h->next();
                valid = valid && ( n <= hds.size_of_halfedges() && n != 0);
            } while ( valid && (h != g));
        }
        ++f;
    }
    if ( valid && f != hds.size_of_facets())
        verr << "counting facets failed." << std::endl;
    if ( valid && level >= 3 && check_tag( Supports_facet_halfedge()) &&
         n + nb  != hds.size_of_halfedges())
        verr << "counting halfedges via facets failed." << std::endl;
    valid = valid && ( f == hds.size_of_facets());
    valid = valid && ( level < 3 ||
                       ! check_tag( Supports_facet_halfedge()) ||
                       n + nb  == hds.size_of_halfedges());

    if ( level >= 4) {
        verr << "level 4: normalized_border_is_valid( verbose = true)"
             << std::endl;
        valid = valid && ( normalized_border_is_valid( hds, verb));
    }
    verr << "end of Halfedge_data_structure_decorator<HDS>::"
            "is_valid(): structure is " << ( valid ? "valid." :
            "NOT VALID.") << std::endl;
    return valid;
}


template < class _HDS >  CGAL_MEDIUM_INLINE
void
Halfedge_data_structure_decorator<_HDS>::
reorient_facet( Halfedge* first) {
    if ( first == NULL)
        return;
    Halfedge_data_structure_decorator<_HDS> decorator;
    Halfedge* last  = first;
    Halfedge* prev  = first;
    Halfedge* start = first;
    first = first->next();
    Vertex*  new_v = decorator.get_vertex( start);
    while (first != last) {
        Vertex*  tmp_v = decorator.get_vertex( first);
        decorator.set_vertex( first, new_v);
        decorator.set_vertex_halfedge( first);
        new_v = tmp_v;
        Halfedge* next = first->next();
        first->set_next( prev);
        decorator.set_prev( first, next);
        prev  = first;
        first = next;
    }
    decorator.set_vertex( start, new_v);
    decorator.set_vertex_halfedge( start);
    Halfedge* next = start->next();
    start->set_next( prev);
    decorator.set_prev( start, next);
}

template < class _HDS >
void                            // Supports: halfedge pointer in facets.
Halfedge_data_structure_decorator<_HDS>::
inside_out( _HDS& hds, Tag_false) {
    typedef typename _HDS::Halfedge_iterator Halfedge_iterator;
    Halfedge_iterator begin = hds.halfedges_begin();
    Halfedge_iterator end   = hds.halfedges_end();
    for( ; begin != end; ++begin) {
        // Define the halfedge with the `smallest' pointer
        // value as the one that is used to reorient the facet.
        Halfedge* h = &*begin;
        bool is_min = true;
        Halfedge* c = h;
        Halfedge* d = c;
        do {
            CGAL_assertion( c != NULL);
            if ( h > c)
                is_min = false;
            c = c->next();
        } while ( c != d && is_min);
        if ( is_min)
            reorient_facet( h);
    }
}

template < class _HDS >
void                            // Supports: halfedge pointer in facets.
Halfedge_data_structure_decorator<_HDS>::
inside_out( _HDS& hds, Tag_true) {
    typedef typename _HDS::Halfedge_iterator Halfedge_iterator;
    typedef typename _HDS::Facet_iterator    Facet_iterator;
    Facet_iterator begin = hds.facets_begin();
    Facet_iterator end   = hds.facets_end();
    for( ; begin != end; ++begin)
        reorient_facet( begin->halfedge());
    // Note: A border edge is now parallel to its opposite edge.
    // We scan all border edges for this property. If it holds, we
    // reorient the associated hole and search again until no border
    // edge with that property exists any longer. Then, all holes are
    // reoriented.
    Halfedge_iterator first = hds.halfedges_begin();
    Halfedge_iterator last  = hds.halfedges_end();
    for( ; first != last; ++first) {
        if ( first->is_border() &&
             first->vertex() == first->opposite()->vertex())
            reorient_facet( &*first);
    }
}

CGAL_END_NAMESPACE
#endif // CGAL_HALFEDGE_DATA_STRUCTURE_DECORATOR_H //
// EOF //
