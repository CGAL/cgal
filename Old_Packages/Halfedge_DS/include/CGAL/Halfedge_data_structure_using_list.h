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
// file          : Halfedge_data_structure_using_list.h
// chapter       : $CGAL_Chapter: Halfedge Data Structures $
// package       : $CGAL_Package: Halfedge_DS 2.8 (13 Sep 2000) $
// source        : hds.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : MPI Saarbruecken (Stefan Schirra <stschirr@mpi-sb.mpg.de>)
//
// Halfedge Data Structure Using a List Implementation.
// ============================================================================

#ifndef CGAL_HALFEDGE_DATA_STRUCTURE_USING_LIST_H
#define CGAL_HALFEDGE_DATA_STRUCTURE_USING_LIST_H 1
#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif // CGAL_BASIC_H
#ifndef CGAL_PROTECT_CSTDDEF
#include <cstddef>
#define CGAL_PROTECT_CSTDDEF
#endif
#ifndef CGAL_CIRCULATOR_H
#include <CGAL/circulator.h>
#endif
#ifndef CGAL_IN_PLACE_LIST_H
#include <CGAL/In_place_list.h>
#endif
#ifndef CGAL_N_STEP_ADAPTOR_H
#include <CGAL/N_step_adaptor.h>
#endif
#ifndef CGAL_HALFEDGE_DATA_STRUCTURE_DECORATOR_H
#include <CGAL/Halfedge_data_structure_decorator.h>
#endif
#ifndef CGAL_PROTECT_MAP
#include <map>
#define CGAL_PROTECT_MAP
#endif
#include <memory>

// Define shorter names to please linker (g++/egcs)
#define _HDS_list_vertex                    _Hlv
#define _HDS_list_halfedge                  _Hlh
#define _HDS_list_facet                     _Hlf

CGAL_BEGIN_NAMESPACE


template <class V, class H, class F> class _HDS_list_vertex;
template <class V, class H, class F> class _HDS_list_halfedge;
template <class V, class H, class F> class _HDS_list_facet;

template <class V, class H, class F>
class _HDS_list_vertex
    : public  V,
      public  In_place_list_base< _HDS_list_vertex<V,H,F> >
{
public:
    typedef V                                   Base;
    typedef _HDS_list_vertex<V,H,F>             Vertex;
    typedef _HDS_list_halfedge<V,H,F>           Halfedge;
    typedef _HDS_list_facet<V,H,F>              Facet;

    // Needed for point constructor, which is needed for efficiency.
    typedef typename V::Point                   Point;

    _HDS_list_vertex() {}
    _HDS_list_vertex( const Point& p) : V(p) {}

    Halfedge*       halfedge()       {return (Halfedge*)(V::halfedge());}
    const Halfedge* halfedge() const {return (const Halfedge*)(V::halfedge());}
    void            set_halfedge( Halfedge* h) { V::set_halfedge(h);}
};

template <class V, class H, class F>
class _HDS_list_halfedge
    : public  H,
      public  In_place_list_base< _HDS_list_halfedge<V,H,F> >
{
public:
    typedef H                                   Base;
    typedef _HDS_list_vertex<V,H,F>             Vertex;
    typedef _HDS_list_halfedge<V,H,F>           Halfedge;
    typedef _HDS_list_facet<V,H,F>              Facet;

    typedef typename H::Supports_halfedge_prev    Supports_halfedge_prev;
    typedef typename H::Supports_halfedge_vertex  Supports_halfedge_vertex;
    typedef typename H::Supports_halfedge_facet   Supports_halfedge_facet;

    Halfedge*       opposite()       {return (Halfedge*)(H::opposite());}
    Halfedge*       next()           {return (Halfedge*)(H::next());}
    Halfedge*       prev()           {return (Halfedge*)(H::prev());}
    Vertex*         vertex()         {return (Vertex*)(H::vertex());}
    Facet*          facet()          {return (Facet*)(H::facet());}

    const Halfedge* opposite() const {return (const Halfedge*)(H::opposite());}
    const Halfedge* next()     const {return (const Halfedge*)(H::next());}
    const Halfedge* prev()     const {return (const Halfedge*)(H::prev());}
    const Vertex*   vertex()   const {return (const Vertex*)(H::vertex());}
    const Facet*    facet()    const {return (const Facet*)(H::facet());}

    void  set_next( Halfedge* h)     { H::set_next(h);}
    void  set_prev( Halfedge* h)     { H::set_prev(h);}
    void  set_vertex( Vertex* ve)    { H::set_vertex(ve);}
    void  set_facet( Facet* facet)   { H::set_facet(facet);}

private:
    void  set_opposite( void* h)     { H::set_opposite(h);}
};


template <class V, class H, class F>
class _HDS_list_facet
    : public  F,
      public  In_place_list_base< _HDS_list_facet<V,H,F> >
{
public:
    typedef F                                   Base;
    typedef _HDS_list_vertex<V,H,F>             Vertex;
    typedef _HDS_list_halfedge<V,H,F>           Halfedge;
    typedef _HDS_list_facet<V,H,F>              Facet;

    Halfedge*       halfedge()       {return (Halfedge*)(F::halfedge());}
    const Halfedge* halfedge() const {return (const Halfedge*)(F::halfedge());}
    void            set_halfedge( Halfedge* h) { F::set_halfedge(h);}
};


// A Halfedge Data Structure Using Lists
// -----------------------------------------
//
// DEFINITION
//
// The class Halfedge_data_structure_using_list<Vertex,Halfedge,Facet>
// is a halfedge data structure parameterized with vertex, halfedge,
// and facet types. The base classes defined in the previous subsection
// could be used therefore. It is sufficient for the parameter classes to
// implement the pointers as `void*'. They do not have to know the types
// of their relatives.
//
// Halfedge_data_structure_using_list<Vertex,Halfedge,Facet> uses a list
// implementation and supports therefore removal, but the iterators are
// only bidirectional iterators. The capacity is not restricted and calls
// to reserve do not invalidate any iterator or circulator.

template < class V, class H, class F>
class Halfedge_data_structure_using_list {
public:
    typedef Halfedge_data_structure_using_list<V,H,F>   Self;
    typedef _HDS_list_vertex<V,H,F>         Vertex;
    typedef _HDS_list_halfedge<V,H,F>       Halfedge;
    typedef _HDS_list_facet<V,H,F>          Facet;

    // Point needed for Vertex constructor for efficiency reasons.
    typedef typename Vertex::Point          Point;

    // Halfegdes are allocated in pairs. Here is the type for that.
    struct Halfedge_pair {
        Halfedge first;
        Halfedge second;
        Halfedge_pair() {}
        Halfedge_pair( const Halfedge& h, const Halfedge& g)
            : first(h), second(g) {}
    };

    typedef std::allocator<Vertex>          Vertex_allocator;
    typedef std::allocator<Halfedge_pair>   Edge_allocator;
    typedef std::allocator<Facet>           Facet_allocator;
protected:


    // Three in-place lists for the elements. They are unmanaged.
    typedef In_place_list<Vertex,  false, Vertex_allocator> Vertex_list;
    typedef In_place_list<Halfedge,false, std::allocator<Halfedge> >
                                                            Halfedge_list;
    typedef In_place_list<Facet,   false, Facet_allocator>  Facet_list;

public:
    typedef typename Halfedge_list::size_type   Size;
    // typedef typename Halfedge_list::size_type   size_type;
    typedef typename Halfedge_list::difference_type
                                                Difference;
protected:
    Vertex_list    vertices;
    Halfedge_list  halfedges;
    Facet_list     facets;
    Size           nb_border_halfedges;
    Size           nb_border_edges;
    Halfedge*      border_halfedges;

    Vertex_allocator vertex_allocator;
    Edge_allocator   edge_allocator;
    Facet_allocator  facet_allocator;

    Vertex* get_vertex_node( const Vertex& t) {
        Vertex* p = vertex_allocator.allocate(1);
        vertex_allocator.construct(p, t);
        return p;
    }
    void put_vertex_node( Vertex* p) {
        vertex_allocator.destroy( p);
        vertex_allocator.deallocate( p, 1);
    }

    Halfedge* get_edge_node( const Halfedge& h, const Halfedge& g) {
        // creates a new pair of opposite border halfedges.
        Halfedge_pair* hpair = edge_allocator.allocate(1);
        edge_allocator.construct(hpair, Halfedge_pair( h, g));
        Halfedge* h2 = &(hpair->first);
        Halfedge* g2 = &(hpair->second);
        CGAL_assertion( h2 == (Halfedge*)hpair);
        CGAL_assertion( g2 == h2 + 1);
        h2->H::set_opposite(g2);
        g2->H::set_opposite(h2);
        return h2;
    }
    void put_edge_node( Halfedge* h) {
        // deletes the pair of opposite halfedges h and h->opposite().
        Halfedge* g = &* (h->opposite());
        Halfedge_pair* hpair = (Halfedge_pair*)h;
        if ( h > g)
            hpair = (Halfedge_pair*)g;
        CGAL_assertion( &(hpair->first) == (Halfedge*)hpair);
        edge_allocator.destroy( hpair);
        edge_allocator.deallocate( hpair, 1);
    }

    Facet* get_facet_node( const Facet& t) {
        Facet* p = facet_allocator.allocate(1);
        facet_allocator.construct(p, t);
        return p;
    }
    void put_facet_node( Facet* p) {
        facet_allocator.destroy( p);
        facet_allocator.deallocate( p, 1);
    }


public:
    typedef typename Vertex::Supports_vertex_halfedge Supports_vertex_halfedge;
    typedef typename Halfedge::Supports_halfedge_prev Supports_halfedge_prev;
    typedef typename Halfedge::Supports_halfedge_vertex
                                                      Supports_halfedge_vertex;
    typedef typename Halfedge::Supports_halfedge_facet
                                                      Supports_halfedge_facet;
    typedef typename Facet::Supports_facet_halfedge   Supports_facet_halfedge;

    typedef typename Vertex::Supports_vertex_point    Supports_vertex_point;
    typedef typename Facet::Supports_facet_plane      Supports_facet_plane;
    typedef typename Facet::Supports_facet_normal     Supports_facet_normal;

    typedef Tag_true                                  Supports_removal;
    typedef std::bidirectional_iterator_tag           iterator_category;

    typedef typename Vertex_list::iterator      Vertex_iterator;
    typedef typename Halfedge_list::iterator    Halfedge_iterator;
    typedef typename Facet_list::iterator       Facet_iterator;
    typedef N_step_adaptor< Halfedge_iterator,
                                2,
                                Halfedge&,
                                Halfedge*,
                                Halfedge,
                                std::ptrdiff_t,
                                iterator_category>
                                                Edge_iterator;

    typedef typename Vertex_list::const_iterator    Vertex_const_iterator;
    typedef typename Halfedge_list::const_iterator  Halfedge_const_iterator;
    typedef typename Facet_list::const_iterator     Facet_const_iterator;
    typedef N_step_adaptor< Halfedge_const_iterator,
                                2,
                                const Halfedge&,
                                const Halfedge*,
                                Halfedge,
                                std::ptrdiff_t,
                                iterator_category>
                                                Edge_const_iterator;


// CREATION

private:
    void pointer_update( const Self& hds);
        // Update own pointers assuming that they lived previously
        // in a halfedge data structure `hds' with lists.

public:
    Halfedge_data_structure_using_list()
        : nb_border_halfedges(0), nb_border_edges(0), border_halfedges(NULL) {}
        // the empty polyhedron `P'.

    Halfedge_data_structure_using_list( Size, Size, Size)
        : nb_border_halfedges(0), nb_border_edges(0), border_halfedges(NULL) {}
        // Parameter order is v,h,f.
        // a polyhedron `P' with storage reserved for v vertices, h
        // halfedges, and f facets. The reservation sizes are a hint for
        // optimizing storage allocation. They are not used here.


    void reserve( Size, Size, Size) {}
        // Parameter order is v,h,f.
        // reserve storage for v vertices, h halfedges, and f facets. The
        // reservation sizes are a hint for optimizing storage allocation.
        // If the `capacity' is already greater than the requested size
        // nothing happens. If the `capacity' changes all iterators and
        // circulators invalidates. Function is void here.

    ~Halfedge_data_structure_using_list() { delete_all(); }

    Halfedge_data_structure_using_list( const Self& hds)
    :  vertices( hds.vertices),
       //halfedges( hds.halfedges): need copying pairs instead
       facets( hds.facets),
       nb_border_halfedges( hds.nb_border_halfedges),
       nb_border_edges( hds.nb_border_edges),
       border_halfedges( hds.border_halfedges)
    {
	// goal is halfedges = hds.halfedges, but we have pairs here
	Halfedge_const_iterator i = hds.halfedges.begin();
	for ( ; i != hds.halfedges.end(); ++ ++ i) {
	    new_edge( & * i);
	}
        pointer_update( hds);
    }

    Self& operator=( const Self& hds)  {
        if ( this != &hds) {
            delete_all();
            vertices            = hds.vertices;
            // goal is halfedges = hds.halfedges, but we have pairs here
	    halfedges = Halfedge_list();
	    Halfedge_const_iterator i = hds.halfedges.begin();
	    for ( ; i != hds.halfedges.end(); ++ ++ i) {
		new_edge( & * i);
	    }
            facets              = hds.facets;
            nb_border_halfedges = hds.nb_border_halfedges;
            nb_border_edges     = hds.nb_border_edges;
            border_halfedges    = hds.border_halfedges;
            pointer_update( hds);
        }
        return *this;
    }

// Access Member Functions

    Size size_of_vertices() const  { return vertices.size();}
    Size size_of_halfedges() const { return halfedges.size();}
        // number of all halfedges (including border halfedges).
    Size size_of_facets() const    { return facets.size();}

    Size capacity_of_vertices() const    { return vertices.max_size();}
    Size capacity_of_halfedges() const   { return halfedges.max_size();}
    Size capacity_of_facets() const      { return facets.max_size();}

    std::size_t bytes() const {
        return sizeof(Self)
               + vertices.size()  * sizeof( Vertex)
               + halfedges.size() * sizeof( Halfedge)
               + facets.size()    * sizeof( Facet);
    }
    std::size_t bytes_reserved() const { return bytes();}

    Vertex_iterator   vertices_begin()   { return vertices.begin();}
    Vertex_iterator   vertices_end()     { return vertices.end();}
    Halfedge_iterator halfedges_begin()  { return halfedges.begin();}
    Halfedge_iterator halfedges_end()    { return halfedges.end();}
    Facet_iterator    facets_begin()     { return facets.begin();}
    Facet_iterator    facets_end()       { return facets.end();}

    Edge_iterator edges_begin() { return Edge_iterator(halfedges_begin());}
        // iterator over all edges. The iterator refers to halfedges, but
        // enumerates only one of the two corresponding opposite
        // halfedges.

    Edge_iterator edges_end() { return Edge_iterator(halfedges_end());}
        // end of the range over all edges.

    // The constant iterators and circulators.

    Vertex_const_iterator   vertices_begin()  const{ return vertices.begin();}
    Vertex_const_iterator   vertices_end()    const{ return vertices.end();}
    Halfedge_const_iterator halfedges_begin() const{ return halfedges.begin();}
    Halfedge_const_iterator halfedges_end()   const{ return halfedges.end();}
    Facet_const_iterator    facets_begin()    const{ return facets.begin();}
    Facet_const_iterator    facets_end()      const{ return facets.end();}

    Edge_const_iterator edges_begin() const {
        return Edge_const_iterator(halfedges_begin());
    }
    Edge_const_iterator edges_end() const {
        return Edge_const_iterator(halfedges_end());
    }

// Insertion
//
// The following operations simply allocate a new element of that type.
// Halfedges are always allocated in pairs of opposite halfedges. The
// opposite pointers are automatically set.

    Vertex* new_vertex() {
                Vertex* v = get_vertex_node( Vertex());
                vertices.push_back( *v);
                return v;
    }

    Vertex* new_vertex( const Vertex* w) {
                Vertex* v = get_vertex_node( *w);
                vertices.push_back( *v);
                return v;
    }

    Vertex* new_vertex( const Point& p) {
                Vertex* v = get_vertex_node( Vertex(p));
                vertices.push_back( *v);
                return v;
    }

    Halfedge* new_edge() {
                // creates a new pair of opposite border halfedges.
                Halfedge* h = get_edge_node( Halfedge(), Halfedge());
                Halfedge* g = h->opposite();
                CGAL_assertion( g == h + 1);
                halfedges.push_back( *h);
                halfedges.push_back( *g);
                return h;
    }

    Halfedge* new_edge( const Halfedge* he) {
                Halfedge* h = get_edge_node(*he,*(he->opposite()));
                Halfedge* g = h->opposite();
                CGAL_assertion( g == h + 1);
                halfedges.push_back( *h);
                halfedges.push_back( *g);
                return h;
    }

    Facet* new_facet(){
                Facet* f = get_facet_node( Facet());
                facets.push_back( *f);
                return f;
    }

    Facet* new_facet( const Facet* g) {
                Facet* f = get_facet_node( *g);
                facets.push_back( *f);
                return f;
    }

// Removal
//
// The following operations erase an element referenced by a pointer.
// Halfedges are always deallocated in pairs of opposite halfedges. Erase
// of single elements is optional. The deletion of all is mandatory.

    void delete_vertex( Vertex* v) {
                vertices.erase(v);
                put_vertex_node( v);
    }

    void delete_edge( Halfedge* h) {
                // deletes the pair of opposite halfedges h.
                Halfedge* g = h->opposite();
                halfedges.erase(h);
                halfedges.erase(g);
                put_edge_node(h);
    }

    void delete_facet( Facet* f) {
                facets.erase(f);
                put_facet_node( f);
    }

    void delete_all() {
                vertices.destroy();
                while ( halfedges.begin() != halfedges.end()) {
                    Halfedge* h = &* (halfedges.begin());
                    halfedges.pop_front();
                    halfedges.pop_front();
                    put_edge_node( h);
                }
                facets.destroy();
                nb_border_halfedges = 0;
                nb_border_edges = 0;
                border_halfedges = NULL;
    }

    void vertex_pop_back() {
                Vertex* v = & (vertices.back());
                vertices.pop_back();
                put_vertex_node(v);
    }

    void edge_pop_back() {
                Halfedge_iterator i = halfedges.end();
                --i;
                --i;
                Halfedge* h = &*i;
                halfedges.pop_back();
                halfedges.pop_back();
                put_edge_node(h);
    }

    void facet_pop_back() {
                Facet* f = & (facets.back());
                facets.pop_back();
                put_facet_node(f);
    }

// Operations with Border Halfedges

    Size size_of_border_halfedges() const { return nb_border_halfedges;}
        // number of border halfedges. An edge with no incident facet
        // counts as two border halfedges. Precondition: `normalize_border()'
        // has been called and no halfedge insertion or removal and no
        // change in border status of the halfedges have occured since
        // then.

    Size size_of_border_edges() const { return nb_border_edges;}
        // number of border edges. If `size_of_border_edges() ==
        // size_of_border_halfedges()' all border edges are incident to a
        // facet on one side and to a hole on the other side.
        // Precondition: `normalize_border()' has been called and no
        // halfedge insertion or removal and no change in border status of
        // the halfedges have occured since then.

    Halfedge_iterator border_halfedges_begin() {
        // halfedge iterator starting with the border edges. The range [
        // `halfedges_begin(), border_halfedges_begin()') denotes all
        // non-border edges. The range [`border_halfedges_begin(),
        // halfedges_end()') denotes all border edges. Precondition:
        // `normalize_border()' has been called and no halfedge insertion
        // or removal and no change in border status of the halfedges have
        // occured since then.
        CGAL_assertion( border_halfedges);
        return Halfedge_iterator( border_halfedges);
    }

    Edge_iterator border_edges_begin()  {
        // ... trial to make Edge_iterator obsolete.
        return Edge_iterator(border_halfedges_begin());
    }

    Halfedge_const_iterator border_halfedges_begin() const {
        return Halfedge_const_iterator( border_halfedges);
    }

    Edge_const_iterator border_edges_begin()  const {
        return Edge_const_iterator(border_halfedges_begin());
    }

    void normalize_border();
        // sorts halfedges such that the non-border edges precedes the
        // border edges. For each border edge that is incident to a facet
        // the halfedge iterator will reference the halfedge incident to
        // the facet right before the halfedge incident to the hole.
};


template < class V, class H, class F>
void
Halfedge_data_structure_using_list<V,H,F>::
pointer_update( const Halfedge_data_structure_using_list<V,H,F>& hds) {
    // Update own pointers assuming that they lived previously
    // in a halfedge data structure `hds' with lists.
    typedef std::map< const Vertex*, Vertex*, std::less<const Vertex*> > V_map;
    typedef std::map< const Halfedge*, Halfedge*,
                                            std::less<const Halfedge*> > H_map;
    typedef std::map< const Facet*, Facet*, std::less<const Facet*> >    F_map;
    V_map v_map;
    H_map h_map;
    F_map f_map;
    // initialize maps.
    Halfedge_iterator ii = halfedges.begin();
    Halfedge_const_iterator i = hds.halfedges.begin();
    for ( ; i != hds.halfedges.end(); ++i) {
        h_map[&*i] = &*ii;
        ++ii;
    }
    h_map[&*i] = &*ii;
    if ( check_tag( Supports_halfedge_vertex())) {
        Vertex_iterator vv = vertices.begin();
        for ( Vertex_const_iterator v = hds.vertices.begin();
              v != hds.vertices.end(); ++v) {
            v_map[&*v] = &*vv;
            ++vv;
        }
    }
    if ( check_tag( Supports_halfedge_facet())) {
        Facet_iterator ff = facets.begin();
        for ( Facet_const_iterator f = hds.facets.begin();
              f != hds.facets.end(); ++f) {
            f_map[&*f] = &*ff;
            ++ff;
        }
    }

    Halfedge_data_structure_decorator<Self> D;
    for ( Halfedge_iterator h = halfedges.begin(); h != halfedges.end(); ++h) {
        h->set_next( h_map[ h->next()]);
	// Superfluous and false: opposite pointer get set upon creation
        //h->H::set_opposite( h_map[ h->opposite()]);
        if ( check_tag( Supports_halfedge_prev()))
            D.set_prev( &*h, h_map[ D.get_prev(&*h)]);
        if ( check_tag( Supports_halfedge_vertex())) {
            D.set_vertex( &*h, v_map[ D.get_vertex(&*h)]);
            D.set_vertex_halfedge( &*h);
        }
        if ( h->is_border())
            D.set_facet( &*h, 0);
        else if ( check_tag( Supports_halfedge_facet())) {
            D.set_facet( &*h, f_map[ D.get_facet(&*h)]);
            D.set_facet_halfedge( &*h);
        }
    }
    border_halfedges = h_map[ border_halfedges];
}

template < class V, class H, class F>
void
Halfedge_data_structure_using_list<V,H,F>::normalize_border() {
    CGAL_assertion_code( std::size_t count = halfedges.size();)
    nb_border_halfedges = 0;
    nb_border_edges = 0;
    Halfedge_list  border;
    Halfedge_iterator i = halfedges_begin();
    while ( i != halfedges_end()) {
        Halfedge_iterator j = i;
        ++i;
        ++i;
        Halfedge_iterator k = j;
        ++k;
        if ( j->is_border()) {
            nb_border_halfedges++;
            nb_border_edges++;
            if (k->is_border())
                nb_border_halfedges++;
            border.splice( border.end(), halfedges, k);
            border.splice( border.end(), halfedges, j);
        } else if ( k->is_border()) {
            nb_border_halfedges++;
            nb_border_edges++;
            border.splice( border.end(), halfedges, j);
            border.splice( border.end(), halfedges, k);
        } else {
            CGAL_assertion_code( count -= 2;)
        }
    }
    CGAL_assertion( count == 2 * nb_border_edges );
    CGAL_assertion( count == border.size());
    if ( i == halfedges_begin()) {
        halfedges.splice( halfedges.end(), border);
        i = halfedges_begin();
    } else {
        --i;
        --i;
        CGAL_assertion( ! (*i).is_border() &&
                   ! (*i).opposite()->is_border());
        halfedges.splice( halfedges.end(), border);
        ++i;
        ++i;
    }
    CGAL_assertion( i == halfedges_end() ||
               (*i).opposite()->is_border());
    border_halfedges = &*i;
}

CGAL_END_NAMESPACE

// Undef shorter names (g++/egcs)
//#undef _HDS_list_vertex
//#undef _HDS_list_halfedge
//#undef _HDS_list_facet
#endif // CGAL_HALFEDGE_DATA_STRUCTURE_USING_LIST_H //
// EOF //
