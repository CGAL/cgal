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
// file          : HalfedgeDS_using_list.h
// chapter       : $CGAL_Chapter: Halfedge Data Structures $
// package       : $CGAL_Package: HalfedgeDS 3.3 (27 Sep 2000) $
// source        : hds_list.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : MPI Saarbruecken (Stefan Schirra <stschirr@mpi-sb.mpg.de>)
//
// Halfedge Data Structure Using a List Implementation.
// ============================================================================

#ifndef CGAL_HALFEDGEDS_USING_LIST_H
#define CGAL_HALFEDGEDS_USING_LIST_H 1
#ifndef CGAL_PROTECT_ALGORITHM
#include <algorithm>
#define CGAL_PROTECT_ALGORITHM
#endif
#ifndef CGAL_PROTECT_LIST
#include <list>
#define CGAL_PROTECT_LIST
#endif
#ifndef CGAL_PROTECT_MAP
#include <map>
#define CGAL_PROTECT_MAP
#endif

#ifndef CGAL_HALFEDGEDS_ITERATOR_ADAPTOR_H
#include <CGAL/HalfedgeDS_iterator_adaptor.h>
#endif
#ifndef CGAL_HALFEDGEDS_ITEMS_DECORATOR_H
#include <CGAL/HalfedgeDS_items_decorator.h>
#endif

CGAL_BEGIN_NAMESPACE

#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
template < class Traits_, class HalfedgeDSItems>
class HalfedgeDS_using_list {
public:
    typedef HalfedgeDS_using_list<Traits_,HalfedgeDSItems> Self;
#else
struct HalfedgeDS_using_list {
template < class Traits_, class HalfedgeDSItems>
class HDS {
public:
    typedef HDS<Traits_,HalfedgeDSItems> Self;
#endif
    typedef Traits_                                       Traits;
    typedef HalfedgeDSItems                               Items;

#ifdef __GNUC__
    typedef typename Items::Vertex_wrapper<Self,Traits>   Vertex_wrapper;
    typedef typename Items::Halfedge_wrapper<Self,Traits> Halfedge_wrapper;
    typedef typename Items::Face_wrapper<Self,Traits>     Face_wrapper;
#else // __GNUC__ //
    // here is the standard conforming way
    typedef Items::template Vertex_wrapper<Self,Traits>   Vertex_wrapper;
    typedef Items::template Halfedge_wrapper<Self,Traits> Halfedge_wrapper;
    typedef Items::template Face_wrapper<Self,Traits>     Face_wrapper;
#endif // __GNUC__ //

    typedef typename Vertex_wrapper::Vertex            Vertex;
    typedef std::list<Vertex>                          Vertex_list;
    typedef typename Vertex_list::iterator             Vertex_I;
    typedef typename Vertex_list::const_iterator       Vertex_CI;
    typedef HalfedgeDS_iterator_adaptor<Vertex_I>      Vertex_iterator;
    typedef HalfedgeDS_iterator_adaptor<Vertex_CI>     Vertex_const_iterator;
    typedef Vertex_iterator                            Vertex_handle;
    typedef Vertex_const_iterator                      Vertex_const_handle;

    typedef typename Halfedge_wrapper::Halfedge        Halfedge;
    typedef std::list<Halfedge>                        Halfedge_list;
    typedef typename Halfedge_list::iterator           Halfedge_I;
    typedef typename Halfedge_list::const_iterator     Halfedge_CI;
    typedef HalfedgeDS_iterator_adaptor<Halfedge_I>    Halfedge_iterator;
    typedef HalfedgeDS_iterator_adaptor<Halfedge_CI>   Halfedge_const_iterator;
    typedef Halfedge_iterator                          Halfedge_handle;
    typedef Halfedge_const_iterator                    Halfedge_const_handle;

    typedef typename Face_wrapper::Face                Face;
    typedef std::list<Face>                            Face_list;
    typedef typename Face_list::iterator               Face_I;
    typedef typename Face_list::const_iterator         Face_CI;
    typedef HalfedgeDS_iterator_adaptor<Face_I>        Face_iterator;
    typedef HalfedgeDS_iterator_adaptor<Face_CI>       Face_const_iterator;
    typedef Face_iterator                              Face_handle;
    typedef Face_const_iterator                        Face_const_handle;

    typedef typename Halfedge_list::size_type          size_type;
    typedef typename Halfedge_list::difference_type    difference_type;
    typedef std::bidirectional_iterator_tag            iterator_category;
    typedef Tag_true                                   Supports_removal;

    typedef typename Vertex::Supports_vertex_halfedge Supports_vertex_halfedge;
    typedef typename Halfedge::Supports_halfedge_prev Supports_halfedge_prev;
    typedef typename Halfedge::Supports_halfedge_vertex
                                                      Supports_halfedge_vertex;
    typedef typename Halfedge::Supports_halfedge_face
                                                      Supports_halfedge_face;
    typedef typename Face::Supports_face_halfedge     Supports_face_halfedge;

    static inline Vertex_handle vertex_handle( Vertex* v) {
        Vertex_I vv = 0;
        vv = (Vertex_I)((char*) v - (ptrdiff_t)(&*vv));
        return vv;
    }
    static inline Vertex_const_handle vertex_handle( const Vertex* v) {
        const Vertex_CI vv = 0;
        vv = (Vertex_CI)((const char*) v - (ptrdiff_t)(&*vv));
        return vv;
    }

    static inline Halfedge_handle halfedge_handle( Halfedge* h) {
        Halfedge_I hh = 0;
        char* ptr = (char*) h - (ptrdiff_t)(&*hh);
        hh = (Halfedge_I&)ptr;
        return hh;
    }
    static inline Halfedge_const_handle halfedge_handle( const Halfedge* h) {
        const Halfedge_CI hh = 0;
        hh = (Halfedge_CI)((const char*) h - (ptrdiff_t)(&*hh));
        return hh;
    }

    static inline Face_handle face_handle( Face* f) {
        Face_I ff = 0;
        ff = (Face_I)((char*) f - (ptrdiff_t)(&*ff));
        return ff;
    }
    static inline Face_const_handle face_handle( const Face* f) {
        const Face_CI ff = 0;
        ff = (Face_CI)((const char*) f - (ptrdiff_t)(&*ff));
        return ff;
    }

protected:
    typedef typename Vertex::Base                      VBase;
    typedef typename Halfedge::Base                    HBase;
    typedef typename Halfedge::Base_base               HBase_base;
    typedef typename Face::Base                        FBase;

    Vertex_list        vertices;
    Halfedge_list      halfedges;
    Face_list          faces;

    size_type          nb_border_halfedges;
    size_type          nb_border_edges;
    Halfedge_iterator  border_halfedges;

// CREATION

private:
    void pointer_update( const Self& hds);
        // Update own pointers assuming that they lived previously
        // in a halfedge data structure `hds' with lists.

public:
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
    HalfedgeDS_using_list() : nb_border_halfedges(0), nb_border_edges(0) {}
        // the empty polyhedron `P'.

    HalfedgeDS_using_list( size_type, size_type, size_type)
        : nb_border_halfedges(0), nb_border_edges(0) {}
        // Parameter order is v,h,f.
        // a polyhedron `P' with storage reserved for v vertices, h
        // halfedges, and f faces. The reservation sizes are a hint for
        // optimizing storage allocation. They are not used here.

    HalfedgeDS_using_list( const Self& hds)
#else
    HDS() : nb_border_halfedges(0), nb_border_edges(0) {}
    HDS( size_type, size_type, size_type)
          : nb_border_halfedges(0), nb_border_edges(0) {}
    HDS( const Self& hds)
#endif // CGAL_CFG_NO_TMPL_IN_TMPL_PARAM //
    :  vertices( hds.vertices),
       halfedges( hds.halfedges),
       faces( hds.faces),
       nb_border_halfedges( hds.nb_border_halfedges),
       nb_border_edges( hds.nb_border_edges),
       border_halfedges( hds.border_halfedges)
    {
        pointer_update( hds);
    }

    Self& operator=( const Self& hds)  {
        if ( this != &hds) {
            clear();
            vertices            = hds.vertices;
            halfedges           = hds.halfedges;
            faces               = hds.faces;
            nb_border_halfedges = hds.nb_border_halfedges;
            nb_border_edges     = hds.nb_border_edges;
            border_halfedges    = hds.border_halfedges;
            pointer_update( hds);
        }
        return *this;
    }

    void reserve( size_type, size_type, size_type) {}
        // Parameter order is v,h,f.
        // reserve storage for v vertices, h halfedges, and f faces. The
        // reservation sizes are a hint for optimizing storage allocation.
        // If the `capacity' is already greater than the requested size
        // nothing happens. If the `capacity' changes all iterators and
        // circulators invalidates. The function is void here.

// Access Member Functions

    size_type size_of_vertices()  const { return vertices.size();}
    size_type size_of_halfedges() const { return halfedges.size();}
        // number of all halfedges (including border halfedges).
    size_type size_of_faces()     const { return faces.size();}

    size_type capacity_of_vertices() const    { return vertices.max_size();}
    size_type capacity_of_halfedges() const   { return halfedges.max_size();}
    size_type capacity_of_faces() const       { return faces.max_size();}

    size_t bytes() const {
        return sizeof(Self)
               + vertices.size()  * (sizeof( Vertex)   + 2 * sizeof(void*))
               + halfedges.size() * (sizeof( Halfedge) + 2 * sizeof(void*))
               + faces.size()     * (sizeof( Face)     + 2 * sizeof(void*));
    }
    size_t bytes_reserved() const { return bytes();}

    Vertex_iterator   vertices_begin()   { return vertices.begin();}
    Vertex_iterator   vertices_end()     { return vertices.end();}
    Halfedge_iterator halfedges_begin()  { return halfedges.begin();}
    Halfedge_iterator halfedges_end()    { return halfedges.end();}
    Face_iterator     faces_begin()      { return faces.begin();}
    Face_iterator     faces_end()        { return faces.end();}

    // The constant iterators and circulators.

    Vertex_const_iterator   vertices_begin()  const{ return vertices.begin();}
    Vertex_const_iterator   vertices_end()    const{ return vertices.end();}
    Halfedge_const_iterator halfedges_begin() const{ return halfedges.begin();}
    Halfedge_const_iterator halfedges_end()   const{ return halfedges.end();}
    Face_const_iterator     faces_begin()     const{ return faces.begin();}
    Face_const_iterator     faces_end()       const{ return faces.end();}

// Insertion
//
// The following operations simply allocate a new element of that type.
// Halfedges are always allocated in pairs of opposite halfedges. The
// opposite pointers are automatically set.

    Vertex_handle vertices_push_back( const Vertex& v) {
        vertices.push_back(v);
        Vertex_handle vv = vertices.end();
        return --vv;
    }

    Halfedge_handle edges_push_back( const Halfedge& h, const Halfedge& g) {
        // creates a new pair of opposite border halfedges.
        halfedges.push_back(h);
        Halfedge_handle hh = halfedges.end();
        --hh;
        halfedges.push_back(g);
        Halfedge_handle gg = hh;
        ++gg;
        hh->HBase_base::set_opposite(gg);
        gg->HBase_base::set_opposite(hh);
        return hh;
    }

    Halfedge_handle edges_push_back( const Halfedge& h) {
        CGAL_precondition( h.opposite() != Halfedge_const_handle());
        return edges_push_back( h, *(h.opposite()));
    }

    Face_handle faces_push_back( const Face& f) {
        faces.push_back(f);
        Face_handle ff = faces.end();
        return --ff;
    }


// Removal
//
// The following operations erase an element referenced by a handle.
// Halfedges are always deallocated in pairs of opposite halfedges. Erase
// of single elements is optional. The deletion of all is mandatory.

    void vertices_pop_front()             { vertices.pop_front(); }
    void vertices_pop_back()              { vertices.pop_back(); }
    void vertices_erase( Vertex_handle v) { vertices.erase(v.iterator()); }
    void vertices_erase( Vertex_iterator first, Vertex_iterator last) {
        vertices.erase( first.iterator(), last.iterator());
    }

    void edges_erase( Halfedge_handle h) {
        // deletes the pair of opposite halfedges h and h->opposite().
        Halfedge_handle g = h->opposite();
        halfedges.erase(h.iterator());
        halfedges.erase(g.iterator());
    }
    void edges_pop_front() { edges_erase( halfedges_begin()); }
    void edges_pop_back()  {
        Halfedge_iterator h = halfedges_end();
        edges_erase( --h);
    }
    void edges_erase( Halfedge_iterator first, Halfedge_iterator last) {
        CGAL_precondition( std::distance( first.iterator(),
                                          last.iterator()) % 2 == 0);
        halfedges.erase( first.iterator(), last.iterator());
    }

    void faces_pop_front()           { faces.pop_front(); }
    void faces_pop_back()            { faces.pop_back(); }
    void faces_erase( Face_handle f) { faces.erase(f.iterator()); }
    void faces_erase( Face_iterator first, Face_iterator last) {
        faces.erase( first.iterator(), last.iterator());
    }

    void vertices_clear() { vertices.erase( vertices.begin(), vertices.end());}
    void edges_clear() {
        halfedges.erase( halfedges.begin(), halfedges.end());
        nb_border_halfedges = 0;
        nb_border_edges = 0;
        border_halfedges = Halfedge_handle();
    }
    void faces_clear() { faces.erase( faces.begin(), faces.end()); }
    void clear() {
        vertices_clear();
        edges_clear();
        faces_clear();
    }

// Operations with Border Halfedges

    size_type size_of_border_halfedges() const { return nb_border_halfedges;}
        // number of border halfedges. An edge with no incident face
        // counts as two border halfedges. Precondition: `normalize_border()'
        // has been called and no halfedge insertion or removal and no
        // change in border status of the halfedges have occured since
        // then.

    size_type size_of_border_edges() const { return nb_border_edges;}
        // number of border edges. If `size_of_border_edges() ==
        // size_of_border_halfedges()' all border edges are incident to a
        // face on one side and to a hole on the other side.
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
        CGAL_postcondition( border_halfedges != Halfedge_iterator());
        return border_halfedges;
    }

    Halfedge_const_iterator border_halfedges_begin() const {
        CGAL_postcondition( border_halfedges != Halfedge_iterator());
        return border_halfedges;
    }

    void normalize_border();
        // sorts halfedges such that the non-border edges precedes the
        // border edges. For each border edge that is incident to a face
        // the halfedge iterator will reference the halfedge incident to
        // the face right before the halfedge incident to the hole.
};
#ifdef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
};
#endif


#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
#define CGAL__HalfedgeDS_using_list HalfedgeDS_using_list
#else
#define CGAL__HalfedgeDS_using_list HalfedgeDS_using_list::HDS
#endif

// A class for comparing handles, used in the maps below.
template < class Handle>
struct CGAL__HDS_Cmp_handle_2 {
    bool operator()( Handle a, Handle b) const { return &*a < &*b; }
};

template < class Traits_, class HalfedgeDSItems>
void
CGAL__HalfedgeDS_using_list<Traits_,HalfedgeDSItems>::
pointer_update( const CGAL__HalfedgeDS_using_list<Traits_,
                    HalfedgeDSItems>& hds) {
    // Update own pointers assuming that they lived previously
    // in a halfedge data structure `hds' with lists.
    typedef std::map< Vertex_const_handle,
                      Vertex_handle,
                      CGAL__HDS_Cmp_handle_2< Vertex_const_handle> >    V_map;
    typedef std::map< Halfedge_const_handle,
                      Halfedge_handle,
                      CGAL__HDS_Cmp_handle_2< Halfedge_const_handle > > H_map;
    typedef std::map< Face_const_handle,
                      Face_handle,
                      CGAL__HDS_Cmp_handle_2< Face_const_handle > >     F_map;
    V_map v_map;
    H_map h_map;
    F_map f_map;
    h_map[Halfedge_const_iterator()] = Halfedge_iterator();
    v_map[Vertex_const_iterator()]   = Vertex_iterator();
    f_map[Face_const_iterator()]     = Face_iterator();
    // initialize maps.
    Halfedge_iterator ii = halfedges_begin();
    Halfedge_const_iterator i = hds.halfedges_begin();
    for ( ; i != hds.halfedges_end(); ++i) {
        h_map[i] = ii;
        ++ii;
    }
    h_map[i] = ii;
    if ( check_tag( Supports_halfedge_vertex())) {
        Vertex_iterator vv = vertices_begin();
        for ( Vertex_const_iterator v = hds.vertices_begin();
              v != hds.vertices_end(); ++v) {
            v_map[v] = vv;
            ++vv;
        }
    }
    if ( check_tag( Supports_halfedge_face())) {
        Face_iterator ff = faces_begin();
        for ( Face_const_iterator f = hds.faces_begin();
              f != hds.faces_end(); ++f) {
            f_map[f] = ff;
            ++ff;
        }
    }
    HalfedgeDS_items_decorator<Self> D;
    for ( Halfedge_iterator h = halfedges_begin(); h != halfedges_end(); ++h) {
        h->HBase::set_next( h_map[ h->next()]);
        h->HBase_base::set_opposite( h_map[ h->opposite()]);
        if ( check_tag( Supports_halfedge_prev()))
            D.set_prev( h, h_map[ D.get_prev(h)]);
        if ( check_tag( Supports_halfedge_vertex())) {
            D.set_vertex( h, v_map[ D.get_vertex(h)]);
            D.set_vertex_halfedge( h);
        }
        if ( h->is_border())
            D.set_face( h, Face_handle());
        else if ( check_tag( Supports_halfedge_face())) {
            D.set_face( h, f_map[ D.get_face(h)]);
            D.set_face_halfedge( h);
        }
    }
    border_halfedges = h_map[ border_halfedges];
}

template < class Traits_, class HalfedgeDSItems>
void
CGAL__HalfedgeDS_using_list<Traits_,HalfedgeDSItems>::
normalize_border() {
    CGAL_assertion_code( size_type count = halfedges.size();)
    nb_border_halfedges = 0;
    nb_border_edges = 0;
    Halfedge_list  border;
    Halfedge_I i = halfedges.begin();
    while ( i != halfedges.end()) {
        Halfedge_I j = i;
        ++i;
        ++i;
        Halfedge_I k = j;
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
    if ( i == halfedges.begin()) {
        halfedges.splice( halfedges.end(), border);
        i = halfedges.begin();
    } else {
        --i;
        --i;
        CGAL_assertion( ! i->is_border() && ! i->opposite()->is_border());
        halfedges.splice( halfedges.end(), border);
        ++i;
        ++i;
    }
    CGAL_assertion( i == halfedges.end() || i->opposite()->is_border());
    border_halfedges = i;
}

#undef CGAL__HalfedgeDS_using_list

CGAL_END_NAMESPACE
#endif // CGAL_HALFEDGEDS_USING_LIST_H //
// EOF //
