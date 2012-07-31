// Copyright (c) 1997  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>

#ifndef CGAL_HALFEDGEDS_LIST_H
#define CGAL_HALFEDGEDS_LIST_H 1

#include <CGAL/In_place_list.h>
#include <CGAL/HalfedgeDS_items_decorator.h>
#include <CGAL/memory.h>
#include <CGAL/Unique_hash_map.h>
#include <cstddef>

namespace CGAL {

template < class VertexBase>
class HalfedgeDS_in_place_list_vertex
    : public VertexBase, public In_place_list_base<
                      HalfedgeDS_in_place_list_vertex< VertexBase> > {
public:
    typedef HalfedgeDS_in_place_list_vertex< VertexBase> Self;
    typedef typename VertexBase::Vertex_handle       Vertex_handle;
    typedef typename VertexBase::Vertex_const_handle Vertex_const_handle;
    HalfedgeDS_in_place_list_vertex() {}
    HalfedgeDS_in_place_list_vertex( const VertexBase& v)   // down cast
        : VertexBase(v) {}
    Self& operator=( const Self& v) {
        // This self written assignment avoids that assigning vertices will
        // overwrite the list linking of the target vertex.
        *((VertexBase*)this) = ((const VertexBase&)v);
        return *this;
    }
};

template < class HalfedgeBase>
class HalfedgeDS_in_place_list_halfedge
    : public HalfedgeBase, public In_place_list_base<
                  HalfedgeDS_in_place_list_halfedge< HalfedgeBase> > {
public:
    typedef HalfedgeDS_in_place_list_halfedge< HalfedgeBase> Self;
    typedef typename HalfedgeBase::Halfedge_handle       Halfedge_handle;
    typedef typename HalfedgeBase::Halfedge_const_handle
                                                    Halfedge_const_handle;
    HalfedgeDS_in_place_list_halfedge() {}                   // down cast
    HalfedgeDS_in_place_list_halfedge( const HalfedgeBase& h)
        : HalfedgeBase(h) {}
    Self& operator=( const Self& h) {
        // This self written assignment avoids that assigning halfedges will
        // overwrite the list linking of the target halfedge.
        *((HalfedgeBase*)this) = ((const HalfedgeBase&)h);
        return *this;
    }
};

template < class FaceBase>
class HalfedgeDS_in_place_list_face
    : public FaceBase, public In_place_list_base<
                            HalfedgeDS_in_place_list_face< FaceBase> > {
public:
    typedef HalfedgeDS_in_place_list_face< FaceBase>  Self;
    typedef typename FaceBase::Face_handle       Face_handle;
    typedef typename FaceBase::Face_const_handle Face_const_handle;
    HalfedgeDS_in_place_list_face() {}                   // down cast
    HalfedgeDS_in_place_list_face( const FaceBase& f) : FaceBase(f) {}
    Self& operator=( const Self& f) {
        // This self written assignment avoids that assigning faces will
        // overwrite the list linking of the target face.
        *((FaceBase*)this) = ((const FaceBase&)f);
        // this->FaceBase::operator=(f); // does not compile on SGI
        return *this;
    }
};

template < class Traits_, class HalfedgeDSItems, class Alloc>
class HalfedgeDS_list_types {
public:
    typedef HalfedgeDS_list_types<Traits_, HalfedgeDSItems, Alloc> Self;
    typedef Traits_                                    Traits;
    typedef HalfedgeDSItems                            Items;
    typedef Alloc                                      Allocator;
    typedef Alloc                                      allocator_type;

    typedef typename Items::template Vertex_wrapper<Self,Traits>
                                                       Vertex_wrapper;
    typedef typename Items::template Halfedge_wrapper<Self,Traits> 
                                                       Halfedge_wrapper;
    typedef typename Items::template Face_wrapper<Self,Traits>
                                                       Face_wrapper;

    typedef typename Vertex_wrapper::Vertex            Vertex_base;
    typedef HalfedgeDS_in_place_list_vertex< Vertex_base> Vertex;
    typedef typename Halfedge_wrapper::Halfedge        Halfedge_base;
    typedef HalfedgeDS_in_place_list_halfedge< Halfedge_base> Halfedge;
    typedef typename Face_wrapper::Face                Face_base;
    typedef HalfedgeDS_in_place_list_face< Face_base>  Face;

    typedef typename Allocator::template rebind< Vertex> Vertex_alloc_rebind;
    typedef typename Vertex_alloc_rebind::other        Vertex_allocator;
    typedef typename Allocator::template rebind< Halfedge>
                                                       Halfedge_alloc_rebind;
    typedef typename Halfedge_alloc_rebind::other      Halfedge_allocator;
    typedef typename Allocator::template rebind< Face> Face_alloc_rebind;
    typedef typename Face_alloc_rebind::other          Face_allocator;

    typedef In_place_list<Vertex,false,Vertex_allocator>  Vertex_list;
    typedef typename Vertex_list::iterator             Vertex_handle;
    typedef typename Vertex_list::const_iterator       Vertex_const_handle;
    typedef typename Vertex_list::iterator             Vertex_iterator;
    typedef typename Vertex_list::const_iterator       Vertex_const_iterator;

    typedef In_place_list<Halfedge,false,Halfedge_allocator>  Halfedge_list;
    typedef typename Halfedge_list::iterator           Halfedge_handle;
    typedef typename Halfedge_list::const_iterator     Halfedge_const_handle;
    typedef typename Halfedge_list::iterator           Halfedge_iterator;
    typedef typename Halfedge_list::const_iterator     Halfedge_const_iterator;

    typedef In_place_list<Face,false,Face_allocator>   Face_list;
    typedef typename Face_list::iterator               Face_handle;
    typedef typename Face_list::const_iterator         Face_const_handle;
    typedef typename Face_list::iterator               Face_iterator;
    typedef typename Face_list::const_iterator         Face_const_iterator;

    typedef typename Halfedge_list::size_type          size_type;
    typedef typename Halfedge_list::difference_type    difference_type;
    typedef std::bidirectional_iterator_tag            iterator_category;
    static inline Vertex_handle vertex_handle( Vertex_base* v) {
        Vertex* vv = 0;
        vv = (Vertex*)((char*) v - (std::ptrdiff_t)((Vertex_base*)vv));
        return vv;
    }
    static inline Vertex_const_handle vertex_handle( const Vertex_base* v) {
        const Vertex* vv = 0;
        vv = (const Vertex*)((const char*) v -
                 (std::ptrdiff_t)((const Vertex_base*)vv));
        return vv;
    }

    static inline Halfedge_handle halfedge_handle( Halfedge_base* h) {
        Halfedge* hh = 0;
        hh = (Halfedge*)((char*) h - (std::ptrdiff_t)((Halfedge_base*)hh));
        return hh;
    }
    static inline
    Halfedge_const_handle halfedge_handle( const Halfedge_base* h) {
        const Halfedge* hh = 0;
        hh = (const Halfedge*)((const char*) h -
                 (std::ptrdiff_t)((const Halfedge_base*)hh));
        return hh;
    }

    static inline Face_handle face_handle( Face_base* f) {
        Face* ff = 0;
        ff = (Face*)((char*) f - (std::ptrdiff_t)((Face_base*)ff));
        return ff;
    }
    static inline Face_const_handle face_handle( const Face_base* f) {
        const Face* ff = 0;
        ff = (const Face*)((const char*)f -
                 (std::ptrdiff_t)((const Face_base*)ff));
        return ff;
    }
};

template < class Traits_, class HalfedgeDSItems, 
           class Alloc = CGAL_ALLOCATOR(int)>
class HalfedgeDS_list
    : public HalfedgeDS_list_types<Traits_, HalfedgeDSItems, Alloc> {
public:
    typedef HalfedgeDS_list<Traits_, HalfedgeDSItems, Alloc> Self;
public:
    typedef HalfedgeDS_list_types<Traits_, HalfedgeDSItems, Alloc> Types;
    typedef typename Types::Traits                     Traits;
    typedef typename Types::Items                      Items;
    typedef typename Types::Allocator                  Allocator;
    typedef typename Types::allocator_type             allocator_type;

    typedef typename Types::Vertex                     Vertex;
    typedef typename Types::Halfedge                   Halfedge;
    typedef typename Types::Face                       Face;

    typedef typename Types::Vertex_allocator           Vertex_allocator;
    typedef typename Types::Vertex_list                Vertex_list;
    typedef typename Types::Vertex_handle              Vertex_handle;
    typedef typename Types::Vertex_const_handle        Vertex_const_handle;
    typedef typename Types::Vertex_iterator            Vertex_iterator;
    typedef typename Types::Vertex_const_iterator      Vertex_const_iterator;

    typedef typename Types::Halfedge_allocator         Halfedge_allocator;
    typedef typename Types::Halfedge_list              Halfedge_list;
    typedef typename Types::Halfedge_handle            Halfedge_handle;
    typedef typename Types::Halfedge_const_handle      Halfedge_const_handle;
    typedef typename Types::Halfedge_iterator          Halfedge_iterator;
    typedef typename Types::Halfedge_const_iterator    Halfedge_const_iterator;

    typedef typename Types::Face_allocator             Face_allocator;
    typedef typename Types::Face_list                  Face_list;
    typedef typename Types::Face_handle                Face_handle;
    typedef typename Types::Face_const_handle          Face_const_handle;
    typedef typename Types::Face_iterator              Face_iterator;
    typedef typename Types::Face_const_iterator        Face_const_iterator;

    typedef typename Types::size_type                  size_type;
    typedef typename Types::difference_type            difference_type;
    typedef typename Types::iterator_category          iterator_category;

    typedef Tag_true                                   Supports_removal;

    typedef typename Vertex::Supports_vertex_halfedge Supports_vertex_halfedge;
    typedef typename Halfedge::Supports_halfedge_prev Supports_halfedge_prev;
    typedef typename Halfedge::Supports_halfedge_vertex
                                                      Supports_halfedge_vertex;
    typedef typename Halfedge::Supports_halfedge_face
                                                      Supports_halfedge_face;
    typedef typename Face::Supports_face_halfedge     Supports_face_halfedge;

    // Halfedges are allocated in pairs. Here is the type for that.
    typedef std::pair<Halfedge,Halfedge>              Halfedge_pair;

    typedef typename Allocator::template rebind< Halfedge_pair>
                                                       Edge_alloc_rebind;
    typedef typename Edge_alloc_rebind::other          Edge_allocator;

protected:
    // Changed from static to local variable
    Vertex_allocator vertex_allocator;
    Edge_allocator   edge_allocator;  // allocates pairs of halfedges
    Face_allocator   face_allocator;
    
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
        h2->HBase_base::set_opposite(g2);
        g2->HBase_base::set_opposite(h2);
        return h2;
    }
    void put_edge_node( Halfedge* h) {
        // deletes the pair of opposite halfedges h and h->opposite().
        Halfedge_handle g = h->opposite();
        Halfedge_pair* hpair = (Halfedge_pair*)(&*h);
        if ( &*h > &*g)
            hpair = (Halfedge_pair*)(&*g);
        CGAL_assertion( &(hpair->first) == (Halfedge*)hpair);
        edge_allocator.destroy( hpair);
        edge_allocator.deallocate( hpair, 1);
    }

    Face* get_face_node( const Face& t) {
        Face* p = face_allocator.allocate(1);
        face_allocator.construct(p, t);
        return p;
    }
    void put_face_node( Face* p) {
        face_allocator.destroy( p);
        face_allocator.deallocate( p, 1);
    }

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
    void pointer_update( const Self& hds) {
        // Update own pointers assuming that they lived previously
        // in a halfedge data structure `hds' with lists.
        // Update own pointers assuming that they lived previously
        // in a halfedge data structure `hds' with lists.
        typedef Unique_hash_map< Vertex_const_iterator, Vertex_iterator> V_map;
        typedef Unique_hash_map< Halfedge_const_iterator, Halfedge_iterator>
                                                                         H_map;
        typedef Unique_hash_map< Face_const_iterator, Face_iterator>     F_map;
        // initialize maps.
        H_map h_map( hds.halfedges_begin(), hds.halfedges_end(),
                     halfedges_begin(), Halfedge_iterator(), 
                     3 * hds.size_of_halfedges() / 2);
        Vertex_iterator vii;
        V_map v_map( vii, 3 * hds.size_of_vertices() / 2);
        Face_iterator fii;
        F_map f_map( fii, 3 * hds.size_of_faces() / 2);
        // some special values
        h_map[Halfedge_const_iterator()] = Halfedge_iterator();
        h_map[hds.halfedges_end()]       = halfedges_end();
        v_map[Vertex_const_iterator()]   = Vertex_iterator();
        v_map[hds.vertices_end()]        = vertices_end();
        f_map[Face_const_iterator()]     = Face_iterator();
        f_map[hds.faces_end()]           = faces_end();
        // vertices and faces are optional
        if ( check_tag( Supports_halfedge_vertex())) {
            v_map.insert( hds.vertices_begin(),
                          hds.vertices_end(),
                          vertices_begin());
        }
        if ( check_tag( Supports_halfedge_face())) {
            f_map.insert( hds.faces_begin(), hds.faces_end(), faces_begin());
        }
        HalfedgeDS_items_decorator<Self> D;
        for (Halfedge_iterator h = halfedges_begin(); h!=halfedges_end(); ++h){
            h->HBase::set_next( h_map[ h->next()]);
            // Superfluous and false: opposite pointer get set upon creation
            // h->HBase_base::set_opposite( h_map[ h->opposite()]);
            if ( check_tag( Supports_halfedge_prev()))
                D.set_prev( h, h_map[ D.get_prev(h)]);
            if ( check_tag( Supports_halfedge_vertex()))
                D.set_vertex( h, v_map[ D.get_vertex(h)]);
            if ( check_tag( Supports_halfedge_face()))
                D.set_face( h, f_map[ D.get_face(h)]);
        }
        border_halfedges = h_map[ border_halfedges];
        if (check_tag( Supports_vertex_halfedge())) {
            for (Vertex_iterator v = vertices_begin(); v != vertices_end();++v)
                D.set_vertex_halfedge(v, h_map[ D.get_vertex_halfedge(v)]);
        }
        if (check_tag( Supports_face_halfedge())) {
            for ( Face_iterator f = faces_begin(); f != faces_end(); ++f)
                D.set_face_halfedge(f, h_map[ D.get_face_halfedge(f)]);
        }
        //h_map.statistics();
        //v_map.statistics();
        //f_map.statistics();
    }

public:
    HalfedgeDS_list()
        : nb_border_halfedges(0), nb_border_edges(0) {}
        // the empty polyhedron `P'.

    HalfedgeDS_list( size_type, size_type, size_type)
        : nb_border_halfedges(0), nb_border_edges(0) {}
        // Parameter order is v,h,f.
        // a polyhedron `P' with storage reserved for v vertices, h
        // halfedges, and f faces. The reservation sizes are a hint for
        // optimizing storage allocation. They are not used here.

    ~HalfedgeDS_list() { clear(); }

    HalfedgeDS_list( const Self& hds)
    :  vertices( hds.vertices),
       //halfedges( hds.halfedges),
       faces( hds.faces),
       nb_border_halfedges( hds.nb_border_halfedges),
       nb_border_edges( hds.nb_border_edges),
       border_halfedges( hds.border_halfedges)
    {
        // goal is halfedges = hds.halfedges, but we have pairs here
        Halfedge_const_iterator i = hds.halfedges_begin();
        for ( ; i != hds.halfedges_end(); ++ ++ i) {
            edges_push_back( *i);
        }
        pointer_update( hds);
    }

    Self& operator=( const Self& hds)  {
        if ( this != &hds) {
            clear();
            vertices            = hds.vertices;
            // goal is halfedges = hds.halfedges, but we have pairs here
            halfedges = Halfedge_list();
            Halfedge_const_iterator i = hds.halfedges_begin();
            for ( ; i != hds.halfedges_end(); ++ ++ i) {
                edges_push_back( *i);
            }
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

    allocator_type  get_allocator() const { return allocator_type(); }

    size_type size_of_vertices() const  { return vertices.size();}
    size_type size_of_halfedges() const { return halfedges.size();}
        // number of all halfedges (including border halfedges).
    size_type size_of_faces() const     { return faces.size();}

    size_type capacity_of_vertices() const    { return vertices.max_size();}
    size_type capacity_of_halfedges() const   { return halfedges.max_size();}
    size_type capacity_of_faces() const       { return faces.max_size();}

    std::size_t bytes() const {
        return sizeof(Self)
               + vertices.size()  * sizeof( Vertex)
               + halfedges.size() * sizeof( Halfedge)
               + faces.size()     * sizeof( Face);
    }
    std::size_t bytes_reserved() const { return bytes();}

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
        vertices.push_back( * get_vertex_node(v));
        Vertex_handle vh = vertices.end();
        return --vh;
    }

    Halfedge_handle edges_push_back( const Halfedge& h, const Halfedge& g) {
        // creates a new pair of opposite border halfedges.
        Halfedge* ptr = get_edge_node( h, g);
        halfedges.push_back( *ptr);
        Halfedge_handle hh = halfedges.end();
        --hh;
        halfedges.push_back( *(ptr->opposite()));
        return hh;
    }

    Halfedge_handle edges_push_back( const Halfedge& h) {
        CGAL_precondition( h.opposite() != Halfedge_const_handle());
        return edges_push_back( h, *(h.opposite()));
    }

    Face_handle faces_push_back( const Face& f) {
        faces.push_back( * get_face_node(f));
        Face_handle fh = faces.end();
        return --fh;
    }


// Removal
//
// The following operations erase an element referenced by a handle.
// Halfedges are always deallocated in pairs of opposite halfedges. Erase
// of single elements is optional. The deletion of all is mandatory.

    void vertices_pop_front() {
        Vertex* v = &(vertices.front());
        vertices.pop_front();
        put_vertex_node( v);
    }
    void vertices_pop_back() {
        Vertex* v = &(vertices.back());
        vertices.pop_back();
        put_vertex_node( v);
    }
    void vertices_erase( Vertex_handle v) {
        Vertex* ptr = &*v;
        vertices.erase(v);
        put_vertex_node( ptr);
    }
    void vertices_erase( Vertex_iterator first, Vertex_iterator last) {
        while (first != last)
            vertices_erase(first++);
    }

    void edges_erase( Halfedge_handle h) {
        // deletes the pair of opposite halfedges h and h->opposite().
        Halfedge_handle g = h->opposite();
        halfedges.erase(h);
        halfedges.erase(g);
        put_edge_node(&*h);
    }
    void edges_pop_front() { edges_erase( halfedges.begin()); }
    void edges_pop_back()  {
        Halfedge_iterator h = halfedges.end();
        edges_erase( --h);
    }
    void edges_erase( Halfedge_iterator first, Halfedge_iterator last) {
        while (first != last) {
            Halfedge_iterator nxt = first;
            ++nxt;
            CGAL_assertion( nxt != last);
            ++nxt;
            edges_erase(first);
            first = nxt;
        }
    }

    void faces_pop_front() {
        Face* f = &(faces.front());
        faces.pop_front();
        put_face_node( f);
    }
    void faces_pop_back() {
        Face* f = &(faces.back());
        faces.pop_back();
        put_face_node( f);
    }
    void faces_erase( Face_handle f) {
        Face* ptr = &*f;
        faces.erase(f);
        put_face_node( ptr);
    }
    void faces_erase( Face_iterator first, Face_iterator last) {
        while (first != last)
            faces_erase(first++);
    }

    void vertices_clear() { vertices.destroy(); }
    void edges_clear() {
        edges_erase( halfedges.begin(), halfedges.end());
        nb_border_halfedges = 0;
        nb_border_edges = 0;
        border_halfedges = Halfedge_handle();
    }
    void faces_clear() { faces.destroy(); }
    void clear() {
        vertices_clear();
        edges_clear();
        faces_clear();
    }

    void vertices_splice( Vertex_iterator target, Self &source,
                          Vertex_iterator begin, Vertex_iterator end) {
        vertices.splice( target, source.vertices, begin, end);
    }

    void halfedges_splice( Halfedge_iterator target, Self &source,
                           Halfedge_iterator begin, Halfedge_iterator end) {
        halfedges.splice( target, source.halfedges, begin, end);
    }

    void faces_splice( Face_iterator target, Self &source,
                       Face_iterator begin, Face_iterator end) {
        faces.splice( target, source.faces, begin, end);
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
        return border_halfedges;
    }

    Halfedge_const_iterator border_halfedges_begin() const {
        return border_halfedges;
    }

    void normalize_border() {
        // sorts halfedges such that the non-border edges precedes the
        // border edges. For each border edge that is incident to a face
        // the halfedge iterator will reference the halfedge incident to
        // the face right before the halfedge incident to the hole.
        CGAL_assertion_code( size_type count = halfedges.size();)
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
            CGAL_assertion( ! i->is_border() && ! i->opposite()->is_border());
            halfedges.splice( halfedges.end(), border);
            ++i;
            ++i;
        }
        CGAL_assertion( i == halfedges_end() || i->opposite()->is_border());
        border_halfedges = i;
    }
};


//  #ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
//  #define CGAL__HDS_IP_List HalfedgeDS_list
//  #else
//  #define CGAL__HDS_IP_List HalfedgeDS_list::HDS
//  #endif

// init static member allocator objects (no longer static)
//template < class Traits_, class HalfedgeDSItems, class Alloc>
//typename CGAL__HDS_IP_List<Traits_, HalfedgeDSItems, Alloc>::Vertex_allocator
//CGAL__HDS_IP_List<Traits_, HalfedgeDSItems, Alloc>::vertex_allocator;
//
//template < class Traits_, class HalfedgeDSItems, class Alloc>
//typename CGAL__HDS_IP_List<Traits_, HalfedgeDSItems, Alloc>::Edge_allocator
//CGAL__HDS_IP_List<Traits_, HalfedgeDSItems, Alloc>::edge_allocator;
//
//template < class Traits_, class HalfedgeDSItems, class Alloc>
//typename CGAL__HDS_IP_List<Traits_, HalfedgeDSItems, Alloc>::Face_allocator
//CGAL__HDS_IP_List<Traits_, HalfedgeDSItems, Alloc>::face_allocator;


//  #undef CGAL__HDS_IP_List

} //namespace CGAL
#endif // CGAL_HALFEDGEDS_LIST_H //
// EOF //
