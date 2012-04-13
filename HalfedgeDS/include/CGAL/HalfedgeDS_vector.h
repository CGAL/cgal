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

#ifndef CGAL_HALFEDGEDS_VECTOR_H
#define CGAL_HALFEDGEDS_VECTOR_H 1

// Define this if HalfedgeDS_vector is based on CGAL internal vector.
#define CGAL__HALFEDGEDS_USE_INTERNAL_VECTOR 1

#include <CGAL/basic.h>
#include <CGAL/memory.h>
#include <CGAL/HalfedgeDS_items_decorator.h>
#include <algorithm>
#include <vector>
#include <map>
#include <cstddef>

#ifdef CGAL__HALFEDGEDS_USE_INTERNAL_VECTOR
#include <CGAL/vector.h>
#else
#include <CGAL/HalfedgeDS_iterator_adaptor.h>
#endif

namespace CGAL {


template < class Traits_, class HalfedgeDSItems, class Alloc>
class HalfedgeDS_vector_types {
public:
    typedef HalfedgeDS_vector_types<Traits_,HalfedgeDSItems,Alloc> Self;
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

    typedef typename Vertex_wrapper::Vertex            Vertex;
    typedef typename Halfedge_wrapper::Halfedge        Halfedge;
    typedef typename Face_wrapper::Face                Face;

    typedef typename Allocator::template rebind< Vertex> Vertex_alloc_rebind;
    typedef typename Vertex_alloc_rebind::other        Vertex_allocator;
    typedef typename Allocator::template rebind< Halfedge>
                                                       Halfedge_alloc_rebind;
    typedef typename Halfedge_alloc_rebind::other      Halfedge_allocator;
    typedef typename Allocator::template rebind< Face> Face_alloc_rebind;
    typedef typename Face_alloc_rebind::other          Face_allocator;

#ifdef CGAL__HALFEDGEDS_USE_INTERNAL_VECTOR
    typedef internal::vector<Vertex, Vertex_allocator>    Vertex_vector;
    typedef typename Vertex_vector::iterator           Vertex_I;
    typedef typename Vertex_vector::const_iterator     Vertex_CI;
    typedef typename Vertex_vector::iterator           Vertex_iterator;
    typedef typename Vertex_vector::const_iterator     Vertex_const_iterator;

    typedef internal::vector<Halfedge, Halfedge_allocator>  Halfedge_vector;
    typedef typename Halfedge_vector::iterator         Halfedge_I;
    typedef typename Halfedge_vector::const_iterator   Halfedge_CI;
    typedef typename Halfedge_vector::iterator         Halfedge_iterator;
    typedef typename Halfedge_vector::const_iterator   Halfedge_const_iterator;

    typedef internal::vector<Face, Face_allocator>        Face_vector;
    typedef typename Face_vector::iterator             Face_I;
    typedef typename Face_vector::const_iterator       Face_CI;
    typedef typename Face_vector::iterator             Face_iterator;
    typedef typename Face_vector::const_iterator       Face_const_iterator;
#else // CGAL__HALFEDGEDS_USE_INTERNAL_VECTOR //
    typedef std::vector<Vertex, Vertex_allocator>      Vertex_vector;
    typedef typename Vertex_vector::iterator           Vertex_I;
    typedef typename Vertex_vector::const_iterator     Vertex_CI;
    typedef HalfedgeDS_iterator_adaptor<Vertex_I>      Vertex_iterator;
    typedef HalfedgeDS_iterator_adaptor<Vertex_CI>     Vertex_const_iterator;

    typedef std::vector<Halfedge, Halfedge_allocator>  Halfedge_vector;
    typedef typename Halfedge_vector::iterator         Halfedge_I;
    typedef typename Halfedge_vector::const_iterator   Halfedge_CI;
    typedef HalfedgeDS_iterator_adaptor<Halfedge_I>    Halfedge_iterator;
    typedef HalfedgeDS_iterator_adaptor<Halfedge_CI>   Halfedge_const_iterator;

    typedef std::vector<Face, Face_allocator>          Face_vector;
    typedef typename Face_vector::iterator             Face_I;
    typedef typename Face_vector::const_iterator       Face_CI;
    typedef HalfedgeDS_iterator_adaptor<Face_I>        Face_iterator;
    typedef HalfedgeDS_iterator_adaptor<Face_CI>       Face_const_iterator;
#endif // CGAL__HALFEDGEDS_USE_INTERNAL_VECTOR //

    typedef Vertex_iterator                            Vertex_handle;
    typedef Vertex_const_iterator                      Vertex_const_handle;
    typedef Halfedge_iterator                          Halfedge_handle;
    typedef Halfedge_const_iterator                    Halfedge_const_handle;
    typedef Face_iterator                              Face_handle;
    typedef Face_const_iterator                        Face_const_handle;

    typedef typename Halfedge_vector::size_type        size_type;
    typedef typename Halfedge_vector::difference_type  difference_type;
    typedef std::random_access_iterator_tag            iterator_category;

    static inline Vertex_handle vertex_handle( Vertex* v) {
        return Vertex_handle( Vertex_I(v));
    }
    static inline Vertex_const_handle vertex_handle( const Vertex* v) {
        return Vertex_const_handle( Vertex_CI(v));
    }
    static inline Halfedge_handle halfedge_handle( Halfedge* h) {
        return Halfedge_handle( Halfedge_I(h));
    }
    static inline Halfedge_const_handle halfedge_handle( const Halfedge* h) {
        return Halfedge_const_handle( Halfedge_CI(h));
    }
    static inline Face_handle face_handle( Face* f) {
        return Face_handle( Face_I(f));
    }
    static inline Face_const_handle face_handle( const Face* f) {
        return Face_const_handle( Face_CI(f));
    }
};


template < class Traits_, class HalfedgeDSItems, 
           class Alloc = CGAL_ALLOCATOR(int)>
class HalfedgeDS_vector
    : public HalfedgeDS_vector_types<Traits_, HalfedgeDSItems, Alloc> {
public:
    typedef HalfedgeDS_vector<Traits_,HalfedgeDSItems,Alloc> Self;
    typedef HalfedgeDS_vector_types<Traits_, HalfedgeDSItems, Alloc> Types;
    typedef typename Types::Traits                     Traits;
    typedef typename Types::Items                      Items;
    typedef typename Types::Allocator                  Allocator;
    typedef typename Types::allocator_type             allocator_type;

    typedef typename Types::Vertex                     Vertex;
    typedef typename Types::Halfedge                   Halfedge;
    typedef typename Types::Face                       Face;

    typedef typename Types::Vertex_vector              Vertex_vector;
    typedef typename Types::Vertex_handle              Vertex_handle;
    typedef typename Types::Vertex_const_handle        Vertex_const_handle;
    typedef typename Types::Vertex_iterator            Vertex_iterator;
    typedef typename Types::Vertex_const_iterator      Vertex_const_iterator;
    typedef typename Types::Vertex_I                   Vertex_I;
    typedef typename Types::Vertex_CI                  Vertex_CI;

    typedef typename Types::Halfedge_vector            Halfedge_vector;
    typedef typename Types::Halfedge_handle            Halfedge_handle;
    typedef typename Types::Halfedge_const_handle      Halfedge_const_handle;
    typedef typename Types::Halfedge_iterator          Halfedge_iterator;
    typedef typename Types::Halfedge_const_iterator    Halfedge_const_iterator;
    typedef typename Types::Halfedge_I                 Halfedge_I;
    typedef typename Types::Halfedge_CI                Halfedge_CI;

    typedef typename Types::Face_vector                Face_vector;
    typedef typename Types::Face_handle                Face_handle;
    typedef typename Types::Face_const_handle          Face_const_handle;
    typedef typename Types::Face_iterator              Face_iterator;
    typedef typename Types::Face_const_iterator        Face_const_iterator;
    typedef typename Types::Face_I                     Face_I;
    typedef typename Types::Face_CI                    Face_CI;

    typedef typename Types::size_type                  size_type;
    typedef typename Types::difference_type            difference_type;
    typedef typename Types::iterator_category          iterator_category;

    typedef Tag_false                                  Supports_removal;

    typedef typename Vertex::Supports_vertex_halfedge Supports_vertex_halfedge;
    typedef typename Halfedge::Supports_halfedge_prev Supports_halfedge_prev;
    typedef typename Halfedge::Supports_halfedge_vertex
                                                      Supports_halfedge_vertex;
    typedef typename Halfedge::Supports_halfedge_face
                                                      Supports_halfedge_face;
    typedef typename Face::Supports_face_halfedge     Supports_face_halfedge;

protected:
    typedef typename Vertex::Base                      VBase;
    typedef typename Halfedge::Base                    HBase;
    typedef typename Halfedge::Base_base               HBase_base;
    typedef typename Face::Base                        FBase;

    Vertex_vector      vertices;
    Halfedge_vector    halfedges;
    Face_vector        faces;

    size_type          nb_border_halfedges;
    size_type          nb_border_edges;
    Halfedge_iterator  border_halfedges;

#ifdef CGAL__HALFEDGEDS_USE_INTERNAL_VECTOR
    Vertex_I    get_v_iter( const Vertex_I&  i)   const { return i; }
    Vertex_CI   get_v_iter( const Vertex_CI& i)   const { return i; }
    Halfedge_I  get_h_iter( const Halfedge_I&  i) const { return i; }
    Halfedge_CI get_h_iter( const Halfedge_CI& i) const { return i; }
    Face_I      get_f_iter( const Face_I&  i)     const { return i; }
    Face_CI     get_f_iter( const Face_CI& i)     const { return i; }
#else // CGAL__HALFEDGEDS_USE_INTERNAL_VECTOR //
    Vertex_I    get_v_iter( const Vertex_iterator&  i) const {
        return i.iterator();
    }
    Vertex_CI   get_v_iter( const Vertex_const_iterator& i)   const {
        return i.iterator();
    }
    Halfedge_I  get_h_iter( const Halfedge_iterator&  i) const {
        return i.iterator();
    }
    Halfedge_CI get_h_iter( const Halfedge_const_iterator& i) const {
        return i.iterator();
    }
    Face_I      get_f_iter( const Face_iterator&  i)     const {
        return i.iterator();
    }
    Face_CI     get_f_iter( const Face_const_iterator& i)     const {
        return i.iterator();
    }
#endif // CGAL__HALFEDGEDS_USE_INTERNAL_VECTOR //

// CREATION

private:
#define CGAL__V_UPDATE(v) (((v) == Vertex_handle()) ? (v) : \
                           (v_new + ( Vertex_CI   (get_v_iter(v)) - v_old)))
#define CGAL__H_UPDATE(h) (((h) == Halfedge_handle()) ? (h) : \
                           (h_new + ( Halfedge_CI (get_h_iter(h)) - h_old)))
#define CGAL__F_UPDATE(f) (((f) == Face_handle()) ? (f) : \
                           (f_new + ( Face_CI     (get_f_iter(f)) - f_old)))

    void pointer_update( Vertex_CI v_old, Halfedge_CI h_old, Face_CI f_old) {
        // Update own pointers assuming that they lived previously
        // in a halfedge data structure with vector starting addresses
        // as given as parameters v_old, h_old, f_old.
        HalfedgeDS_items_decorator<Self> D;
        Vertex_iterator   v_new = vertices.begin();
        Halfedge_iterator h_new = halfedges.begin();
        Face_iterator     f_new = faces.begin();
        for (Halfedge_iterator h = halfedges_begin();h != halfedges_end();++h){
            h->HBase::set_next( CGAL__H_UPDATE( h->next()));
            h->HBase_base::set_opposite( CGAL__H_UPDATE( h->opposite()));
            D.set_prev(   h, CGAL__H_UPDATE( D.get_prev(h)));
            D.set_vertex( h, CGAL__V_UPDATE( D.get_vertex(h)));
            D.set_face(   h, CGAL__F_UPDATE( D.get_face(h)));
        }
        border_halfedges = CGAL__H_UPDATE( border_halfedges);
        if (check_tag( Supports_vertex_halfedge())) {
            for (Vertex_iterator v = vertices_begin();v != vertices_end();++v){
                D.set_vertex_halfedge(v, CGAL__H_UPDATE(
                    D.get_vertex_halfedge(v)));
            }
        }
        if (check_tag( Supports_face_halfedge())) {
            for ( Face_iterator f = faces_begin(); f != faces_end(); ++f) {
                D.set_face_halfedge(f,CGAL__H_UPDATE( D.get_face_halfedge(f)));
            }
        }
    }
#undef CGAL__V_UPDATE
#undef CGAL__H_UPDATE
#undef CGAL__F_UPDATE

public:

    HalfedgeDS_vector()
        : nb_border_halfedges(0), nb_border_edges(0) {}
        // empty halfedge data structure.

    HalfedgeDS_vector( size_type v, size_type h, size_type f)
        : nb_border_halfedges(0), nb_border_edges(0) {
        // halfedge data structure with storage reserved for v vertices, h
        // halfedges, and f faces. The reservation sizes are a hint for
        // optimizing storage allocation. They are not used here.
        vertices.reserve(v);
        halfedges.reserve(h);
        faces.reserve(f);
    }

    HalfedgeDS_vector( const Self& hds)
    :  vertices( hds.vertices),
       halfedges( hds.halfedges),
       faces( hds.faces),
       nb_border_halfedges( hds.nb_border_halfedges),
       nb_border_edges( hds.nb_border_edges),
       border_halfedges( hds.border_halfedges)
    {
        pointer_update( hds.vertices.begin(),
                        hds.halfedges.begin(),
                        hds.faces.begin());
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
            pointer_update( hds.vertices.begin(),
                            hds.halfedges.begin(),
                            hds.faces.begin());
        }
        return *this;
    }

    void reserve( size_type v, size_type h, size_type f) {
        // reserve storage for v vertices, h halfedges, and f faces. The
        // reservation sizes are a hint for optimizing storage allocation.
        // If the `capacity' is already greater than the requested size
        // nothing happens. If the `capacity' changes all iterators and
        // circulators invalidates. Function is void here.
        if ( (check_tag( Supports_halfedge_vertex())
                && v > capacity_of_vertices())
             || h > capacity_of_halfedges()
             || (check_tag( Supports_halfedge_face())
                && f > capacity_of_faces())) {
            Vertex_CI   v_old = vertices.begin();
            Halfedge_CI h_old = halfedges.begin();
            Face_CI     f_old = faces.begin();
            if ( check_tag( Supports_halfedge_vertex()))
                vertices.reserve(v);
            halfedges.reserve(h);
            if ( check_tag( Supports_halfedge_face()))
                faces.reserve(f);
            pointer_update( v_old, h_old, f_old);
        }
    }

// Access Member Functions

    allocator_type  get_allocator() const { return allocator_type(); }

    size_type size_of_vertices() const  { return vertices.size();}
    size_type size_of_halfedges() const { return halfedges.size();}
        // number of all halfedges (including border halfedges).
    size_type size_of_faces() const     { return faces.size();}

    size_type capacity_of_vertices() const    { return vertices.capacity();}
    size_type capacity_of_halfedges() const   { return halfedges.capacity();}
    size_type capacity_of_faces() const       { return faces.capacity();}

    std::size_t bytes() const {
        return sizeof(Self)
               + vertices.size()  * sizeof( Vertex)
               + halfedges.size() * sizeof( Halfedge)
               + faces.size()     * sizeof( Face);
    }
    std::size_t bytes_reserved() const {
        return sizeof(Self)
               + vertices.capacity()  * sizeof( Vertex)
               + halfedges.capacity() * sizeof( Halfedge)
               + faces.capacity()     * sizeof( Face);
    }

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
        CGAL_precondition( 1+size_of_vertices() <= capacity_of_vertices());
        vertices.push_back(v);
        Vertex_handle vv = vertices.end();
        return --vv;
    }

    Halfedge_handle edges_push_back( const Halfedge& h, const Halfedge& g) {
        // creates a new pair of opposite border halfedges.
        CGAL_precondition( 1 + size_of_halfedges() < capacity_of_halfedges());
        halfedges.push_back(h);
        Halfedge_handle hh = halfedges.end();
        --hh;
        halfedges.push_back(g);
        Halfedge_handle gg = halfedges.end();
        --gg;
        CGAL_assertion( hh + 1 == gg);
        CGAL_assertion( (char*)(&*gg) - (char*)(&*hh) == sizeof( Halfedge));
        hh->HBase_base::set_opposite(gg);
        gg->HBase_base::set_opposite(hh);
        return hh;
    }

    Halfedge_handle edges_push_back( const Halfedge& h) {
        CGAL_precondition( h.opposite() != Halfedge_const_handle());
        return edges_push_back( h, *(h.opposite()));
    }

    Face_handle faces_push_back( const Face& f) {
        CGAL_precondition( 1+size_of_faces() <= capacity_of_faces());
        faces.push_back(f);
        Face_handle ff = faces.end();
        return --ff;
    }


// Removal
//
// The following operations erase an element referenced by a handle.
// Halfedges are always deallocated in pairs of opposite halfedges. Erase
// of single elements is optional. The deletion of all is mandatory.

    void vertices_pop_back() { vertices.pop_back(); }
    void edges_pop_back()    {
        CGAL_precondition(( halfedges_end()-1)->opposite() ==
                          ( halfedges_end()-2));
        halfedges.pop_back();
        halfedges.pop_back();
    }
    void faces_pop_back()    { faces.pop_back(); }

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

// Special Operations on Polyhedral Surfaces
protected:
    // Update operation used in pointer_update(...).
    void update_opposite( Halfedge_I h) {
        Halfedge_I g = h + 1;
        h->HBase_base::set_opposite(g);
        g->HBase_base::set_opposite(h);
    }

// Operations with Border Halfedges
public:
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
        nb_border_halfedges = 0;
        nb_border_edges = 0;
        border_halfedges = halfedges_end();
        // Lets run one partition step over the array of halfedges.
        // First find a pivot -- that means a border edge.
        Halfedge_I ll = halfedges.begin();
        while ( ll != halfedges.end() && (! ll->is_border()) &&
                (! ll->opposite()->is_border() ))
            ++ ++ll;
        if ( ll == halfedges.end()) // Done. No border edges found.
            return;

        // An array of pointers to update the changed halfedge pointers.
        typedef typename Allocator::template rebind< Halfedge_I>
                                                      HI_alloc_rebind;
        typedef typename HI_alloc_rebind::other       HI_allocator;

        typedef std::vector<Halfedge_I, HI_allocator> HVector;
        typedef typename HVector::iterator Hiterator;
        HVector hvector;
        // Initialize it.
        hvector.reserve( halfedges.size());
        for ( Halfedge_I i = halfedges.begin(); i != halfedges.end(); ++i) {
            hvector.push_back(i);
        }
        Hiterator llhv = hvector.begin() + (ll-halfedges.begin());
        // Start with the partitioning.
        Halfedge_I rr = halfedges.end();
        -- --rr;
        Hiterator rrhv = hvector.end();
        -- --rrhv;
        // The comments proove the invariant of the partitioning step.
        // Note that + 1 or - 1 denotes plus one edge or minus one edge,
        // so they mean actually + 2 and - 2.
                              // Pivot is in *ll
                              // Elements in [rr+1..end) >= pivot (border)
                              // Elements in [begin..ll) <  pivot (non border)
        while (ll < rr) {
                              // Pivot is in *ll, ll <= rr.
            while ( rr > ll && (rr->is_border() 
                              || rr->opposite()->is_border())) {
                if ( ! rr->opposite()->is_border()) {
                    CGAL_assertion( rr + 1 == get_h_iter(rr->opposite()));
                    std::swap( *rr, *(rr+1));
                    update_opposite( rr);
                    std::swap( *rrhv, *(rrhv+1));
                }
                -- --rr;
                -- --rrhv;
            }
                              // Elements in [rr+1..end) >= pivot (border)
                              // *rr <= pivot, ll <= rr.
            CGAL_assertion( rr + 1 == get_h_iter( rr->opposite()));
            CGAL_assertion( ll + 1 == get_h_iter( ll->opposite()));
            std::swap( *(ll+1), *(rr+1));
            std::swap( *ll, *rr);
            update_opposite( ll);
            update_opposite( rr);
            std::swap( *(llhv+1), *(rrhv+1));
            std::swap( *llhv, *rrhv);
                              // Elements in [begin..ll) < pivot
                              // Pivot is in *rr, ll <= rr.
            while ( !ll->is_border() && ! ll->opposite()->is_border()) {
                ++ ++ll;
                ++ ++llhv;
            }
                              // Elements in [begin..ll) < pivot
                              // *ll >= pivot
                              // ll <= rr (since *rr is pivot.)
            CGAL_assertion( ll <= rr);
            CGAL_assertion( llhv <= rrhv);
            CGAL_assertion( rr + 1 == get_h_iter( rr->opposite()));
            CGAL_assertion( ll + 1 == get_h_iter( ll->opposite()));
            std::swap( *(ll+1), *(rr+1));
            std::swap( *ll, *rr);
            update_opposite( ll);
            update_opposite( rr);
            std::swap( *(llhv+1), *(rrhv+1));
            std::swap( *llhv, *rrhv);
            if ( ! rr->opposite()->is_border()) {
                CGAL_assertion( rr + 1 == get_h_iter( rr->opposite()));
                std::swap( *rr, *(rr+1));
                update_opposite( rr);
                std::swap( *rrhv, *(rrhv+1));
            }

            // This guard is needed here because, rr==ll==begin, might be true
            // at this point, causing the decrement to result in undefined
            // behaviour.
            // [Fernando Cacciola]
            if ( ll < rr )
            {
              -- --rr;
              -- --rrhv;
            }
                              // Elements in [rr+1..end) >= pivot
                              // Pivot is in *ll
        }
        CGAL_assertion( llhv >= rrhv);
                              // rr + 1 >= ll >= rr
                              // Elements in [rr+1..end) >= pivot
                              // Elemente in [begin..ll) <  pivot
                              // Pivot is in a[ll]
        if ( ll == rr) {
            // Check for the possibly missed swap.
            if ( rr->is_border() && ! rr->opposite()->is_border()) {
                CGAL_assertion( rr + 1 == get_h_iter (rr->opposite()));
                std::swap( *rr, *(rr+1));
                update_opposite( rr);
                std::swap( *rrhv, *(rrhv+1));
            }
        }
        CGAL_assertion( ll->opposite()->is_border());
        CGAL_assertion( ll == halfedges.begin() || ! (ll-2)->is_border());
        CGAL_assertion( ll == halfedges.begin() || ! (ll-1)->is_border());
        border_halfedges = ll;
        nb_border_edges = (halfedges.end() - ll) / 2;
        nb_border_halfedges = 0;

        HVector inv_vector( halfedges.size());
        // Initialize inverse index.
        for ( Hiterator k = hvector.begin(); k != hvector.end(); ++k){
            inv_vector[*k - halfedges.begin()] =
                halfedges.begin() + (k - hvector.begin());
        }
        // Update halfedge pointers.
        HalfedgeDS_items_decorator<Self> D;
        for (Halfedge_iterator h = halfedges_begin();h != halfedges_end();++h){
            h->HBase::set_next( inv_vector[ h->next() - halfedges_begin()]);
            D.set_vertex_halfedge(h);
            if ( D.get_prev( h) == Halfedge_iterator())
                D.set_prev( h, Halfedge_iterator());
            else
                D.set_prev( h, inv_vector[ D.get_prev(h) - halfedges_begin()]);
            if ( h->is_border())
                nb_border_halfedges++;
            else
                D.set_face_halfedge(h);
        }
    }
};


} //namespace CGAL
#endif // CGAL_HALFEDGEDS_VECTOR_H //
// EOF //
