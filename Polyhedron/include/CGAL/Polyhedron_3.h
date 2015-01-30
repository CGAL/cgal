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

#ifndef CGAL_POLYHEDRON_3_H
#define CGAL_POLYHEDRON_3_H 1

#include <CGAL/basic.h>
#include <algorithm>
#include <cstddef>

#include <CGAL/HalfedgeDS_iterator.h>
#include <CGAL/Iterator_project.h>
#include <CGAL/function_objects.h>
#include <CGAL/N_step_adaptor_derived.h>
#include <CGAL/Polyhedron_items_3.h>
#include <CGAL/HalfedgeDS_default.h>
#include <CGAL/HalfedgeDS_const_decorator.h>
#include <CGAL/HalfedgeDS_decorator.h>
#include <CGAL/Modifier_base.h>
#include <CGAL/IO/Verbose_ostream.h>
#include <CGAL/Polyhedron_traits_3.h>

namespace CGAL {

template <class VertexBase>
class I_Polyhedron_vertex  : public VertexBase  {
public:
    typedef VertexBase                            Base;
    //typedef typename Base::HalfedgeDS              HDS;
    typedef typename Base::Point                   Point;
    typedef Point                                  Point_3;

    // Handles have to explicitly repeated, although they are derived
    typedef typename Base::Vertex_handle           Vertex_handle;
    typedef typename Base::Halfedge_handle         Halfedge_handle;
    typedef typename Base::Face_handle             Face_handle;
    typedef typename Base::Face_handle             Facet_handle;
    typedef typename Base::Vertex_const_handle     Vertex_const_handle;
    typedef typename Base::Halfedge_const_handle   Halfedge_const_handle;
    typedef typename Base::Face_const_handle       Face_const_handle;
    typedef typename Base::Face_const_handle       Facet_const_handle;
    typedef typename Base::Halfedge                Halfedge;
    typedef typename Base::Face                    Face;
    typedef typename Base::Face                    Facet;

    // Supported options by HDS.
    typedef typename Base::Supports_vertex_halfedge
                                                  Supports_vertex_halfedge;
    typedef typename Base::Supports_vertex_point  Supports_vertex_point;

    // Circulator category.
    typedef typename Halfedge::Supports_halfedge_prev  Supports_prev;

public:
    // Circulator category.
    typedef HalfedgeDS_circulator_traits<Supports_prev> Ctr;
    typedef typename Ctr::iterator_category circulator_category;

    // Circulators around a vertex and around a facet.
    typedef I_HalfedgeDS_facet_circ< Halfedge_handle, circulator_category>
                                         Halfedge_around_facet_circulator;

    typedef I_HalfedgeDS_vertex_circ< Halfedge_handle, circulator_category>
                                        Halfedge_around_vertex_circulator;

    typedef I_HalfedgeDS_facet_circ<
        Halfedge_const_handle,
        circulator_category>       Halfedge_around_facet_const_circulator;

    typedef I_HalfedgeDS_vertex_circ<
        Halfedge_const_handle,
        circulator_category>      Halfedge_around_vertex_const_circulator;



    typedef typename Halfedge_around_vertex_circulator::size_type
        size_type;
    typedef typename Halfedge_around_vertex_circulator::difference_type
        difference_type;

public:
    // We need to repeat the constructors here.
    I_Polyhedron_vertex() {}
    I_Polyhedron_vertex( const VertexBase& b) : VertexBase(b) {}
    I_Polyhedron_vertex( const Point_3& p) : VertexBase(p) {}

// New Access Functions (not provided in VertexBase).

    Halfedge_around_vertex_circulator vertex_begin() {
        // a circulator of halfedges around the vertex (clockwise).
        return Halfedge_around_vertex_circulator( this->halfedge());
    }
    Halfedge_around_vertex_const_circulator vertex_begin() const {
        // a circulator of halfedges around the vertex (clockwise).
        return Halfedge_around_vertex_const_circulator( this->halfedge());
    }

    // the degree of the vertex, i.e., edges emanating from this vertex
    std::size_t vertex_degree() const {
        return this->halfedge()->vertex_degree();
    }
    size_type degree() const { return vertex_degree(); } //backwards compatible

    // returns true if the vertex has exactly two incident edges
    bool is_bivalent() const { return  this->halfedge()->is_bivalent(); }

    // returns true if the vertex has exactly three incident edges
    bool is_trivalent() const { return  this->halfedge()->is_trivalent(); }

    // No longer hidden. Now the restricted version with precondition.
    // sets incident halfedge to h. Precondition: h is incident, i.e.,
    // h->vertex() == v.
    void  set_halfedge( Halfedge_handle hh) {
        CGAL_assertion( &*(hh->vertex()) == this);
        Base::set_halfedge(hh);
    }
};

// A halfedge is an oriented edge. Both orientations exist, i.e.
// an edge is represented by two opposite halfedges. The geometric
// relations are as follows:
//
//              _ _ _   .
//             /        |\.
//                      | \.
//           /             \ next half
//                          \ edge
//         /                 \.
//
//        |                   O  incident vertex
//                facet      ,
//        |                 /| |
//                         / | | opposite
//         \                 | | half edge
//                      half | |
//           \          edge | | /
//                           | |/
//             \_ _ _ _ _ _    '
//

template <class HalfedgeBase>
class I_Polyhedron_halfedge : public HalfedgeBase {
public:
    typedef HalfedgeBase                          Base;
    typedef typename Base::HalfedgeDS              HDS;

    // Handles have to explicitly repeated, although they are derived
    typedef typename Base::Vertex_handle           Vertex_handle;
    typedef typename Base::Halfedge_handle         Halfedge_handle;
    typedef typename Base::Face_handle             Face_handle;
    typedef typename Base::Face_handle             Facet_handle;
    typedef typename Base::Vertex_const_handle     Vertex_const_handle;
    typedef typename Base::Halfedge_const_handle   Halfedge_const_handle;
    typedef typename Base::Face_const_handle       Face_const_handle;
    typedef typename Base::Face_const_handle       Facet_const_handle;

    typedef typename Base::Vertex                  Vertex;
    typedef typename Base::Face                    Face;
    typedef typename Base::Face                    Facet;

    // Supported options by HDS.
    typedef typename Base::Supports_halfedge_prev Supports_halfedge_prev;
    typedef typename Base::Supports_halfedge_vertex
                                                  Supports_halfedge_vertex;
    typedef typename Base::Supports_halfedge_face Supports_halfedge_face;

    // Circulator category.
    typedef typename Base::Supports_halfedge_prev Supports_prev;

public:
    // Circulator category.
    typedef HalfedgeDS_circulator_traits<Supports_prev> Ctr;
    typedef typename Ctr::iterator_category circulator_category;

    // Circulators around a vertex and around a facet.
    typedef I_HalfedgeDS_facet_circ< Halfedge_handle, circulator_category>
                                         Halfedge_around_facet_circulator;

    typedef I_HalfedgeDS_vertex_circ< Halfedge_handle, circulator_category>
                                        Halfedge_around_vertex_circulator;

    typedef I_HalfedgeDS_facet_circ<
        Halfedge_const_handle,
        circulator_category>       Halfedge_around_facet_const_circulator;

    typedef I_HalfedgeDS_vertex_circ<
        Halfedge_const_handle,
        circulator_category>      Halfedge_around_vertex_const_circulator;



public:
    I_Polyhedron_halfedge() {}
    I_Polyhedron_halfedge( const HalfedgeBase& b) : HalfedgeBase(b) {}

// New Access Functions (not provided in HalfedgeBase).

    // Change semantic of prev: it is always available at this level.
    // If the HDS does not provide a prev-function, the previous
    // halfedge will be searched around the incident facet.
private:
    Halfedge_handle       find_prev( Halfedge_handle,       Tag_true) {
        return Base::prev();
    }
    Halfedge_const_handle find_prev( Halfedge_const_handle, Tag_true) const {
        return Base::prev();
    }
    Halfedge_handle find_prev( Halfedge_handle h, Tag_false) const {
        CGAL_precondition( &*h != this); // we have at least 2-gons
        while ( &*(h->next()) != this)
            h = h->next();
        return h;
    }
    Halfedge_const_handle find_prev( Halfedge_const_handle h, Tag_false) const{
        CGAL_precondition( &*h != this); // we have at least 2-gons
        while ( &*(h->next()) != this)
            h = h->next();
        return h;
    }

public:
    Halfedge_handle       prev() {
        return find_prev( this->next(), Supports_halfedge_prev());
    }
    Halfedge_const_handle prev() const {
        return find_prev( this->next(), Supports_halfedge_prev());
    }

    // Make face-functions also available as facet-functions.
    Face_handle           facet()       { return this->face();}
    Face_const_handle     facet() const { return this->face();}


    // the next halfedge around the vertex (clockwise). This is equal to
    // `h.next()->opposite()'.
    Halfedge_handle       next_on_vertex() { return this->next()->opposite(); }
    Halfedge_const_handle next_on_vertex() const {
        return this->next()->opposite();
    }

    // the previous halfedge around the vertex (counterclockwise). Is
    // equal to `h.opposite()->prev()'.
    Halfedge_handle       prev_on_vertex() { return this->opposite()->prev(); }
    Halfedge_const_handle prev_on_vertex() const {
        return this->opposite()->prev();
    }

    bool is_border_edge() const {
        // is true if `h' or `h.opposite()' is a border halfedge.
        return (this->opposite()->is_border() || this->is_border());
    }

    // a circulator of halfedges around the vertex (clockwise).
    Halfedge_around_vertex_circulator vertex_begin() {
        return Halfedge_around_vertex_circulator(
            HDS::halfedge_handle(this));
    }
    Halfedge_around_vertex_const_circulator vertex_begin() const {
        return Halfedge_around_vertex_const_circulator(
            HDS::halfedge_handle(this));
    }

    // a circulator of halfedges around the facet (counterclockwise).
    Halfedge_around_facet_circulator facet_begin() {
        return Halfedge_around_facet_circulator(
            HDS::halfedge_handle(this));
    }
    Halfedge_around_facet_const_circulator facet_begin() const {
        return Halfedge_around_facet_const_circulator(
            HDS::halfedge_handle(this));
    }

    // the degree of the incident vertex, i.e., edges emanating from this
    // vertex
    std::size_t vertex_degree() const {
        return circulator_size( vertex_begin());
    }

    // the degree of the incident facet, i.e., edges on the boundary of this
    // facet
    std::size_t facet_degree() const {
        return circulator_size( facet_begin());
    }

    // returns true if the incident vertex has exactly two incident edges
    bool is_bivalent() const {
        CGAL_precondition( this != &* (this->next()->opposite()));
        return  (this == &* (this->next()->opposite()->next()->opposite()));
    }

    // returns true if the incident vertex has exactly three incident edges
    bool is_trivalent() const {
        CGAL_precondition( this != &* (this->next()->opposite()));
        return  (   this != &* (this->next()->opposite()->next()->opposite())
                 && this == &* (this->next()->opposite()->next()->opposite()
                                ->next()->opposite()));
    }

    // returns true if the incident facet is a triangle.
    bool is_triangle() const {
        CGAL_precondition( this != &* (this->next()));
        CGAL_precondition( this != &* (this->next()->next()));
        return  (this == &* (this->next()->next()->next()));
    }

    // returns true if the incident facet is a quadrilateral.
    bool is_quad()     const {
        CGAL_precondition( this != &* (this->next()));
        CGAL_precondition( this != &* (this->next()->next()));
        return  (this == &* (this->next()->next()->next()->next()));
    }


private:
    // Hide some other functions of H.
    void  set_next( Halfedge_handle hh)  { Base::set_next(hh);}
    void  set_prev( Halfedge_handle hh)  { Base::set_prev(hh);}
    void  set_vertex( Vertex_handle vv)  { Base::set_vertex(vv);}
    void  set_face( Face_handle ff)      { Base::set_face(ff);}
    void  set_facet( Face_handle ff)     { set_face(ff);}
};


template <class FacetBase>
class I_Polyhedron_facet  : public FacetBase  {
public:
    typedef FacetBase                             Base;
    //typedef typename Base::HalfedgeDS              HDS;
    typedef typename Base::Plane                   Plane;
    typedef Plane                                  Plane_3;

    // Handles have to explicitly repeated, although they are derived
    typedef typename Base::Vertex_handle           Vertex_handle;
    typedef typename Base::Halfedge_handle         Halfedge_handle;
    typedef typename Base::Face_handle             Face_handle;
    typedef typename Base::Face_handle             Facet_handle;
    typedef typename Base::Vertex_const_handle     Vertex_const_handle;
    typedef typename Base::Halfedge_const_handle   Halfedge_const_handle;
    typedef typename Base::Face_const_handle       Face_const_handle;
    typedef typename Base::Face_const_handle       Facet_const_handle;
    typedef typename Base::Vertex                  Vertex;
    typedef typename Base::Halfedge                Halfedge;

    // Supported options by HDS.
    typedef typename Base::Supports_face_halfedge Supports_face_halfedge;
    typedef typename Base::Supports_face_plane    Supports_face_plane;

    // No longer required.
    typedef Tag_false                             Supports_face_normal;

    // Circulator category.
    typedef typename Halfedge::Supports_halfedge_prev  Supports_prev;

public:
    // Circulator category.
    typedef HalfedgeDS_circulator_traits<Supports_prev> Ctr;
    typedef typename Ctr::iterator_category circulator_category;

    // Circulators around a vertex and around a facet.
    typedef I_HalfedgeDS_facet_circ< Halfedge_handle, circulator_category>
                                         Halfedge_around_facet_circulator;

    typedef I_HalfedgeDS_vertex_circ< Halfedge_handle, circulator_category>
                                        Halfedge_around_vertex_circulator;

    typedef I_HalfedgeDS_facet_circ<
        Halfedge_const_handle,
        circulator_category>       Halfedge_around_facet_const_circulator;

    typedef I_HalfedgeDS_vertex_circ<
        Halfedge_const_handle,
        circulator_category>      Halfedge_around_vertex_const_circulator;



    typedef typename Halfedge_around_vertex_circulator::size_type
        size_type;
    typedef typename Halfedge_around_vertex_circulator::difference_type
        difference_type;

public:
    // We need to repeat the constructors here.
    I_Polyhedron_facet() {}
    I_Polyhedron_facet( const FacetBase& b) : FacetBase(b) {}

// New Access Functions (not provided in FacetBase).

    Halfedge_around_facet_circulator facet_begin() {
        // a circulator of halfedges around the facet (counterclockwise).
        return Halfedge_around_facet_circulator( this->halfedge());
    }
    Halfedge_around_facet_const_circulator facet_begin() const {
        // a circulator of halfedges around the facet (counterclockwise).
        return Halfedge_around_facet_const_circulator( this->halfedge());
    }

    // the degree of the incident facet, i.e., edges on the boundary of this
    // facet
    std::size_t facet_degree() const {return this->halfedge()->facet_degree();}
    size_type size() const { return facet_degree(); } // backwards compatible

    // returns true if the facet is a triangle.
    bool is_triangle() const { return this->halfedge()->is_triangle(); }

    // returns true if the facet is a quadrilateral.
    bool is_quad()     const { return this->halfedge()->is_quad(); }

    // No longer hidden. Now the restricted version with precondition.
    // sets incident halfedge to h. Precondition: h is incident, i.e.,
    // h->face() == v.
    void  set_halfedge( Halfedge_handle hh) {
        CGAL_assertion( &*(hh->facet()) == this);
        Base::set_halfedge(hh);
    }
};


template < class Items>
class I_Polyhedron_derived_items_3 {
public:
    template < class Refs, class Traits>
    class Vertex_wrapper {
    public:
        typedef typename Items::template Vertex_wrapper<Refs,Traits> VWrap;
        typedef typename VWrap::Vertex Vertex_base;
        typedef I_Polyhedron_vertex< Vertex_base> Vertex;
    };
    template < class Refs, class Traits>
    class Halfedge_wrapper {
    public:
        typedef typename Items::template Halfedge_wrapper<Refs,Traits> HWrap;
        typedef typename HWrap::Halfedge Halfedge_base;
        typedef I_Polyhedron_halfedge< Halfedge_base> Halfedge;
    };
    template < class Refs, class Traits>
    class Face_wrapper {
    public:
        typedef typename Items::template Face_wrapper<Refs,Traits> FWrap;
        typedef typename FWrap::Face Face_base;
        typedef I_Polyhedron_facet< Face_base> Face;
    };
};


template < class PolyhedronTraits_3,
           class PolyhedronItems_3 = Polyhedron_items_3,
           template < class T, class I, class A>
           class T_HDS = HalfedgeDS_default,
           class Alloc = CGAL_ALLOCATOR(int)>
class Polyhedron_3 {
    //
    // DEFINITION
    //
    // The boundary representation of a 3d-polyhedron P of the type
    // Polyhedron consists of vertices, edges and facets. The
    // vertices are points in space. The edges are straight line
    // segments. The facets are planar polygons. We restrict here
    // the facets to be simple planar polygons without holes and the
    // boundary of the polyhedron to be an oriented 2-manifold. Thus
    // facets are consistently oriented and an edge is incident to
    // exactly two facets. We restrict the representation further
    // that an edge has two distinct incident endpoints and
    // following duality that an edge has two distinct incident
    // facets. The class Polyhedron is able to guarantee
    // the combinatorial properties, but not all geometric
    // properties. Support functions are provided for testing
    // geometric properties, e.g. test for self intersections which
    // is  too expensive to be guaranteed as a class invariant.
public:
    typedef Polyhedron_3< PolyhedronTraits_3, PolyhedronItems_3, T_HDS, Alloc>
                                                  Self;
    typedef PolyhedronTraits_3                    Traits;
    typedef PolyhedronItems_3                     Items;
    typedef I_Polyhedron_derived_items_3<Items>   Derived_items;
    typedef T_HDS< Traits, Derived_items, Alloc>  HDS;
    typedef HDS                                   HalfedgeDS;

    // portability with older CGAL release
    typedef HDS                                   Halfedge_data_structure;

    typedef Alloc                                 Allocator;
    typedef Alloc                                 allocator_type; // STL name

    // Container stuff.
    typedef typename HDS::size_type               size_type;
    typedef typename HDS::difference_type         difference_type;
    typedef typename HDS::iterator_category       iterator_category;
    typedef typename HDS::Supports_removal        Supports_removal;

    // Geometry
    typedef typename Traits::Point_3              Point_3;
    typedef Point_3                               Point;
    typedef typename Traits::Plane_3              Plane_3;
    // No longer required.
    //typedef typename Traits::Normal               Normal;

    // Items
    typedef typename HDS::Vertex                  Vertex;
    typedef typename HDS::Halfedge                Halfedge;
    typedef typename HDS::Face                    Face;

    typedef typename Vertex::Base                 VBase;
    typedef typename Halfedge::Base               HBase;
    typedef typename Face::Base                   FBase;

    // Handles and Iterators
    typedef typename HDS::Vertex_handle           Vertex_handle;
    typedef typename HDS::Halfedge_handle         Halfedge_handle;
    typedef typename HDS::Face_handle             Face_handle;
    typedef typename HDS::Vertex_iterator         Vertex_iterator;
    typedef typename HDS::Halfedge_iterator       Halfedge_iterator;
    typedef typename HDS::Face_iterator           Face_iterator;

    typedef typename HDS::Vertex_const_handle     Vertex_const_handle;
    typedef typename HDS::Halfedge_const_handle   Halfedge_const_handle;
    typedef typename HDS::Face_const_handle       Face_const_handle;
    typedef typename HDS::Vertex_const_iterator   Vertex_const_iterator;
    typedef typename HDS::Halfedge_const_iterator Halfedge_const_iterator;
    typedef typename HDS::Face_const_iterator     Face_const_iterator;

    // Auxiliary iterators for convenience
    typedef Project_point<Vertex>                 Proj_point;
    typedef Iterator_project<Vertex_iterator, Proj_point>
                                                  Point_iterator;
    typedef Iterator_project<Vertex_const_iterator, Proj_point,
        const Point_3&, const Point_3*>           Point_const_iterator;

    typedef Project_plane<Face>                   Proj_plane;
    typedef Iterator_project<Face_iterator, Proj_plane>
                                                  Plane_iterator;
    typedef Iterator_project<Face_const_iterator, Proj_plane,
        const Plane_3&, const Plane_3*>           Plane_const_iterator;

    typedef N_step_adaptor_derived<Halfedge_iterator, 2>
                                                  Edge_iterator;
    typedef N_step_adaptor_derived<Halfedge_const_iterator, 2>
                                                  Edge_const_iterator;

    // All face related types get a related facet type name.
    typedef Face                                  Facet;
    typedef Face_handle                           Facet_handle;
    typedef Face_iterator                         Facet_iterator;
    typedef Face_const_handle                     Facet_const_handle;
    typedef Face_const_iterator                   Facet_const_iterator;

    // Supported options by HDS.
    typedef typename VBase::Supports_vertex_halfedge
                                                  Supports_vertex_halfedge;
    typedef typename HBase::Supports_halfedge_prev  Supports_halfedge_prev;
    typedef typename HBase::Supports_halfedge_prev  Supports_prev;
    typedef typename HBase::Supports_halfedge_vertex
                                                  Supports_halfedge_vertex;
    typedef typename HBase::Supports_halfedge_face  Supports_halfedge_face;
    typedef typename FBase::Supports_face_halfedge  Supports_face_halfedge;

    // Supported options especially for Polyhedron_3.
    typedef typename VBase::Supports_vertex_point   Supports_vertex_point;
    typedef typename FBase::Supports_face_plane     Supports_face_plane;

    // No longer required.
    typedef Tag_false                               Supports_face_normal;

    // Renamed versions for facet
    typedef Supports_halfedge_face  Supports_halfedge_facet;
    typedef Supports_face_halfedge  Supports_facet_halfedge;
    typedef Supports_face_plane     Supports_facet_plane;
    // No longer required.
    typedef Supports_face_normal    Supports_facet_normal;

public:
    // Circulator category.
    typedef HalfedgeDS_circulator_traits<Supports_prev> Ctr;
    typedef typename Ctr::iterator_category circulator_category;

    // Circulators around a vertex and around a facet.
    typedef I_HalfedgeDS_facet_circ< Halfedge_handle, circulator_category>
                                         Halfedge_around_facet_circulator;

    typedef I_HalfedgeDS_vertex_circ< Halfedge_handle, circulator_category>
                                        Halfedge_around_vertex_circulator;

    typedef I_HalfedgeDS_facet_circ<
        Halfedge_const_handle,
        circulator_category>       Halfedge_around_facet_const_circulator;

    typedef I_HalfedgeDS_vertex_circ<
        Halfedge_const_handle,
        circulator_category>      Halfedge_around_vertex_const_circulator;



protected:
    HDS     hds_;  // the boundary representation.
    Traits  m_traits;

public:
    HDS& hds() { return hds_; }
    const HDS& hds() const { return hds_; }

// CREATION
public:

    Polyhedron_3( const Traits& trts = Traits()) : m_traits(trts) {}
        // the empty polyhedron `P'.

    Polyhedron_3( size_type v, size_type h, size_type f,
                  const Traits& traits = Traits())
    : hds_(v,h,f), m_traits(traits) {}
        // a polyhedron `P' with storage reserved for v vertices, h
        // halfedges, and f facets. The reservation sizes are a hint for
        // optimizing storage allocation.

    void reserve( size_type v, size_type h, size_type f) {
        // reserve storage for v vertices, h halfedges, and f facets. The
        // reservation sizes are a hint for optimizing storage allocation.
        // If the `capacity' is already greater than the requested size
        // nothing happens. If the `capacity' changes all iterators and
        // circulators invalidates.
        hds_.reserve(v,h,f);
    }

protected:
    Halfedge_handle
    make_triangle( Vertex_handle v1, Vertex_handle v2, Vertex_handle v3) {
        HalfedgeDS_decorator<HDS> decorator(hds_);
        Halfedge_handle h  = hds_.edges_push_back( Halfedge(), Halfedge());
        h->HBase::set_next( hds_.edges_push_back( Halfedge(), Halfedge()));
        h->next()->HBase::set_next( hds_.edges_push_back( Halfedge(),
                                                         Halfedge()));
        h->next()->next()->HBase::set_next( h);
        decorator.set_prev( h, h->next()->next());
        decorator.set_prev( h->next(), h);
        decorator.set_prev( h->next()->next(), h->next());
        h->opposite()->HBase::set_next( h->next()->next()->opposite());
        h->next()->opposite()->HBase::set_next( h->opposite());
        h->next()->next()->opposite()->HBase::set_next(
            h->next()->opposite());
        decorator.set_prev( h->opposite(), h->next()->opposite());
        decorator.set_prev( h->next()->opposite(),
                            h->next()->next()->opposite());
        decorator.set_prev( h->next()->next()->opposite(), h->opposite());
        // the vertices
        decorator.set_vertex( h, v1);
        decorator.set_vertex( h->next(), v2);
        decorator.set_vertex( h->next()->next(), v3);
        decorator.set_vertex( h->opposite(), v3);
        decorator.set_vertex( h->next()->opposite(), v1);
        decorator.set_vertex( h->next()->next()->opposite(), v2);
        decorator.set_vertex_halfedge( h);
        decorator.set_vertex_halfedge( h->next());
        decorator.set_vertex_halfedge( h->next()->next());
        // the facet
        Facet_handle f = decorator.faces_push_back( Facet());
        decorator.set_face( h, f);
        decorator.set_face( h->next(), f);
        decorator.set_face( h->next()->next(), f);
        decorator.set_face_halfedge( h);
        return h;
    }

    Halfedge_handle
    make_tetrahedron( Vertex_handle v1,
                      Vertex_handle v2,
                      Vertex_handle v3,
                      Vertex_handle v4) {
        HalfedgeDS_decorator<HDS> decorator(hds_);
        Halfedge_handle h  = make_triangle(v1,v2,v3);
        // The remaining tip.
        Halfedge_handle g  = hds_.edges_push_back( Halfedge(), Halfedge());
        decorator.insert_tip( g->opposite(), h->opposite());
        decorator.close_tip( g);
        decorator.set_vertex( g, v4);
        Halfedge_handle e  = hds_.edges_push_back( Halfedge(), Halfedge());
        Halfedge_handle d  = hds_.edges_push_back( Halfedge(), Halfedge());
        decorator.insert_tip( e->opposite(), h->next()->opposite());
        decorator.insert_tip( e, g);
        decorator.insert_tip( d->opposite(),h->next()->next()->opposite());
        decorator.insert_tip( d, e);
        decorator.set_vertex_halfedge( g);
        // facets
        Facet_handle f = decorator.faces_push_back( Facet());
        decorator.set_face( h->opposite(), f);
        decorator.set_face( g, f);
        decorator.set_face( e->opposite(), f);
        decorator.set_face_halfedge( g);
        f = decorator.faces_push_back( Facet());
        decorator.set_face( h->next()->opposite(), f);
        decorator.set_face( e, f);
        decorator.set_face( d->opposite(), f);
        decorator.set_face_halfedge( e);
        f = decorator.faces_push_back( Facet());
        decorator.set_face( h->next()->next()->opposite(), f);
        decorator.set_face( d, f);
        decorator.set_face( g->opposite(), f);
        decorator.set_face_halfedge( d);
        return h;
    }

public:
    Halfedge_handle make_tetrahedron() {
        // the combinatorial structure of a tetrahedron is added to the
        // actual polyhedral surface. Returns an arbitrary halfedge of
        // this structure.
        reserve( 4 + size_of_vertices(),
                12 + size_of_halfedges(),
                 4 + size_of_facets());
        HalfedgeDS_decorator<HDS> decorator(hds_);
        return make_tetrahedron( decorator.vertices_push_back( Vertex()),
                                 decorator.vertices_push_back( Vertex()),
                                 decorator.vertices_push_back( Vertex()),
                                 decorator.vertices_push_back( Vertex()));
    }

    Halfedge_handle make_tetrahedron( const Point_3& p1,
                                      const Point_3& p2,
                                      const Point_3& p3,
                                      const Point_3& p4) {
        reserve( 4 + size_of_vertices(),
                12 + size_of_halfedges(),
                 4 + size_of_facets());
        HalfedgeDS_decorator<HDS> decorator(hds_);
        return make_tetrahedron( decorator.vertices_push_back( Vertex(p1)),
                                 decorator.vertices_push_back( Vertex(p2)),
                                 decorator.vertices_push_back( Vertex(p3)),
                                 decorator.vertices_push_back( Vertex(p4)));

    }

    Halfedge_handle make_triangle() {
        // the combinatorial structure of a single triangle with border
        // edges is added to the actual polyhedral surface. Returns an
        // arbitrary halfedge of this structure.
        reserve( 3 + size_of_vertices(),
                 6 + size_of_halfedges(),
                 1 + size_of_facets());
        HalfedgeDS_decorator<HDS> decorator(hds_);
        return make_triangle( decorator.vertices_push_back( Vertex()),
                              decorator.vertices_push_back( Vertex()),
                              decorator.vertices_push_back( Vertex()));
    }

    Halfedge_handle make_triangle( const Point_3& p1,
                                   const Point_3& p2,
                                   const Point_3& p3) {
        // the single triangle p_1, p_2, p_3 with border edges is added to
        // the actual polyhedral surface. Returns an arbitrary halfedge of
        // this structure.
        reserve( 3 + size_of_vertices(),
                 6 + size_of_halfedges(),
                 1 + size_of_facets());
        HalfedgeDS_decorator<HDS> decorator(hds_);
        return make_triangle( decorator.vertices_push_back( Vertex(p1)),
                              decorator.vertices_push_back( Vertex(p2)),
                              decorator.vertices_push_back( Vertex(p3)));
    }

// Access Member Functions

    allocator_type get_allocator() const { return hds_.get_allocator(); }

    size_type size_of_vertices() const { return hds_.size_of_vertices();}
        // number of vertices.

    size_type size_of_halfedges() const { return hds_.size_of_halfedges();}
        // number of all halfedges (including border halfedges).

    size_type size_of_facets() const { return hds_.size_of_faces();}
        // number of facets.

    bool empty() const { return size_of_halfedges() == 0; }

    bool is_empty() const { return size_of_halfedges() == 0; }

    size_type capacity_of_vertices() const {
        // space reserved for vertices.
        return hds_.capacity_of_vertices();
    }

    size_type capacity_of_halfedges() const {
        // space reserved for halfedges.
        return hds_.capacity_of_halfedges();
    }

    size_type capacity_of_facets() const {
        // space reserved for facets.
        return hds_.capacity_of_faces();
    }

    std::size_t bytes() const {
        // bytes used for the polyhedron.
        return sizeof(Self) - sizeof(HDS) + hds_.bytes();
    }

    std::size_t bytes_reserved() const {
        // bytes reserved for the polyhedron.
        return sizeof(Self) - sizeof(HDS) + hds_.bytes_reserved();
    }

    Vertex_iterator vertices_begin() { return hds_.vertices_begin();}
        // iterator over all vertices.

    Vertex_iterator vertices_end() { return hds_.vertices_end();}

    Halfedge_iterator halfedges_begin() { return hds_.halfedges_begin();}
        // iterator over all halfedges

    Halfedge_iterator halfedges_end() { return hds_.halfedges_end();}

    Facet_iterator facets_begin() { return hds_.faces_begin();}
        // iterator over all facets

    Facet_iterator facets_end() { return hds_.faces_end();}

    // The constant iterators and circulators.

    Vertex_const_iterator vertices_begin() const {
        return hds_.vertices_begin();
    }
    Vertex_const_iterator vertices_end() const {
        return hds_.vertices_end();
    }

    Halfedge_const_iterator halfedges_begin() const {
      return hds_.halfedges_begin();
    }
    Halfedge_const_iterator halfedges_end() const {
        return hds_.halfedges_end();
    }
    Facet_const_iterator facets_begin() const { return hds_.faces_begin();}
    Facet_const_iterator facets_end()   const { return hds_.faces_end();}

    // Auxiliary iterators for convinience
    Point_iterator       points_begin()       { return vertices_begin();}
    Point_iterator       points_end()         { return vertices_end();}

    Point_const_iterator points_begin() const { return vertices_begin();}
    Point_const_iterator points_end()   const { return vertices_end();}

    Plane_iterator       planes_begin()       { return facets_begin();}
    Plane_iterator       planes_end()         { return facets_end();}

    Plane_const_iterator planes_begin() const { return facets_begin();}
    Plane_const_iterator planes_end()   const { return facets_end();}

    Edge_iterator        edges_begin()        { return halfedges_begin();}
        // iterator over all edges. The iterator refers to halfedges, but
        // enumerates only one of the two corresponding opposite
        // halfedges.
    Edge_iterator        edges_end()          { return halfedges_end();}
        // end of the range over all edges.

    Edge_const_iterator  edges_begin()  const { return halfedges_begin();}
    Edge_const_iterator  edges_end()    const { return halfedges_end();}

    Traits&       traits()       { return m_traits; }
    const Traits& traits() const { return m_traits; }


// Combinatorial Predicates

    bool is_closed() const {
        for ( Halfedge_const_iterator i = halfedges_begin();
              i != halfedges_end(); ++i) {
            if ( i->is_border())
                return false;
        }
        return true;
    }

private:
    bool is_pure_bivalent( Tag_true) const {
        for ( Vertex_const_iterator i = vertices_begin();
              i != vertices_end(); ++i)
            if ( ! i->is_bivalent())
                return false;
        return true;
    }
    bool is_pure_bivalent( Tag_false) const {
        for ( Halfedge_const_iterator i = halfedges_begin();
              i != halfedges_end(); ++i)
            if ( ! i->is_bivalent())
                return false;
        return true;
    }

public:
    // returns true if all vertices have exactly two incident edges
    bool is_pure_bivalent() const {
        return is_pure_bivalent( Supports_vertex_halfedge());
    }

private:
    bool is_pure_trivalent( Tag_true) const {
        for ( Vertex_const_iterator i = vertices_begin();
              i != vertices_end(); ++i)
            if ( ! i->is_trivalent())
                return false;
        return true;
    }
    bool is_pure_trivalent( Tag_false) const {
        for ( Halfedge_const_iterator i = halfedges_begin();
              i != halfedges_end(); ++i)
            if ( ! i->is_trivalent())
                return false;
        return true;
    }

public:
    // returns true if all vertices have exactly three incident edges
    bool is_pure_trivalent() const {
        return is_pure_trivalent( Supports_vertex_halfedge());
    }

private:
    bool is_pure_triangle( Tag_true) const {
        for ( Facet_const_iterator i = facets_begin();
              i != facets_end(); ++i)
            if ( ! i->is_triangle())
                return false;
        return true;
    }
    bool is_pure_triangle( Tag_false) const {
        for ( Halfedge_const_iterator i = halfedges_begin();
              i != halfedges_end(); ++i)
            if ( ! i->is_border() && ! i->is_triangle())
                return false;
        return true;
    }

public:
    // returns true if all facets are triangles
    bool is_pure_triangle() const {
        return is_pure_triangle( Supports_facet_halfedge());
    }

private:
    bool is_pure_quad( Tag_true) const {
        for ( Facet_const_iterator i = facets_begin();
              i != facets_end(); ++i)
            if ( ! i->is_quad())
                return false;
        return true;
    }
    bool is_pure_quad( Tag_false) const {
        for ( Halfedge_const_iterator i = halfedges_begin();
              i != halfedges_end(); ++i)
            if ( ! i->is_border() && ! i->is_quad())
                return false;
        return true;
    }

public:
    // returns true if all facets are quadrilaterals
    bool is_pure_quad() const {
        return is_pure_quad( Supports_facet_halfedge());
    }


// Geometric Predicates

    bool
    is_triangle( Halfedge_const_handle h1) const {
        Halfedge_const_handle h2 = h1->next();
        Halfedge_const_handle h3 = h1->next()->next();
        // check halfedge combinatorics.
        // exact two edges at vertices 1, 2, 3.
        if ( h1->opposite()->next() != h3->opposite())    return false;
        if ( h2->opposite()->next() != h1->opposite())    return false;
        if ( h3->opposite()->next() != h2->opposite())    return false;
        // The facet is a triangle.
        if ( h1->next()->next()->next() != h1) return false;

        if ( check_tag( Supports_halfedge_face())
             &&  ! h1->is_border_edge())
            return false;  // implies h2 and h3
        CGAL_assertion( ! h1->is_border() || ! h1->opposite()->is_border());

        // Assert consistency.
        CGAL_assertion( h1 != h2);
        CGAL_assertion( h1 != h3);
        CGAL_assertion( h3 != h2);

        // check prev pointer.
        CGAL_assertion_code( HalfedgeDS_items_decorator<HDS> D;)
        CGAL_assertion(D.get_prev(h1) == Halfedge_handle() ||
                       D.get_prev(h1) == h3);
        CGAL_assertion(D.get_prev(h2) == Halfedge_handle() ||
                       D.get_prev(h2) == h1);
        CGAL_assertion(D.get_prev(h3) == Halfedge_handle() ||
                       D.get_prev(h3) == h2);

        // check vertices.
        CGAL_assertion( D.get_vertex(h1) == D.get_vertex( h2->opposite()));
        CGAL_assertion( D.get_vertex(h2) == D.get_vertex( h3->opposite()));
        CGAL_assertion( D.get_vertex(h3) == D.get_vertex( h1->opposite()));

        CGAL_assertion( ! check_tag( Supports_halfedge_vertex()) ||
                   D.get_vertex(h1) != D.get_vertex(h2));
        CGAL_assertion( ! check_tag( Supports_halfedge_vertex()) ||
                   D.get_vertex(h1) != D.get_vertex(h3));
        CGAL_assertion( ! check_tag( Supports_halfedge_vertex()) ||
                   D.get_vertex(h2) != D.get_vertex(h3));

        // check facets.
        CGAL_assertion( D.get_face(h1) == D.get_face(h2));
        CGAL_assertion( D.get_face(h1) == D.get_face(h3));

        return true;
    }

    bool
    is_tetrahedron( Halfedge_const_handle h1) const {
        Halfedge_const_handle h2 = h1->next();
        Halfedge_const_handle h3 = h1->next()->next();
        Halfedge_const_handle h4 = h1->opposite()->next();
        Halfedge_const_handle h5 = h2->opposite()->next();
        Halfedge_const_handle h6 = h3->opposite()->next();
        // check halfedge combinatorics.
        // at least three edges at vertices 1, 2, 3.
        if ( h4 == h3->opposite())    return false;
        if ( h5 == h1->opposite())    return false;
        if ( h6 == h2->opposite())    return false;
        // exact three edges at vertices 1, 2, 3.
        if ( h4->opposite()->next() != h3->opposite())    return false;
        if ( h5->opposite()->next() != h1->opposite())    return false;
        if ( h6->opposite()->next() != h2->opposite())    return false;
        // three edges at v4.
        if ( h4->next()->opposite() != h5) return false;
        if ( h5->next()->opposite() != h6) return false;
        if ( h6->next()->opposite() != h4) return false;
        // All facets are triangles.
        if ( h1->next()->next()->next() != h1) return false;
        if ( h4->next()->next()->next() != h4) return false;
        if ( h5->next()->next()->next() != h5) return false;
        if ( h6->next()->next()->next() != h6) return false;
        // all edges are non-border edges.
        if ( h1->is_border()) return false;  // implies h2 and h3
        CGAL_assertion( ! h2->is_border());
        CGAL_assertion( ! h3->is_border());
        if ( h4->is_border()) return false;
        if ( h5->is_border()) return false;
        if ( h6->is_border()) return false;

        // Assert consistency.
        CGAL_assertion( h1 != h2);
        CGAL_assertion( h1 != h3);
        CGAL_assertion( h3 != h2);
        CGAL_assertion( h4 != h5);
        CGAL_assertion( h5 != h6);
        CGAL_assertion( h6 != h4);

        // check prev pointer.
        CGAL_assertion_code( HalfedgeDS_items_decorator<HDS> D;)
        CGAL_assertion(D.get_prev(h1) == Halfedge_handle() ||
                       D.get_prev(h1) == h3);
        CGAL_assertion(D.get_prev(h2) == Halfedge_handle() ||
                       D.get_prev(h2) == h1);
        CGAL_assertion(D.get_prev(h3) == Halfedge_handle() ||
                       D.get_prev(h3) == h2);
        CGAL_assertion(D.get_prev(h4) == Halfedge_handle() ||
                  D.get_prev(h4) == h1->opposite());
        CGAL_assertion(D.get_prev(h5) == Halfedge_handle() ||
                  D.get_prev(h5) == h2->opposite());
        CGAL_assertion(D.get_prev(h6) == Halfedge_handle() ||
                  D.get_prev(h6) == h3->opposite());

        // check vertices.
        CGAL_assertion( D.get_vertex(h1) == D.get_vertex( h2->opposite()));
        CGAL_assertion( D.get_vertex(h1) == D.get_vertex( h5->opposite()));
        CGAL_assertion( D.get_vertex(h2) == D.get_vertex( h3->opposite()));
        CGAL_assertion( D.get_vertex(h2) == D.get_vertex( h6->opposite()));
        CGAL_assertion( D.get_vertex(h3) == D.get_vertex( h1->opposite()));
        CGAL_assertion( D.get_vertex(h3) == D.get_vertex( h4->opposite()));
        CGAL_assertion( D.get_vertex(h4) == D.get_vertex( h5));
        CGAL_assertion( D.get_vertex(h4) == D.get_vertex( h6));

        CGAL_assertion( ! check_tag( Supports_halfedge_vertex()) ||
                   D.get_vertex(h1) != D.get_vertex(h2));
        CGAL_assertion( ! check_tag( Supports_halfedge_vertex()) ||
                   D.get_vertex(h1) != D.get_vertex(h3));
        CGAL_assertion( ! check_tag( Supports_halfedge_vertex()) ||
                   D.get_vertex(h1) != D.get_vertex(h4));
        CGAL_assertion( ! check_tag( Supports_halfedge_vertex()) ||
                   D.get_vertex(h2) != D.get_vertex(h3));
        CGAL_assertion( ! check_tag( Supports_halfedge_vertex()) ||
                   D.get_vertex(h2) != D.get_vertex(h4));
        CGAL_assertion( ! check_tag( Supports_halfedge_vertex()) ||
                   D.get_vertex(h3) != D.get_vertex(h4));

        // check facets.
        CGAL_assertion( D.get_face(h1) == D.get_face(h2));
        CGAL_assertion( D.get_face(h1) == D.get_face(h3));
        CGAL_assertion( D.get_face(h4) == D.get_face(h4->next()));
        CGAL_assertion( D.get_face(h4) == D.get_face(h1->opposite()));
        CGAL_assertion( D.get_face(h5) == D.get_face(h5->next()));
        CGAL_assertion( D.get_face(h5) == D.get_face(h2->opposite()));
        CGAL_assertion( D.get_face(h6) == D.get_face(h6->next()));
        CGAL_assertion( D.get_face(h6) == D.get_face(h3->opposite()));

        CGAL_assertion( ! check_tag( Supports_halfedge_face()) ||
                   D.get_face(h1) != D.get_face(h4));
        CGAL_assertion( ! check_tag( Supports_halfedge_face()) ||
                   D.get_face(h1) != D.get_face(h5));
        CGAL_assertion( ! check_tag( Supports_halfedge_face()) ||
                   D.get_face(h1) != D.get_face(h6));
        CGAL_assertion( ! check_tag( Supports_halfedge_face()) ||
                   D.get_face(h4) != D.get_face(h5));
        CGAL_assertion( ! check_tag( Supports_halfedge_face()) ||
                   D.get_face(h4) != D.get_face(h6));
        CGAL_assertion( ! check_tag( Supports_halfedge_face()) ||
                   D.get_face(h5) != D.get_face(h6));

        return true;
    }

// Euler Operators (Combinatorial Modifications)
//
// The following Euler operations modify consistently the combinatorial
// structure of the polyhedral surface. The geometry remains unchanged.

    Halfedge_handle split_facet( Halfedge_handle h, Halfedge_handle g) {
        // split the facet incident to `h' and `g' into two facets with
        // new diagonal between the two vertices denoted by `h' and `g'
        // respectively. The second (new) facet is a copy of the first
        // facet. It returns the new diagonal. The time is proportional to
        // the distance from `h' to `g' around the facet. Precondition:
        // `h' and `g' are incident to the same facet. `h != g' (no
        // loops). `h->next() != g' and `g->next() != h' (no multi-edges).
        reserve( size_of_vertices(),
                 2 + size_of_halfedges(),
                 1 + size_of_facets());
        HalfedgeDS_decorator<HDS> D(hds_);
        CGAL_precondition( D.get_face(h) == D.get_face(g));
        CGAL_precondition( h != g);
        CGAL_precondition( h != g->next());
        CGAL_precondition( h->next() != g);
        return D.split_face( h, g);
    }

    Halfedge_handle join_facet( Halfedge_handle h) {
        // join the two facets incident to h. The facet incident to
        // `h->opposite()' gets removed. Both facets might be holes.
        // Returns the predecessor of h. The invariant `join_facet(
        // split_facet( h, g))' returns h and keeps the polyhedron
        // unchanged. The time is proportional to the size of the facet
        // removed and the time to compute `h.prev()'. Precondition:
        // `HDS' supports removal of facets. The degree of both
        // vertices incident to h is at least three (no antennas).
        HalfedgeDS_decorator<HDS> D(hds_);
        CGAL_precondition( circulator_size(h->vertex_begin())
                           >= size_type(3));
        CGAL_precondition( circulator_size(h->opposite()->vertex_begin())
                           >= size_type(3));
        return D.join_face(h);
    }

    Halfedge_handle split_vertex( Halfedge_handle h, Halfedge_handle g) {
        // split the vertex incident to `h' and `g' into two vertices and
        // connects them with a new edge. The second (new) vertex is a
        // copy of the first vertex. It returns the new edge. The time is
        // proportional to the distance from `h' to `g' around the vertex.
        // Precondition: `h' and `g' are incident to the same vertex. `h
        // != g' (no antennas). `h->next() != g' and `g->next() != h'.
        reserve( 1 + size_of_vertices(),
                 2 + size_of_halfedges(),
                 size_of_facets());
        HalfedgeDS_decorator<HDS> D(hds_);
        CGAL_precondition( D.get_vertex(h) == D.get_vertex(g));
        CGAL_precondition( h != g);
        return D.split_vertex( h, g);
    }

    Halfedge_handle join_vertex( Halfedge_handle h) {
        // join the two vertices incident to h. The vertex denoted by
        // `h->opposite()' gets removed. Returns the predecessor of h. The
        // invariant `join_vertex( split_vertex( h, g))' returns h and
        // keeps the polyhedron unchanged. The time is proportional to
        // the degree of the vertex removed and the time to compute
        // `h.prev()'.
        // Precondition: `HDS' supports removal of vertices. The size of
        // both facets incident to h is at least four (no multi-edges)
        HalfedgeDS_decorator<HDS> D(hds_);
        CGAL_precondition( circulator_size( h->facet_begin())
                           >= size_type(4));
        CGAL_precondition( circulator_size( h->opposite()->facet_begin())
                           >= size_type(4));
        return D.join_vertex(h);
    }

    Halfedge_handle split_edge( Halfedge_handle h) {
        return split_vertex( h->prev(), h->opposite())->opposite();
    }

    Halfedge_handle flip_edge( Halfedge_handle h) {
        HalfedgeDS_items_decorator<HDS> D;
        return D.flip_edge(h);
    }

    Halfedge_handle create_center_vertex( Halfedge_handle h) {
        HalfedgeDS_decorator<HDS> D(hds_);
        CGAL_assertion( circulator_size( h->facet_begin())
                        >= size_type(3));
        return D.create_center_vertex(h);
    }

    Halfedge_handle erase_center_vertex( Halfedge_handle h) {
        HalfedgeDS_decorator<HDS> D(hds_);
        return D.erase_center_vertex(h);
    }

// Euler Operators Modifying Genus

    Halfedge_handle split_loop( Halfedge_handle h,
                                Halfedge_handle i,
                                Halfedge_handle j) {
        // cut the polyhedron into two parts along the cycle (h,i,j).
        // Three copies of the vertices and two new triangles will be
        // created. h,i,j will be incident to the first new triangle. The
        // returnvalue will be an halfedge iterator denoting the new
        // halfegdes of the second new triangle which was h beforehand.
        // Precondition: h,i,j are distinct, consecutive vertices of the
        // polyhedron and form a cycle: i.e. `h->vertex() == i->opposite()
        // ->vertex()', ..., `j->vertex() == h->opposite()->vertex()'. The
        // six facets incident to h,i,j are all distinct.
        reserve( 3 + size_of_vertices(),
                 6 + size_of_halfedges(),
                 2 + size_of_facets());
        HalfedgeDS_decorator<HDS> D(hds_);
        CGAL_precondition( h != i);
        CGAL_precondition( h != j);
        CGAL_precondition( i != j);
        CGAL_precondition( D.get_vertex(h) == D.get_vertex(i->opposite()));
        CGAL_precondition( D.get_vertex(i) == D.get_vertex(j->opposite()));
        CGAL_precondition( D.get_vertex(j) == D.get_vertex(h->opposite()));
        CGAL_precondition( D.get_face(h) == Facet_handle() ||
                           D.get_face(h) != D.get_face(i));
        CGAL_precondition( D.get_face(h) == Facet_handle() ||
                           D.get_face(h) != D.get_face(j));
        CGAL_precondition( D.get_face(i) == Facet_handle() ||
                           D.get_face(i) != D.get_face(j));
        CGAL_precondition( D.get_face(h) == Facet_handle() ||
                           D.get_face(h) != D.get_face(h->opposite()));
        CGAL_precondition( D.get_face(h) == Facet_handle() ||
                           D.get_face(h) != D.get_face(i->opposite()));
        CGAL_precondition( D.get_face(h) == Facet_handle() ||
                           D.get_face(h) != D.get_face(j->opposite()));
        CGAL_precondition( D.get_face(i) == Facet_handle() ||
                           D.get_face(i) != D.get_face(h->opposite()));
        CGAL_precondition( D.get_face(i) == Facet_handle() ||
                           D.get_face(i) != D.get_face(i->opposite()));
        CGAL_precondition( D.get_face(i) == Facet_handle() ||
                           D.get_face(i) != D.get_face(j->opposite()));
        CGAL_precondition( D.get_face(j) == Facet_handle() ||
                           D.get_face(j) != D.get_face(h->opposite()));
        CGAL_precondition( D.get_face(j) == Facet_handle() ||
                           D.get_face(j) != D.get_face(i->opposite()));
        CGAL_precondition( D.get_face(j) == Facet_handle() ||
                           D.get_face(j) != D.get_face(j->opposite()));
        CGAL_precondition( D.get_face(h->opposite()) == Facet_handle() ||
            D.get_face(h->opposite()) != D.get_face(i->opposite()));
        CGAL_precondition( D.get_face(h->opposite()) == Facet_handle() ||
            D.get_face(h->opposite()) != D.get_face(j->opposite()));
        CGAL_precondition( D.get_face(i->opposite()) == Facet_handle() ||
            D.get_face(i->opposite()) != D.get_face(j->opposite()));
        return D.split_loop( h, i, j);
    }

    Halfedge_handle join_loop( Halfedge_handle h, Halfedge_handle g) {
        // glues the boundary of two facets together. Both facets and the
        // vertices of g gets removed. Returns an halfedge iterator for h.
        // The invariant `join_loop( h, split_loop( h, i, j))' returns h
        // and keeps the polyhedron unchanged. Precondition: `HDS'
        // supports removal of vertices and facets. The facets denoted by
        // h and g have equal size.
        HalfedgeDS_decorator<HDS> D(hds_);
        CGAL_precondition( D.get_face(h) == Facet_handle() ||
                           D.get_face(h) != D.get_face(g));
        CGAL_precondition( circulator_size( h->facet_begin())
                           >= size_type(3));
        CGAL_precondition( circulator_size( h->facet_begin())
                           == circulator_size( g->facet_begin()));
        return D.join_loop( h, g);
    }

// Modifying Facets and Holes

    Halfedge_handle make_hole( Halfedge_handle h) {
        // removes incident facet and makes all halfedges incident to the
        // facet to border edges. Returns h. Precondition: `HDS'
        // supports removal of facets. `! h.is_border()'.
        HalfedgeDS_decorator<HDS> D(hds_);
        return D.make_hole(h);
    }

    Halfedge_handle fill_hole( Halfedge_handle h) {
        // fill a hole with a new created facet. Makes all border
        // halfedges of the hole denoted by h incident to the new facet.
        // Returns h. Precondition: `h.is_border()'.
        reserve( size_of_vertices(),
                 size_of_halfedges(),
                 1 + size_of_facets());
        HalfedgeDS_decorator<HDS> D(hds_);
        return D.fill_hole(h);
    }

    Halfedge_handle add_vertex_and_facet_to_border( Halfedge_handle h,
                                                    Halfedge_handle g) {
        // creates a new facet within the hole incident to h and g by
        // connecting the tip of g with the tip of h with two new
        // halfedges and a new vertex and filling this separated part of
        // the hole with a new facet. Returns the new halfedge incident to
        // the new facet and the new vertex. Precondition: `h->is_border(
        // )', `g->is_border()', `h != g', and g can be reached along the
        // same hole starting with h.
        CGAL_precondition( h != g);
        reserve( 1 + size_of_vertices(),
                 4 + size_of_halfedges(),
                 1 + size_of_facets());
        HalfedgeDS_decorator<HDS> D(hds_);
        Halfedge_handle hh = D.add_face_to_border( h, g);
        CGAL_assertion( hh == g->next());
        D.split_vertex( g, hh->opposite());
        return g->next();
    }

    Halfedge_handle add_facet_to_border( Halfedge_handle h,
                                         Halfedge_handle g) {
        // creates a new facet within the hole incident to h and g by
        // connecting the tip of g with the tip of h with a new halfedge
        // and filling this separated part of the hole with a new facet.
        // Returns the new halfedge incident to the new facet.
        // Precondition: `h->is_border()', `g->is_border()', `h != g',
        // `h->next() != g', and g can be reached along the same hole
        // starting with h.
        CGAL_precondition( h != g);
        CGAL_precondition( h->next() != g);
        reserve( size_of_vertices(),
                 2 + size_of_halfedges(),
                 1 + size_of_facets());
        HalfedgeDS_decorator<HDS> D(hds_);
        return D.add_face_to_border( h, g);
    }

// Erasing

    void erase_facet( Halfedge_handle h) {
        // removes the incident facet of h and changes all halfedges
        // incident to the facet into border edges or removes them from
        // the polyhedral surface if they were already border edges. See
        // `make_hole(h)' for a more specialized variant. Precondition:
        // `Traits' supports removal.
        HalfedgeDS_decorator<HDS> D(hds_);
        D.erase_face(h);
    }

    void erase_connected_component( Halfedge_handle h) {
        // removes the vertices, halfedges, and facets that belong to the
        // connected component of h. Precondition: `Traits' supports
        // removal.
        HalfedgeDS_decorator<HDS> D(hds_);
        D.erase_connected_component(h);
    }

    /// Erases the small connected components and the isolated vertices.
    ///
    /// @commentheading Preconditions:
    /// supports vertices, halfedges, and removal operation.
    ///
    /// @commentheading Template Parameters:
    /// @param nb_components_to_keep the number of large connected components to keep.
    ///
    /// @return the number of connected components erased (ignoring isolated vertices).
    unsigned int keep_largest_connected_components(unsigned int nb_components_to_keep)
    {
        HalfedgeDS_decorator<HDS> D(hds_);
        return D.keep_largest_connected_components(nb_components_to_keep);
    }

    void clear() { hds_.clear(); }
        // removes all vertices, halfedges, and facets.

    void erase_all() { clear(); }
        // equivalent to `clear()'. Depricated.

// Special Operations on Polyhedral Surfaces

    void delegate( Modifier_base<HDS>& modifier) {
        // calls the `operator()' of the `modifier'. Precondition: The
        // `modifier' returns a consistent representation.
        modifier( hds_);
        CGAL_expensive_postcondition( is_valid());
    }

// Operations with Border Halfedges

    size_type size_of_border_halfedges() const {
        // number of border halfedges. An edge with no incident facet
        // counts as two border halfedges. Precondition: `normalize_border
        // ()' has been called and no halfedge insertion or removal and no
        // change in border status of the halfedges have occured since
        // then.
        return hds_.size_of_border_halfedges();
    }

    size_type size_of_border_edges() const {
        // number of border edges. If `size_of_border_edges() ==
        // size_of_border_halfedges()' all border edges are incident to a
        // facet on one side and to a hole on the other side.
        // Precondition: `normalize_border()' has been called and no
        // halfedge insertion or removal and no change in border status of
        // the halfedges have occured since then.
        return hds_.size_of_border_edges();
    }

    Halfedge_iterator border_halfedges_begin() {
        // halfedge iterator starting with the border edges. The range [
        // `halfedges_begin(), border_halfedges_begin()') denotes all
        // non-border edges. The range [`border_halfedges_begin(),
        // halfedges_end()') denotes all border edges. Precondition:
        // `normalize_border()' has been called and no halfedge insertion
        // or removal and no change in border status of the halfedges have
        // occured since then.
        return hds_.border_halfedges_begin();
    }
    Halfedge_const_iterator border_halfedges_begin() const {
        return hds_.border_halfedges_begin();
    }

    // Convenient edge iterator
    Edge_iterator border_edges_begin() { return border_halfedges_begin(); }
    Edge_const_iterator border_edges_begin() const {
        return border_halfedges_begin();
    }

    bool normalized_border_is_valid( bool verbose = false) const {
        // checks whether all non-border edges precedes the border edges.
        HalfedgeDS_const_decorator<HDS> decorator(hds_);
        bool valid = decorator.normalized_border_is_valid( verbose);
        for ( Halfedge_const_iterator i = border_halfedges_begin();
              valid && (i != halfedges_end()); (++i, ++i)) {
            if ( i->is_border()) {
                Verbose_ostream verr(verbose);
                verr << "    both halfedges of an edge are border "
                        "halfedges." << std::endl;
                valid = false;
            }
        }
        return valid;
    }

    void normalize_border() {
        // sorts halfedges such that the non-border edges precedes the
        // border edges.
        hds_.normalize_border();
        CGAL_postcondition( normalized_border_is_valid());
    }

protected:            // Supports_face_plane
    void inside_out_geometry( Tag_false) {}
    void inside_out_geometry( Tag_true) {
        typename Traits::Construct_opposite_plane_3 opp
            = traits().construct_opposite_plane_3_object();
        std::transform( planes_begin(), planes_end(), planes_begin(), opp);
    }

public:
    void inside_out() {
        // reverse facet orientation.
        HalfedgeDS_decorator<HDS> decorator(hds_);
        decorator.inside_out();
        inside_out_geometry( Supports_face_plane());
    }

    bool is_valid( bool verb = false, int level = 0) const {
        // checks the combinatorial consistency.
        Verbose_ostream verr(verb);
        verr << "begin CGAL::Polyhedron_3<...>::is_valid( verb=true, "
                          "level = " << level << "):" << std::endl;
        HalfedgeDS_const_decorator<HDS> D(hds_);
        bool valid = D.is_valid( verb, level + 3);
        // All halfedges.
        Halfedge_const_iterator i   = halfedges_begin();
        Halfedge_const_iterator end = halfedges_end();
        size_type  n = 0;
        for( ; valid && (i != end); ++i) {
            verr << "halfedge " << n << std::endl;
            // At least triangular facets and distinct geometry.
            valid = valid && ( i->next() != i);
            valid = valid && ( i->next()->next() != i);
            valid = valid && ( ! check_tag( Supports_halfedge_vertex()) ||
                               D.get_vertex(i) != D.get_vertex(i->opposite()));
            valid = valid && ( ! check_tag( Supports_halfedge_vertex()) ||
                               D.get_vertex(i) != D.get_vertex(i->next()));
            valid = valid && ( ! check_tag( Supports_halfedge_vertex()) ||
                        D.get_vertex(i) != D.get_vertex(i->next()->next()));
            if ( ! valid) {
                verr << "    incident facet is not at least a triangle."
                     << std::endl;
                break;
            }
            // Distinct facets on each side of an halfegde.
            valid = valid && ( ! check_tag( Supports_halfedge_face()) ||
                               D.get_face(i) != D.get_face(i->opposite()));
            if ( ! valid) {
                verr << "    both incident facets are equal." << std::endl;
                break;
            }
            ++n;
        }
        valid = valid && (n == size_of_halfedges());
        if ( n != size_of_halfedges())
            verr << "counting halfedges failed." << std::endl;

        verr << "end of CGAL::Polyhedron_3<...>::is_valid(): structure is "
             << ( valid ? "valid." : "NOT VALID.") << std::endl;
        return valid;
    }
};

} //namespace CGAL

#ifndef CGAL_NO_DEPRECATED_CODE
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#endif

#endif // CGAL_POLYHEDRON_3_H //
