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
// file          : Halfedge_data_structure_bases.h
// chapter       : $CGAL_Chapter: Halfedge Data Structures $
// package       : $CGAL_Package: Halfedge_DS 2.8 (13 Sep 2000) $
// source        : hds.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : MPI Saarbruecken (Stefan Schirra <stschirr@mpi-sb.mpg.de>)
//
// Halfedge Data Structure Base Classes for Vertices, Edges, Facets.
// ============================================================================

#ifndef CGAL_HALFEDGE_DATA_STRUCTURE_BASES_H
#define CGAL_HALFEDGE_DATA_STRUCTURE_BASES_H 1

#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

class Vertex_min_base {
    // defines the minimal vertex functionality with no data at all.
public:
    typedef Tag_false Supports_vertex_point;
    typedef Tag_false Supports_vertex_halfedge;
    typedef void*     Point;
};

class Halfedge_min_base {
    // defines the minimal halfedge functionality with next and opposite
    // pointers.
protected:
    void* opp;
    void* nxt;
public:
    typedef Tag_false Supports_halfedge_prev;
    typedef Tag_false Supports_halfedge_vertex;
    typedef Tag_false Supports_halfedge_facet;

    void*       opposite()       { return opp;}
    const void* opposite() const { return opp;}
    void*       next()           { return nxt;}
    const void* next() const     { return nxt;}
        // the next halfedge along the facet.

    bool is_border() const { return false;}
        // is always false since facets are not supported here.

    void  set_opposite( void* h)  { opp = h;}
    void  set_next( void* h)      { nxt = h;}
};

class Facet_min_base {
    // defines the minimal facet functionality with no data at all.
public:
    typedef Tag_false Supports_facet_plane;
    typedef Tag_false Supports_facet_normal;
    typedef Tag_false Supports_facet_halfedge;
    typedef void*     Plane;
    typedef void*     Normal;
};

template <class Pt>
class Vertex_max_base {
    // defines the maximal vertex functionality including halfedge pointer
    // and a template parameter for the point.
protected:
    void* hdg;
    Pt    pt;
public:
    typedef Tag_true  Supports_vertex_point;
    typedef Tag_true  Supports_vertex_halfedge;
    typedef Pt        Point;

    Vertex_max_base() {}
    Vertex_max_base( const Pt& p) : pt(p) {}

    void*       halfedge()               { return hdg;}
    const void* halfedge() const         { return hdg;}
    void        set_halfedge( void* h)   { hdg = h;}
        // an incident halfedge pointing to `v'.

    Point&       point()       { return pt;}
    const Point& point() const { return pt;}
};

class Halfedge_max_base : public  Halfedge_min_base{
    // defines the maximal halfedge functionality including previous,
    // vertex and facet pointers.
protected:
    void* prv;
    void* v;
    void* f;
public:
    typedef Halfedge_min_base Base;
    typedef Tag_true  Supports_halfedge_prev;
    typedef Tag_true  Supports_halfedge_vertex;
    typedef Tag_true  Supports_halfedge_facet;

    Halfedge_max_base() : f(NULL) {}

    void*       prev()       { return prv;}
    const void* prev() const { return prv;}
        // the previous halfedge along the facet.

    void*       vertex()       { return v;}
    const void* vertex() const { return v;}
        // the incident vertex.

    void*       facet()       { return f;}
    const void* facet() const { return f;}
        // the facet to the left.

    bool is_border() const { return f == NULL;}
        // is true if `h' is a border halfedge).

    void  set_prev( void* h)        { prv = h;}
    void  set_vertex( void* _v)     { v = _v;}
    void  set_facet( void* _f)      { f = _f;}

    // Avoids unnecessary matchings with base class. (g++ 2.7.2 bug)
    void*       opposite()       { return Base::opposite();}
    const void* opposite() const { return Base::opposite();}
    void*       next()           { return Base::next();}
    const void* next() const     { return Base::next();}
};

class Facet_max_base {
    // defines the maximal facet functionality including halfedge pointer
    // but no geometry.
protected:
    void* hdg;
public:
    typedef Tag_false Supports_facet_plane;
    typedef Tag_false Supports_facet_normal;
    typedef Tag_true  Supports_facet_halfedge;
    typedef void*     Plane;
    typedef void*     Normal;

    void*       halfedge()       { return hdg;}
    const void* halfedge() const { return hdg;}
        // an incident halfedge pointing to `f'.

    void set_halfedge( void* h)  { hdg = h;}
};

template < class Kernel_>
class Polyhedron_facet_base_3 : public Facet_max_base {
    // defines the maximal facet functionality including halfedge pointer,
    // plane equation and normal vector.
public:
    typedef Kernel_                    Kernel;
    // typedef Kernel            R;  // maybe for backwards compatibility
    typedef typename Kernel::Vector_3  Vector_3;
    typedef typename Kernel::Plane_3   Plane_3;
    typedef Tag_true                   Supports_facet_plane;
    typedef Tag_true                   Supports_facet_normal;
protected:
    Plane_3 pln;
public:
    Vector_3        normal() const { return pln.orthogonal_vector();}
    Plane_3&        plane()        { return pln;}
    const Plane_3&  plane() const  { return pln;}
};

CGAL_END_NAMESPACE
#endif // CGAL_HALFEDGE_DATA_STRUCTURE_BASES_H //
// EOF //
