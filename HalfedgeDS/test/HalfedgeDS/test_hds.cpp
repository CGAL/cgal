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
// file          : test/HalfedgeDS/test_hds.C
// package       : HalfedgeDS 3.3 (27 Sep 2000)
// chapter       : Halfedge Data Structures
//
// revision      : $Id$
// revision_date : $Date$
//
// author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>
// maintainer    :
// coordinator   : MPI Saarbruecken
//
// Test hds.
// ============================================================================


#include <CGAL/Cartesian.h>
#include <CGAL/HalfedgeDS_default.h>
#include <CGAL/HalfedgeDS_min_items.h>
#include <CGAL/Polyhedron_items_3.h>
#include <CGAL/HalfedgeDS_decorator.h>
#include <cstddef>
#include <cassert>

typedef CGAL::Cartesian<double>   Kernel;
struct Dummy_traits_3 {
    typedef Kernel::Point_3    Point_3;
    typedef Kernel::Plane_3    Plane_3;
};
struct Dummy_traits_2 {
    typedef Kernel::Point_2    Point_2;
};

struct Empty_traits {};

void test_HalfedgeDS_polyhedron_default() {
    // Simple instantiation of the default for polyhedrons.
    typedef CGAL_HALFEDGEDS_DEFAULT< Dummy_traits_3,
        CGAL::Polyhedron_items_3> HDS;
    typedef HDS::Vertex           Vertex;
    typedef HDS::Halfedge         Halfedge;
    typedef HDS::Face             Face;
    typedef Halfedge::Base        HBase;

    HDS hds;
    hds.vertices_push_back( Vertex());
    hds.edges_push_back( Halfedge(), Halfedge());
    hds.faces_push_back( Face());
    hds.halfedges_begin()->HBase::set_face( hds.faces_begin());
    assert( hds.size_of_vertices() == 1);
    assert( hds.size_of_halfedges() == 2);
    assert( hds.size_of_faces() == 1);
    hds.normalize_border();
    assert( hds.size_of_border_halfedges() == 1);
    assert( hds.size_of_border_edges() == 1);
    typedef HDS::Halfedge_iterator Halfedge_iterator;
    Halfedge_iterator i = hds.halfedges_begin();
    assert(!(i == hds.halfedges_end()));
}
void test_HalfedgeDS_default() {
    // Simple instantiation of the default halfedge data structure.
    typedef CGAL_HALFEDGEDS_DEFAULT<Dummy_traits_2> HDS;
    typedef CGAL::HalfedgeDS_decorator<HDS> Decorator;
    HDS hds;
    Decorator D(hds);

    D.create_loop();
    assert( hds.size_of_vertices() == 1);
    assert( hds.size_of_halfedges() == 2);
    assert( hds.size_of_faces() == 2);
    assert( D.is_valid( false, 3));
    D.make_hole( hds.halfedges_begin()->opposite());
    hds.normalize_border();
    assert( hds.size_of_border_halfedges() == 1);
    assert( hds.size_of_border_edges() == 1);
    assert( D.is_valid( false, 4));
    D.fill_hole( hds.halfedges_begin()->opposite());
    hds.normalize_border();
    assert( hds.size_of_border_halfedges() == 0);
    assert( hds.size_of_border_edges() == 0);
    assert( D.is_valid( false, 4));

    HDS hds2(hds); // copy constructor.
    Decorator D2(hds2);
    assert( D2.is_valid( false, 4));

    hds = hds2; // assignment.
    assert( D.is_valid( false, 4));

    typedef HDS::Halfedge_iterator Halfedge_iterator;
    Halfedge_iterator i = hds.halfedges_begin();
    assert(!(i == hds.halfedges_end()));
}

#include <CGAL/HalfedgeDS_vector.h>

void test_HalfedgeDS_vector() {
    // Instantiation of the halfedge data structure using vector
    // and maximal bases for polyhedral surfaces.
    typedef CGAL::HalfedgeDS_vector< Dummy_traits_3,
                                           CGAL::Polyhedron_items_3> HDS;
    typedef HDS::Halfedge          Halfedge;
    typedef Halfedge::Base         HBase;
    typedef HDS::Face_handle       Face_handle;
    typedef CGAL::HalfedgeDS_decorator<HDS> Decorator;

    HDS hds(1,2,2);
    Decorator D(hds);

    D.create_loop();
    assert( hds.size_of_vertices() == 1);
    assert( hds.size_of_halfedges() == 2);
    assert( hds.size_of_faces() == 2);
    assert( D.is_valid( false, 3));
    (hds.halfedges_begin()+1)->HBase::set_face( Face_handle());
    hds.normalize_border();
    assert( hds.size_of_border_halfedges() == 1);
    assert( hds.size_of_border_edges() == 1);
    (hds.halfedges_begin()+1)->HBase::set_face( hds.faces_begin() + 1);
    hds.normalize_border();
    assert( hds.size_of_border_halfedges() == 0);
    assert( hds.size_of_border_edges() == 0);
    assert( D.is_valid( false, 4));
    hds.reserve(2,4,4);
    assert( hds.capacity_of_faces() == 4);
    assert( D.is_valid( false, 4));

    HDS hds2(hds); // copy constructor.
    Decorator D2(hds2);
    assert( D2.is_valid( false, 4));

    hds = hds2; // assignment.
    assert( D.is_valid( false, 4));
}

void test_HalfedgeDS_vector_min() {
    // Instantiation of the halfedge data structure using vector
    // and minimal bases as for an undirected graph.
    typedef CGAL::HalfedgeDS_vector< Empty_traits,
                                           CGAL::HalfedgeDS_min_items> HDS;
    typedef HDS::Halfedge   Halfedge;
    typedef Halfedge::Base  HBase;

    HDS hds(1,2,1);
    hds.edges_push_back( Halfedge(), Halfedge());
    hds.halfedges_begin()->HBase::set_next(hds.halfedges_begin()+1);
    (hds.halfedges_begin()+1)->HBase::set_next(hds.halfedges_begin());
    assert( hds.size_of_vertices() == 0);
    assert( hds.size_of_halfedges() == 2);
    assert( hds.size_of_faces() == 0);
    hds.normalize_border();
    assert( hds.size_of_border_halfedges() == 0);
    assert( hds.size_of_border_edges() == 0);
}


int main() {
    test_HalfedgeDS_polyhedron_default();
    test_HalfedgeDS_default();
    test_HalfedgeDS_vector();
    test_HalfedgeDS_vector_min();
    return 0;
}
// EOF //
