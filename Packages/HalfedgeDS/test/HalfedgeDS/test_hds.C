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
// file          : test_hds.C
// chapter       : $CGAL_Chapter: Halfedge Data Structures $
// package       : $CGAL_Package: HalfedgeDS 3.3 (27 Sep 2000) $
// source        : hds.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : MPI Saarbruecken (Stefan Schirra <stschirr@mpi-sb.mpg.de>)
//
// Test hds.
// ============================================================================


#include <CGAL/basic.h>
#include <cstddef>

#include <CGAL/Cartesian.h>
#include <CGAL/HalfedgeDS_default.h>
#include <CGAL/HalfedgeDS_min_items.h>
#include <CGAL/Polyhedron_items_3.h>
#include <CGAL/HalfedgeDS_decorator.h>
#include <CGAL/Point_2.h>
#include <CGAL/Point_3.h>
#include <CGAL/Plane_3.h>

typedef CGAL::Cartesian<double>   Rep;
struct Dummy_traits_3 {
    typedef CGAL::Point_3<Rep>    Point;
    typedef CGAL::Vector_3<Rep>   Normal;
    typedef CGAL::Plane_3<Rep>    Plane;
};
struct Dummy_traits_2 {
    typedef CGAL::Point_2<Rep>    Point;
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
    CGAL_assertion( hds.size_of_vertices() == 1);
    CGAL_assertion( hds.size_of_halfedges() == 2);
    CGAL_assertion( hds.size_of_faces() == 1);
    hds.normalize_border();
    CGAL_assertion( hds.size_of_border_halfedges() == 1);
    CGAL_assertion( hds.size_of_border_edges() == 1);
    typedef HDS::Halfedge_iterator Halfedge_iterator;
    Halfedge_iterator i = hds.halfedges_begin();
    CGAL_assertion(!(i == hds.halfedges_end()));
}
void test_HalfedgeDS_default() {
    // Simple instantiation of the default halfedge data structure.
    typedef CGAL_HALFEDGEDS_DEFAULT<Dummy_traits_2> HDS;
    typedef CGAL::HalfedgeDS_decorator<HDS> Decorator;
    HDS hds;
    Decorator D(hds);

    D.create_loop();
    CGAL_assertion( hds.size_of_vertices() == 1);
    CGAL_assertion( hds.size_of_halfedges() == 2);
    CGAL_assertion( hds.size_of_faces() == 2);
    CGAL_assertion( D.is_valid( false, 3));
    D.make_hole( hds.halfedges_begin()->opposite());
    hds.normalize_border();
    CGAL_assertion( hds.size_of_border_halfedges() == 1);
    CGAL_assertion( hds.size_of_border_edges() == 1);
    CGAL_assertion( D.is_valid( false, 4));
    D.fill_hole( hds.halfedges_begin()->opposite());
    hds.normalize_border();
    CGAL_assertion( hds.size_of_border_halfedges() == 0);
    CGAL_assertion( hds.size_of_border_edges() == 0);
    CGAL_assertion( D.is_valid( false, 4));

    HDS hds2(hds); // copy constructor.
    Decorator D2(hds2);
    CGAL_assertion( D2.is_valid( false, 4));

    hds = hds2; // assignment.
    CGAL_assertion( D.is_valid( false, 4));

    typedef HDS::Halfedge_iterator Halfedge_iterator;
    Halfedge_iterator i = hds.halfedges_begin();
    CGAL_assertion(!(i == hds.halfedges_end()));
}
#ifndef CGAL_HALFEDGEDS_USING_LIST_H
#include <CGAL/HalfedgeDS_using_list.h>
#endif

void test_HalfedgeDS_using_list() {
    // Instantiation of the halfedge data structure using list
    // and maximal bases for polyhedral surfaces.
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
    typedef CGAL::HalfedgeDS_using_list< Dummy_traits_3,
                                         CGAL::Polyhedron_items_3> HDS;
#else
    typedef CGAL::HalfedgeDS_using_list::HDS< Dummy_traits_3,
                                         CGAL::Polyhedron_items_3> HDS;
#endif
    typedef HDS::Halfedge          Halfedge;
    typedef Halfedge::Base         HBase;
    typedef HDS::Face_handle       Face_handle;
    typedef CGAL::HalfedgeDS_decorator<HDS> Decorator;

    HDS hds;
    Decorator D(hds);

    D.create_loop();
    CGAL_assertion( hds.size_of_vertices() == 1);
    CGAL_assertion( hds.size_of_halfedges() == 2);
    CGAL_assertion( hds.size_of_faces() == 2);
    CGAL_assertion( D.is_valid( false, 3));
    D.make_hole( hds.halfedges_begin()->opposite());
    hds.normalize_border();
    CGAL_assertion( hds.size_of_border_halfedges() == 1);
    CGAL_assertion( hds.size_of_border_edges() == 1);
    CGAL_assertion( D.is_valid( false, 4));
    D.fill_hole( hds.halfedges_begin()->opposite());
    hds.normalize_border();
    CGAL_assertion( hds.size_of_border_halfedges() == 0);
    CGAL_assertion( hds.size_of_border_edges() == 0);
    CGAL_assertion( D.is_valid( false, 4));

    HDS hds2(hds); // copy constructor.
    Decorator D2(hds2);
    CGAL_assertion( D2.is_valid( false, 4));

    hds = hds2; // assignment.
    CGAL_assertion( D.is_valid( false, 4));

    typedef HDS::Halfedge_iterator Halfedge_iterator;
    Halfedge_iterator i = hds.halfedges_begin();
    CGAL_assertion(!(i == hds.halfedges_end()));
}

void test_HalfedgeDS_using_list_min() {
    // Instantiation of the halfedge data structure using vector
    // and minimal bases as for an undirected graph.
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
    typedef CGAL::HalfedgeDS_using_list< Empty_traits,
                                         CGAL::HalfedgeDS_min_items> HDS;
#else
    typedef CGAL::HalfedgeDS_using_list::HDS< Empty_traits,
                                         CGAL::HalfedgeDS_min_items> HDS;
#endif
    typedef HDS::Halfedge            Halfedge;
    typedef HDS::Halfedge_iterator   Halfedge_iterator;
    typedef Halfedge::Base           HBase;

    HDS hds;
    hds.edges_push_back( Halfedge(), Halfedge());
    Halfedge_iterator i = hds.halfedges_begin();
    Halfedge_iterator j = i;
    ++j;
    i->HBase::set_next(j);
    j->HBase::set_next(i);
    CGAL_assertion( hds.size_of_vertices() == 0);
    CGAL_assertion( hds.size_of_halfedges() == 2);
    CGAL_assertion( hds.size_of_faces() == 0);
    hds.normalize_border();
    CGAL_assertion( hds.size_of_border_halfedges() == 0);
    CGAL_assertion( hds.size_of_border_edges() == 0);
}
#ifndef CGAL_HALFEDGEDS_USING_VECTOR_H
#include <CGAL/HalfedgeDS_using_vector.h>
#endif

void test_HalfedgeDS_using_vector() {
    // Instantiation of the halfedge data structure using vector
    // and maximal bases for polyhedral surfaces.
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
    typedef CGAL::HalfedgeDS_using_vector< Dummy_traits_3,
                                           CGAL::Polyhedron_items_3> HDS;
#else
    typedef CGAL::HalfedgeDS_using_vector::HDS< Dummy_traits_3,
                                           CGAL::Polyhedron_items_3> HDS;
#endif
    typedef HDS::Halfedge          Halfedge;
    typedef Halfedge::Base         HBase;
    typedef HDS::Face_handle       Face_handle;
    typedef CGAL::HalfedgeDS_decorator<HDS> Decorator;

    HDS hds(1,2,2);
    Decorator D(hds);

    D.create_loop();
    CGAL_assertion( hds.size_of_vertices() == 1);
    CGAL_assertion( hds.size_of_halfedges() == 2);
    CGAL_assertion( hds.size_of_faces() == 2);
    CGAL_assertion( D.is_valid( false, 3));
    (hds.halfedges_begin()+1)->HBase::set_face( Face_handle());
    hds.normalize_border();
    CGAL_assertion( hds.size_of_border_halfedges() == 1);
    CGAL_assertion( hds.size_of_border_edges() == 1);
    (hds.halfedges_begin()+1)->HBase::set_face( hds.faces_begin() + 1);
    hds.normalize_border();
    CGAL_assertion( hds.size_of_border_halfedges() == 0);
    CGAL_assertion( hds.size_of_border_edges() == 0);
    CGAL_assertion( D.is_valid( false, 4));
    hds.reserve(2,4,4);
    CGAL_assertion( hds.capacity_of_faces() == 4);
    CGAL_assertion( D.is_valid( false, 4));

    HDS hds2(hds); // copy constructor.
    Decorator D2(hds2);
    CGAL_assertion( D2.is_valid( false, 4));

    hds = hds2; // assignment.
    CGAL_assertion( D.is_valid( false, 4));
}

void test_HalfedgeDS_using_vector_min() {
    // Instantiation of the halfedge data structure using vector
    // and minimal bases as for an undirected graph.
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
    typedef CGAL::HalfedgeDS_using_vector< Empty_traits,
                                           CGAL::HalfedgeDS_min_items> HDS;
#else
    typedef CGAL::HalfedgeDS_using_vector::HDS< Empty_traits,
                                           CGAL::HalfedgeDS_min_items> HDS;
#endif
    typedef HDS::Halfedge   Halfedge;
    typedef Halfedge::Base  HBase;

    HDS hds(1,2,1);
    hds.edges_push_back( Halfedge(), Halfedge());
    hds.halfedges_begin()->HBase::set_next(hds.halfedges_begin()+1);
    (hds.halfedges_begin()+1)->HBase::set_next(hds.halfedges_begin());
    CGAL_assertion( hds.size_of_vertices() == 0);
    CGAL_assertion( hds.size_of_halfedges() == 2);
    CGAL_assertion( hds.size_of_faces() == 0);
    hds.normalize_border();
    CGAL_assertion( hds.size_of_border_halfedges() == 0);
    CGAL_assertion( hds.size_of_border_edges() == 0);
}


int main() {
    test_HalfedgeDS_polyhedron_default();
    test_HalfedgeDS_default();
    test_HalfedgeDS_using_list();
    test_HalfedgeDS_using_list_min();
    test_HalfedgeDS_using_vector();
    test_HalfedgeDS_using_vector_min();
    return 0;
}
// EOF //
