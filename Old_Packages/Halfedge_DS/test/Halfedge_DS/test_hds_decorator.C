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
// file          : test_hds_decorator.C
// chapter       : $CGAL_Chapter: Halfedge Data Structures $
// package       : $CGAL_Package: Halfedge_DS 2.8 (13 Sep 2000) $
// source        : hds_decorator.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : MPI Saarbruecken (Stefan Schirra <stschirr@mpi-sb.mpg.de>)
//
// Test Halfedge Data Structure Decorator.
// ============================================================================


#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif
#ifndef CGAL_CARTESIAN_H
#include <CGAL/Cartesian.h>
#endif // CGAL_CARTESIAN_H
#ifndef CGAL_POINT_2_H
#include <CGAL/Point_2.h>
#endif // CGAL_POINT_2_H
#ifndef CGAL_POINT_3_H
#include <CGAL/Point_3.h>
#endif // CGAL_POINT_3_H
#ifndef CGAL_PLANE_3_H
#include <CGAL/Plane_3.h>
#endif // CGAL_PLANE_3_H
#ifndef CGAL_HALFEDGE_DATA_STRUCTURE_DEFAULT_H
#include <CGAL/Halfedge_data_structure_default.h>
#endif // CGAL_HALFEDGE_DATA_STRUCTURE_DEFAULT_H
#ifndef CGAL_HALFEDGE_DATA_STRUCTURE_DECORATOR_H
#include <CGAL/Halfedge_data_structure_decorator.h>
#endif // CGAL_HALFEDGE_DATA_STRUCTURE_DECORATOR_H
#ifndef CGAL_HALFEDGE_DATA_STRUCTURE_USING_VECTOR_H
#include <CGAL/Halfedge_data_structure_using_vector.h>
#endif

// MSVC fixes:
CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(CGAL::Halfedge_max_base)
CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(CGAL::Halfedge_max_base*)
CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(CGAL::Vertex_min_base)
CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(CGAL::Halfedge_min_base)
CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(CGAL::Halfedge_min_base*)
CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(CGAL::Facet_min_base)
CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(CGAL::Facet_min_base*)

typedef CGAL::Vertex_max_base<CGAL::Point_3<CGAL::Cartesian<double> > > A0;
typedef CGAL::Polyhedron_facet_base_3<CGAL::Cartesian<double> > A1;
typedef CGAL::_HDS_vector_halfedge<A0, CGAL::Halfedge_max_base, A1> A2;
typedef CGAL::_HDS_vector_vertex<A0, CGAL::Halfedge_max_base, A1> A2v;
typedef CGAL::_HDS_vector_facet<A0, CGAL::Halfedge_max_base, A1> A2f;

CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(A0)
CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(A1)
CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(A2)
CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(A2v)
CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(A2f)

typedef CGAL::Facet_min_base A1m;
typedef CGAL::_HDS_vector_halfedge<
    CGAL::Vertex_min_base,CGAL::Halfedge_min_base, A1m> A2m;
typedef CGAL::_HDS_vector_vertex<
    CGAL::Vertex_min_base,CGAL::Halfedge_min_base, A1m> A2mv;
typedef CGAL::_HDS_vector_facet<
    CGAL::Vertex_min_base,CGAL::Halfedge_min_base, A1m> A2mf;

CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(A2m)
CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(A2m*)
CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(A2mv)
CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(A2mv*)
CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(A2mf)
CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC(A2mf*)

using namespace CGAL;

void test_Halfedge_data_structure_decorator() {
    // Simple instantiation of the default halfedge data structure.
    typedef Cartesian<double>                      Rep;
    typedef Point_2<Rep>                           Point;
    typedef Halfedge_data_structure_default<Point> HDS;
    typedef Halfedge_data_structure_decorator<HDS> Decorator;
    typedef HDS::Vertex                            Vertex;
    typedef HDS::Halfedge                          Halfedge;
    typedef HDS::Facet                             Facet;
    HDS hds;
    Decorator  decorator;
    // Check create single loop.
    Halfedge* h = decorator.create_loop(hds);
    (void)h; // suppress compiler warning
    hds.normalize_border();
    CGAL_assertion( hds.size_of_vertices() == 1);
    CGAL_assertion( hds.size_of_halfedges() == 2);
    CGAL_assertion( hds.size_of_facets() == 2);
    CGAL_assertion( decorator.is_valid( hds, false, 4));

    // Restart with open segment.
    hds.delete_all();
    hds.normalize_border();
    CGAL_assertion( decorator.is_valid( hds, false, 4));
    h = decorator.create_segment(hds);
    CGAL_assertion( hds.size_of_vertices() == 2);
    CGAL_assertion( hds.size_of_halfedges() == 2);
    CGAL_assertion( hds.size_of_facets() == 1);
    CGAL_assertion( decorator.is_valid( hds, false, 4));

    // Create border edge and check normalization.
    decorator.set_facet( h->opposite(), NULL);
    hds.normalize_border();
    CGAL_assertion( hds.size_of_border_halfedges() == 1);
    CGAL_assertion( hds.size_of_border_edges() == 1);
    CGAL_assertion( decorator.normalized_border_is_valid( hds));
    decorator.set_facet( h->opposite(), h->facet());
    hds.normalize_border();
    CGAL_assertion( hds.size_of_border_halfedges() == 0);
    CGAL_assertion( hds.size_of_border_edges() == 0);
    CGAL_assertion( decorator.is_valid( hds, false, 4));

    // Extend edge to two triangles.
    Halfedge* g = decorator.split_vertex( hds, h, h);
    CGAL_assertion( decorator.is_valid( hds, false, 4));
    CGAL_assertion( h != g);
    CGAL_assertion( h->next()->next() == g);
    CGAL_assertion( h == g->next()->next());
    CGAL_assertion( h->opposite() == g->next());
    CGAL_assertion( g->opposite() == h->next());
    Halfedge* g2 = decorator.split_facet( hds, h->opposite(),
                                          g->opposite());
    CGAL_assertion( decorator.is_valid( hds, false, 4));
    CGAL_assertion( h->opposite()->next() == g2);
    CGAL_assertion( g2->next() == g);
    decorator.split_vertex( hds, g2, g->opposite());
    CGAL_assertion( decorator.is_valid( hds, false, 4));
    CGAL_assertion( g->next()->next()->next()->next() == g);
    CGAL_assertion_code( Halfedge* g3 = 
        decorator.split_facet( hds, g2->next()->opposite(), h); )
    CGAL_assertion( decorator.is_valid( hds, false, 4));
    CGAL_assertion( g->next()->next()->next()->next() == g);
    CGAL_assertion( h->next()->next()->next() == h);
    CGAL_assertion( g3->next()->next()->next() == g3);
    CGAL_assertion( g3->next() == g->opposite());
    CGAL_assertion( g3->opposite()->next() == g2->opposite());
    CGAL_assertion( g3->opposite() == h->next());

    // Edge flip within the triangle.
    CGAL_assertion_code( Halfedge* g4 = decorator.flip_edge( g3);)
    CGAL_assertion( decorator.is_valid( hds, false, 4));
    CGAL_assertion( g4 == g3);
    CGAL_assertion( g3->next()->next() == g2->opposite());
    CGAL_assertion( g3->opposite()->next() == h);
    CGAL_assertion( g->next()->next()->next()->next() == g);
    CGAL_assertion( h->next()->next()->next() == h);
    CGAL_assertion( g3->next()->next()->next() == g3);

    // Reverse facet orientation.
    decorator.inside_out(hds);
    CGAL_assertion( decorator.is_valid( hds, false, 4));
    decorator.inside_out(hds);
    CGAL_assertion( decorator.is_valid( hds, false, 4));

    // Check hole manipulations.
    decorator.make_hole( hds, g);
    hds.normalize_border();
    CGAL_assertion( hds.size_of_border_halfedges() == 4);
    CGAL_assertion( hds.size_of_border_edges() == 4);
    CGAL_assertion( decorator.is_valid( hds, false, 4));

    // Reverse facet orientation, deal also with the hole..
    decorator.inside_out(hds);
    CGAL_assertion( decorator.is_valid( hds, false, 3));
    hds.normalize_border();
    CGAL_assertion( decorator.is_valid( hds, false, 4));

    // Check add_facet_to_border.
    hds.delete_all();
    h = decorator.create_loop(hds);
    decorator.make_hole( hds, h->opposite());
    hds.normalize_border();
    CGAL_assertion( decorator.is_valid( hds, false, 4));
    decorator.add_facet_to_border( hds, h->opposite(), h->opposite());
    CGAL_assertion( hds.size_of_halfedges() == 4);
    CGAL_assertion( hds.size_of_facets() == 2);
    CGAL_assertion( decorator.is_valid( hds, false, 3));
}

void test_Halfedge_data_structure_decorator2() {
    // Instantiation of the halfedge data structure using vector
    // with max-bases for a polyhedral surface.
    typedef Cartesian<double>                       Rep;
    typedef Point_3<Rep>                            Point;
    typedef Halfedge_data_structure_using_vector<
                Vertex_max_base<Point>,
                Halfedge_max_base,
                Polyhedron_facet_base_3<Rep> >      HDS;
    typedef Halfedge_data_structure_decorator<HDS>  Decorator;
    typedef HDS::Vertex                             Vertex;
    typedef HDS::Halfedge                           Halfedge;
    typedef HDS::Facet                              Facet;
    HDS hds(4,10,3);
    Decorator  decorator;
    // Check create single loop.
    Halfedge* h = decorator.create_loop(hds);
    (void)h; // suppress compiler warning
    hds.normalize_border();
    CGAL_assertion( hds.size_of_vertices() == 1);
    CGAL_assertion( hds.size_of_halfedges() == 2);
    CGAL_assertion( hds.size_of_facets() == 2);
    CGAL_assertion( decorator.is_valid( hds, false, 4));

    // Restart with open segment.
    hds.delete_all();
    hds.normalize_border();
    CGAL_assertion( decorator.is_valid( hds, false, 4));
    h = decorator.create_segment(hds);
    CGAL_assertion( hds.size_of_vertices() == 2);
    CGAL_assertion( hds.size_of_halfedges() == 2);
    CGAL_assertion( hds.size_of_facets() == 1);
    CGAL_assertion( decorator.is_valid( hds, false, 3));
    hds.normalize_border();
    CGAL_assertion( decorator.is_valid( hds, false, 4));

    // Create border edge and check normalization.
    decorator.set_facet( h->opposite(), NULL);
    hds.normalize_border();
    CGAL_assertion( hds.size_of_border_halfedges() == 1);
    CGAL_assertion( hds.size_of_border_edges() == 1);
    CGAL_assertion( decorator.normalized_border_is_valid( hds));
    decorator.set_facet( h->opposite(),  h->facet());
    hds.normalize_border();
    CGAL_assertion( hds.size_of_border_halfedges() == 0);
    CGAL_assertion( hds.size_of_border_edges() == 0);
    CGAL_assertion( decorator.is_valid( hds, false, 4));

    // Extend edge to two triangles.
    Halfedge* g = decorator.split_vertex( hds, h, h);
    CGAL_assertion( decorator.is_valid( hds, false, 3));
    hds.normalize_border();
    CGAL_assertion( decorator.is_valid( hds, false, 4));
    CGAL_assertion( h != g);
    CGAL_assertion( h->next()->next() == g);
    CGAL_assertion( h == g->next()->next());
    CGAL_assertion( h->opposite() == g->next());
    CGAL_assertion( g->opposite() == h->next());
    Halfedge* g2 = decorator.split_facet( hds, h->opposite(),
                                          g->opposite());
    hds.normalize_border();
    CGAL_assertion( decorator.is_valid( hds, false, 4));
    CGAL_assertion( h->opposite()->next() == g2);
    CGAL_assertion( g2->next() == g);
    decorator.split_vertex( hds, g2, g->opposite());
    CGAL_assertion( decorator.is_valid( hds, false, 3));
    hds.normalize_border();
    CGAL_assertion( decorator.is_valid( hds, false, 4));
    CGAL_assertion( g->next()->next()->next()->next() == g);
    CGAL_assertion_code( Halfedge* g3 = 
        decorator.split_facet( hds, g2->next()->opposite(), h);)
    hds.normalize_border();
    CGAL_assertion( decorator.is_valid( hds, false, 4));
    CGAL_assertion( g->next()->next()->next()->next() == g);
    CGAL_assertion( h->next()->next()->next() == h);
    CGAL_assertion( g3->next()->next()->next() == g3);
    CGAL_assertion( g3->next() == g->opposite());
    CGAL_assertion( g3->opposite()->next() == g2->opposite());
    CGAL_assertion( g3->opposite() == h->next());

    // Edge flip within the triangle.
    CGAL_assertion_code( Halfedge* g4 = decorator.flip_edge( g3);)
    CGAL_assertion( decorator.is_valid( hds, false, 3));
    hds.normalize_border();
    CGAL_assertion( decorator.is_valid( hds, false, 4));
    CGAL_assertion( g4 == g3);
    CGAL_assertion( g3->next()->next() == g2->opposite());
    CGAL_assertion( g3->opposite()->next() == h);
    CGAL_assertion( g->next()->next()->next()->next() == g);
    CGAL_assertion( h->next()->next()->next() == h);
    CGAL_assertion( g3->next()->next()->next() == g3);

    // Reverse facet orientation.
    decorator.inside_out(hds);
    CGAL_assertion( decorator.is_valid( hds, false, 4));
}

int main() {
    test_Halfedge_data_structure_decorator();
    test_Halfedge_data_structure_decorator2();
    return 0;
}
// EOF //
