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
// file          : test_hds_2.C
// chapter       : $CGAL_Chapter: Halfedge Data Structures $
// package       : $CGAL_Package: Halfedge_DS 2.8 (13 Sep 2000) $
// source        : hds.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : MPI Saarbruecken (Stefan Schirra <stschirr@mpi-sb.mpg.de>)
//
// Test hds.
// ============================================================================


#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif // CGAL_BASIC_H
#ifndef CGAL_PROTECT_CSTDDEF
#include <cstddef>
#define CGAL_PROTECT_CSTDDEF
#endif
#ifndef CGAL_CIRCULATOR_H
#include <CGAL/circulator.h>
#endif // CGAL_CIRCULATOR_H
#ifndef CGAL_CIRCULATOR_IMPL_H
#include <CGAL/circulator_impl.h>
#endif // CGAL_CIRCULATOR_IMPL_H
#ifndef CGAL_PROTECT_VECTOR
#include <vector>
#define CGAL_PROTECT_VECTOR
#endif
#ifndef CGAL_PROTECT_LIST
#include <list>
#define CGAL_PROTECT_LIST
#endif
#ifndef CGAL_CARTESIAN_H
#include <CGAL/Cartesian.h>
#endif // CGAL_CARTESIAN_H
#ifndef CGAL_HALFEDGE_DATA_STRUCTURE_POLYHEDRON_DEFAULT_3_H
#include <CGAL/Halfedge_data_structure_polyhedron_default_3.h>
#endif
#ifndef CGAL_POINT_3_H
#include <CGAL/Point_3.h>
#endif // CGAL_POINT_3_H
#ifndef CGAL_PLANE_3_H
#include <CGAL/Plane_3.h>
#endif // CGAL_PLANE_3_H
#ifndef CGAL_HALFEDGE_DATA_STRUCTURE_DEFAULT_H
#include <CGAL/Halfedge_data_structure_default.h>
#endif
#ifndef CGAL_POINT_2_H
#include <CGAL/Point_2.h>
#endif // CGAL_POINT_2_H
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

void test_Halfedge_data_structure_default() {
    // Simple instantiation of the default halfedge data structure.
    typedef Cartesian<double>                      Rep;
    typedef Point_2<Rep>                           Point;
    typedef Halfedge_data_structure_default<Point> HDS;
    Halfedge_data_structure_decorator<HDS> D;
    HDS hds;
    D.create_loop(hds);
    CGAL_assertion( hds.size_of_vertices() == 1);
    CGAL_assertion( hds.size_of_halfedges() == 2);
    CGAL_assertion( hds.size_of_facets() == 2);
    CGAL_assertion( D.is_valid( hds, false, 3));
    D.make_hole( hds, hds.halfedges_begin()->opposite());
    hds.normalize_border();
    CGAL_assertion( hds.size_of_border_halfedges() == 1);
    CGAL_assertion( hds.size_of_border_edges() == 1);
    CGAL_assertion( D.is_valid( hds, false, 4));
    D.fill_hole( hds, hds.halfedges_begin()->opposite());
    hds.normalize_border();
    CGAL_assertion( hds.size_of_border_halfedges() == 0);
    CGAL_assertion( hds.size_of_border_edges() == 0);
    CGAL_assertion( D.is_valid( hds, false, 4));

    HDS hds2(hds); // copy constructor.
    CGAL_assertion( D.is_valid( hds2, false, 4));
    hds = hds2; // assignment.
    CGAL_assertion( D.is_valid( hds, false, 4));

    // g++ 2.8 and egcs 2.90 linker relocation error:
    // This here is quite similar code, but it works out nicely, no error.
    typedef HDS::Edge_iterator Edge_iterator;
    CGAL_assertion_code( Edge_iterator i = hds.edges_begin();)
    // This operator== causes NO error here.
    CGAL_assertion(!(i == hds.edges_end()));

    D.erase_connected_component( hds, D.create_loop(hds));
    CGAL_assertion( hds.size_of_vertices() == 1);
    CGAL_assertion( hds.size_of_halfedges() == 2);
    CGAL_assertion( hds.size_of_facets() == 2);
    CGAL_assertion( D.is_valid( hds, false, 4));
    D.erase_connected_component( hds, D.create_segment(hds));
    CGAL_assertion( hds.size_of_vertices() == 1);
    CGAL_assertion( hds.size_of_halfedges() == 2);
    CGAL_assertion( hds.size_of_facets() == 2);
    CGAL_assertion( D.is_valid( hds, false, 4));
}

int main() {
    test_Halfedge_data_structure_default();
    return 0;
}
// EOF //
