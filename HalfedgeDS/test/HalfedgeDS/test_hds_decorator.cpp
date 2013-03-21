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
// file          : test/HalfedgeDS/test_hds_decorator.C
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
// Test Halfedge Data Structure Decorator.
// ============================================================================


#include <CGAL/Cartesian.h>
#include <CGAL/HalfedgeDS_default.h>
#include <CGAL/Polyhedron_items_3.h>
#include <CGAL/HalfedgeDS_decorator.h>

typedef CGAL::Cartesian<double> Kernel;
struct Dummy_traits_2 {
    typedef Kernel::Point_2 Point_2;
};

void test_HalfedgeDS_decorator() {
    // Simple instantiation of the default halfedge data structure.
    typedef CGAL_HALFEDGEDS_DEFAULT<Dummy_traits_2>  HDS;
    typedef CGAL::HalfedgeDS_decorator<HDS>          Decorator;
    typedef HDS::Halfedge_handle                     Halfedge_handle;
    typedef HDS::Face_handle                         Face_handle;
    HDS hds;
    Decorator  decorator(hds);
    // Check create single loop.
    Halfedge_handle h = decorator.create_loop();
    hds.normalize_border();
    CGAL_assertion( hds.size_of_vertices() == 1);
    CGAL_assertion( hds.size_of_halfedges() == 2);
    CGAL_assertion( hds.size_of_faces() == 2);
    CGAL_assertion( decorator.is_valid( false, 4));

    // Restart with open segment.
    hds.clear();
    hds.normalize_border();
    CGAL_assertion( decorator.is_valid( false, 4));
    h = decorator.create_segment();
    CGAL_assertion( hds.size_of_vertices() == 2);
    CGAL_assertion( hds.size_of_halfedges() == 2);
    CGAL_assertion( hds.size_of_faces() == 1);
    CGAL_assertion( decorator.is_valid( false, 4));

    // Create border edge and check normalization.
    decorator.set_face( h->opposite(), Face_handle());
    hds.normalize_border();
    CGAL_assertion( hds.size_of_border_halfedges() == 1);
    CGAL_assertion( hds.size_of_border_edges() == 1);
    CGAL_assertion( decorator.normalized_border_is_valid());
    decorator.set_face( h->opposite(), h->face());
    hds.normalize_border();
    CGAL_assertion( hds.size_of_border_halfedges() == 0);
    CGAL_assertion( hds.size_of_border_edges() == 0);
    CGAL_assertion( decorator.is_valid( false, 4));

    // Extend edge to two triangles.
    Halfedge_handle g = decorator.split_vertex( h, h);
    CGAL_assertion( decorator.is_valid( false, 4));
    CGAL_assertion( h != g);
    CGAL_assertion( h->next()->next() == g);
    CGAL_assertion( h == g->next()->next());
    CGAL_assertion( h->opposite() == g->next());
    CGAL_assertion( g->opposite() == h->next());
    Halfedge_handle g2 = decorator.split_face(h->opposite(),g->opposite());
    CGAL_assertion( decorator.is_valid( false, 4));
    CGAL_assertion( h->opposite()->next() == g2);
    CGAL_assertion( g2->next() == g);
    decorator.split_vertex( g2, g->opposite());
    CGAL_assertion( decorator.is_valid( false, 4));
    CGAL_assertion( g->next()->next()->next()->next() == g);
    Halfedge_handle g3 = 
        decorator.split_face( g2->next()->opposite(), h);
    CGAL_assertion( decorator.is_valid( false, 4));
    CGAL_assertion( g->next()->next()->next()->next() == g);
    CGAL_assertion( h->next()->next()->next() == h);
    CGAL_assertion( g3->next()->next()->next() == g3);
    CGAL_assertion( g3->next() == g->opposite());
    CGAL_assertion( g3->opposite()->next() == g2->opposite());
    CGAL_assertion( g3->opposite() == h->next());

    // Edge flip within the triangle.
    Halfedge_handle g4 = decorator.flip_edge( g3);
    CGAL_assertion( decorator.is_valid( false, 4));
    CGAL_assertion( g4 == g3);
    CGAL_assertion( g3->next()->next() == g2->opposite());
    CGAL_assertion( g3->opposite()->next() == h);
    CGAL_assertion( g->next()->next()->next()->next() == g);
    CGAL_assertion( h->next()->next()->next() == h);
    CGAL_assertion( g3->next()->next()->next() == g3);

    // Reverse face orientation.
    decorator.inside_out();
    CGAL_assertion( decorator.is_valid( false, 4));
    decorator.inside_out();
    CGAL_assertion( decorator.is_valid( false, 4));

    // Check hole manipulations.
    decorator.make_hole(g);
    hds.normalize_border();
    CGAL_assertion( hds.size_of_border_halfedges() == 4);
    CGAL_assertion( hds.size_of_border_edges() == 4);
    CGAL_assertion( decorator.is_valid( false, 4));

    // Reverse face orientation, deal also with the hole..
    decorator.inside_out();
    CGAL_assertion( decorator.is_valid( false, 3));
    hds.normalize_border();
    CGAL_assertion( decorator.is_valid( false, 4));

    // Check add_face_to_border.
    hds.clear();
    h = decorator.create_loop();
    decorator.make_hole( h->opposite());
    hds.normalize_border();
    CGAL_assertion( decorator.is_valid( false, 4));
    decorator.add_face_to_border( h->opposite(), h->opposite());
    CGAL_assertion( hds.size_of_halfedges() == 4);
    CGAL_assertion( hds.size_of_faces() == 2);
    CGAL_assertion( decorator.is_valid( false, 3));
}

#include <CGAL/HalfedgeDS_vector.h>

struct Dummy_traits_3 {
    typedef Kernel::Point_3    Point_3;
    typedef Kernel::Plane_3    Plane_3;
};

void test_HalfedgeDS_decorator2() {
    // Instantiation of the halfedge data structure using vector
    // with max-bases for a polyhedral surface.
    typedef CGAL::HalfedgeDS_vector< Dummy_traits_3,
                                           CGAL::Polyhedron_items_3> HDS;
    typedef CGAL::HalfedgeDS_decorator<HDS> Decorator;
    typedef HDS::Halfedge_handle            Halfedge_handle;
    typedef HDS::Face_handle                Face_handle;
    HDS hds(4,10,3);
    Decorator  decorator(hds);
    // Check create single loop.
    Halfedge_handle h = decorator.create_loop();
    hds.normalize_border();
    CGAL_assertion( hds.size_of_vertices() == 1);
    CGAL_assertion( hds.size_of_halfedges() == 2);
    CGAL_assertion( hds.size_of_faces() == 2);
    CGAL_assertion( decorator.is_valid( false, 4));

    // Restart with open segment.
    hds.clear();
    hds.normalize_border();
    CGAL_assertion( decorator.is_valid( false, 4));
    h = decorator.create_segment();
    CGAL_assertion( hds.size_of_vertices() == 2);
    CGAL_assertion( hds.size_of_halfedges() == 2);
    CGAL_assertion( hds.size_of_faces() == 1);
    CGAL_assertion( decorator.is_valid( false, 3));
    hds.normalize_border();
    CGAL_assertion( decorator.is_valid( false, 4));

    // Create border edge and check normalization.
    decorator.set_face( h->opposite(), Face_handle());
    hds.normalize_border();
    CGAL_assertion( hds.size_of_border_halfedges() == 1);
    CGAL_assertion( hds.size_of_border_edges() == 1);
    CGAL_assertion( decorator.normalized_border_is_valid());
    decorator.set_face( h->opposite(), h->face());
    hds.normalize_border();
    CGAL_assertion( hds.size_of_border_halfedges() == 0);
    CGAL_assertion( hds.size_of_border_edges() == 0);
    CGAL_assertion( decorator.is_valid( false, 4));

    // Extend edge to two triangles.
    Halfedge_handle g = decorator.split_vertex( h, h);
    CGAL_assertion( decorator.is_valid( false, 3));
    hds.normalize_border();
    CGAL_assertion( decorator.is_valid( false, 4));
    CGAL_assertion( h != g);
    CGAL_assertion( h->next()->next() == g);
    CGAL_assertion( h == g->next()->next());
    CGAL_assertion( h->opposite() == g->next());
    CGAL_assertion( g->opposite() == h->next());
    Halfedge_handle g2 = decorator.split_face(h->opposite(),g->opposite());
    hds.normalize_border();
    CGAL_assertion( decorator.is_valid( false, 4));
    CGAL_assertion( h->opposite()->next() == g2);
    CGAL_assertion( g2->next() == g);
    decorator.split_vertex( g2, g->opposite());
    CGAL_assertion( decorator.is_valid( false, 3));
    hds.normalize_border();
    CGAL_assertion( decorator.is_valid( false, 4));
    CGAL_assertion( g->next()->next()->next()->next() == g);
    Halfedge_handle g3 = 
        decorator.split_face( g2->next()->opposite(), h);
    hds.normalize_border();
    CGAL_assertion( decorator.is_valid( false, 4));
    CGAL_assertion( g->next()->next()->next()->next() == g);
    CGAL_assertion( h->next()->next()->next() == h);
    CGAL_assertion( g3->next()->next()->next() == g3);
    CGAL_assertion( g3->next() == g->opposite());
    CGAL_assertion( g3->opposite()->next() == g2->opposite());
    CGAL_assertion( g3->opposite() == h->next());

    // Edge flip within the triangle.
    Halfedge_handle g4 = decorator.flip_edge( g3);
    CGAL_assertion( decorator.is_valid( false, 3));
    hds.normalize_border();
    CGAL_assertion( decorator.is_valid( false, 4));
    CGAL_assertion( g4 == g3);
    CGAL_assertion( g3->next()->next() == g2->opposite());
    CGAL_assertion( g3->opposite()->next() == h);
    CGAL_assertion( g->next()->next()->next()->next() == g);
    CGAL_assertion( h->next()->next()->next() == h);
    CGAL_assertion( g3->next()->next()->next() == g3);

    // Reverse face orientation.
    decorator.inside_out();
    CGAL_assertion( decorator.is_valid( false, 4));
}


#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Modifier_base.h>

template <class HDS>
struct Build_triangle : public CGAL::Modifier_base<HDS> {
  double x_,y_; //origin

  Build_triangle(double x, double y):x_(x),y_(y) {}
  void operator()( HDS& hds) {
      // Postcondition: `hds' is a valid polyhedral surface.
      CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
      B.begin_surface( 4, 2, 8);
      typedef typename HDS::Vertex   Vertex;
      typedef typename Vertex::Point Point;
      B.add_vertex( Point( x_, y_));
      B.add_vertex( Point( x_+1, y_+0));
      B.add_vertex( Point( x_+0, y_+1));
      B.add_vertex( Point( x_+0, y_+-1));
      B.begin_facet();
      B.add_vertex_to_facet( 0);
      B.add_vertex_to_facet( 1);
      B.add_vertex_to_facet( 2);
      B.end_facet();
      B.begin_facet();
      B.add_vertex_to_facet( 1);
      B.add_vertex_to_facet( 0);
      B.add_vertex_to_facet( 3);
      B.end_facet();
      B.end_surface();
  }
};

#include <CGAL/HalfedgeDS_vertex_base.h>
#include <CGAL/HalfedgeDS_halfedge_base.h>
#include <CGAL/HalfedgeDS_face_base.h>

void test_HalfedgeDS_decorator3() {
    // Simple instantiation of the default halfedge data structure.
    typedef CGAL::HalfedgeDS_default<Dummy_traits_2>  HDS;
    typedef CGAL::HalfedgeDS_decorator<HDS>          Decorator;

    CGAL_USE_TYPE(HDS::Halfedge_handle);
    CGAL_USE_TYPE(HDS::Face_handle);

    HDS hds;

    Build_triangle<HDS> modifier(0,0);
    modifier(hds);
    modifier.x_=33;
    modifier.y_=33;
    modifier(hds);


    Decorator decorator(hds);
    decorator.keep_largest_connected_components(1);

    hds.normalize_border();

    CGAL_assertion( hds.size_of_vertices() == 4);
    CGAL_assertion( hds.size_of_halfedges() == 10);
    CGAL_assertion( hds.size_of_faces() == 2);
    CGAL_assertion( decorator.is_valid( false, 4));

    decorator.vertices_erase(hds.vertices_begin());
    decorator.vertices_erase(hds.vertices_begin(), hds.vertices_end());

    decorator.faces_erase(hds.faces_begin());
    decorator.faces_erase(hds.faces_begin(), hds.faces_end());

    CGAL_assertion( hds.size_of_vertices() == 0);
    CGAL_assertion( hds.size_of_halfedges() == 10);
    CGAL_assertion( hds.size_of_faces() == 0);

}
int main() {
    test_HalfedgeDS_decorator();
    test_HalfedgeDS_decorator2();
    test_HalfedgeDS_decorator3();
    return 0;
}
// EOF //
