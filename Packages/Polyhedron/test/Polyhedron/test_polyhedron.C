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
// file          : test_polyhedron.C
// chapter       : $CGAL_Chapter: 3D-Polyhedral Surfaces $
// package       : $CGAL_Package: Polyhedron 2.9 (13 Sep 2000) $
// source        : polyhedron.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : MPI Saarbruecken (Stefan Schirra <stschirr@mpi-sb.mpg.de>)
//
// Test Polyhedral Surfaces.
// ============================================================================


#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>

// This is the test file for the new design. Skip new design test for 
// disabled compilers.
#ifndef CGAL_USE_POLYHEDRON_DESIGN_ONE

#include <CGAL/HalfedgeDS_vector.h>
#include <CGAL/HalfedgeDS_list.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_traits_3.h>
#include <CGAL/Polyhedron_traits_with_normals_3.h>
#include <CGAL/Iterator_project.h>
#include <CGAL/function_objects.h>

// A polyhedron modifier that creates a tetrahedron using the
// incremental builder.
template < class HDS >
class Build_tetrahedron : public CGAL::Modifier_base<HDS> {
public:
    Build_tetrahedron() {}
        // creates the modifier.
    void operator()( HDS& target);
        // builds a tetrahedron.
        // Postcondition: `target' is a valid polyhedral surface.
};

template < class HDS >
void
Build_tetrahedron<HDS>:: operator()( HDS& target) {
    CGAL::Polyhedron_incremental_builder_3<HDS> B( target, true);
    B.begin_surface( 4, 4, 12);
    // Point coordinates suitable for integer coordinates.
    typedef typename HDS::Vertex Vertex;
    typedef typename Vertex::Point Point;
    B.add_vertex( Point( 0, 0, 1));
    B.add_vertex( Point( 1, 1, 1));
    B.add_vertex( Point( 0, 1, 0));
    B.add_vertex( Point( 1, 0, 0));
    B.begin_facet();
    B.add_vertex_to_facet( 1);
    B.add_vertex_to_facet( 3);
    B.add_vertex_to_facet( 0);
    B.end_facet();
    B.begin_facet();
    B.add_vertex_to_facet( 2);
    B.add_vertex_to_facet( 1);
    B.add_vertex_to_facet( 0);
    B.end_facet();
    B.begin_facet();
    B.add_vertex_to_facet( 3);
    B.add_vertex_to_facet( 2);
    B.add_vertex_to_facet( 0);
    B.end_facet();
    B.begin_facet();
    B.add_vertex_to_facet( 2);
    B.add_vertex_to_facet( 3);
    B.add_vertex_to_facet( 1);
    B.end_facet();
    B.end_surface();
}


void test_Polyhedron() {
    typedef CGAL::Cartesian<double>                     Kernel;
    typedef CGAL::Cartesian<int>                        KernelI;
    typedef CGAL::Point_3<Kernel>                       Point;
    typedef CGAL::Plane_3<Kernel>                       Plane;
    typedef CGAL::Polyhedron_traits_3<Kernel>           Traits;
    typedef CGAL::Polyhedron_traits_3<KernelI>          TraitsI;
    typedef CGAL::Polyhedron_traits_with_normals_3<Kernel>  TraitsN;

    // Using traits for testing of requirements minimality
    typedef CGAL::Polyhedron_3<Traits>                  Polyhedron;
    typedef CGAL::Polyhedron_3<TraitsI>                 PolyhedronI;

    // Test use of kernel as well
    typedef CGAL::Polyhedron_3<Kernel>                  PolyhedronK;

    // Test the use of normal vectors instead of plane equations
    typedef CGAL::Polyhedron_3<TraitsN>                 PolyhedronN;

    typedef CGAL::Polyhedron_3<Traits,
                             CGAL::Polyhedron_items_3,
                             CGAL::HalfedgeDS_vector>   PolyhedronV;
    typedef CGAL::Polyhedron_3<Traits,
                             CGAL::Polyhedron_items_3,
                             CGAL::HalfedgeDS_list>     PolyhedronL;

    typedef Polyhedron::HDS                     HDS;
    typedef PolyhedronI::HDS                    HDSI;
    typedef PolyhedronK::HDS                    HDSK;
    typedef PolyhedronV::HDS                    HDSV;
    typedef PolyhedronL::HDS                    HDSL;

    typedef Polyhedron::Vertex                  Vertex;
    typedef Polyhedron::Halfedge                Halfedge;
    typedef Polyhedron::Facet                   Facet;

    typedef Polyhedron::Vertex_iterator         Vertex_iterator;
    typedef Polyhedron::Facet_iterator          Facet_iterator;
    typedef Polyhedron::Halfedge_handle         Halfedge_handle;
    typedef Polyhedron::Halfedge_iterator       Halfedge_iterator;
    typedef Polyhedron::Edge_iterator           Edge_iterator;
    typedef Polyhedron::Halfedge_const_handle   Halfedge_const_handle;
    typedef Polyhedron::Halfedge_const_iterator Halfedge_const_iterator;
    typedef Polyhedron::Edge_const_iterator     Edge_const_iterator;

    typedef Polyhedron::Halfedge_around_vertex_circulator
                                   Halfedge_around_vertex_circulator;
    typedef Polyhedron::Halfedge_around_vertex_const_circulator
                                   Halfedge_around_vertex_const_circulator;
    typedef Polyhedron::Halfedge_around_facet_circulator
                                   Halfedge_around_facet_circulator;
    typedef Polyhedron::Halfedge_around_facet_const_circulator
                                   Halfedge_around_facet_const_circulator;

    typedef PolyhedronK::Halfedge_handle         HalfedgeK_handle;

    typedef PolyhedronN::Halfedge_handle         HalfedgeN_handle;

    typedef PolyhedronV::Halfedge_handle         HalfedgeV_handle;
    typedef PolyhedronV::Halfedge_iterator       HalfedgeV_iterator;
    typedef PolyhedronV::Edge_iterator           EdgeV_iterator;
    typedef PolyhedronV::Halfedge_const_handle   HalfedgeV_const_handle;
    typedef PolyhedronV::Halfedge_const_iterator HalfedgeV_const_iterator;
    typedef PolyhedronV::Edge_const_iterator     EdgeV_const_iterator;

    typedef PolyhedronV::Halfedge_around_vertex_circulator
                                  HalfedgeV_around_vertex_circulator;
    typedef PolyhedronV::Halfedge_around_vertex_const_circulator
                                  HalfedgeV_around_vertex_const_circulator;
    typedef PolyhedronV::Halfedge_around_facet_circulator
                                  HalfedgeV_around_facet_circulator;
    typedef PolyhedronV::Halfedge_around_facet_const_circulator
                                  HalfedgeV_around_facet_const_circulator;

    typedef PolyhedronL::Halfedge_handle        HalfedgeL_handle;

    {
        // Check if all automatic conversions work as promised.
        // First the conversions to Halfedge_const_iterator.
        // Simultaneously the conversion of an Halfedge to its handle
        // with the handle() member function is checked.
        Polyhedron P;

        Halfedge_handle  hh = P.make_triangle();  // Check handle.
        CGAL_assertion_code( Halfedge_const_handle  hch(hh);)
        CGAL_assertion( P.is_triangle( hh));
        CGAL_assertion( P.is_triangle( hch));
        CGAL_assertion( P.is_triangle( HDS::halfedge_handle(&*hh)));
        CGAL_assertion( P.is_triangle( HDS::halfedge_handle(&*hch)));
        //CGAL_assertion( P.is_triangle( hh->handle()));
        //CGAL_assertion( P.is_triangle( hch->handle()));

        Halfedge_iterator  hi(hh);          // Check Iterator.
        CGAL_assertion( P.is_triangle( hi));
        CGAL_assertion_code( Halfedge_const_iterator  hci(hh);)
        CGAL_assertion( P.is_triangle( hci));

        CGAL_assertion_code( Edge_iterator ei(hh);) // Check Edge Iterator.
        CGAL_assertion( P.is_triangle( ei));
        CGAL_assertion_code( Edge_const_iterator eci(hh);)
        CGAL_assertion( P.is_triangle( eci));

        Halfedge_around_facet_circulator  hfc(hh);  // Check circulator 1.
        CGAL_assertion( P.is_triangle( hfc));
        CGAL_assertion_code( Halfedge_around_facet_const_circulator hfcc(hh);)
        CGAL_assertion( P.is_triangle( hfcc));

        Halfedge_around_vertex_circulator  hvc(hh);  // Check circulator 2.
        CGAL_assertion( P.is_triangle( hvc));
        CGAL_assertion_code( Halfedge_around_vertex_const_circulator hvcc(hh);)
        CGAL_assertion( P.is_triangle( hvcc));


        // Next the conversions to Halfedge_iterator.
        hh = hh->opposite();
        P.fill_hole(hh);
        P.make_hole(hi);
        P.fill_hole(hfc);
        CGAL_assertion( ! P.is_triangle( hh->opposite()));
        P.make_hole(hvc);
        CGAL_assertion( P.is_triangle( hh->opposite()));
    }
    {
        // Same conversion checks for polyhedron with vector.
        PolyhedronV P;

        // Check handle.
        CGAL_assertion_code( HalfedgeV_handle  hh = P.make_triangle();)
        CGAL_assertion_code( HalfedgeV_const_handle  hch(hh);)
        CGAL_assertion( P.is_triangle( hh));
        CGAL_assertion( P.is_triangle( hch));
        CGAL_assertion( P.is_triangle( HDSV::halfedge_handle(&*hh)));
        CGAL_assertion( P.is_triangle( HDSV::halfedge_handle(&*hch)));

        CGAL_assertion_code( HalfedgeV_iterator  hi(hh);)  // Check Iterator.
        CGAL_assertion( P.is_triangle( hi));
        CGAL_assertion_code( HalfedgeV_const_iterator  hci(hh);)
        CGAL_assertion( P.is_triangle( hci));

        CGAL_assertion_code( EdgeV_iterator  ei(hh);) // Check Edge Iterator.
        CGAL_assertion( P.is_triangle( ei));
        CGAL_assertion_code( EdgeV_const_iterator  eci(hh);)
        CGAL_assertion( P.is_triangle( eci));

        // Check circulator 1.
        CGAL_assertion_code( HalfedgeV_around_facet_circulator hfc(hh);)
        CGAL_assertion( P.is_triangle( hfc));
        CGAL_assertion_code( HalfedgeV_around_facet_const_circulator hfcc(hh);)
        CGAL_assertion( P.is_triangle( hfcc));

        // Check circulator 2.
        CGAL_assertion_code( HalfedgeV_around_vertex_circulator  hvc(hh);)
        CGAL_assertion( P.is_triangle( hvc));
        CGAL_assertion_code(HalfedgeV_around_vertex_const_circulator hvcc(hh);)
        CGAL_assertion( P.is_triangle( hvcc));
    }
    {
        // The first check that the polyhedron and its normalization works.
        Polyhedron P;
        Halfedge_handle h = P.make_triangle();
        CGAL_assertion( P.is_valid());
        CGAL_assertion( P.is_triangle( h));
        CGAL_assertion( ! P.is_tetrahedron( h));
        P.normalize_border();
        CGAL_assertion( P.is_valid( false, 1));
        P.inside_out();
        CGAL_assertion( P.is_valid( false, 1));
        CGAL_assertion( P.is_triangle( h));
        P.normalize_border();
        CGAL_assertion( P.is_valid( false, 1));
        CGAL_assertion( P.is_triangle( h));

        h = P.make_tetrahedron();
        CGAL_assertion( P.is_valid());
        CGAL_assertion( P.is_tetrahedron( h));
        CGAL_assertion( ! P.is_triangle( h));
        P.normalize_border();
        CGAL_assertion( P.is_valid( false, 1));
        P.inside_out();
        CGAL_assertion( P.is_valid( false, 1));
        CGAL_assertion( P.is_tetrahedron( h));
        P.normalize_border();
        CGAL_assertion( P.is_valid( false, 1));
        CGAL_assertion( P.is_tetrahedron( h));
        P.make_hole(h);
        CGAL_assertion( ! P.is_tetrahedron( h));
        P.fill_hole(h);
        CGAL_assertion( P.is_tetrahedron( h));

        Polyhedron P2;
        Build_tetrahedron<HDS> modifier;
        P2.delegate( modifier);
        CGAL_assertion( P2.is_tetrahedron(P2.halfedges_begin()));
        P2.normalize_border();
        CGAL_assertion( P2.is_valid( false, 1));

        Polyhedron P3(P2);
        CGAL_assertion( P3.is_tetrahedron(P3.halfedges_begin()));
        P3.inside_out();
        P3.inside_out();
        P2 = P3;
        CGAL_assertion( P2.is_tetrahedron(P2.halfedges_begin()));
        P2.inside_out();
        CGAL_assertion( P2.is_tetrahedron(P2.halfedges_begin()));
    }
    {
        // Check the predefined point iterator
        Polyhedron P;
        CGAL_assertion_code( Halfedge_handle h = P.make_tetrahedron();)
        CGAL_assertion( P.is_valid());
        CGAL_assertion( P.is_tetrahedron( h));
        typedef Polyhedron::Point_iterator Point_iterator;
        Point_iterator begin( P.points_begin());
        Point_iterator end( P.points_end());
        Vertex_iterator i = P.vertices_begin();
        while( begin != end) {
            CGAL_assertion( i->point() == *begin);
            ++begin;
            ++i;
        }
        CGAL_assertion( i == P.vertices_end());
    }
    {
        // Check the predefined point const_iterator
        Polyhedron P;
        CGAL_assertion_code( Halfedge_handle h = P.make_tetrahedron();)
        CGAL_assertion( P.is_valid());
        CGAL_assertion( P.is_tetrahedron( h));
        typedef Polyhedron::Point_const_iterator Point_const_iterator;
        const Polyhedron& P2(P);
        Point_const_iterator begin( P2.points_begin());
        Point_const_iterator end( P2.points_end());
        Vertex_iterator i = P.vertices_begin();
        while( begin != end) {
            CGAL_assertion( i->point() == *begin);
            ++begin;
            ++i;
        }
        CGAL_assertion( i == P.vertices_end());
    }
    {
        // Check the predefined plane iterator
        Polyhedron P;
        P.make_tetrahedron();
        typedef Polyhedron::Plane_iterator Plane_iterator;
        Plane_iterator begin( P.planes_begin());
        Plane_iterator end( P.planes_end());
        Facet_iterator i = P.facets_begin();
        while( begin != end) {
            CGAL_assertion( i->plane() == *begin);
            ++begin;
            ++i;
        }
        CGAL_assertion( i == P.facets_end());
    }
    {
        // Check the predefined plane const_iterator
        Polyhedron P;
        CGAL_assertion_code( Halfedge_handle h = P.make_tetrahedron();)
        CGAL_assertion( P.is_valid());
        CGAL_assertion( P.is_tetrahedron( h));
        typedef Polyhedron::Plane_const_iterator Plane_const_iterator;
        const Polyhedron& P2(P);
        Plane_const_iterator begin( P2.planes_begin());
        Plane_const_iterator end( P2.planes_end());
        Facet_iterator i = P.facets_begin();
        while( begin != end) {
            CGAL_assertion( i->plane() == *begin);
            ++begin;
            ++i;
        }
        CGAL_assertion( i == P.facets_end());
    }
    {
        // Check the easy creation of a point iterator.
        Polyhedron P;
        CGAL_assertion_code( Halfedge_handle h = P.make_tetrahedron();)
        CGAL_assertion( P.is_valid());
        CGAL_assertion( P.is_tetrahedron( h));

        typedef CGAL::Project_point<Vertex>        Proj_point;
        typedef CGAL::Iterator_project<Vertex_iterator, Proj_point>
                                                   Point_iterator;

        Point_iterator begin( P.vertices_begin());
        Point_iterator end( P.vertices_end());
        Vertex_iterator i = P.vertices_begin();
        while( begin != end) {
            CGAL_assertion( i->point() == *begin);
            ++begin;
            ++i;
        }
        CGAL_assertion( i == P.vertices_end());
    }
    {
        // Check border facet generation.
        Polyhedron P;
        CGAL_assertion_code( Halfedge_handle h = P.make_triangle();)
        CGAL_assertion( P.is_valid());
        CGAL_assertion( P.is_triangle( h));

        CGAL_assertion_code( Halfedge_handle g = 
                             P.add_vertex_and_facet_to_border(
                                h->next()->opposite(),
                                h->opposite());)
        CGAL_assertion( P.is_valid());
        CGAL_assertion( h->next()->next()->next() == h);
        CGAL_assertion( g->next()->next()->next() == g);
        CGAL_assertion( h->opposite()->next() == g);

        CGAL_assertion_code( Halfedge_handle gg = P.add_facet_to_border(
                                 h->next()->next()->opposite(),
                                 g->next()->opposite());)
        CGAL_assertion( P.is_valid());
        CGAL_assertion(  h->next()->next()->next() == h);
        CGAL_assertion(  g->next()->next()->next() == g);
        CGAL_assertion( gg->next()->next()->next() == gg);
        CGAL_assertion( h->opposite()->next() == g);
        CGAL_assertion( gg->next()->opposite() == h->next());
    }
    {
        // Check erasing of facets and connected components.
        Polyhedron P;
        Halfedge_handle h = P.make_tetrahedron();
        Halfedge_handle g = P.make_tetrahedron();
        CGAL_assertion( P.is_valid());
        CGAL_assertion( P.is_tetrahedron( h));
        CGAL_assertion( P.is_tetrahedron( g));
        P.erase_connected_component(h);
        CGAL_assertion( P.is_valid());
        CGAL_assertion( P.is_tetrahedron( g));
        CGAL_assertion( P.size_of_vertices()  == 4);
        CGAL_assertion( P.size_of_halfedges() == 12);
        CGAL_assertion( P.size_of_facets()    == 4);
        P.erase_connected_component( P.make_triangle());
        CGAL_assertion( P.is_valid());
        CGAL_assertion( P.is_tetrahedron( g));
        CGAL_assertion( P.size_of_vertices()  == 4);
        CGAL_assertion( P.size_of_halfedges() == 12);
        CGAL_assertion( P.size_of_facets()    == 4);
        P.erase_facet( g);
        CGAL_assertion( P.is_valid());
        CGAL_assertion( P.size_of_vertices()  == 4);
        CGAL_assertion( P.size_of_halfedges() == 12);
        CGAL_assertion( P.size_of_facets()    == 3);
        P.erase_facet( g->opposite());
        CGAL_assertion( P.is_valid());
        CGAL_assertion( P.size_of_vertices()  == 4);
        CGAL_assertion( P.size_of_halfedges() == 10);
        CGAL_assertion( P.size_of_facets()    == 2);
        P.clear();
        CGAL_assertion( P.is_valid());
        CGAL_assertion( P.size_of_vertices()  == 0);
        CGAL_assertion( P.size_of_halfedges() == 0);
        CGAL_assertion( P.size_of_facets()    == 0);
    }
    {
        // Check bug-fix in join_vertex with respect to border edges.
        Polyhedron P;
        Halfedge_handle h = P.make_triangle();
        P.split_vertex( h, h->next()->opposite());
        P.join_vertex( h);
        CGAL_assertion( P.is_valid());
    }
    {
        // Check set_halfedge() for vertices and facets
        Polyhedron P;
        Halfedge_handle h = P.make_tetrahedron();
        h->vertex()->set_halfedge(h);
        CGAL_assertion( h->vertex()->halfedge() == h);
        CGAL_assertion( h->vertex() == h->next()->opposite()->vertex());
        h->vertex()->set_halfedge(h->next()->opposite());
        CGAL_assertion( h->vertex()->halfedge() == h->next()->opposite());
        CGAL_assertion( P.is_valid());

        h->facet()->set_halfedge(h);
        CGAL_assertion( h->facet()->halfedge() == h);
        CGAL_assertion( h->facet() == h->next()->facet());
        h->facet()->set_halfedge(h->next());
        CGAL_assertion( h->facet()->halfedge() == h->next());
        CGAL_assertion( P.is_valid());
    }
    {
        // A first check that integer coordinates are supported.
        PolyhedronI P;
        Build_tetrahedron<HDSI> modifierI;
        P.delegate( modifierI);
        CGAL_assertion( P.is_tetrahedron(P.halfedges_begin()));
        P.normalize_border();
        CGAL_assertion( P.is_valid());

        PolyhedronI P2(P);
        CGAL_assertion( P2.is_tetrahedron(P2.halfedges_begin()));
    }
    {
        // Check normalization as in polyhedron_prog5.C.
        Point p( 0.0, 0.0, 0.0);
        Point q( 1.0, 0.0, 0.0);
        Point r( 0.0, 1.0, 0.0);
        Point s( 0.0, 0.0, 1.0);

        Polyhedron P;
        P.make_tetrahedron( p, q, r, s);
        /* Remove the first facet, disturbing the border edge order. */
        P.make_hole( P.halfedges_begin());
        /* Reastablish border edge order. */
        P.normalize_border();
        /* Check it out with an halfedge iterator. */
        Halfedge_iterator h = P.halfedges_begin();
        /* The first three edges must be non-border edges. */
        CGAL_assertion( ! h->is_border_edge());
        ++ ++h;
        CGAL_assertion( ! h->is_border_edge());
        ++ ++h;
        CGAL_assertion( ! h->is_border_edge());
        ++ ++h;
        /* Here the three border edges should start. */
        CGAL_assertion( h == P.border_halfedges_begin());
        CGAL_assertion( h->is_border_edge());
        CGAL_assertion( h->opposite()->is_border());
        ++ ++h;
        CGAL_assertion( h->is_border_edge());
        CGAL_assertion( h->opposite()->is_border());
        ++ ++h;
        CGAL_assertion( h->is_border_edge());
        CGAL_assertion( h->opposite()->is_border());
        ++ ++h;
        CGAL_assertion( h == P.halfedges_end());
    }
    {
        // Check invariants of Euler operations.
        Polyhedron P;
        Halfedge_handle h = P.make_tetrahedron();
        CGAL_assertion( P.is_tetrahedron( h));
        Halfedge_handle g = P.split_vertex( h, h->next()->opposite());
        CGAL_assertion( g == h->next()->opposite());
        CGAL_assertion( P.is_valid());
        Halfedge_handle i = P.join_vertex( g);
        CGAL_assertion( i == h);
        CGAL_assertion( P.is_valid());
        CGAL_assertion( P.is_tetrahedron( h));

        g = P.split_vertex( h, h->next()->opposite());
        CGAL_assertion( g == h->next()->opposite());
        CGAL_assertion( P.is_valid());
        // Create a diagonal in the quadrangle.
        i = P.split_facet( g, g->next()->next());
        CGAL_assertion( i == g->next());
        CGAL_assertion( P.is_valid());
        Halfedge_handle j = P.join_facet( i);
        CGAL_assertion( j == g);
        CGAL_assertion( P.is_valid());
        j = P.join_vertex( g);
        CGAL_assertion( j == h);
        CGAL_assertion( P.is_valid());
        CGAL_assertion( P.is_tetrahedron( h));
        g = P.make_tetrahedron();
        CGAL_assertion_code( Halfedge_handle gg = g->opposite()->next();)
        CGAL_assertion( P.size_of_vertices()  == 8);
        CGAL_assertion( P.size_of_halfedges() == 24);
        CGAL_assertion( P.size_of_facets()    == 8);
        CGAL_assertion( P.is_valid());
        i = P.join_loop( h, g);
        CGAL_assertion( i == h);
        CGAL_assertion( P.size_of_vertices()  == 5);
        CGAL_assertion( P.size_of_halfedges() == 18);
        CGAL_assertion( P.size_of_facets()    == 6);
        // The unchanged facets.
        CGAL_assertion( h->opposite()->facet()->halfedge()->facet() ==
                   h->opposite()->facet());
        CGAL_assertion( h->opposite()->next_on_vertex()->facet()->
                        halfedge()->facet() ==
                        h->opposite()->next_on_vertex()->facet());
        CGAL_assertion( h->opposite()->next_on_vertex()->next_on_vertex()->
                        facet()->halfedge()->facet() ==
                        h->opposite()->next_on_vertex()->
                        next_on_vertex()->facet());
        // The changed facets.
        CGAL_assertion( h->facet()->halfedge()->facet() == h->facet());
        CGAL_assertion( h->next_on_vertex()->facet()->halfedge()->facet()
                        == h->next_on_vertex()->facet());
        CGAL_assertion( h->next_on_vertex()->next_on_vertex()->facet()->
                   halfedge()->facet() ==
                   h->next_on_vertex()->next_on_vertex()->facet());
        CGAL_assertion( P.is_valid());
        i = h->next()->opposite()->next();
        j = i->next()->opposite()->next();
        i = P.split_loop( h, i, j);
        CGAL_assertion( i->opposite()->next() == gg);
        CGAL_assertion( P.size_of_vertices()  == 8);
        CGAL_assertion( P.size_of_halfedges() == 24);
        CGAL_assertion( P.size_of_facets()    == 8);
        CGAL_assertion( P.is_valid());
        CGAL_assertion( P.is_tetrahedron( h));
        CGAL_assertion( P.is_tetrahedron( i));

        // create and erase center_vertex
        i = P.create_center_vertex( h);
        CGAL_assertion( i == h->next());
        CGAL_assertion( P.is_valid());
        j = P.erase_center_vertex( i);
        CGAL_assertion( j == h);
        CGAL_assertion( P.is_valid());
        CGAL_assertion( P.is_tetrahedron( h));
    }
    {
        // The first check that the polyhedron and its normalization works.
        PolyhedronL P;
        HalfedgeL_handle h = P.make_triangle();
        CGAL_assertion( P.is_valid());
        CGAL_assertion( P.is_triangle( h));
        CGAL_assertion( ! P.is_tetrahedron( h));
        P.normalize_border();
        CGAL_assertion( P.is_valid( false, 1));
        P.inside_out();
        CGAL_assertion( P.is_valid( false, 1));
        CGAL_assertion( P.is_triangle( h));
        P.normalize_border();
        CGAL_assertion( P.is_valid( false, 1));
        CGAL_assertion( P.is_triangle( h));

        h = P.make_tetrahedron();
        CGAL_assertion( P.is_valid());
        CGAL_assertion( P.is_tetrahedron( h));
        CGAL_assertion( ! P.is_triangle( h));
        P.normalize_border();
        CGAL_assertion( P.is_valid( false, 1));
        P.inside_out();
        CGAL_assertion( P.is_valid( false, 1));
        CGAL_assertion( P.is_tetrahedron( h));
        P.normalize_border();
        CGAL_assertion( P.is_valid( false, 1));
        CGAL_assertion( P.is_tetrahedron( h));
        P.make_hole(h);
        CGAL_assertion( ! P.is_tetrahedron( h));
        P.fill_hole(h);
        CGAL_assertion( P.is_tetrahedron( h));

        PolyhedronL P2;
        Build_tetrahedron<HDSL> modifier;
        P2.delegate( modifier);
        CGAL_assertion( P2.is_tetrahedron(P2.halfedges_begin()));
        P2.normalize_border();
        CGAL_assertion( P2.is_valid( false, 1));

        PolyhedronL P3(P2);
        CGAL_assertion( P3.is_tetrahedron(P3.halfedges_begin()));
        P3.inside_out();
        P3.inside_out();
        P2 = P3;
        CGAL_assertion( P2.is_tetrahedron(P2.halfedges_begin()));
        P2.inside_out();
        CGAL_assertion( P2.is_tetrahedron(P2.halfedges_begin()));
    }
    {
        // Check invariants of Euler operations.
        PolyhedronL P;
        HalfedgeL_handle h = P.make_tetrahedron();
        CGAL_assertion( P.is_tetrahedron( h));
        HalfedgeL_handle g = P.split_vertex( h, h->next()->opposite());
        CGAL_assertion( g == h->next()->opposite());
        CGAL_assertion( P.is_valid());
        HalfedgeL_handle i = P.join_vertex( g);
        CGAL_assertion( i == h);
        CGAL_assertion( P.is_valid());
        CGAL_assertion( P.is_tetrahedron( h));

        g = P.split_vertex( h, h->next()->opposite());
        CGAL_assertion( g == h->next()->opposite());
        CGAL_assertion( P.is_valid());
        // Create a diagonal in the quadrangle.
        i = P.split_facet( g, g->next()->next());
        CGAL_assertion( i == g->next());
        CGAL_assertion( P.is_valid());
        HalfedgeL_handle j = P.join_facet( i);
        CGAL_assertion( j == g);
        CGAL_assertion( P.is_valid());
        j = P.join_vertex( g);
        CGAL_assertion( j == h);
        CGAL_assertion( P.is_valid());
        CGAL_assertion( P.is_tetrahedron( h));
        g = P.make_tetrahedron();
        CGAL_assertion_code( HalfedgeL_handle gg = g->opposite()->next();)
        CGAL_assertion( P.size_of_vertices()  == 8);
        CGAL_assertion( P.size_of_halfedges() == 24);
        CGAL_assertion( P.size_of_facets()    == 8);
        CGAL_assertion( P.is_valid());
        i = P.join_loop( h, g);
        CGAL_assertion( i == h);
        CGAL_assertion( P.size_of_vertices()  == 5);
        CGAL_assertion( P.size_of_halfedges() == 18);
        CGAL_assertion( P.size_of_facets()    == 6);
        // The unchanged facets.
        CGAL_assertion( h->opposite()->facet()->halfedge()->facet() ==
                   h->opposite()->facet());
        CGAL_assertion( h->opposite()->next_on_vertex()->facet()->
                        halfedge()->facet() ==
                        h->opposite()->next_on_vertex()->facet());
        CGAL_assertion( h->opposite()->next_on_vertex()->next_on_vertex()->
                        facet()->halfedge()->facet() ==
                        h->opposite()->next_on_vertex()->
                        next_on_vertex()->facet());
        // The changed facets.
        CGAL_assertion( h->facet()->halfedge()->facet() == h->facet());
        CGAL_assertion( h->next_on_vertex()->facet()->halfedge()->facet()
                        == h->next_on_vertex()->facet());
        CGAL_assertion( h->next_on_vertex()->next_on_vertex()->facet()->
                   halfedge()->facet() ==
                   h->next_on_vertex()->next_on_vertex()->facet());
        CGAL_assertion( P.is_valid());
        i = h->next()->opposite()->next();
        j = i->next()->opposite()->next();
        i = P.split_loop( h, i, j);
        CGAL_assertion( i->opposite()->next() == gg);
        CGAL_assertion( P.size_of_vertices()  == 8);
        CGAL_assertion( P.size_of_halfedges() == 24);
        CGAL_assertion( P.size_of_facets()    == 8);
        CGAL_assertion( P.is_valid());
        CGAL_assertion( P.is_tetrahedron( h));
        CGAL_assertion( P.is_tetrahedron( i));
    }
    {
        // The first check that the polyhedron and its normalization works.
        PolyhedronV P(7,18,5);  // 1 tetra + 1 triangle.
        HalfedgeV_handle h = P.make_triangle();
        CGAL_assertion( P.is_valid());
        CGAL_assertion( P.is_triangle( h));
        CGAL_assertion( ! P.is_tetrahedron( h));
        P.normalize_border();
        CGAL_assertion( P.is_valid( false, 1));
        P.inside_out();
        CGAL_assertion( P.is_valid( false, 1));
        CGAL_assertion( P.is_triangle( h));
        P.normalize_border();
        CGAL_assertion( P.is_valid( false, 1));
        CGAL_assertion( P.is_triangle( h));

        h = P.make_tetrahedron();
        CGAL_assertion( P.is_valid());
        CGAL_assertion( P.is_tetrahedron( h));
        CGAL_assertion( ! P.is_triangle( h));
        P.normalize_border();
        CGAL_assertion( P.is_valid( false, 1));
        P.inside_out();
        CGAL_assertion( P.is_valid( false, 1));
        CGAL_assertion( P.is_tetrahedron( h));
        P.normalize_border();
        CGAL_assertion( P.is_valid( false, 1));
        CGAL_assertion( P.is_tetrahedron( h));

        PolyhedronV P2(4,12,4);  // 1 tetra.
        Build_tetrahedron<HDSV> modifier;
        P2.delegate( modifier);
        CGAL_assertion( P2.is_tetrahedron( P2.halfedges_begin()));
        P2.normalize_border();
        CGAL_assertion( P2.is_valid( false, 1));

        PolyhedronV P3(P2);
        CGAL_assertion( P3.is_tetrahedron( P3.halfedges_begin()));
        P3.inside_out();
        P3.inside_out();
        P2 = P3;
        CGAL_assertion( P2.is_tetrahedron( P2.halfedges_begin()));
        P2.inside_out();
        CGAL_assertion( P2.is_tetrahedron( P2.halfedges_begin()));
    }
    {
        // Check use of kernel as traits class.
        PolyhedronK P;
        HalfedgeK_handle h = P.make_triangle();
        CGAL_assertion( P.is_valid());
        CGAL_assertion( P.is_triangle( h));
        CGAL_assertion( ! P.is_tetrahedron( h));
        P.normalize_border();
        CGAL_assertion( P.is_valid( false, 1));
        P.inside_out();
        CGAL_assertion( P.is_valid( false, 1));
        CGAL_assertion( P.is_triangle( h));
        P.normalize_border();
        CGAL_assertion( P.is_valid( false, 1));
        CGAL_assertion( P.is_triangle( h));

        h = P.make_tetrahedron();
        CGAL_assertion( P.is_valid());
        CGAL_assertion( P.is_tetrahedron( h));
        CGAL_assertion( ! P.is_triangle( h));
        P.normalize_border();
        CGAL_assertion( P.is_valid( false, 1));
        P.inside_out();
        CGAL_assertion( P.is_valid( false, 1));
        CGAL_assertion( P.is_tetrahedron( h));
        P.normalize_border();
        CGAL_assertion( P.is_valid( false, 1));
        CGAL_assertion( P.is_tetrahedron( h));
        P.make_hole(h);
        CGAL_assertion( ! P.is_tetrahedron( h));
        P.fill_hole(h);
        CGAL_assertion( P.is_tetrahedron( h));

        Polyhedron P2;
        Build_tetrahedron<HDS> modifier;
        P2.delegate( modifier);
        CGAL_assertion( P2.is_tetrahedron(P2.halfedges_begin()));
        P2.normalize_border();
        CGAL_assertion( P2.is_valid( false, 1));

        Polyhedron P3(P2);
        CGAL_assertion( P3.is_tetrahedron(P3.halfedges_begin()));
        P3.inside_out();
        P3.inside_out();
        P2 = P3;
        CGAL_assertion( P2.is_tetrahedron(P2.halfedges_begin()));
        P2.inside_out();
        CGAL_assertion( P2.is_tetrahedron(P2.halfedges_begin()));
    }
    {
        // Check use of kernel as traits class.
        PolyhedronN P;
        CGAL_assertion_code( HalfedgeN_handle h =) P.make_triangle();
        CGAL_assertion( P.is_valid());
        CGAL_assertion( P.is_triangle( h));
        CGAL_assertion( ! P.is_tetrahedron( h));
        P.normalize_border();
        CGAL_assertion( P.is_valid( false, 1));
        P.inside_out();
        CGAL_assertion( P.is_valid( false, 1));
        CGAL_assertion( P.is_triangle( h));
        P.normalize_border();
        CGAL_assertion( P.is_valid( false, 1));
        CGAL_assertion( P.is_triangle( h));
    }
}
#endif // CGAL_USE_POLYHEDRON_DESIGN_ONE //



int main() {
// This is the test file for the new design. Skip new design test for 
// disabled compilers.
#ifndef CGAL_USE_POLYHEDRON_DESIGN_ONE
    test_Polyhedron();
#endif
    return 0;
}
// EOF //
