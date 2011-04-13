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
#ifdef CGAL_USE_POLYHEDRON_DESIGN_TWO
#undef CGAL_USE_POLYHEDRON_DESIGN_TWO
#endif
#define CGAL_USE_POLYHEDRON_DESIGN_ONE 1

#include <CGAL/circulator.h>
#include <CGAL/circulator_impl.h>
#include <list>
#include <vector>

#include <CGAL/Halfedge_data_structure_polyhedron_default_3.h>
#include <CGAL/Halfedge_data_structure_using_vector.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Iterator_project.h>
#include <CGAL/function_objects.h>

using namespace CGAL;

// A polyhedron modifier that creates a tetrahedron using the
// incremental builder.
template < class HDS >
class Build_tetrahedron : public Modifier_base<HDS> {
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
    Polyhedron_incremental_builder_3<HDS> B( target, true);
    B.begin_surface( 4, 4, 12);
    // Point coordinates suitable for integer coordinates.
    typedef typename HDS::Point Point;
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
    typedef Cartesian<double>                   Rep;
    typedef Cartesian<int>                      RepI;
    typedef Point_3<Rep>                        Point;
    typedef Plane_3<Rep>                        Plane;
    typedef Halfedge_data_structure_polyhedron_default_3<Rep>  HDS;
    typedef Halfedge_data_structure_polyhedron_default_3<RepI> HDSI;
    typedef Halfedge_data_structure_using_vector<
                Vertex_max_base<Point>,
                Halfedge_max_base,
                Polyhedron_facet_base_3<Rep> >  HDSV;
    typedef Polyhedron_default_traits_3<Rep>    Traits;
    typedef Polyhedron_default_traits_3<RepI>   TraitsI;

    typedef Polyhedron_3<Traits, HDS>           Polyhedron;
    typedef Polyhedron_3<TraitsI,HDSI>          PolyhedronI;
    typedef Polyhedron_3<Traits, HDSV>          PolyhedronV;

    typedef Polyhedron::Vertex                  Vertex;
    typedef Polyhedron::Halfedge                Halfedge;
    typedef Polyhedron::Facet                   Facet;

    typedef Polyhedron::Edge_iterator           Edge_iterator;
    typedef Polyhedron::Vertex_iterator         Vertex_iterator;
    typedef Polyhedron::Halfedge_handle         Halfedge_handle;
    typedef Polyhedron::Halfedge_iterator       Halfedge_iterator;
    typedef Polyhedron::Halfedge_const_handle   Halfedge_const_handle;
    typedef Polyhedron::Halfedge_const_iterator Halfedge_const_iterator;
    typedef Polyhedron::Halfedge_around_vertex_circulator
                                   Halfedge_around_vertex_circulator;
    typedef Polyhedron::Halfedge_around_vertex_const_circulator
                                   Halfedge_around_vertex_const_circulator;
    typedef Polyhedron::Halfedge_around_facet_circulator
                                   Halfedge_around_facet_circulator;
    typedef Polyhedron::Halfedge_around_facet_const_circulator
                                   Halfedge_around_facet_const_circulator;

    typedef PolyhedronV::Halfedge_handle        HalfedgeV_handle;

    {
        // Check if all automatic conversions work as promised.
        // First the conversions to Halfedge_const_iterator.
        // Simultaneously all iterators/circulators are checked for
        // initialization with a pointer to Halfedge.
        Polyhedron P;
        // Halfedge* h = P.make_triangle().ptr();
        // CGAL_assertion( P.is_valid());
        // CGAL_assertion( P.is_triangle( h));  // Accept pointer.
        // const Halfedge* ch = h;
        // CGAL_assertion( P.is_triangle( ch));

        Halfedge_handle  hh = P.make_triangle();  // Check handle.
        CGAL_assertion( P.is_triangle( hh));
        CGAL_assertion_code( Halfedge_const_handle  hch(hh);)
        CGAL_assertion( P.is_triangle( hch));

        Halfedge_iterator  hi(hh);          // Check Iterator.
        CGAL_assertion( P.is_triangle( hi));
        CGAL_assertion_code( Halfedge_const_iterator  hci(hh);)
        CGAL_assertion( P.is_triangle( hci));

        Halfedge_around_facet_circulator  hfc(hh);  // Check circulator 1.
        CGAL_assertion( P.is_triangle( hfc));  //SunPro
        CGAL_assertion_code( Halfedge_around_facet_const_circulator  hfcc(hh);)
        CGAL_assertion( P.is_triangle( hfcc));

        Halfedge_around_vertex_circulator  hvc(hh);  // Check circulator 2.
        CGAL_assertion( P.is_triangle( hvc));  // SunPro
        CGAL_assertion_code(Halfedge_around_vertex_const_circulator  hvcc(hh);)
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
        P2 = P3;
        CGAL_assertion( P2.is_tetrahedron(P2.halfedges_begin()));
        P2.inside_out();
        CGAL_assertion( P2.is_tetrahedron(P2.halfedges_begin()));
    }
    {
        // Check the easy creation of a point iterator.
        Polyhedron P;
        CGAL_assertion_code( Halfedge_handle h = P.make_tetrahedron();)
        CGAL_assertion( P.is_valid());
        CGAL_assertion( P.is_tetrahedron( h));

        typedef Project_point<Vertex>        Project_point;
        typedef Polyhedron::Difference           Difference;
        typedef Polyhedron::iterator_category    iterator_category;
        typedef Iterator_project<Vertex_iterator, Project_point,
            Point&, Point*, Difference, iterator_category>
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
        CGAL_assertion_code( Halfedge_handle h = P.make_triangle(); )
        CGAL_assertion( P.is_valid());
        CGAL_assertion( P.is_triangle( h));

        CGAL_assertion_code( Halfedge_handle g = 
                                P.add_vertex_and_facet_to_border(
                                h->next()->opposite(),
                                h->opposite()); )
        CGAL_assertion( P.is_valid());
        CGAL_assertion( h->next()->next()->next() == h);
        CGAL_assertion( g->next()->next()->next() == g);
        CGAL_assertion( h->opposite()->next() == g);

        CGAL_assertion_code( Halfedge_handle gg = P.add_facet_to_border(
                                 h->next()->next()->opposite(),
                                 g->next()->opposite()); )
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
        if ( P.size_of_vertices() > 28)  // fool SGI compiler
            std::cerr << P.size_of_vertices() << ' '
                      << P.size_of_halfedges() << ' '
                      << P.size_of_facets() << std::endl;
        CGAL_assertion( P.is_valid());
        CGAL_assertion( P.size_of_vertices()  == 4);
        CGAL_assertion( P.size_of_halfedges() == 10);
        CGAL_assertion( P.size_of_facets()    == 2);
        P.clear();
        P.erase_all();
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
        /* Check it out with an edge iterator. */
        Edge_iterator h = P.edges_begin();
        /* The first three edges must be non-border edges. */
        CGAL_assertion( ! h->is_border_edge());
        ++h;
        CGAL_assertion( ! h->is_border_edge());
        ++h;
        CGAL_assertion( ! h->is_border_edge());
        ++h;
        /* Here the three border edges should start. */
        CGAL_assertion( h == P.border_edges_begin());
        CGAL_assertion( h->is_border_edge());
        CGAL_assertion( h->opposite()->is_border());
        ++h;
        CGAL_assertion( h->is_border_edge());
        CGAL_assertion( h->opposite()->is_border());
        ++h;
        CGAL_assertion( h->is_border_edge());
        CGAL_assertion( h->opposite()->is_border());
        ++h;
        CGAL_assertion( h == P.edges_end());
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
                   facet()->halfedge()->
                   facet() == h->opposite()->next_on_vertex()->
                   next_on_vertex()->facet());
        // The changed facets.
        CGAL_assertion( h->facet()->halfedge()->facet() ==
                   h->facet());
        CGAL_assertion( h->next_on_vertex()->facet()->halfedge()->facet()==
                   h->next_on_vertex()->facet());
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
        P2 = P3;
        CGAL_assertion( P2.is_tetrahedron( P2.halfedges_begin()));
        P2.inside_out();
        CGAL_assertion( P2.is_tetrahedron( P2.halfedges_begin()));
    }
}


int main(){
    test_Polyhedron();
    return 0;
}
// EOF //
