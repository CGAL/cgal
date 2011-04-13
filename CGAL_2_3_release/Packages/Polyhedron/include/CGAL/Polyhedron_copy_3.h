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
// file          : Polyhedron_copy_3.h
// chapter       : $CGAL_Chapter: 3D-Polyhedral Surfaces $
// package       : $CGAL_Package: Polyhedron 2.9 (13 Sep 2000) $
// source        : polyhedron_builder.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : MPI Saarbruecken (Stefan Schirra <stschirr@mpi-sb.mpg.de>)
//
// Copy of Polyhedral Surfaces.
// ============================================================================

#ifndef CGAL_POLYHEDRON_COPY_3_H
#define CGAL_POLYHEDRON_COPY_3_H 1

#include <CGAL/basic.h>
// MS Visual C++ 6.0 does not work with the new design.
#if defined( _MSC_VER) && (_MSC_VER <= 1200)
#ifndef CGAL_USE_POLYHEDRON_DESIGN_TWO
#define CGAL_USE_POLYHEDRON_DESIGN_ONE 1
#endif
#endif
#ifndef CGAL_USE_POLYHEDRON_DESIGN_ONE
#define CGAL_USE_POLYHEDRON_DESIGN_TWO 1
#endif

#include <CGAL/Modifier_base.h>
#include <CGAL/Inverse_index.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

CGAL_BEGIN_NAMESPACE

template < class Poly, class HDS >
class Polyhedron_copy_3 : public Modifier_base<HDS> {
protected:
    const Poly& source;
public:
    typedef Poly Polyhedron;
    typedef HDS  Halfedge_data_structure;
    Polyhedron_copy_3( const Polyhedron& poly) : source(poly) {}
        // creates the copy modifier and stores the `source' in its
        // internal state.
    void operator()( HDS& target);
        // copies the `source' known from creation time into the `target'.
        // Postcondition: the `target' is a valid polyhedral surface.
};

template < class Poly, class HDS>
void
Polyhedron_copy_3<Poly,HDS>:: operator()( HDS& target) {
    typedef typename Poly::Vertex_const_iterator  Vertex_const_iterator;
    typedef typename Poly::Facet_const_iterator   Facet_const_iterator;
    typedef Inverse_index< Vertex_const_iterator> Index;
    typedef typename HDS::Point                   Point;

#ifdef CGAL_USE_POLYHEDRON_DESIGN_ONE
    target.delete_all();
#else
    target.clear();
#endif
    Polyhedron_incremental_builder_3<HDS> B( target);
    B.begin_surface( source.size_of_vertices(),
                     source.size_of_facets(),
                     source.size_of_halfedges());
    for ( Vertex_const_iterator vi = source.vertices_begin();
          vi != source.vertices_end();
          ++vi) {
        B.add_vertex( Point( vi->point()));
    }
    Index index( source.vertices_begin(), source.vertices_end());

    for ( Facet_const_iterator fi = source.facets_begin();
          fi != source.facets_end();
          ++fi) {
        B.begin_facet();
        typedef typename Poly::Halfedge_around_facet_const_circulator
            Halfedge_around_facet_const_circulator;
        Halfedge_around_facet_const_circulator hc = fi->facet_begin();
        Halfedge_around_facet_const_circulator hc_end = hc;
        CGAL_assertion( hc != NULL);
        do {
#ifdef CGAL_USE_POLYHEDRON_DESIGN_ONE
            B.add_vertex_to_facet( index[ Vertex_const_iterator(
                hc->vertex().ptr())]);
#else
            B.add_vertex_to_facet( index[ hc->vertex()]);
#endif
            ++hc;
        } while( hc != hc_end);
        B.end_facet();
    }
    B.end_surface();
    target.normalize_border();
}

CGAL_END_NAMESPACE
#endif // CGAL_POLYHEDRON_COPY_3_H //
// EOF //
