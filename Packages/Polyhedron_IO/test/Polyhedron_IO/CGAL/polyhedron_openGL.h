#line 12 "cgal_header.fw"
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
// file          : polyhedron_openGL.h
// source        : polyhedron_io.fw
#line 26 "cgal_header.fw"
                                   
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : MPI Saarbruecken (Stefan Schirra <stschirr@mpi-sb.mpg.de>)
//
// Render a Polyhedron<Traits> using the openGL library.
// ============================================================================
#line 279 "polyhedron_io.fw"

#ifndef CGAL_POLYHEDRON_OPENGL_H
#define CGAL_POLYHEDRON_OPENGL_H 1
#line 2635 "polyhedron_io.fw"
#ifndef __gl_h_
#include <GL/gl.h>
#endif
#ifndef CGAL_PRINT_OPENGL_H
#include <CGAL/print_openGL.h>
#endif
#ifndef CGAL_POLYHEDRON_H
#include <CGAL/Polyhedron.h>
#endif

template <class Traits>
void
CGAL_print_openGL( const CGAL_Polyhedron<Traits>& P) {
    // Renders P into the current openGL window. Draws only polygons.
    typedef CGAL_Polyhedron<Traits>                    Poly;
    typedef Poly::Facet_const_iterator             FCI;
    typedef Poly::Halfedge_facet_const_circulator  HFCC;
    for( FCI fi = P.facets_begin(); fi != P.facets_end(); ++fi) {
        HFCC hc = (*fi).facet_begin();
        HFCC hc_end = hc;
        CGAL_assertion( hc != NULL);
        glBegin( GL_POLYGON);
        do {
            CGAL_openGL_vertex((*hc).vertex()->point());
            ++hc;
        } while( hc != hc_end);
        glEnd();
    }
}
#line 282 "polyhedron_io.fw"
 
#line 2667 "polyhedron_io.fw"
template <class Traits>
void
CGAL_print_openGL_vertices( const CGAL_Polyhedron<Traits>& P) {
    // Renders P into the current openGL window.
    // Draws polygon vertices.
    typedef CGAL_Polyhedron<Traits>                 Poly;
    typedef Poly::Vertex_const_iterator            VCI;
    glBegin( GL_POINTS);
    for ( VCI i = P.vertices_begin(); i != P.vertices_end(); ++i) {
        CGAL_openGL_vertex((*i).point());
    }
    glEnd();
}
#line 283 "polyhedron_io.fw"
 
#line 2683 "polyhedron_io.fw"
template <class Traits>
void
CGAL_print_openGL_normal( const CGAL_Polyhedron<Traits>& P) {
    // Renders P into the current openGL window.
    // Draws polygons with normals.
    typedef CGAL_Polyhedron<Traits>                    Poly;
    typedef Poly::Facet_const_iterator             FCI;
    typedef Poly::Halfedge_facet_const_circulator  HFCC;
    for( FCI fi = P.facets_begin(); fi != P.facets_end(); ++fi) {
        HFCC hc = (*fi).facet_begin();
        HFCC hc_end = hc;
        CGAL_assertion( hc != NULL);
        CGAL_openGL_normal((*fi).normal());
        glBegin( GL_POLYGON);
        do {
            CGAL_openGL_vertex((*hc).vertex()->point());
            ++hc;
        } while( hc != hc_end);
        glEnd();
    }
}
#line 284 "polyhedron_io.fw"
 
#line 2707 "polyhedron_io.fw"
#ifndef CGAL_CENTER_POINTS_H
#include <CGAL/center_point.h>
#endif

template <class Traits>
void
CGAL_print_normals_openGL( const CGAL_Polyhedron<Traits>& P,
                          double length,
                          double tip_ratio = 0.3) {
    // Renders the normal vectors of P into the current openGL window.
    typedef CGAL_Polyhedron<Traits>                    Poly;
    typedef Poly::Facet_const_iterator             FCI;
    typedef Poly::Halfedge_facet_const_circulator  HFCC;
    typedef Poly::Point                            Point;
    typedef Poly::Normal                           Normal;
    for( FCI fi = P.facets_begin(); fi != P.facets_end(); ++fi) {
        HFCC hc = (*fi).facet_begin();
        HFCC hc_end = hc;
        CGAL_assertion( hc != NULL);
        Point   foot;
        CGAL_center_point( (*fi), foot);
        Normal  ortho( 0.0, 0.0, 0.0);
        do {
            HFCC next_vertex = hc;
            ++ next_vertex;
            ortho = (*next_vertex).vertex()->point() -
                    (*hc).vertex()->point();
            if ( ortho != Normal( 0.0, 0.0, 0.0)) {
                CGAL_print_openGL( foot, (*fi).normal(),
                                 length, tip_ratio, ortho);
                hc = hc_end;
            } else
                ++hc;
        } while( hc != hc_end);
    }
}
#line 285 "polyhedron_io.fw"
 
#line 2746 "polyhedron_io.fw"
template <class Traits>
int
CGAL_print_zero_facet_openGL( const CGAL_Polyhedron<Traits>& P) {
    // Counts and renders facets with zero normals as a closed line.
    typedef CGAL_Polyhedron<Traits>                    Poly;
    typedef Poly::Facet_const_iterator             FCI;
    typedef Poly::Halfedge_facet_const_circulator  HFCC;
    typedef Poly::Point                            Point;
    typedef Poly::Normal                           Normal;
    int count = 0;
    for( FCI fi = P.facets_begin(); fi != P.facets_end(); ++fi) {
        if ( (*fi).normal().sqr_length() == 0) {
            count++;
            HFCC hc = (*fi).facet_begin();
            HFCC hc_end = hc;
            CGAL_assertion( hc != NULL);
            glBegin( GL_LINE_LOOP);
            do {
                CGAL_openGL_vertex( (*hc).vertex()->point());
                ++hc;
            } while( hc != hc_end);
            glEnd();
        }
    }
    return count;
}
#line 286 "polyhedron_io.fw"
 
#line 2775 "polyhedron_io.fw"
// VP could be either a viewing direction (orthogonal projection),
// or a view point (perspective projection). It is internally
// solved nicely by overloading the "is_contour" member function
// for halfedges.
template <class Traits, class VP>
int
CGAL_print_contour_openGL( const CGAL_Polyhedron<Traits>& P, const VP& vp) {
    // Counts and renders contour edges with respect to the
    // viewing direction v.
    typedef CGAL_Polyhedron<Traits>                    Poly;
    typedef Poly::Edge_const_iterator              ECI;
    typedef Poly::Point                            Point;
    typedef Poly::Normal                           Normal;
    int count = 0;
    glBegin( GL_LINES);
    ECI ei = P.edges_begin();
    for( ; ei != P.border_edges_begin(); ++ei) {
        // These are the non-border edges. Check against viewing vector.
        CGAL_assertion( ! (*ei).is_border_edge());
        if ( (*ei).is_contour_non_border( vp)) {
            count++;
            CGAL_openGL_vertex( (*ei).vertex()->point());
            CGAL_openGL_vertex( (*ei).opposite()->vertex()->point());
        }
    }
    for( ; ei != P.edges_end(); ++ei) {
        // The border edges are always contour edges.
        count++;
        CGAL_assertion( (*ei).is_border_edge());
        CGAL_openGL_vertex( (*ei).vertex()->point());
        CGAL_openGL_vertex( (*ei).opposite()->vertex()->point());
    }
    glEnd();
    return count;
}
#line 287 "polyhedron_io.fw"
 
#endif // CGAL_POLYHEDRON_OPENGL_H //
// EOF //
