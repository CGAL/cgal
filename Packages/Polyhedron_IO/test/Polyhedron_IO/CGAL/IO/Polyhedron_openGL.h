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
// file          : Polyhedron_openGL.h
// package       : $CGAL_Package: $
// source        : polyhedron_io.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : Herve Bronnimann  <Herve.Bronnimann@sophia.inria.fr>
//
// Render a Polyhedron_3 using the openGL library.
// ============================================================================

#ifndef CGAL_IO_POLYHEDRON_OPENGL_H
#define CGAL_IO_POLYHEDRON_OPENGL_H 1
#ifndef CGAL_PROTECT_GL_GL_H
#include <GL/gl.h>
#define CGAL_PROTECT_GL_GL_H
#endif // CGAL_PROTECT_GL_GL_H
#ifndef CGAL_PRINT_OPENGL_H
#include <CGAL/print_openGL.h>
#endif
#ifndef CGAL_POLYHEDRON_3_H
#include <CGAL/Polyhedron_3.h>
#endif

template <class Traits, class HDS>
void
CGAL_print_openGL( const CGAL_Polyhedron_3<Traits,HDS>& P) {
    // Renders P into the current openGL window. Draws only polygons.
    typedef CGAL_Polyhedron_3<Traits,HDS>                    Poly;
    typedef typename Poly::Facet_const_iterator             FCI;
    typedef typename Poly::Halfedge_around_facet_const_circulator  HFCC;
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
template <class Traits, class HDS>
void
CGAL_print_openGL_vertices( const CGAL_Polyhedron_3<Traits,HDS>& P) {
    // Renders P into the current openGL window.
    // Draws polygon vertices.
    typedef CGAL_Polyhedron_3<Traits,HDS>                    Poly;
    typedef typename Poly::Vertex_const_iterator            VCI;
    glBegin( GL_POINTS);
    for ( VCI i = P.vertices_begin(); i != P.vertices_end(); ++i) {
        CGAL_openGL_vertex((*i).point());
    }
    glEnd();
}
template <class Traits, class HDS>
void
CGAL_print_openGL_normal( const CGAL_Polyhedron_3<Traits,HDS>& P) {
    // Renders P into the current openGL window.
    // Draws polygons with normals.
    typedef CGAL_Polyhedron_3<Traits,HDS>                    Poly;
    typedef typename Poly::Facet_const_iterator             FCI;
    typedef typename Poly::Halfedge_around_facet_const_circulator  HFCC;
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
#ifndef CGAL_CENTER_POINT_3_H
#include <CGAL/center_point_3.h>
#endif // CGAL_CENTER_POINT_3_H

template <class Traits, class HDS>
void
CGAL_print_normals_openGL( const CGAL_Polyhedron_3<Traits,HDS>& P,
                          double length,
                          double tip_ratio = 0.3) {
    // Renders the normal vectors of P into the current openGL window.
    typedef CGAL_Polyhedron_3<Traits,HDS>                    Poly;
    typedef typename Poly::Facet_const_iterator             FCI;
    typedef typename Poly::Halfedge_around_facet_const_circulator  HFCC;
    typedef typename Poly::Point                            Point;
    typedef typename Poly::Normal                           Normal;
    for( FCI fi = P.facets_begin(); fi != P.facets_end(); ++fi) {
        HFCC hc = (*fi).facet_begin();
        HFCC hc_end = hc;
        CGAL_assertion( hc != NULL);
        Point   foot;
        CGAL_center_point_3( (*fi), foot);
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
template <class Traits, class HDS>
int
CGAL_print_zero_facet_openGL( const CGAL_Polyhedron_3<Traits,HDS>& P) {
    // Counts and renders facets with zero normals as a closed line.
    typedef CGAL_Polyhedron_3<Traits,HDS>                    Poly;
    typedef typename Poly::Facet_const_iterator             FCI;
    typedef typename Poly::Halfedge_around_facet_const_circulator  HFCC;
    typedef typename Poly::Point                            Point;
    typedef typename Poly::Normal                           Normal;
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
template < class Halfedge_handle, class NT, class Traits>
bool CGAL__is_contour_edge( Halfedge_handle h, const CGAL_Point_3<NT>& p,
                           const Traits& traits) {
    return traits.is_contour( h->facet()->plane(),
                              h->opposite()->facet()->plane(), p);
}

template < class Halfedge_handle, class NT, class Traits>
bool CGAL__is_contour_edge( Halfedge_handle h, const CGAL_Vector_3<NT>& v,
                           const Traits& traits) {
    return traits.is_contour( h->facet()->normal(),
                              h->opposite()->facet()->normal(), v);
}


// VP could be either a viewing direction (orthogonal projection),
// or a view point (perspective projection). It is internally
// solved nicely by overloading the "is_contour_edge" function
// above.
template <class Traits, class HDS, class VP>
int
CGAL_print_contour_openGL( const CGAL_Polyhedron_3<Traits,HDS>& P,
                          const VP& vp) {
    // Counts and renders contour edges with respect to the
    // viewing direction v.
    typedef CGAL_Polyhedron_3<Traits,HDS>                    Poly;
    typedef typename Poly::Halfedge_const_iterator          HECI;
    typedef typename Poly::Point                            Point;
    typedef typename Poly::Normal                           Normal;
    int count = 0;
    glBegin( GL_LINES);
    HECI h = P.halfedges_begin();
    for( ; h != P.border_halfedges_begin(); ++h,++h) {
        // These are the non-border edges. Check against viewing vector.
        CGAL_assertion( ! h->is_border_edge());
        if ( CGAL__is_contour_edge( h, vp, P.traits())) {
            count++;
            CGAL_openGL_vertex( h->vertex()->point());
            CGAL_openGL_vertex( h->opposite()->vertex()->point());
        }
    }
    for( ; h != P.halfedges_end(); ++h, ++h) {
        // The border edges are always contour edges.
        count++;
        CGAL_assertion( h->is_border_edge());
        CGAL_openGL_vertex( h->vertex()->point());
        CGAL_openGL_vertex( h->opposite()->vertex()->point());
    }
    glEnd();
    return count;
}
#endif // CGAL_IO_POLYHEDRON_OPENGL_H //
// EOF //
