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
// file          : HalfedgeDS_min_items.h
// chapter       : $CGAL_Chapter: Halfedge Data Structures $
// package       : $CGAL_Package: HalfedgeDS 3.3 (27 Sep 2000) $
// source        : hds_bases.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : MPI Saarbruecken (Stefan Schirra <stschirr@mpi-sb.mpg.de>)
//
// Minimal items class for HDS (Vertex, Halfedge, Face).
// ============================================================================

#ifndef CGAL_HALFEDGEDS_MIN_ITEMS_H
#define CGAL_HALFEDGEDS_MIN_ITEMS_H 1
#ifndef CGAL_HALFEDGEDS_VERTEX_MIN_BASE_H
#include <CGAL/HalfedgeDS_vertex_min_base.h>
#endif
#ifndef CGAL_HALFEDGEDS_HALFEDGE_MIN_BASE_H
#include <CGAL/HalfedgeDS_halfedge_min_base.h>
#endif
#ifndef CGAL_HALFEDGEDS_FACE_MIN_BASE_H
#include <CGAL/HalfedgeDS_face_min_base.h>
#endif

CGAL_BEGIN_NAMESPACE

class HalfedgeDS_min_items {
public:
    template < class Refs, class Traits>
    struct Vertex_wrapper {
        typedef HalfedgeDS_vertex_min_base< Refs>   Vertex;
    };
    template < class Refs, class Traits>
    struct Halfedge_wrapper {
        typedef HalfedgeDS_halfedge_min_base< Refs> Halfedge;
    };
    template < class Refs, class Traits>
    struct Face_wrapper {
        typedef HalfedgeDS_face_min_base< Refs>     Face;
    };
};

CGAL_END_NAMESPACE
#endif // CGAL_HALFEDGEDS_MIN_ITEMS_H //
// EOF //
