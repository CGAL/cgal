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
// file          : Halfedge_data_structure_default.h
// chapter       : $CGAL_Chapter: Halfedge Data Structures $
// package       : $CGAL_Package: Halfedge_DS 2.8 (13 Sep 2000) $
// source        : hds.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : MPI Saarbruecken (Stefan Schirra <stschirr@mpi-sb.mpg.de>)
//
// Halfedge Data Structure Default Implementation for CGAL.
// ============================================================================

#ifndef CGAL_HALFEDGE_DATA_STRUCTURE_DEFAULT_H
#define CGAL_HALFEDGE_DATA_STRUCTURE_DEFAULT_H 1
#ifndef CGAL_HALFEDGE_DATA_STRUCTURE_BASES_H
#include <CGAL/Halfedge_data_structure_bases.h>
#endif

#ifndef CGAL_HALFEDGE_DATA_STRUCTURE_USING_LIST_H
#include <CGAL/Halfedge_data_structure_using_list.h>
#endif

CGAL_BEGIN_NAMESPACE

template <class Pt>
class Halfedge_data_structure_default
    : public Halfedge_data_structure_using_list<
          Vertex_max_base<Pt>, Halfedge_max_base, Facet_max_base> {
public:  // CREATION
    typedef typename   Halfedge_data_structure_using_list<
        Vertex_max_base<Pt>, Halfedge_max_base, Facet_max_base>::Size Size;
    Halfedge_data_structure_default() {}
    Halfedge_data_structure_default( Size v, Size h, Size f)
    : Halfedge_data_structure_using_list<
        Vertex_max_base<Pt>, Halfedge_max_base, Facet_max_base> (v,h,f) {}
};

CGAL_END_NAMESPACE
#endif // CGAL_HALFEDGE_DATA_STRUCTURE_DEFAULT_H //
// EOF //
