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
// file          : Halfedge_data_structure_polyhedron_default_3.h
// chapter       : $CGAL_Chapter: Halfedge Data Structures $
// package       : $CGAL_Package: Halfedge_DS 2.8 (13 Sep 2000) $
// source        : hds.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : MPI Saarbruecken (Stefan Schirra <stschirr@mpi-sb.mpg.de>)
//
// Halfedge Data Structure Default for Polyhedral Surfaces.
// ============================================================================

#ifndef CGAL_HALFEDGE_DATA_STRUCTURE_POLYHEDRON_DEFAULT_3_H
#define CGAL_HALFEDGE_DATA_STRUCTURE_POLYHEDRON_DEFAULT_3_H 1
#ifndef CGAL_HALFEDGE_DATA_STRUCTURE_BASES_H
#include <CGAL/Halfedge_data_structure_bases.h>
#endif
#ifndef CGAL_HALFEDGE_DATA_STRUCTURE_USING_LIST_H
#include <CGAL/Halfedge_data_structure_using_list.h>
#endif
#ifndef CGAL_POINT_3_H
#include <CGAL/Point_3.h>
#endif
#ifndef CGAL_PLANE_3_H
#include <CGAL/Plane_3.h>
#endif
#ifndef CGAL_VECTOR_3_H
#include <CGAL/Vector_3.h>
#endif

CGAL_BEGIN_NAMESPACE

template <class Traits_>
class Halfedge_data_structure_polyhedron_default_3
    : public Halfedge_data_structure_using_list<
          Vertex_max_base< CGAL_TYPENAME_MSVC_NULL Traits_::Point_3 >,
          Halfedge_max_base,
          Polyhedron_facet_base_3<Traits_>
      > {
  public:  // CREATION
    typedef Traits_                   Traits;
    typedef typename Traits::Point_3  Point_3;
    typedef typename Halfedge_data_structure_using_list<
          Vertex_max_base< Point_3>,
          Halfedge_max_base,
          Polyhedron_facet_base_3<Traits>
      >::Size Size;
    Halfedge_data_structure_polyhedron_default_3() {}
    Halfedge_data_structure_polyhedron_default_3( Size v, Size h, Size f)
        : Halfedge_data_structure_using_list<
              Vertex_max_base< Point_3>,
              Halfedge_max_base,
              Polyhedron_facet_base_3<Traits>
          > (v,h,f) {}
};

CGAL_END_NAMESPACE
#endif // CGAL_HALFEDGE_DATA_STRUCTURE_POLYHEDRON_DEFAULT_3_H //
// EOF //
