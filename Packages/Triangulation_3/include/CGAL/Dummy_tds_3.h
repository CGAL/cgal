// ============================================================================
//
// Copyright (c) 2003 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Dummy_tds_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Sylvain Pion
//
// coordinator   : INRIA Sophia Antipolis
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_DUMMY_TDS_3_H
#define CGAL_TRIANGULATION_DUMMY_TDS_3_H

#include <CGAL/basic.h>
#include <CGAL/Triangulation_short_names_3.h>

CGAL_BEGIN_NAMESPACE

// Dummy TDS which provides all types that a vertex_base or cell_base can use.
struct Dummy_tds_3 {
  struct Vertex {};
  struct Cell {};
  struct Facet {};
  struct Edge {};

  struct Vertex_handle {};
  struct Cell_handle {};

  struct Vertex_iterator {};
  struct Cell_iterator {};
  struct Facet_iterator {};
  struct Edge_iterator {};

  struct Cell_circulator {};
  struct Facet_circulator {};
};

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGULATION_DUMMY_TDS_3_H
