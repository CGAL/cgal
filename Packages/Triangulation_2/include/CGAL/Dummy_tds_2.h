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
// file          : include/CGAL/Dummy_tds_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec
//
// coordinator   : INRIA Sophia Antipolis
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_DUMMY_TDS_2_H
#define CGAL_TRIANGULATION_DUMMY_TDS_2_H

#include <CGAL/basic.h>
#include <CGAL/Triangulation_short_names_2.h>

CGAL_BEGIN_NAMESPACE

// Dummy TDS which provides all types that a vertex_base or cell_base can use.
struct Dummy_tds_2 {
  struct Vertex {};
  struct Face {};
  struct Edge {};

  struct Vertex_handle {};
  struct Face_handle {};

  struct Vertex_iterator {};
  struct Face_iterator {};
  struct Edge_iterator {};

  struct Edge_circulator {};
  struct Facet_circulator {};
  struct Vertex_circulator {};
};

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGULATION_DUMMY_TDS_2_H
