// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
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
// source        : 
// file          : include/CGAL/_test_cls_tds_vertex.C
// revision      : 
// revision_date : 
// author(s)     : Rebufat Francois (Francois.Rebufat@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================

#include <cassert>

template <class Vertex>
void
_test_vertex_tds_3(const Vertex &)
{
  typedef typename Vertex::Point             Point;
  typedef typename Vertex::Cell              Cell;

  // Build object 
  Cell c1;
  Point p;
  Vertex v1(p,&c1);
  assert(v1.cell()==&c1);
  c1.set_vertex(0,&v1);
  assert(v1.is_valid());
  Cell c2;
  v1.set_cell(&c2);
  c2.set_vertex(0,&v1);
  assert(v1.is_valid());
}
