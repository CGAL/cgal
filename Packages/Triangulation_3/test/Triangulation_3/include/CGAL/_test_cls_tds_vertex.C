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
// author(s)     : Rebufat Francois
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================

#include <cassert>

template <class Vertex>
void
_test_vertex_tds_3(const Vertex &)
{
  typedef typename Vertex::Triangulation_data_structure  Tds;
  typedef typename Tds::Cell_handle                      Cell_handle;
  typedef typename Tds::Vertex_handle                    Vertex_handle;

  Tds tds;

  Cell_handle c1 = tds.create_cell();
  Vertex_handle v1 = tds.create_vertex();
  v1->set_cell(c1);
  assert(v1->cell() == c1);
  c1->set_vertex(0, v1);
  assert(v1->is_valid());

  Cell_handle c2 = tds.create_cell();
  v1->set_cell(c2);
  c2->set_vertex(0, v1);
  assert(v1->is_valid());

  // Compatibility of handles with NULL.
  Vertex_handle v = NULL;
  v = NULL;
  assert(v == NULL);
  assert(v1 != NULL);

  Cell_handle c = NULL;
  c = NULL;
  assert(c == NULL);
  assert(c1 != NULL);

  // We want the following comparisons to work for use in std::set<>...
  bool b1 = v<v1;  (void) b1;
  bool b2 = c<c1;  (void) b2;
}
