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
// file          : include/CGAL/_test_triangulation_iterators.C
// revision      : 
// revision_date : 
// author(s)     : Herve Bronnimann (Herve.Bronnimann@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================


#include <CGAL/_test_cls_vertex_iterator.C>
#include <CGAL/_test_cls_edge_iterator.C>
#include <CGAL/_test_cls_face_iterator.C>

template < class Triangulation >
void
CGAL__test_iterators( const Triangulation &T )
{
  int nv = CGAL__test_cls_vertex_iterator(T);
  int ne = CGAL__test_cls_edge_iterator(T);
  int nf = CGAL__test_cls_face_iterator(T);
  // cout << "Euler's relation: " << nv -ne + nf << endl;
  assert( nv -ne + nf == 1); // Euler's relation
}
