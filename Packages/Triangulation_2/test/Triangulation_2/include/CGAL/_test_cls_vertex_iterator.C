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
// file          : include/CGAL/_test_cls_vertex_iterator.C
// revision      : 
// revision_date : 
// author(s)     : Herve Bronnimann (Herve.Bronnimann@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================


template < class Triangulation >
int
CGAL__test_cls_vertex_iterator( const Triangulation &T )
{
  typedef typename Triangulation::Vertex_iterator Vertex_iterator;

  int n = 0;
  Vertex_iterator vit;
  for (vit = T.vertices_begin(); vit != T.vertices_end(); ++vit)
    n++;
  assert( n == T.number_of_vertices() );

  return n;
}
