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
// file          : include/CGAL/_test_cls_edge_iterator.C
// revision      : 
// revision_date : 
// author(s)     : Herve Bronnimann (Herve.Bronnimann@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================

template < class Triangulation >
int
_test_cls_edge_iterator( const Triangulation &T )
{
  typedef typename Triangulation::Edge_iterator   Edge_iterator;

  int n = 0;
  Edge_iterator eit;
  for (eit = T.edges_begin(); eit != T.edges_end(); ++eit)
    n++;

  return n;
}
