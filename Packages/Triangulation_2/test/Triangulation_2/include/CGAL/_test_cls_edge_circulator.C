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
// file          : include/CGAL/_test_cls_edge_circulator.C
// revision      : 
// revision_date : 
// author(s)     : Herve Bronnimann (Herve.Bronnimann@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================

template < class Triangulation >
int
CGAL__test_cls_edge_circulator( const Triangulation &T )
{
  typedef typename Triangulation::Vertex_iterator   Vertex_iterator;
  typedef typename Triangulation::Edge_circulator Edge_circulator;

  int n = 0;
  Vertex_iterator vit;
  Edge_circulator ec, ec0;
  for (vit = T.vertices_begin(); vit != T.vertices_end(); ++vit)
    {
      ec0 = ec = vit->incident_edges( vit->face() );
      //n++;
      do {
	ec++; n++;
      } while (ec != ec0);
    }

  return n;
}
