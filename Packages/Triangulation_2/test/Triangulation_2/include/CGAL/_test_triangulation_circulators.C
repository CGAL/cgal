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
// file          : include/CGAL/_test_triangulation_circulators.C
// revision      : 
// revision_date : 
// author(s)     : Herve Bronnimann (Herve.Bronnimann@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================


#include <CGAL/_test_cls_vertex_circulator.C>
#include <CGAL/_test_cls_edge_circulator.C>
#include <CGAL/_test_cls_face_circulator.C>

template < class Triangulation >
void
_test_circulators( const Triangulation &T )
{
  int nvi = _test_cls_vertex_circulator(T); 
  int nei = _test_cls_edge_circulator(T); 
  int nfi = _test_cls_face_circulator(T); 
  assert( nvi == nei );
  assert( nvi == nfi );

  // test the circulators provided by the Triangulation class 
  typedef typename Triangulation::Vertex_iterator   Vertex_iterator;
  typedef typename Triangulation::Vertex_circulator Vertex_circulator;
  typedef typename Triangulation::Face_iterator     Face_iterator;
  typedef typename Triangulation::Face_circulator   Face_circulator;
  typedef typename Triangulation::Edge_iterator     Edge_iterator;
  typedef typename Triangulation::Edge_circulator   Edge_circulator;

  int n = 0;
  Vertex_iterator vit;
  Vertex_circulator vc, vc0;
  for (vit = T.vertices_begin(); vit != T.vertices_end(); ++vit)
    {
      vc0 = vc = T.incident_vertices( vit, vit->face() );
      do {
	vc++; n++;
      } while (vc != vc0);
    }
  assert(nvi==n);

  n = 0;
  Face_circulator fc, fc0;
  for (vit = T.vertices_begin(); vit != T.vertices_end(); ++vit)
    {
      fc0 = fc = T.incident_faces( vit, vit->face() );
      do {
	fc++; n++;
      } while (fc != fc0);
    }
  assert(nfi==n);

  n = 0;
  Edge_circulator ec, ec0;
  for (vit = T.vertices_begin(); vit != T.vertices_end(); ++vit)
    {
      ec0 = ec = T.incident_edges( vit, vit->face() );
      do {
	ec++; n++;
      } while (ec != ec0);
    }
  assert(nei==n);

  //Traverse convex_hull- this count infinite
  n = 0;
  fc = fc0 = T.incident_faces(T.infinite_vertex());
  do {
    fc++; n++;
  } while (fc != fc0);
  
  //Count finite faces
  int m=0;
  Face_iterator fit;
  for(fit = T.faces_begin(); fit != T.faces_end(); ++fit) {
    m++;
  }
  //Check Euler formula
  assert( n+m == 2*(T.number_of_vertices() +1) -4);
}


