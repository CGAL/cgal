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


//#include <CGAL/_test_cls_vertex_circulator.C>
//#include <CGAL/_test_cls_edge_circulator.C>
//#include <CGAL/_test_cls_face_circulator.C>

template < class Triangulation >
void
_test_circulators( const Triangulation &T )
{
  // test the circulators provided by the Triangulation class 
  typedef typename Triangulation::All_vertices_iterator All_vertices_iterator;
  typedef typename Triangulation::All_faces_iterator    All_faces_iterator;
  typedef typename Triangulation::All_edges_iterator    All_edges_iterator;
  typedef typename Triangulation::Face_circulator   Face_circulator;
  typedef typename Triangulation::Vertex_circulator Vertex_circulator;
  typedef typename Triangulation::Edge_circulator   Edge_circulator;

  int nvi = 0;
  All_vertices_iterator vit;
  Vertex_circulator vc, vc0;
  for (vit = T.all_vertices_begin(); vit != T.all_vertices_end(); ++vit)
    {
      vc0 = vc = T.incident_vertices( vit, vit->face() );
      if( !vc.is_empty()){
	do {
	  vc++; nvi++;
	} while (vc != vc0);
      }
    }
  
  int nfi = 0;
  Face_circulator fc, fc0;
  for (vit = T.all_vertices_begin(); vit != T.all_vertices_end(); ++vit)
    {
      fc0 = fc = T.incident_faces( vit, vit->face() );
      if( !fc.is_empty()){
	do {
	  fc++; nfi++;
	} while (fc != fc0);
      }
    }
  

  int nei = 0;
  Edge_circulator ec, ec0;
  for (vit = T.all_vertices_begin(); vit != T.all_vertices_end(); ++vit)
    {
      ec0 = ec = T.incident_edges( vit, vit->face() );
       if( !ec.is_empty()){
	 do {
	   ec++; nei++;
	 } while (ec != ec0);
       }
    }
  
  //Traverse convex_hull- this count infinite
  int nch = 0;
  fc = fc0 = T.incident_faces(T.infinite_vertex());
  if( !fc.is_empty()){
    do {
      fc++; nch++;
    } while (fc != fc0);
  }
  
  //Check Total 
  int mf = T.number_of_faces() + nch;
  int mv = T.number_of_vertices() +1;
  int me ;
  if (T.dimension() <= 0) me=0;
  else me = T.dimension() == 1 ? mv : 3*mv-6 ;
    
  assert ( nvi ==  nei);
  assert ( nvi == 2*me);
  assert ( nfi == 3*mf);
 
}


