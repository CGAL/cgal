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
// file          : include/CGAL/_test_triangulation_circulators.h
// revision      : 
// revision_date : 
// author(s)     : Herve Bronnimann (Herve.Bronnimann@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================


//#include <CGAL/_test_cls_vertex_circulator.h>
//#include <CGAL/_test_cls_edge_circulator.h>
//#include <CGAL/_test_cls_face_circulator.h>

template < class Tr>
void
_test_circulators( const Tr &T )
{
  // test the circulators provided by the Triangulation class 
  typedef typename Tr::All_vertices_iterator All_vertices_iterator;
  typedef typename Tr::All_faces_iterator    All_faces_iterator;
  typedef typename Tr::All_edges_iterator    All_edges_iterator;
  typedef typename Tr::Face_circulator   Face_circulator;
  typedef typename Tr::Vertex_circulator Vertex_circulator;
  typedef typename Tr::Edge_circulator   Edge_circulator;

  CGAL_USE_TYPE(All_faces_iterator);
  CGAL_USE_TYPE(All_edges_iterator);

  typename Tr::size_type nvi = 0;
  typename Tr::size_type nvi_r = 0;
  All_vertices_iterator vit;
  Vertex_circulator vc, vc0;
  for (vit = T.all_vertices_begin(); vit != T.all_vertices_end(); ++vit)
    {
      vc0 = vc = T.incident_vertices( vit, vit->face() );
      if( !vc.is_empty()){
	if( vc != NULL){
	  do {
	    vc++; nvi++;
	  } while (vc != vc0);
	}
      }
      //test operator --()
       vc0 = vc = T.incident_vertices( vit, vit->face() );
      if( !vc.is_empty()){
	if( vc != NULL){
	  do {
	    vc--; nvi_r++;
	  } while (vc != vc0);
	}
      }
      assert(nvi_r == nvi);
    }
  
  typename Tr::size_type nfi = 0;
  typename Tr::size_type nfi_r = 0;
  Face_circulator fc, fc0;
  for (vit = T.all_vertices_begin(); vit != T.all_vertices_end(); ++vit)
    {
      fc0 = fc = T.incident_faces( vit, vit->face() );
      if( !fc.is_empty()){
	do {
	  fc++; nfi++;
	} while (fc != fc0);
      }
      //test operator --()
      fc0 = fc = T.incident_faces( vit, vit->face() );
      if( !fc.is_empty()){
	do {
	  fc--; nfi_r++;
	} while (fc != fc0);
      }
      assert(nfi_r == nfi);
    }
  

  typename Tr::size_type nei = 0;
  typename Tr::size_type nei_r = 0;
  Edge_circulator ec, ec0;
  for (vit = T.all_vertices_begin(); vit != T.all_vertices_end(); ++vit)
    {
      ec0 = ec = T.incident_edges( vit, vit->face() );
       if( !ec.is_empty()){
	 do {
	   ec++; nei++;
	 } while (ec != ec0);
       }
       //test operator --()
       ec0 = ec = T.incident_edges( vit, vit->face() );
       if( !ec.is_empty()){
	 do {
	   ec--; nei_r++;
	 } while (ec != ec0);
       }
       assert(nei_r == nei);
    }
  
  //Traverse convex_hull- this count infinite
  //using pre incrementation to test it
  int nch = 0;
  fc = fc0 = T.incident_faces(T.infinite_vertex());
  if( !fc.is_empty()){
    do {
      ++fc; nch++;
    } while (fc != fc0);
  }
  
  //Check Total 
  typename Tr::size_type mf = T.number_of_faces() + nch;
  typename Tr::size_type mv = T.number_of_vertices() +1;
  typename Tr::size_type me ;
  if (T.dimension() <= 0) me=0;
  else me = T.dimension() == 1 ? mv : 3*mv-6 ;
    
  assert ( nvi ==  nei);
  assert ( nvi == 2*me);
  assert ( nfi == 3*mf);

}


