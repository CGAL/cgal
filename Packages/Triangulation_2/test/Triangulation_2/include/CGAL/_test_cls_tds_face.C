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
// file          : include/CGAL/_test_cls_tds_face.C
// revision      : 
// revision_date : 
// author(s)     : Herve Bronnimann (Herve.Bronnimann@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================


#include <cassert>

CGAL_BEGIN_NAMESPACE

template <class Ftds>
void
_test_cls_tds_face( const Ftds )
{
  std::cout << "    face" << std::endl;

  typedef typename Ftds::Tds                  Tds;
  typedef typename Ftds::Vertex               Vertex;
  typedef typename Ftds::Face                 Face;
  typedef typename Ftds::Face_handle          Face_handle;
  typedef typename Ftds::Vertex_handle        Vertex_handle;

  // Build a few objects
  int i;
  Vertex v1, v2, v3, v4;
  Vertex_handle vh1 = v1.handle();
  Vertex_handle vh2 = v2.handle();
  Vertex_handle vh3 = v3.handle();
  Vertex_handle vh4 = v4.handle();

  // Test constructors
  Face f1,f2,f3;
  Face_handle fh1 = f1.handle();
  Face_handle fh2 = f2.handle();
  Face_handle fh3 = f3.handle();
  Face f4(vh1,vh2,vh3);
  Face f5(vh1,vh2,vh3,fh1,fh2,fh3);
  Face_handle fh4 = f4.handle();
  
  // Test vertex()
  assert( f4.vertex(0) == vh1 );
  assert( f4.vertex(1) == vh2 );
  assert( f4.vertex(2) == vh3 );
  assert( f5.vertex(0) == vh1 );
  assert( f5.vertex(1) == vh2 );
  assert( f5.vertex(2) == vh3 );
  
  // Test index(Vertex *)
  assert( f4.index(vh1) == 0 );
  assert( f4.index(vh2) == 1 );
  assert( f4.index(vh3) == 2 );
  assert( f5.index(vh1) == 0 );
  assert( f5.index(vh2) == 1 );
  assert( f5.index(vh3) == 2 );
      
  // Test has_vertex()
  assert( f4.has_vertex(vh1) );
  assert( ! f4.has_vertex(vh4) );
  assert( f4.has_vertex(vh2,i) && (i==1) );
  assert( f5.has_vertex(vh3,i) && (i==2) );

  // Test set_vertex()
  f1.set_vertex(0,vh4);
  f1.set_vertex(1,vh3);
  f1.set_vertex(2,vh2);
  f2.set_vertex(0,vh3);
  f2.set_vertex(1,vh4);
  f2.set_vertex(2,vh1);
  f3.set_vertex(0,vh2);
  f3.set_vertex(1,vh1);
  f3.set_vertex(2,vh4);
  assert( f1.vertex(0) == vh4 );
  assert( f1.vertex(1) == vh3 );
  assert( f1.vertex(2) == vh2 );
  assert( f2.vertex(0) == vh3 );
  assert( f2.vertex(1) == vh4 );
  assert( f2.vertex(2) == vh1 );
  assert( f3.vertex(0) == vh2 );
  assert( f3.vertex(1) == vh1 );
  assert( f3.vertex(2) == vh4 );
  
  // Test set_vertices()
  f5.set_vertices();
  assert( f5.vertex(0) == NULL );
  assert( f5.vertex(1) == NULL );
  assert( f5.vertex(2) == NULL);
  f5.set_vertices( vh1, vh2, vh3 );
  assert( f5.vertex(0) == vh1 );
  assert( f5.vertex(1) == vh2 );
  assert( f5.vertex(2) == vh3 );

  // Test neighbor()
  assert( f5.neighbor(0) == fh1 );
  assert( f5.neighbor(1) == fh2 );
  assert( f5.neighbor(2) == fh3 );
  
  // Test index(Face *)
  assert( f5.index(fh1) == 0 );
  assert( f5.index(fh2) == 1 );
  assert( f5.index(fh3) == 2 );
   
  // Test has_neighbor()
  assert( f5.has_neighbor(fh1) );
  assert( ! f5.has_neighbor(fh4) );
  assert( f5.has_neighbor(fh2,i) && (i==1) );
  assert( f5.has_neighbor(fh3,i) && (i==2) );

  // Test set_neighbor()
  f1.set_neighbor(0,fh4);
  f2.set_neighbor(1,fh4);
  f3.set_neighbor(2,fh4);
  assert( f1.neighbor(0) == fh4 );
  assert( f2.neighbor(1) == fh4 );
  assert( f3.neighbor(2) == fh4 );
  
  // Test set_neighbors()
  f5.set_neighbors();
  assert( f5.neighbor(0) == NULL );
  assert( f5.neighbor(1) == NULL );
  assert( f5.neighbor(2) == NULL );
  f5.set_neighbors( fh1, fh2, fh3 );
  assert( f5.neighbor(0) == fh1 );
  assert( f5.neighbor(1) == fh2 );
  assert( f5.neighbor(2) == fh3 );
  f4.set_neighbors( fh1, fh2, fh3 );
  assert( f4.neighbor(0) == fh1 );
  assert( f4.neighbor(1) == fh2 );
  assert( f4.neighbor(2) == fh3 );

  //Test mirror_vertex() mirror_index()
  assert( f4.mirror_vertex(0) == vh4);
  assert( f1.mirror_vertex(0) == vh1);
  assert( f2.mirror_vertex(1) == vh2);
  assert( f3.mirror_vertex(2) == vh3);
  assert( f4.mirror_index(0) == 0);
  assert( f1.mirror_index(0) == 0);
  assert( f2.mirror_index(1) == 1);
  assert( f3.mirror_index(2) == 2);

  // Test is_valid
  assert( f4.is_valid() );
  
  // Test ccw() and cw()
  assert( f1.ccw(0) == 1 );
  assert( f1.ccw(1) == 2 );
  assert( f1.ccw(2) == 0 );
  assert( f1.cw(0) == 2 );
  assert( f1.cw(1) == 0 );
  assert( f1.cw(2) == 1 );

  // Test dimension
  assert(f4.dimension() == 2);

  //Test low dimensional faces
  Face g1(vh2,vh3,NULL);
  Face g2(vh3,vh1, NULL);
  Face g3(vh1,vh2, NULL); 
  Face_handle gh1 = g1.handle();
  Face_handle gh2 = g2.handle();
  Face_handle gh3 = g3.handle();
  g1.set_neighbors( gh2, gh3, NULL);
  g2.set_neighbors( gh3, gh1, NULL);
  g3.set_neighbors( gh1, gh2, NULL);

  assert(g1.dimension() == 1);
  assert (g1.is_valid());

  Face h1(vh1, NULL, NULL);
  Face_handle hh1=h1.handle();
  Face h2(vh2, NULL, NULL, hh1, NULL, NULL);
  Face_handle hh2=h2.handle();
  h1.set_neighbor(0, hh2);
  assert (h1.dimension() == 0);
  assert (h1.is_valid()); 

  assert(g1.mirror_vertex(0) == vh1);
  assert(g1.mirror_vertex(1) == vh1);
  assert(g1.mirror_index(0) == 1);
  assert(g1.mirror_index(1) == 0);

  return;
}

CGAL_END_NAMESPACE
