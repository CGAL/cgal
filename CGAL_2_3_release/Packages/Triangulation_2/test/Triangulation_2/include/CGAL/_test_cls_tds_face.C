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

template <class Face, class Gt>
void
_test_cls_tds_face( const Face &, const Gt & )
{
  std::cout << "    face" << std::endl;

  typedef typename Face::Vertex               Vertex;
  // typedef typename Face::Edge                 Edge;

  // Build a few objects
  int i;
  Vertex v1, v2, v3, v4;
  
  // Test constructors
  Face f1,f2,f3;
  Face f4(&v1,&v2,&v3);
  Face f5(&v1,&v2,&v3,&f1,&f2,&f3);
  
  // Test vertex()
  assert( f4.vertex(0) == &v1 );
  assert( f4.vertex(1) == &v2 );
  assert( f4.vertex(2) == &v3 );
  assert( f5.vertex(0) == &v1 );
  assert( f5.vertex(1) == &v2 );
  assert( f5.vertex(2) == &v3 );
  
  // Test index(Vertex *)
  assert( f4.index(&v1) == 0 );
  assert( f4.index(&v2) == 1 );
  assert( f4.index(&v3) == 2 );
  assert( f5.index(&v1) == 0 );
  assert( f5.index(&v2) == 1 );
  assert( f5.index(&v3) == 2 );
      
  // Test has_vertex()
  assert( f4.has_vertex(&v1) );
  assert( ! f4.has_vertex(&v4) );
  assert( f4.has_vertex(&v2,i) && (i==1) );
  assert( f5.has_vertex(&v3,i) && (i==2) );

  // Test set_vertex()
  f1.set_vertex(0,&v4);
  f1.set_vertex(1,&v3);
  f1.set_vertex(2,&v2);
  f2.set_vertex(0,&v3);
  f2.set_vertex(1,&v4);
  f2.set_vertex(2,&v1);
  f3.set_vertex(0,&v2);
  f3.set_vertex(1,&v1);
  f3.set_vertex(2,&v4);
  assert( f1.vertex(0) == &v4 );
  assert( f1.vertex(1) == &v3 );
  assert( f1.vertex(2) == &v2 );
  assert( f2.vertex(0) == &v3 );
  assert( f2.vertex(1) == &v4 );
  assert( f2.vertex(2) == &v1 );
  assert( f3.vertex(0) == &v2 );
  assert( f3.vertex(1) == &v1 );
  assert( f3.vertex(2) == &v4 );
  
  // Test set_vertices()
  f5.set_vertices();
  assert( f5.vertex(0) == NULL );
  assert( f5.vertex(1) == NULL );
  assert( f5.vertex(2) == NULL);
  f5.set_vertices( &v1, &v2, &v3 );
  assert( f5.vertex(0) == &v1 );
  assert( f5.vertex(1) == &v2 );
  assert( f5.vertex(2) == &v3 );

  // Test neighbor()
  assert( f5.neighbor(0) == &f1 );
  assert( f5.neighbor(1) == &f2 );
  assert( f5.neighbor(2) == &f3 );
  
  // Test index(Face *)
  assert( f5.index(&f1) == 0 );
  assert( f5.index(&f2) == 1 );
  assert( f5.index(&f3) == 2 );
   
  // Test has_neighbor()
  assert( f5.has_neighbor(&f1) );
  assert( ! f5.has_neighbor(&f4) );
  assert( f5.has_neighbor(&f2,i) && (i==1) );
  assert( f5.has_neighbor(&f3,i) && (i==2) );

  // Test set_neighbor()
  f1.set_neighbor(0,&f4);
  f2.set_neighbor(1,&f4);
  f3.set_neighbor(2,&f4);
  assert( f1.neighbor(0) == &f4 );
  assert( f2.neighbor(1) == &f4 );
  assert( f3.neighbor(2) == &f4 );
  
  // Test set_neighbors()
  f5.set_neighbors();
  assert( f5.neighbor(0) == NULL );
  assert( f5.neighbor(1) == NULL );
  assert( f5.neighbor(2) == NULL );
  f5.set_neighbors( &f1, &f2, &f3 );
  assert( f5.neighbor(0) == &f1 );
  assert( f5.neighbor(1) == &f2 );
  assert( f5.neighbor(2) == &f3 );
  f4.set_neighbors( &f1, &f2, &f3 );
  assert( f4.neighbor(0) == &f1 );
  assert( f4.neighbor(1) == &f2 );
  assert( f4.neighbor(2) == &f3 );

  //Test mirror_vertex() mirror_index()
  assert( f4.mirror_vertex(0) == &v4);
  assert( f1.mirror_vertex(0) == &v1);
  assert( f2.mirror_vertex(1) == &v2);
  assert( f3.mirror_vertex(2) == &v3);
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
  Face g1(&v2,&v3,NULL);
  Face g2(&v3,&v1, NULL);
  Face g3(&v1,&v2, NULL);
  g1.set_neighbors( &g2, &g3, NULL);
  g2.set_neighbors( &g3, &g1, NULL);
  g3.set_neighbors( &g1, &g2, NULL);
  assert(g1.dimension() == 1);
  assert (g1.is_valid());

  Face h1(&v1, NULL, NULL);
  Face h2(&v2, NULL, NULL, &h1, NULL, NULL);
  h1.set_neighbor(0, &h2);
  assert (h1.dimension() == 0);
  assert (h1.is_valid()); 

  assert(g1.mirror_vertex(0) == &v1);
  assert(g1.mirror_vertex(1) == &v1);
  assert(g1.mirror_index(0) == 1);
  assert(g1.mirror_index(1) == 0);

  return;
}

CGAL_END_NAMESPACE
