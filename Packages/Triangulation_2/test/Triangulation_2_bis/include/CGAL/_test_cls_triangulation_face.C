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
// file          : include/CGAL/_test_cls_triangulation_face.C
// revision      : 
// revision_date : 
// author(s)     : Herve Bronnimann (Herve.Bronnimann@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================


#include <cassert>

template <class Face>
void
_test_cls_triangulation_face( const Face & )
{
  std::cout << "    face" << std::endl;
  typedef typename Face::Geom_traits          Gt;

  typedef typename Face::Point                Point;
  typedef typename Face::Segment              Segment;
  typedef typename Face::Triangle             Triangle;

  typedef typename Face::Vertex               Vertex;
  // typedef typename Face::Edge                 Edge;

  typedef typename Face::Vertex_handle        Vertex_handle;
  typedef typename Face::Face_handle          Face_handle;

  // Build a few objects
  int i;
  Point p1(5,6,1);
  Point p2(6,14,1);
  Vertex v1, v2, v3(p1), v4(p2);
  
  // Test constructors
  Face f1,f2,f3;
  Face f4(v1.handle(),v2.handle(),v3.handle());
  Face f5(v1.handle(),v2.handle(),v3.handle(),f1.handle(),f2.handle(),f3.handle());
  
  // Test vertex()
  assert( f4.vertex(0) == v1.handle() );
  assert( f4.vertex(1) == v2.handle() );
  assert( f4.vertex(2) == v3.handle() );
  assert( f5.vertex(0) == v1.handle() );
  assert( f5.vertex(1) == v2.handle() );
  assert( f5.vertex(2) == v3.handle() );
  
  // Test index()
  assert( f4.index(v1.handle()) == 0 );
  assert( f4.index(v2.handle()) == 1 );
  assert( f4.index(v3.handle()) == 2 );
  assert( f5.index(v1.handle()) == 0 );
  assert( f5.index(v2.handle()) == 1 );
  assert( f5.index(v3.handle()) == 2 );
      
  // Test has_vertex()
  assert( f4.has_vertex(v1.handle()) );
  assert( ! f4.has_vertex(v4.handle()) );
  assert( f4.has_vertex(v2.handle(),i) && (i==1) );
  assert( f5.has_vertex(v3.handle(),i) && (i==2) );

  // Test set_vertex()
  f1.set_vertex(1,v3.handle());
  f1.set_vertex(2,v2.handle());
  f2.set_vertex(2,v1.handle());
  f2.set_vertex(0,v3.handle());
  f3.set_vertex(0,v2.handle());
  f3.set_vertex(1,v1.handle());
  assert( f1.vertex(1) == v3.handle() );
  assert( f1.vertex(2) == v2.handle() );
  assert( f2.vertex(0) == v3.handle() );
  assert( f2.vertex(2) == v1.handle() );
  assert( f3.vertex(0) == v2.handle() );
  assert( f3.vertex(1) == v1.handle() );
  
  // Test set_vertices()
  f5.set_vertices();
  assert( f5.vertex(0).is_null() );
  assert( f5.vertex(1).is_null() );
  assert( f5.vertex(2).is_null() );
  f5.set_vertices( v1.handle(), v2.handle(), v3.handle() );
  assert( f5.vertex(0) == v1.handle() );
  assert( f5.vertex(1) == v2.handle() );
  assert( f5.vertex(2) == v3.handle() );

  // Test neighbor()
  assert( f5.neighbor(0) == f1.handle() );
  assert( f5.neighbor(1) == f2.handle() );
  assert( f5.neighbor(2) == f3.handle() );
  
  // Test index()
  assert( f5.index(f1.handle()) == 0 );
  assert( f5.index(f2.handle()) == 1 );
  assert( f5.index(f3.handle()) == 2 );
   
  // Test has_neighbor()
  assert( f5.has_neighbor(f1.handle()) );
  assert( ! f5.has_neighbor(f4.handle()) );
  assert( f5.has_neighbor(f2.handle(),i) && (i==1) );
  assert( f5.has_neighbor(f3.handle(),i) && (i==2) );

  // Test set_neighbor()
  f1.set_neighbor(0,f4.handle());
  f2.set_neighbor(1,f4.handle());
  f3.set_neighbor(2,f4.handle());
  assert( f1.neighbor(0) == f4.handle() );
  assert( f2.neighbor(1) == f4.handle() );
  assert( f3.neighbor(2) == f4.handle() );
  
  // Test set_neighbors()
  f5.set_neighbors();
  assert( f5.neighbor(0).is_null() );
  assert( f5.neighbor(1).is_null() );
  assert( f5.neighbor(2).is_null() );
  f5.set_neighbors( f1.handle(), f2.handle(), f3.handle() );
  assert( f5.neighbor(0) == f1.handle() );
  assert( f5.neighbor(1) == f2.handle() );
  assert( f5.neighbor(2) == f3.handle() );
  f4.set_neighbors( f1.handle(), f2.handle(), f3.handle() );
  assert( f4.neighbor(0) == f1.handle() );
  assert( f4.neighbor(1) == f2.handle() );
  assert( f4.neighbor(2) == f3.handle() );

  // Test is_valid
  assert( f4.is_valid() );
  //assert( !f5.is_valid() );

  // Test ccw() and cw()
  assert( f1.ccw(0) == 1 );
  assert( f1.ccw(1) == 2 );
  assert( f1.ccw(2) == 0 );
  assert( f1.cw(0) == 2 );
  assert( f1.cw(1) == 0 );
  assert( f1.cw(2) == 1 );
}
