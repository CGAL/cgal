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
// file          : include/CGAL/_test_cls_triangulation_vertex.C
// revision      : 
// revision_date : 
// author(s)     : Herve Bronnimann (Herve.Bronnimann@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================


#include <cassert>

template <class Vertex>
void
_test_cls_triangulation_vertex( const Vertex & )
{
  std::cout << "    vertex" << std::endl;

  typedef typename Vertex::Geom_traits          Gt;

  typedef typename Vertex::Point                Point;
  typedef typename Vertex::Segment              Segment;
  typedef typename Vertex::Triangle             Triangle;

  typedef typename Vertex::Face                 Face;
  typedef typename Vertex::Edge                 Edge;

  typedef typename Vertex::Vertex_handle        Vertex_handle;
  typedef typename Vertex::Face_handle          Face_handle;

  typedef typename Vertex::Vertex_circulator    Vertex_circulator;
  typedef typename Vertex::Face_circulator      Face_circulator;
  typedef typename Vertex::Edge_circulator      Edge_circulator;

  // Build a few objects
  Point p2(5,6,1);
  Point p3(6,14,1);
  Face f3;
  
  // Test constructors
  Vertex v1;
  Vertex v2(p2);
  Vertex v3(p3,f3.handle());
  
  // test is_valid
  assert( v3.face() == f3.handle() );
  f3.set_vertex(1,v3.handle());
  assert(v3.is_valid());

  // Test point()
  assert( Gt().compare_x_2_object()( v2.point(), p2) == CGAL::EQUAL &&
	  Gt().compare_y_2_object()( v2.point(), p2) == CGAL::EQUAL);
  assert( Gt().compare_x_2_object()( v3.point(), p3) == CGAL::EQUAL &&
	  Gt().compare_y_2_object()( v3.point(), p3) == CGAL::EQUAL);
    
  // Test face()
  assert( v3.face() == f3.handle() );
      
  // Test set_face()
  v2.set_face(f3.handle());
  assert( v2.face() == f3.handle() );
  
  // Test ccw() and cw()
  assert( v1.ccw(0) == 1 );
  assert( v1.ccw(1) == 2 );
  assert( v1.ccw(2) == 0 );
  assert( v1.cw(0) == 2 );
  assert( v1.cw(1) == 0 );
  assert( v1.cw(2) == 1 );
  
  // The functions degree(), incident_faces(), incident_vertices(),
  // incident_edges() need a vertex in some triangulation,
  // They are tested in _test_cls_triangulation_2.h
  
  // The following are here to test the presence of the types
  // otherwise the compiler might not instantiate them
  // the asert() is to avoid an unused variable warning
  Face_circulator   fc = v1.incident_faces(); assert( &fc == &fc );
  Edge_circulator   ec = v1.incident_edges(); assert( &ec == &ec );
  Vertex_circulator vc = v1.incident_vertices(); assert( &vc == &vc );
}
