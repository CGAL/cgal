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
// file          : include/CGAL/_test_cls_tds_vertex.C
// revision      : 
// revision_date : 
// author(s)     : Herve Bronnimann (Herve.Bronnimann@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================


#include <cassert>
CGAL_BEGIN_NAMESPACE

template <class Vertex, class Gt>
void
_test_cls_tds_vertex( const Vertex &, const Gt & )
{
  std::cout << "    vertex" << std::endl;

  typedef typename Vertex::Point                Point;
  typedef typename Vertex::Face                 Face;
  typedef typename Vertex::Edge                 Edge;

  typedef typename Vertex::Vertex_circulator    Vertex_circulator;
  typedef typename Vertex::Face_circulator      Face_circulator;
  typedef typename Vertex::Edge_circulator      Edge_circulator;


  // Build a few objects
  // Build a few objects
  Point p2(1,2);
  Point p3(2,3);
  Face f3;
  
  // Test constructors
  Vertex v1;
  Vertex v2(p2);
  Vertex v3(p3,&f3);
  
  // Test point()
  assert( Gt().compare_x_2_object()(v2.point(),p2) == CGAL::EQUAL &&
	  Gt().compare_y_2_object()(v2.point(),p2) == CGAL::EQUAL) ;
  assert( Gt().compare_x_2_object()(v3.point(),p3) == CGAL::EQUAL &&
	  Gt().compare_y_2_object()(v3.point(),p3) == CGAL::EQUAL) ;
  
  // Test face()
  assert( v3.face() == &f3 );
  // to avoid "unused variable warning
  v3.set_face(&f3);
      
  // Test set_face()
  v2.set_face(&f3);
  assert( v2.face() == &f3 );
  
  // Test set_point()
  v1.set_point(p3);
  assert( Gt().compare_x_2_object()(v1.point(),p3) == CGAL::EQUAL &&
	  Gt().compare_y_2_object()(v1.point(),p3) == CGAL::EQUAL) ;
  
  // Test ccw() and cw()
  assert( v1.ccw(0) == 1 );
  assert( v1.ccw(1) == 2 );
  assert( v1.ccw(2) == 0 );
  assert( v1.cw(0) == 2 );
  assert( v1.cw(1) == 0 );
  assert( v1.cw(2) == 1 );
  
  // The functions degree(), incident_faces(), incident_vertices(),
  // incident_edges() and is_valid() need a vertex in some
  // triangulation,
  // idem  for circulators
  // They are tested in _test_cls_triangulation_2.h
  
  // The following are here to test the presence of the types
  // otherwise the compiler might not instantiate them
  // v1 is chosen because its face() is NULL, 
  // the constructors calls function that crashes with those
  // incomplete data
     Face_circulator   fc = v1.incident_faces(); 
     if (fc != 0) fc++;
     Edge_circulator   ec = v1.incident_edges(); 
     if (ec != 0) ec++;
     Vertex_circulator vc = v1.incident_vertices(); 
     if (vc != 0) vc++;
  // Test degree()
  assert( v1.degree() == 0 );
}

CGAL_END_NAMESPACE
