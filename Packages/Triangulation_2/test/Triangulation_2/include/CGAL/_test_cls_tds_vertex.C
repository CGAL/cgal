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

template <class Vtds>
void
_test_cls_tds_vertex( const Vtds&)
{
  std::cout << "    vertex" << std::endl;

  typedef typename Vtds::Face_handle          Face_handle;
  typedef typename Vtds::Vertex_handle        Vertex_handle;

  typedef typename Vtds::Vertex_circulator    Vertex_circulator;
  typedef typename Vtds::Face_circulator      Face_circulator;
  typedef typename Vtds::Edge_circulator      Edge_circulator;

  typedef typename Vtds::Triangulation_data_structure Tds;
  typedef typename Tds::Face                  Face;
  typedef typename Vtds::Vertex               Vertex;

  Face f3;
  Face_handle fh3 = f3.handle();
  
  // Test constructors
  Vertex v1;
  Vertex v2, v3;
  v3.set_face(f3.handle());
  

  // Test face()
  assert( v3.face() == fh3 );
  // to avoid "unused variable warning
  v3.set_face(fh3);
      
  // Test set_face()
  v2.set_face(fh3);
  assert( v2.face() == fh3 );
  
//   // Test ccw() and cw()
//   assert( v1.ccw(0) == 1 );
//   assert( v1.ccw(1) == 2 );
//   assert( v1.ccw(2) == 0 );
//   assert( v1.cw(0) == 2 );
//   assert( v1.cw(1) == 0 );
//   assert( v1.cw(2) == 1 );
  
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
     assert( v1.face() == NULL);
     Face_circulator   fc = v1.incident_faces(); 
     if (fc != 0) fc++;
     Edge_circulator   ec = v1.incident_edges(); 
     if (ec != 0) ec++;
     Vertex_circulator vc = v1.incident_vertices(); 
     if (vc != 0) vc++;
  // Test degree()
     assert (v1.face() == NULL);
     assert( v1.degree() == 0 );
}

CGAL_END_NAMESPACE
