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
// file          : include/CGAL/_test_cls_tds_2.C
// revision      : 
// revision_date : 
// author(s)     : Herve Bronnimann (Herve.Bronnimann@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================


#include <cassert>

#include <CGAL/_test_cls_tds_vertex.C>
#include <CGAL/_test_cls_tds_face.C>

template <class Tds, class Gt>
void
CGAL::_test_cls_tds_2( const Tds &, const Gt & )
{
  // Since tds is extensively used by Triangulation_2,
  // there is no real need to test the functionality
  // we simply test for the presence of the types and the
  // functions

  typedef typename Tds::Vertex            Vertex;
  typedef typename Tds::Face              Face;
  // missing in the documentation, it's a bug
  // typedef typename Tds::Edge              Edge;

  typedef typename Tds::Vertex_iterator   Vertex_iterator;
  typedef typename Tds::Face_iterator     Face_iterator;
  typedef typename Tds::Edge_iterator     Edge_iterator;
  
  typedef typename Tds::Vertex_circulator Vertex_circulator;
  typedef typename Tds::Face_circulator   Face_circulator;
  typedef typename Tds::Edge_circulator   Edge_circulator;

  // Test subclasses
  CGAL::_test_cls_tds_vertex( Vertex(), Gt() );
  CGAL::_test_cls_tds_face( Face(), Gt() );

  // Test constructors
  cout << "    constructors" << endl;
  Tds tds1;
  Tds tds2(new Vertex);
  Tds tds3(tds2);
  Tds tds4 = tds2;
  tds4.swap(tds2);

  // Setting functions
  cout << "    setting functions" << endl;
  Face   *f1 = new Face;
  Vertex *vt1 = new Vertex; vt1->set_face(f1);
  f1->set_vertices(NULL,vt1,NULL);
  tds1.set_finite_vertex(vt1);
  tds1.set_number_of_vertices(1);
  assert( tds1.number_of_vertices() == 1 );
  Face   *f2 = new Face;
  Vertex *vt2 = new Vertex; vt2->set_face(f2);
  f2->set_vertices(NULL,NULL,vt2);
  tds1.set_infinite_vertex(vt2);

  // Finite and infinite vertices and faces
  cout << "    finite/infinite faces and vertices" << endl;
  assert( !tds1.is_infinite(vt1) );
  assert( tds1.is_infinite(vt2) );
  assert( !tds1.is_infinite(f1) );
  assert( tds1.is_infinite(f2) );
  assert( tds1.is_infinite(f2,1) );

  // assert( tds1.infinite_face() == f2 );
  assert( tds1.infinite_vertex() == vt2 );
  assert( tds1.finite_vertex() == vt1 );

  // The other functions are tested in CGAL::_test_cls_triangulation_2()
  cout << "    insert... and remove... are tested by test_triangulation_2" << endl;
  cout << "    iterators and circulators are tested by test_triangulation_2" << endl;
  
  // misc.
  cout << "    miscellaneous" << endl;
  assert( tds1.ccw(0) == 1 );
  assert( tds1.ccw(1) == 2 );
  assert( tds1.ccw(2) == 0 );
  assert( tds1.cw(0) == 2 );
  assert( tds1.cw(1) == 0 );
  assert( tds1.cw(2) == 1 );
  
  // make tds1 and tds2 valid in order to allow clear() to work
  cout << "    resetting all tds for clear to work" << endl;
  tds1.set_number_of_vertices(0);
  tds1.set_infinite_vertex(vt2); vt2->set_face(NULL);
  tds1.set_finite_vertex(NULL);
  assert( tds1.is_valid() );

  // THE FOLLOWING IS NOT NEEDED AND CAUSE A SEGMENTATION FAULT
  // ON LINUX WHEN LEAVING THE PROCEDURE
//   tds2.set_number_of_vertices(0);
//   tds2.set_infinite_vertex(vt2); vt2->set_face(NULL);
//   tds2.set_finite_vertex(NULL);
//   assert( tds2.is_valid() );
// 
//   tds3.set_number_of_vertices(0);
//   tds3.set_infinite_vertex(vt2); vt2->set_face(NULL);
//   tds3.set_finite_vertex(NULL);
//   assert( tds3.is_valid() );
// 
//   tds4.set_number_of_vertices(0);
//   tds4.set_infinite_vertex(vt2); vt2->set_face(NULL);
//   tds4.set_finite_vertex(NULL);
//   assert( tds4.is_valid() );

  // test destructor and return
  cout << "    test destructors and return" << endl;
  tds1.clear();
  tds2.clear();
  tds3.clear();
  tds4.clear(); 
    
}
