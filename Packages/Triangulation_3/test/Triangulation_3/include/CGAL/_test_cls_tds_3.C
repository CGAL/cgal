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
// file          : include/CGAL/_test_cls_tds_3.C
// revision      : 
// revision_date : 
// author(s)     : Herve Bronnimann (Herve.Bronnimann@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================
#include <cassert>
#include "_test_cls_tds_vertex.C"
#include "_test_cls_tds_cell.C"

template <class Tds>
void
_test_cls_tds_3( const Tds &)
{

  typedef typename Tds::Vertex            Vertex;
  typedef typename Tds::Cell              Cell;
  typedef typename Tds::Edge              Edge;
  typedef typename Tds::Facet             Facet;

  typedef typename Tds::Vertex_iterator   Vertex_iterator;
  typedef typename Tds::Facet_iterator    Facet_iterator;
  typedef typename Tds::Edge_iterator     Edge_iterator;
  typedef typename Tds::Cell_iterator     Cell_iterator;


  // test Vertex and cell :
  cout << "    Test Vertex " << endl;
  _test_vertex_tds_3(Vertex());

  cout << "    Test Cell " << endl;
  
  _test_cell_tds_3(Cell());


  cout << "   Testing TDS " << endl;
  
  // Test constructors
  cout << "    constructors" << endl;
  Tds tds1;
  Tds tds2;
  Vertex_iterator vit;
  Vertex* v1 = new Vertex;
  tds2.insert_outside_affine_hull(*v1, NULL);
  Tds tds3(tds2);
  Vertex* v2 = new Vertex;

  vit=tds3.vertices_begin();
  tds3.insert_outside_affine_hull(*v2,&*vit);
  cout << "ok" << endl;
  Tds tds4 = tds3;
  Vertex* v3 = new Vertex;
  vit=tds4.vertices_begin();
  tds4.insert_outside_affine_hull(*v3,&*vit);
  cout << "ok" << endl;
  Tds tds5;
  tds5.swap(tds4);
  tds4=tds5;
  Vertex* v4 = new Vertex;
  vit=tds5.vertices_begin();
  tds5.insert_outside_affine_hull(*v4,&*vit);
  cout << "ok" << endl;
  Tds tds6;
  tds6.swap(tds5);
  tds5=tds6;
  Vertex* v5 = new Vertex;
  vit=tds6.vertices_begin();
  tds6.insert_outside_affine_hull(*v5,&*vit);
  cout << "ok" << endl;

  // Setting functions
  cout << "    setting functions" << endl;
  tds1.set_number_of_vertices(1);
  tds2.set_number_of_vertices(5);

  
  cout << "  Insert are tested in test_triangulation_3  " << endl;

  cout << "  Iterator and circulator are tested in test_triangulation_3  " << endl;

  // Access functions

  assert(tds1.dimension()==-2);
  assert(tds2.dimension()==-1);
  assert(tds3.dimension()==0);
  assert(tds4.dimension()==1);
  assert(tds5.dimension()==2);
  assert(tds6.dimension()==3);

  assert(tds3.number_of_vertices()==2);
  assert(tds2.number_of_vertices()==5);
  assert(tds1.number_of_vertices()==1);


  // test destructor and return
  cout << "    test destructors and return" << endl;
  tds1.clear();
  tds2.clear();
  tds3.clear();
}
