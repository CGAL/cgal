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
// file          : include/CGAL/_test_cls_delaunay_triangulation_2.C
// revision      : 
// revision_date : 
// author(s)     : Herve Bronnimann (Herve.Bronnimann@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================

#include <iostream.h>
#include <fstream.h>

#include <pair.h>
#include <list.h>
#include <vector.h>

#include <CGAL/_test_fct_delaunay_duality.C>
#include <CGAL/_test_cls_triangulation_2.C>

template <class Triangulation>
void
CGAL__test_cls_delaunay_triangulation_2( const Triangulation & )
{
  typedef Triangulation                      Cls;

  // We assume the traits class has been tested already
  // actually, any traits is good if it has been tested
  typedef typename Cls::Geom_traits          Gt;

  typedef typename Cls::Point                Point;
  typedef typename Cls::Segment              Segment;
  typedef typename Cls::Triangle             Triangle;

  typedef typename Cls::Distance             Distance;

  typedef typename Cls::Line                 Line;
  typedef typename Cls::Direction            Direction;
  typedef typename Cls::Ray                  Ray;

  typedef typename Cls::Vertex               Vertex;
  typedef typename Cls::Face                 Face;

  typedef typename Cls::Vertex_handle        Vertex_handle;
  typedef typename Cls::Face_handle          Face_handle;

  typedef pair<Face_handle,int>              Edge;

  typedef typename Cls::Vertex_iterator      Vertex_iterator;
  typedef typename Cls::Face_iterator        Face_iterator;
  typedef typename Cls::Edge_iterator        Edge_iterator;

  typedef typename Cls::Vertex_circulator    Vertex_circulator;
  typedef typename Cls::Face_circulator      Face_circulator;
  typedef typename Cls::Edge_circulator      Edge_circulator;
  typedef typename Cls::Line_face_circulator Line_face_circulator;

  typedef typename Cls::Locate_type          Locate_type;
  
  /***********************/
  /***** SUBCLASSES ******/
  CGAL__test_cls_triangulation_2( Cls() );

  // Constructors
  cout << "    constructors(3)" << endl;

  // Build dummy delaunay triangulations, 1- and 2-dimensional
  Cls T1;
  int m,p;
  for (m=0; m<20; m++)
      T1.insert( Point(3*m, 2*m) );
  assert( T1.is_valid() );
   
  Cls T2;
  for (m=0; m<20; m++)
      for (p=0; p<20; p++)
	  T2.insert( Point(3*m+p, m-2*p) );
  assert( T2.is_valid() );
  
  Cls T3;
  // All these points are on a circle of radius 325
  Point pt[28] = {
      Point(36,323), Point(80,315), Point(91,-312),
      Point(125,300), Point(165,280), Point(195,-260), Point(204,253),
      Point(36,-323), Point(80,-315), Point(91,-312),
      Point(125,-300), Point(165,-280), Point(195,-260), Point(204,-253),
      Point(-36,-323), Point(-80,-315), Point(-91,-312),
      Point(-125,-300), Point(-165,-280), Point(-195,-260), Point(-204,-253),
      Point(-36,323), Point(-80,315), Point(-91,312),
      Point(-125,300), Point(-165,280), Point(-195,260), Point(-204,253)
  };
  for (m=0; m<28; m++) 
      T3.insert( Point(pt[m]) );
  assert( T3.is_valid() );
 
  /********************/
  /***** Duality ******/
  cout << "    duality" << endl;
  CGAL__test_delaunay_duality(T1);
  CGAL__test_delaunay_duality(T2);
  CGAL__test_delaunay_duality(T3);
}
