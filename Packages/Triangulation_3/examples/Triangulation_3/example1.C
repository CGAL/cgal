// ============================================================================
//
// Copyright (c) 1998-1999 The CGAL Consortium
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
// file          : example/Triangulation3/example1.C
// revision      : $Revision$
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis (Mariette Yvinec)
//
// ============================================================================
#include <CGAL/basic.h>

#include <iostream>
#include <fstream>

#include <assert.h>
#include <list>
#include <vector>

#include <CGAL/Cartesian.h>

#include <CGAL/Triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulation_geom_traits_3.h>
#include <CGAL/Triangulation_3.h>

// for this simple example, using doubles is ok
//
typedef double NT;

// for more complicated examples with degenerate configurations,
// using Filtered_exact number type is advised :
// 
// #include <CGAL/Arithmetic_filter.h>
// #include <CGAL/MP_Float.h>
// 
// typedef CGAL::Filtered_exact<double, CGAL::MP_Float> NT;

typedef CGAL::Cartesian<NT> Repr;
typedef CGAL::Triangulation_geom_traits_3<Repr> Gt;

typedef CGAL::Triangulation_vertex_base_pointer_3<Gt> Vb;
typedef CGAL::Triangulation_cell_base_3<Gt>  Cb;

typedef CGAL::Triangulation_data_structure_3<Vb,Cb> TDS;
typedef CGAL::Triangulation_3<Gt,TDS> Triangulation;

typedef Triangulation::Cell_handle Cell_handle;
typedef Triangulation::Vertex_handle Vertex_handle;
typedef Triangulation::Locate_type Locate_type;

typedef Gt::Point_3 Point;

int main()
{
  Triangulation T;

  // insertion from a list :
  std::list<Point> L;
  L.push_front(Point(0,0,0));
  L.push_front(Point(1,0,0));
  L.push_front(Point(0,1,0));

  int n = T.insert(L.begin(), L.end());

  // insertion from a vector :
  std::vector<Point> V(3);
  V[0] = Point(0,0,1);
  V[1] = Point(1,1,1);
  V[2] = Point(2,2,2);

  n = n + T.insert(V.begin(), V.end());

  // 6 points have been inserted :
  assert( n == 6 );

  // checking validity of T :
  assert( T.is_valid(false) );

  Locate_type lt;
  int li, lj;
  Point p(0,0,0);
  Cell_handle c = T.locate(p, lt, li, lj);
  // p is the vertex of c of index li :
  assert( lt == Triangulation::VERTEX );
  assert(  c->vertex(li)->point() == p );

  Vertex_handle v = c->vertex( (li+1)&3 );
  // v is another vertex of c
  Cell_handle nc = c->neighbor(li);
  // nc = neighbor of c opposite to the vertex associated with p
  // nc must have vertex v :
  int nli;
  assert(  nc->has_vertex( v, nli ) );
  // nli is the index of v in nc

  std::ofstream oFileT("output",std::ios::out);
  // writing file output; 
  // this file is meant to be read only by the operator >>
  oFileT << T; 

  Triangulation T1;
  std::ifstream iFileT("output",std::ios::in);
  // reading file output; 
  iFileT >> T1; 
  assert( T1.is_valid() );
  assert( T1.number_of_vertices() == T.number_of_vertices() );
  assert( T1.number_of_cells() == T.number_of_cells() );

  return 0;
}
