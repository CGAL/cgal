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
// file          : demo/Triangulation3/blind_demo.C
// revision      : $Revision$
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis (Mariette Yvinec)
//
// ============================================================================
#include <CGAL/basic.h>

#include <cstring>
#include <iostream>
#include <fstream>
#include <strstream.h>

#include <list>

#include <CGAL/Cartesian.h>
#include <CGAL/Point_3.h>

#include <CGAL/Triangulation_iterators_3.h>
#include <CGAL/Triangulation_circulators_3.h>
#include <CGAL/Triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulation_geom_traits_3.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Delaunay_triangulation_3.h>

typedef CGAL::Cartesian<double>  Rep;

typedef CGAL::Triangulation_geom_traits_3<Rep> Gt;
typedef CGAL::Triangulation_vertex_base_3<Gt> Vb;
typedef CGAL::Triangulation_cell_base_3<Gt>  Cb;

typedef CGAL::Triangulation_data_structure_3<Vb,Cb> TDS;
typedef CGAL::Triangulation_3<Gt,TDS> Triangulation;
typedef CGAL::Delaunay_triangulation_3<Gt,TDS> Delaunay;

typedef CGAL::Triangulation_vertex_iterator_3<Gt,TDS> Vertex_iterator;
typedef CGAL::Triangulation_edge_iterator_3<Gt,TDS> Edge_iterator;
typedef CGAL::Triangulation_cell_iterator_3<Gt,TDS> Cell_iterator;
typedef CGAL::Triangulation_facet_iterator_3<Gt,TDS> Facet_iterator;
typedef CGAL::Triangulation_cell_circulator_3<Gt,TDS> Cell_circulator;

typedef typename Triangulation::Cell Cell;
typedef typename Triangulation::Vertex Vertex;
typedef typename Triangulation::Cell_handle Cell_handle;
typedef typename Triangulation::Vertex_handle Vertex_handle;
typedef typename Triangulation::Locate_type Locate_type;

typedef Gt::Point Point;
//typedef CGAL::Point_3<Rep>  Point;

int main(int argc, char* argv[])
{

  Delaunay T;

  ifstream iFile("data",ios::in);

  if (! iFile) {
    cout <<"A file named data containing points should be provided," << endl
	 <<"see README"<<endl;
    return 1;
  }

  cout <<"                              reading file data" << endl ;
  Point nouv;
  while ( iFile >> nouv ) {
    T.insert(nouv);
  }

  T.is_valid(true);

  ofstream oFileT("output",ios::out);
  cout <<"                              writing file output" << endl << flush;
  oFileT << T;

  char ch;
  cout << "enter any character to quit" << endl;
  cin >> ch;

  return 1;
}
