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

#include <iostream>
#include <fstream>

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

typedef Triangulation::Cell Cell;
typedef Triangulation::Vertex Vertex;
typedef Triangulation::Cell_handle Cell_handle;
typedef Triangulation::Vertex_handle Vertex_handle;
typedef Triangulation::Locate_type Locate_type;

typedef Gt::Point_3 Point;
//typedef CGAL::Point_3<Rep>  Point;

int main()
{

  Delaunay T;

  std::ifstream iFile("data/points",std::ios::in);

  if (! iFile) {
    std::cout <<"A file named points in directory data"
	      <<" containing points should be provided," << std::endl
	      <<"see README"<<std::endl;
    return 1;
  }

  std::cout <<"                reading file data/points" << std::endl ;
  Point nouv;
  while ( iFile >> nouv ) {
    T.insert(nouv);
  }

  T.is_valid(true);

  std::cout <<"                locating point (1,1,1) :" << std::endl;
  Point p(1,1,1);
  Locate_type lt;
  int li, lj;
  T.locate(p,lt,li,lj);
  if ( lt == Triangulation::CELL ) 
    std::cout <<"                     CELL" << std::endl;
  if ( lt == Triangulation::FACET )
    std::cout <<"                     FACET" << std::endl;
  if ( lt == Triangulation::EDGE ) 
    std::cout <<"                     EDGE" << std::endl;
  if ( lt == Triangulation::VERTEX )
    std::cout <<"                     VERTEX" << std::endl;
  if ( lt == Triangulation::OUTSIDE_CONVEX_HULL ) 
    std::cout <<"                     OUTSIDE_CONVEX_HULL" << std::endl;
  if ( lt == Triangulation::OUTSIDE_AFFINE_HULL ) 
    std::cout <<"                     OUTSIDE_AFFINE_HULL" << std::endl;

  std::ofstream oFileT("data/output",std::ios::out);
  std::cout <<"                writing file data/output" << std::endl 
	    << std::flush;
  oFileT << T;

  return 1;
}
