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
// file          : demo/Triangulation3/demo.C
// revision      : $Revision$
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis (Mariette Yvinec)
//
// ============================================================================

// Geomview doesn't work on M$ at the moment, so we don't compile this file.
#if defined(__BORLANDC__) || defined(_MSC_VER)
#include <iostream>
int main()
{
  std::cerr << "Geomview doesn't work on Windows, so this demo doesn't work"
            << std::endl;
  return 0;
}
#else

#include <CGAL/basic.h>

#include <iostream>
#include <fstream>
#include <unistd.h>

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

#include <CGAL/IO/Geomview_stream.h>

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

////////////////////// 
// VISU GEOMVIEW
////////////////////// 
template<class TRIANGULATION>
void visu_cells(CGAL::Geomview_stream & os, const TRIANGULATION & T)
{
  Cell_iterator cit = T.finite_cells_begin();
  Cell_iterator cdone = T.cells_end();
  
  if ( cit == cdone ) { std::cout << "no cell" << std::endl ;}
  else {
    while(cit != cdone) {
      os << T.tetrahedron(&(*cit));
    ++cit;
    }
  }
}
template<class TRIANGULATION>
void visu_cell(CGAL::Geomview_stream & os, const TRIANGULATION & T,
	       Cell_handle c)
{
  if ( ! T.is_infinite(c) )
    os << T.tetrahedron(c);
  else 
    os << T.triangle(c,c->index(T.infinite_vertex()));
}
template<class TRIANGULATION>
void visu_facets(CGAL::Geomview_stream & os, const TRIANGULATION & T)
{
  Facet_iterator fit = T.finite_facets_begin();
  Facet_iterator fdone = T.facets_end();
  
  if ( fit == fdone ) { std::cout << "no facet" << std::endl ;}
  else {
    while(fit != fdone) {
      os << T.triangle(*fit);
      ++fit;
    }
  }
}
template<class TRIANGULATION>
void visu_facet(CGAL::Geomview_stream & os, const TRIANGULATION & T,
	       Cell_handle c, int i)
{
  if ( ! T.is_infinite(c,i) )
    os << T.triangle(c,i);
}
template<class TRIANGULATION>
void visu_edges(CGAL::Geomview_stream & os, const TRIANGULATION & T)
{
  Edge_iterator eit = T.finite_edges_begin();
  Edge_iterator edone = T.edges_end();
  
  if ( eit == edone ) { std::cout << "no edge" << std::endl ;}
  else {
    while(eit != edone) {
      os << T.segment(*eit);
      ++eit;
    }
  }
}
template<class TRIANGULATION>
void visu_edge(CGAL::Geomview_stream & os, const TRIANGULATION & T,
	       Cell_handle c, int i, int j)
{
  if ( ! T.is_infinite(c,i,j) )
    os << T.segment(c,i,j);
}
template<class TRIANGULATION>
void visu_vertices(CGAL::Geomview_stream & os, const TRIANGULATION & T)
{
  Vertex_iterator vit = T.finite_vertices_begin();
  Vertex_iterator vdone = T.vertices_end();
  
  if ( vit == vdone ) { std::cout << "no vertex" << std::endl ;}
  else {
    while(vit != vdone) {
      os << vit->point();
      ++vit;
    }
  }
}
template<class TRIANGULATION>
void visu_vertex(CGAL::Geomview_stream & os, const TRIANGULATION & T,
	       Cell_handle c, int i)
{
  if ( ! T.is_infinite(c->vertex(i)) )
    os << c->vertex(i)->point();
}

int main()
{
  CGAL::Geomview_stream gv(CGAL::Bbox_3(0,0,0, 2, 2, 2));

  gv.set_line_width(4);
  gv.set_trace(false);
  gv.set_bg_color(CGAL::Color(0, 200, 200));

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

  std::cout <<"                visualizing vertices and edges" << std::endl;
  gv << CGAL::GREEN;
  visu_vertices(gv,T);
  visu_edges(gv,T);

  sleep(3);

  std::cout <<"                locating point (1,1,1) :" << std::endl;
  Point p(1,1,1);
  gv.set_vertex_color(CGAL::ORANGE);
  gv << p;
  Locate_type lt;
  int li, lj;
  Cell_handle c = T.locate(p,lt,li,lj);

  sleep(3);

  gv << CGAL::VIOLET;
  if ( lt == Triangulation::CELL ) {
    std::cout <<"                     CELL" << std::endl;
    visu_cell(gv,T,c);
  }
  if ( lt == Triangulation::FACET ) {
    std::cout <<"                     FACET" << std::endl;
    visu_facet(gv,T,c,li);
  }
  if ( lt == Triangulation::EDGE ) {
    std::cout <<"                     EDGE" << std::endl;
    visu_edge(gv,T,c,li,lj);
  }
  if ( lt == Triangulation::VERTEX ) {
    std::cout <<"                     VERTEX" << std::endl;
    visu_vertex(gv,T,c,li);
  }
  if ( lt == Triangulation::OUTSIDE_CONVEX_HULL ) {
    std::cout <<"                     OUTSIDE_CONVEX_HULL" << std::endl;
  }
  if ( lt == Triangulation::OUTSIDE_AFFINE_HULL ) {
    std::cout <<"                     OUTSIDE_AFFINE_HULL" << std::endl;
  }

  sleep(6);

  gv.clear();
  gv.set_face_color(CGAL::BLUE);
  gv.set_edge_color(CGAL::GREEN);
  gv.set_vertex_color(CGAL::RED);

  std::cout <<"                visualizing T" << std::endl;
  visu_cells(gv,T);
  visu_vertices(gv,T);
  visu_edges(gv,T);

  std::ofstream oFileT("data/output",std::ios::out);
  std::cout <<"                writing file data/output" << std::endl ;
  oFileT << T;


  char ch;
  std::cout << "enter any character to quit" << std::endl;
  std::cin >> ch;

  return 1;
}

#endif // if defined(__BORLANDC__) || defined(_MSC_VER)
