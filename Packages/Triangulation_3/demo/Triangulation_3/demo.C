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

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>

#include <CGAL/Delaunay_triangulation_3.h>

#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/IO/Triangulation_geomview_ostream_3.h>

#include <iostream>
#include <fstream>
#include <unistd.h>
#include <list>

typedef CGAL::Filtered_kernel<CGAL::Simple_cartesian<double> > K;

typedef CGAL::Triangulation_3<K> Triangulation;
typedef CGAL::Delaunay_triangulation_3<K> Delaunay;

typedef Triangulation::Vertex_iterator Vertex_iterator;
typedef Triangulation::Edge_iterator Edge_iterator;
typedef Triangulation::Cell_iterator Cell_iterator;
typedef Triangulation::Facet_iterator Facet_iterator;
typedef Triangulation::Cell_circulator Cell_circulator;

typedef Triangulation::Cell Cell;
typedef Triangulation::Vertex Vertex;
typedef Triangulation::Cell_handle Cell_handle;
typedef Triangulation::Vertex_handle Vertex_handle;
typedef Triangulation::Locate_type Locate_type;

typedef K::Point_3 Point;

////////////////////// 
// VISU GEOMVIEW
////////////////////// 
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
void visu_facet(CGAL::Geomview_stream & os, const TRIANGULATION & T,
	       Cell_handle c, int i)
{
  if ( ! T.is_infinite(c,i) )
    os << T.triangle(c,i);
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

////////////////////// 

int main()
{
  CGAL::Geomview_stream gv(CGAL::Bbox_3(0,0,0, 2, 2, 2));
  gv.set_bg_color(CGAL::Color(0, 200, 200));
  gv.clear();

  Delaunay T;

  std::ifstream iFile("data/points",std::ios::in);

  if (! iFile) {
    std::cout <<"A file named points in directory data" 
              <<" containing points should be provided," << std::endl 
	      <<"see README"<<std::endl;
    return 1;
  }

  std::cout <<"          Reading file data/points" << std::endl ;
  Point nouv;
  while ( iFile >> nouv ) 
    T.insert(nouv);

  T.is_valid(true);

  std::cout <<"          Visualizing vertices and edges" << std::endl;
  visu_vertices(gv,T);
  gv.set_wired(true);
  gv << T;
  gv.set_wired(false);

  sleep(3);

  std::cout <<"          Locating point (1,1,1) :" << std::endl;
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
  if ( lt == Triangulation::OUTSIDE_CONVEX_HULL ) 
    std::cout <<"                     OUTSIDE_CONVEX_HULL" << std::endl;
  if ( lt == Triangulation::OUTSIDE_AFFINE_HULL ) 
    std::cout <<"                     OUTSIDE_AFFINE_HULL" << std::endl;

  sleep(6);

  std::cout <<"          Visualizing T" << std::endl;
  gv.clear();
  std::cout <<"                - facets" << std::endl;
  gv << T;
  std::cout <<"                - edges only" << std::endl;
  gv.set_wired(true);
  gv << T;
  gv.set_wired(false);
  std::cout <<"          You can move one of the" <<std::endl
	    <<"          two triangulations by selecting it"    <<std::endl
	    <<"          in the Geomview targets" <<std::endl;

  char ch;
  std::cout << "Enter any character to quit" << std::endl;
  std::cin >> ch;

  return 1;
}

#endif // if defined(__BORLANDC__) || defined(_MSC_VER)
