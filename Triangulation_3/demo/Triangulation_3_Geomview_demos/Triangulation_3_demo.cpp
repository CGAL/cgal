// Copyright (c) 1999, 2000, 2001, 2002, 2003, 2004, 2005  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>

#include <CGAL/basic.h>

#ifndef CGAL_USE_GEOMVIEW
#include <iostream>
int main()
{
  std::cerr << "Geomview doesn't work on this platform,"
               " so this demo doesn't work" << std::endl;
  return 0;
}
#else

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_3.h>

#include <CGAL/IO/Triangulation_geomview_ostream_3.h>

#include <iostream>
#include <fstream>
#include <iterator>
#include <unistd.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Delaunay_triangulation_3<K>  Triangulation;

typedef Triangulation::Finite_vertices_iterator Finite_vertices_iterator;
typedef Triangulation::Cell_handle              Cell_handle;
typedef Triangulation::Locate_type              Locate_type;
typedef Triangulation::Point                    Point;

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
  Finite_vertices_iterator vit = T.finite_vertices_begin();
  Finite_vertices_iterator vdone = T.finite_vertices_end();

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

  Triangulation  T;

  std::ifstream iFile("data/points",std::ios::in);

  if (! iFile) {
    std::cout <<"A file named points in directory data"
              <<" containing points should be provided," << std::endl
              <<"see README"<<std::endl;
    return 1;
  }

  std::cout <<"          Reading file data/points" << std::endl ;
  {
      std::istream_iterator<Point> begin (iFile), end;
      T.insert (begin, end);
  }

  T.is_valid(true);

  std::cout <<"          Visualizing vertices and edges" << std::endl;
  visu_vertices(gv,T);
  gv.set_wired(true);
  gv << T;
  gv.set_wired(false);

  sleep(3);

  std::cout <<"          Locating point (1,1,1) :" << std::endl;
  Point p(1,1,1);
  gv.set_vertex_color(CGAL::orange());
  gv << p;
  Locate_type lt;
  int li, lj;
  Cell_handle c = T.locate(p,lt,li,lj);

  sleep(3);

  gv << CGAL::violet();
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

  return 0;
}

#endif // CGAL_USE_GEOMVIEW
