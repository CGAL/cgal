// ======================================================================
//
// Copyright (c) 1999 The GALIA Consortium
//
// This software and related documentation is part of the
// Computational Geometry Algorithms Library (CGAL).
//
// Every use of CGAL requires a license. Licenses come in three kinds:
//
// - For academic research and teaching purposes, permission to use and
//   copy the software and its documentation is hereby granted free of  
//   charge, provided that
//   (1) it is not a component of a commercial product, and
//   (2) this notice appears in all copies of the software and
//       related documentation.
// - Development licenses grant access to the source code of the library 
//   to develop programs. These programs may be sold to other parties as 
//   executable code. To obtain a development license, please contact
//   the GALIA Consortium (at cgal@cs.uu.nl).
// - Commercialization licenses grant access to the source code and the
//   right to sell development licenses. To obtain a commercialization 
//   license, please contact the GALIA Consortium (at cgal@cs.uu.nl).
//
// This software and documentation is provided "as-is" and without
// warranty of any kind. In no event shall the CGAL Consortium be
// liable for any damage of any kind.
//
// The GALIA Consortium consists of Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Free University of Berlin (Germany),
// INRIA Sophia-Antipolis (France), Trier University
// (Germany), Max-Planck-Institute Saarbrucken (Germany),
// and Tel-Aviv University (Israel).
//
// ----------------------------------------------------------------------
//
// file          : demo/GeoWin/CGAL_alpha_shape_2.C
//
// ======================================================================

#include <CGAL/basic.h>

#if !defined(CGAL_USE_LEDA) || (__LEDA__ < 400)
#include <iostream>

int main(int argc, char *argv[])
{
 std::cout << "No LEDA 4.0 or higher installed!\n";
 std::cout << "A LEDA version >= 4.0 is required to run GeoWin!\n";
 return 0;
}
#else 

#define CGAL_ALPHA_WINDOW_STREAM

#include <CGAL/Cartesian.h>
#include <CGAL/squared_distance_2.h> 
#include <CGAL/Point_2.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Alpha_shape_euclidean_traits_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Alpha_shape_2.h>


typedef double coord_type;
typedef CGAL::Cartesian<coord_type>  CRep;
typedef CGAL::Point_2<CRep>  Point;
typedef CGAL::Segment_2<CRep>  Segment;
typedef CGAL::Line_2<CRep>  Line;
typedef CGAL::Triangle_2<CRep>  Triangle;

typedef CGAL::Alpha_shape_euclidean_traits_2<CRep> Gt;

typedef CGAL::Alpha_shape_vertex_base_2<Gt> Vb;

typedef CGAL::Triangulation_face_base_2<Gt> Df;
typedef CGAL::Alpha_shape_face_base_2<Gt, Df>  Fb;

typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;
typedef CGAL::Delaunay_triangulation_2<Gt,Tds> Triangulation_2;

typedef CGAL::Alpha_shape_2<Triangulation_2>  Alpha_shape;

typedef Alpha_shape::Face  Face;
typedef Alpha_shape::Vertex Vertex;
typedef Alpha_shape::Edge Edge;
typedef Alpha_shape::Face_handle  Face_handle;
typedef Alpha_shape::Vertex_handle Vertex_handle;

typedef Alpha_shape::Face_circulator  Face_circulator;
typedef Alpha_shape::Vertex_circulator  Vertex_circulator;

typedef Alpha_shape::Locate_type Locate_type;

typedef Alpha_shape::Face_iterator  Face_iterator;
typedef Alpha_shape::Vertex_iterator  Vertex_iterator;
typedef Alpha_shape::Edge_iterator  Edge_iterator;
typedef Alpha_shape::Edge_circulator  Edge_circulator;
typedef Alpha_shape::Coord_type Coord_type;
typedef Alpha_shape::Alpha_iterator Alpha_iterator;

#include <CGAL/geowin_support.h>


class geo_delaunay_triang : public geowin_update<std::list<CGALPoint>, std::list<CGALSegment> >,
                            public geowin_redraw
{
 Triangulation_2 dt;  

public:
 virtual ~geo_delaunay_triang() {}


 bool insert(const CGALPoint& p)
 {
  std::cout << "insert:" << p << "\n";
  dt.insert(p);
  return true; 
 }
 
 bool del(const CGALPoint& p)
 {
  std::cout << "del:" << p << "\n";
  Triangulation_2::Vertex_handle vh = dt.nearest_vertex(p);
  dt.remove(vh);
  return true;  
 } 

 void draw(leda_window& W,leda_color c1,leda_color c2,double x1,double y1,double x2,double y2)
 {
  Edge_iterator eit = dt.edges_begin();
  Edge_iterator beyond = dt.edges_end();   
  Edge eact;

  while (eit != beyond) {
       eact = *eit;   
       W.draw_segment(convert_to_leda(dt.segment(eact)),c1);                    
       ++eit;  
  }    
 }

 void update(const CGALPointlist& L, CGALSegmentlist&)
 {
  dt.clear();       
  dt.insert(L.begin(),L.end());
 }
};


double alpha_index;
bool reg;

class geo_alpha : public geowin_update<std::list<CGALPoint>, std::list<CGALSegment> >,
                  public geowin_redraw
{
public:
 Alpha_shape A;
 
 virtual ~geo_alpha() { }

 void draw(leda_window& W, leda_color c1, leda_color c2,double x1,double y1,double x2,double y2)
 {  W.set_color(c1); 
    W << A; 
 }  
 
 void update(const CGALPointlist& L, CGALSegmentlist& Sl)
 {
  A.make_alpha_shape(L.begin(), L.end());
  A.set_alpha(alpha_index);
  if (reg) A.set_mode(Alpha_shape::REGULARIZED);
  else A.set_mode(Alpha_shape::GENERAL);
 }
};

int main()
{
  geowin_init_default_type((CGALPointlist*)0, leda_string("CGALPointList"));

  CGALPointlist L;

  GeoWin GW("Alpha shapes");
  GW.add_help_text(leda_string("CGAL_alpha_shape_2"));

  geo_scene my_scene= GW.new_scene(L);  
  GW.set_point_style(my_scene, leda_disc_point);
 
  geo_delaunay_triang deltria;  
  geo_scene res2 = GW.new_scene(deltria, deltria, my_scene, leda_string("Delaunay Triangulation"));
  GW.set_color(res2, leda_yellow);

  geo_alpha alpha;  
  geo_scene res3 = GW.new_scene(alpha, alpha, my_scene, leda_string("Alpha shape"));
  GW.set_line_width(res3, 3);
  GW.set_color(res3, leda_green2);

  GW.set_all_visible(true);
  
  GW.init_menu();
  // add a double item ...
  alpha_index = 1000;
  GW.get_window().double_item("Alpha:",alpha_index);
  reg = false;
  GW.get_window().bool_item("Regularized:",reg);
  GW.edit(my_scene);
  
  return 0;  
}

#endif
