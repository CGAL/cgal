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
// file          : demo/GeoWin/CGAL_geometry_2.C
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

#include <CGAL/Cartesian.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Min_circle_2.h>
#include <CGAL/Min_circle_2_traits_2.h>
#include <CGAL/Min_ellipse_2.h>
#include <CGAL/Min_ellipse_2_traits_2.h>
#include <CGAL/IO/Window_stream.h>
#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/geowin_support.h>
#include <set>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

typedef CGAL::Cartesian<double>                                   K;
typedef K::Point_2                                                Point;
typedef K::Segment_2                                              Segment;
typedef K::Circle_2                                               Circle;
typedef K::Ray_2                                                  Ray;

typedef  CGAL::Min_circle_2_traits_2<K>                           Circ_Traits;
typedef  CGAL::Min_circle_2<Circ_Traits>                          Min_circle;
typedef  Min_circle::Circle                                       OptCircle;

typedef  CGAL::Min_ellipse_2_traits_2<K>                          Ell_Traits;
typedef  CGAL::Min_ellipse_2<Ell_Traits>                          Min_ellipse;

typedef CGAL::Triangulation_euclidean_traits_2<K>                 Gt;
typedef CGAL::Triangulation_vertex_base_2<Gt>                     Vb;
typedef CGAL::Triangulation_face_base_2<Gt>                       Fb;
typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb>    Tds;
typedef CGAL::Triangulation_2<Gt,Tds>                             Triangulation_2;
typedef CGAL::Delaunay_triangulation_2<Gt,Tds>                    Delaunay_triangulation_2;

typedef Triangulation_2::Edge                                     Edge;
typedef Triangulation_2::Locate_type                              Locate_type;
typedef Triangulation_2::Face_handle                              Face_handle;
typedef Triangulation_2::Edge_iterator                            Edge_iterator;
typedef Triangulation_2::Face_iterator                            Face_iterator;

typedef K::Less_xy_2                                              Point_compare;


class geo_hull : public geowin_update<std::list<Point>, std::list<Segment> >
{
public:
 void update(const std::list<Point>& L, std::list<Segment>& Slst)
 {
  Slst.clear();
  std::list<Point> out;
  CGAL::convex_hull_points_2(L.begin(),L.end(), std::back_inserter(out));   

  // build the segment list ...
  if( out.size() > 1 ) {
    Point pakt,prev,pstart;

    std::list<Point>::const_iterator it;
    it=out.begin();
    prev= *it; pstart=prev;
    it++;

    for(; it != out.end(); ++it) {
       pakt= *it;
       Slst.push_back(Segment(prev,pakt));
       prev=pakt;
    }
    Slst.push_back(Segment(pakt,pstart));
  }
 }
 
};

class geo_triang : public geowin_update<std::list<Point>, std::list<Segment> >
{
public:
 void update(const std::list<Point>& L, std::list<Segment>& Sl)
 {
  Triangulation_2 tr;    
  Sl.clear();      
                   
  std::list<Point>::const_iterator it;
  it= L.begin();
  Point pakt;
 
  for (; it != L.end() ; ++it) {
        pakt= *it;
        tr.push_back(pakt);                            
  }

  Edge_iterator eit = tr.edges_begin();
  Edge_iterator beyond = tr.edges_end();   
  Edge eact;

  while (eit != beyond) {
       eact = *eit;         
       Sl.push_back(tr.segment(eact));               
       ++eit;  
  }         
 }
};

class geo_crust : public geowin_update<std::list<Point>, std::list<Segment> >
{
public:
 void update(const std::list<Point>& L, std::list<Segment>& Sl)
 {
  Delaunay_triangulation_2 dt;    
  Sl.clear();      
  dt.insert(L.begin(),L.end());

  Face_iterator f;

  typedef std::set<Point, Point_compare> point_set;
  point_set s;
  for (f=dt.faces_begin(); f!=dt.faces_end(); ++f)
    s.insert( dt.dual(f) );

  point_set::iterator p;
  for (p=s.begin(); p!=s.end(); ++p) dt.insert(*p);  
  
  Edge_iterator e;
  for (e=dt.edges_begin(); e!=dt.edges_end(); ++e) {
    Face_handle f=(*e).first; int i=(*e).second;
    bool s1=s.find(f->vertex(f->ccw(i))->point())==s.end();
    bool s2=s.find(f->vertex(f->cw(i))->point())==s.end();
    if (s1&&s2) {
      Sl.push_back(dt.segment(e));
    }
  }
  
 } 
};

#if !defined(_MSC_VER)
class geo_ellipse : public geowin_redraw, public geowin_update<std::list<Point>, std::list<Point> >
{
public:
  virtual ~geo_ellipse() {}

  Min_ellipse min_ell;

  void draw(leda_window& W, leda_color c1, leda_color c2,double x1,double y1,double x2,double y2)
  {  W.set_color(c1); W << min_ell.ellipse(); }  

  void update(const std::list<Point>& L, std::list<Point>& LP)
  {  min_ell.clear(); min_ell.insert(L.begin(),L.end()); }
};
#endif

class geo_circ : public geowin_update<std::list<Point>,std::list<Circle> >
{
public:
 void update(const std::list<Point>& L,std::list<Circle>& Cl)
 {
   Cl.clear();
   if (L.size() < 2) return;

   Min_circle  mc1( L.begin(), L.end(), true);
   OptCircle ci= mc1.circle();

   Point ctp=ci.center();
   Circle conv(ctp,ci.squared_radius());
   Cl.push_back(conv); 
 }
};

Delaunay_triangulation_2 mydt;

class geo_delaunay_triang : public geowin_update<std::list<Point>, std::list<Segment> >,
                            public geowin_redraw
{
public:

 bool insert(const Point& p)
 {
  //std::cout << "insert:" << p << "\n";
  mydt.insert(p);
  return true; 
 }
 
 bool del(const Point& p)
 {
  //std::cout << "del:" << p << "\n";
  Delaunay_triangulation_2::Vertex_handle vh = mydt.nearest_vertex(p);
  mydt.remove(vh);
  return true;  
 }

 void draw(leda_window& W,leda_color c1,leda_color c2,double x1,double y1,double x2,double y2)
 {
  Edge_iterator eit = mydt.edges_begin();
  Edge_iterator beyond = mydt.edges_end();   
  Edge eact;

  while (eit != beyond) {
       eact = *eit;   
       W.draw_segment(convert_to_leda(mydt.segment(eact)),c1);                    
       ++eit;  
  }    
 }

 void update(const std::list<Point>& L, std::list<Segment>&)
 {
  //std::cout << "recompute!\n";
  mydt.clear();
  mydt.insert(L.begin(),L.end());  
 }
};

class geo_locate : public geowin_update<std::list<Point>, std::list<Circle> >
{
public:
 void update(const std::list<Point>& L, std::list<Circle>& Sl)
 {
    Sl.clear();
    Locate_type lt;
    int li;

    typedef Delaunay_triangulation_2::Face_handle   DT_Face_handle;
 
    std::list<Point>::const_iterator it;
    it= L.begin();
    Point pakt;
 
    for (; it != L.end() ; ++it) {
      pakt = *it;
      DT_Face_handle f = mydt.locate(pakt,lt, li);
      Circle c;
      if (lt==Delaunay_triangulation_2::FACE) {
        Point p0=f->vertex(0)->point();
        Point p1=f->vertex(1)->point();
        Point p2=f->vertex(2)->point();
        c=Circle(p0,p1,p2);
	Sl.push_back(c);
      }
    }
 }
};

geo_scene res2;

class geo_voro1 : public geowin_update<std::list<Point>, std::list<Segment> >
{
public:
 void update(const std::list<Point>& L, std::list<Segment>& Sl)
 {
  GeoWin* gw = GeoWin::get_call_geowin();
 
  bool visi = gw->get_visible(res2);
 
  Delaunay_triangulation_2 dt;
  Delaunay_triangulation_2* delptr;
  
  if (! visi) { // use own DT         
    dt.insert(L.begin(),L.end()); 
    delptr = &dt;
  } 
  else delptr = &mydt; // use other DT
  
  Sl.clear();  

  Edge_iterator eit, est=delptr->edges_begin(), eend=delptr->edges_end();
  for (eit=est; eit != eend; ++eit){
    Segment s;
    CGAL::Object o = delptr->dual(eit);
    if (CGAL::assign(s,o)) Sl.push_back(s);
  }
 }
};

class geo_voro2 : public geowin_update<std::list<Point>, std::list<Ray> >
{
public:
 void update(const std::list<Point>& L, std::list<Ray>& Sl)
 {
  GeoWin* gw = GeoWin::get_call_geowin();
 
  bool visi = gw->get_visible(res2);
 
  Delaunay_triangulation_2 dt;
  Delaunay_triangulation_2* delptr;
  
  if (! visi) { // use own DT         
    dt.insert(L.begin(),L.end()); 
    delptr = &dt;
  } 
  else delptr = &mydt; // use other DT
  
  Sl.clear();  

  Edge_iterator eit, est=delptr->edges_begin(), eend=delptr->edges_end();
  for (eit=est; eit != eend; ++eit){
    Ray s;
    CGAL::Object o = delptr->dual(eit);
    if (CGAL::assign(s,o)) Sl.push_back(s);
  }
 }
};



int main()
{
  geowin_init_default_type((std::list<Point>*)0, leda_string("CGALPointList"));

  std::list<Point> L, L2;

  GeoWin GW("CGAL - 2d geometry demo");
  // help for the user ...
  GW.add_help_text(leda_string("CGAL_geometry_2"));
  
  // black background
  GW.get_window().set_bg_color(leda_black);

  // create input scene for DT, Triangulation, VD, ...
  geo_scene my_scene= GW.new_scene(L);  
  GW.set_point_style(my_scene, leda_box_point);
  GW.set_color(my_scene, leda_yellow); 
  
  // second input scene for location of points
  geo_scene loc_input= GW.new_scene(L2);
  GW.set_point_style(loc_input, leda_disc_point);
  GW.set_color(loc_input, leda_green);   

  // create result scenes ...
  // for traingulation
  geo_triang triangulate;
  geo_scene res1 = GW.new_scene(triangulate ,my_scene , leda_string("Triangulation"));
  GW.set_color(res1, leda_orange);
 
  // for delaunay triangulation
  geo_delaunay_triang deltria;  
  res2 = GW.new_scene(deltria, deltria, my_scene, leda_string("Delaunay Triangulation"));
  GW.set_line_width(res2, 2);
  GW.set_color(res2, leda_red);

  // for voronoi diagram
  geo_voro1 voro_edges;
  geo_scene res3 = GW.new_scene(voro_edges, my_scene, leda_string("Voronoi Edges"));
  GW.set_color(res3, leda_blue2);

  geo_voro2 voro_rays;
  geo_scene res4 = GW.new_scene(voro_rays, my_scene, leda_string("Voronoi Edges"));
  GW.set_color(res4, leda_blue2);

  // for the convex hull
  geo_hull hull_update;  
  geo_scene res5 = GW.new_scene(hull_update,my_scene,leda_string("Convex Hull"));  
  GW.set_line_width(res5, 3);
  GW.set_color(res5, leda_blue2);
  
  // for the minimum enclosing circle
  geo_circ min_circ;
  geo_scene res6  = GW.new_scene(min_circ ,my_scene , leda_string("Minimal enclosing circle"));
  GW.set_color(res6,leda_green2);
  GW.set_fill_color(res6,leda_invisible);
  GW.set_line_width(res6, 3);

#if !defined(_MSC_VER)
  // for the minimum enclosing ellipse
  geo_ellipse EL;
  geo_scene res7  = GW.new_scene(EL, EL , my_scene, leda_string("Minimal enclosing ellipse"));
  GW.set_color(res7, leda_white);
#endif
  
  // for the crust
  geo_crust crust_update;  
  geo_scene res8 = GW.new_scene(crust_update,my_scene,leda_string("Crust"));  
  GW.set_color(res8, leda_cyan);  
  GW.set_line_width(res8,3);
  
  // for the point location
  // (this is the circles scene in "delaunay and circles" in the Open GL demo)
  geo_locate pt_loc;
  geo_scene res9 = GW.new_scene(pt_loc,loc_input,leda_string("Locate on DT"));  
  GW.set_color(res9, leda_blue2); 
  GW.set_fill_color(res9, leda_invisible); 
  GW.set_line_width(res9, 2);    
  
  GW.set_all_visible(true);
  
  GW.edit(my_scene);
  
  return 0;  
}

#endif
