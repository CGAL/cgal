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
// file          : demo/GeoWin/CGAL_triangulation_2.C
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
#include <CGAL/squared_distance_2.h> 
#include <CGAL/Point_2.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

typedef CGAL::Cartesian<double>  Rep;
typedef CGAL::Point_2<Rep>  Point;
typedef CGAL::Segment_2<Rep>  Segment;
typedef CGAL::Ray_2<Rep> Ray;

typedef CGAL::Triangulation_euclidean_traits_2<Rep> Gt;
typedef CGAL::Triangulation_vertex_base_2<Gt> Vb;
typedef CGAL::Triangulation_face_base_2<Gt>  Fb;
typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;
typedef CGAL::Triangulation_2<Gt,Tds>  Triangulation_2;
typedef CGAL::Delaunay_triangulation_2<Gt,Tds> Delaunay_triangulation_2;

typedef Triangulation_2::Edge Edge;
typedef Triangulation_2::Locate_type Locate_type;
typedef Triangulation_2::Edge_iterator  Edge_iterator;

#include <CGAL/Object.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/geowin_support.h>

class geo_triang : public geowin_update<std::list<CGALPoint>, std::list<CGALSegment> >
{
public:
 void update(const CGALPointlist& L, CGALSegmentlist& Sl)
 {
  Triangulation_2 tr;    
  Sl.clear();      
                   
  std::list<CGALPoint>::const_iterator it;
  it= L.begin();
  CGALPoint pakt;
 
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

Delaunay_triangulation_2 dt;

class geo_delaunay_triang : public geowin_update<std::list<CGALPoint>, std::list<CGALSegment> >
{
public:
 void update(const CGALPointlist& L, CGALSegmentlist& Sl)
 {
  dt.clear();    
  Sl.clear();      

  dt.insert(L.begin(),L.end());

  Edge_iterator eit = dt.edges_begin();
  Edge_iterator beyond = dt.edges_end();   
  Edge eact;

  while (eit != beyond) {
       eact = *eit;         
       Sl.push_back(dt.segment(eact));               
       ++eit;  
  }       
 }
};

geo_scene res2;


class geo_voro1 : public geowin_update<std::list<CGALPoint>, std::list<CGALSegment> >
{
public:
 void update(const CGALPointlist& L, CGALSegmentlist& Sl)
 {
  GeoWin* gw = GeoWin::get_call_geowin();
 
  bool visi = gw->get_visible(res2);
 
  Delaunay_triangulation_2 delt;
  Delaunay_triangulation_2* delptr;
  
  if (! visi) { // use own DT         
    delt.insert(L.begin(),L.end()); 
    delptr = &delt;
  } 
  else delptr = &dt; // use other DT
  
  Sl.clear();  

  Edge_iterator eit, est=delptr->edges_begin(), eend=delptr->edges_end();
  for (eit=est; eit != eend; ++eit){
    Segment s;
    CGAL::Object o = delptr->dual(eit);
    if (CGAL::assign(s,o)) Sl.push_back(s);
  }
 }
};

class geo_voro2 : public geowin_update<std::list<CGALPoint>, std::list<CGALRay> >
{
public:
 void update(const CGALPointlist& L, CGALRaylist& Sl)
 {
  GeoWin* gw = GeoWin::get_call_geowin();
 
  bool visi = gw->get_visible(res2);
 
  Delaunay_triangulation_2 delt;
  Delaunay_triangulation_2* delptr;
  
  if (! visi) { // use own DT         
    delt.insert(L.begin(),L.end()); 
    delptr = &delt;
  } 
  else delptr = &dt; // use other DT
  
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
  geowin_init_default_type((CGALPointlist*)0, leda_string("CGALPointList"));

  CGALPointlist L;

  GeoWin GW("CGAL - Triangulation demo");
  GW.add_help_text(leda_string("CGAL_triangulation_2"));

  geo_scene my_scene= GW.new_scene(L);  
  GW.set_point_style(my_scene, leda_disc_point);

  geo_triang triangulate;
  geo_scene res1 = GW.new_scene(triangulate ,my_scene , leda_string("Triangulation"));
  GW.set_color(res1, leda_blue2);
 
  geo_delaunay_triang deltria;  
  res2 = GW.new_scene(deltria, my_scene, leda_string("Delaunay Triangulation"));
  GW.set_line_width(res2, 2);
  GW.set_color(res2, leda_red);

  geo_voro1 voro_edges;
  geo_scene res3 = GW.new_scene(voro_edges, my_scene, leda_string("Voronoi Edges"));
  GW.set_color(res3, leda_blue);

  geo_voro2 voro_rays;
  geo_scene res4 = GW.new_scene(voro_rays, my_scene, leda_string("Voronoi Edges"));
  GW.set_color(res4, leda_blue);
  GW.set_all_visible(true);
  
  GW.edit(my_scene);
  
  return 0;  
}

#endif
