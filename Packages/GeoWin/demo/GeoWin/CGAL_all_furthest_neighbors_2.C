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
// file          : demo/GeoWin/CGAL_all_furthest_neighbors_2.C
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

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/random_convex_set_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/all_furthest_neighbors_2.h>
#include <list>
#include <vector>
#include <CGAL/geowin_support.h>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

typedef CGAL::Cartesian<double>           K;
typedef CGAL::Polygon_traits_2<K>         Traits;
typedef Traits::Point_2                   Point;
typedef Traits::Segment_2                 Segment;
typedef std::vector<Point>                Container;
typedef CGAL::Polygon_2<Traits,Container> Polygon;

class geo_all_furthest_nb : public geowin_update<std::list<Point>,std::list<Segment> >,
                            public geowin_redraw
{
public:
 std::list<Segment> ST;
 
 virtual ~geo_all_furthest_nb() {}

 void draw(leda_window& W,leda_color c1,leda_color c2,double x1,double y1,double x2,double y2)
 {
  std::list<Segment>::const_iterator it;
  
  for (it=ST.begin(); it != ST.end(); it++){
    Segment seg = *it;
    W.draw_arrow(seg.source().x(), seg.source().y(), seg.target().x(), seg.target().y(), c1);
  } 
 }

 void update(const std::list<Point>& L, std::list<Segment>& Sl)
 {
  ST.clear();
  Polygon P;

  // building convex polygon...
  CGAL::convex_hull_points_2(L.begin(),L.end(), std::back_inserter(P));   
  
  if (P.size()<2) return;

  std::list<int> il;

  CGAL::all_furthest_neighbors(P.vertices_begin(), P.vertices_end(), std::back_inserter(il));

  std::list<int>::const_iterator lit= il.begin();
  int z=0;
  
  std::vector<Point> CT = P.container();

  for(; lit != il.end(); ++lit) {
     Point p1= CT[z];
     Point p2= CT[*lit]; 

     ST.push_back(Segment(p1,p2));
     z++;
  }
 }
};

class conv_hull_seg : public geowin_update<std::list<Point>,std::list<Segment> >
{
public:

 void update(const std::list<Point>& L, std::list<Segment>& Sl)
 {
  Sl.clear();
  std::list<Point> out;

  CGAL::convex_hull_points_2(L.begin(),L.end(), std::back_inserter(out));   

  // building the segment list ...
  if( out.size() > 1 ) {
    Point pakt,prev,pstart;

    std::list<Point>::const_iterator it;
    it=out.begin();
    prev= *it; pstart=prev;
    it++;

    for(; it != out.end(); ++it) {
       pakt= *it;
       Sl.push_back(Segment(prev,pakt));
       prev=pakt;
    }

    Sl.push_back(Segment(pakt,pstart));
  }
 }
};

int main()
{
  geowin_init_default_type((std::list<Point>*)0, leda_string("CGALPointList"));
 
  std::list<Point> L;

  GeoWin GW("CGAL - All furthest neighbors");

  geo_scene my_scene= GW.new_scene(L);  
  GW.set_point_style(my_scene, leda_disc_point); 
  
  conv_hull_seg CHS;
  geo_scene hull_scene  = GW.new_scene(CHS,my_scene,leda_string("Convex Hull"));
  GW.set_visible(hull_scene,true);

  geo_all_furthest_nb AFN;
  geo_scene result  = GW.new_scene(AFN, AFN, my_scene,leda_string("All furthest neighbors"));
  GW.set_color(result,leda_red);
  GW.set_visible(result,true);
 
  GW.message("The furthest neighbors of the vertices of the convex hull are displayed");
  GW.edit(my_scene);
  return 0;  
}

#endif
