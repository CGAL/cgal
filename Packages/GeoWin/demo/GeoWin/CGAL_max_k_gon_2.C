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
// file          : demo/GeoWin/CGAL_max_k_gon_2.C
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
#include <CGAL/extremal_polygon_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <LEDA/list.h>
#include<list>
#include<vector>
#include <CGAL/geowin_support.h>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

typedef CGAL::Cartesian<double>           K;
typedef CGAL::Polygon_traits_2<K>         Traits;
typedef K::Segment_2                      Segment;
typedef Traits::Point_2                   Point;
typedef std::vector<Point>                Container;
typedef CGAL::Polygon_2<Traits,Container> Polygonvec;


class k_gon_5 : public geowin_update<std::list<Point>, std::list<leda_polygon> >
{
public:

 void update(const std::list<Point>& L, std::list<leda_polygon>& Sl)
 {
  Sl.clear();
  Polygonvec p;

  // building convex polygon...
  CGAL::convex_hull_points_2(L.begin(),L.end(), std::back_inserter(p));   

  Polygonvec kg;
  if (p.size()>4) {
      CGAL::maximum_area_inscribed_k_gon(p.vertices_begin(),p.vertices_end(),5, std::back_inserter(kg));

     leda_list<leda_point> HLP;
     int i;

     std::vector<Point> CT = kg.container();

     for (i=0;i<kg.size();i++){
       leda_point lakt(convert_to_leda(CT[i]));
       HLP.append(lakt);
     }
     leda_polygon back(HLP); 
     Sl.push_back(back);
  }
 }
};

class k_gon_3 : public geowin_update<std::list<Point>, std::list<leda_polygon> >
{
public:

 void update(const std::list<Point>& L, std::list<leda_polygon>& Sl)
 { 
  Sl.clear();
  Polygonvec p;

  // building convex polygon...
  CGAL::convex_hull_points_2(L.begin(),L.end(), std::back_inserter(p));   

  Polygonvec kg;
  if (p.size()>2) {
     CGAL::maximum_area_inscribed_k_gon(p.vertices_begin(),p.vertices_end(),3, std::back_inserter(kg));

     leda_list<leda_point> HLP;
     int i;

     std::vector<Point> CT = kg.container();

     for (i=0;i<kg.size();i++){
       leda_point lakt(convert_to_leda(CT[i]));
       HLP.append(lakt);
     }
     leda_polygon back(HLP); 
     Sl.push_back(back);
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

  GeoWin GW("CGAL - Maximum Area inscribed k-gon of a convex polygon");
 
  // build a new edit scene
  geo_scene my_scene= GW.new_scene(L);
  GW.set_point_style(my_scene, leda_disc_point);  

  conv_hull_seg CHS;
  geo_scene hull_scene  = GW.new_scene(CHS,my_scene,leda_string("Convex Hull"));
  GW.set_visible(hull_scene,true);

  k_gon_5 K5;
  geo_scene result  = GW.new_scene(K5,my_scene,leda_string("maximum inscribed 5-gon"));
  GW.set_visible(result,true);
 
  k_gon_3 K3;
  geo_scene result2  = GW.new_scene(K3,my_scene,leda_string("maximum inscribed triangle"));
  GW.set_fill_color(result2, leda_green);
  GW.set_visible(result2,true);

  GW.message("We compute the maximum inscribed 3/5-gon of the convex hull of the input");
  GW.edit(my_scene);
  return 0;  
}

#endif
