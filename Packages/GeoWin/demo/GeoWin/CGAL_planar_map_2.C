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
// file          : demo/GeoWin/CGAL_planar_map_2.C
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
#include <CGAL/Pm_segment_epsilon_traits.h>
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include<list>

#include <CGAL/geowin_support.h>
#include <CGAL/distance_predicates_2.h>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

typedef double                                            number_type;
typedef CGAL::Cartesian<number_type>                      K;
typedef CGAL::Pm_segment_epsilon_traits<K>                pmtraits;
typedef K::Point_2                                        Point;
typedef K::Segment_2                                      Segment;
typedef CGAL::Polygon_traits_2<K>                         PTraits;
typedef CGAL::Polygon_2<PTraits,std::list<Point> >        Polygon;
typedef pmtraits::X_curve                                 Curve;
typedef CGAL::Pm_default_dcel<pmtraits>                   pmdcel;

typedef CGAL::Planar_map_2<pmdcel, pmtraits>	          planar_map;
typedef planar_map::Locate_type                           loc_type;
typedef planar_map::Halfedge_handle                       halfedge_handle;

std::list<Point>   L;
std::list<Point>   shoot;
std::list<Segment> Lseg;
std::list<Point>   Loc;


bool snap_to_points(Segment& seg)
{
 if (L.size() < 2) return false;
 Point ps = seg.source();
 Point pt = seg.target();

 std::list<Point>::const_iterator it = L.begin();
 Point psnew=*it, ptnew=*it;

 for(; it != L.end(); ++it) {
  Point pakt= *it;
  if (CGAL::has_smaller_dist_to_point(ps,pakt,psnew)) psnew=pakt;
 }
 it=L.begin();
 for(; it != L.end(); ++it) {
  Point pakt= *it;
  if (CGAL::has_smaller_dist_to_point(pt,pakt,ptnew)) ptnew=pakt;
 }
 
 seg=Segment(psnew,ptnew);
 if (seg.is_degenerate()) return false; 
 else return true;
}


bool segment_add_changer(GeoWin& gw, const Segment& segold, Segment& segnew)
{
 segnew=segold;
 bool b = snap_to_points(segnew);
 return b;
}

bool segment_start_change(GeoWin& gw, const Segment& segold)
{ return false; }
bool point_start_change(GeoWin& gw, const Point& p)
{ return false; }

CGAL::Planar_map_2<pmdcel, pmtraits>* pm;

class geo_shoot : public geowin_redraw, public geowin_update<std::list<Point>, std::list<Point> >
{
public:
  virtual ~geo_shoot() {}
  
  std::list<Point>     PT;
  std::list<Segment>   PSEG;

  virtual void draw(leda_window& W, leda_color c1, leda_color c2,double x1,double y1,double x2,double y2)
  {
    leda_color cold = W.set_color(leda_blue);
    int lwold = W.set_line_width(3);
    
    std::list<Segment>::const_iterator it = PSEG.begin();
    std::list<Point>::const_iterator pit  = PT.begin();
    
    for(; it != PSEG.end(); it++, pit++){
      W << *it;
      //draw shooting arrow ...
      Point src = *pit;
      double x1,y1,x2,y2,m,n;
      x1 = (*it).source().x(); x2 = (*it).target().x();
      y1 = (*it).source().y(); y2 = (*it).target().y();
      m = (y2-y1)/(x2-x1);
      n = y2-m*x2;
      
      W.set_color(leda_red);
      W.draw_arrow((*pit).x(),(*pit).y(),(*pit).x(),m*((*pit).x())+n);
      W.set_color(leda_blue);
    }
    
    W.set_color(cold);
    W.set_line_width(lwold);    
  }

  virtual void update(const std::list<Point>& L, std::list<Point>&)
  {
    PT.clear(); PSEG.clear(); 
    std::list<Point>::const_iterator it = L.begin();
    loc_type lt;
    
    for(;it != L.end(); it++){
      halfedge_handle he = pm->vertical_ray_shoot(*it,lt,true);
      switch(lt){
        case 1: {   //Vertex
	  break;
	}
	case 2: {   //Edge
          Point pt1 = (*he).source()->point();
	  Point pt2 = (*he).target()->point();
	  PSEG.push_back(Segment(pt1,pt2));
	  PT.push_back(*it);
	  break;
	}
	default: { break; }
      }
    }
  }
};

class geo_plmap : public geowin_redraw, public geowin_update<std::list<Segment>, std::list<Point> >
{
public:

  std::list<Polygon> CPL;
  std::list<Point>  red;    //for location...
  std::list<Point>  green;
  std::list<Point>  white;

  virtual ~geo_plmap() {}

  virtual void draw(leda_window& W, leda_color c1, leda_color c2,double x1,double y1,double x2,double y2)
  {     
    std::list<Polygon>::const_iterator it = CPL.begin();
    int cnt = 0;
    for(; it != CPL.end(); ++it, cnt ++) {    
       W.set_fill_color(cnt % 7 + 6);
       Polygon pol= *it;
       if (pol.is_simple()) W << convert_to_leda(pol); 
    }
    // output the lists for location ...
    std::list<Point>::const_iterator rt;

    W.set_color(leda_red);
    for (rt=red.begin(); rt != red.end(); rt++) W.draw_disc((*rt).x(),(*rt).y(),5);
    W.set_color(leda_green);
    for (rt=green.begin(); rt != green.end(); rt++) W.draw_disc((*rt).x(),(*rt).y(),5);
    W.set_color(leda_white);
    for (rt=white.begin(); rt != white.end(); rt++) W.draw_disc((*rt).x(),(*rt).y(),5);
  }
  
  virtual bool write_postscript(ps_file& PS, leda_color c1, leda_color c2)
  {     
    std::list<Polygon>::const_iterator it = CPL.begin();
    int cnt = 0;
    for(; it != CPL.end(); ++it, cnt ++) {    
       PS.set_fill_color(cnt % 7 + 6);
       Polygon pol= *it;
       if (pol.is_simple()) PS << convert_to_leda(pol); 
    }
    // output the lists for location ...
    std::list<CGALPoint>::const_iterator rt;

    PS.set_color(leda_red);
    for (rt=red.begin(); rt != red.end(); rt++) PS.draw_disc((*rt).x(),(*rt).y(),5);
    PS.set_color(leda_green);
    for (rt=green.begin(); rt != green.end(); rt++) PS.draw_disc((*rt).x(),(*rt).y(),5);
    PS.set_color(leda_white);
    for (rt=white.begin(); rt != white.end(); rt++) PS.draw_disc((*rt).x(),(*rt).y(),5);
    return false;
  }


  virtual void update(const std::list<Segment>& L, std::list<Point>&)
  {
     CPL.clear();
     pm->clear();
     
     std::list<Segment>::const_iterator it = L.begin();

     for(; it != L.end(); ++it) {
        Segment sakt= *it;
        pm->insert(Curve(sakt.source(),sakt.target()));
     }

     if (! pm->is_valid()) std::cout << "not a valid map!\n";

     int cnt = 0;
     CGAL::Planar_map_2<pmdcel, pmtraits>::Face_iterator fit= pm->faces_begin();
     for(; fit != pm->faces_end(); ++fit) {  
        CGAL::Planar_map_2<pmdcel, pmtraits>::Face fc = *fit;

        // now lets do something with the face ...
        if (! fc.is_unbounded()) {

           CGAL::Planar_map_2<pmdcel, pmtraits>::Halfedge_handle he, he_next;
           he= fc.halfedge_on_outer_ccb();
	   he_next = he;
	   std::list<Point> pcon;

	   do {
             Point pt = (*he_next).source()->point();
             pcon.push_back(pt);
	     he_next = (*he_next).next_halfedge();
             
	   } while (he!=he_next);
           
           Polygon poly(pcon.begin(), pcon.end());
           CPL.push_back(poly);
        }

        cnt++;
     }      

     green.clear(); red.clear(); white.clear();

     std::list<Point>::const_iterator loc_it = Loc.begin();
     loc_type lt;

     for(; loc_it != Loc.end(); ++loc_it) {
           pm->locate(*loc_it,lt);        
           if (lt==4) green.push_back(*loc_it);
           if (lt==3) red.push_back(*loc_it);
           if (lt==1 || lt==2) white.push_back(*loc_it);
     }  
  }
};


int main()
{
  geowin_init_default_type((std::list<Point>*)0, leda_string("CGALPointList"));
  geowin_init_default_type((std::list<Segment>*)0, leda_string("CGALSegmentList"));
  
  CGAL::Planar_map_2<pmdcel, pmtraits> pmap;
  pm = &pmap;
 
  GeoWin GW("CGAL - Planar map demo");
  GW.add_help_text(leda_string("CGAL_planar_map_2"));

  geo_scene point_scene   = GW.new_scene(L);  

  GeoEditScene<std::list<Point> >* PTR1 =
       (GeoEditScene<std::list<Point> >*)point_scene;
  GW.set_start_change_handler(PTR1, point_start_change);
  
  // vertical ray shooting scene ...
  geo_scene shoot_scene  = GW.new_scene(shoot);
  GW.set_color(shoot_scene, leda_blue);
  GW.set_point_style(shoot_scene, leda_disc_point);
 
  geo_scene segment_scene = GW.new_scene(Lseg); 

  GeoEditScene<std::list<Segment> >* PTR2 =
       (GeoEditScene<std::list<Segment> >*)segment_scene;

  GW.set_pre_add_change_handler(PTR2,segment_add_changer);
  GW.set_start_change_handler(PTR2,segment_start_change);

  // second point scene for locating ...
  
  geo_scene ploc_scene   = GW.new_scene(Loc);  
  GW.set_color(ploc_scene,leda_black);

  geo_plmap PM;
 
  // Result Scenes ...
  geo_scene sc1 = GW.new_scene( PM, PM, segment_scene, "Planar Map"); 
  GW.set_color(sc1, leda_blue);
  
  geo_shoot SHT;
  geo_scene sc2 = GW.new_scene( SHT, SHT, shoot_scene, "Vertical Ray shooting");   
  
  GW.set_all_visible(true);

  GW.add_dependence(segment_scene,sc2);
  GW.add_dependence(ploc_scene,sc1);
  GW.edit(point_scene);
  
  return 0;  
}

#endif
