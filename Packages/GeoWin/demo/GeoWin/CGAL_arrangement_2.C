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
// file          : demo/GeoWin/CGAL_arrangement_2.C
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
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arr_leda_segment_exact_traits.h>
#include <CGAL/Pm_walk_along_line_point_location.h>

typedef CGAL::Arr_leda_segment_exact_traits            Traits;

typedef Traits::Point                                  IPoint;
typedef CGAL::Cartesian<double>                        K;
typedef K::Point_2                                     Point;
typedef K::Segment_2                                   Segment;
typedef Traits::X_curve                                X_curve;

typedef CGAL::Arr_base_node<X_curve>                   Base_node;
typedef CGAL::Arr_2_default_dcel<Traits>               Dcel;
typedef CGAL::Arrangement_2<Dcel,Traits,Base_node >    Arr_2;

typedef CGAL::Arrangement_2<Dcel,Traits,Base_node >::Curve_iterator  Curve_iterator;
typedef CGAL::Arrangement_2<Dcel,Traits,Base_node >::Vertex_iterator Vertex_iterator;
typedef CGAL::Arrangement_2<Dcel,Traits,Base_node >::Vertex          Vertex;

#include <CGAL/geowin_support.h>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

Arr_2 arr;

class geo_locate : public geowin_update<std::list<Point>,std::list<Segment> >
{
public:
 void update(const std::list<Point>& LP, std::list<Segment>& Sl)
 {
  Sl.clear();
  
  // empty ?
  Curve_iterator cit = arr.curve_node_begin();
  
  if (cit==arr.curve_node_end()) return; //empty
 
  // perform locate operations (and compute output)...
  Arr_2::Locate_type lt;
  Arr_2::Halfedge_handle e;  
    
  std::list<Point>::const_iterator pit= LP.begin();
  
  for(;pit != LP.end(); ++pit){
    Point pt = *pit;
    
    e = arr.locate(IPoint(pt.x(), pt.y()),lt);

    //color the face on the screen
    Arr_2::Face_handle f=e->face();
    if (f->does_outer_ccb_exist()) {
      Arr_2::Ccb_halfedge_circulator cc=f->outer_ccb();
      do {
        leda_segment seg = (cc->curve()).to_float();
	Segment cseg(Point(seg.xcoord1(),seg.ycoord1()), Point(seg.xcoord2(),seg.ycoord2()) );   
      
        Sl.push_back(cseg);
      } while (++cc != f->outer_ccb());
      
    }
    
    Arr_2::Holes_iterator hit=f->holes_begin(),eit=f->holes_end();
    for (;hit!=eit; ++hit) {
      Arr_2::Ccb_halfedge_circulator cc=*hit; 
      do {
        leda_segment seg = (cc->curve()).to_float();
	Segment cseg(Point(seg.xcoord1(),seg.ycoord1()), Point(seg.xcoord2(),seg.ycoord2()) );
	
        Sl.push_back(cseg);
      } while (++cc != *hit);
    }
   }     
 }
};


class geo_seg_arr : public geowin_update<std::list<Segment>,std::list<Segment> >,
                    public geowin_redraw
{
public:

 bool insert(const Segment& seg)
 {
  leda_rat_segment ls(convert_to_leda(seg));
  arr.insert(X_curve(ls.start(), ls.end() ));
  //std::cout << "insert:" << seg << "\n";
  return true; 
 } 

 void draw(leda_window& W, leda_color c1, leda_color c2,double x1,double y1,double x2,double y2)
 {  leda_color cold = W.set_color(c1); 
    // draw vertices of arrangement
    Vertex_iterator vit = arr.vertices_begin();
    
    for(;vit != arr.vertices_end();vit++){
     Vertex vpm = *vit;
     W << vpm.point();
    }
    
    W.set_color(cold);
 } 

 void update(const std::list<Segment>& L, std::list<Segment>& Sl)
 {  
  // clear Arrangement
  arr.clear();
 
  // insert segments ...
  std::list<Segment>::const_iterator it1= L.begin();
  
  for(;it1 != L.end(); ++it1){
    const Segment& s1= *it1;
    leda_rat_segment ls(convert_to_leda(s1));
    arr.insert(X_curve(ls.start(), ls.end() ));
  }
 }
};

int main()
{
  geowin_init_default_type((std::list<Segment>*)0, leda_string("CGALSegmentList"));
  geowin_init_default_type((std::list<Point>*)0, leda_string("CGALPointList"));
 
  std::list<Segment> L;
  std::list<Point>   LP;

  GeoWin GW("CGAL - Segment Arrangements and locate operations");
 
  geo_scene seg_scene = GW.new_scene(L);  
  GW.set_name(seg_scene,"Input for line segment arrangement");
  GW.set_color( seg_scene, leda_green );
  GW.set_line_width( seg_scene, 3 );  
  GW.set_active_line_width( seg_scene, 3 );
  
  geo_scene loc_input = GW.new_scene(LP);
  GW.set_name(loc_input, "Input for point location");
  GW.set_point_style(loc_input, leda_disc_point);

  geo_seg_arr arr;
  //geo_scene res  = 
  GW.new_scene(arr,arr,seg_scene,leda_string("Segment arrangement"));
  
  geo_locate locate;
  geo_scene loc  = GW.new_scene(locate,loc_input,leda_string("Arrangement locate operations"));
  GW.set_color( loc, leda_red );
  GW.set_line_width( loc, 3 );
  GW.set_fill_color( loc, leda_red);  
   
  GW.add_dependence(seg_scene,loc);
 
  GW.set_all_visible(true);
  
  GW.message("To perform point locations, activate in the Scenes menu the input scene for point location.");
  GW.edit(seg_scene);
  
  return 0;  
}

#endif
