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

/*

Ideen:
- point location einstellbar machen !
- Nutzerinteraktion einstellbar
- counter einblenden 
- Farben

*/

//#define DEBUG
#define CGAL_PROVIDE_LEDA_RAT_KERNEL_TRAITS_3
#define CGAL_NO_DEPRECATED_CODE
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
#include <CGAL/Kernel_special.h>
#include <CGAL/kernel_event_support.h>
#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>
#include <CEP/Leda_rat_kernel/geowin_leda_rat_kernel.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/geowin_support.h>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif


typedef CGAL::leda_rat_kernel_traits                   LEDA_KERNEL;
typedef CGAL::kernel_event<LEDA_KERNEL>                KEV;
typedef CGAL::kernel_event<int>                        KRES;
typedef CGAL::Kernel_special<LEDA_KERNEL, KEV, KRES>   K;

typedef CGAL::Arr_segment_traits_2<K>                  Traits;

typedef K::Point_2                                     IPoint;
typedef K::Point_2                                     Point;
typedef K::Segment_2                                   Segment;
typedef Traits::X_curve                                X_curve;

typedef CGAL::Arr_base_node<X_curve>                   Base_node;
typedef CGAL::Arr_2_default_dcel<Traits>               Dcel;
typedef CGAL::Arrangement_2<Dcel,Traits,Base_node >    Arr_2;

typedef CGAL::Arrangement_2<Dcel,Traits,Base_node >::Curve_iterator  Curve_iterator;
typedef CGAL::Arrangement_2<Dcel,Traits,Base_node >::Vertex_iterator Vertex_iterator;
typedef CGAL::Arrangement_2<Dcel,Traits,Base_node >::Vertex          Vertex;

// functors of the adapted kernel ...
typedef LEDA_KERNEL::Is_vertical_2                      Is_vertical_2;
typedef LEDA_KERNEL::Construct_vertex_2                 Construct_vertex_2;
typedef LEDA_KERNEL::Less_x_2                           Less_x_2;
typedef LEDA_KERNEL::Construct_segment_2                Construct_segment_2;
typedef LEDA_KERNEL::Compare_xy_2                       Compare_xy_2;
		      
typedef LEDA_KERNEL::Compare_x_2                        Compare_x_2;
typedef LEDA_KERNEL::Compare_y_2                        Compare_y_2;
typedef LEDA_KERNEL::Less_y_2                           Less_y_2;
typedef LEDA_KERNEL::Equal_2                            Equal_2;
typedef LEDA_KERNEL::Compare_y_at_x_2                   Compare_y_at_x_2;
typedef LEDA_KERNEL::Compare_slope_2                    Compare_slope_2;
typedef LEDA_KERNEL::Counterclockwise_in_between_2      Counterclockwise_in_between_2; 
typedef LEDA_KERNEL::Construct_direction_2              Construct_direction_2; 
typedef LEDA_KERNEL::Construct_opposite_direction_2     Construct_opposite_direction_2;       



Arr_2 arr;
bool  graphical_output = false;

class geo_locate : public geowin_update<std::list<Point>,std::list<Segment> >
{
public:

 int is_vertical_counter;
 int construct_vertex_counter;
 int less_x_counter;
 int construct_segment_counter;
 int compare_xy_counter;

 CGAL::event_item Is_vertical_2_it; 
 CGAL::event_item Construct_vertex_2_it;
 CGAL::event_item Less_x_2_it;
 CGAL::event_item Construct_segment_2_it;
 CGAL::event_item Compare_xy_2_it;   
 
 // check some others ...
 CGAL::event_item Compare_x_2_it;
 CGAL::event_item Compare_y_2_it;
 CGAL::event_item Less_y_2_it;
 CGAL::event_item Equal_2_point_it;
 CGAL::event_item Equal_2_segment_it;
 CGAL::event_item Compare_y_at_x_2_point_segment_segment_it;
 CGAL::event_item Compare_y_at_x_2_point_segment_it;
 CGAL::event_item Compare_slope_2_it;
 CGAL::event_item Counterclockwise_in_between_2_it;
 CGAL::event_item Construct_direction_2_it;
 CGAL::event_item Construct_opposite_direction_2_it; 
 
 int equal_2_op_point;
 int equal_2_op_segment;
 int compare_x_2_op;
 int compare_y_2_op;
 int less_y_2_op;
 int compare_y_at_x_2_op_point_segment_segment;
 int compare_y_at_x_2_op_point_segment;
 int compare_slope_2_op;
 int counterclockwise_in_between_2_op;
 int construct_direction_2_op;
 int construct_opposite_direction_2_op;
 
 Point   loc_point;
 GeoWin& gw; 
 leda_window& win; 
 
 geo_locate(GeoWin& g) : gw(g), win(g.get_window()) 
 { } 
 
 // ---------------------------------------------------------------------------------------------------------
 // other stuff 
 // --------------------------------------------------------------------------------------------------------- 
 
 void equal_2_op_point_counter(const LEDA_KERNEL::Equal_2&, const Point& p1, const Point& p2)
 { equal_2_op_point++; }
 void equal_2_op_segment_counter(const LEDA_KERNEL::Equal_2&, const Segment& s1, const Segment& s2)
 { equal_2_op_segment++; }
 void compare_x_2_op_counter(const LEDA_KERNEL::Compare_x_2&, const Point& p1, const Point& p2)
 { compare_x_2_op++; }
 void compare_y_2_op_counter(const LEDA_KERNEL::Compare_y_2&, const Point& p1, const Point& p2)
 { compare_y_2_op++; }
 void less_y_2_op_counter(const LEDA_KERNEL::Less_y_2&, const Point& p1, const Point& p2)
 { less_y_2_op++; }
 void compare_y_at_x_2_op_point_segment_segment_counter(const LEDA_KERNEL::Compare_y_at_x_2&, const Point& p, const Point& s1, 
                                                        const Segment& s2)
 { compare_y_at_x_2_op_point_segment_segment++; }
 void compare_y_at_x_2_op_point_segment_counter(const LEDA_KERNEL::Compare_y_at_x_2&, const Point& p, const Segment& s1)
 { compare_y_at_x_2_op_point_segment++; }

 void compare_slope_2_op_counter(const LEDA_KERNEL::Compare_slope_2&,const Segment& s1, const Segment& s2)
 { compare_slope_2_op++; }

 void counterclockwise_in_between_2_op_counter(const LEDA_KERNEL::Counterclockwise_in_between_2&,
                                               const leda::rat_direction & d,
                                               const leda::rat_direction & d1, 
					       const leda::rat_direction & d2)
 { counterclockwise_in_between_2_op++; }
 void construct_direction_2_op_counter(const LEDA_KERNEL::Construct_direction_2&, const Segment& s)
 { construct_direction_2_op++; }
 void construct_opposite_direction_2_op_counter(const LEDA_KERNEL::Construct_opposite_direction_2&, const leda::rat_direction & d)
 { construct_opposite_direction_2_op++; }

 // ---------------------------------------------------------------------------------------------------------
 
 void is_vertical_occurence(const LEDA_KERNEL::Is_vertical_2&, const Segment& obj)
 { std::cout << "is_vertical:" << obj <<  "\n"; 
   is_vertical_counter++;
 }
 
 void construct_vertex_occurence(const LEDA_KERNEL::Construct_vertex_2&, const Segment& obj, const int& i)
 { std::cout << "construct_vertex:" << obj << " " << i << "\n";
   construct_vertex_counter++; 
 }
 
 void less_x_occurence(const LEDA_KERNEL::Less_x_2&, const Point& p1, const Point& p2, const bool&)
 { std::cout << "less_x_2:" << p1 << " " << p2 << "\n"; 
   less_x_counter++;
   
   if (graphical_output){
     Point other;
     other = leda::identical(p1,loc_point) ? p2 : p1;
     win.draw_arrow(loc_point.to_float(), other.to_float(), blue);
     leda_wait(0.5);
   }   
 }
 
 void construct_segment_occurence(const LEDA_KERNEL::Construct_segment_2&, const Point& p1, const Point& p2)
 { std::cout << "construct_segment:" << p1 << " " << p2 << "\n";
   construct_segment_counter++; 
 }  
 
 void compare_xy_occurence(const LEDA_KERNEL::Compare_xy_2&, const Point& p1, const Point& p2)
 { std::cout << "compare_xy_2:" << p1 << " " << p2 << "\n"; 
   compare_xy_counter++;
 } 
 
/*
Problem: nur das eine Resultevent geht - warum ???
*/ 
 
 void init_visualization(const Point& locate)
 {
   loc_point = locate;
 
   // set counters ...
   std::cout << "\n\n";
   is_vertical_counter = 0;
   construct_vertex_counter = 0;
   less_x_counter = 0;
   construct_segment_counter = 0;
   compare_xy_counter = 0;
 
   // others ...
   equal_2_op_point = 0;
   equal_2_op_segment = 0;
   compare_x_2_op = 0;
   compare_y_2_op = 0;
   less_y_2_op = 0;
   compare_y_at_x_2_op_point_segment_segment = 0;
   compare_y_at_x_2_op_point_segment = 0;
   compare_slope_2_op = 0;
   counterclockwise_in_between_2_op = 0;
   construct_direction_2_op = 0;
   construct_opposite_direction_2_op = 0;
       
   Is_vertical_2_it = CGAL::attach(KEV::EVENT, *this, &geo_locate::is_vertical_occurence); 
   Construct_vertex_2_it = CGAL::attach(KEV::EVENT, *this, &geo_locate::construct_vertex_occurence);
   Less_x_2_it = CGAL::attach(KRES::EVENT, *this, &geo_locate::less_x_occurence);
   Construct_segment_2_it = CGAL::attach(KEV::EVENT, *this, &geo_locate::construct_segment_occurence);
   Compare_xy_2_it = CGAL::attach(KEV::EVENT, *this, &geo_locate::compare_xy_occurence);
		      
   Compare_x_2_it = CGAL::attach(KEV::EVENT, *this, &geo_locate::compare_x_2_op_counter);   
   Compare_y_2_it = CGAL::attach(KEV::EVENT, *this, &geo_locate::compare_y_2_op_counter);   
   Less_y_2_it    = CGAL::attach(KEV::EVENT, *this, &geo_locate::less_y_2_op_counter);   
   Equal_2_point_it  = CGAL::attach(KEV::EVENT, *this, &geo_locate::equal_2_op_point_counter);   
   Equal_2_segment_it  = CGAL::attach(KEV::EVENT, *this, &geo_locate::equal_2_op_segment_counter);   
   Compare_y_at_x_2_point_segment_segment_it  = CGAL::attach(KEV::EVENT, *this, &geo_locate::compare_y_at_x_2_op_point_segment_segment_counter);  
   Compare_y_at_x_2_point_segment_it  = CGAL::attach(KEV::EVENT, *this, &geo_locate::compare_y_at_x_2_op_point_segment_counter);   
   Compare_slope_2_it  = CGAL::attach(KEV::EVENT, *this, &geo_locate::compare_slope_2_op_counter);   
   Counterclockwise_in_between_2_it  = CGAL::attach(KEV::EVENT, *this, &geo_locate::counterclockwise_in_between_2_op_counter);   
   Construct_direction_2_it  = CGAL::attach(KEV::EVENT, *this, &geo_locate::construct_direction_2_op_counter);  
   Construct_opposite_direction_2_it  = CGAL::attach(KEV::EVENT, *this, &geo_locate::construct_opposite_direction_2_op_counter);          
 }
 
 void reset_visualization() {
   // output counters ...
   std::cout << "\ncounters:\n";
   std::cout << "---------------\n";     
   std::cout << "is_vertical_counter:" << is_vertical_counter  << "\n";
   std::cout << "construct_vertex_counter:" << construct_vertex_counter << "\n";
   std::cout << "less_x_counter:" << less_x_counter << "\n";
   std::cout << "construct_segment_counter:" << construct_segment_counter << "\n";
   std::cout << "compare_xy_counter:" << compare_xy_counter << "\n";
   std::cout << "\nother counters:\n";
   std::cout << "---------------\n";   
   std::cout << "point equality:" << equal_2_op_point << "\n";
   std::cout << "segment equality:" << equal_2_op_segment << "\n";
   std::cout << "point compare x:" << compare_x_2_op << "\n";
   std::cout << "point compare y:" << compare_y_2_op << "\n";
   std::cout << "point less y:" << less_y_2_op << "\n";
   std::cout << "compare_y_at_x_2 (point/segment/segment):" << compare_y_at_x_2_op_point_segment_segment << "\n";
   std::cout << "compare_y_at_x_2 (point/segment):" << compare_y_at_x_2_op_point_segment << "\n";
   std::cout << "compare_slope:" << compare_slope_2_op << "\n";
   std::cout << "counterclockwise_in_between:" << counterclockwise_in_between_2_op << "\n";
   std::cout << "construct_direction:" << construct_direction_2_op << "\n";
   std::cout << "construct_opposite_direction:" << construct_opposite_direction_2_op << "\n";   
   std::cout << "\n\n";
   
   gw.msg_open(leda::string("construct_vertex_counter:%d",construct_vertex_counter));
   gw.msg_open(leda::string("less_x_counter:%d",less_x_counter));   
    
   CGAL::detach(Is_vertical_2_it);
   CGAL::detach(Construct_vertex_2_it);
   CGAL::detach(Less_x_2_it);
   CGAL::detach(Construct_segment_2_it);
   CGAL::detach(Compare_xy_2_it);  
   
   // others ...
   CGAL::detach(Compare_x_2_it);
   CGAL::detach(Compare_y_2_it);
   CGAL::detach(Less_y_2_it);
   CGAL::detach(Equal_2_point_it);
   CGAL::detach(Equal_2_segment_it);
   CGAL::detach(Compare_y_at_x_2_point_segment_segment_it);
   CGAL::detach(Compare_y_at_x_2_point_segment_it);
   CGAL::detach(Compare_slope_2_it);
   CGAL::detach(Counterclockwise_in_between_2_it);
   CGAL::detach(Construct_direction_2_it);
   CGAL::detach(Construct_opposite_direction_2_it);    
 } 

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
    
    // switch events on ...
    init_visualization(pt);
    e = arr.locate(pt,lt);
    reset_visualization(); 
    // ... and off   

    //color the face on the screen
    Arr_2::Face_handle f=e->face();
    if (f->does_outer_ccb_exist()) {
      Arr_2::Ccb_halfedge_circulator cc=f->outer_ccb();
      do {
        Segment seg = cc->curve();
        Sl.push_back(seg);
      } while (++cc != f->outer_ccb());
      
    }
    
    Arr_2::Holes_iterator hit=f->holes_begin(),eit=f->holes_end();
    for (;hit!=eit; ++hit) {
      Arr_2::Ccb_halfedge_circulator cc=*hit; 
      do {
        Segment seg = cc->curve();
        Sl.push_back(seg);
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
  arr.insert(X_curve(seg));
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
    arr.insert(X_curve(s1));
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
  
  geo_locate locate(GW);
  geo_scene loc  = GW.new_scene(locate,loc_input,leda_string("Arrangement locate operations"));
  GW.set_color( loc, leda_red );
  GW.set_line_width( loc, 3 );
  GW.set_fill_color( loc, leda_red);  
   
  GW.add_dependence(seg_scene,loc);
 
  GW.set_all_visible(true);
  
  GW.message("To perform point locations, activate in the Scenes menu the input scene for point location.");

  GW.init_menu();
  GW.get_window().bool_item(" Output in algorithm:",graphical_output);
  GW.edit(seg_scene);
  
  return 0;  
}

#endif
