//Achtung - LEDA kernel macht Probleme ohne audgeschaltete POSTCONDITIONS
//rausfinden, warum ! - Problem war Fehler in Less_rotate_ccw_2, fixed
//Loesung: Kernel_traits fuer LEDA spezialisieren ? - nein, diese Version wird nicht verwendet
//
//noch kleines Problem: bei w.clear() wird GeoWin-redraw Funktion gerufen, was
//in diesem Fall Unsinn ist -> sollte man aendern !!
//
// Loesung: GeoWin::redraw sollte man mehr beeinflussen koennen ...

#define CGAL_CH_NO_POSTCONDITIONS
#define CGAL_PROVIDE_LEDA_RAT_KERNEL_TRAITS_3
#define CGAL_GEOMETRY_EVENTS

#include <CGAL/basic.h>

#if !defined(CGAL_USE_LEDA) || (__LEDA__ < 420)
#include <iostream>

int main(int argc, char *argv[])
{
 std::cout << "No LEDA 4.2 or higher installed!\n";
 std::cout << "A LEDA version >= 4.2 is required !\n";
 return 0;
}
#else 

#include <CGAL/Cartesian.h>
#include <CGAL/leda_rational.h>
#include <CGAL/Kernel_checker.h>
#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>
#include <CEP/Leda_rat_kernel/geowin_leda_rat_kernel.h>
#include <CGAL/ch_selected_extreme_points_2.h>
#include <CGAL/ch_jarvis.h>
#include <CGAL/geowin_support.h>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

//typedef CGAL::Cartesian<leda_rational>    K;

typedef CGAL::leda_rat_kernel_traits      K;
//typedef CGAL::Homogeneous<leda_integer>   K2;
//typedef CGAL::Kernel_checker<K, K2, CGAL::leda_to_cgal_2 > Kernel;

typedef K::Less_rotate_ccw_2              Less_rotate_ccw_2;

typedef K::Point_2                        Point;
typedef K::Segment_2                      Segment;


void new_redraw(leda_window* wp, double x0, double y0, double x1, double y1)
{ }

class geo_hull : public geowin_update<std::list<Point>, std::list<Segment> >
{
public:
 CGAL::event_item less_rotate_ccw_it;
 CGAL::event_item less_xy_it; 
 int less_rotate_ccw_counter;
 int less_xy_counter;
 GeoWin& gw;
 leda_window& w;
 
 std::list<Point>         current_hull;
 leda_point_style              pold;
 const std::list<Point>*  input_set;

 geo_hull(GeoWin& g) : gw(g), w(g.get_window()) 
 { less_rotate_ccw_counter=0; less_xy_counter=0; input_set = NULL; }
 
 virtual ~geo_hull() { }

 // --------------------------------------------------------------------------------------------------
 // event handling ...
 // -------------------------------------------------------------------------------------------------- 

 void user_interaction() { w.read_mouse(); }
 
 void draw_points()
 { leda_point_style ps = w.set_point_style(leda_cross_point);
   std::list<Point>::const_iterator it= input_set->begin();
   for(;it != input_set->end(); it++) w.draw_point(it->to_float());
   w.set_point_style(ps);
 }
 
 void draw_current_hull()
 {
   w.clear();
   draw_points();
   std::list<Point>::const_iterator cit = current_hull.begin();
   Point plast;
   bool first = true;
   for(;cit != current_hull.end(); cit++){
     w.draw_point((*cit).to_float(), leda_green);
     if (! first) w.draw_segment((*cit).to_float(), plast.to_float(), leda_green);
     plast = *cit;
     first= false;
   }
 }
 
 bool less_rotate_ccw_2(const Point& p1, const Point& p2, const Point& p3)
 {
   // switch off the event
   CGAL::disable(less_rotate_ccw_it);
   // compute result ...
   Less_rotate_ccw_2 functor;
   bool result = functor(p1,p2,p3);
   // switch it on again ...
   CGAL::enable(less_rotate_ccw_it);
   return result;
 }
 
 void less_rotate_ccw_occurence(const Point& p1, const Point& p2, const Point& p3)
 {
   less_rotate_ccw_counter++;
   //std::cout << "less_rotate_ccw:" << p1 << " " << p2 << " " << p3 << "\n";
   //store hull ...
   if (current_hull.empty() || p1 != current_hull.front()) current_hull.push_front(p1);
   
   bool result = less_rotate_ccw_2(p1,p2,p3);  
   
   draw_current_hull();
    
   if (result) {
    gw.msg_clear();
    gw.msg_open(leda_string("less_rotate_ccw - new point found ..."));
    w.draw_ray(p1.to_float(), p2.to_float(), leda_blue);
    w.draw_arrow(p2.to_float(), p3.to_float(), leda_black);
    user_interaction();
   }  
 }
 
 void less_xy_occurence(const Point& p1, const Point& p2)
 {
   //std::cout << "less_xy:" << p1 << " " << p2 << "\n";
   bool result = (leda_rat_point::cmp_xy(p1,p2)  <  0);
   
   if (less_xy_counter == 0) // draw start point ...
     w.draw_point(p1.to_float(), leda_red);
   if (result) {
     w.draw_point(p1.to_float(), leda_red);
     gw.msg_clear();
     gw.msg_open(leda_string("new minimum point (xy-lexicographically) found ..."));
     user_interaction();
   }
   less_xy_counter++;
 } 
 
 void init_visualization(const std::list<Point>& L)
 {
   less_rotate_ccw_counter = 0; less_rotate_ccw_counter = 0;   
   less_rotate_ccw_it = CGAL::attach(CGAL::Predicate_leda_rat_less_rotate_ccw_2<K>::ev_leda_rat_point, \
                                     *this, &geo_hull::less_rotate_ccw_occurence);
   less_xy_it         = CGAL::attach(CGAL::Predicate_leda_rat_less_xy_2<K>::ev_leda_rat_point, \
                                     *this, &geo_hull::less_xy_occurence);				     
   w.clear();
   w.set_redraw(new_redraw);
   input_set = &L;
   draw_points();
   current_hull.clear(); 
   pold = w.set_point_style(leda_disc_point);
 }
 
 void reset_visualization() {
   std::cout << "less_rotate_ccw_counter:" << less_rotate_ccw_counter << "\n";
   std::cout << "less_xy_counter        :" << less_xy_counter << "\n";      
   CGAL::detach(less_rotate_ccw_it);
   CGAL::detach(less_xy_it);
   w.set_point_style(pold);
   w.set_redraw(GeoWin::redraw_geowin);
 }
 // --------------------------------------------------------------------------------------------------

 void update(const std::list<Point>& L, std::list<Segment>& Sl)
 {
  Sl.clear();
  std::list<Point> out;  
  //Kernel traits; 
  K traits;
  init_visualization(L);
  CGAL::ch_jarvis(L.begin(),L.end(), std::back_inserter(out), traits);   
  reset_visualization();

  if( out.size() > 1 ) {
    Point pakt,prev,pstart;
    std::list<Point>::const_iterator it=out.begin();
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
  geowin_init_default_type((std::list<Point>*)0, leda_string("point"));
  std::list<Point> L;

  GeoWin GW("Jarvis march"); 
  geo_hull update_obj(GW);
  
  geo_scene my_scene= GW.new_scene(L);  
  geo_scene result  = GW.new_scene(update_obj,my_scene,leda_string("Convex Hull")); 
  GW.set_visible(result,true);
 
  GW.message("Giftwrapping - 1. phase: find point with minimal x-value, 2. phase: wrap the input point set");
  GW.edit(my_scene);
  
  return 0;  
}
#endif
