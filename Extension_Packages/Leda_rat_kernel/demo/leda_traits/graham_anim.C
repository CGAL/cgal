// Graham scan - wir haben nach der Sortierung zwei Phasen:
// 1. untere Huelle wird gebaut
// 2. obere Huelle wird gebaut

#define CGAL_CH_NO_POSTCONDITIONS
#define CGAL_PROVIDE_LEDA_RAT_KERNEL_TRAITS_3
#define CGAL_GEOMETRY_EVENTS

#include <CGAL/basic.h>

#if !defined(CGAL_USE_LEDA) || (__LEDA__ < 430)
#include <iostream>

int main(int argc, char *argv[])
{
 std::cout << "No LEDA 4.3 or higher installed!\n";
 std::cout << "A LEDA version >= 4.3 is required !\n";
 return 0;
}
#else 

#include <CGAL/Cartesian.h>
#include <CGAL/leda_rational.h>
#include <CGAL/Kernel_checker.h>
#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>
#include <CEP/Leda_rat_kernel/geowin_leda_rat_kernel.h>
#include <CGAL/ch_graham_andrew.h>
#include <CGAL/geowin_support.h>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

//typedef CGAL::Cartesian<leda_rational>    K;

typedef CGAL::leda_rat_kernel_traits      K;
//typedef CGAL::Homogeneous<leda_integer>   K2;
//typedef CGAL::Kernel_checker<K, K2, CGAL::leda_to_cgal_2 > Kernel;

typedef K::Point_2                        Point;
typedef K::Segment_2                      Segment;


void new_redraw(window* wp, double x0, double y0, double x1, double y1)
{ }

class geo_hull : public geowin_update<std::list<Point>, std::list<Segment> >
{
public:
 CGAL::event_item left_turn_it;
 CGAL::event_item less_xy_it; 
 GeoWin& gw;
 window& w;
 
 std::list<Point>         current_hull;
 Point                    p_left, p_right;
 std::list<Point>         out;  
 
 point_style              pold;
 const std::list<Point>*  input_set;
 bool                     first_scan; 
 int                      algorithm_phase; // 0 ... sorting; 1 ... lower hull; 2 - upper hull

 geo_hull(GeoWin& g) : gw(g), w(g.get_window()) 
 { input_set = NULL; }
 
 virtual ~geo_hull() { }

 // --------------------------------------------------------------------------------------------------
 // event handling ...
 // -------------------------------------------------------------------------------------------------- 

 void user_interaction() { w.read_mouse(); }
 
 void draw_points()
 { point_style ps = w.set_point_style(cross_point);
   std::list<Point>::const_iterator it= input_set->begin();
   for(;it != input_set->end(); it++) w.draw_point(it->to_float());
   w.set_point_style(ps);
 }
 
 void draw_current_hull()
 {
   w.clear();
   draw_points();
   std::list<Point>::const_iterator cit = out.begin();
   Point plast;
   bool first = true;
   for(;cit != out.end(); cit++){
     w.draw_point((*cit).to_float(), green);
     if (! first) w.draw_segment((*cit).to_float(), plast.to_float(), black);
     plast = *cit;
     first= false;
   }
   w.draw_segment(p_right.to_float(), plast.to_float(), black);
   w.draw_line(p_left.to_float(), p_right.to_float(), blue);
 }
 
 void left_turn_occurence(const Point& p1, const Point& p2, const Point& p3)
 {
   // compute result ...
   bool result = leda::left_turn(p1,p2,p3);   
   std::cout << "left_turn:" << p1 << " " << p2 << " " << p3 << " - result:" << result << "\n";
   
   if (first_scan) { // first left_turn call ...
     first_scan = false;
     algorithm_phase = 1;
     std::cout << "construct lower hull ...\n";
     p_right = p1;
     p_left  = p2;
     
     // line separating upper and lower hull ...
     w.draw_line(p1.to_float(), p2.to_float(), blue);
   }
   else {
     if (algorithm_phase == 1 && identical(p2, p_right)) {
        algorithm_phase =2;
	draw_current_hull();
	std::cout << "construct upper hull ...\n";
     }
   }
   
   if (result) {
    w.draw_segment(p1.to_float(), p2.to_float(), black);
    w.draw_arrow(p2.to_float(), p3.to_float(), black);
   }
   else {
    w.draw_arrow(p1.to_float(), p2.to_float(), green);
    w.draw_arrow(p2.to_float(), p3.to_float(), green);
    w.draw_point(p3.to_float(), red);
   }
   user_interaction();
 }
 
 void less_xy_occurence(const Point& p1, const Point& p2)
 { } 
 
 void init_visualization(const std::list<Point>& L)
 {
   left_turn_it       = CGAL::attach(CGAL::Predicate_leda_rat_leftturn_2::ev_leda_rat_point, \
                                     *this, &geo_hull::left_turn_occurence);
   less_xy_it         = CGAL::attach(CGAL::Predicate_leda_rat_less_xy_2::ev_leda_rat_point, \
                                     *this, &geo_hull::less_xy_occurence);				     
   w.clear();
   w.set_redraw(new_redraw);
   input_set = &L;
   draw_points();
   first_scan = true;
   algorithm_phase = 0;
   current_hull.clear(); 
   pold = w.set_point_style(disc_point);
 }
 
 void reset_visualization() {
   CGAL::detach(left_turn_it);
   CGAL::detach(less_xy_it);
   w.set_point_style(pold);
   w.set_redraw(GeoWin::redraw_geowin);
 }
 // --------------------------------------------------------------------------------------------------

 void update(const std::list<Point>& L, std::list<Segment>& Sl)
 {
  Sl.clear(); 
  out.clear();
  //Kernel traits; 
  K traits;
  init_visualization(L);
  CGAL::ch_graham_andrew(L.begin(),L.end(), std::back_inserter(out), traits);   
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

  GeoWin GW; 
  geo_hull update_obj(GW);
  
  geo_scene my_scene= GW.new_scene(L);  
  geo_scene result  = GW.new_scene(update_obj,my_scene,leda_string("Convex Hull")); 
  GW.set_visible(result,true);
 
  GW.message("Graham-Andrew algorithm");
  GW.edit(my_scene);
  
  return 0;  
}
#endif
