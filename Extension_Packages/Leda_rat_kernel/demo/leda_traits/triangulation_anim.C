#define CGAL_GEOMETRY_EVENTS

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

#include <CGAL/Triangulation_2.h>
#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>
#include <CEP/Leda_rat_kernel/geowin_leda_rat_kernel.h>
#include <CGAL/geowin_support.h>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

typedef CGAL::leda_rat_kernel_traits    K;
typedef K::Point_2                      Point;
typedef K::Segment_2                    Segment;
typedef K::Ray_2                        Ray;

typedef CGAL::Triangulation_2<K>        Triangulation_2;
typedef Triangulation_2::Edge           Edge;
typedef Triangulation_2::Edge_iterator  Edge_iterator;

void new_redraw(window* wp, double x0, double y0, double x1, double y1)
{ }

struct geo_triang : public geowin_update<std::list<Point>, std::list<Segment> >
{
 Triangulation_2 tr; 
 const std::list<Point>*  input_set;

 CGAL::event_item orientation_it;
 
 GeoWin& gw;
 window& w; 

 geo_triang(GeoWin& g) : gw(g), w(g.get_window()) 
 { } 
 
 virtual ~geo_triang() { }
 
 void user_interaction() { w.read_mouse(); } 
 
 void draw_points()
 { 
   std::list<Point>::const_iterator it= input_set->begin();
   for(;it != input_set->end(); it++) w.draw_point(it->to_float());
 }
 
 void draw_triangulation(const Point& p1, const Point& p2, const Point& p3)
 {
   w.clear();
   draw_points();
   Edge_iterator eit = tr.edges_begin();
   for(;eit != tr.edges_end();eit++) w.draw_segment(tr.segment(*eit).to_float(),leda_blue2);      
   w.draw_arrow(p1.to_float(),p2.to_float(),green);
   w.draw_arrow(p2.to_float(),p3.to_float(),green);
   user_interaction();           
 }
 
 void orientation_occurence(const Point& p1, const Point& p2, const Point& p3)
 {
   CGAL::disable(orientation_it);
   std::cout << "orientation ...\n"; 
   std::cout << p1 << " " << p2 << " " << p3 << "\n";
   draw_triangulation(p1,p2,p3);
   CGAL::enable(orientation_it);
 }

 // --------------------------------------------------------------------------------------------------  
 void init_visualization(const std::list<Point>& L)
 {
   orientation_it   = CGAL::attach(CGAL::Predicate_leda_rat_orientation_2<K>::ev_leda_rat_point, \
                                   *this, &geo_triang::orientation_occurence);	
				     				     			     
   w.clear();
   input_set = &L;
   w.set_redraw(new_redraw);
 }
 
 void reset_visualization() {
   CGAL::detach(orientation_it);
   w.set_redraw(GeoWin::redraw_geowin);
 }
 // -------------------------------------------------------------------------------------------------- 

 void update(const std::list<Point>& L, std::list<Segment>& Sl)
 {
  tr.clear();    
  Sl.clear();             
    
  init_visualization(L);            
  tr.insert(L.begin(), L.end());
  reset_visualization();

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

int main()
{
  geowin_init_default_type((std::list<Point>*)0, leda_string("LEDA-rat_point"));

  std::list<Point> L;
  GeoWin GW("CGAL - Triangulation demo using LEDA kernel");

  geo_scene my_scene= GW.new_scene(L);  
  GW.set_point_style(my_scene, leda_disc_point);

  geo_triang triangulate(GW);
  geo_scene res1 = GW.new_scene(triangulate ,my_scene , leda_string("Triangulation"));
  GW.set_color(res1, leda_blue2);
  GW.set_all_visible(true);
  
  GW.edit(my_scene);
  
  return 0;  
}

#endif
