#define LEDA_NO_MIN_MAX_TEMPL

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

#define CGAL_ALPHA_WINDOW_STREAM

#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>
#include <CEP/Leda_rat_kernel/geowin_leda_rat_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/geowin_support.h>

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

typedef CGAL::leda_rat_kernel_traits                K;

typedef K::Point_2                                  Point;
typedef K::Segment_2                                Segment;
typedef CGAL::Alpha_shape_vertex_base_2<K>          Vb;

typedef CGAL::Triangulation_face_base_2<K>          Df;
typedef CGAL::Alpha_shape_face_base_2<K, Df>        Fb;
typedef CGAL::Triangulation_default_data_structure_2<K,Vb,Fb> Tds;
typedef CGAL::Delaunay_triangulation_2<K,Tds>       Triangulation_2;
typedef CGAL::Alpha_shape_2<Triangulation_2>        Alpha_shape_2;

typedef Alpha_shape_2::Edge                          Edge;
typedef Alpha_shape_2::Vertex_handle                 Vertex_handle;
typedef Alpha_shape_2::Edge_iterator                 Edge_iterator;


class geo_delaunay_triang : public geowin_update<std::list<Point>, std::list<Segment> >,
                            public geowin_redraw
{
 Triangulation_2 dt;  

public:
 virtual ~geo_delaunay_triang() {}


 bool insert(const Point& p)
 {
  dt.insert(p);
  return true; 
 }
 
 bool del(const Point& p)
 {
  Triangulation_2::Vertex_handle vh = dt.nearest_vertex(p);
  dt.remove(vh);
  return true;  
 } 

 void draw(leda_window& W,leda_color c1,leda_color c2,double x1,double y1,double x2,double y2)
 {
  Edge_iterator eit = dt.edges_begin();
  Edge_iterator beyond = dt.edges_end();   
  Edge eact;

  while (eit != beyond) {
       eact = *eit;   
       W.draw_segment(dt.segment(eact).to_float(),c1);                    
       ++eit;  
  }    
 }

 void update(const std::list<Point>& L, std::list<Segment>&)
 {
  dt.clear();       
  dt.insert(L.begin(),L.end());
 }
};


int alpha_index;
bool   reg;

class geo_alpha : public geowin_update<std::list<Point>, std::list<Segment> >,
                  public geowin_redraw
{
public:
 Alpha_shape_2 A;
 
 virtual ~geo_alpha() { }

 void draw(leda_window& W, leda_color c1, leda_color c2,double x1,double y1,double x2,double y2)
 {  W.set_color(c1); 
    W << A; 
 }  
 
 void update(const std::list<Point>& L, std::list<Segment>& Sl)
 {
  A.make_alpha_shape(L.begin(), L.end());
  A.set_alpha((double)alpha_index);
  
  if (reg) A.set_mode(Alpha_shape_2::REGULARIZED);
  else     A.set_mode(Alpha_shape_2::GENERAL);
 }
};

geo_scene res2, res3;
GeoWin* gwin;

void fcn(int val){
   res2->update(); res3->update();
   gwin->redraw();
}

int main()
{
  geowin_init_default_type((std::list<Point>*)0, leda_string("LEDA-rat_point"));

  std::list<Point> L;

  GeoWin GW("Alpha shapes");
  GW.add_help_text(leda_string("CGAL_alpha_shape_2"));
  
  gwin = &GW;

  geo_scene my_scene= GW.new_scene(L);  
  GW.set_point_style(my_scene, leda_disc_point);
 
  geo_delaunay_triang deltria;  
  res2 = GW.new_scene(deltria, deltria, my_scene, leda_string("Delaunay Triangulation"));
  GW.set_color(res2, leda_yellow);

  geo_alpha alpha;  
  res3 = GW.new_scene(alpha, alpha, my_scene, leda_string("Alpha shape"));
  GW.set_line_width(res3, 3);
  GW.set_color(res3, leda_green2);

  GW.set_all_visible(true);
  
  GW.init_menu();
  alpha_index = 1000;
  GW.get_window().int_item("Alpha:",alpha_index,1,10000,fcn); 
  reg = false;
  GW.get_window().bool_item("Regularized:",reg);
  GW.edit(my_scene);
  
  return 0;  
}

#endif
