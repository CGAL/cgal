#include <CGAL/config.h>
#include <list>
#include <LEDA/window.h>
#include <LEDA/rat_window.h>
#include <CGAL/Cartesian.h>

#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_euclidean_leda_traits_2.h>

#include <CGAL/point_set_leda_traits_2.h>
#include <CGAL/Point_set_2_tb.h>
#include <CGAL/IO/Window_stream.h>
#include <LEDA/point.h>

//typedef CGAL::Triangulation_euclidean_leda_rat_traits_2 Gt;
typedef CGAL::Triangulation_euclidean_leda_float_traits_2 Gt;
typedef CGAL::Triangulation_vertex_base_2<Gt> Vb;
typedef CGAL::Triangulation_face_base_2<Gt>  Fb;
typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;

//typedef CGAL::point_set_leda_rat_traits_2      TRAITS;
typedef CGAL::point_set_leda_float_traits_2      TRAITS;
typedef CGAL::Point_set_2_tb<TRAITS,Gt,Tds>::Edge    Edge;
typedef CGAL::Point_set_2_tb<TRAITS,Gt,Tds>::Edge_iterator  Edge_iterator;
typedef CGAL::Point_set_2_tb<TRAITS,Gt,Tds>::Vertex_handle  Vertex_handle;
typedef CGAL::Point_set_2_tb<TRAITS,Gt,Tds>::Vertex  Vertex;
typedef Gt::Segment_2 Segment;
typedef Gt::Point_2 Point;
//typedef leda_rat_circle Circle;
typedef leda_circle Circle;


void output(leda_window& W, const CGAL::Point_set_2_tb<TRAITS,Gt,Tds>& PSet)
{
  W.clear();
  Edge e;
  Edge_iterator eit = PSet.finite_edges_begin();
  
  for(;eit != PSet.finite_edges_end(); eit++) {
    e = *eit;
    Segment s= PSet.seg(e);
    W << s.to_float();
  }
}  

int main()
{
  int i;
  CGAL::Point_set_2_tb<TRAITS,Gt,Tds> PSR_rat;
  leda_window W(500,400);
  W.init(-500,500,-400);
  W.display();

  for (i=0; i<30; i++) {
    Point pnew;
    W >> pnew;
    PSR_rat.insert(pnew);
    output(W,PSR_rat);    
  }

  for (i=0; i<3; i++) {
    Point pnew;
    W >> pnew;
    Vertex_handle v = PSR_rat.nearest_neighbor(pnew);
    PSR_rat.del(v);
    output(W,PSR_rat);    
  }  

  Circle rc;
  W >> rc; W << rc;
  
#if (__LEDA__ >= 400)
  W.set_point_style(leda_disc_point);
#endif

  std::list<Vertex_handle> LV;
  std::list<Vertex_handle>::const_iterator vit;
  
  PSR_rat.range_search(rc,std::back_inserter(LV));
  
  W.set_color(leda_red);
  for(vit=LV.begin(); vit!=LV.end(); vit++){
    W << PSR_rat.pos(*vit);
  }
  W.set_color(leda_black);
  
  std::cout << "circular range_search - found " << LV.size() << " vertices!\n";

  std::cout << "range search for triangle !\n";    
  Point pt1,pt2,pt3,pt4;
  W >> pt1; W << pt1;
  W >> pt2; W << pt2;
  W >> pt3; W << pt3;
  
  LV.clear();
  PSR_rat.range_search(pt1,pt2,pt3,std::back_inserter(LV));
  std::list<Vertex_handle>::const_iterator it;
  W.set_color(leda_green);
  
  for (it=LV.begin();it != LV.end(); it++)
      W <<  PSR_rat.pos(*it);
      
  std::cout << "triangular range_search - found " << LV.size() << " vertices!\n";
  
  LV.clear();
  
  W >> pt1; W << pt1;
  W >> pt3; W << pt3;
  
  pt2 = Point(pt3.xcoord(),pt1.ycoord());
  pt4 = Point(pt1.xcoord(),pt3.ycoord());
  
  W.set_color(leda_orange);
  PSR_rat.range_search(pt1,pt2,pt3,pt4,std::back_inserter(LV));
  for (it=LV.begin();it != LV.end(); it++)
      W << PSR_rat.pos(*it);     
      
  std::cout << "rectangular range_search - found " << LV.size() << " vertices!\n";

  W.read_mouse();
  return 0;
}

