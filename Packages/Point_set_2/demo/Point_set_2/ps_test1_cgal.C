#include <CGAL/basic.h>

#if !defined(CGAL_USE_LEDA)
#include <iostream>

int main(int argc, char *argv[])
{
 std::cout << "No LEDA installed!\n";
 return 0;
}
#else 

#include <CGAL/config.h>
#include <list>
#include <LEDA/window.h>
#include <CGAL/Point_set_2.h>
#include <CGAL/IO/Window_stream.h>

typedef CGAL::Cartesian<double>          REP;
typedef CGAL::point_set_traits_2<REP>      TRAITS;
typedef CGAL::Point_set_2<TRAITS>::Edge    Edge;
typedef CGAL::Point_set_2<TRAITS>::Vertex  Vertex;

leda_segment get_seg(const CGAL::Segment_2<REP>& s)
{
  CGAL::Point_2<REP> ps=s.source();
  CGAL::Point_2<REP> pt=s.target();

  leda_point p1(ps.x(),ps.y()), p2(pt.x(),pt.y());
  return leda_segment(p1,p2);
}

void output(leda_window& W, const CGAL::Point_set_2<TRAITS>& PSR_rat)
{
  W.clear();
  leda_edge e;
  forall_edges(e,PSR_rat) {
    CGAL::Segment_2<REP> s= PSR_rat.seg(e);
    W << get_seg(s);
  }
}

int main()
{
  int i;
  CGAL::Point_set_2<TRAITS> PSR_rat;

  leda_window W(600,500);  
  CGAL::cgalize( W);

  W.init(-500,500,-400);
  W.display(100,100);
  
  output(W,PSR_rat);

  for (i=0; i<25; i++) {
    CGAL::Point_2<REP> pnew;
    W >> pnew;
    PSR_rat.insert(pnew);
    output(W,PSR_rat);    
  }

  for (i=0; i<2; i++) {
    CGAL::Point_2<REP> pnew;
    W >> pnew;
    Vertex v = PSR_rat.nearest_neighbor(pnew);
    PSR_rat.del(v);
    output(W,PSR_rat);    
  }  

  std::cout << "circular range search !\n";
  
  CGAL::Circle_2<REP> rc;
  W >> rc;

  std::list<Vertex> LV;
  std::list<Vertex>::const_iterator vit;
  
  PSR_rat.range_search(rc,std::back_inserter(LV));
  
  W.set_color(leda_red);
  for(vit=LV.begin(); vit!=LV.end(); vit++){
    W << PSR_rat.pos(*vit);
  }
  W.set_color(leda_black);
   
  std::cout << LV.size() << "\n";
 
  std::cout << "triangular range search !\n";
    
  CGAL::Point_2<REP> pt1,pt2,pt3,pt4;
  W >> pt1;
  W >> pt2;
  W >> pt3;
  
  std::list<Vertex>  outlist;
  
  PSR_rat.range_search(pt1,pt2,pt3,std::back_inserter(outlist));
  W.set_color(leda_green);
  for(vit=outlist.begin(); vit!=outlist.end(); vit++){
    W << PSR_rat.pos(*vit);
  }  
  W.set_color(leda_black);
  
  std::cout << outlist.size() << "\n";  
  
  outlist.clear();
 
  std::cout << "rectangular range search !\n";
  W >> pt1; W >> pt3; 
 
  pt2 = CGAL::Point_2<REP>(pt3.x(),pt1.y());
  pt4 = CGAL::Point_2<REP>(pt1.x(),pt3.y());
  
  PSR_rat.range_search(pt1,pt2,pt3,pt4,std::back_inserter(outlist));
  W.set_color(leda_yellow);
  for(vit=outlist.begin(); vit!=outlist.end(); vit++){
    W << PSR_rat.pos(*vit);
  }  
  W.set_color(leda_black);
  std::cout << outlist.size() << "\n";    

  leda_list<Edge> El = PSR_rat.minimum_spanning_tree();
  std::cout << "MST:\n" << El.size() << " vertices\n";

  W.read_mouse();

  return 0;
}

#endif
