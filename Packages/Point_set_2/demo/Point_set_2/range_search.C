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
#include <CGAL/leda_rational.h>
#include <CGAL/IO/Window_stream.h>

typedef CGAL::Cartesian<leda_rational>     REP;
typedef CGAL::point_set_traits_2<REP>      TRAITS;
typedef CGAL::Point_set_2<TRAITS>::Edge    Edge;
typedef CGAL::Point_set_2<TRAITS>::Vertex  Vertex;

leda_segment get_seg(const CGAL::Segment_2<REP>& s)
{
  CGAL::Point_2<REP> ps=s.source();
  CGAL::Point_2<REP> pt=s.target();

  leda_rat_point p1(ps.x(),ps.y()), p2(pt.x(),pt.y());
  
  return leda_segment(p1.to_point(),p2.to_point());
}

void output(leda_window& W, const CGAL::Point_set_2<TRAITS>& PSet)
{
  W.clear();
  leda_edge e;
  forall_edges(e,PSet) {
    CGAL::Segment_2<REP> s= PSet.seg(e);
    W << get_seg(s);
  }
}

int main()
{
  int i;
  CGAL::Point_set_2<TRAITS> PSet;
  CGAL::Point_2<REP> pnew;
  leda_window W(600,500,leda_string("Range searches"));  
  CGAL::cgalize( W);

  W.init(-500,500,-400);
  W.display(100,100);

  for (i=0; i<30; i++) {
    W >> pnew;
    PSet.insert(pnew);
    output(W,PSet);    
  }

  for (i=0; i<5; i++) {
    W >> pnew;
    Vertex v = PSet.nearest_neighbor(pnew);
    std::cout << "delete!\n"; std::cout.flush();
    PSet.del(v);
    output(W,PSet);    
  }  

  std::cout << "range search for circle !\n";  
  CGAL::Circle_2<REP> rc;
  W >> rc;

  std::list<Vertex> LV;
  PSet.range_search(rc,std::back_inserter(LV));

  std::list<Vertex>::const_iterator it;
  for (it=LV.begin();it != LV.end(); it++)
      W << CGAL::GREEN << PSet.pos(*it) << CGAL::BLACK;      
 
  std::cout << "range search for triangle !\n";    
  CGAL::Point_2<REP> pt1,pt2,pt3,pt4;
  W >> pt1;
  W >> pt2;
  W >> pt3;
  
  LV.clear();
  PSet.range_search(pt1,pt2,pt3,std::back_inserter(LV));
  for (it=LV.begin();it != LV.end(); it++)
      W << CGAL::BLUE << PSet.pos(*it) << CGAL::BLACK;
        
  LV.clear();
 
  std::cout << "range search for iso rectangle !\n";
  W >> pt1; // lower left
  W >> pt3; // upper right 
  pt2 = CGAL::Point_2<REP>(pt3.x(),pt1.y());
  pt4 = CGAL::Point_2<REP>(pt1.x(),pt3.y());
  
  PSet.range_search(pt1,pt2,pt3,pt4,std::back_inserter(LV));
  for (it=LV.begin();it != LV.end(); it++)
      W << CGAL::RED << PSet.pos(*it) << CGAL::BLACK;   

  leda_list<Edge> El = PSet.minimum_spanning_tree();
  Edge eakt;
  
  std::cout << "minimum spanning tree !\n";
  forall(eakt,El) W << CGAL::VIOLET << PSet.seg(eakt);  

  W.read_mouse();

  return 0;
}

#endif
