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

typedef CGAL::Cartesian<leda_rational>    REP;
typedef CGAL::point_set_traits_2<REP>     TRAITS;
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
  CGAL::Point_set_2<TRAITS> PSet;

  leda_window W(600,500, leda_string("Finding nearest neighbors"));  
  CGAL::cgalize( W);

  W.init(-500,500,-400);
  W.display(100,100);
  
  CGAL::Point_2<REP> actual;
  int i;
  
  for (i=0; i<25; i++){
   W >> actual;
   PSet.insert(actual);
   output(W,PSet);
  }  

  // nearest neighbor ...  
  for (i=0; i<5; i++){
    W >> actual;
    Vertex v = PSet.nearest_neighbor(actual);
    
    CGAL::Segment_2<REP> my_seg(actual,PSet.pos(v));
    
    W << CGAL::RED << PSet.pos(v) << CGAL::BLACK;
    W << CGAL::BLUE << my_seg << CGAL::BLACK;
  }
  
  // k nearest neighbors ...
  std::list<Vertex> L;
  std::list<Vertex>::const_iterator it;
  
  for (i=0; i<5; i++){
    L.clear();
    W >> actual;
    PSet.nearest_neighbors(actual,5, std::back_inserter(L));
    std::cout << "actual point: " << actual << "\n";
    
    W.clear();
    output(W,PSet);
    W << CGAL::RED << actual << CGAL::BLACK;
    
    for (it=L.begin();it != L.end(); it++){
      W << CGAL::GREEN << PSet.pos(*it) << CGAL::BLACK;      
      std::cout << PSet.pos(*it) << "\n";
    }
    std::cout << "\n";
  }  

  W.read_mouse();

  return 0;
}
#endif
