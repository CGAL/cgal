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
#include <CGAL/leda_integer.h>
#include <CGAL/IO/Window_stream.h>

typedef CGAL::Homogeneous<leda_integer>      REP;
typedef CGAL::point_set_traits_2<REP>        TRAITS;
typedef CGAL::Point_set_2<TRAITS>::Edge       Edge;
typedef CGAL::Point_set_2<TRAITS>::Vertex     Vertex;

leda_segment get_seg(const CGAL::Segment_2<REP>& s)
{
  CGAL::Point_2<REP> ps=s.source();
  CGAL::Point_2<REP> pt=s.target();
  
  leda_rat_point p1(ps.hx()/ps.hw(),ps.hy()/ps.hw());
  leda_rat_point p2(pt.hx()/ps.hw(),pt.hy()/ps.hw());
  
  return leda_segment(p1.to_point(),p2.to_point());
}

void output(leda_window& W, const CGAL::Point_set_2<TRAITS>& PSet)
{
  W.clear();
  leda_edge e;
  forall_edges(e,PSet) {
    CGAL::Segment_2<REP> s = PSet.seg(e);
    W << get_seg(s);
  }
}

int main()
{
  int i;
  CGAL::Point_set_2<TRAITS> PSet;
  CGAL::Point_2<REP> actual;
  
  leda_window W(600,500,leda_string("Input vertices of the Delaunay Triangulation!"));  
  CGAL::cgalize( W);

  W.init(-500,500,-400);
  W.display(100,100);

  // read in 20 points ...
  for (i=0; i<20; i++){
   W >> actual;
   PSet.insert(actual);
   output(W,PSet);
  }

  // locate ...  
  for (i=0; i<5; i++){
    W >> actual;
    Edge e = PSet.locate(actual);
    
    CGAL::Segment_2<REP> my_seg = PSet.seg(e);
 
    W << CGAL::BLUE << my_seg << CGAL::BLACK;
  }

  W.read_mouse();

  return 0;
}

#endif
