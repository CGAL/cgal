#include <CGAL/config.h>
#include <list>
#include <LEDA/window.h>
#include <CGAL/Cartesian.h>

#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_euclidean_traits_2.h>

#include <CGAL/Point_set_2_tb.h>
#include <CGAL/IO/Window_stream.h>

typedef CGAL::Cartesian<double>    REP;
typedef CGAL::Triangulation_euclidean_traits_2<REP> Gt;
typedef CGAL::Triangulation_vertex_base_2<Gt> Vb;
typedef CGAL::Triangulation_face_base_2<Gt>  Fb;
typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;

typedef CGAL::point_set_traits_2<REP>     TRAITS;
typedef CGAL::Point_set_2_tb<TRAITS,Gt,Tds>::Edge    Edge;
typedef CGAL::Point_set_2_tb<TRAITS,Gt,Tds>::Edge_iterator  Edge_iterator;
typedef CGAL::Point_set_2_tb<TRAITS,Gt,Tds>::Vertex_handle  Vertex_handle;
typedef CGAL::Point_set_2_tb<TRAITS,Gt,Tds>::Vertex  Vertex;

leda_segment get_seg(const CGAL::Segment_2<REP>& s)
{
  CGAL::Point_2<REP> ps=s.source();
  CGAL::Point_2<REP> pt=s.target();

  leda_point p1(ps.x(),ps.y()), p2(pt.x(),pt.y());
  
  return leda_segment(p1, p2);
}

void output(leda_window& W, const CGAL::Point_set_2_tb<TRAITS,Gt,Tds>& PSet)
{
  W.clear();
  Edge e;
  Edge_iterator eit = PSet.finite_edges_begin();
  
  for(;eit != PSet.finite_edges_end(); eit++) {
    e = *eit;
    CGAL::Segment_2<REP> s= PSet.seg(e);
    W << get_seg(s);
  }
}

int main()
{
  CGAL::Point_set_2_tb<TRAITS,Gt,Tds> PSet;

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
    Vertex_handle v = PSet.nearest_neighbor(actual);
    
    CGAL::Segment_2<REP> my_seg(actual,PSet.pos(v));
    
    W << CGAL::RED << PSet.pos(v) << CGAL::BLACK;
    W << CGAL::BLUE << my_seg << CGAL::BLACK;
  }
  
  // k nearest neighbors ...
  std::list<Vertex_handle> L;
  std::list<Vertex_handle>::const_iterator it;
  
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

  return 1;
}

