#include <CGAL/config.h>
#include <list>
#include <LEDA/window.h>
#include <CGAL/Cartesian.h>

#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_euclidean_traits_2.h>

#include <CGAL/Point_set_2_tb.h>
#include <CGAL/IO/Window_stream.h>

typedef CGAL::Cartesian<double>     REP;
typedef CGAL::Triangulation_euclidean_traits_2<REP> Gt;
typedef CGAL::Triangulation_vertex_base_2<Gt> Vb;
typedef CGAL::Triangulation_face_base_2<Gt>  Fb;
typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;


typedef CGAL::point_set_traits_2<REP>      TRAITS;
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
  int i;
  CGAL::Point_set_2_tb<TRAITS,Gt,Tds> PSet;
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
    Vertex_handle v = PSet.nearest_neighbor(pnew);
    std::cout << "delete!\n"; std::cout.flush();
    PSet.del(v);
    output(W,PSet);    
  }  


  std::cout << "range search for circle !\n";  
  CGAL::Circle_2<REP> rc;
  W >> rc;

  std::list<Vertex_handle> LV;
  PSet.range_search(rc,std::back_inserter(LV));

  std::list<Vertex_handle>::const_iterator it;
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

  W.read_mouse();

  return 1;
}

