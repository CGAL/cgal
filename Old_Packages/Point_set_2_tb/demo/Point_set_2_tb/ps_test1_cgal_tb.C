#include <CGAL/config.h>
#include <list>
#include <LEDA/window.h>
#include <CGAL/Cartesian.h>

#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_euclidean_traits_2.h>

#include <CGAL/Point_set_2_tb.h>
#include <CGAL/IO/Window_stream.h>
#include <CGAL/leda_rational.h>
#include <LEDA/rat_point.h>

//typedef CGAL::Cartesian<double>          REP;
typedef CGAL::Cartesian<leda_rational>          REP;
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

  leda_rat_point p1(ps.x(),ps.y()), p2(pt.x(),pt.y());
  
  return leda_segment(p1.to_float(), p2.to_float());
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
  CGAL::Point_set_2_tb<TRAITS,Gt,Tds> PSR_rat;

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
    Vertex_handle v = PSR_rat.nearest_neighbor(pnew);
    PSR_rat.del(v);
    output(W,PSR_rat);    
  }  

  std::cout << "circular range search !\n";
  
  CGAL::Circle_2<REP> rc;
  W >> rc;

  std::list<Vertex_handle> LV;
  std::list<Vertex_handle>::const_iterator vit;
  
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
  
  std::list<Vertex_handle>  outlist;
  
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

  W.read_mouse();

  return 1;
}

