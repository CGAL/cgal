#include <CGAL/basic.h>

#if !defined(CGAL_USE_LEDA) || (__LEDA__ < 400)
#include <iostream>

int main(int argc, char *argv[])
{
 std::cout << "No LEDA 4.0 or higher installed!\n";
 std::cout << "A LEDA version >= 4.0 is required !\n";
 return 0;
}
#else 


// Pointset demo using LEDA lists as containers;
// LEDA 4.0 or higher should be used

#include <CGAL/config.h>
#include <list>
#include <LEDA/window.h>
#include <LEDA/rat_window.h>
#include <CGAL/Point_set_2.h>


typedef CGAL::point_set_leda_rat_traits_2  TRAITS;
typedef CGAL::Point_set_2<TRAITS>::Edge    Edge;
typedef CGAL::Point_set_2<TRAITS>::Vertex  Vertex;

void output(leda_window& W, const CGAL::Point_set_2<TRAITS>& PSR_rat)
{
  W.clear();
  leda_edge e;
  forall_edges(e,PSR_rat) W << leda_segment(PSR_rat.seg(e).to_segment());
}

int main()
{
  int i;
  CGAL::Point_set_2<TRAITS> PSR_rat;
  leda_window W(500,400);
  W.init(-500,500,-400);
  W.display();

  for (i=0; i<30; i++) {
    leda_rat_point pnew;
    W >> pnew;
    PSR_rat.insert(pnew);
    output(W,PSR_rat);    
  }

  for (i=0; i<3; i++) {
    leda_rat_point pnew;
    W >> pnew;
    Vertex v = PSR_rat.nearest_neighbor(pnew);
    PSR_rat.del(v);
    output(W,PSR_rat);    
  }  

  leda_rat_circle rc;
  W >> rc; W << rc;
  
#if (__LEDA__ >= 400)
  W.set_point_style(leda_disc_point);
#endif

  leda_list<Vertex> LV = PSR_rat.range_search(rc);
  Vertex v;  
  W.set_color(leda_red);
  forall(v,LV) W << PSR_rat.pos(v);
  
  std::cout << "circular range_search - found " << LV.size() << " vertices!\n";
  
  std::cout << "range search for triangle !\n";    
  leda_rat_point pt1,pt2,pt3,pt4;
  W >> pt1; W << pt1;
  W >> pt2; W << pt2;
  W >> pt3; W << pt3;
  
  LV.clear();
  PSR_rat.range_search(pt1,pt2,pt3,std::back_inserter(LV));
  leda_list<Vertex>::iterator it;
  W.set_color(leda_green);
  
  for (it=LV.begin();it != LV.end(); it++)
      W <<  PSR_rat.pos(*it);
      
  std::cout << "triangular range_search - found " << LV.size() << " vertices!\n";
  
  LV.clear();
  
  W >> pt1; W << pt1;
  W >> pt3; W << pt3;
  
  pt2 = leda_rat_point(pt3.xcoord(),pt1.ycoord());
  pt4 = leda_rat_point(pt1.xcoord(),pt3.ycoord());
  
  W.set_color(leda_orange);
  PSR_rat.range_search(pt1,pt2,pt3,pt4,std::back_inserter(LV));
  for (it=LV.begin();it != LV.end(); it++)
      W << PSR_rat.pos(*it);     
      
  std::cout << "rectangular range_search - found " << LV.size() << " vertices!\n";
          
  leda_list<Edge> El = PSR_rat.minimum_spanning_tree();
  
  std::cout << "\nMST: " << El.size() << " vertices\n";

  W.read_mouse();
  return 0;
}
#endif
