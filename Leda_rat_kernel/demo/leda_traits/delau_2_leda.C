
#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/IO/Window_stream.h>

#include <LEDA/window.h>
#include <LEDA/rat_window.h>
#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>

#include <list> 

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif


typedef CGAL::leda_rat_kernel_traits                           K;
typedef K::Point_2                                             Point;

typedef CGAL::Delaunay_triangulation_2<K>                      Delaunay_triang_2;
typedef Delaunay_triang_2::Edge                                Edge;
typedef Delaunay_triang_2::Edge_iterator                       Edge_iterator;

Delaunay_triang_2 dt;



void draw(leda_window& W, Delaunay_triang_2& D ) 
{
  W.clear();
  W.set_color(leda_red);

  Edge_iterator eit = D.edges_begin();
  Edge_iterator beyond = D.edges_end();   
  Edge eact;

  while (eit != beyond) {
       eact = *eit;   
       W << (D.segment(eact));                    
       ++eit;  
  }    
}


int main()
{
  std::list<Point> input;

  leda_window win(500,500,"Delaunay triangulation of points");
  win.init(0,500,0);
  
  win.display();
  
  Point act;
  
  while(win >> act){
   win << act;
   input.push_back(act);
   dt.insert(act);
   draw(win, dt);
  }
  return 0;  
}
