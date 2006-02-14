
#include <CGAL/basic.h>
#include <CGAL/Triangulation_2.h>
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


typedef CGAL::Triangulation_2<K>                               Triang_2;
typedef Triang_2::Edge                                         Edge;
typedef Triang_2::Edge_iterator                                Edge_iterator;


void draw(leda_window& W, Triang_2& Tr ) 
{
  W.clear();
  W.set_color(leda_red);

  Edge_iterator eit = Tr.edges_begin();
  Edge_iterator beyond = Tr.edges_end();   
  Edge eact;

  while (eit != beyond) {
       eact = *eit;   
       W << (Tr.segment(eact));                    
       ++eit;  
  }    
}


int main()
{
  Triang_2 T;
  std::list<Point> input;

  leda_window win(500,500,"Delaunay triangulation of points");
  win.init(0,500,0);
  
  win.display();
  
  Point act;
  
  while(win >> act){
   win << act;
   input.push_back(act);
   T.insert(act);
   draw(win, T);
  }
  return 0;  
}
