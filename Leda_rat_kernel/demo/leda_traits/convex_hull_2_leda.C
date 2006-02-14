#include <CGAL/basic.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Polygon_2.h>
#include <LEDA/window.h>
#include <LEDA/rat_window.h>
#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>
#include <list>

#undef list

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

typedef leda_rat_point                               Point_2;

int main()
{
  std::list<Point_2> input;

  leda_window win(500,500,"convex hull of points");
  win.init(0,500,0);
  
  win.display();
  
  Point_2 act, prev, first;
  
  while(win >> act){
   win << act;
   input.push_back(act);
  }

  std::list<Point_2> output;

  CGAL::leda_rat_kernel_traits  tr;

  CGAL::convex_hull_2(input.begin(), input.end(), std::back_inserter(output), tr );
		       
  // draw the output
  if (output.size() != 0) {
   std::list<Point_2>::const_iterator it = output.begin();
   prev = *it;
   first = prev;
    
   for(;it != output.end(); it++){
     win.draw_segment(prev.to_float(), (*it).to_float(), leda_red);
     prev = *it;
   }
   
   win.draw_segment(first.to_float(), prev.to_float(), leda_red);
  }		       
		    
  win.read_mouse();          // wait for mouse click in window
  
  return 0;
}
