
#include <CGAL/basic.h>
#include <CGAL/convex_hull_2.h>
#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>
#include <LEDA/random_rat_point.h>

#include <list>
#include <iostream> 

#if defined(LEDA_NAMESPACE)
using namespace leda;
#endif

typedef CGAL::leda_rat_kernel_traits                           K;
typedef K::Point_2                                             Point;



int main()
{
  std::list<Point>          input;
  
  // insert some points  
  int i;
  Point p;
  
  for(i=0;i<10000;i++){
     random_point_in_square(p, 5000);
     input.push_back(p);
  }
  
  K  tr;
  std::list<Point> output;
  CGAL::convex_hull_2(input.begin(), input.end(), std::back_inserter(output), tr ); 
  
  return 0;
}
