#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Timer.h>
#include <iostream>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = Kernel::Point_3;
using Segment_3 = Kernel::Segment_3;
using Triangle_3 = Kernel::Triangle_3;

int main()
{
  Point_3 p(0, 0, 0);
  Point_3 q(10, 0, 0);
  Point_3 r(0, 10, 0);

  Point_3 s(1, 1, -1);
  Point_3 t(1, 1, 1);
  Point_3 u(1, 1, -2);

  Segment_3 st(s, t);
  Segment_3 tu(t, u);
  Triangle_3 tri(p, q, r);

  int does = 0;
  int does_not = 0;

   CGAL::Timer timer;
   timer.start();


  for (int i = 0; i < 1000000000; ++i) {
    int j = i/3;
    if(i == 3*j){
      if(CGAL::do_intersect(st, tri)){
        ++does;
      }
    }else{
      if(CGAL::do_intersect(tu, tri)){
        ++does_not;
      }
    }
  }

   timer.stop();
   std::cout << "Does intersect: " << does << std::endl;
   std::cout << "Does not intersect: " << does_not << std::endl;
   std::cout << "Time: " << timer.time() << " seconds." << std::endl;

  return 0;
}
