#include "short_names.h"

#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/Pm_segment_traits_2.h>
#include <CGAL/IO/Pm_iostream.h>
#include <iostream>
#include <CGAL/Timer.h>

typedef CGAL::Quotient<int>                     NT;
typedef CGAL::Cartesian<NT>                     Kernel;
typedef CGAL::Pm_segment_traits_2<Kernel>       Traits;
typedef CGAL::Pm_default_dcel<Traits>           Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>         Planar_map;

int main()
{
  
    // X_curve_2 cv[5];
  
  CGAL::Timer timer;
  timer.start();
  int i;
  for (i = 0; i < 1000; i++) {
    Planar_map pm;
    std::cin >> pm;
  }
  timer.stop();

  std::cout << "Total insertion time: " <<  timer.time() << std::endl;
  
  return 0;
}
