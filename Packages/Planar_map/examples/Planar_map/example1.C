// examples/Planar_map/example1.C
// ------------------------------

#include "short_names.h"

#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/Pm_segment_traits_2.h>
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>
#include <iostream>
#include <iterator>
#include <algorithm>

typedef CGAL::Quotient<long>              Number_type;
typedef CGAL::Cartesian<Number_type  >    Kernel;
typedef CGAL::Pm_segment_traits_2<Kernel> Traits;
typedef Traits::Point_2                   Point_2;
typedef Traits::X_curve_2                 X_curve_2;
typedef CGAL::Pm_default_dcel<Traits>     Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>   Planar_map;

int main()
{
  // Create an instance of a Planar_map:
  Planar_map pm;
  X_curve_2 cv[5];

  Point_2 a1(100, 0), a2(20, 50), a3(180, 50), a4(100, 100);

  // Create the curves:
  cv[0] = X_curve_2(a1, a2);
  cv[1] = X_curve_2(a1, a3);
  cv[2] = X_curve_2(a2, a3);
  cv[3] = X_curve_2(a2, a4);
  cv[4] = X_curve_2(a3, a4);

  std::cout << "The curves of the map :" << std::endl;
  std::copy(&cv[0], &cv[5], std::ostream_iterator<X_curve_2>(std::cout, "\n"));
  std::cout << std::endl;

  // Insert the curves into the Planar_map:
  std::cout << "Inserting the curves to the map ... ";
  pm.insert(&cv[0], &cv[5]);
  std::cout << ((pm.is_valid()) ? "map valid!" : "map invalid!") << std::endl
            << std::endl;
  
  // Shoot a vertical ray upward from p:
  Point_2 p(95, 30);
  Planar_map::Locate_type lt;

  std::cout << "Upward vertical ray shooting from " << p << std::endl; 
  Planar_map::Halfedge_handle e = pm.vertical_ray_shoot(p, lt, true);
  std::cout << "returned the curve " << e->curve() <<  ", oriented toward " 
            << e->target()->point() << std::endl; 
  return 0;
}

