// see examples/Planar_map/example1.C
// for original version ...
// ------------------------------

// needed in Pm_segment_traits_2 - otherwise
// we use x() / y()
#define NBUG
#define LEDA_NO_MIN_MAX_TEMPL

#define CGAL_PROVIDE_LEDA_RAT_KERNEL_TRAITS_3

#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>
#include <CGAL/Kernel_checker.h>
#include <CGAL/Pm_segment_traits_2.h>
#include <CGAL/Pm_segment_traits_leda_kernel_2.h>

#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>
#include <iostream>
#include <iterator>
#include <algorithm>

//typedef CGAL::leda_rat_kernel_traits      K1;
//typedef CGAL::Homogeneous<leda_integer>   K2;

typedef CGAL::leda_rat_kernel_traits      Kernel;
//typedef CGAL::Pm_segment_traits_leda_kernel_2<leda_rational> Kernel;
//typedef CGAL::Kernel_checker<K1, K2, CGAL::leda_to_cgal_2 > Kernel;

typedef CGAL::Pm_segment_traits_2<Kernel> Traits;
typedef Traits::Point_2                   Point_2;
typedef Traits::X_curve_2                 X_curve_2; // segment ...
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
  std::cout << "returned the curve " << e->curve() << ", oriented toward " 
            << e->target()->point() << std::endl; 
  return 0;
}

