
//#if CATEGORY == 1  //left category
//#else
#if CATEGORY == 2  //reflect category
#define HAS_LEFT_NOT 
#define HAS_REFLECT
#else
#if CATEGORY == 3  //no category
#define HAS_LEFT_NOT
#endif
#endif
//#endif

#include <iostream>
#include <CGAL/basic.h>
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Pm_segment_traits_2.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/assertions.h>
#include <CGAL/Polygon_2.h>
#include <set> 
#include <CGAL/Planar_map_2/Pm_traits_wrap_2.h>

typedef CGAL::Quotient<CGAL::MP_Float>          NT;
typedef CGAL::Cartesian<NT>                     R;
typedef CGAL::Pm_segment_traits_2<R>             Traits;
typedef CGAL::Pm_default_dcel<Traits>           Dcel;
typedef CGAL::Planar_map_2< Dcel, Traits >      Planar_map;
typedef CGAL::Pm_traits_wrap_2<Traits>          Traits_wrap;
typedef Traits::Point_2                         Point_2;
typedef Traits::X_monotone_curve_2              X_monotone_curve_2;

int main()
{
  std::cout << "Test categories" << std::endl;

#ifdef HAS_LEFT_NOT
  std::cout << "HAS_LEFT_NOT" << std::endl;
#endif

#ifdef HAS_REFLECT
  std::cout << "HAS_REFLECT" << std::endl;
#endif
  
  //create wrapper 
  Traits_wrap Tw;

  //create planar map
  Planar_map Pm;

  //create points:
  Point_2 a1(1, 3), a2(4, 3), a3(3, 5);

  // Create the curves:
  X_monotone_curve_2 cv[2];

  cv[0]  = X_monotone_curve_2(a1, a2);
  cv[1]  = X_monotone_curve_2(a2, a3);

  //curves_compare_y_at_x_left test - for categoris
  CGAL::Comparison_result l_res = Tw.curves_compare_y_at_x_left (cv[0], cv[1], a2);
  std::cout << "l-res = " << l_res <<  std::endl;  
  l_res = Tw.curves_compare_y_at_x_left (cv[1], cv[0], a2);
  std::cout << "l-res = " << l_res <<  std::endl;  
  l_res = Tw.curves_compare_y_at_x_left (cv[0], cv[0], a2);
  std::cout << "l-res = " << l_res <<  std::endl; 

  //check if the results are o.k.
  

  return 0;
}
