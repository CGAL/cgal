// examples/Arrangement_2/example3.C   
// ---------------------------------

#include "short_names.h"

#include <CGAL/basic.h>

#ifndef CGAL_USE_LEDA
// To enable compilation without leda:
int main ()
{
  return (0);
}

#else

#include <CGAL/Cartesian.h>
#include <CGAL/leda_real.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arr_conic_traits_2.h> 
#include <CGAL/Arrangement_2.h>

typedef leda_real                                       NT;
typedef CGAL::Cartesian<NT>                             Kernel;
typedef CGAL::Arr_conic_traits_2<Kernel>                Traits;
typedef Traits::Point_2                                 Point_2;
typedef Traits::Circle_2                                Circle_2;
typedef Traits::Curve_2                                 Curve_2;
typedef Traits::X_monotone_curve_2                      X_monotone_curve_2;
typedef CGAL::Arr_2_default_dcel<Traits>                Dcel;
typedef CGAL::Arrangement_2<Dcel,Traits>                Arr_2;

int main()
{
  Arr_2 arr;  
              
  // 2 ccw circles with radius 5 and center (0,0) and (6,0) resp.
  Circle_2  c1 (Point_2(0,0), 5*5, CGAL::COUNTERCLOCKWISE);
  Circle_2  c2 (Point_2(6,0), 5*5, CGAL::COUNTERCLOCKWISE);

  Arr_2::Curve_iterator cit = arr.insert(c1);
  cit = arr.insert(c2); 

  // upward vertical ray shooting
  Arr_2::Locate_type lt;
  Arr_2::Halfedge_handle e = arr.vertical_ray_shoot(Point_2(-1, 0), lt, true);

  CGAL_assertion(e->source()->point() == Point_2(3, 4)); 
  CGAL_assertion(e->target()->point() == Point_2(-5, 0));

  return 0;
}

#endif
