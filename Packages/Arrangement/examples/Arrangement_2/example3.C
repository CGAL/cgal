// examples/Arrangement_2/example3.C   
// ---------------------------------

#include "short_names.h"

#include <CGAL/basic.h>
#include <CGAL/Arr_2_bases.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arr_circles_real_traits.h> 
#include <CGAL/Arrangement_2.h>

// We use here double instead of leda_real to enable compilation without LEDA.
// This is not recommended generally.
// Read more in the README file or in the manual.

typedef double                                          NT;
typedef CGAL::Arr_circles_real_traits<NT>               Traits;
typedef Traits::Point_2                                 Point_2;
typedef Traits::X_curve_2                               X_curve_2;
typedef Traits::Curve_2                                 Curve_2;
typedef CGAL::Arr_base_node<Curve_2>                    Base_node;
typedef CGAL::Arr_2_default_dcel<Traits>                Dcel;
typedef CGAL::Arrangement_2<Dcel,Traits,Base_node>      Arr_2;

int main()
{
  Arr_2 arr;  

  // 2 ccw circles with radius 5 and center (0,0) and (6,0) resp.
  Arr_2::Curve_iterator cit = arr.insert(Curve_2(0, 0, 25));
  cit=arr.insert(Curve_2(6, 0, 25)); 

  // upward vertical ray shooting
  Arr_2::Locate_type lt;
  Arr_2::Halfedge_handle e = arr.vertical_ray_shoot(Point_2(-1, 0), lt, true);

  CGAL_assertion(e->source()->point() == Point_2(3, 4)); 
  CGAL_assertion(e->target()->point() == Point_2(-5, 0));

  return 0;
}
